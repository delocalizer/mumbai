#!/usr/bin/env python3
"""
Query a region across multiple bams.

Default log level is 'INFO' — set to something else with the LOG_LEVEL
environment variable
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import importlib.resources as pkg_resources
import logging
from math import floor
import os
from pathlib import Path
import pathlib
import sqlite3
import sys
import traceback
import yaml

import psycopg2
from psycopg2.extras import LoggingConnection
from pysam import AlignedSegment, AlignmentFile, AlignmentHeader

from mumbai.mumbai_db import valid_dbid, ERR_UNEXPECTED, MSG_DSN_FORMATS
import mumbai.resources

LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(os.environ.get('LOG_LEVEL', 'INFO'))
LOGGER.addHandler(logging.StreamHandler())

BIN_MAX_RNG = (2**29) - 1
BIN_ID_STARTS = (1, 9, 73, 585, 4681)
MODE = ('count', 'sam', 'tview')
WINDOW_SIZE = 2**14

DEFAULT_MAX_WORKERS = 20

MSG_NBAMS = '%s bams queried'
MSG_NRECORDS = '%s records'
MSG_REGION = '%s:%s-%s'

with pkg_resources.open_text(mumbai.resources, 'config.yaml') as fh:
    CONF = yaml.full_load(fh)
SQL = CONF['SQL']


# BAI0 is a well-formatted but empty index file. pysam.AlignmentFile
# constructor is 'helpful' and by default tries to open any canonically named
# bai file adjacent to the specified bam. We wish to avoid this (all our bams
# have adjacent indexes, and not reading them is the whole point of this
# exercise) so we use the index_filename constructor option to specify a
# trivial index file which is quickly read but never used. Longer term
# solution would be to patch AlignmentFile with an option to force it not to
# use an index at all.
BAI0 = pkg_resources.path(mumbai.resources, '0_index')


def aligned_view(alignedseg):
    """
    Returns a 'tview' style representation of the AlignedSegment sequence as
    a list of chars inferred from the CIGAR and the query sequence
    """
    rep, pos, seq = [], 0, alignedseg.query_sequence
    for (op, num) in alignedseg.cigartuples:
        if op == 0:    # M
            if alignedseg.is_reverse:
                rep.extend(map(str.lower, seq[pos:pos+num]))
            else:
                rep.extend(seq[pos:pos+num])
            pos += num
        elif op == 1:  # I
            rep[-1] = '|'
            pos += num
        elif op == 2:  # D
            rep.extend(['*'] * num)
        elif op == 3:  # N
            rep.extend(['-'] * num)
        elif op == 4:  # S
            rep.extend(['.'] * num)
            pos += num
        elif op == 5:  # H
            rep.extend([' '] * num)
            pos += num
        else:
            raise NotImplementedError(op)
    return rep


def bams_offsets(dbid, ref, start, stop):
    """
    Query the bam index database at the region 'ref:start-stop' and return a
    list of (bam, offset) tuples, where offset is the bgzip offset to start
    reading from the bam to obtain records that overlap the region.

    Args:
        dbid: a dictionary of postgres database connection parameters OR
            the path to a sqlite3 database
        ref: Query interval reference.
        start: Query interval start position (1-based inclusive).
        stop: Query interval stop position (1-based inclusive).
    """
    dbtype = 'postgres' if isinstance(dbid, dict) else 'sqlite3'

    window = floor((start-1)/WINDOW_SIZE)
    bins = tuple(reg2bins(start, stop))
    query = SQL[dbtype]['select_offset']

    if dbtype == 'postgres':
        with psycopg2.connect(connection_factory=LoggingConnection,
                              **dbid) as con:
            con.initialize(LOGGER)
            with con.cursor() as cur:
                cur.execute(query, (ref, window, bins))
                offsets = cur.fetchall()
    else:
        with sqlite3.connect(dbid) as con:
            con.set_trace_callback(LOGGER.debug)
            cur = con.cursor()
            # workaround sqlite3 DB-API not supporting placeholder for arrays
            cur.execute(query % ','.join('?' for b in bins),
                        (ref, window, *bins))
            offsets = cur.fetchall()
    return offsets


def bam_query(bam, offset, refnames, ref, start, stop):
    """
    Return mapped records that overlap the query region.

    Reference names are passed explicitly to avoid a disk seek to read the bam
    header but still end up with reference names and not numeric ids in the
    SAM record RNAME field.

    Args:
        bam: Path to the coordinate-sorted bam.
        offset: Start reading the bam from here.
        refnames: Names of references in the header
        ref: Query interval reference.
        start: Query interval start position (1-based inclusive).
        stop: Query interval stop position (1-based inclusive).
    Returns:
        serialized SAM format records.
    """
    # NB: pysam uses 0-based coords so convert
    start -= 1
    stop -= 1

    # Don't read the header as that requires a seek
    bamf = AlignmentFile(bam, mode='rb', check_header=False, check_sq=False,
                         index_filename=BAI0, reference_names=refnames)
    bamf.seek(offset)
    records = []
    for read in bamf:

        read_start, read_stop = read.pos, read.pos + len(read.seq) - 1

        # Read on the wrong contig => stop
        if read.reference_name != ref:
            break

        # Read too far right => stop
        if read_start > stop:
            break

        # Read is too far left => skip
        if read_stop < start:
            continue

        # Read is unmapped or has empty sequence (thanks, cutadapt) => skip
        if read.is_unmapped or not read.seq:
            continue

        # Read originates outside but overlaps or read originates inside region
        assert read_start < start <= read_stop or start <= read_start <= stop

        # serialize AlignedSegment for pickling
        # https://github.com/pysam-developers/pysam/issues/950#issuecomment-737361642
        records.append(read.to_string())

    return records


def bams_query(dbid, ref, start, stop, max_workers):
    """
    Return mapped records that overlap the query region from all bams in the
    index database.

    Args:
        dbid: a dictionary of postgres database connection parameters OR
            the path to a sqlite3 database
        ref: Query interval reference.
        start: Query interval start position (1-based inclusive).
        stop: Query interval stop position (1-based inclusive).
        max_workers: Maximum number of processes to spawn for multiprocessing.
    Returns:
        (header, [SAM format records])
    """
    offsets = bams_offsets(dbid, ref, start, stop)
    LOGGER.info(MSG_NBAMS, len(offsets))
    if not offsets:
        return (None, [])

    # header template from the first bam we see
    h0 = AlignmentFile(offsets[0][0], index_filename=BAI0).header.to_dict()
    header = AlignmentHeader.from_dict({
        'HD': h0['HD'],
        'SQ': h0['SQ'],
        'PG': [{'ID': Path(sys.argv[0]).name, 'CL': ' '.join(sys.argv[1:])}]
    })

    batch_args = [
        (bam, offset, header.references, ref, start, stop)
        for bam, offset in offsets]
    records = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        for batch in executor.map(bam_query, *zip(*batch_args)):
            records.extend(batch)
    return (header, records)


def parse_cmdargs(args):
    """
    Returns: parsed arguments
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        'mode',
        choices=MODE,
        help='How to display the results.')
    parser.add_argument(
        'dbid',
        metavar='DSN|PATH',
        type=valid_dbid,
        help='The bam index database to query. If a valid libpq connection '
        f'string is supplied e.g. {MSG_DSN_FORMATS} then a postgres database '
        'is inferred, otherwise the value is interpreted as the path to a '
        'sqlite3 database file.')
    parser.add_argument(
        'ref',
        help='Query region reference sequence name e.g. chr1.')
    parser.add_argument(
        'start',
        type=int,
        help='Query region start pos (1-based, inclusive).')
    parser.add_argument(
        'stop',
        type=int,
        help='Query region stop pos (1-based, inclusive).')
    parser.add_argument(
        '--max-workers',
        type=int,
        default=DEFAULT_MAX_WORKERS,
        help=f'max processes to use (default={DEFAULT_MAX_WORKERS})')

    return parser.parse_args(args)


def reg2bins(rbeg, rend):
    """
    Generate bin ids which overlap the specified region.

    Args:
        rbeg (int): 1-based inclusive beginning position of region
        rend (int): 1-based inclusive end position of region
    Yields:
        (int): bin IDs for overlapping bins of region
    Raises:
        AssertionError (Exception): if the range is malformed or invalid

    Credit: https://github.com/betteridiot/bamnostic
    """
    # 1-based inputs to 0-based calculation
    rbeg -= 1
    rend -= 1

    assert 0 <= rbeg <= rend <= BIN_MAX_RNG, 'Invalid region {}, {}'.format(
        rbeg, rend)

    yield 0
    for start, shift in zip(BIN_ID_STARTS, range(26, 13, -3)):
        i = rbeg >> shift if rbeg > 0 else 0
        j = rend >> shift if rend < BIN_MAX_RNG else BIN_MAX_RNG >> shift

        for bin_id_offset in range(i, j + 1):
            yield start + bin_id_offset


def tview(records, header, histart=0, histop=0):
    """
    Display 'tview' style text view of the SAM records.

    Args:
        records: SAM format text records.
        header: AlignmentHeader to reconstitute AlignedSegment records.
        histart: Start highlighting at this position (1-based inclusive).
        histop: Stop higlighting at this position (1-based inclusive).

    Returns:
        text view
    """
    segs = sorted(
        [AlignedSegment.fromstring(r, header) for r in records],
        key=lambda x: x.reference_start)
    leftmost = segs[0].reference_start
    lines = []
    for seg in segs:
        leftpad = (seg.reference_start - leftmost) * ' '
        av = aligned_view(seg)
        if histart and histop:
            # pysam reference_start (POS) is 0-based so subtract 1
            left = histart - seg.reference_start - 1
            right = histop - seg.reference_start - 1
            for pos in range(left, right+1):
                if 0 <= pos < len(av):
                    av[pos] = '\033[1m' + av[pos] + '\033[0m'
        lines.append(leftpad + ''.join(av))
    return '\n'.join(lines)


def main():
    args = parse_cmdargs(sys.argv[1:])
    LOGGER.info(MSG_REGION, args.ref, args.start, args.stop)

    try:
        header, records = bams_query(args.dbid, args.ref,
                                     args.start, args.stop, args.max_workers)
        if args.mode == 'count':
            print(MSG_NRECORDS % len(records))
        elif args.mode == 'sam':
            recs = sorted(records, key=lambda x: int(x.split('\t')[3]))
            print(str(header) + '\n'.join(recs))
        elif args.mode == 'tview':
            print(tview(records, header, args.start, args.stop))
        else:
            raise NotImplementedError(args.mode)

    except Exception as ex:
        LOGGER.debug(traceback.format_exc())
        LOGGER.error(ERR_UNEXPECTED, ex)


if __name__ == '__main__':
    main()
