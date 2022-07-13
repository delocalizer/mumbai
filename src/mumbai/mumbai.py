#!/usr/bin/env python3
"""
Describe a region across multiple bams.

Default log level is 'INFO' — set to something else with the LOG_LEVEL
environment variable
"""
import argparse
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from importlib import metadata
import importlib.resources as pkg_resources
from itertools import chain
import logging
from math import floor
from operator import itemgetter
import os
from pathlib import Path
import sqlite3
from struct import calcsize, unpack, unpack_from
import struct
import sys
import traceback
import yaml

import psycopg2
from psycopg2.extras import LoggingConnection
from pysam import AlignedSegment, AlignmentFile, AlignmentHeader, BGZFile

from mumbai.mumbai_db import valid_dbid, ERR_UNEXPECTED, MSG_DSN_FORMATS
import mumbai.resources

LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(os.environ.get('LOG_LEVEL', 'INFO'))
LOGGER.addHandler(logging.StreamHandler())

# See: Sequence Alignment/Map Format Specification|4.2.3
BASE_ENC = tuple('=ACMGRSVTWYHKDBN')

# See: Sequence Alignment/Map Format Specification|4.2.4
VALTYPE_ENC = {
    'A': '<c', 'c': '<b', 'C': '<B', 's': '<h', 'S': '<H', 'i': '<i',
    'I': '<I', 'f': '<f'}

BIN_MAX_RNG = (2**29) - 1
BIN_ID_STARTS = (1, 9, 73, 585, 4681)
MODE = ('count', 'pileup', 'sam', 'tview')
WINDOW_SIZE = 2**14

DEFAULT_MAX_WORKERS = 20

ERR_NOTPOSINT = '%s cannot be understood as a positive integer'
MSG_NBAMS = '%s bams queried'
MSG_NRECORDS = '%s records'

with pkg_resources.open_text(mumbai.resources, 'config.yaml') as fh:
    CONF = yaml.full_load(fh)
SQL = CONF['SQL']


# BAI0 is a well-formatted but empty index file. pysam.AlignmentFile
# constructor is 'helpful' and by default tries to open any canonically named
# bai file adjacent to the specified bam. We wish to avoid this (all our bams
# have adjacent indexes, and not reading them is the whole point of this
# exercise) so we use the index_filename constructor option to specify a
# trivial index file which is quickly read but never used. Longer term TODO
# would be to patch AlignmentFile with an option to force it not to use an
# index at all.
BAI0 = pkg_resources.path(mumbai.resources, '0_index')


@dataclass(frozen=True)
class Region:
    """
    A genomic region i.e. contig:[start-stop]

    The coordinates are 1-based and boundaries are included.
    """
    contig: str
    start: int  # 1-based, inclusive
    stop: int   # 1-based, inclusive

    @property
    def start0(self):
        """
        0-based for calculations in pysam.
        """
        return self.start - 1

    @property
    def stop0(self):
        """
        0-based for calculations in pysam.
        """
        return self.stop - 1

    def __str__(self):
        return f'{self.contig}:[{self.start}-{self.stop}]'


def aligned_view(alignedseg):
    """
    Generates a unpadded 'tview' style representation of the aligned portion
    of an AlignedSegment sequence as inferred from the CIGAR and query strings.

    Returns: An unpadded text representation of the alignment (one character
             per reference base).
    """
    aln, qpos, seq = [], 0, alignedseg.query_sequence
    for (op, num) in alignedseg.cigartuples:
        if op == 0:    # M
            if alignedseg.is_reverse:
                aln.extend(map(str.lower, seq[qpos:qpos+num]))
            else:
                aln.extend(seq[qpos:qpos+num])
            qpos += num
        elif op == 1:  # I
            aln[-1] = '|'
            qpos += num
        elif op == 2:  # D
            aln.extend(['-'] * num)
        elif op == 4:  # S
            qpos += num
        elif op == 5:  # H
            pass
        else:
            raise NotImplementedError(op)
    return aln


def bam_fetch(bam, ref_names, ref_lengths, offset, reg, mode):
    """
    Return the records that overlap the region in the bam.

    Args:
        bam: Path to the coordinate-sorted bam.
        ref_names: List of reference names for constructing a header.
        ref_lengths: List of reference lengths for constructing a header.
        offset: Start reading the bam from here.
        reg: Region (contig, start, stop).
        mode: Type of results to return (count, pileup, sam, tview).
    Returns:
        results of the form:
        mode=='count'  => n                    : int
        mode=='pileup' => {pileup}             : dict
        mode=='sam'    => [(pos, SAM record)]  : list
        mode=='tview': => [(pos, alignedview)] : list
    """
    # Implementation note:
    # The price to pay for multiprocessing is that this function has to return
    # pickle-able results; in particular it can't return pysam.AlignedRead
    # instances. We could serialize using .to_string() but then we'd have to
    # reconstitute the SAM records elsewhere with .fromstring() to do
    # calculations for tview and pileup that use reads' reference_positions
    # and reference_sequence. To avoid that round trip cost we have this
    # mildly ugly construction where we do all those calculations here and
    # return different types of results for the different modes.

    bgzf = BGZFile(bam, mode='rb')
    bgzf.seek(offset)
    header = AlignmentHeader.from_references(ref_names, ref_lengths)

    # intialize accumulators
    count = 0
    pileup = {(reg.contig, reg.start + i): defaultdict(int)
              for i in range(reg.stop - reg.start + 1)}
    sam = []
    tview = []

    # See SAM spec section 4.2 "The BAM format"
    while (next_block := bgzf.read(4)):

        block_size = unpack('I', next_block)[0]
        read = parse_bam_record(bgzf.read(block_size), header)

        # Read on the wrong contig => stop walking
        if read.reference_name != reg.contig:
            break

        # Read too far right => stop walking
        if read.pos > reg.stop0:
            break

        # Read is unmapped or has empty sequence (thanks, cutadapt) => skip
        if read.is_unmapped or not read.seq:
            continue

        # Account for any clipping at ends
        aln_start, aln_stop = read.reference_start, read.reference_end - 1

        # alignment has no overlap => skip
        if aln_stop < reg.start0 or aln_start > reg.stop0:
            continue

        if mode == 'count':
            count += 1
        elif mode == 'pileup':
            update_pileup(pileup, reg, read)
        elif mode == 'sam':
            sam.append((read.pos, read.to_string()))
        elif mode == 'tview':
            tview.append((aln_start, aligned_view(read)))
        else:
            raise NotImplementedError(mode)

    if mode == 'count':
        return count
    if mode == 'pileup':
        return pileup
    if mode == 'tview':
        return tview
    if mode == 'sam':
        return sam


def bams_fetch(dbid, reg, mode, max_workers=1):
    """
    Return the records that overlap the region from all bams in the index
    database.

    Args:
        dbid: A dictionary of postgres database connection parameters OR
            the path to a sqlite3 database.
        reg: Region (contig, start, stop).
        mode: Type of results to return (count, pileup, sam, tview).
        max_workers: Maximum number of processes to spawn for multiprocessing.
    Returns:
        mode=='count'  =>  n                       : int
        mode=='pileup' =>  pileup                  : dict
        mode=='sam'    =>  (header, [SAM records])
        mode=='tview'  =>  [tview]                 : list
    """
    offsets = bams_offsets(dbid, reg)
    LOGGER.info(MSG_NBAMS, len(offsets))
    if not offsets:
        return None

    # header template from the first bam we see
    h0 = AlignmentFile(offsets[0][0], index_filename=BAI0).header.to_dict()
    header = AlignmentHeader.from_dict({
        'HD': h0['HD'],
        'SQ': h0['SQ'],
        'PG': [{'ID': Path(sys.argv[0]).name, 'CL': ' '.join(sys.argv[1:])}]
    })

    batch_args = [(bam, header.references, header.lengths, offset, reg, mode)
                  for bam, offset in offsets]
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        batches = executor.map(bam_fetch, *zip(*batch_args))
    return _agg_results(batches, header, mode)


def bams_offsets(dbid, reg):
    """
    Query the bam index database at the region and return a list of
    (bam, offset) tuples, where offset is the bgzip offset to begin reading
    from the bam to obtain records that overlap the region.

    Args:
        dbid: A dictionary of postgres database connection parameters OR
            the path to a sqlite3 database.
        reg: Region (contig, start, stop).
    """
    dbtype = 'postgres' if isinstance(dbid, dict) else 'sqlite3'

    window = floor((reg.start0)/WINDOW_SIZE)
    bins = tuple(reg2bins(reg.start0, reg.stop0))
    query = SQL[dbtype]['select_offset']

    if dbtype == 'postgres':
        with psycopg2.connect(connection_factory=LoggingConnection,
                              **dbid) as con:
            con.initialize(LOGGER)
            with con.cursor() as cur:
                cur.execute(query, (reg.contig, window, bins))
                offsets = cur.fetchall()
    else:
        with sqlite3.connect(dbid) as con:
            con.set_trace_callback(LOGGER.debug)
            cur = con.cursor()
            # workaround sqlite3 DB-API not supporting placeholder for arrays
            cur.execute(query % ','.join('?' for b in bins),
                        (reg.contig, window, *bins))
            offsets = cur.fetchall()
    return offsets


def parse_cmdargs(args):
    """
    Returns: Parsed arguments.
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
        help='Genomic region reference sequence name e.g. chr1.')
    parser.add_argument(
        'start',
        type=valid_pos_int,
        help='Genomic region start pos (1-based, inclusive).')
    parser.add_argument(
        'stop',
        type=valid_pos_int,
        help='Genomic region stop pos (1-based, inclusive).')
    parser.add_argument(
        '--max-workers',
        type=valid_pos_int,
        default=DEFAULT_MAX_WORKERS,
        help=f'Max processes to use (default={DEFAULT_MAX_WORKERS})')
    parser.add_argument(
        '--version',
        action='version',
        version=metadata.version('mumbai'))
    parsed = parser.parse_args(args)
    parsed.reg = Region(parsed.ref, parsed.start, parsed.stop)
    return parsed


def parse_bam_record(byts, header):
    """
    Return AlignedSegment parsed from binary BAM record.

    byts: The BAM record.
    header: AlignmentHeader containing reference sequence names.

    See: Sequence Alignment/Map Format Specification|4.2
    """
    read = AlignedSegment(header)
    fmt, offset = '<l', 0
    read.reference_id = unpack_from(prev := fmt, byts, offset)[0]

    fmt, offset = '<l', offset + calcsize(prev)
    read.reference_start = unpack_from(prev := fmt, byts, offset)[0]

    fmt, offset = '<B', offset + calcsize(prev)
    l_read_name = unpack_from(prev := fmt, byts, offset)[0]

    fmt, offset = '<B', offset + calcsize(prev)
    read.mapping_quality = unpack_from(prev := fmt, byts, offset)[0]

    fmt, offset = '<H', offset + calcsize(prev)
    bin_ = unpack_from(prev := fmt, byts, offset)[0]

    fmt, offset = '<H', offset + calcsize(prev)
    n_cigar_op = unpack_from(prev := fmt, byts, offset)[0]

    fmt, offset = '<H', offset + calcsize(prev)
    read.flag = unpack_from(prev := fmt, byts, offset)[0]

    fmt, offset = '<L', offset + calcsize(prev)
    l_seq = unpack_from(prev := fmt, byts, offset)[0]

    fmt, offset = '<l', offset + calcsize(prev)
    read.next_reference_id = unpack_from(prev := fmt, byts, offset)[0]

    fmt, offset = '<l', offset + calcsize(prev)
    read.next_reference_start = unpack_from(prev := fmt, byts, offset)[0]

    fmt, offset = '<l', offset + calcsize(prev)
    read.template_length = unpack_from(prev := fmt, byts, offset)[0]

    fmt, offset = f'<{l_read_name}s', offset + calcsize(prev)
    read.query_name = unpack_from(prev := fmt, byts, offset)[0][:-1].decode('ascii')

    fmt, offset = f'<{n_cigar_op}L', offset + calcsize(prev)
    cigar_enc = unpack_from(prev := fmt, byts, offset)
    read.cigartuples = [(c & 15, c >> 4) for c in cigar_enc]

    fmt, offset = f'<{(l_seq+1)//2}B', offset + calcsize(prev)
    seq_enc = unpack_from(prev := fmt, byts, offset)
    seq = ''.join(chain.from_iterable(
        (BASE_ENC[i >> 4], BASE_ENC[i & 15]) for i in seq_enc))
    # omit undefined last value when l_seq is odd
    read.query_sequence = seq[:-1] if l_seq % 2 else seq

    fmt, offset = f'<{l_seq}s', offset + calcsize(prev)
    read.query_qualities = unpack_from(prev := fmt, byts, offset)[0]

    offmax = len(byts)
    fmt, offset = '<2s', offset + calcsize(prev)
    while offset < offmax :
        tag = unpack_from(prev := fmt, byts, offset)[0].decode('ascii')
        fmt, offset = '<c', offset + calcsize(prev)
        val_type = unpack_from(prev := fmt, byts, offset)[0].decode('ascii')
        offset += calcsize(prev)
        if val_type in VALTYPE_ENC:
            fmt = f'{VALTYPE_ENC[val_type]}'
            val = unpack_from(prev := fmt, byts, offset)[0]
        elif val_type == 'Z':
            remain = byts[offset:]
            nul = remain.index(b'\x00')
            val = ''.join(chr(b) for b in remain[:nul])
            prev = f'<{len(val)+1}c'
        else:
            raise NotImplementedError(val_type)
        read.set_tag(tag, val, val_type)
        fmt, offset = '<2s', offset + calcsize(prev)

    return read


def reg2bins(rbeg, rend):
    """
    Generate bin ids which overlap the specified region.

    Args:
        rbeg (int): 0-based inclusive beginning position of region
        rend (int): 0-based inclusive end position of region
    Yields:
        (int): bin IDs for overlapping bins of region
    Raises:
        AssertionError (Exception): if the range is malformed or invalid

    Credit: https://github.com/betteridiot/bamnostic
    """
    assert 0 <= rbeg <= rend <= BIN_MAX_RNG, f'Invalid region {rbeg}, {rend}'

    yield 0
    for start, shift in zip(BIN_ID_STARTS, range(26, 13, -3)):
        i = rbeg >> shift if rbeg > 0 else 0
        j = rend >> shift if rend < BIN_MAX_RNG else BIN_MAX_RNG >> shift

        for bin_id_offset in range(i, j + 1):
            yield start + bin_id_offset


def update_pileup(pileup, reg, read):
    """
    Update the region pileup dictionary from the read.

    Args:
        pileup: Pileup dictionary to be updated in-place.
        reg: Genomic region.
        read: AlignedSegment that overlaps the region.
    Returns:
        the updated pileup.
    """
    ref_bases = read.get_reference_sequence()  # requires MD tag
    aln_start, aln_stop = read.reference_start, read.reference_end - 1
    aln = aligned_view(read)
    aln_0 = 0 if reg.start0 <= aln_start else reg.start0 - aln_start
    aln_1 = len(aln) if reg.stop0 >= aln_stop else reg.stop0 - aln_stop + len(aln)
    for aln_i in range(aln_0, aln_1):
        pos = aln_start + aln_i + 1
        base = aln[aln_i].upper()
        p_loc = pileup[(reg.contig, pos)]
        p_loc['ref'] = ref_bases[aln_i]
        p_loc[base] += 1
    return pileup


def valid_pos_int(arg):
    """
    Validate arg as +ve integer.
    """
    try:
        ival = int(arg)
        if ival <= 0:
            raise ValueError
    except ValueError as valerr:
        raise argparse.ArgumentTypeError(ERR_NOTPOSINT % arg) from valerr
    return ival


def _agg_results(results, header, mode):
    """
    Aggregate results according to their type.

    Args:
        results: Iterable of outputs from bam_fetch.
        header: AlignmentHeader
        mode: Type of results (count, pileup, sam, tview).
    Returns:
        results of the form:
        mode=='count'  => n                        : int
        mode=='pileup' => {pileup}                 : dict
        mode=='sam'    => (header, [SAM records])  : list
        mode=='tview': => [(pos, alignedview)]     : list
    """
    if mode == 'count':
        return sum(results)
    if mode == 'pileup':
        agg = defaultdict(lambda: defaultdict(int))
        for pileup in results:
            for loc in pileup:
                agg[loc]['ref'] = pileup[loc]['ref']
                for base in ('A', 'C', 'G', 'T', '-', '|'):
                    agg[loc][base] += pileup[loc].get(base, 0)
        return agg
    if mode == 'sam':
        records = [sam for pos, sam in sorted(chain.from_iterable(results),
                                              key=itemgetter(0))]
        return (header, records)
    if mode == 'tview':
        return list(sorted(chain.from_iterable(results), key=itemgetter(0)))


def _fmt_results(results, reg, mode):
    """
    Format results for output according to their type.

    Args:
        results: Aggregated outputs from bams_fetch.
        reg: Region (contig, start, stop).
        mode: Type of results to format (count, pileup, sam, tview).
    Returns:
        str
    """
    if mode == 'count':
        total = results
        return MSG_NRECORDS % total
    if mode == 'pileup':
        pileup = results
        lines = []
        for (contig, pos), pup in pileup.items():
            counts = [(key, pup.get(key, 0))
                    for key in ('A', 'C', 'G', 'T', '-', '|')]
            # order by counts desc
            freqs = sorted(counts, key=lambda x: x[1], reverse=True)
            lines.append('{}\t{}\t{}\t{:>9}{:>9}{:>9}{:>9}{:>9}{:>9}'.format(
                contig, pos, pup['ref'], *[f'{f[0]}:{f[1]}' for f in freqs]))
        return '\n'.join(lines)
    if mode == 'sam':
        header, records = results
        return str(header) + '\n'.join(records)
    if mode == 'tview':
        alignments = results
        leftmost = alignments[0][0]
        lines = []
        for ref_start, aln in alignments:
            leftpad = (ref_start - leftmost) * ' '
            ref_stop = ref_start + len(aln) - 1
            # region highlighting
            aln_0 = 0 if reg.start0 <= ref_start else reg.start0 - ref_start
            aln_1 = -1 if reg.stop0 >= ref_stop else reg.stop0 - ref_stop + len(aln) - 1
            aln[aln_0] = '\033[1m' + aln[aln_0]
            aln[aln_1] = aln[aln_1] + '\033[0m'
            lines.append(leftpad + ''.join(aln))
        return '\n'.join(lines)
    raise NotImplementedError(mode)


def main():
    args = parse_cmdargs(sys.argv[1:])
    LOGGER.info(args.reg)
    try:
        results = bams_fetch(args.dbid, args.reg, args.mode, args.max_workers)
        print(_fmt_results(results, args.reg, args.mode))
    except Exception as ex:
        LOGGER.debug(traceback.format_exc())
        LOGGER.error(ERR_UNEXPECTED, ex)


if __name__ == '__main__':
    main()
