#!/usr/bin/env python3
"""
Create and load a bam index database.
"""
import argparse
import logging
from importlib import metadata
import importlib.resources as pkg_resources
import io
import os
from os.path import abspath, basename, exists, expanduser
import sqlite3
import string
import struct
import sys
import traceback
import uuid
import yaml

import psycopg2
from psycopg2.extensions import parse_dsn, ISOLATION_LEVEL_READ_UNCOMMITTED
from psycopg2.sql import Identifier, SQL as pgSQL
import pysam

from . import resources


LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(os.environ.get('LOG_LEVEL', 'INFO'))
LOGGER.addHandler(logging.StreamHandler())

MAGICBAI = 'BAI\1'
MAGICSQLITE = 'SQLite format 3\0'

ERR_BAD_MAGIC = 'the file %s is not a sqlite3 database'
ERR_DB_EXISTS = 'the file %s already exists.'
ERR_INVALID_INPUTFILE = 'the file %s does not exist or is not readable'
ERR_INVALID_DBID = ('"%s" is neither a valid libpq connection string '
                    'nor a sqlite3 database file.')
ERR_METADATA_CHUNKS = 'Expected 2 chunks in metadata pseudo-bin, got %s'
ERR_NOT_EOF = 'not at EOF'
ERR_UNEXPECTED = 'Unexpected error: %s'

MSG_BAM_DUPLICATE = '%s already processed.'
MSG_COUNTS = '%s bam indexes processed.'
MSG_DB_CREATED = '%s database %s created'
MSG_DSN_FORMATS = ('\n"host=myhost dbname=mydb user=myuser"\n'
                   'or\n"postgresql://myuser@myhost/mydb"')
MSG_PROCESSING = 'processing %s...'
MSG_REF_MISMATCH = ('@SQ SN LN values for %s are not a subset of the values '
                    'from the database ref table:\n%s\nWas it '
                    'aligned to a different reference?')

# resources from adjacent files
with pkg_resources.open_text(resources, 'schema.sql') as fh:
    SCHEMA_TEMPLATE = string.Template(fh.read())
with pkg_resources.open_text(resources, 'config.yaml') as fh:
    CONF = yaml.full_load(fh)
SQL = CONF['SQL']
REFERENCE = CONF['reference']


def copy_stringio(cur, tablename, data):
    """
    Efficiently copy data into a postgres table

    See https://hakibenita.com/fast-load-data-python-postgresql#copy

    Args:
        con: Database connection.
        tablename: Name of the destination table.
        data: Rows (iterable of iterables).
    """
    sep = '|'
    newline = '\n'
    csv = io.StringIO()
    for row in data:
        csv.write(sep.join(map(str, row)) + newline)
    csv.seek(0)
    cur.copy_from(csv, tablename, sep=sep)


def create_baidb(newdb, reference, dsn=None):
    """
    Create a new bam index database.

    Args:
        newdb: Name for the new bam index database. If dsn is not supplied the
            value is interpreted as the path to the new sqlite3 database file.
        reference: {GRCh37, GRCh38} — reference genome used by indexed bams.
        dsn: A dictionary of connection parameters. If this argument is
            supplied a postgres database will be created, otherwise a sqlite3
            file database.
    """
    dbtype = 'postgres' if dsn else 'sqlite3'

    # explicitly create postgres db
    if dbtype == 'postgres':
        con = psycopg2.connect(**dsn)
        con.autocommit = True
        cur = con.cursor()
        cur.execute(pgSQL(SQL[dbtype]['create_db']).format(
                    Identifier(newdb),
                    Identifier(dsn['user'])))
        con.close()
    else:
        newdb = abspath(expanduser(newdb))

    # DDL as per dbtype
    create_schema = SCHEMA_TEMPLATE.substitute(SQL[dbtype]['coltypes'])
    LOGGER.debug(create_schema)

    # ref table contents as per reference
    refs = [(i, sn, ln) for i, (sn, ln) in enumerate(REFERENCE[reference])]

    if dbtype == 'postgres':
        dsn['dbname'] = newdb
        with psycopg2.connect(**dsn) as con:
            with con.cursor() as cur:
                cur.execute(create_schema)
                cur.executemany(SQL[dbtype]['insert_ref'], refs)
    else:
        with sqlite3.connect(newdb) as con:
            cur = con.cursor()
            cur.executescript(create_schema)
            cur.executemany(SQL[dbtype]['insert_ref'], refs)

    con.close()
    LOGGER.info(MSG_DB_CREATED, dbtype, newdb)


def load_baidb(dbid, bams):
    """
    Parse and load bai files into a bam index database

    Args:
        dbid: a dictionary of postgres database connection parameters OR
            the path to a sqlite3 database
        bams: iterable of bam file paths with adjacent .bai files
    """
    dbtype = 'postgres' if isinstance(dbid, dict) else 'sqlite3'
    if dbtype == 'postgres':
        con = psycopg2.connect(**dbid)
        con.set_isolation_level(ISOLATION_LEVEL_READ_UNCOMMITTED)
        cur = con.cursor()
    else:
        con = sqlite3.connect(f'file:{dbid}?mode=rw', uri=True)
        con.isolation_level = None
        cur = con.cursor()
        # speed up batch inserts
        cur.execute('PRAGMA cache_size=-1048576')  # 1M pages ~ 1G
        cur.execute('PRAGMA journal_mode = OFF;')
        cur.execute('PRAGMA synchronous = 0;')
        cur.execute('PRAGMA cache_size = 1000000;')
        cur.execute('PRAGMA locking_mode = EXCLUSIVE;')
        cur.execute('PRAGMA temp_store = MEMORY;')

    # what reference sequences does the db expect
    cur.execute(SQL['common']['select_ref'])
    db_refs = set((r[1], r[2]) for r in cur.fetchall())

    loaded = 0
    cur.execute('BEGIN;')
    cur.execute(SQL['common']['drop_chunk_idx'])
    cur.execute(SQL['common']['drop_intv_idx'])
    for bam in bams:
        bai = f'{bam}.bai'
        LOGGER.info(MSG_PROCESSING, bai)
        bamuuid = str(uuid.UUID(basename(bai).split('.')[0]))

        af = pysam.AlignmentFile(bam)
        bam_refs = set(zip(af.references, af.lengths))
        if not bam_refs <= db_refs:
            LOGGER.warning(MSG_REF_MISMATCH, bam, sorted(bam_refs - db_refs))
            continue
        cur.execute(SQL[dbtype]['insert_bam'], (bamuuid, bam))
        if cur.rowcount == 0:
            LOGGER.warning(MSG_BAM_DUPLICATE, bam)
            continue
        cur.execute(SQL[dbtype]['last_inserted'])
        bam_id = cur.fetchone()[0]

        chunks, intvs = parse_bai(bam_id, bai)
        if dbtype == 'postgres':
            copy_stringio(cur, 'chunk', chunks)
            copy_stringio(cur, 'intv', intvs)
        elif dbtype == 'sqlite3':
            cur.executemany(SQL[dbtype]['insert_chunk'], chunks)
            cur.executemany(SQL[dbtype]['insert_intv'], intvs)
        loaded += 1

    cur.execute(SQL['common']['create_chunk_idx'])
    cur.execute(SQL['common']['create_intv_idx'])
    cur.execute('COMMIT;')
    con.close()
    LOGGER.info(MSG_COUNTS, loaded)


def parse_bai(bam_id, bai):
    """
    Read chunk and intv data from a bam index file.

    Args:
        bam_id: identifier for the indexed bam
        bai: Path to a bam index file

    Returns:
        (chunks, intvs)
        chunks is a list of records:
            (bam_id, ref_id, bin_id, chunk_beg, chunk_end)
        intvs is a list of records:
            (bam_id, ref_id, window, ioffset)
    """
    chunks = []
    intvs = []
    with open(bai, 'rb') as fh:
        fh.seek(0)
        magic = fh.read(4).decode('ascii')
        assert magic == MAGICBAI
        n_ref = struct.unpack('i', fh.read(4))[0]
        for ref_id in range(n_ref):
            n_bin = struct.unpack('i', fh.read(4))[0]
            for _ in range(n_bin):
                bin_id = struct.unpack('i', fh.read(4))[0]
                n_chunk = struct.unpack('i', fh.read(4))[0]
                if bin_id != 37450:
                    for __ in range(n_chunk):
                        chunk_beg = struct.unpack('l', fh.read(8))[0]
                        chunk_end = struct.unpack('l', fh.read(8))[0]
                        chunks.append(
                            (bam_id, ref_id, bin_id, chunk_beg, chunk_end))
                else:
                    # read but don't load the metadata pseudo-bin
                    assert n_chunk == 2, ERR_METADATA_CHUNKS % n_chunk
                    # 4x8 bytes for ref_[beg|end] and n_[mapped|unmapped]
                    fh.read(32)
            n_intv = struct.unpack('i', fh.read(4))[0]
            for window in range(n_intv):
                ioffset = struct.unpack('l', fh.read(8))[0]
                intvs.append((bam_id, ref_id, window, ioffset))
        # 8 bytes for n_no_coor
        fh.read(8)
        assert fh.read() == b'', ERR_NOT_EOF
    return(chunks, intvs)


# Design note: do as much input validation up-front as reasonably possible.
# For the CLI user this means earlier and friendlier errors. For the author
# this means freedom to write downstream functions with less explicit error
# checking and handling.
def parse_cmdargs(arglist):
    """
    Returns: valid parsed arguments as Namespace.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    subs = parser.add_subparsers(
        required=True,
        help=f'{sys.argv[0]} SUBCOMMAND --help for help on each subcommand')
    pcreate = subs.add_parser(
        'create',
        description='Create a new bam index database. If the DSN argument is '
        'supplied then a postgres server hosted database is inferred, '
        'otherwise a sqlite3 file database. ')
    pcreate.set_defaults(cmd=create_baidb)
    pcreate.add_argument(
        'newdb',
        help='Name of the database to create. If no DSN is supplied the value '
        'is interpreted as the path to the new sqlite3 database file.')
    pcreate.add_argument(
        'reference',
        choices=REFERENCE,
        help='Reference genome for all indexed bams.')
    pcreate.add_argument(
        'dsn',
        metavar='DSN',
        nargs='?',
        type=valid_dsn,
        help='If this argument is supplied it must be a valid libpq '
        'connection string for an existing postgres maintenance database and '
        f'user with CREATE DATABASE permissions. Examples: {MSG_DSN_FORMATS}. '
        'A corresponding credentials entry in the user\'s ~/.pgpass file is '
        'required.')
    pload = subs.add_parser(
        'load',
        description='Parse bai index file(s) into a bam index database.')
    pload.set_defaults(cmd=load_baidb)
    pload.add_argument(
        'dbid',
        metavar='DSN|PATH',
        type=valid_dbid,
        help='The database to load. If a valid libpq connection string is '
        f'supplied e.g. {MSG_DSN_FORMATS} then a postgres database is '
        'inferred, otherwise the value is interpreted as the path to an '
        'existing sqlite3 database file.')
    pload.add_argument(
        'bams',
        metavar='BAM',
        nargs='+',
        type=valid_bam_with_index,
        help='path to bam with adjacent .bai file.')
    parser.add_argument(
        '--version',
        action='version',
        version=metadata.version('mumbai'))

    parsed = parser.parse_args(arglist)

    if parsed.cmd is create_baidb and not parsed.dsn:
        parsed.newdb = abspath(expanduser(parsed.newdb))
        if exists(parsed.newdb):
            parser.error(ERR_DB_EXISTS % parsed.newdb)

    return parsed


def valid_bam_with_index(bam):
    """
    Validate bam is readable and has adjacent .bai file
    """
    for file_ in (bam, f'{bam}.bai'):
        try:
            with open(file_, 'rb'):
                pass
        except (FileNotFoundError, PermissionError) as ex:
            LOGGER.debug(ex)
            raise argparse.ArgumentTypeError(ERR_INVALID_INPUTFILE % file_)
    return bam


def valid_dbid(dbid):
    """
    Validate dbid as either a libpq connection string or path to a sqlite3 db

    Returns:
        If dbid is a libpq connection string, the DSN is returned as a
        dictionary of connection parameters.
        If dbid is the path to a sqlite3 db, the absolute path is returned.
    Raises:
        argparse.ArgumentTypeError if dbid is neither a valid connection
        string nor the path to a sqlite3 db.
    """
    try:
        return parse_dsn(dbid)
    except psycopg2.ProgrammingError as pgpe:
        LOGGER.debug(pgpe)
    path = abspath(expanduser(dbid))
    try:
        with open(path, 'rb') as fh:
            magic = fh.read(16).decode('utf-8')
            assert magic == MAGICSQLITE, ERR_BAD_MAGIC % path
        return path
    except Exception as ex:
        LOGGER.debug(ex)
    raise argparse.ArgumentTypeError(ERR_INVALID_DBID % dbid)


def valid_dsn(dsn):
    """
    Validate dsn as a libpq connection string

    Returns:
        dsn parsed into a dictionary of connection parameters.
    Raises:
        argparse.ArgumentTypeError if dsn is not a valid connection string.
    """
    try:
        return parse_dsn(dsn)
    except psycopg2.ProgrammingError as pgpe:
        LOGGER.debug(pgpe)
        raise argparse.ArgumentTypeError('Required format: ' + MSG_DSN_FORMATS)


def main():
    args = parse_cmdargs(sys.argv[1:])
    argsd = vars(args)
    cmd = argsd.pop('cmd')
    try:
        cmd(**argsd)
    except Exception as ex:
        LOGGER.debug(traceback.format_exc())
        LOGGER.error(ERR_UNEXPECTED, ex)


if __name__ == '__main__':
    main()
