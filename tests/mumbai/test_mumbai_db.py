"""
Test mumbai_db module functions.
"""
import importlib.resources as pkg_resources
import sqlite3
import sys

from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase
from unittest.mock import patch

from . import resources

from mumbai import mumbai_db

TESTBAI = pkg_resources.path(resources,
                             'a05168de-0e6d-4867-82b1-22d330dac0f8.bam.bai')
TESTBAM = pkg_resources.path(resources,
                             'a05168de-0e6d-4867-82b1-22d330dac0f8.bam')

TESTBAI_CHUNKS = [
    (1, 0, 4744, 111740733, 111741122),
    (1, 0, 4681, 111738880, 111739179),
    # reads in bin 9 span the first 2**20 == 1048576 bp boundary
    (1, 0, 9, 111741122, 191102976),
    (1, 0, 4687, 111739179, 111740733)
]

TESTBAI_INTVS = [
    (1, 0, 0, 111738880), (1, 0, 1, 111738880), (1, 0, 2, 111738880),
    (1, 0, 3, 111738880), (1, 0, 4, 111738880), (1, 0, 5, 111738880),
    (1, 0, 6, 111739179), (1, 0, 7, 111739179), (1, 0, 8, 111739179),
    (1, 0, 9, 111739179), (1, 0, 10, 111739179), (1, 0, 11, 111739179),
    (1, 0, 12, 111739179), (1, 0, 13, 111739179), (1, 0, 14, 111739179),
    (1, 0, 15, 111739179), (1, 0, 16, 111739179), (1, 0, 17, 111739179),
    (1, 0, 18, 111739179), (1, 0, 19, 111739179), (1, 0, 20, 111739179),
    (1, 0, 21, 111739179), (1, 0, 22, 111739179), (1, 0, 23, 111739179),
    (1, 0, 24, 111739179), (1, 0, 25, 111739179), (1, 0, 26, 111739179),
    (1, 0, 27, 111739179), (1, 0, 28, 111739179), (1, 0, 29, 111739179),
    (1, 0, 30, 111739179), (1, 0, 31, 111739179), (1, 0, 32, 111739179),
    (1, 0, 33, 111739179), (1, 0, 34, 111739179), (1, 0, 35, 111739179),
    (1, 0, 36, 111739179), (1, 0, 37, 111739179), (1, 0, 38, 111739179),
    (1, 0, 39, 111739179), (1, 0, 40, 111739179), (1, 0, 41, 111739179),
    (1, 0, 42, 111739179), (1, 0, 43, 111739179), (1, 0, 44, 111739179),
    (1, 0, 45, 111739179), (1, 0, 46, 111739179), (1, 0, 47, 111739179),
    (1, 0, 48, 111739179), (1, 0, 49, 111739179), (1, 0, 50, 111739179),
    (1, 0, 51, 111739179), (1, 0, 52, 111739179), (1, 0, 53, 111739179),
    (1, 0, 54, 111739179), (1, 0, 55, 111739179), (1, 0, 56, 111739179),
    (1, 0, 57, 111739179), (1, 0, 58, 111739179), (1, 0, 59, 111739179),
    (1, 0, 60, 111739179), (1, 0, 61, 111739179), (1, 0, 62, 111739179),
    (1, 0, 63, 111740733), (1, 0, 64, 111741122), ]

class TestMumbaiDb(TestCase):
    """
    Test mumbai_db module functions
    """

    def test_parse_bai(self):
        """
        Check that bai data is read as expected.
        """
        chunks, intvs = mumbai_db.parse_bai(1, TESTBAI)
        self.assertEqual(chunks, TESTBAI_CHUNKS)
        self.assertEqual(intvs, TESTBAI_INTVS)

    def test_create_baidb(self):
        """
        Check that sqlite3 database is created as expected.
        (No test for postgres as that would require either local postgresql
        server or heroic mocking to no real gain)
        """
        with TemporaryDirectory() as tmpdir:
            createdb = Path(tmpdir) / 'test.db'
            testargs = list(map(str,
                            ('mumbai_db', 'create', createdb, 'GRCh37')))
            with patch.object(sys, 'argv', testargs):
                mumbai_db.main()
                con = sqlite3.connect(createdb)
                cur = con.cursor()
                cur.execute('SELECT sql FROM sqlite_master WHERE sql NOT NULL')
                schema = '\n'.join(f'{row[0]};' for row in cur.fetchall())+'\n'
                expected = pkg_resources.read_text(resources,
                                                   'sqlite3_schema.sql')
                self.assertEqual(schema, expected)
                cur.close()

    def test_load_baidb(self):
        """
        Check that sqlite3 database is loaded as expected.
        (No test for postgres as that would require either local postgresql
        server or heroic mocking to no real gain)
        """
        with TemporaryDirectory() as tmpdir:
            createdb = Path(tmpdir) / 'test.db'
            testargs = list(map(str,
                            ('mumbai_db', 'create', createdb, 'GRCh37')))
            with patch.object(sys, 'argv', testargs):
                mumbai_db.main()
            testargs = list(map(str,
                            ('mumbai_db', 'load', createdb, TESTBAM)))
            with patch.object(sys, 'argv', testargs):
                mumbai_db.main()
            con = sqlite3.connect(createdb)
            cur = con.cursor()
            cur.execute('SELECT * FROM bam')
            bams = cur.fetchall()
            cur.execute('SELECT * FROM chunk')
            chunks = cur.fetchall()
            cur.execute('SELECT * FROM intv')
            intvs = cur.fetchall()
            self.assertEqual(bams, [(1, f'{TESTBAM.stem}', f'{TESTBAM}')])
            self.assertEqual(chunks, TESTBAI_CHUNKS)
            self.assertEqual(intvs, TESTBAI_INTVS)
            cur.close()
