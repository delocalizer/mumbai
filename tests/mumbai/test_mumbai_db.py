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
    (1, 0, 4744, 111543828, 111544217),
    (1, 0, 4681, 111542272, 111542571),
    # reads in bin 9 span the first 2**20 == 1048576 bp boundary
    (1, 0, 9, 111544217, 186777600),
    (1, 0, 4687, 111542571, 111543828)
]

TESTBAI_INTVS = [
    (1, 0, 0, 111542571), (1, 0, 1, 111542571), (1, 0, 2, 111542571),
    (1, 0, 3, 111542571), (1, 0, 4, 111542571), (1, 0, 5, 111542571),
    (1, 0, 6, 111542571), (1, 0, 7, 111543828), (1, 0, 8, 111543828),
    (1, 0, 9, 111543828), (1, 0, 10, 111543828), (1, 0, 11, 111543828),
    (1, 0, 12, 111543828), (1, 0, 13, 111543828), (1, 0, 14, 111543828),
    (1, 0, 15, 111543828), (1, 0, 16, 111543828), (1, 0, 17, 111543828),
    (1, 0, 18, 111543828), (1, 0, 19, 111543828), (1, 0, 20, 111543828),
    (1, 0, 21, 111543828), (1, 0, 22, 111543828), (1, 0, 23, 111543828),
    (1, 0, 24, 111543828), (1, 0, 25, 111543828), (1, 0, 26, 111543828),
    (1, 0, 27, 111543828), (1, 0, 28, 111543828), (1, 0, 29, 111543828),
    (1, 0, 30, 111543828), (1, 0, 31, 111543828), (1, 0, 32, 111543828),
    (1, 0, 33, 111543828), (1, 0, 34, 111543828), (1, 0, 35, 111543828),
    (1, 0, 36, 111543828), (1, 0, 37, 111543828), (1, 0, 38, 111543828),
    (1, 0, 39, 111543828), (1, 0, 40, 111543828), (1, 0, 41, 111543828),
    (1, 0, 42, 111543828), (1, 0, 43, 111543828), (1, 0, 44, 111543828),
    (1, 0, 45, 111543828), (1, 0, 46, 111543828), (1, 0, 47, 111543828),
    (1, 0, 48, 111543828), (1, 0, 49, 111543828), (1, 0, 50, 111543828),
    (1, 0, 51, 111543828), (1, 0, 52, 111543828), (1, 0, 53, 111543828),
    (1, 0, 54, 111543828), (1, 0, 55, 111543828), (1, 0, 56, 111543828),
    (1, 0, 57, 111543828), (1, 0, 58, 111543828), (1, 0, 59, 111543828),
    (1, 0, 60, 111543828), (1, 0, 61, 111543828), (1, 0, 62, 111543828),
    (1, 0, 63, 111543828), (1, 0, 64, 111544217)
]

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
