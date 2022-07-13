"""
Test mumbai module functions.

Note that the test resource 'test.db' was created and loaded using mumbai_db
so many of the tests here implicitly test that module's logic for creating and
loading a bam index database.
"""
import contextlib
import importlib.resources as pkg_resources
import io
from pathlib import Path
import sys
from unittest import TestCase
from unittest.mock import patch

from pysam import AlignmentFile

from . import resources

from mumbai import mumbai


TESTDB = pkg_resources.path(resources, 'test.db')
TESTBAM = pkg_resources.path(resources,
                             'a05168de-0e6d-4867-82b1-22d330dac0f8.bam')
# The virtual file offset where SAM records start in TESTBAM
TESTBAM_RECORDSTART = 111738880
TESTBAM_HEADER = AlignmentFile(TESTBAM).header

UNPATCHED_bams_offsets = mumbai.bams_offsets


def patched_paths(dbid, reg):
    """
    Patch bam paths returned from the test database so they refer to test bams
    as they existxs locally.
    """
    return [
        (TESTBAM.parent / Path(bam).name, offset)
        for bam, offset in UNPATCHED_bams_offsets(dbid, reg)
    ]


class TestMumbai(TestCase):
    """
    Test mumbai module functions.
    """

    def test_reg2bins_1(self):
        """
        Check that 1-based inclusive region [1,1] gives expected bins.
        """
        expected = (0, 1, 9, 73, 585, 4681)
        self.assertEqual(tuple(mumbai.reg2bins(0, 0)), expected)

    def test_reg2bins_2(self):
        """
        Check that 1-based inclusive region [16383,16383] gives expected bins.
        """
        expected = (0, 1, 9, 73, 585, 4681)
        self.assertEqual(tuple(mumbai.reg2bins(16383, 16383)), expected)

    def test_reg2bins_3(self):
        """
        Check that 1-based inclusive region [16383,16384] gives expected bins.
        """
        expected = (0, 1, 9, 73, 585, 4681, 4682)
        self.assertEqual(tuple(mumbai.reg2bins(16383, 16384)), expected)

    def test_bam_fetch_1(self):
        """
        Check that expected reads are returned from the bam.
        Inspect a05168de-0e6d-4867-82b1-22d330dac0f8.sam for details.
        The reads are paired 101bp.
        """
        reg = mumbai.Region('chr1', 10379, 10379)
        qnames = [sam.split('\t')[0] for pos, sam in mumbai.bam_fetch(
                  TESTBAM, ['chr1'], [249250621], TESTBAM_RECORDSTART,
                  reg, 'sam')]
        # 0 mapped qnames overlap chr1:10379-10379; the one placed at
        # POS = 10379 is unmapped so not returned.
        expected = []
        self.assertEqual(qnames, expected)

    def test_bam_fetch_2(self):
        """
        Check that expected reads are returned from the bam.
        Inspect a05168de-0e6d-4867-82b1-22d330dac0f8.sam for details.
        The reads are paired 101bp.
        """
        reg = mumbai.Region('chr1', 100000, 100000)
        qnames = [sam.split('\t')[0] for pos, sam in mumbai.bam_fetch(
                  TESTBAM, ['chr1'], [249250621], TESTBAM_RECORDSTART,
                  reg, 'sam')]
        # 1 mapped record overlaps chr1:100000-100000; the one that starts at
        # POS = 100000
        expected = ['D8QSB6V1:72:C1K64ACXX:4:2311:13531:97563']
        self.assertEqual(qnames, expected)

    def test_bam_fetch_3(self):
        """
        Check that expected reads are returned from the bam.
        Inspect a05168de-0e6d-4867-82b1-22d330dac0f8.sam for details.
        The reads are paired 101bp.
        """
        reg = mumbai.Region('chr1', 100010, 100010)
        qnames = [sam.split('\t')[0] for pos, sam in mumbai.bam_fetch(
                  TESTBAM, ['chr1'], [249250621], TESTBAM_RECORDSTART,
                  reg, 'sam')]
        # 3 mapped qnames overlap chr1:100010-100010; the one that starts at
        # POS = 100000 and both that start at 100005
        expected = ['D8QSB6V1:72:C1K64ACXX:4:2311:13531:97563',
                    'HWI-ST1213:151:C1DTBACXX:2:1210:12782:5555',
                    'artefact_soft_clip_1']
        self.assertEqual(qnames, expected)

    def test_bam_fetch_4(self):
        """
        Check that expected reads are returned from the bam.
        Inspect a05168de-0e6d-4867-82b1-22d330dac0f8.sam for details.
        The reads are paired 101bp.
        """
        reg = mumbai.Region('chr1', 100100, 100100)
        qnames = [sam.split('\t')[0] for pos, sam in mumbai.bam_fetch(
                  TESTBAM, ['chr1'], [249250621], TESTBAM_RECORDSTART,
                  reg, 'sam')]
        # 2 mapped qnames overlap chr1:100100-100100; the one that starts at
        # POS = 100000 and the unclipped one that starts at POS = 100005
        expected = ['D8QSB6V1:72:C1K64ACXX:4:2311:13531:97563',
                    'HWI-ST1213:151:C1DTBACXX:2:1210:12782:5555']
        self.assertEqual(qnames, expected)

    def test_bam_fetch_5(self):
        """
        Check that expected reads are returned from the bam.
        Inspect a05168de-0e6d-4867-82b1-22d330dac0f8.sam for details.
        The reads are paired 101bp.
        """
        reg = mumbai.Region('chr1', 100101, 100101)
        qnames = [sam.split('\t')[0] for pos, sam in mumbai.bam_fetch(
                  TESTBAM, ['chr1'], [249250621], TESTBAM_RECORDSTART,
                  reg, 'sam')]
        # 2 mapped qnames overlap chr1:100101-100101; the unclipped one that
        # starts at POS = 100005 and the one that starts at POS = 100101
        expected = ['HWI-ST1213:151:C1DTBACXX:2:1210:12782:5555',
                    'D8QSB6V1:72:C1K64ACXX:4:2206:20773:63822']
        self.assertEqual(qnames, expected)

    def test_bams_fetch_1(self):
        """
        Check that querying the db for a region with no reads returns empty.
        """
        reg = mumbai.Region('chr1', 20000, 20000)
        records = mumbai.bams_fetch(TESTDB, reg, 'sam')
        expected = None
        self.assertEqual(records, expected)

    def test_bams_fetch_2(self):
        """
        Check that querying the db for a region with mapped reads returns them.
        """
        reg = mumbai.Region('chr1', 100000, 100000)
        with patch('mumbai.mumbai.bams_offsets', side_effect=patched_paths):
            header, records = mumbai.bams_fetch(TESTDB, reg, 'sam')
            expected = ['D8QSB6V1:72:C1K64ACXX:4:2311:13531:97563\t163\tchr1\t100000\t0\t101M\t=\t100191\t292\tCACTAAGCACACAGAGAATAATGTCTAGAATCTGAGTGCCATGTTATCAAATTGTACTGAGACTCTTGCAGTCACACAGGCTGACATGTAAGCATCGCCAT\tCCCFFFFFHHHHHJJJJJJJJIJIJJIJJJJJJJGHCHIIJIJHHJJJJJFHJJHIJJJIJJJJJJIJIIJIJJIJJIJHHGHFFFFFCEECEEEDDDDDD\tZC:i:2\tMD:Z:101\tPG:Z:MarkDuplicates\tRG:Z:9f46438b-c352-4153-b46e-5c057ea58c90\tNM:i:0\tAS:i:101\tXS:i:101']
            self.assertEqual(records, expected)

    def test_count_mode(self):
        """
        Integration test; implicitly test parse_cmdargs, bams_offsets and
        bams_fetch.
        """
        testargs = list(map(str,
                        ('mumbai', 'count', TESTDB, 'chr1', 1048575, 1048576)))
        expected = '3 records\n'
        with patch('mumbai.mumbai.bams_offsets', side_effect=patched_paths):
            with patch.object(sys, 'argv', testargs):
                outstr = io.StringIO()
                with contextlib.redirect_stdout(outstr):
                    mumbai.main()
                    self.assertEqual(outstr.getvalue(), expected)

    def test_pileup_mode_1(self):
        """
        Integration test.
        """
        testargs = list(map(str,
                        ('mumbai', 'pileup', TESTDB, 'chr1', 100050, 100050)))
        expected = """chr1\t100050\tA\t      A:3      C:0      G:0      T:0      -:0      |:0\n"""
        with patch('mumbai.mumbai.bams_offsets', side_effect=patched_paths):
            with patch.object(sys, 'argv', testargs):
                outstr = io.StringIO()
                with contextlib.redirect_stdout(outstr):
                    mumbai.main()
                    self.assertEqual(outstr.getvalue(), expected)


    def test_pileup_mode_2(self):
        """
        Integration test.
        """
        testargs = list(map(str,
                        ('mumbai', 'pileup', TESTDB, 'chr1', 100101, 100101)))
        expected = """chr1\t100101\tG\t      G:2      A:0      C:0      T:0      -:0      |:0\n"""
        with patch('mumbai.mumbai.bams_offsets', side_effect=patched_paths):
            with patch.object(sys, 'argv', testargs):
                outstr = io.StringIO()
                with contextlib.redirect_stdout(outstr):
                    mumbai.main()
                    self.assertEqual(outstr.getvalue(), expected)


    def test_sam_mode(self):
        """
        Integration test.
        """
        testargs = list(map(str,
                        ('mumbai', 'sam', TESTDB, 'chr1', 1048575, 1048576)))
        expected = pkg_resources.read_text(
            resources, 'chr1:1048575-1048576.sam').replace(
                    'tests/mumbai/resources/test.db',
                    str(TESTDB))
        with patch('mumbai.mumbai.bams_offsets', side_effect=patched_paths):
            with patch.object(sys, 'argv', testargs):
                outstr = io.StringIO()
                with contextlib.redirect_stdout(outstr):
                    mumbai.main()
                    self.assertEqual(outstr.getvalue(), expected)

    def test_tview_mode(self):
        """
        Integration test.
        """
        testargs = list(map(str,
                        ('mumbai', 'tview', TESTDB, 'chr1', 1048575, 1048576)))
        expected = """CTGAGGCAGGAGAATGACATGAACCCGGGAGGCAGAGCTTGCAGTGAGCCGAGATTGTGCCACTGCACTCCAGCCTGGGCGACAGAGACAGACTCCATCT\x1b[1mC\x1b[0m
                                                  GAGATTGTGCCACTGCACTCCAGCCTGGGCGACAGAGACAGACTCCATCT\x1b[1mCA\x1b[0mAAAAAAAAAAAAAAA
                                                                                                     \x1b[1ma\x1b[0maaaaaaaaaaaaaaacttgcctcagtctctccttcagcctcatgtccttccatcaaattctttcttctgaggcagcgagaatcgaggccgctgctgacat\n"""
        with patch('mumbai.mumbai.bams_offsets', side_effect=patched_paths):
            with patch.object(sys, 'argv', testargs):
                outstr = io.StringIO()
                with contextlib.redirect_stdout(outstr):
                    mumbai.main()
                    self.assertEqual(outstr.getvalue(), expected)

    def test_count_pileup_consistency(self):
        """
        Sanity check that count at a position == pileup total.
        """
        regions = (
            mumbai.Region('chr1', 100000, 100000),
            mumbai.Region('chr1', 100050, 100050),
            mumbai.Region('chr1', 100101, 100101),
            mumbai.Region('chr1', 100151, 100151),
        )
        for reg in regions:
            with patch('mumbai.mumbai.bams_offsets',
                    side_effect=patched_paths):
                count = mumbai.bams_fetch(TESTDB, reg, 'count')
                pileup = mumbai.bams_fetch(TESTDB, reg, 'pileup')
                putotal = sum(pileup[(reg.contig, reg.start)][base]
                              for base in ('A', 'C', 'G', 'T', '-' , '|'))
                self.assertEqual(count, putotal)
