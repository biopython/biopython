# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO BlatIO parsers."""

import os
import unittest

from Bio.SearchIO import parse


# test case files are in the Blast directory
TEST_DIR = "Blat"
FMT = "blat-psl"


def get_file(filename):
    """Return the path of a test file."""
    return os.path.join(TEST_DIR, filename)


class BlatPslCases(unittest.TestCase):
    def test_psl_34_001(self, testf="psl_34_001.psl", pslx=False):
        """Test parsing blat output (psl_34_001.psl)."""
        blat_file = get_file(testf)
        self.qresults = list(parse(blat_file, FMT, pslx=pslx))
        self.assertEqual(2, len(self.qresults))
        # check common attributes
        for qresult in self.qresults:
            for hit in qresult:
                self.assertEqual(qresult.id, hit.query_id)
                for hsp in hit:
                    self.assertEqual(hit.id, hsp.hit_id)
                    self.assertEqual(qresult.id, hsp.query_id)

        # test first qresult
        qresult = self.qresults[0]
        self.assertEqual("hg18_dna", qresult.id)
        self.assertEqual("blat", qresult.program)
        self.assertEqual(33, qresult.seq_len)
        self.assertEqual(3, len(qresult))
        # first qresult, first hit
        hit = qresult[0]
        self.assertEqual("chr4", hit.id)
        self.assertEqual(191154276, hit.seq_len)
        self.assertEqual(1, len(hit.hsps))
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(16, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(0, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(11, hsp.query_start)
        self.assertEqual(61646095, hsp.hit_start)
        self.assertEqual(27, hsp.query_end)
        self.assertEqual(61646111, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([16], hsp.query_span_all)
        self.assertEqual([16], hsp.hit_span_all)
        self.assertEqual([(11, 27)], hsp.query_range_all)
        self.assertEqual([(61646095, 61646111)], hsp.hit_range_all)
        # first qresult, second hit
        hit = qresult[1]
        self.assertEqual("chr1", hit.id)
        self.assertEqual(249250621, hit.seq_len)
        self.assertEqual(1, len(hit.hsps))
        # first qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(33, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(0, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(10271783, hsp.hit_start)
        self.assertEqual(33, hsp.query_end)
        self.assertEqual(10271816, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([33], hsp.query_span_all)
        self.assertEqual([33], hsp.hit_span_all)
        self.assertEqual([(0, 33)], hsp.query_range_all)
        self.assertEqual([(10271783, 10271816)], hsp.hit_range_all)
        # first qresult, third hit
        hit = qresult[2]
        self.assertEqual("chr2", hit.id)
        self.assertEqual(243199373, hit.seq_len)
        self.assertEqual(1, len(hit.hsps))
        # first qresult, third hit, first hsp
        hsp = qresult[2].hsps[0]
        self.assertEqual(17, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(0, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(-1, hsp[0].query_strand)
        self.assertEqual(8, hsp.query_start)
        self.assertEqual(53575980, hsp.hit_start)
        self.assertEqual(25, hsp.query_end)
        self.assertEqual(53575997, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([17], hsp.query_span_all)
        self.assertEqual([17], hsp.hit_span_all)
        self.assertEqual([(8, 25)], hsp.query_range_all)
        self.assertEqual([(53575980, 53575997)], hsp.hit_range_all)

        # test second qresult
        qresult = self.qresults[1]
        self.assertEqual("hg19_dna", qresult.id)
        self.assertEqual("blat", qresult.program)
        self.assertEqual(50, qresult.seq_len)
        self.assertEqual(10, len(qresult))
        # second qresult, first hit
        hit = qresult[0]
        self.assertEqual("chr9", hit.id)
        self.assertEqual(141213431, hit.seq_len)
        self.assertEqual(1, len(hit.hsps))
        # second qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(38, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(3, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(9, hsp.query_start)
        self.assertEqual(85737865, hsp.hit_start)
        self.assertEqual(50, hsp.query_end)
        self.assertEqual(85737906, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([41], hsp.query_span_all)
        self.assertEqual([41], hsp.hit_span_all)
        self.assertEqual([(9, 50)], hsp.query_range_all)
        self.assertEqual([(85737865, 85737906)], hsp.hit_range_all)
        # second qresult, second hit
        hit = qresult[1]
        self.assertEqual("chr8", hit.id)
        self.assertEqual(146364022, hit.seq_len)
        self.assertEqual(1, len(hit.hsps))
        # second qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(41, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(0, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(8, hsp.query_start)
        self.assertEqual(95160479, hsp.hit_start)
        self.assertEqual(49, hsp.query_end)
        self.assertEqual(95160520, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([41], hsp.query_span_all)
        self.assertEqual([41], hsp.hit_span_all)
        self.assertEqual([(8, 49)], hsp.query_range_all)
        self.assertEqual([(95160479, 95160520)], hsp.hit_range_all)
        # second qresult, third hit
        hit = qresult[2]
        self.assertEqual("chr22", hit.id)
        self.assertEqual(51304566, hit.seq_len)
        self.assertEqual(2, len(hit.hsps))
        # second qresult, third hit, first hsp
        hsp = qresult[2].hsps[0]
        self.assertEqual(33, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(3, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(11, hsp.query_start)
        self.assertEqual(42144400, hsp.hit_start)
        self.assertEqual(47, hsp.query_end)
        self.assertEqual(42144436, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([36], hsp.query_span_all)
        self.assertEqual([36], hsp.hit_span_all)
        self.assertEqual([(11, 47)], hsp.query_range_all)
        self.assertEqual([(42144400, 42144436)], hsp.hit_range_all)
        # second qresult, third hit, second hsp
        hsp = qresult[2].hsps[1]
        self.assertEqual(35, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(2, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(-1, hsp[0].query_strand)
        self.assertEqual(12, hsp.query_start)
        self.assertEqual(48997405, hsp.hit_start)
        self.assertEqual(49, hsp.query_end)
        self.assertEqual(48997442, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([37], hsp.query_span_all)
        self.assertEqual([37], hsp.hit_span_all)
        self.assertEqual([(12, 49)], hsp.query_range_all)
        self.assertEqual([(48997405, 48997442)], hsp.hit_range_all)
        # second qresult, fourth hit
        hit = qresult[3]
        self.assertEqual("chr2", hit.id)
        self.assertEqual(243199373, hit.seq_len)
        self.assertEqual(2, len(hit.hsps))
        # second qresult, fourth hit, first hsp
        hsp = qresult[3].hsps[0]
        self.assertEqual(43, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(1, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(1, hsp.query_gapopen_num)
        self.assertEqual(4, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(1, hsp.query_start)
        self.assertEqual(183925984, hsp.hit_start)
        self.assertEqual(49, hsp.query_end)
        self.assertEqual(183926028, hsp.hit_end)
        self.assertEqual(2, len(hsp))
        self.assertEqual([6, 38], hsp.query_span_all)
        self.assertEqual([6, 38], hsp.hit_span_all)
        self.assertEqual([(1, 7), (11, 49)], hsp.query_range_all)
        self.assertEqual([1, 11], hsp.query_start_all)
        self.assertEqual([7, 49], hsp.query_end_all)
        self.assertEqual(
            [(183925984, 183925990), (183925990, 183926028)], hsp.hit_range_all
        )
        self.assertEqual([183925984, 183925990], hsp.hit_start_all)
        self.assertEqual([183925990, 183926028], hsp.hit_end_all)
        # second qresult, fourth hit, second hsp
        hsp = qresult[3].hsps[1]
        self.assertEqual(35, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(1, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(-1, hsp[0].query_strand)
        self.assertEqual(13, hsp.query_start)
        self.assertEqual(120641740, hsp.hit_start)
        self.assertEqual(49, hsp.query_end)
        self.assertEqual(120641776, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([36], hsp.query_span_all)
        self.assertEqual([36], hsp.hit_span_all)
        self.assertEqual([(13, 49)], hsp.query_range_all)
        self.assertEqual([(120641740, 120641776)], hsp.hit_range_all)
        # second qresult, fifth hit
        hit = qresult[4]
        self.assertEqual("chr19", hit.id)
        self.assertEqual(59128983, hit.seq_len)
        self.assertEqual(3, len(hit.hsps))
        # second qresult, fifth hit, first hsp
        hsp = qresult[4].hsps[0]
        self.assertEqual(34, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(2, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(1, hsp.hit_gapopen_num)
        self.assertEqual(134, hsp.hit_gap_num)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(10, hsp.query_start)
        self.assertEqual(35483340, hsp.hit_start)
        self.assertEqual(46, hsp.query_end)
        self.assertEqual(35483510, hsp.hit_end)
        self.assertEqual(2, len(hsp))
        self.assertEqual([25, 11], hsp.query_span_all)
        self.assertEqual([25, 11], hsp.hit_span_all)
        self.assertEqual([(10, 35), (35, 46)], hsp.query_range_all)
        self.assertEqual(
            [(35483340, 35483365), (35483499, 35483510)], hsp.hit_range_all
        )
        # second qresult, fifth hit, second hsp
        hsp = qresult[4].hsps[1]
        self.assertEqual(39, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(0, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(-1, hsp[0].query_strand)
        self.assertEqual(10, hsp.query_start)
        self.assertEqual(54017130, hsp.hit_start)
        self.assertEqual(49, hsp.query_end)
        self.assertEqual(54017169, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([39], hsp.query_span_all)
        self.assertEqual([39], hsp.hit_span_all)
        self.assertEqual([(10, 49)], hsp.query_range_all)
        self.assertEqual([(54017130, 54017169)], hsp.hit_range_all)
        # second qresult, fifth hit, third hsp
        hsp = qresult[4].hsps[2]
        self.assertEqual(36, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(3, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(-1, hsp[0].query_strand)
        self.assertEqual(10, hsp.query_start)
        self.assertEqual(553742, hsp.hit_start)
        self.assertEqual(49, hsp.query_end)
        self.assertEqual(553781, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([39], hsp.query_span_all)
        self.assertEqual([39], hsp.hit_span_all)
        self.assertEqual([(10, 49)], hsp.query_range_all)
        self.assertEqual([(553742, 553781)], hsp.hit_range_all)

    def test_psl_34_002(self, testf="psl_34_002.psl", pslx=False):
        """Test parsing blat output (psl_34_001.psl)."""
        blat_file = get_file(testf)
        self.qresults = list(parse(blat_file, FMT, pslx=pslx))
        self.assertEqual(0, len(self.qresults))

    def test_psl_34_003(self, testf="psl_34_003.psl", pslx=False):
        """Test parsing blat output (psl_34_003.psl)."""
        blat_file = get_file(testf)
        self.qresults = list(parse(blat_file, FMT, pslx=pslx))
        self.assertEqual(1, len(self.qresults))
        # check common attributes
        for qresult in self.qresults:
            for hit in qresult:
                self.assertEqual(qresult.id, hit.query_id)
                for hsp in hit:
                    self.assertEqual(hit.id, hsp.hit_id)
                    self.assertEqual(qresult.id, hsp.query_id)

        # test first qresult
        qresult = self.qresults[0]
        self.assertEqual("hg18_dna", qresult.id)
        self.assertEqual("blat", qresult.program)
        self.assertEqual(33, qresult.seq_len)
        self.assertEqual(3, len(qresult))
        # first qresult, first hit
        hit = qresult[0]
        self.assertEqual("chr4", hit.id)
        self.assertEqual(191154276, hit.seq_len)
        self.assertEqual(1, len(hit.hsps))
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(16, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(0, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(11, hsp.query_start)
        self.assertEqual(61646095, hsp.hit_start)
        self.assertEqual(27, hsp.query_end)
        self.assertEqual(61646111, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([16], hsp.query_span_all)
        self.assertEqual([16], hsp.hit_span_all)
        self.assertEqual([(11, 27)], hsp.query_range_all)
        self.assertEqual([(61646095, 61646111)], hsp.hit_range_all)
        # first qresult, second hit
        hit = qresult[1]
        self.assertEqual("chr1", hit.id)
        self.assertEqual(249250621, hit.seq_len)
        self.assertEqual(1, len(hit.hsps))
        # first qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(33, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(0, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(10271783, hsp.hit_start)
        self.assertEqual(33, hsp.query_end)
        self.assertEqual(10271816, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([33], hsp.query_span_all)
        self.assertEqual([33], hsp.hit_span_all)
        self.assertEqual([(0, 33)], hsp.query_range_all)
        self.assertEqual([(10271783, 10271816)], hsp.hit_range_all)
        # first qresult, third hit
        hit = qresult[2]
        self.assertEqual("chr2", hit.id)
        self.assertEqual(243199373, hit.seq_len)
        self.assertEqual(1, len(hit.hsps))
        # first qresult, third hit, first hsp
        hsp = qresult[2].hsps[0]
        self.assertEqual(17, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(0, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(-1, hsp[0].query_strand)
        self.assertEqual(8, hsp.query_start)
        self.assertEqual(53575980, hsp.hit_start)
        self.assertEqual(25, hsp.query_end)
        self.assertEqual(53575997, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([17], hsp.query_span_all)
        self.assertEqual([17], hsp.hit_span_all)
        self.assertEqual([(8, 25)], hsp.query_range_all)
        self.assertEqual([(53575980, 53575997)], hsp.hit_range_all)

    def test_psl_34_004(self, testf="psl_34_004.psl", pslx=False):
        """Test parsing blat output (psl_34_004.psl)."""
        blat_file = get_file(testf)
        self.qresults = list(parse(blat_file, FMT, pslx=pslx))
        self.assertEqual(1, len(self.qresults))
        # check common attributes
        for qresult in self.qresults:
            for hit in qresult:
                self.assertEqual(qresult.id, hit.query_id)
                for hsp in hit:
                    self.assertEqual(hit.id, hsp.hit_id)
                    self.assertEqual(qresult.id, hsp.query_id)

        # test first qresult
        qresult = self.qresults[0]
        self.assertEqual("hg19_dna", qresult.id)
        self.assertEqual("blat", qresult.program)
        self.assertEqual(50, qresult.seq_len)
        self.assertEqual(10, len(qresult))
        # first qresult, first hit
        hit = qresult[0]
        self.assertEqual("chr9", hit.id)
        self.assertEqual(141213431, hit.seq_len)
        self.assertEqual(1, len(hit.hsps))
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(38, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(3, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(9, hsp.query_start)
        self.assertEqual(85737865, hsp.hit_start)
        self.assertEqual(50, hsp.query_end)
        self.assertEqual(85737906, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([41], hsp.query_span_all)
        self.assertEqual([41], hsp.hit_span_all)
        self.assertEqual([(9, 50)], hsp.query_range_all)
        self.assertEqual([(85737865, 85737906)], hsp.hit_range_all)
        # first qresult, second hit
        hit = qresult[1]
        self.assertEqual("chr8", hit.id)
        self.assertEqual(146364022, hit.seq_len)
        self.assertEqual(1, len(hit.hsps))
        # first qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(41, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(0, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(8, hsp.query_start)
        self.assertEqual(95160479, hsp.hit_start)
        self.assertEqual(49, hsp.query_end)
        self.assertEqual(95160520, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([41], hsp.query_span_all)
        self.assertEqual([41], hsp.hit_span_all)
        self.assertEqual([(8, 49)], hsp.query_range_all)
        self.assertEqual([(95160479, 95160520)], hsp.hit_range_all)
        # first qresult, third hit
        hit = qresult[2]
        self.assertEqual("chr22", hit.id)
        self.assertEqual(51304566, hit.seq_len)
        self.assertEqual(2, len(hit.hsps))
        # first qresult, third hit, first hsp
        hsp = qresult[2].hsps[0]
        self.assertEqual(33, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(3, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(11, hsp.query_start)
        self.assertEqual(42144400, hsp.hit_start)
        self.assertEqual(47, hsp.query_end)
        self.assertEqual(42144436, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([36], hsp.query_span_all)
        self.assertEqual([36], hsp.hit_span_all)
        self.assertEqual([(11, 47)], hsp.query_range_all)
        self.assertEqual([(42144400, 42144436)], hsp.hit_range_all)
        # first qresult, third hit, second hsp
        hsp = qresult[2].hsps[1]
        self.assertEqual(35, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(2, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(-1, hsp[0].query_strand)
        self.assertEqual(12, hsp.query_start)
        self.assertEqual(48997405, hsp.hit_start)
        self.assertEqual(49, hsp.query_end)
        self.assertEqual(48997442, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([37], hsp.query_span_all)
        self.assertEqual([37], hsp.hit_span_all)
        self.assertEqual([(12, 49)], hsp.query_range_all)
        self.assertEqual([(48997405, 48997442)], hsp.hit_range_all)
        # first qresult, fourth hit
        hit = qresult[3]
        self.assertEqual("chr2", hit.id)
        self.assertEqual(243199373, hit.seq_len)
        self.assertEqual(2, len(hit.hsps))
        # first qresult, fourth hit, first hsp
        hsp = qresult[3].hsps[0]
        self.assertEqual(43, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(1, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(1, hsp.query_gapopen_num)
        self.assertEqual(4, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(1, hsp.query_start)
        self.assertEqual(183925984, hsp.hit_start)
        self.assertEqual(49, hsp.query_end)
        self.assertEqual(183926028, hsp.hit_end)
        self.assertEqual(2, len(hsp))
        self.assertEqual([6, 38], hsp.query_span_all)
        self.assertEqual([6, 38], hsp.hit_span_all)
        self.assertEqual([(1, 7), (11, 49)], hsp.query_range_all)
        self.assertEqual(
            [(183925984, 183925990), (183925990, 183926028)], hsp.hit_range_all
        )
        # first qresult, fourth hit, second hsp
        hsp = qresult[3].hsps[1]
        self.assertEqual(35, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(1, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(-1, hsp[0].query_strand)
        self.assertEqual(13, hsp.query_start)
        self.assertEqual(120641740, hsp.hit_start)
        self.assertEqual(49, hsp.query_end)
        self.assertEqual(120641776, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([36], hsp.query_span_all)
        self.assertEqual([36], hsp.hit_span_all)
        self.assertEqual([(13, 49)], hsp.query_range_all)
        self.assertEqual([(120641740, 120641776)], hsp.hit_range_all)
        # first qresult, fifth hit
        hit = qresult[4]
        self.assertEqual("chr19", hit.id)
        self.assertEqual(59128983, hit.seq_len)
        self.assertEqual(3, len(hit.hsps))
        # first qresult, fifth hit, first hsp
        hsp = qresult[4].hsps[0]
        self.assertEqual(34, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(2, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(1, hsp.hit_gapopen_num)
        self.assertEqual(134, hsp.hit_gap_num)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(10, hsp.query_start)
        self.assertEqual(35483340, hsp.hit_start)
        self.assertEqual(46, hsp.query_end)
        self.assertEqual(35483510, hsp.hit_end)
        self.assertEqual(2, len(hsp))
        self.assertEqual([25, 11], hsp.query_span_all)
        self.assertEqual([25, 11], hsp.hit_span_all)
        self.assertEqual([(10, 35), (35, 46)], hsp.query_range_all)
        self.assertEqual(
            [(35483340, 35483365), (35483499, 35483510)], hsp.hit_range_all
        )
        # first qresult, fifth hit, second hsp
        hsp = qresult[4].hsps[1]
        self.assertEqual(39, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(0, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(-1, hsp[0].query_strand)
        self.assertEqual(10, hsp.query_start)
        self.assertEqual(54017130, hsp.hit_start)
        self.assertEqual(49, hsp.query_end)
        self.assertEqual(54017169, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([39], hsp.query_span_all)
        self.assertEqual([39], hsp.hit_span_all)
        self.assertEqual([(10, 49)], hsp.query_range_all)
        self.assertEqual([(54017130, 54017169)], hsp.hit_range_all)
        # first qresult, fifth hit, third hsp
        hsp = qresult[4].hsps[2]
        self.assertEqual(36, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(3, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(-1, hsp[0].query_strand)
        self.assertEqual(10, hsp.query_start)
        self.assertEqual(553742, hsp.hit_start)
        self.assertEqual(49, hsp.query_end)
        self.assertEqual(553781, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([39], hsp.query_span_all)
        self.assertEqual([39], hsp.hit_span_all)
        self.assertEqual([(10, 49)], hsp.query_range_all)
        self.assertEqual([(553742, 553781)], hsp.hit_range_all)

    def test_psl_34_005(self, testf="psl_34_005.psl", pslx=False):
        """Test parsing blat output (psl_34_005.psl)."""
        blat_file = get_file(testf)
        self.qresults = list(parse(blat_file, FMT, pslx=pslx))
        self.assertEqual(2, len(self.qresults))
        # check common attributes
        for qresult in self.qresults:
            for hit in qresult:
                self.assertEqual(qresult.id, hit.query_id)
                for hsp in hit:
                    self.assertEqual(hit.id, hsp.hit_id)
                    self.assertEqual(qresult.id, hsp.query_id)

        # test first qresult
        qresult = self.qresults[0]
        self.assertEqual("hg18_dna", qresult.id)
        self.assertEqual("blat", qresult.program)
        self.assertEqual(33, qresult.seq_len)
        self.assertEqual(3, len(qresult))
        # first qresult, first hit
        hit = qresult[0]
        self.assertEqual("chr4", hit.id)
        self.assertEqual(191154276, hit.seq_len)
        self.assertEqual(1, len(hit.hsps))
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(16, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(0, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(11, hsp.query_start)
        self.assertEqual(61646095, hsp.hit_start)
        self.assertEqual(27, hsp.query_end)
        self.assertEqual(61646111, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([16], hsp.query_span_all)
        self.assertEqual([16], hsp.hit_span_all)
        self.assertEqual([(11, 27)], hsp.query_range_all)
        self.assertEqual([(61646095, 61646111)], hsp.hit_range_all)
        # first qresult, second hit
        hit = qresult[1]
        self.assertEqual("chr1", hit.id)
        self.assertEqual(249250621, hit.seq_len)
        self.assertEqual(1, len(hit.hsps))
        # first qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(33, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(0, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(10271783, hsp.hit_start)
        self.assertEqual(33, hsp.query_end)
        self.assertEqual(10271816, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([33], hsp.query_span_all)
        self.assertEqual([33], hsp.hit_span_all)
        self.assertEqual([(0, 33)], hsp.query_range_all)
        self.assertEqual([(10271783, 10271816)], hsp.hit_range_all)
        # first qresult, third hit
        hit = qresult[2]
        self.assertEqual("chr2", hit.id)
        self.assertEqual(243199373, hit.seq_len)
        self.assertEqual(1, len(hit.hsps))
        # first qresult, third hit, first hsp
        hsp = qresult[2].hsps[0]
        self.assertEqual(17, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(0, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(-1, hsp[0].query_strand)
        self.assertEqual(8, hsp.query_start)
        self.assertEqual(53575980, hsp.hit_start)
        self.assertEqual(25, hsp.query_end)
        self.assertEqual(53575997, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([17], hsp.query_span_all)
        self.assertEqual([17], hsp.hit_span_all)
        self.assertEqual([(8, 25)], hsp.query_range_all)
        self.assertEqual([(53575980, 53575997)], hsp.hit_range_all)

        # test second qresult
        qresult = self.qresults[1]
        self.assertEqual("hg19_dna", qresult.id)
        self.assertEqual("blat", qresult.program)
        self.assertEqual(50, qresult.seq_len)
        self.assertEqual(10, len(qresult))
        # second qresult, first hit
        hit = qresult[0]
        self.assertEqual("chr9", hit.id)
        self.assertEqual(141213431, hit.seq_len)
        self.assertEqual(1, len(hit.hsps))
        # second qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(38, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(3, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(9, hsp.query_start)
        self.assertEqual(85737865, hsp.hit_start)
        self.assertEqual(50, hsp.query_end)
        self.assertEqual(85737906, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([41], hsp.query_span_all)
        self.assertEqual([41], hsp.hit_span_all)
        self.assertEqual([(9, 50)], hsp.query_range_all)
        self.assertEqual([(85737865, 85737906)], hsp.hit_range_all)
        # second qresult, second hit
        hit = qresult[1]
        self.assertEqual("chr8", hit.id)
        self.assertEqual(146364022, hit.seq_len)
        self.assertEqual(1, len(hit.hsps))
        # second qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(41, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(0, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(8, hsp.query_start)
        self.assertEqual(95160479, hsp.hit_start)
        self.assertEqual(49, hsp.query_end)
        self.assertEqual(95160520, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([41], hsp.query_span_all)
        self.assertEqual([41], hsp.hit_span_all)
        self.assertEqual([(8, 49)], hsp.query_range_all)
        self.assertEqual([(95160479, 95160520)], hsp.hit_range_all)
        # second qresult, third hit
        hit = qresult[2]
        self.assertEqual("chr22", hit.id)
        self.assertEqual(51304566, hit.seq_len)
        self.assertEqual(2, len(hit.hsps))
        # second qresult, third hit, first hsp
        hsp = qresult[2].hsps[0]
        self.assertEqual(33, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(3, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(11, hsp.query_start)
        self.assertEqual(42144400, hsp.hit_start)
        self.assertEqual(47, hsp.query_end)
        self.assertEqual(42144436, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([36], hsp.query_span_all)
        self.assertEqual([36], hsp.hit_span_all)
        self.assertEqual([(11, 47)], hsp.query_range_all)
        self.assertEqual([(42144400, 42144436)], hsp.hit_range_all)
        # second qresult, third hit, second hsp
        hsp = qresult[2].hsps[1]
        self.assertEqual(35, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(2, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(-1, hsp[0].query_strand)
        self.assertEqual(12, hsp.query_start)
        self.assertEqual(48997405, hsp.hit_start)
        self.assertEqual(49, hsp.query_end)
        self.assertEqual(48997442, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([37], hsp.query_span_all)
        self.assertEqual([37], hsp.hit_span_all)
        self.assertEqual([(12, 49)], hsp.query_range_all)
        self.assertEqual([(48997405, 48997442)], hsp.hit_range_all)
        # second qresult, fourth hit
        hit = qresult[3]
        self.assertEqual("chr2", hit.id)
        self.assertEqual(243199373, hit.seq_len)
        self.assertEqual(2, len(hit.hsps))
        # second qresult, fourth hit, first hsp
        hsp = qresult[3].hsps[0]
        self.assertEqual(43, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(1, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(1, hsp.query_gapopen_num)
        self.assertEqual(4, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(1, hsp.query_start)
        self.assertEqual(183925984, hsp.hit_start)
        self.assertEqual(49, hsp.query_end)
        self.assertEqual(183926028, hsp.hit_end)
        self.assertEqual(2, len(hsp))
        self.assertEqual([6, 38], hsp.query_span_all)
        self.assertEqual([6, 38], hsp.hit_span_all)
        self.assertEqual([(1, 7), (11, 49)], hsp.query_range_all)
        self.assertEqual(
            [(183925984, 183925990), (183925990, 183926028)], hsp.hit_range_all
        )
        # second qresult, fourth hit, second hsp
        hsp = qresult[3].hsps[1]
        self.assertEqual(35, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(1, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(-1, hsp[0].query_strand)
        self.assertEqual(13, hsp.query_start)
        self.assertEqual(120641740, hsp.hit_start)
        self.assertEqual(49, hsp.query_end)
        self.assertEqual(120641776, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([36], hsp.query_span_all)
        self.assertEqual([36], hsp.hit_span_all)
        self.assertEqual([(13, 49)], hsp.query_range_all)
        self.assertEqual([(120641740, 120641776)], hsp.hit_range_all)
        # second qresult, fifth hit
        hit = qresult[4]
        self.assertEqual("chr19", hit.id)
        self.assertEqual(59128983, hit.seq_len)
        self.assertEqual(3, len(hit.hsps))
        # second qresult, fifth hit, first hsp
        hsp = qresult[4].hsps[0]
        self.assertEqual(34, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(2, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(1, hsp.hit_gapopen_num)
        self.assertEqual(134, hsp.hit_gap_num)
        self.assertEqual(1, hsp[0].query_strand)
        self.assertEqual(10, hsp.query_start)
        self.assertEqual(35483340, hsp.hit_start)
        self.assertEqual(46, hsp.query_end)
        self.assertEqual(35483510, hsp.hit_end)
        self.assertEqual(2, len(hsp))
        self.assertEqual([25, 11], hsp.query_span_all)
        self.assertEqual([25, 11], hsp.hit_span_all)
        self.assertEqual([(10, 35), (35, 46)], hsp.query_range_all)
        self.assertEqual(
            [(35483340, 35483365), (35483499, 35483510)], hsp.hit_range_all
        )
        # second qresult, fifth hit, second hsp
        hsp = qresult[4].hsps[1]
        self.assertEqual(39, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(0, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(-1, hsp[0].query_strand)
        self.assertEqual(10, hsp.query_start)
        self.assertEqual(54017130, hsp.hit_start)
        self.assertEqual(49, hsp.query_end)
        self.assertEqual(54017169, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([39], hsp.query_span_all)
        self.assertEqual([39], hsp.hit_span_all)
        self.assertEqual([(10, 49)], hsp.query_range_all)
        self.assertEqual([(54017130, 54017169)], hsp.hit_range_all)
        # second qresult, fifth hit, third hsp
        hsp = qresult[4].hsps[2]
        self.assertEqual(36, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(3, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(-1, hsp[0].query_strand)
        self.assertEqual(10, hsp.query_start)
        self.assertEqual(553742, hsp.hit_start)
        self.assertEqual(49, hsp.query_end)
        self.assertEqual(553781, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([39], hsp.query_span_all)
        self.assertEqual([39], hsp.hit_span_all)
        self.assertEqual([(10, 49)], hsp.query_range_all)
        self.assertEqual([(553742, 553781)], hsp.hit_range_all)

    def test_psl_35_001(self, testf="psl_35_001.psl", pslx=False):
        """Test parsing blat output (psl_35_001.psl)."""
        blat_file = get_file(testf)
        self.qresults = list(parse(blat_file, FMT, pslx=pslx))
        self.assertEqual(1, len(self.qresults))
        # check common attributes
        for qresult in self.qresults:
            for hit in qresult:
                self.assertEqual(qresult.id, hit.query_id)
                for hsp in hit:
                    self.assertEqual(hit.id, hsp.hit_id)
                    self.assertEqual(qresult.id, hsp.query_id)

        # test first qresult
        qresult = self.qresults[0]
        self.assertEqual("CAG33136.1", qresult.id)
        self.assertEqual("blat", qresult.program)
        self.assertEqual(230, qresult.seq_len)
        self.assertEqual(2, len(qresult))
        # first qresult, first hit
        hit = qresult[0]
        self.assertEqual("chr13", hit.id)
        self.assertEqual(114364328, hit.seq_len)
        self.assertEqual(6, len(hit.hsps))
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(52, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(0, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(0, hsp[0].query_strand)
        self.assertEqual(61, hsp.query_start)
        self.assertEqual(75566694, hsp.hit_start)
        self.assertEqual(113, hsp.query_end)
        self.assertEqual(75566850, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([52], hsp.query_span_all)
        self.assertEqual([156], hsp.hit_span_all)
        self.assertEqual([(61, 113)], hsp.query_range_all)
        self.assertEqual([(75566694, 75566850)], hsp.hit_range_all)

    def test_psl_35_002(self, testf="psl_35_002.psl", pslx=False):
        """Test parsing blat output (psl_35_002.psl)."""
        blat_file = get_file(testf)
        self.qresults = list(parse(blat_file, FMT, pslx=pslx))
        self.assertEqual(1, len(self.qresults))
        # check common attributes
        for qresult in self.qresults:
            for hit in qresult:
                self.assertEqual(qresult.id, hit.query_id)
                for hsp in hit:
                    self.assertEqual(hit.id, hsp.hit_id)
                    self.assertEqual(qresult.id, hsp.query_id)

        # test first qresult
        qresult = self.qresults[0]
        self.assertEqual("CAG33136.1", qresult.id)
        self.assertEqual("blat", qresult.program)
        self.assertEqual(230, qresult.seq_len)
        self.assertEqual(3, len(qresult))
        # first qresult, last hit
        hit = qresult[-1]
        self.assertEqual("KI537194", hit.id)
        self.assertEqual(37111980, hit.seq_len)
        self.assertEqual(1, len(hit.hsps))
        # # first qresult, last hit, first hsp
        hsp = hit.hsps[-1]
        self.assertEqual(204, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(6, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(1, hsp.query_gapopen_num)
        self.assertEqual(20, hsp.query_gap_num)
        self.assertEqual(1, hsp.hit_gapopen_num)
        self.assertEqual(1, hsp.hit_gap_num)
        self.assertEqual(0, hsp[0].query_strand)
        self.assertEqual(-1, hsp[0].hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(20872390, hsp.hit_start)
        self.assertEqual(230, hsp.query_end)
        self.assertEqual(20873021, hsp.hit_end)
        self.assertEqual(2, len(hsp))
        self.assertEqual([183, 27], hsp.query_span_all)
        self.assertEqual([549, 81], hsp.hit_span_all)
        self.assertEqual([(0, 183), (203, 230)], hsp.query_range_all)
        self.assertEqual(
            [(20872472, 20873021), (20872390, 20872471)], hsp.hit_range_all
        )


class BlatPslxCases(BlatPslCases):
    def test_pslx_34_001(self, testf="pslx_34_001.pslx"):
        """Test parsing blat output (pslx_34_001.pslx)."""
        BlatPslCases.test_psl_34_001(self, "pslx_34_001.pslx", pslx=True)

        # test first qresult
        qresult = self.qresults[0]
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual("aggtaaactgccttca", hsp.query_all[0].seq)
        self.assertEqual("aggtaaactgccttca", hsp.hit_all[0].seq)
        # first qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual("atgagcttccaaggtaaactgccttcaagattc", hsp.query_all[0].seq)
        self.assertEqual("atgagcttccaaggtaaactgccttcaagattc", hsp.hit_all[0].seq)
        # first qresult, third hit, first hsp
        hsp = qresult[2].hsps[0]
        self.assertEqual("aaggcagtttaccttgg", hsp.query_all[0].seq)
        self.assertEqual("aaggcagtttaccttgg", hsp.hit_all[0].seq)

        # test second qresult
        qresult = self.qresults[1]
        # second qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(
            "acaaaggggctgggcgtggtggctcacacctgtaatcccaa", hsp.query_all[0].seq
        )
        self.assertEqual(
            "acaaaggggctgggcgcagtggctcacgcctgtaatcccaa", hsp.hit_all[0].seq
        )
        # second qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(
            "cacaaaggggctgggcgtggtggctcacacctgtaatccca", hsp.query_all[0].seq
        )
        self.assertEqual(
            "cacaaaggggctgggcgtggtggctcacacctgtaatccca", hsp.hit_all[0].seq
        )
        # second qresult, third hit, first hsp
        hsp = qresult[2].hsps[0]
        self.assertEqual("aaaggggctgggcgtggtggctcacacctgtaatcc", hsp.query_all[0].seq)
        self.assertEqual("aaaggggctgggcgtggtagctcatgcctgtaatcc", hsp.hit_all[0].seq)
        # second qresult, third hit, second hsp
        hsp = qresult[2].hsps[1]
        self.assertEqual("tgggattacaggtgtgagccaccacgcccagcccctt", hsp.query_all[0].seq)
        self.assertEqual("tgggattacaggcgggagccaccacgcccagcccctt", hsp.hit_all[0].seq)
        # second qresult, fourth hit, first hsp
        hsp = qresult[3].hsps[0]
        self.assertEqual("aaaaat", hsp.query_all[0].seq)
        self.assertEqual("aaaaat", hsp.hit_all[0].seq)
        self.assertEqual("aaaggggctgggcgtggtggctcacacctgtaatccca", hsp.query_all[1].seq)
        self.assertEqual("aaaggggctgggcgtggtggctcacgcctgtaatccca", hsp.hit_all[1].seq)
        # second qresult, fourth hit, second hsp
        hsp = qresult[3].hsps[1]
        self.assertEqual("tgggattacaggtgtgagccaccacgcccagcccct", hsp.query_all[0].seq)
        self.assertEqual("tgggattacaggcgtgagccaccacgcccagcccct", hsp.hit_all[0].seq)
        # second qresult, fifth hit, first hsp
        hsp = qresult[4].hsps[0]
        self.assertEqual("caaaggggctgggcgtggtggctca", hsp.query_all[0].seq)
        self.assertEqual("caaaggggctgggcgtagtggctga", hsp.hit_all[0].seq)
        self.assertEqual("cacctgtaatc", hsp.query_all[1].seq)
        self.assertEqual("cacctgtaatc", hsp.hit_all[1].seq)
        # second qresult, fifth hit, second hsp
        hsp = qresult[4].hsps[1]
        self.assertEqual(
            "tgggattacaggtgtgagccaccacgcccagcccctttg", hsp.query_all[0].seq
        )
        self.assertEqual("tgggattacaggtgtgagccaccacgcccagcccctttg", hsp.hit_all[0].seq)
        # second qresult, fifth hit, third hsp
        hsp = qresult[4].hsps[2]
        self.assertEqual(
            "tgggattacaggtgtgagccaccacgcccagcccctttg", hsp.query_all[0].seq
        )
        self.assertEqual("tgggatgacaggggtgaggcaccacgcccagcccctttg", hsp.hit_all[0].seq)

    def test_pslx_34_002(self, testf="pslx_34_002.pslx"):
        """Test parsing blat output (pslx_34_002.pslx)."""
        BlatPslCases.test_psl_34_002(self, "pslx_34_002.pslx", pslx=True)

    def test_pslx_34_003(self, testf="pslx_34_003.pslx"):
        """Test parsing blat output (pslx_34_003.pslx)."""
        BlatPslCases.test_psl_34_003(self, "pslx_34_003.pslx", pslx=True)

        # test first qresult
        qresult = self.qresults[0]
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual("aggtaaactgccttca", hsp.query_all[0].seq)
        self.assertEqual("aggtaaactgccttca", hsp.hit_all[0].seq)
        # first qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual("atgagcttccaaggtaaactgccttcaagattc", hsp.query_all[0].seq)
        self.assertEqual("atgagcttccaaggtaaactgccttcaagattc", hsp.hit_all[0].seq)
        # first qresult, third hit, first hsp
        hsp = qresult[2].hsps[0]
        self.assertEqual("aaggcagtttaccttgg", hsp.query_all[0].seq)
        self.assertEqual("aaggcagtttaccttgg", hsp.hit_all[0].seq)

    def test_pslx_34_004(self, testf="pslx_34_004.pslx"):
        """Test parsing blat output (pslx_34_004.pslx)."""
        BlatPslCases.test_psl_34_004(self, "pslx_34_004.pslx", pslx=True)

        # test first qresult
        qresult = self.qresults[0]
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(
            "acaaaggggctgggcgtggtggctcacacctgtaatcccaa", hsp.query_all[0].seq
        )
        self.assertEqual(
            "acaaaggggctgggcgcagtggctcacgcctgtaatcccaa", hsp.hit_all[0].seq
        )
        # first qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(
            "cacaaaggggctgggcgtggtggctcacacctgtaatccca", hsp.query_all[0].seq
        )
        self.assertEqual(
            "cacaaaggggctgggcgtggtggctcacacctgtaatccca", hsp.hit_all[0].seq
        )
        # first qresult, third hit, first hsp
        hsp = qresult[2].hsps[0]
        self.assertEqual("aaaggggctgggcgtggtggctcacacctgtaatcc", hsp.query_all[0].seq)
        self.assertEqual("aaaggggctgggcgtggtagctcatgcctgtaatcc", hsp.hit_all[0].seq)
        # first qresult, third hit, second hsp
        hsp = qresult[2].hsps[1]
        self.assertEqual("tgggattacaggtgtgagccaccacgcccagcccctt", hsp.query_all[0].seq)
        self.assertEqual("tgggattacaggcgggagccaccacgcccagcccctt", hsp.hit_all[0].seq)
        # first qresult, fourth hit, first hsp
        hsp = qresult[3].hsps[0]
        self.assertEqual("aaaaat", hsp.query_all[0].seq)
        self.assertEqual("aaaaat", hsp.hit_all[0].seq)
        self.assertEqual("aaaggggctgggcgtggtggctcacacctgtaatccca", hsp.query_all[1].seq)
        self.assertEqual("aaaggggctgggcgtggtggctcacgcctgtaatccca", hsp.hit_all[1].seq)
        # first qresult, fourth hit, second hsp
        hsp = qresult[3].hsps[1]
        self.assertEqual("tgggattacaggtgtgagccaccacgcccagcccct", hsp.query_all[0].seq)
        self.assertEqual("tgggattacaggcgtgagccaccacgcccagcccct", hsp.hit_all[0].seq)
        # first qresult, fifth hit, first hsp
        hsp = qresult[4].hsps[0]
        self.assertEqual("caaaggggctgggcgtggtggctca", hsp.query_all[0].seq)
        self.assertEqual("caaaggggctgggcgtagtggctga", hsp.hit_all[0].seq)
        self.assertEqual("cacctgtaatc", hsp.query_all[1].seq)
        self.assertEqual("cacctgtaatc", hsp.hit_all[1].seq)
        # first qresult, fifth hit, second hsp
        hsp = qresult[4].hsps[1]
        self.assertEqual(
            "tgggattacaggtgtgagccaccacgcccagcccctttg", hsp.query_all[0].seq
        )
        self.assertEqual("tgggattacaggtgtgagccaccacgcccagcccctttg", hsp.hit_all[0].seq)
        # first qresult, fifth hit, third hsp
        hsp = qresult[4].hsps[2]
        self.assertEqual(
            "tgggattacaggtgtgagccaccacgcccagcccctttg", hsp.query_all[0].seq
        )
        self.assertEqual("tgggatgacaggggtgaggcaccacgcccagcccctttg", hsp.hit_all[0].seq)

    def test_pslx_34_005(self, testf="pslx_34_005.pslx"):
        """Test parsing blat output (pslx_34_005.pslx)."""
        BlatPslCases.test_psl_34_005(self, "pslx_34_005.pslx", pslx=True)

        # test first qresult
        qresult = self.qresults[0]
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual("aggtaaactgccttca", hsp.query_all[0].seq)
        self.assertEqual("aggtaaactgccttca", hsp.hit_all[0].seq)
        # first qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual("atgagcttccaaggtaaactgccttcaagattc", hsp.query_all[0].seq)
        self.assertEqual("atgagcttccaaggtaaactgccttcaagattc", hsp.hit_all[0].seq)
        # first qresult, third hit, first hsp
        hsp = qresult[2].hsps[0]
        self.assertEqual("aaggcagtttaccttgg", hsp.query_all[0].seq)
        self.assertEqual("aaggcagtttaccttgg", hsp.hit_all[0].seq)

        # test second qresult
        qresult = self.qresults[1]
        # second qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(
            "acaaaggggctgggcgtggtggctcacacctgtaatcccaa", hsp.query_all[0].seq
        )
        self.assertEqual(
            "acaaaggggctgggcgcagtggctcacgcctgtaatcccaa", hsp.hit_all[0].seq
        )
        # second qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(
            "cacaaaggggctgggcgtggtggctcacacctgtaatccca", hsp.query_all[0].seq
        )
        self.assertEqual(
            "cacaaaggggctgggcgtggtggctcacacctgtaatccca", hsp.hit_all[0].seq
        )
        # second qresult, third hit, first hsp
        hsp = qresult[2].hsps[0]
        self.assertEqual("aaaggggctgggcgtggtggctcacacctgtaatcc", hsp.query_all[0].seq)
        self.assertEqual("aaaggggctgggcgtggtagctcatgcctgtaatcc", hsp.hit_all[0].seq)
        # second qresult, third hit, second hsp
        hsp = qresult[2].hsps[1]
        self.assertEqual("tgggattacaggtgtgagccaccacgcccagcccctt", hsp.query_all[0].seq)
        self.assertEqual("tgggattacaggcgggagccaccacgcccagcccctt", hsp.hit_all[0].seq)
        # second qresult, fourth hit, first hsp
        hsp = qresult[3].hsps[0]
        self.assertEqual("aaaaat", hsp.query_all[0].seq)
        self.assertEqual("aaaaat", hsp.hit_all[0].seq)
        self.assertEqual("aaaggggctgggcgtggtggctcacacctgtaatccca", hsp.query_all[1].seq)
        self.assertEqual("aaaggggctgggcgtggtggctcacgcctgtaatccca", hsp.hit_all[1].seq)
        # second qresult, fourth hit, second hsp
        hsp = qresult[3].hsps[1]
        self.assertEqual("tgggattacaggtgtgagccaccacgcccagcccct", hsp.query_all[0].seq)
        self.assertEqual("tgggattacaggcgtgagccaccacgcccagcccct", hsp.hit_all[0].seq)
        # second qresult, fifth hit, first hsp
        hsp = qresult[4].hsps[0]
        self.assertEqual("caaaggggctgggcgtggtggctca", hsp.query_all[0].seq)
        self.assertEqual("caaaggggctgggcgtagtggctga", hsp.hit_all[0].seq)
        self.assertEqual("cacctgtaatc", hsp.query_all[1].seq)
        self.assertEqual("cacctgtaatc", hsp.hit_all[1].seq)
        # second qresult, fifth hit, second hsp
        hsp = qresult[4].hsps[1]
        self.assertEqual(
            "tgggattacaggtgtgagccaccacgcccagcccctttg", hsp.query_all[0].seq
        )
        self.assertEqual("tgggattacaggtgtgagccaccacgcccagcccctttg", hsp.hit_all[0].seq)
        # second qresult, fifth hit, third hsp
        hsp = qresult[4].hsps[2]
        self.assertEqual(
            "tgggattacaggtgtgagccaccacgcccagcccctttg", hsp.query_all[0].seq
        )
        self.assertEqual("tgggatgacaggggtgaggcaccacgcccagcccctttg", hsp.hit_all[0].seq)

    def test_pslx_35_002(self, testf="pslx_35_002.pslx"):
        """Test parsing blat output (pslx_35_002.pslx)."""
        BlatPslCases.test_psl_35_002(self, "pslx_35_002.pslx", pslx=True)

        # first qresult, last hit, first hsp
        qresult = self.qresults[0]
        hsp = qresult[-1].hsps[0]

        self.assertEqual(
            "MEGQRWLPLEANPEVTNQFLKQLGLHPNWQFVDVY", hsp.query_all[0].seq[:35]
        )
        self.assertEqual(
            "ETSAHEGQTEAPSIDEKVDLHFIALVHVDGHLYEL", hsp.query_all[0].seq[-35:]
        )
        self.assertEqual("DAIEVCKKFMERDPDELRFNAIALSAA", hsp.query_all[1].seq)

        self.assertEqual("MESQRWLPLEANPEVTNQFLKQLGLHPNWQCVDVY", hsp.hit_all[0].seq[:35])
        self.assertEqual(
            "ETSAHEGQTEAPNIDEKVDLHFIALVHVDGHLYEL", hsp.hit_all[0].seq[-35:]
        )
        self.assertEqual("DAIEVCKKFMERDPDELRFNAIALSAA", hsp.hit_all[1].seq)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
