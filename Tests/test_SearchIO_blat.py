# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO BlatIO parsers."""

import os
import unittest

from Bio.SearchIO import parse

# test case files are in the Blast directory
TEST_DIR = 'Blat'


def get_file(filename):
    """Returns the path of a test file."""
    return os.path.join(TEST_DIR, filename)


class BlatPslCases(unittest.TestCase):

    def test_psl_34_001(self, testf='psl_34_001.psl', fmt='blat-psl'):
        """Test parsing blat output (psl_34_001.psl)"""
        blat_file = get_file(testf)
        self.qresults = list(parse(blat_file, fmt))
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
        self.assertEqual('hg18_dna', qresult.id)
        self.assertEqual('blat', qresult.program)
        self.assertEqual(33, qresult.seq_len)
        self.assertEqual(3, len(qresult))
        # first qresult, first hit
        hit = qresult[0]
        self.assertEqual('chr4', hit.id)
        self.assertEqual(191154276, hit.seq_len)
        self.assertEqual(1, len(hit.batch_hsps))
        # first qresult, first hit, first hsp
        hsp = qresult[0].batch_hsps[0]
        self.assertEqual(16, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(0, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(11, hsp.query_start)
        self.assertEqual(61646095, hsp.hit_start)
        self.assertEqual(27, hsp.query_end)
        self.assertEqual(61646111, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([16], hsp.query_spans)
        self.assertEqual([16], hsp.hit_spans)
        self.assertEqual([(11, 27)], hsp.query_ranges)
        self.assertEqual([(61646095, 61646111)], hsp.hit_ranges)
        # first qresult, second hit
        hit = qresult[1]
        self.assertEqual('chr1', hit.id)
        self.assertEqual(249250621, hit.seq_len)
        self.assertEqual(1, len(hit.batch_hsps))
        # first qresult, second hit, first hsp
        hsp = qresult[1].batch_hsps[0]
        self.assertEqual(33, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(0, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(10271783, hsp.hit_start)
        self.assertEqual(33, hsp.query_end)
        self.assertEqual(10271816, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([33], hsp.query_spans)
        self.assertEqual([33], hsp.hit_spans)
        self.assertEqual([(0, 33)], hsp.query_ranges)
        self.assertEqual([(10271783, 10271816)], hsp.hit_ranges)
        # first qresult, third hit
        hit = qresult[2]
        self.assertEqual('chr2', hit.id)
        self.assertEqual(243199373, hit.seq_len)
        self.assertEqual(1, len(hit.batch_hsps))
        # first qresult, third hit, first hsp
        hsp = qresult[2].batch_hsps[0]
        self.assertEqual(17, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(0, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(-1, hsp.query_strand)
        self.assertEqual(8, hsp.query_start)
        self.assertEqual(53575980, hsp.hit_start)
        self.assertEqual(25, hsp.query_end)
        self.assertEqual(53575997, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([17], hsp.query_spans)
        self.assertEqual([17], hsp.hit_spans)
        self.assertEqual([(8, 25)], hsp.query_ranges)
        self.assertEqual([(53575980, 53575997)], hsp.hit_ranges)

        # test second qresult
        qresult = self.qresults[1]
        self.assertEqual('hg19_dna', qresult.id)
        self.assertEqual('blat', qresult.program)
        self.assertEqual(50, qresult.seq_len)
        self.assertEqual(10, len(qresult))
        # second qresult, first hit
        hit = qresult[0]
        self.assertEqual('chr9', hit.id)
        self.assertEqual(141213431, hit.seq_len)
        self.assertEqual(1, len(hit.batch_hsps))
        # second qresult, first hit, first hsp
        hsp = qresult[0].batch_hsps[0]
        self.assertEqual(38, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(3, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(9, hsp.query_start)
        self.assertEqual(85737865, hsp.hit_start)
        self.assertEqual(50, hsp.query_end)
        self.assertEqual(85737906, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([41], hsp.query_spans)
        self.assertEqual([41], hsp.hit_spans)
        self.assertEqual([(9, 50)], hsp.query_ranges)
        self.assertEqual([(85737865, 85737906)], hsp.hit_ranges)
        # second qresult, second hit
        hit = qresult[1]
        self.assertEqual('chr8', hit.id)
        self.assertEqual(146364022, hit.seq_len)
        self.assertEqual(1, len(hit.batch_hsps))
        # second qresult, second hit, first hsp
        hsp = qresult[1].batch_hsps[0]
        self.assertEqual(41, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(0, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(8, hsp.query_start)
        self.assertEqual(95160479, hsp.hit_start)
        self.assertEqual(49, hsp.query_end)
        self.assertEqual(95160520, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([41], hsp.query_spans)
        self.assertEqual([41], hsp.hit_spans)
        self.assertEqual([(8, 49)], hsp.query_ranges)
        self.assertEqual([(95160479, 95160520)], hsp.hit_ranges)
        # second qresult, third hit
        hit = qresult[2]
        self.assertEqual('chr22', hit.id)
        self.assertEqual(51304566, hit.seq_len)
        self.assertEqual(2, len(hit.batch_hsps))
        # second qresult, third hit, first hsp
        hsp = qresult[2].batch_hsps[0]
        self.assertEqual(33, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(3, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(11, hsp.query_start)
        self.assertEqual(42144400, hsp.hit_start)
        self.assertEqual(47, hsp.query_end)
        self.assertEqual(42144436, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([36], hsp.query_spans)
        self.assertEqual([36], hsp.hit_spans)
        self.assertEqual([(11, 47)], hsp.query_ranges)
        self.assertEqual([(42144400, 42144436)], hsp.hit_ranges)
        # second qresult, third hit, second hsp
        hsp = qresult[2].batch_hsps[1]
        self.assertEqual(35, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(2, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(-1, hsp.query_strand)
        self.assertEqual(12, hsp.query_start)
        self.assertEqual(48997405, hsp.hit_start)
        self.assertEqual(49, hsp.query_end)
        self.assertEqual(48997442, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([37], hsp.query_spans)
        self.assertEqual([37], hsp.hit_spans)
        self.assertEqual([(12, 49)], hsp.query_ranges)
        self.assertEqual([(48997405, 48997442)], hsp.hit_ranges)
        # second qresult, fourth hit
        hit = qresult[3]
        self.assertEqual('chr2', hit.id)
        self.assertEqual(243199373, hit.seq_len)
        self.assertEqual(2, len(hit.batch_hsps))
        # second qresult, fourth hit, first hsp
        hsp = qresult[3].batch_hsps[0]
        self.assertEqual(43, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(1, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(1, hsp.query_gapopen_num)
        self.assertEqual(4, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(1, hsp.query_start)
        self.assertEqual(183925984, hsp.hit_start)
        self.assertEqual(49, hsp.query_end)
        self.assertEqual(183926028, hsp.hit_end)
        self.assertEqual(2, len(hsp))
        self.assertEqual([6, 38], hsp.query_spans)
        self.assertEqual([6, 38], hsp.hit_spans)
        self.assertEqual([(1, 7), (11, 49)], hsp.query_ranges)
        self.assertEqual([(183925984, 183925990), (183925990, 183926028)], hsp.hit_ranges)
        # second qresult, fourth hit, second hsp
        hsp = qresult[3].batch_hsps[1]
        self.assertEqual(35, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(1, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(-1, hsp.query_strand)
        self.assertEqual(13, hsp.query_start)
        self.assertEqual(120641740, hsp.hit_start)
        self.assertEqual(49, hsp.query_end)
        self.assertEqual(120641776, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([36], hsp.query_spans)
        self.assertEqual([36], hsp.hit_spans)
        self.assertEqual([(13, 49)], hsp.query_ranges)
        self.assertEqual([(120641740, 120641776)], hsp.hit_ranges)
        # second qresult, fifth hit
        hit = qresult[4]
        self.assertEqual('chr19', hit.id)
        self.assertEqual(59128983, hit.seq_len)
        self.assertEqual(3, len(hit.batch_hsps))
        # second qresult, fifth hit, first hsp
        hsp = qresult[4].batch_hsps[0]
        self.assertEqual(34, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(2, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(1, hsp.hit_gapopen_num)
        self.assertEqual(134, hsp.hit_gap_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(10, hsp.query_start)
        self.assertEqual(35483340, hsp.hit_start)
        self.assertEqual(46, hsp.query_end)
        self.assertEqual(35483510, hsp.hit_end)
        self.assertEqual(2, len(hsp))
        self.assertEqual([25, 11], hsp.query_spans)
        self.assertEqual([25, 11], hsp.hit_spans)
        self.assertEqual([(10, 35), (35, 46)], hsp.query_ranges)
        self.assertEqual([(35483340, 35483365), (35483499, 35483510)], hsp.hit_ranges)
        # second qresult, fifth hit, second hsp
        hsp = qresult[4].batch_hsps[1]
        self.assertEqual(39, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(0, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(-1, hsp.query_strand)
        self.assertEqual(10, hsp.query_start)
        self.assertEqual(54017130, hsp.hit_start)
        self.assertEqual(49, hsp.query_end)
        self.assertEqual(54017169, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([39], hsp.query_spans)
        self.assertEqual([39], hsp.hit_spans)
        self.assertEqual([(10, 49)], hsp.query_ranges)
        self.assertEqual([(54017130, 54017169)], hsp.hit_ranges)
        # second qresult, fifth hit, third hsp
        hsp = qresult[4].batch_hsps[2]
        self.assertEqual(36, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(3, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(-1, hsp.query_strand)
        self.assertEqual(10, hsp.query_start)
        self.assertEqual(553742, hsp.hit_start)
        self.assertEqual(49, hsp.query_end)
        self.assertEqual(553781, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([39], hsp.query_spans)
        self.assertEqual([39], hsp.hit_spans)
        self.assertEqual([(10, 49)], hsp.query_ranges)
        self.assertEqual([(553742, 553781)], hsp.hit_ranges)

    def test_psl_34_002(self, testf='psl_34_002.psl', fmt='blat-psl'):
        """Test parsing blat output (psl_34_001.psl)"""
        blat_file = get_file(testf)
        self.qresults = list(parse(blat_file, fmt))
        self.assertEqual(0, len(self.qresults))

    def test_psl_34_003(self, testf='psl_34_003.psl', fmt='blat-psl'):
        """Test parsing blat output (psl_34_003.psl)"""
        blat_file = get_file(testf)
        self.qresults = list(parse(blat_file, fmt))
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
        self.assertEqual('hg18_dna', qresult.id)
        self.assertEqual('blat', qresult.program)
        self.assertEqual(33, qresult.seq_len)
        self.assertEqual(3, len(qresult))
        # first qresult, first hit
        hit = qresult[0]
        self.assertEqual('chr4', hit.id)
        self.assertEqual(191154276, hit.seq_len)
        self.assertEqual(1, len(hit.batch_hsps))
        # first qresult, first hit, first hsp
        hsp = qresult[0].batch_hsps[0]
        self.assertEqual(16, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(0, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(11, hsp.query_start)
        self.assertEqual(61646095, hsp.hit_start)
        self.assertEqual(27, hsp.query_end)
        self.assertEqual(61646111, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([16], hsp.query_spans)
        self.assertEqual([16], hsp.hit_spans)
        self.assertEqual([(11, 27)], hsp.query_ranges)
        self.assertEqual([(61646095, 61646111)], hsp.hit_ranges)
        # first qresult, second hit
        hit = qresult[1]
        self.assertEqual('chr1', hit.id)
        self.assertEqual(249250621, hit.seq_len)
        self.assertEqual(1, len(hit.batch_hsps))
        # first qresult, second hit, first hsp
        hsp = qresult[1].batch_hsps[0]
        self.assertEqual(33, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(0, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(10271783, hsp.hit_start)
        self.assertEqual(33, hsp.query_end)
        self.assertEqual(10271816, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([33], hsp.query_spans)
        self.assertEqual([33], hsp.hit_spans)
        self.assertEqual([(0, 33)], hsp.query_ranges)
        self.assertEqual([(10271783, 10271816)], hsp.hit_ranges)
        # first qresult, third hit
        hit = qresult[2]
        self.assertEqual('chr2', hit.id)
        self.assertEqual(243199373, hit.seq_len)
        self.assertEqual(1, len(hit.batch_hsps))
        # first qresult, third hit, first hsp
        hsp = qresult[2].batch_hsps[0]
        self.assertEqual(17, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(0, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(-1, hsp.query_strand)
        self.assertEqual(8, hsp.query_start)
        self.assertEqual(53575980, hsp.hit_start)
        self.assertEqual(25, hsp.query_end)
        self.assertEqual(53575997, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([17], hsp.query_spans)
        self.assertEqual([17], hsp.hit_spans)
        self.assertEqual([(8, 25)], hsp.query_ranges)
        self.assertEqual([(53575980, 53575997)], hsp.hit_ranges)

    def test_psl_34_004(self, testf='psl_34_004.psl', fmt='blat-psl'):
        """Test parsing blat output (psl_34_004.psl)"""
        blat_file = get_file(testf)
        self.qresults = list(parse(blat_file, fmt))
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
        self.assertEqual('hg19_dna', qresult.id)
        self.assertEqual('blat', qresult.program)
        self.assertEqual(50, qresult.seq_len)
        self.assertEqual(10, len(qresult))
        # first qresult, first hit
        hit = qresult[0]
        self.assertEqual('chr9', hit.id)
        self.assertEqual(141213431, hit.seq_len)
        self.assertEqual(1, len(hit.batch_hsps))
        # first qresult, first hit, first hsp
        hsp = qresult[0].batch_hsps[0]
        self.assertEqual(38, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(3, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(9, hsp.query_start)
        self.assertEqual(85737865, hsp.hit_start)
        self.assertEqual(50, hsp.query_end)
        self.assertEqual(85737906, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([41], hsp.query_spans)
        self.assertEqual([41], hsp.hit_spans)
        self.assertEqual([(9, 50)], hsp.query_ranges)
        self.assertEqual([(85737865, 85737906)], hsp.hit_ranges)
        # first qresult, second hit
        hit = qresult[1]
        self.assertEqual('chr8', hit.id)
        self.assertEqual(146364022, hit.seq_len)
        self.assertEqual(1, len(hit.batch_hsps))
        # first qresult, second hit, first hsp
        hsp = qresult[1].batch_hsps[0]
        self.assertEqual(41, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(0, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(8, hsp.query_start)
        self.assertEqual(95160479, hsp.hit_start)
        self.assertEqual(49, hsp.query_end)
        self.assertEqual(95160520, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([41], hsp.query_spans)
        self.assertEqual([41], hsp.hit_spans)
        self.assertEqual([(8, 49)], hsp.query_ranges)
        self.assertEqual([(95160479, 95160520)], hsp.hit_ranges)
        # first qresult, third hit
        hit = qresult[2]
        self.assertEqual('chr22', hit.id)
        self.assertEqual(51304566, hit.seq_len)
        self.assertEqual(2, len(hit.batch_hsps))
        # first qresult, third hit, first hsp
        hsp = qresult[2].batch_hsps[0]
        self.assertEqual(33, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(3, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(11, hsp.query_start)
        self.assertEqual(42144400, hsp.hit_start)
        self.assertEqual(47, hsp.query_end)
        self.assertEqual(42144436, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([36], hsp.query_spans)
        self.assertEqual([36], hsp.hit_spans)
        self.assertEqual([(11, 47)], hsp.query_ranges)
        self.assertEqual([(42144400, 42144436)], hsp.hit_ranges)
        # first qresult, third hit, second hsp
        hsp = qresult[2].batch_hsps[1]
        self.assertEqual(35, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(2, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(-1, hsp.query_strand)
        self.assertEqual(12, hsp.query_start)
        self.assertEqual(48997405, hsp.hit_start)
        self.assertEqual(49, hsp.query_end)
        self.assertEqual(48997442, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([37], hsp.query_spans)
        self.assertEqual([37], hsp.hit_spans)
        self.assertEqual([(12, 49)], hsp.query_ranges)
        self.assertEqual([(48997405, 48997442)], hsp.hit_ranges)
        # first qresult, fourth hit
        hit = qresult[3]
        self.assertEqual('chr2', hit.id)
        self.assertEqual(243199373, hit.seq_len)
        self.assertEqual(2, len(hit.batch_hsps))
        # first qresult, fourth hit, first hsp
        hsp = qresult[3].batch_hsps[0]
        self.assertEqual(43, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(1, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(1, hsp.query_gapopen_num)
        self.assertEqual(4, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(1, hsp.query_start)
        self.assertEqual(183925984, hsp.hit_start)
        self.assertEqual(49, hsp.query_end)
        self.assertEqual(183926028, hsp.hit_end)
        self.assertEqual(2, len(hsp))
        self.assertEqual([6, 38], hsp.query_spans)
        self.assertEqual([6, 38], hsp.hit_spans)
        self.assertEqual([(1, 7), (11, 49)], hsp.query_ranges)
        self.assertEqual([(183925984, 183925990), (183925990, 183926028)], hsp.hit_ranges)
        # first qresult, fourth hit, second hsp
        hsp = qresult[3].batch_hsps[1]
        self.assertEqual(35, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(1, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(-1, hsp.query_strand)
        self.assertEqual(13, hsp.query_start)
        self.assertEqual(120641740, hsp.hit_start)
        self.assertEqual(49, hsp.query_end)
        self.assertEqual(120641776, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([36], hsp.query_spans)
        self.assertEqual([36], hsp.hit_spans)
        self.assertEqual([(13, 49)], hsp.query_ranges)
        self.assertEqual([(120641740, 120641776)], hsp.hit_ranges)
        # first qresult, fifth hit
        hit = qresult[4]
        self.assertEqual('chr19', hit.id)
        self.assertEqual(59128983, hit.seq_len)
        self.assertEqual(3, len(hit.batch_hsps))
        # first qresult, fifth hit, first hsp
        hsp = qresult[4].batch_hsps[0]
        self.assertEqual(34, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(2, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(1, hsp.hit_gapopen_num)
        self.assertEqual(134, hsp.hit_gap_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(10, hsp.query_start)
        self.assertEqual(35483340, hsp.hit_start)
        self.assertEqual(46, hsp.query_end)
        self.assertEqual(35483510, hsp.hit_end)
        self.assertEqual(2, len(hsp))
        self.assertEqual([25, 11], hsp.query_spans)
        self.assertEqual([25, 11], hsp.hit_spans)
        self.assertEqual([(10, 35), (35, 46)], hsp.query_ranges)
        self.assertEqual([(35483340, 35483365), (35483499, 35483510)], hsp.hit_ranges)
        # first qresult, fifth hit, second hsp
        hsp = qresult[4].batch_hsps[1]
        self.assertEqual(39, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(0, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(-1, hsp.query_strand)
        self.assertEqual(10, hsp.query_start)
        self.assertEqual(54017130, hsp.hit_start)
        self.assertEqual(49, hsp.query_end)
        self.assertEqual(54017169, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([39], hsp.query_spans)
        self.assertEqual([39], hsp.hit_spans)
        self.assertEqual([(10, 49)], hsp.query_ranges)
        self.assertEqual([(54017130, 54017169)], hsp.hit_ranges)
        # first qresult, fifth hit, third hsp
        hsp = qresult[4].batch_hsps[2]
        self.assertEqual(36, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(3, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(-1, hsp.query_strand)
        self.assertEqual(10, hsp.query_start)
        self.assertEqual(553742, hsp.hit_start)
        self.assertEqual(49, hsp.query_end)
        self.assertEqual(553781, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([39], hsp.query_spans)
        self.assertEqual([39], hsp.hit_spans)
        self.assertEqual([(10, 49)], hsp.query_ranges)
        self.assertEqual([(553742, 553781)], hsp.hit_ranges)

    def test_psl_34_005(self, testf='psl_34_005.psl', fmt='blat-psl'):
        """Test parsing blat output (psl_34_005.psl)"""
        blat_file = get_file(testf)
        self.qresults = list(parse(blat_file, fmt))
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
        self.assertEqual('hg18_dna', qresult.id)
        self.assertEqual('blat', qresult.program)
        self.assertEqual(33, qresult.seq_len)
        self.assertEqual(3, len(qresult))
        # first qresult, first hit
        hit = qresult[0]
        self.assertEqual('chr4', hit.id)
        self.assertEqual(191154276, hit.seq_len)
        self.assertEqual(1, len(hit.batch_hsps))
        # first qresult, first hit, first hsp
        hsp = qresult[0].batch_hsps[0]
        self.assertEqual(16, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(0, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(11, hsp.query_start)
        self.assertEqual(61646095, hsp.hit_start)
        self.assertEqual(27, hsp.query_end)
        self.assertEqual(61646111, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([16], hsp.query_spans)
        self.assertEqual([16], hsp.hit_spans)
        self.assertEqual([(11, 27)], hsp.query_ranges)
        self.assertEqual([(61646095, 61646111)], hsp.hit_ranges)
        # first qresult, second hit
        hit = qresult[1]
        self.assertEqual('chr1', hit.id)
        self.assertEqual(249250621, hit.seq_len)
        self.assertEqual(1, len(hit.batch_hsps))
        # first qresult, second hit, first hsp
        hsp = qresult[1].batch_hsps[0]
        self.assertEqual(33, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(0, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(10271783, hsp.hit_start)
        self.assertEqual(33, hsp.query_end)
        self.assertEqual(10271816, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([33], hsp.query_spans)
        self.assertEqual([33], hsp.hit_spans)
        self.assertEqual([(0, 33)], hsp.query_ranges)
        self.assertEqual([(10271783, 10271816)], hsp.hit_ranges)
        # first qresult, third hit
        hit = qresult[2]
        self.assertEqual('chr2', hit.id)
        self.assertEqual(243199373, hit.seq_len)
        self.assertEqual(1, len(hit.batch_hsps))
        # first qresult, third hit, first hsp
        hsp = qresult[2].batch_hsps[0]
        self.assertEqual(17, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(0, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(-1, hsp.query_strand)
        self.assertEqual(8, hsp.query_start)
        self.assertEqual(53575980, hsp.hit_start)
        self.assertEqual(25, hsp.query_end)
        self.assertEqual(53575997, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([17], hsp.query_spans)
        self.assertEqual([17], hsp.hit_spans)
        self.assertEqual([(8, 25)], hsp.query_ranges)
        self.assertEqual([(53575980, 53575997)], hsp.hit_ranges)

        # test second qresult
        qresult = self.qresults[1]
        self.assertEqual('hg19_dna', qresult.id)
        self.assertEqual('blat', qresult.program)
        self.assertEqual(50, qresult.seq_len)
        self.assertEqual(10, len(qresult))
        # second qresult, first hit
        hit = qresult[0]
        self.assertEqual('chr9', hit.id)
        self.assertEqual(141213431, hit.seq_len)
        self.assertEqual(1, len(hit.batch_hsps))
        # second qresult, first hit, first hsp
        hsp = qresult[0].batch_hsps[0]
        self.assertEqual(38, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(3, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(9, hsp.query_start)
        self.assertEqual(85737865, hsp.hit_start)
        self.assertEqual(50, hsp.query_end)
        self.assertEqual(85737906, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([41], hsp.query_spans)
        self.assertEqual([41], hsp.hit_spans)
        self.assertEqual([(9, 50)], hsp.query_ranges)
        self.assertEqual([(85737865, 85737906)], hsp.hit_ranges)
        # second qresult, second hit
        hit = qresult[1]
        self.assertEqual('chr8', hit.id)
        self.assertEqual(146364022, hit.seq_len)
        self.assertEqual(1, len(hit.batch_hsps))
        # second qresult, second hit, first hsp
        hsp = qresult[1].batch_hsps[0]
        self.assertEqual(41, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(0, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(8, hsp.query_start)
        self.assertEqual(95160479, hsp.hit_start)
        self.assertEqual(49, hsp.query_end)
        self.assertEqual(95160520, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([41], hsp.query_spans)
        self.assertEqual([41], hsp.hit_spans)
        self.assertEqual([(8, 49)], hsp.query_ranges)
        self.assertEqual([(95160479, 95160520)], hsp.hit_ranges)
        # second qresult, third hit
        hit = qresult[2]
        self.assertEqual('chr22', hit.id)
        self.assertEqual(51304566, hit.seq_len)
        self.assertEqual(2, len(hit.batch_hsps))
        # second qresult, third hit, first hsp
        hsp = qresult[2].batch_hsps[0]
        self.assertEqual(33, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(3, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(11, hsp.query_start)
        self.assertEqual(42144400, hsp.hit_start)
        self.assertEqual(47, hsp.query_end)
        self.assertEqual(42144436, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([36], hsp.query_spans)
        self.assertEqual([36], hsp.hit_spans)
        self.assertEqual([(11, 47)], hsp.query_ranges)
        self.assertEqual([(42144400, 42144436)], hsp.hit_ranges)
        # second qresult, third hit, second hsp
        hsp = qresult[2].batch_hsps[1]
        self.assertEqual(35, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(2, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(-1, hsp.query_strand)
        self.assertEqual(12, hsp.query_start)
        self.assertEqual(48997405, hsp.hit_start)
        self.assertEqual(49, hsp.query_end)
        self.assertEqual(48997442, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([37], hsp.query_spans)
        self.assertEqual([37], hsp.hit_spans)
        self.assertEqual([(12, 49)], hsp.query_ranges)
        self.assertEqual([(48997405, 48997442)], hsp.hit_ranges)
        # second qresult, fourth hit
        hit = qresult[3]
        self.assertEqual('chr2', hit.id)
        self.assertEqual(243199373, hit.seq_len)
        self.assertEqual(2, len(hit.batch_hsps))
        # second qresult, fourth hit, first hsp
        hsp = qresult[3].batch_hsps[0]
        self.assertEqual(43, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(1, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(1, hsp.query_gapopen_num)
        self.assertEqual(4, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(1, hsp.query_start)
        self.assertEqual(183925984, hsp.hit_start)
        self.assertEqual(49, hsp.query_end)
        self.assertEqual(183926028, hsp.hit_end)
        self.assertEqual(2, len(hsp))
        self.assertEqual([6, 38], hsp.query_spans)
        self.assertEqual([6, 38], hsp.hit_spans)
        self.assertEqual([(1, 7), (11, 49)], hsp.query_ranges)
        self.assertEqual([(183925984, 183925990), (183925990, 183926028)], hsp.hit_ranges)
        # second qresult, fourth hit, second hsp
        hsp = qresult[3].batch_hsps[1]
        self.assertEqual(35, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(1, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(-1, hsp.query_strand)
        self.assertEqual(13, hsp.query_start)
        self.assertEqual(120641740, hsp.hit_start)
        self.assertEqual(49, hsp.query_end)
        self.assertEqual(120641776, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([36], hsp.query_spans)
        self.assertEqual([36], hsp.hit_spans)
        self.assertEqual([(13, 49)], hsp.query_ranges)
        self.assertEqual([(120641740, 120641776)], hsp.hit_ranges)
        # second qresult, fifth hit
        hit = qresult[4]
        self.assertEqual('chr19', hit.id)
        self.assertEqual(59128983, hit.seq_len)
        self.assertEqual(3, len(hit.batch_hsps))
        # second qresult, fifth hit, first hsp
        hsp = qresult[4].batch_hsps[0]
        self.assertEqual(34, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(2, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(1, hsp.hit_gapopen_num)
        self.assertEqual(134, hsp.hit_gap_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(10, hsp.query_start)
        self.assertEqual(35483340, hsp.hit_start)
        self.assertEqual(46, hsp.query_end)
        self.assertEqual(35483510, hsp.hit_end)
        self.assertEqual(2, len(hsp))
        self.assertEqual([25, 11], hsp.query_spans)
        self.assertEqual([25, 11], hsp.hit_spans)
        self.assertEqual([(10, 35), (35, 46)], hsp.query_ranges)
        self.assertEqual([(35483340, 35483365), (35483499, 35483510)], hsp.hit_ranges)
        # second qresult, fifth hit, second hsp
        hsp = qresult[4].batch_hsps[1]
        self.assertEqual(39, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(0, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(-1, hsp.query_strand)
        self.assertEqual(10, hsp.query_start)
        self.assertEqual(54017130, hsp.hit_start)
        self.assertEqual(49, hsp.query_end)
        self.assertEqual(54017169, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([39], hsp.query_spans)
        self.assertEqual([39], hsp.hit_spans)
        self.assertEqual([(10, 49)], hsp.query_ranges)
        self.assertEqual([(54017130, 54017169)], hsp.hit_ranges)
        # second qresult, fifth hit, third hsp
        hsp = qresult[4].batch_hsps[2]
        self.assertEqual(36, hsp.match_num)
        self.assertEqual(0, hsp.match_rep_num)
        self.assertEqual(3, hsp.mismatch_num)
        self.assertEqual(0, hsp.n_num)
        self.assertEqual(0, hsp.query_gapopen_num)
        self.assertEqual(0, hsp.query_gap_num)
        self.assertEqual(0, hsp.hit_gapopen_num)
        self.assertEqual(0, hsp.hit_gap_num)
        self.assertEqual(-1, hsp.query_strand)
        self.assertEqual(10, hsp.query_start)
        self.assertEqual(553742, hsp.hit_start)
        self.assertEqual(49, hsp.query_end)
        self.assertEqual(553781, hsp.hit_end)
        self.assertEqual(1, len(hsp))
        self.assertEqual([39], hsp.query_spans)
        self.assertEqual([39], hsp.hit_spans)
        self.assertEqual([(10, 49)], hsp.query_ranges)
        self.assertEqual([(553742, 553781)], hsp.hit_ranges)


class BlatPslxCases(BlatPslCases):

    def test_pslx_34_001(self, testf='pslx_34_001.pslx', fmt='blat-pslx'):
        """Test parsing blat output (pslx_34_001.pslx)"""
        BlatPslCases.test_psl_34_001(self, 'pslx_34_001.pslx', 'blat-pslx')

        # test first qresult
        qresult = self.qresults[0]
        # first qresult, first hit, first hsp
        hsp = qresult[0].batch_hsps[0]
        self.assertEqual('aggtaaactgccttca', str(hsp.queries[0].seq))
        self.assertEqual('aggtaaactgccttca', str(hsp.hits[0].seq))
        # first qresult, second hit, first hsp
        hsp = qresult[1].batch_hsps[0]
        self.assertEqual('atgagcttccaaggtaaactgccttcaagattc', str(hsp.queries[0].seq))
        self.assertEqual('atgagcttccaaggtaaactgccttcaagattc', str(hsp.hits[0].seq))
        # first qresult, third hit, first hsp
        hsp = qresult[2].batch_hsps[0]
        self.assertEqual('aaggcagtttaccttgg', str(hsp.queries[0].seq))
        self.assertEqual('aaggcagtttaccttgg', str(hsp.hits[0].seq))

        # test second qresult
        qresult = self.qresults[1]
        # second qresult, first hit, first hsp
        hsp = qresult[0].batch_hsps[0]
        self.assertEqual('acaaaggggctgggcgtggtggctcacacctgtaatcccaa', str(hsp.queries[0].seq))
        self.assertEqual('acaaaggggctgggcgcagtggctcacgcctgtaatcccaa', str(hsp.hits[0].seq))
        # second qresult, second hit, first hsp
        hsp = qresult[1].batch_hsps[0]
        self.assertEqual('cacaaaggggctgggcgtggtggctcacacctgtaatccca', str(hsp.queries[0].seq))
        self.assertEqual('cacaaaggggctgggcgtggtggctcacacctgtaatccca', str(hsp.hits[0].seq))
        # second qresult, third hit, first hsp
        hsp = qresult[2].batch_hsps[0]
        self.assertEqual('aaaggggctgggcgtggtggctcacacctgtaatcc', str(hsp.queries[0].seq))
        self.assertEqual('aaaggggctgggcgtggtagctcatgcctgtaatcc', str(hsp.hits[0].seq))
        # second qresult, third hit, second hsp
        hsp = qresult[2].batch_hsps[1]
        self.assertEqual('tgggattacaggtgtgagccaccacgcccagcccctt', str(hsp.queries[0].seq))
        self.assertEqual('tgggattacaggcgggagccaccacgcccagcccctt', str(hsp.hits[0].seq))
        # second qresult, fourth hit, first hsp
        hsp = qresult[3].batch_hsps[0]
        self.assertEqual('aaaaat', str(hsp.queries[0].seq))
        self.assertEqual('aaaaat', str(hsp.hits[0].seq))
        self.assertEqual('aaaggggctgggcgtggtggctcacacctgtaatccca', str(hsp.queries[1].seq))
        self.assertEqual('aaaggggctgggcgtggtggctcacgcctgtaatccca', str(hsp.hits[1].seq))
        # second qresult, fourth hit, second hsp
        hsp = qresult[3].batch_hsps[1]
        self.assertEqual('tgggattacaggtgtgagccaccacgcccagcccct', str(hsp.queries[0].seq))
        self.assertEqual('tgggattacaggcgtgagccaccacgcccagcccct', str(hsp.hits[0].seq))
        # second qresult, fifth hit, first hsp
        hsp = qresult[4].batch_hsps[0]
        self.assertEqual('caaaggggctgggcgtggtggctca', str(hsp.queries[0].seq))
        self.assertEqual('caaaggggctgggcgtagtggctga', str(hsp.hits[0].seq))
        self.assertEqual('cacctgtaatc', str(hsp.queries[1].seq))
        self.assertEqual('cacctgtaatc', str(hsp.hits[1].seq))
        # second qresult, fifth hit, second hsp
        hsp = qresult[4].batch_hsps[1]
        self.assertEqual('tgggattacaggtgtgagccaccacgcccagcccctttg', str(hsp.queries[0].seq))
        self.assertEqual('tgggattacaggtgtgagccaccacgcccagcccctttg', str(hsp.hits[0].seq))
        # second qresult, fifth hit, third hsp
        hsp = qresult[4].batch_hsps[2]
        self.assertEqual('tgggattacaggtgtgagccaccacgcccagcccctttg', str(hsp.queries[0].seq))
        self.assertEqual('tgggatgacaggggtgaggcaccacgcccagcccctttg', str(hsp.hits[0].seq))

    def test_pslx_34_002(self, testf='pslx_34_002.pslx', fmt='blat-pslx'):
        """Test parsing blat output (pslx_34_002.pslx)"""
        BlatPslCases.test_psl_34_002(self, 'pslx_34_002.pslx', 'blat-pslx')

    def test_pslx_34_003(self, testf='pslx_34_003.pslx', fmt='blat-pslx'):
        """Test parsing blat output (pslx_34_003.pslx)"""
        BlatPslCases.test_psl_34_003(self, 'pslx_34_003.pslx', 'blat-pslx')

        # test first qresult
        qresult = self.qresults[0]
        # first qresult, first hit, first hsp
        hsp = qresult[0].batch_hsps[0]
        self.assertEqual('aggtaaactgccttca', str(hsp.queries[0].seq))
        self.assertEqual('aggtaaactgccttca', str(hsp.hits[0].seq))
        # first qresult, second hit, first hsp
        hsp = qresult[1].batch_hsps[0]
        self.assertEqual('atgagcttccaaggtaaactgccttcaagattc', str(hsp.queries[0].seq))
        self.assertEqual('atgagcttccaaggtaaactgccttcaagattc', str(hsp.hits[0].seq))
        # first qresult, third hit, first hsp
        hsp = qresult[2].batch_hsps[0]
        self.assertEqual('aaggcagtttaccttgg', str(hsp.queries[0].seq))
        self.assertEqual('aaggcagtttaccttgg', str(hsp.hits[0].seq))

    def test_pslx_34_004(self, testf='pslx_34_004.pslx', fmt='blat-pslx'):
        """Test parsing blat output (pslx_34_004.pslx)"""
        BlatPslCases.test_psl_34_004(self, 'pslx_34_004.pslx', 'blat-pslx')

        # test first qresult
        qresult = self.qresults[0]
        # first qresult, first hit, first hsp
        hsp = qresult[0].batch_hsps[0]
        self.assertEqual('acaaaggggctgggcgtggtggctcacacctgtaatcccaa', str(hsp.queries[0].seq))
        self.assertEqual('acaaaggggctgggcgcagtggctcacgcctgtaatcccaa', str(hsp.hits[0].seq))
        # first qresult, second hit, first hsp
        hsp = qresult[1].batch_hsps[0]
        self.assertEqual('cacaaaggggctgggcgtggtggctcacacctgtaatccca', str(hsp.queries[0].seq))
        self.assertEqual('cacaaaggggctgggcgtggtggctcacacctgtaatccca', str(hsp.hits[0].seq))
        # first qresult, third hit, first hsp
        hsp = qresult[2].batch_hsps[0]
        self.assertEqual('aaaggggctgggcgtggtggctcacacctgtaatcc', str(hsp.queries[0].seq))
        self.assertEqual('aaaggggctgggcgtggtagctcatgcctgtaatcc', str(hsp.hits[0].seq))
        # first qresult, third hit, second hsp
        hsp = qresult[2].batch_hsps[1]
        self.assertEqual('tgggattacaggtgtgagccaccacgcccagcccctt', str(hsp.queries[0].seq))
        self.assertEqual('tgggattacaggcgggagccaccacgcccagcccctt', str(hsp.hits[0].seq))
        # first qresult, fourth hit, first hsp
        hsp = qresult[3].batch_hsps[0]
        self.assertEqual('aaaaat', str(hsp.queries[0].seq))
        self.assertEqual('aaaaat', str(hsp.hits[0].seq))
        self.assertEqual('aaaggggctgggcgtggtggctcacacctgtaatccca', str(hsp.queries[1].seq))
        self.assertEqual('aaaggggctgggcgtggtggctcacgcctgtaatccca', str(hsp.hits[1].seq))
        # first qresult, fourth hit, second hsp
        hsp = qresult[3].batch_hsps[1]
        self.assertEqual('tgggattacaggtgtgagccaccacgcccagcccct', str(hsp.queries[0].seq))
        self.assertEqual('tgggattacaggcgtgagccaccacgcccagcccct', str(hsp.hits[0].seq))
        # first qresult, fifth hit, first hsp
        hsp = qresult[4].batch_hsps[0]
        self.assertEqual('caaaggggctgggcgtggtggctca', str(hsp.queries[0].seq))
        self.assertEqual('caaaggggctgggcgtagtggctga', str(hsp.hits[0].seq))
        self.assertEqual('cacctgtaatc', str(hsp.queries[1].seq))
        self.assertEqual('cacctgtaatc', str(hsp.hits[1].seq))
        # first qresult, fifth hit, second hsp
        hsp = qresult[4].batch_hsps[1]
        self.assertEqual('tgggattacaggtgtgagccaccacgcccagcccctttg', str(hsp.queries[0].seq))
        self.assertEqual('tgggattacaggtgtgagccaccacgcccagcccctttg', str(hsp.hits[0].seq))
        # first qresult, fifth hit, third hsp
        hsp = qresult[4].batch_hsps[2]
        self.assertEqual('tgggattacaggtgtgagccaccacgcccagcccctttg', str(hsp.queries[0].seq))
        self.assertEqual('tgggatgacaggggtgaggcaccacgcccagcccctttg', str(hsp.hits[0].seq))

    def test_pslx_34_005(self, testf='pslx_34_005.pslx', fmt='blat-pslx'):
        """Test parsing blat output (pslx_34_005.pslx)"""
        BlatPslCases.test_psl_34_005(self, 'pslx_34_005.pslx', 'blat-pslx')

        # test first qresult
        qresult = self.qresults[0]
        # first qresult, first hit, first hsp
        hsp = qresult[0].batch_hsps[0]
        self.assertEqual('aggtaaactgccttca', str(hsp.queries[0].seq))
        self.assertEqual('aggtaaactgccttca', str(hsp.hits[0].seq))
        # first qresult, second hit, first hsp
        hsp = qresult[1].batch_hsps[0]
        self.assertEqual('atgagcttccaaggtaaactgccttcaagattc', str(hsp.queries[0].seq))
        self.assertEqual('atgagcttccaaggtaaactgccttcaagattc', str(hsp.hits[0].seq))
        # first qresult, third hit, first hsp
        hsp = qresult[2].batch_hsps[0]
        self.assertEqual('aaggcagtttaccttgg', str(hsp.queries[0].seq))
        self.assertEqual('aaggcagtttaccttgg', str(hsp.hits[0].seq))

        # test second qresult
        qresult = self.qresults[1]
        # second qresult, first hit, first hsp
        hsp = qresult[0].batch_hsps[0]
        self.assertEqual('acaaaggggctgggcgtggtggctcacacctgtaatcccaa', str(hsp.queries[0].seq))
        self.assertEqual('acaaaggggctgggcgcagtggctcacgcctgtaatcccaa', str(hsp.hits[0].seq))
        # second qresult, second hit, first hsp
        hsp = qresult[1].batch_hsps[0]
        self.assertEqual('cacaaaggggctgggcgtggtggctcacacctgtaatccca', str(hsp.queries[0].seq))
        self.assertEqual('cacaaaggggctgggcgtggtggctcacacctgtaatccca', str(hsp.hits[0].seq))
        # second qresult, third hit, first hsp
        hsp = qresult[2].batch_hsps[0]
        self.assertEqual('aaaggggctgggcgtggtggctcacacctgtaatcc', str(hsp.queries[0].seq))
        self.assertEqual('aaaggggctgggcgtggtagctcatgcctgtaatcc', str(hsp.hits[0].seq))
        # second qresult, third hit, second hsp
        hsp = qresult[2].batch_hsps[1]
        self.assertEqual('tgggattacaggtgtgagccaccacgcccagcccctt', str(hsp.queries[0].seq))
        self.assertEqual('tgggattacaggcgggagccaccacgcccagcccctt', str(hsp.hits[0].seq))
        # second qresult, fourth hit, first hsp
        hsp = qresult[3].batch_hsps[0]
        self.assertEqual('aaaaat', str(hsp.queries[0].seq))
        self.assertEqual('aaaaat', str(hsp.hits[0].seq))
        self.assertEqual('aaaggggctgggcgtggtggctcacacctgtaatccca', str(hsp.queries[1].seq))
        self.assertEqual('aaaggggctgggcgtggtggctcacgcctgtaatccca', str(hsp.hits[1].seq))
        # second qresult, fourth hit, second hsp
        hsp = qresult[3].batch_hsps[1]
        self.assertEqual('tgggattacaggtgtgagccaccacgcccagcccct', str(hsp.queries[0].seq))
        self.assertEqual('tgggattacaggcgtgagccaccacgcccagcccct', str(hsp.hits[0].seq))
        # second qresult, fifth hit, first hsp
        hsp = qresult[4].batch_hsps[0]
        self.assertEqual('caaaggggctgggcgtggtggctca', str(hsp.queries[0].seq))
        self.assertEqual('caaaggggctgggcgtagtggctga', str(hsp.hits[0].seq))
        self.assertEqual('cacctgtaatc', str(hsp.queries[1].seq))
        self.assertEqual('cacctgtaatc', str(hsp.hits[1].seq))
        # second qresult, fifth hit, second hsp
        hsp = qresult[4].batch_hsps[1]
        self.assertEqual('tgggattacaggtgtgagccaccacgcccagcccctttg', str(hsp.queries[0].seq))
        self.assertEqual('tgggattacaggtgtgagccaccacgcccagcccctttg', str(hsp.hits[0].seq))
        # second qresult, fifth hit, third hsp
        hsp = qresult[4].batch_hsps[2]
        self.assertEqual('tgggattacaggtgtgagccaccacgcccagcccctttg', str(hsp.queries[0].seq))
        self.assertEqual('tgggatgacaggggtgaggcaccacgcccagcccctttg', str(hsp.hits[0].seq))


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
