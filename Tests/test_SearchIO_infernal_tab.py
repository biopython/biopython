# Copyright 2024 by Samuel Prince. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.


"""Tests for SearchIO InfernalIO infernal-tab parser."""

import os
import unittest
import itertools

from Bio.SearchIO import parse

# test case files are in the Infernal directory
TEST_DIR = "Infernal"
FMT = "infernal-tab"


def get_file(filename):
    """Return the path of a test file."""
    return os.path.join(TEST_DIR, filename)


def next_result(qresults, counter):
    """Iterate over the results and counter."""
    return next(qresults), next(counter)


class CmscanCases(unittest.TestCase):
    """Test parsing cmscan output."""

    def test_cmscan_mq_mm(self):
        """Test parsing infernal-tab, cmscan, multiple queries, multiple hit, one hsp, default format (IRES_5S_U2_Yeast)"""
        tab_file = get_file("cmscan_115_IRES_5S_U2_Yeast.tbl")
        qresults = parse(tab_file, FMT)
        counter = itertools.count(start=1)

        # first qresult
        qresult, count = next_result(qresults, counter)
        self.assertEqual(len(qresult), 1)
        self.assertEqual(qresult.id, "ENA|BK006935|BK006935.2")
        self.assertEqual(qresult.accession, "-")
        hit = qresult[0]
        self.assertEqual(len(hit), 2)
        self.assertEqual(hit.id, "U2")
        self.assertEqual(hit.accession, "RF00004")
        self.assertEqual(hit.description, "U2 spliceosomal RNA")
        # first hsp
        hsp = hit[0]
        self.assertEqual(len(hsp), 1)
        self.assertEqual(hsp.evalue, 0.91)
        self.assertEqual(hsp.bitscore, 13.5)
        self.assertEqual(hsp.bias, 0.0)
        self.assertEqual(hsp.gc, 0.44)
        self.assertEqual(hsp.truncated, "no")
        self.assertEqual(hsp.model, "cm")
        self.assertEqual(hsp.pipeline_pass, 1)
        self.assertFalse(hsp.is_included)
        frag = hsp[0]
        self.assertEqual(frag.query_start, 1)
        self.assertEqual(frag.query_end, 193)
        self.assertEqual(52929, frag.hit_start)
        self.assertEqual(53083, frag.hit_end)
        self.assertEqual(frag.hit_strand, 0)
        # second hsp
        hsp = hit[1]
        self.assertEqual(len(hsp), 1)
        self.assertEqual(1.3, hsp.evalue)
        self.assertEqual(12.8, hsp.bitscore)
        self.assertEqual(5.3, hsp.bias)
        self.assertEqual(hsp.gc, 0.33)
        self.assertEqual(hsp.truncated, "no")
        self.assertEqual(hsp.pipeline_pass, 1)
        self.assertFalse(hsp.is_included)
        frag = hsp[0]
        self.assertEqual(frag.query_start, 1)
        self.assertEqual(frag.query_end, 193)
        self.assertEqual(frag.hit_start, 196389)
        self.assertEqual(frag.hit_end, 196571)
        self.assertEqual(frag.hit_strand, -1)

        # second qresult
        qresult, count = next_result(qresults, counter)
        self.assertEqual(len(qresult), 1)
        self.assertEqual("ENA|BK006936|BK006936.2", qresult.id)
        self.assertEqual(qresult.accession, "-")
        hit = qresult[0]
        self.assertEqual(len(hit), 1)
        self.assertEqual(hit.id, "U2")
        self.assertEqual(hit.accession, "RF00004")
        self.assertEqual(hit.description, "U2 spliceosomal RNA")
        hsp = hit[0]
        self.assertEqual(len(hsp), 1)
        self.assertEqual(hsp.evalue, 1.2e-20)
        self.assertEqual(hsp.bitscore, 98.7)
        self.assertEqual(hsp.bias, 0.1)
        self.assertEqual(hsp.gc, 0.33)
        self.assertEqual(hsp.truncated, "no")
        self.assertEqual(hsp.model, "cm")
        self.assertEqual(hsp.pipeline_pass, 1)
        self.assertTrue(hsp.is_included)
        frag = hsp[0]
        self.assertEqual(frag.query_start, 1)
        self.assertEqual(frag.query_end, 193)
        self.assertEqual(frag.hit_start, 681747)
        self.assertEqual(frag.hit_end, 681858)
        self.assertEqual(frag.hit_strand, -1)

        # third qresult
        qresult, count = next_result(qresults, counter)
        self.assertEqual(len(qresult), 2)
        self.assertEqual(qresult.id, "ENA|BK006937|BK006937.2")
        self.assertEqual(qresult.accession, "-")
        # first hit
        hit = qresult[0]
        self.assertEqual(len(hit), 1)
        self.assertEqual(hit.id, "5S_rRNA")
        self.assertEqual(hit.accession, "RF00001")
        self.assertEqual(hit.description, "5S ribosomal RNA")
        hsp = hit[0]
        self.assertEqual(len(hsp), 1)
        self.assertEqual(hsp.evalue, 2.4)
        self.assertEqual(hsp.bitscore, 14.1)
        self.assertEqual(hsp.bias, 0.3)
        self.assertEqual(hsp.gc, 0.41)
        self.assertEqual(hsp.truncated, "no")
        self.assertEqual(hsp.model, "cm")
        self.assertEqual(hsp.pipeline_pass, 1)
        self.assertFalse(hsp.is_included)
        frag = hsp[0]
        self.assertEqual(frag.query_start, 1)
        self.assertEqual(frag.query_end, 119)
        self.assertEqual(frag.hit_start, 644)
        self.assertEqual(frag.hit_end, 761)
        self.assertEqual(frag.hit_strand, -1)
        # second hit
        hit = qresult[1]
        self.assertEqual(len(hit), 1)
        self.assertEqual(hit.id, "U2")
        self.assertEqual(hit.accession, "RF00004")
        self.assertEqual(hit.description, "U2 spliceosomal RNA")
        hsp = hit[0]
        self.assertEqual(len(hsp), 1)
        self.assertEqual(hsp.evalue, 4.7)
        self.assertEqual(hsp.bitscore, 11.1)
        self.assertEqual(hsp.bias, 0.1)
        self.assertEqual(hsp.gc, 0.32)
        self.assertEqual(hsp.truncated, "no")
        self.assertEqual(hsp.pipeline_pass, 1)
        self.assertFalse(hsp.is_included)
        frag = hsp[0]
        self.assertEqual(frag.query_start, 1)
        self.assertEqual(frag.query_end, 193)
        self.assertEqual(frag.hit_start, 229885)
        self.assertEqual(frag.hit_end, 229986)
        self.assertEqual(frag.hit_strand, -1)

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(count, 3)

    def test_cmscan_mq_mm_fmt2(self):
        """Test parsing infernal-tab, cmscan, multiple queries, multiple hit, one hsp, fmt 2 (IRES_5S_U2_Yeast_fmt_2)"""
        tab_file = get_file("cmscan_115_IRES_5S_U2_Yeast_fmt_2.tbl")
        qresults = parse(tab_file, FMT, fmt=2)
        counter = itertools.count(start=1)

        # first qresult
        qresult, count = next_result(qresults, counter)
        self.assertEqual(len(qresult), 1)
        self.assertEqual(qresult.id, "ENA|BK006936|BK006936.2")
        self.assertEqual(qresult.accession, "-")
        self.assertEqual(qresult.clan, "-")
        self.assertEqual(qresult.seq_len, 813184)
        hit = qresult[0]
        self.assertEqual(len(hit), 1)
        self.assertEqual(hit.id, "U2")
        self.assertEqual(hit.accession, "RF00004")
        self.assertEqual(hit.description, "U2 spliceosomal RNA")
        self.assertEqual(hit.seq_len, 193)
        # first hsp
        hsp = hit[0]
        self.assertEqual(len(hsp), 1)
        self.assertEqual(hsp.evalue, 1.2e-20)
        self.assertEqual(hsp.bitscore, 98.7)
        self.assertEqual(hsp.bias, 0.1)
        self.assertEqual(hsp.gc, 0.33)
        self.assertEqual(hsp.truncated, "no")
        self.assertEqual(hsp.model, "cm")
        self.assertEqual(hsp.pipeline_pass, 1)
        self.assertTrue(hsp.is_included)
        self.assertEqual("*", hsp.olp)
        self.assertEqual("-", hsp.anyidx)
        self.assertEqual("-", hsp.afrct1)
        self.assertEqual("-", hsp.afrct2)
        self.assertEqual("-", hsp.winidx)
        self.assertEqual("-", hsp.wfrct1)
        self.assertEqual("-", hsp.wfrct2)
        frag = hsp[0]
        self.assertEqual(frag.query_start, 1)
        self.assertEqual(frag.query_end, 193)
        self.assertEqual(frag.hit_start, 681747)
        self.assertEqual(frag.hit_end, 681858)
        self.assertEqual(frag.hit_strand, -1)

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(count, 1)

    def test_cmscan_mq_mm_fmt3(self):
        """Test parsing infernal-tab, cmscan, multiple queries, multiple hit, one hsp, fmt 3 (IRES_5S_U2_Yeast_fmt_3)"""
        tab_file = get_file("cmscan_115_IRES_5S_U2_Yeast_fmt_3.tbl")
        qresults = parse(tab_file, FMT, fmt=3)
        counter = itertools.count(start=1)

        # first qresult
        qresult, count = next_result(qresults, counter)
        self.assertEqual(len(qresult), 1)
        self.assertEqual("ENA|BK006936|BK006936.2", qresult.id)
        self.assertEqual(qresult.accession, "-")
        self.assertEqual(813184, qresult.seq_len)
        hit = qresult[0]
        self.assertEqual(len(hit), 1)
        self.assertEqual(hit.id, "U2")
        self.assertEqual(hit.accession, "RF00004")
        self.assertEqual(hit.description, "U2 spliceosomal RNA")
        self.assertEqual(193, hit.seq_len)
        # first hsp
        hsp = hit[0]
        self.assertEqual(len(hsp), 1)
        self.assertEqual(hsp.evalue, 1.2e-20)
        self.assertEqual(hsp.bitscore, 98.7)
        self.assertEqual(hsp.bias, 0.1)
        self.assertEqual(hsp.gc, 0.33)
        self.assertEqual(hsp.truncated, "no")
        self.assertEqual(hsp.model, "cm")
        self.assertEqual(hsp.pipeline_pass, 1)
        self.assertTrue(hsp.is_included)
        frag = hsp[0]
        self.assertEqual(frag.query_start, 1)
        self.assertEqual(frag.query_end, 193)
        self.assertEqual(frag.hit_start, 681747)
        self.assertEqual(frag.hit_end, 681858)
        self.assertEqual(frag.hit_strand, -1)

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(count, 1)


class CmsearchCases(unittest.TestCase):
    """Test parsing cmsearch output."""

    def test_1q_0m(self):
        """Test parsing infernal-tab, cmsearch, one query, no hits (IRES_Yeast)"""
        tab_file = get_file("cmsearch_114_IRES_Yeast.tbl")
        qresults = parse(tab_file, FMT)

        self.assertRaises(StopIteration, next, qresults)

    def test_cmsearch_1q_1m(self):
        """Test parsing infernal-tab, cmsearch, one queries, one hit, one hsp (U2_Yeast)"""
        tab_file = get_file("cmsearch_114_U2_Yeast.tbl")
        qresults = parse(tab_file, FMT)
        counter = itertools.count(start=1)

        qresult, count = next_result(qresults, counter)
        self.assertEqual(len(qresult), 1)
        self.assertEqual(qresult.id, "U2")
        self.assertEqual(qresult.accession, "RF00004")
        hit = qresult[0]
        self.assertEqual(len(hit), 1)
        self.assertEqual(hit.id, "ENA|BK006936|BK006936.2")
        self.assertEqual(hit.accession, "-")
        self.assertEqual(
            hit.description,
            "TPA_inf: Saccharomyces cerevisiae S288C chromosome II, complete sequence.",
        )

        hsp = hit[0]
        self.assertEqual(len(hsp), 1)
        self.assertEqual(hsp.evalue, 5.9e-20)
        self.assertEqual(hsp.bitscore, 98.7)
        self.assertEqual(hsp.bias, 0.1)
        self.assertTrue(hsp.is_included)
        self.assertEqual(hsp.gc, 0.33)
        self.assertEqual(hsp.truncated, "no")
        self.assertEqual(hsp.model, "cm")
        self.assertEqual(hsp.pipeline_pass, 1)
        frag = hsp[0]
        self.assertEqual(frag.query_start, 1)
        self.assertEqual(frag.query_end, 193)
        self.assertEqual(frag.hit_start, 681747)
        self.assertEqual(frag.hit_end, 681858)
        self.assertEqual(frag.hit_strand, -1)

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(count, 1)

    def test_cmsearch_1q_mm(self):
        """Test parsing infernal-tab, cmsearch, one queries, multiple hit, one hsp (5S_Yeast)"""
        tab_file = get_file("cmsearch_114_5S_Yeast.tbl")
        qresults = parse(tab_file, FMT)
        counter = itertools.count(start=1)

        qresult, count = next_result(qresults, counter)
        self.assertEqual(len(qresult), 1)
        self.assertEqual(qresult.id, "5S_rRNA")
        self.assertEqual(qresult.accession, "RF00001")
        # first hit
        hit = qresult[0]
        self.assertEqual(6, len(hit))
        self.assertEqual(hit.id, "ENA|BK006945|BK006945.2")
        self.assertEqual(hit.accession, "-")
        self.assertEqual(
            hit.description,
            "TPA_inf: Saccharomyces cerevisiae S288C chromosome XII, complete sequence.",
        )
        hsp = hit[0]
        self.assertEqual(len(hsp), 1)
        self.assertEqual(hsp.evalue, 1.6e-18)
        self.assertEqual(hsp.bitscore, 88.8)
        self.assertEqual(hsp.bias, 0.0)
        self.assertEqual(hsp.gc, 0.52)
        self.assertEqual(hsp.truncated, "no")
        self.assertEqual(hsp.model, "cm")
        self.assertEqual(hsp.pipeline_pass, 1)
        self.assertTrue(hsp.is_included)
        frag = hsp[0]
        self.assertEqual(frag.query_start, 1)
        self.assertEqual(frag.query_end, 119)
        self.assertEqual(frag.hit_start, 459676)
        self.assertEqual(frag.hit_end, 459796)
        self.assertEqual(frag.hit_strand, 0)
        # last hit
        hsp = hit[-1]
        self.assertEqual(len(hsp), 1)
        self.assertEqual(hsp.evalue, 4.4e-17)
        self.assertEqual(hsp.bitscore, 83.2)
        self.assertEqual(hsp.bias, 0.0)
        self.assertEqual(hsp.gc, 0.53)
        self.assertEqual(hsp.truncated, "no")
        self.assertEqual(hsp.pipeline_pass, 1)
        self.assertTrue(hsp.is_included)
        frag = hsp[0]
        self.assertEqual(frag.query_start, 1)
        self.assertEqual(frag.query_end, 119)
        self.assertEqual(frag.hit_start, 485697)
        self.assertEqual(frag.hit_end, 485817)
        self.assertEqual(frag.hit_strand, 0)

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(count, 1)

    def test_cmsearch_1q_mm_shuf(self):
        """Test parsing infernal-tab, cmsearch, one queries, multiple non-consecutive hits, one hsp (U2_Yeast_full_shuffled)"""
        tab_file = get_file("cmsearch_114_U2_Yeast_full_shuffled.tbl")
        qresults = parse(tab_file, FMT)
        counter = itertools.count(start=1)

        qresult, count = next_result(qresults, counter)
        self.assertEqual(2, len(qresult))
        self.assertEqual(qresult.id, "U2")
        self.assertEqual(qresult.accession, "RF00004")
        # first hit
        # first hit (3 hsps at rank 1,3 and 4)
        hit = qresult[0]
        self.assertEqual(3, len(hit))
        self.assertEqual(hit.id, "ENA|BK006936|BK006936.2")
        self.assertEqual(hit.description, "-")
        self.assertEqual(hit.query_id, "U2")
        # first hsp (rank 1)
        hsp = hit[0]
        self.assertEqual(hsp.query_start, 1)
        self.assertEqual(hsp.query_end, 193)
        self.assertEqual(hsp.hit_start, 681747)
        self.assertEqual(hsp.hit_end, 681858)
        # second hsp (rank 3)
        hsp = hit[1]
        self.assertEqual(hsp.query_start, 1)
        self.assertEqual(hsp.query_end, 193)
        self.assertEqual(hsp.hit_start, 1370418)
        self.assertEqual(hsp.hit_end, 1370563)
        # last hsp (rank 4)
        hsp = hit[2]
        self.assertEqual(hsp.query_start, 1)
        self.assertEqual(hsp.query_end, 193)
        self.assertEqual(hsp.hit_start, 1079243)
        self.assertEqual(hsp.hit_end, 1079392)
        # second hit
        hit = qresult[1]
        self.assertEqual(3, len(hit))
        self.assertEqual(hit.id, "ENA|BK006948|BK006948.2")
        self.assertEqual(hit.description, "-")
        self.assertEqual(hit.query_id, "U2")
        # first hsp (rank 2)
        hsp = hit[0]
        self.assertEqual(hsp.query_start, 1)
        self.assertEqual(hsp.query_end, 193)
        self.assertEqual(hsp.hit_start, 737324)
        self.assertEqual(hsp.hit_end, 737498)
        # second hsp (rank 5)
        hsp = hit[1]
        self.assertEqual(hsp.query_start, 1)
        self.assertEqual(hsp.query_end, 193)
        self.assertEqual(hsp.hit_start, 425490)
        self.assertEqual(hsp.hit_end, 425693)
        # last hsp (rank 6)
        hsp = hit[2]
        self.assertEqual(hsp.query_start, 1)
        self.assertEqual(hsp.query_end, 193)
        self.assertEqual(hsp.hit_start, 1073786)
        self.assertEqual(hsp.hit_end, 1073950)

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(count, 1)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
