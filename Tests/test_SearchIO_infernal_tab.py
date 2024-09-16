# Copyright 2024 by Samuel Prince. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Tests for SearchIO InfernalIO infernal-tab parser."""

import os
import unittest

from Bio.SearchIO import parse

# test case files are in the Blast directory
TEST_DIR = "Infernal"
FMT = "infernal-tab"


def get_file(filename):
    """Return the path of a test file."""
    return os.path.join(TEST_DIR, filename)


class CmscanCases(unittest.TestCase):
    """Test parsing cmscan output."""

    def test_cmscan_mq_mm(self):
        """Test parsing infernal-tab, cmscan, multiple queries, multiple match, one hsp, default format"""
        tab_file = get_file("IRES_5S_U2_Yeast-cmscan.tbl")
        qresults = parse(tab_file, FMT)
        counter = 0

        # first qresult
        qresult = next(qresults)
        counter += 1
        self.assertEqual(1, len(qresult))
        self.assertEqual("ENA|BK006935|BK006935.2", qresult.id)
        self.assertEqual("-", qresult.accession)
        self.assertEqual("cm", qresult.model)
        hit = qresult[0]
        self.assertEqual(2, len(hit))
        self.assertEqual("U2", hit.id)
        self.assertEqual("RF00004", hit.accession)
        self.assertEqual("U2 spliceosomal RNA", hit.description)
        # first hsp
        hsp = hit[0]
        self.assertEqual(1, len(hsp))
        self.assertEqual(0.91, hsp.evalue)
        self.assertEqual(13.5, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(0.44, hsp.gc)
        self.assertEqual("no", hsp.truncated)
        self.assertEqual(1, hsp.pipeline_pass)
        self.assertEqual(False, hsp.is_included)
        frag = hsp[0]
        self.assertEqual(1, frag.query_start)
        self.assertEqual(193, frag.query_end)
        self.assertEqual(52929, frag.hit_start)
        self.assertEqual(53083, frag.hit_end)
        self.assertEqual(0, frag.hit_strand)
        # second hsp
        hsp = hit[1]
        self.assertEqual(1, len(hsp))
        self.assertEqual(1.3, hsp.evalue)
        self.assertEqual(12.8, hsp.bitscore)
        self.assertEqual(5.3, hsp.bias)
        self.assertEqual(0.33, hsp.gc)
        self.assertEqual("no", hsp.truncated)
        self.assertEqual(1, hsp.pipeline_pass)
        self.assertEqual(False, hsp.is_included)
        frag = hsp[0]
        self.assertEqual(1, frag.query_start)
        self.assertEqual(193, frag.query_end)
        self.assertEqual(196389, frag.hit_start)
        self.assertEqual(196571, frag.hit_end)
        self.assertEqual(-1, frag.hit_strand)

        # second qresult
        qresult = next(qresults)
        counter += 1
        self.assertEqual(1, len(qresult))
        self.assertEqual("ENA|BK006936|BK006936.2", qresult.id)
        self.assertEqual("-", qresult.accession)
        self.assertEqual("cm", qresult.model)
        hit = qresult[0]
        self.assertEqual(1, len(hit))
        self.assertEqual("U2", hit.id)
        self.assertEqual("RF00004", hit.accession)
        self.assertEqual("U2 spliceosomal RNA", hit.description)
        hsp = hit[0]
        self.assertEqual(1, len(hsp))
        self.assertEqual(1.2e-20, hsp.evalue)
        self.assertEqual(98.7, hsp.bitscore)
        self.assertEqual(0.1, hsp.bias)
        self.assertEqual(0.33, hsp.gc)
        self.assertEqual("no", hsp.truncated)
        self.assertEqual(1, hsp.pipeline_pass)
        self.assertEqual(True, hsp.is_included)
        frag = hsp[0]
        self.assertEqual(1, frag.query_start)
        self.assertEqual(193, frag.query_end)
        self.assertEqual(681747, frag.hit_start)
        self.assertEqual(681858, frag.hit_end)
        self.assertEqual(-1, frag.hit_strand)

        # third qresult
        qresult = next(qresults)
        counter += 1
        self.assertEqual(2, len(qresult))
        self.assertEqual("ENA|BK006937|BK006937.2", qresult.id)
        self.assertEqual("-", qresult.accession)
        self.assertEqual("cm", qresult.model)
        # first hit
        hit = qresult[0]
        self.assertEqual(1, len(hit))
        self.assertEqual("5S_rRNA", hit.id)
        self.assertEqual("RF00001", hit.accession)
        self.assertEqual("5S ribosomal RNA", hit.description)
        hsp = hit[0]
        self.assertEqual(1, len(hsp))
        self.assertEqual(2.4, hsp.evalue)
        self.assertEqual(14.1, hsp.bitscore)
        self.assertEqual(0.3, hsp.bias)
        self.assertEqual(0.41, hsp.gc)
        self.assertEqual("no", hsp.truncated)
        self.assertEqual(1, hsp.pipeline_pass)
        self.assertEqual(False, hsp.is_included)
        frag = hsp[0]
        self.assertEqual(1, frag.query_start)
        self.assertEqual(119, frag.query_end)
        self.assertEqual(644, frag.hit_start)
        self.assertEqual(761, frag.hit_end)
        self.assertEqual(-1, frag.hit_strand)
        # second hit
        hit = qresult[1]
        self.assertEqual(1, len(hit))
        self.assertEqual("U2", hit.id)
        self.assertEqual("RF00004", hit.accession)
        self.assertEqual("U2 spliceosomal RNA", hit.description)
        hsp = hit[0]
        self.assertEqual(1, len(hsp))
        self.assertEqual(4.7, hsp.evalue)
        self.assertEqual(11.1, hsp.bitscore)
        self.assertEqual(0.1, hsp.bias)
        self.assertEqual(0.32, hsp.gc)
        self.assertEqual("no", hsp.truncated)
        self.assertEqual(1, hsp.pipeline_pass)
        self.assertEqual(False, hsp.is_included)
        frag = hsp[0]
        self.assertEqual(1, frag.query_start)
        self.assertEqual(193, frag.query_end)
        self.assertEqual(229885, frag.hit_start)
        self.assertEqual(229986, frag.hit_end)
        self.assertEqual(-1, frag.hit_strand)


    def test_cmscan_mq_mm_fmt2(self):
        """Test parsing infernal-tab, cmscan, multiple queries, multiple match, one hsp, fmt 2"""
        tab_file = get_file("IRES_5S_U2_Yeast-cmscan-fmt_2.tbl")
        qresults = parse(tab_file, FMT)
        counter = 0

        # first qresult
        qresult = next(qresults)
        counter += 1
        self.assertEqual(1, len(qresult))
        self.assertEqual("ENA|BK006936|BK006936.2", qresult.id)
        self.assertEqual("-", qresult.accession)
        self.assertEqual("cm", qresult.model)
        self.assertEqual("-", qresult.clan)
        self.assertEqual(813184, qresult.seq_len)
        hit = qresult[0]
        self.assertEqual(1, len(hit))
        self.assertEqual("U2", hit.id)
        self.assertEqual("RF00004", hit.accession)
        self.assertEqual("U2 spliceosomal RNA", hit.description)
        self.assertEqual(193, hit.seq_len)
        # first hsp
        hsp = hit[0]
        self.assertEqual(1, len(hsp))
        self.assertEqual(1.2e-20, hsp.evalue)
        self.assertEqual(98.7, hsp.bitscore)
        self.assertEqual(0.1, hsp.bias)
        self.assertEqual(0.33, hsp.gc)
        self.assertEqual("no", hsp.truncated)
        self.assertEqual(1, hsp.pipeline_pass)
        self.assertEqual(True, hsp.is_included)
        self.assertEqual("*", hsp.olp)
        self.assertEqual("-", hsp.anyidx)
        self.assertEqual("-", hsp.afrct1)
        self.assertEqual("-", hsp.afrct2)
        self.assertEqual("-", hsp.winidx)
        self.assertEqual("-", hsp.wfrct1)
        self.assertEqual("-", hsp.wfrct2)
        frag = hsp[0]
        self.assertEqual(1, frag.query_start)
        self.assertEqual(193, frag.query_end)
        self.assertEqual(681747, frag.hit_start)
        self.assertEqual(681858, frag.hit_end)
        self.assertEqual(-1, frag.hit_strand)


    def test_cmscan_mq_mm_fmt3(self):
        """Test parsing infernal-tab, cmscan, multiple queries, multiple match, one hsp, fmt 3"""
        tab_file = get_file("IRES_5S_U2_Yeast-cmscan-fmt_3.tbl")
        qresults = parse(tab_file, FMT)
        counter = 0

        # first qresult
        qresult = next(qresults)
        counter += 1
        self.assertEqual(1, len(qresult))
        self.assertEqual("ENA|BK006936|BK006936.2", qresult.id)
        self.assertEqual("-", qresult.accession)
        self.assertEqual("cm", qresult.model)
        self.assertEqual(813184, qresult.seq_len)
        hit = qresult[0]
        self.assertEqual(1, len(hit))
        self.assertEqual("U2", hit.id)
        self.assertEqual("RF00004", hit.accession)
        self.assertEqual("U2 spliceosomal RNA", hit.description)
        self.assertEqual(193, hit.seq_len)
        # first hsp
        hsp = hit[0]
        self.assertEqual(1, len(hsp))
        self.assertEqual(1.2e-20, hsp.evalue)
        self.assertEqual(98.7, hsp.bitscore)
        self.assertEqual(0.1, hsp.bias)
        self.assertEqual(0.33, hsp.gc)
        self.assertEqual("no", hsp.truncated)
        self.assertEqual(1, hsp.pipeline_pass)
        self.assertEqual(True, hsp.is_included)
        frag = hsp[0]
        self.assertEqual(1, frag.query_start)
        self.assertEqual(193, frag.query_end)
        self.assertEqual(681747, frag.hit_start)
        self.assertEqual(681858, frag.hit_end)
        self.assertEqual(-1, frag.hit_strand)


class CmsearchCases(unittest.TestCase):
    """Test parsing cmsearch output."""


    def test_1q_0m(self):
        """Test parsing infernal-tab, cmsearch, single query, no hits"""
        tab_file = get_file("IRES_Yeast.tbl")
        qresults = parse(tab_file, FMT)

        self.assertRaises(StopIteration, next, qresults)


    def test_cmsearch_1q_1m(self):
        """Test parsing infernal-tab, cmsearch, one queries, one match, one hsp"""
        tab_file = get_file("U2_Yeast-threshold.tbl")
        qresults = parse(tab_file, FMT)
        counter = 0

        qresult = next(qresults)

        counter += 1
        self.assertEqual(1, len(qresult))
        self.assertEqual("U2", qresult.id)
        self.assertEqual("RF00004", qresult.accession)
        self.assertEqual("cm", qresult.model)
        hit = qresult[0]
        self.assertEqual(1, len(hit))
        self.assertEqual("ENA|BK006936|BK006936.2", hit.id)
        self.assertEqual("-", hit.accession)
        self.assertEqual("TPA_inf: Saccharomyces cerevisiae S288C chromosome II, complete sequence.", hit.description)
        hsp = hit[0]
        self.assertEqual(1, len(hsp))
        self.assertEqual(5.9e-20, hsp.evalue)
        self.assertEqual(98.7, hsp.bitscore)
        self.assertEqual(0.1, hsp.bias)
        self.assertEqual(0.33, hsp.gc)
        self.assertEqual("no", hsp.truncated)
        self.assertEqual(1, hsp.pipeline_pass)
        self.assertEqual(True, hsp.is_included)
        frag = hsp[0]
        self.assertEqual(1, frag.query_start)
        self.assertEqual(193, frag.query_end)
        self.assertEqual(681747, frag.hit_start)
        self.assertEqual(681858, frag.hit_end)
        self.assertEqual(-1, frag.hit_strand)

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(1, counter)


    def test_cmsearch_1q_mm(self):
        """Test parsing infernal-tab, cmsearch, one queries, multiple match, one hsp"""
        tab_file = get_file("5S_Yeast.tbl")
        qresults = parse(tab_file, FMT)
        counter = 0

        qresult = next(qresults)
        counter += 1
        self.assertEqual(1, len(qresult))
        self.assertEqual("5S_rRNA", qresult.id)
        self.assertEqual("RF00001", qresult.accession)
        self.assertEqual("cm", qresult.model)
        # first hit
        hit = qresult[0]
        self.assertEqual(6, len(hit))
        self.assertEqual("ENA|BK006945|BK006945.2", hit.id)
        self.assertEqual("-", hit.accession)
        self.assertEqual("TPA_inf: Saccharomyces cerevisiae S288C chromosome XII, complete sequence.", hit.description)
        hsp = hit[0]
        self.assertEqual(1, len(hsp))
        self.assertEqual(1.6e-18, hsp.evalue)
        self.assertEqual(88.8, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(0.52, hsp.gc)
        self.assertEqual("no", hsp.truncated)
        self.assertEqual(1, hsp.pipeline_pass)
        self.assertEqual(True, hsp.is_included)
        frag = hsp[0]
        self.assertEqual(1, frag.query_start)
        self.assertEqual(119, frag.query_end)
        self.assertEqual(459676, frag.hit_start)
        self.assertEqual(459796, frag.hit_end)
        self.assertEqual(0, frag.hit_strand)
        # last hit
        hsp = hit[-1]
        self.assertEqual(1, len(hsp))
        self.assertEqual(4.4e-17, hsp.evalue)
        self.assertEqual(83.2, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(0.53, hsp.gc)
        self.assertEqual("no", hsp.truncated)
        self.assertEqual(1, hsp.pipeline_pass)
        self.assertEqual(True, hsp.is_included)
        frag = hsp[0]
        self.assertEqual(1, frag.query_start)
        self.assertEqual(119, frag.query_end)
        self.assertEqual(485697, frag.hit_start)
        self.assertEqual(485817, frag.hit_end)
        self.assertEqual(0, frag.hit_strand)

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(1, counter)



if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)