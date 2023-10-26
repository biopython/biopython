# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO HmmerIO hmmer3-tab parser."""


import os
import unittest

from Bio.SearchIO import parse


# test case files are in the Blast directory
TEST_DIR = "Hmmer"
FMT = "hmmer3-tab"


def get_file(filename):
    """Return the path of a test file."""
    return os.path.join(TEST_DIR, filename)


class HmmscanCases(unittest.TestCase):
    """Test parsing hmmscan output."""

    def test_31b1_hmmscan_001(self):
        """Test parsing hmmer3-tab, hmmscan 3.1b1, multiple queries (tab_31b1_hmmscan_001)."""
        tab_file = get_file("tab_31b1_hmmscan_001.out")
        qresults = list(parse(tab_file, FMT))
        self.assertEqual(4, len(qresults))

        # first qresult, first hit, first hsp
        qresult = qresults[0]
        self.assertEqual(1, len(qresult))
        self.assertEqual("gi|4885477|ref|NP_005359.1|", qresult.id)
        self.assertEqual("-", qresult.accession)
        hit = qresult[0]
        self.assertEqual(1, len(hit))
        self.assertEqual("Globin", hit.id)
        self.assertEqual("PF00042.17", hit.accession)
        self.assertEqual(1e-22, hit.evalue)
        self.assertEqual(80.5, hit.bitscore)
        self.assertEqual(0.3, hit.bias)
        self.assertEqual(1.3, hit.domain_exp_num)
        self.assertEqual(1, hit.region_num)
        self.assertEqual(0, hit.cluster_num)
        self.assertEqual(0, hit.overlap_num)
        self.assertEqual(1, hit.env_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, hit.domain_reported_num)
        self.assertEqual(1, hit.domain_included_num)
        self.assertEqual("Globin", hit.description)
        hsp = hit.hsps[0]
        self.assertEqual(1.6e-22, hsp.evalue)
        self.assertEqual(79.8, hsp.bitscore)
        self.assertEqual(0.3, hsp.bias)

        # last qresult, last hit, last hsp
        qresult = qresults[-1]
        self.assertEqual(5, len(qresult))
        self.assertEqual("gi|125490392|ref|NP_038661.2|", qresult.id)
        self.assertEqual("-", qresult.accession)
        hit = qresult[-1]
        self.assertEqual(1, len(hit))
        self.assertEqual("DUF521", hit.id)
        self.assertEqual("PF04412.8", hit.accession)
        self.assertEqual(0.15, hit.evalue)
        self.assertEqual(10.5, hit.bitscore)
        self.assertEqual(0.1, hit.bias)
        self.assertEqual(1.4, hit.domain_exp_num)
        self.assertEqual(1, hit.region_num)
        self.assertEqual(0, hit.cluster_num)
        self.assertEqual(0, hit.overlap_num)
        self.assertEqual(1, hit.env_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, hit.domain_reported_num)
        self.assertEqual(0, hit.domain_included_num)
        self.assertEqual("Protein of unknown function (DUF521)", hit.description)
        hsp = hit.hsps[0]
        self.assertEqual(0.28, hsp.evalue)
        self.assertEqual(9.6, hsp.bitscore)
        self.assertEqual(0.1, hsp.bias)

    def test_30_hmmscan_001(self):
        """Test parsing hmmer3-tab, hmmscan 3.0, multiple queries (tab_30_hmmscan_001)."""
        tab_file = get_file("tab_30_hmmscan_001.out")
        qresults = parse(tab_file, FMT)
        counter = 0

        # first qresult
        qresult = next(qresults)
        counter += 1
        self.assertEqual(1, len(qresult))
        self.assertEqual("gi|4885477|ref|NP_005359.1|", qresult.id)
        self.assertEqual("-", qresult.accession)
        hit = qresult[0]
        self.assertEqual(1, len(hit))
        self.assertEqual("Globin", hit.id)
        self.assertEqual("PF00042.17", hit.accession)
        self.assertEqual(6e-21, hit.evalue)
        self.assertEqual(74.6, hit.bitscore)
        self.assertEqual(0.3, hit.bias)
        self.assertEqual(1.3, hit.domain_exp_num)
        self.assertEqual(1, hit.region_num)
        self.assertEqual(0, hit.cluster_num)
        self.assertEqual(0, hit.overlap_num)
        self.assertEqual(1, hit.env_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, hit.domain_reported_num)
        self.assertEqual(1, hit.domain_included_num)
        self.assertEqual("Globin", hit.description)
        hsp = hit.hsps[0]
        self.assertEqual(9.2e-21, hsp.evalue)
        self.assertEqual(74.0, hsp.bitscore)
        self.assertEqual(0.2, hsp.bias)

        # second qresult
        qresult = next(qresults)
        counter += 1
        self.assertEqual(2, len(qresult))
        self.assertEqual("gi|126362951:116-221", qresult.id)
        self.assertEqual("-", qresult.accession)
        hit = qresult[0]
        self.assertEqual(1, len(hit))
        self.assertEqual("Ig_3", hit.id)
        self.assertEqual("PF13927.1", hit.accession)
        self.assertEqual(1.4e-09, hit.evalue)
        self.assertEqual(38.2, hit.bitscore)
        self.assertEqual(0.4, hit.bias)
        self.assertEqual(1.3, hit.domain_exp_num)
        self.assertEqual(1, hit.region_num)
        self.assertEqual(0, hit.cluster_num)
        self.assertEqual(0, hit.overlap_num)
        self.assertEqual(1, hit.env_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, hit.domain_reported_num)
        self.assertEqual(1, hit.domain_included_num)
        self.assertEqual("Immunoglobulin domain", hit.description)
        hsp = hit.hsps[0]
        self.assertEqual(2.1e-09, hsp.evalue)
        self.assertEqual(37.6, hsp.bitscore)
        self.assertEqual(0.3, hsp.bias)
        hit = qresult[1]
        self.assertEqual(1, len(hit))
        self.assertEqual("Ig_2", hit.id)
        self.assertEqual("PF13895.1", hit.accession)
        self.assertEqual(3.5e-05, hit.evalue)
        self.assertEqual(23.7, hit.bitscore)
        self.assertEqual(0.1, hit.bias)
        self.assertEqual(1.1, hit.domain_exp_num)
        self.assertEqual(1, hit.region_num)
        self.assertEqual(0, hit.cluster_num)
        self.assertEqual(0, hit.overlap_num)
        self.assertEqual(1, hit.env_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, hit.domain_reported_num)
        self.assertEqual(1, hit.domain_included_num)
        self.assertEqual("Immunoglobulin domain", hit.description)
        hsp = hit.hsps[0]
        self.assertEqual(4.3e-05, hsp.evalue)
        self.assertEqual(23.4, hsp.bitscore)
        self.assertEqual(0.1, hsp.bias)

        # third qresult
        qresult = next(qresults)
        counter += 1
        self.assertEqual(2, len(qresult))
        self.assertEqual("gi|22748937|ref|NP_065801.1|", qresult.id)
        self.assertEqual("-", qresult.accession)
        hit = qresult[0]
        self.assertEqual(1, len(hit))
        self.assertEqual("Xpo1", hit.id)
        self.assertEqual("PF08389.7", hit.accession)
        self.assertEqual(7.8e-34, hit.evalue)
        self.assertEqual(116.6, hit.bitscore)
        self.assertEqual(7.8, hit.bias)
        self.assertEqual(2.8, hit.domain_exp_num)
        self.assertEqual(2, hit.region_num)
        self.assertEqual(0, hit.cluster_num)
        self.assertEqual(0, hit.overlap_num)
        self.assertEqual(2, hit.env_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, hit.domain_reported_num)
        self.assertEqual(1, hit.domain_included_num)
        self.assertEqual("Exportin 1-like protein", hit.description)
        hsp = hit.hsps[0]
        self.assertEqual(1.1e-33, hsp.evalue)
        self.assertEqual(116.1, hsp.bitscore)
        self.assertEqual(3.4, hsp.bias)
        hit = qresult[1]
        self.assertEqual(1, len(hit))
        self.assertEqual("IBN_N", hit.id)
        self.assertEqual("PF03810.14", hit.accession)
        self.assertEqual(0.0039, hit.evalue)
        self.assertEqual(16.9, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual(2.7, hit.domain_exp_num)
        self.assertEqual(2, hit.region_num)
        self.assertEqual(0, hit.cluster_num)
        self.assertEqual(0, hit.overlap_num)
        self.assertEqual(2, hit.env_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, hit.domain_reported_num)
        self.assertEqual(1, hit.domain_included_num)
        self.assertEqual("Importin-beta N-terminal domain", hit.description)
        hsp = hit.hsps[0]
        self.assertEqual(0.033, hsp.evalue)
        self.assertEqual(14.0, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)

        # last qresult
        qresult = next(qresults)
        counter += 1
        self.assertEqual(5, len(qresult))
        self.assertEqual("gi|125490392|ref|NP_038661.2|", qresult.id)
        self.assertEqual("-", qresult.accession)
        # first hit
        hit = qresult[0]
        self.assertEqual(1, len(hit))
        self.assertEqual("Pou", hit.id)
        self.assertEqual("PF00157.12", hit.accession)
        self.assertEqual(7e-37, hit.evalue)
        self.assertEqual(124.8, hit.bitscore)
        self.assertEqual(0.5, hit.bias)
        self.assertEqual(1.5, hit.domain_exp_num)
        self.assertEqual(1, hit.region_num)
        self.assertEqual(0, hit.cluster_num)
        self.assertEqual(0, hit.overlap_num)
        self.assertEqual(1, hit.env_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, hit.domain_reported_num)
        self.assertEqual(1, hit.domain_included_num)
        self.assertEqual("Pou domain - N-terminal to homeobox domain", hit.description)
        hsp = hit.hsps[0]
        self.assertEqual(1.4e-36, hsp.evalue)
        self.assertEqual(123.9, hsp.bitscore)
        self.assertEqual(0.3, hsp.bias)
        # second hit
        hit = qresult[1]
        self.assertEqual(1, len(hit))
        self.assertEqual("Homeobox", hit.id)
        self.assertEqual("PF00046.24", hit.accession)
        self.assertEqual(2.1e-18, hit.evalue)
        self.assertEqual(65.5, hit.bitscore)
        self.assertEqual(1.1, hit.bias)
        self.assertEqual(1.5, hit.domain_exp_num)
        self.assertEqual(1, hit.region_num)
        self.assertEqual(0, hit.cluster_num)
        self.assertEqual(0, hit.overlap_num)
        self.assertEqual(1, hit.env_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, hit.domain_reported_num)
        self.assertEqual(1, hit.domain_included_num)
        self.assertEqual("Homeobox domain", hit.description)
        hsp = hit.hsps[0]
        self.assertEqual(4.1e-18, hsp.evalue)
        self.assertEqual(64.6, hsp.bitscore)
        self.assertEqual(0.7, hsp.bias)
        # third hit
        hit = qresult[2]
        self.assertEqual(1, len(hit))
        self.assertEqual("HTH_31", hit.id)
        self.assertEqual("PF13560.1", hit.accession)
        self.assertEqual(0.012, hit.evalue)
        self.assertEqual(15.6, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual(2.2, hit.domain_exp_num)
        self.assertEqual(2, hit.region_num)
        self.assertEqual(0, hit.cluster_num)
        self.assertEqual(0, hit.overlap_num)
        self.assertEqual(2, hit.env_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, hit.domain_reported_num)
        self.assertEqual(0, hit.domain_included_num)
        self.assertEqual("Helix-turn-helix domain", hit.description)
        hsp = hit.hsps[0]
        self.assertEqual(0.16, hsp.evalue)
        self.assertEqual(12.0, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        # fourth hit
        hit = qresult[3]
        self.assertEqual(1, len(hit))
        self.assertEqual("Homeobox_KN", hit.id)
        self.assertEqual("PF05920.6", hit.accession)
        self.assertEqual(0.039, hit.evalue)
        self.assertEqual(13.5, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual(1.6, hit.domain_exp_num)
        self.assertEqual(1, hit.region_num)
        self.assertEqual(0, hit.cluster_num)
        self.assertEqual(0, hit.overlap_num)
        self.assertEqual(1, hit.env_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, hit.domain_reported_num)
        self.assertEqual(0, hit.domain_included_num)
        self.assertEqual("Homeobox KN domain", hit.description)
        hsp = hit.hsps[0]
        self.assertEqual(0.095, hsp.evalue)
        self.assertEqual(12.3, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        # fifth hit
        hit = qresult[4]
        self.assertEqual(1, len(hit))
        self.assertEqual("DUF521", hit.id)
        self.assertEqual("PF04412.8", hit.accession)
        self.assertEqual(0.14, hit.evalue)
        self.assertEqual(10.5, hit.bitscore)
        self.assertEqual(0.1, hit.bias)
        self.assertEqual(1.4, hit.domain_exp_num)
        self.assertEqual(1, hit.region_num)
        self.assertEqual(0, hit.cluster_num)
        self.assertEqual(0, hit.overlap_num)
        self.assertEqual(1, hit.env_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, hit.domain_reported_num)
        self.assertEqual(0, hit.domain_included_num)
        self.assertEqual("Protein of unknown function (DUF521)", hit.description)
        hsp = hit.hsps[0]
        self.assertEqual(0.26, hsp.evalue)
        self.assertEqual(9.6, hsp.bitscore)
        self.assertEqual(0.1, hsp.bias)

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(4, counter)

    def test_30_hmmscan_002(self):
        """Test parsing hmmer3-tab, hmmscan 3.0, single query, no hits (tab_30_hmmscan_002)."""
        tab_file = get_file("tab_30_hmmscan_002.out")
        qresults = parse(tab_file, FMT)

        self.assertRaises(StopIteration, next, qresults)

    def test_30_hmmscan_003(self):
        """Test parsing hmmer3-tab, hmmscan 3.0, single query, single hit, single hsp (tab_30_hmmscan_003)."""
        tab_file = get_file("tab_30_hmmscan_003.out")
        qresults = parse(tab_file, FMT)
        counter = 0

        qresult = next(qresults)
        counter += 1
        self.assertEqual(1, len(qresult))
        self.assertEqual("gi|4885477|ref|NP_005359.1|", qresult.id)
        self.assertEqual("-", qresult.accession)
        hit = qresult[0]
        self.assertEqual(1, len(hit))
        self.assertEqual("Globin", hit.id)
        self.assertEqual("PF00042.17", hit.accession)
        self.assertEqual(6e-21, hit.evalue)
        self.assertEqual(74.6, hit.bitscore)
        self.assertEqual(0.3, hit.bias)
        self.assertEqual(1.3, hit.domain_exp_num)
        self.assertEqual(1, hit.region_num)
        self.assertEqual(0, hit.cluster_num)
        self.assertEqual(0, hit.overlap_num)
        self.assertEqual(1, hit.env_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, hit.domain_reported_num)
        self.assertEqual(1, hit.domain_included_num)
        self.assertEqual("Globin", hit.description)
        hsp = hit.hsps[0]
        self.assertEqual(9.2e-21, hsp.evalue)
        self.assertEqual(74.0, hsp.bitscore)
        self.assertEqual(0.2, hsp.bias)

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(1, counter)

    def test_30_hmmscan_004(self):
        """Test parsing hmmer3-tab, hmmscan 3.0, single query, multiple hits (tab_30_hmmscan_004)."""
        tab_file = get_file("tab_30_hmmscan_004.out")
        qresults = parse(tab_file, FMT)
        counter = 0

        qresult = next(qresults)
        counter += 1
        self.assertEqual(2, len(qresult))
        self.assertEqual("gi|126362951:116-221", qresult.id)
        self.assertEqual("-", qresult.accession)
        hit = qresult[0]
        self.assertEqual(1, len(hit))
        self.assertEqual("Ig_3", hit.id)
        self.assertEqual("PF13927.1", hit.accession)
        self.assertEqual(1.4e-09, hit.evalue)
        self.assertEqual(38.2, hit.bitscore)
        self.assertEqual(0.4, hit.bias)
        self.assertEqual(1.3, hit.domain_exp_num)
        self.assertEqual(1, hit.region_num)
        self.assertEqual(0, hit.cluster_num)
        self.assertEqual(0, hit.overlap_num)
        self.assertEqual(1, hit.env_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, hit.domain_reported_num)
        self.assertEqual(1, hit.domain_included_num)
        self.assertEqual("Immunoglobulin domain", hit.description)
        hsp = hit.hsps[0]
        self.assertEqual(2.1e-09, hsp.evalue)
        self.assertEqual(37.6, hsp.bitscore)
        self.assertEqual(0.3, hsp.bias)
        hit = qresult[1]
        self.assertEqual(1, len(hit))
        self.assertEqual("Ig_2", hit.id)
        self.assertEqual("PF13895.1", hit.accession)
        self.assertEqual(3.5e-05, hit.evalue)
        self.assertEqual(23.7, hit.bitscore)
        self.assertEqual(0.1, hit.bias)
        self.assertEqual(1.1, hit.domain_exp_num)
        self.assertEqual(1, hit.region_num)
        self.assertEqual(0, hit.cluster_num)
        self.assertEqual(0, hit.overlap_num)
        self.assertEqual(1, hit.env_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, hit.domain_reported_num)
        self.assertEqual(1, hit.domain_included_num)
        self.assertEqual("Immunoglobulin domain", hit.description)
        hsp = hit.hsps[0]
        self.assertEqual(4.3e-05, hsp.evalue)
        self.assertEqual(23.4, hsp.bitscore)
        self.assertEqual(0.1, hsp.bias)

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(1, counter)


class HmmsearchCases(unittest.TestCase):
    """Tests for hmmsearch output."""

    def test_31b1_hmmsearch_001(self):
        """Test parsing hmmer3-tab, hmmsearch 3.1b1, multiple queries (tab_31b1_hmmscan_001)."""
        tab_file = get_file("tab_31b1_hmmsearch_001.out")
        qresults = list(parse(tab_file, FMT))
        self.assertEqual(1, len(qresults))

        # first qresult
        qresult = qresults[0]
        self.assertEqual(4, len(qresult))
        self.assertEqual("Pkinase", qresult.id)
        self.assertEqual("PF00069.17", qresult.accession)

        # first hit, first hsp
        hit = qresult[0]
        self.assertEqual(1, len(hit))
        self.assertEqual("sp|Q9WUT3|KS6A2_MOUSE", hit.id)
        self.assertEqual("-", hit.accession)
        self.assertEqual(8.5e-147, hit.evalue)
        self.assertEqual(492.3, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual(2.1, hit.domain_exp_num)
        self.assertEqual(2, hit.region_num)
        self.assertEqual(0, hit.cluster_num)
        self.assertEqual(0, hit.overlap_num)
        self.assertEqual(2, hit.env_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, hit.domain_reported_num)
        self.assertEqual(2, hit.domain_included_num)
        self.assertEqual(
            "Ribosomal protein S6 kinase alpha-2 OS=Mus musculus GN=Rps6ka2 PE=1 SV=1",
            hit.description,
        )
        hsp = hit.hsps[0]
        self.assertEqual(1.2e-72, hsp.evalue)
        self.assertEqual(249.3, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)

        # last hit, last hsp
        hit = qresult[-1]
        self.assertEqual(1, len(hit))
        self.assertEqual("sp|P18652|KS6AA_CHICK", hit.id)
        self.assertEqual("-", hit.accession)
        self.assertEqual(2.6e-145, hit.evalue)
        self.assertEqual(487.5, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual(2.1, hit.domain_exp_num)
        self.assertEqual(2, hit.region_num)
        self.assertEqual(0, hit.cluster_num)
        self.assertEqual(0, hit.overlap_num)
        self.assertEqual(2, hit.env_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, hit.domain_reported_num)
        self.assertEqual(2, hit.domain_included_num)
        self.assertEqual(
            "Ribosomal protein S6 kinase 2 alpha OS=Gallus gallus GN=RPS6KA PE=2 SV=1",
            hit.description,
        )
        hsp = hit.hsps[-1]
        self.assertEqual(7.6e-72, hsp.evalue)
        self.assertEqual(246.7, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
