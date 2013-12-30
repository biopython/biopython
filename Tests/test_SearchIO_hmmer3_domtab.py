# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO HmmerIO hmmer3-domtab parsers."""


import os
import unittest

from Bio.SearchIO import parse

# test case files are in the Blast directory
TEST_DIR = 'Hmmer'


def get_file(filename):
    """Returns the path of a test file."""
    return os.path.join(TEST_DIR, filename)


class HmmerscanCases(unittest.TestCase):

    fmt = 'hmmscan3-domtab'

    def test_domtab_30_hmmscan_001(self):
        "Test parsing hmmscan-domtab, hmmscan 3.0, multiple queries (domtab_30_hmmscan_001)"

        tab_file = get_file('domtab_30_hmmscan_001.out')
        qresults = parse(tab_file, self.fmt)
        counter = 0

        # first qresult
        qresult = next(qresults)
        counter += 1
        self.assertEqual(1, len(qresult))
        self.assertEqual('gi|4885477|ref|NP_005359.1|', qresult.id)
        self.assertEqual('-', qresult.accession)
        self.assertEqual(154, qresult.seq_len)
        hit = qresult[0]
        self.assertEqual(1, len(hit))
        self.assertEqual('Globin', hit.id)
        self.assertEqual('gi|4885477|ref|NP_005359.1|', hit.query_id)
        self.assertEqual('PF00042.17', hit.accession)
        self.assertEqual(108, hit.seq_len)
        self.assertEqual(6e-21, hit.evalue)
        self.assertEqual(74.6, hit.bitscore)
        self.assertEqual(0.3, hit.bias)
        self.assertEqual('Globin', hit.description)
        hsp = hit.hsps[0]
        self.assertEqual('Globin', hsp.hit_id)
        self.assertEqual('gi|4885477|ref|NP_005359.1|', hsp.query_id)
        self.assertEqual(1, hsp.domain_index)
        self.assertEqual(6.7e-25, hsp.evalue_cond)
        self.assertEqual(9.2e-21, hsp.evalue)
        self.assertEqual(74.0, hsp.bitscore)
        self.assertEqual(0.2, hsp.bias)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(107, hsp.hit_end)
        self.assertEqual(6, hsp.query_start)
        self.assertEqual(112, hsp.query_end)
        self.assertEqual(6, hsp.env_start)
        self.assertEqual(113, hsp.env_end)
        self.assertEqual(0.97, hsp.acc_avg)

        # second qresult
        qresult = next(qresults)
        counter += 1
        self.assertEqual(2, len(qresult))
        self.assertEqual('gi|126362951:116-221', qresult.id)
        self.assertEqual('-', qresult.accession)
        self.assertEqual(106, qresult.seq_len)
        hit = qresult[0]
        self.assertEqual(1, len(hit))
        self.assertEqual('Ig_3', hit.id)
        self.assertEqual('gi|126362951:116-221', hit.query_id)
        self.assertEqual('PF13927.1', hit.accession)
        self.assertEqual(75, hit.seq_len)
        self.assertEqual(1.4e-09, hit.evalue)
        self.assertEqual(38.2, hit.bitscore)
        self.assertEqual(0.4, hit.bias)
        self.assertEqual('Immunoglobulin domain', hit.description)
        hsp = hit.hsps[0]
        self.assertEqual('Ig_3', hsp.hit_id)
        self.assertEqual('gi|126362951:116-221', hsp.query_id)
        self.assertEqual(1, hsp.domain_index)
        self.assertEqual(3e-13, hsp.evalue_cond)
        self.assertEqual(2.1e-09, hsp.evalue)
        self.assertEqual(37.6, hsp.bitscore)
        self.assertEqual(0.3, hsp.bias)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(73, hsp.hit_end)
        self.assertEqual(8, hsp.query_start)
        self.assertEqual(84, hsp.query_end)
        self.assertEqual(8, hsp.env_start)
        self.assertEqual(88, hsp.env_end)
        self.assertEqual(0.94, hsp.acc_avg)
        hit = qresult[1]
        self.assertEqual(1, len(hit))
        self.assertEqual('Ig_2', hit.id)
        self.assertEqual('gi|126362951:116-221', hit.query_id)
        self.assertEqual('PF13895.1', hit.accession)
        self.assertEqual(80, hit.seq_len)
        self.assertEqual(3.5e-05, hit.evalue)
        self.assertEqual(23.7, hit.bitscore)
        self.assertEqual(0.1, hit.bias)
        self.assertEqual('Immunoglobulin domain', hit.description)
        hsp = hit.hsps[0]
        self.assertEqual('Ig_2', hsp.hit_id)
        self.assertEqual('gi|126362951:116-221', hsp.query_id)
        self.assertEqual(1, hsp.domain_index)
        self.assertEqual(6.2e-09, hsp.evalue_cond)
        self.assertEqual(4.3e-05, hsp.evalue)
        self.assertEqual(23.4, hsp.bitscore)
        self.assertEqual(0.1, hsp.bias)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(80, hsp.hit_end)
        self.assertEqual(8, hsp.query_start)
        self.assertEqual(104, hsp.query_end)
        self.assertEqual(8, hsp.env_start)
        self.assertEqual(104, hsp.env_end)
        self.assertEqual(0.71, hsp.acc_avg)

        # third qresult
        qresult = next(qresults)
        counter += 1
        self.assertEqual(2, len(qresult))
        self.assertEqual('gi|22748937|ref|NP_065801.1|', qresult.id)
        self.assertEqual('-', qresult.accession)
        self.assertEqual(1204, qresult.seq_len)
        hit = qresult[0]
        self.assertEqual(2, len(hit))
        self.assertEqual('Xpo1', hit.id)
        self.assertEqual('gi|22748937|ref|NP_065801.1|', hit.query_id)
        self.assertEqual('PF08389.7', hit.accession)
        self.assertEqual(148, hit.seq_len)
        self.assertEqual(7.8e-34, hit.evalue)
        self.assertEqual(116.6, hit.bitscore)
        self.assertEqual(7.8, hit.bias)
        self.assertEqual('Exportin 1-like protein', hit.description)
        hsp = hit.hsps[0]
        self.assertEqual('Xpo1', hsp.hit_id)
        self.assertEqual('gi|22748937|ref|NP_065801.1|', hsp.query_id)
        self.assertEqual(1, hsp.domain_index)
        self.assertEqual(1.6e-37, hsp.evalue_cond)
        self.assertEqual(1.1e-33, hsp.evalue)
        self.assertEqual(116.1, hsp.bitscore)
        self.assertEqual(3.4, hsp.bias)
        self.assertEqual(1, hsp.hit_start)
        self.assertEqual(148, hsp.hit_end)
        self.assertEqual(109, hsp.query_start)
        self.assertEqual(271, hsp.query_end)
        self.assertEqual(108, hsp.env_start)
        self.assertEqual(271, hsp.env_end)
        self.assertEqual(0.98, hsp.acc_avg)
        hsp = hit.hsps[1]
        self.assertEqual('Xpo1', hsp.hit_id)
        self.assertEqual('gi|22748937|ref|NP_065801.1|', hsp.query_id)
        self.assertEqual(2, hsp.domain_index)
        self.assertEqual(0.35, hsp.evalue_cond)
        self.assertEqual(2.4e+03, hsp.evalue)
        self.assertEqual(-1.8, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(111, hsp.hit_start)
        self.assertEqual(139, hsp.hit_end)
        self.assertEqual(498, hsp.query_start)
        self.assertEqual(525, hsp.query_end)
        self.assertEqual(495, hsp.env_start)
        self.assertEqual(529, hsp.env_end)
        self.assertEqual(0.86, hsp.acc_avg)
        # next hit in the third qresult
        hit = qresult[1]
        self.assertEqual(2, len(hit))
        self.assertEqual('IBN_N', hit.id)
        self.assertEqual('gi|22748937|ref|NP_065801.1|', hit.query_id)
        self.assertEqual('PF03810.14', hit.accession)
        self.assertEqual(77, hit.seq_len)
        self.assertEqual(0.0039, hit.evalue)
        self.assertEqual(16.9, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual('Importin-beta N-terminal domain', hit.description)
        hsp = hit.hsps[0]
        self.assertEqual('IBN_N', hsp.hit_id)
        self.assertEqual('gi|22748937|ref|NP_065801.1|', hsp.query_id)
        self.assertEqual(1, hsp.domain_index)
        self.assertEqual(4.8e-06, hsp.evalue_cond)
        self.assertEqual(0.033, hsp.evalue)
        self.assertEqual(14.0, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(3, hsp.hit_start)
        self.assertEqual(75, hsp.hit_end)
        self.assertEqual(35, hsp.query_start)
        self.assertEqual(98, hsp.query_end)
        self.assertEqual(32, hsp.env_start)
        self.assertEqual(100, hsp.env_end)
        self.assertEqual(0.87, hsp.acc_avg)
        hsp = hit.hsps[1]
        self.assertEqual('IBN_N', hsp.hit_id)
        self.assertEqual('gi|22748937|ref|NP_065801.1|', hsp.query_id)
        self.assertEqual(2, hsp.domain_index)
        self.assertEqual(1.2, hsp.evalue_cond)
        self.assertEqual(8e+03, hsp.evalue)
        self.assertEqual(-3.3, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(56, hsp.hit_start)
        self.assertEqual(75, hsp.hit_end)
        self.assertEqual(167, hsp.query_start)
        self.assertEqual(186, hsp.query_end)
        self.assertEqual(164, hsp.env_start)
        self.assertEqual(187, hsp.env_end)
        self.assertEqual(0.85, hsp.acc_avg)

        # fourth qresult
        qresult = next(qresults)
        counter += 1
        self.assertEqual(5, len(qresult))
        self.assertEqual('gi|125490392|ref|NP_038661.2|', qresult.id)
        self.assertEqual('-', qresult.accession)
        self.assertEqual(352, qresult.seq_len)
        # first hit
        hit = qresult[0]
        self.assertEqual(1, len(hit))
        self.assertEqual('Pou', hit.id)
        self.assertEqual('gi|125490392|ref|NP_038661.2|', hit.query_id)
        self.assertEqual('PF00157.12', hit.accession)
        self.assertEqual(75, hit.seq_len)
        self.assertEqual(7e-37, hit.evalue)
        self.assertEqual(124.8, hit.bitscore)
        self.assertEqual(0.5, hit.bias)
        self.assertEqual('Pou domain - N-terminal to homeobox domain', hit.description)
        hsp = hit.hsps[0]
        self.assertEqual('Pou', hsp.hit_id)
        self.assertEqual('gi|125490392|ref|NP_038661.2|', hsp.query_id)
        self.assertEqual(1, hsp.domain_index)
        self.assertEqual(5e-40, hsp.evalue_cond)
        self.assertEqual(1.4e-36, hsp.evalue)
        self.assertEqual(123.9, hsp.bitscore)
        self.assertEqual(0.3, hsp.bias)
        self.assertEqual(2, hsp.hit_start)
        self.assertEqual(75, hsp.hit_end)
        self.assertEqual(132, hsp.query_start)
        self.assertEqual(205, hsp.query_end)
        self.assertEqual(130, hsp.env_start)
        self.assertEqual(205, hsp.env_end)
        self.assertEqual(0.97, hsp.acc_avg)
        # second hit
        hit = qresult[1]
        self.assertEqual(1, len(hit))
        self.assertEqual('Homeobox', hit.id)
        self.assertEqual('gi|125490392|ref|NP_038661.2|', hit.query_id)
        self.assertEqual('PF00046.24', hit.accession)
        self.assertEqual(57, hit.seq_len)
        self.assertEqual(2.1e-18, hit.evalue)
        self.assertEqual(65.5, hit.bitscore)
        self.assertEqual(1.1, hit.bias)
        self.assertEqual('Homeobox domain', hit.description)
        hsp = hit.hsps[0]
        self.assertEqual('Homeobox', hsp.hit_id)
        self.assertEqual('gi|125490392|ref|NP_038661.2|', hsp.query_id)
        self.assertEqual(1, hsp.domain_index)
        self.assertEqual(1.5e-21, hsp.evalue_cond)
        self.assertEqual(4.1e-18, hsp.evalue)
        self.assertEqual(64.6, hsp.bitscore)
        self.assertEqual(0.7, hsp.bias)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(57, hsp.hit_end)
        self.assertEqual(223, hsp.query_start)
        self.assertEqual(280, hsp.query_end)
        self.assertEqual(223, hsp.env_start)
        self.assertEqual(280, hsp.env_end)
        self.assertEqual(0.98, hsp.acc_avg)
        # third hit
        hit = qresult[2]
        self.assertEqual(2, len(hit))
        self.assertEqual('HTH_31', hit.id)
        self.assertEqual('gi|125490392|ref|NP_038661.2|', hit.query_id)
        self.assertEqual('PF13560.1', hit.accession)
        self.assertEqual(64, hit.seq_len)
        self.assertEqual(0.012, hit.evalue)
        self.assertEqual(15.6, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual('Helix-turn-helix domain', hit.description)
        hsp = hit.hsps[0]
        self.assertEqual('HTH_31', hsp.hit_id)
        self.assertEqual('gi|125490392|ref|NP_038661.2|', hsp.query_id)
        self.assertEqual(1, hsp.domain_index)
        self.assertEqual(5.7e-05, hsp.evalue_cond)
        self.assertEqual(0.16, hsp.evalue)
        self.assertEqual(12.0, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(35, hsp.hit_end)
        self.assertEqual(140, hsp.query_start)
        self.assertEqual(181, hsp.query_end)
        self.assertEqual(140, hsp.env_start)
        self.assertEqual(184, hsp.env_end)
        self.assertEqual(0.96, hsp.acc_avg)
        hsp = hit.hsps[1]
        self.assertEqual('HTH_31', hsp.hit_id)
        self.assertEqual('gi|125490392|ref|NP_038661.2|', hsp.query_id)
        self.assertEqual(2, hsp.domain_index)
        self.assertEqual(0.19, hsp.evalue_cond)
        self.assertEqual(5.2e+02, hsp.evalue)
        self.assertEqual(0.8, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(38, hsp.hit_start)
        self.assertEqual(62, hsp.hit_end)
        self.assertEqual(244, hsp.query_start)
        self.assertEqual(268, hsp.query_end)
        self.assertEqual(242, hsp.env_start)
        self.assertEqual(270, hsp.env_end)
        self.assertEqual(0.86, hsp.acc_avg)
        # fourth hit
        hit = qresult[3]
        self.assertEqual(1, len(hit))
        self.assertEqual('Homeobox_KN', hit.id)
        self.assertEqual('gi|125490392|ref|NP_038661.2|', hit.query_id)
        self.assertEqual('PF05920.6', hit.accession)
        self.assertEqual(40, hit.seq_len)
        self.assertEqual(0.039, hit.evalue)
        self.assertEqual(13.5, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual('Homeobox KN domain', hit.description)
        hsp = hit.hsps[0]
        self.assertEqual('Homeobox_KN', hsp.hit_id)
        self.assertEqual('gi|125490392|ref|NP_038661.2|', hsp.query_id)
        self.assertEqual(1, hsp.domain_index)
        self.assertEqual(3.5e-05, hsp.evalue_cond)
        self.assertEqual(0.095, hsp.evalue)
        self.assertEqual(12.3, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(6, hsp.hit_start)
        self.assertEqual(39, hsp.hit_end)
        self.assertEqual(243, hsp.query_start)
        self.assertEqual(276, hsp.query_end)
        self.assertEqual(240, hsp.env_start)
        self.assertEqual(277, hsp.env_end)
        self.assertEqual(0.91, hsp.acc_avg)
        # fifth hit
        hit = qresult[4]
        self.assertEqual(1, len(hit))
        self.assertEqual('DUF521', hit.id)
        self.assertEqual('gi|125490392|ref|NP_038661.2|', hit.query_id)
        self.assertEqual('PF04412.8', hit.accession)
        self.assertEqual(400, hit.seq_len)
        self.assertEqual(0.14, hit.evalue)
        self.assertEqual(10.5, hit.bitscore)
        self.assertEqual(0.1, hit.bias)
        self.assertEqual('Protein of unknown function (DUF521)', hit.description)
        hsp = hit.hsps[0]
        self.assertEqual('DUF521', hsp.hit_id)
        self.assertEqual('gi|125490392|ref|NP_038661.2|', hsp.query_id)
        self.assertEqual(1, hsp.domain_index)
        self.assertEqual(9.4e-05, hsp.evalue_cond)
        self.assertEqual(0.26, hsp.evalue)
        self.assertEqual(9.6, hsp.bitscore)
        self.assertEqual(0.1, hsp.bias)
        self.assertEqual(272, hsp.hit_start)
        self.assertEqual(334, hsp.hit_end)
        self.assertEqual(220, hsp.query_start)
        self.assertEqual(280, hsp.query_end)
        self.assertEqual(196, hsp.env_start)
        self.assertEqual(294, hsp.env_end)
        self.assertEqual(0.77, hsp.acc_avg)

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(4, counter)

    def test_domtab_30_hmmscan_002(self):
        "Test parsing hmmscan-domtab, hmmscan 3.0, single query, no hits (domtab_30_hmmscan_002)"

        tab_file = get_file('domtab_30_hmmscan_002.out')
        qresults = parse(tab_file, self.fmt)

        self.assertRaises(StopIteration, next, qresults)

    def test_domtab_30_hmmscan_003(self):
        "Test parsing hmmscan-domtab, hmmscan 3.0, multiple queries (domtab_30_hmmscan_003)"

        tab_file = get_file('domtab_30_hmmscan_003.out')
        qresults = parse(tab_file, self.fmt)
        counter = 0

        qresult = next(qresults)
        counter += 1
        self.assertEqual(1, len(qresult))
        self.assertEqual('gi|4885477|ref|NP_005359.1|', qresult.id)
        self.assertEqual('-', qresult.accession)
        self.assertEqual(154, qresult.seq_len)
        hit = qresult[0]
        self.assertEqual(1, len(hit))
        self.assertEqual('Globin', hit.id)
        self.assertEqual('gi|4885477|ref|NP_005359.1|', hit.query_id)
        self.assertEqual('PF00042.17', hit.accession)
        self.assertEqual(108, hit.seq_len)
        self.assertEqual(6e-21, hit.evalue)
        self.assertEqual(74.6, hit.bitscore)
        self.assertEqual(0.3, hit.bias)
        self.assertEqual('Globin', hit.description)
        hsp = hit.hsps[0]
        self.assertEqual('Globin', hsp.hit_id)
        self.assertEqual('gi|4885477|ref|NP_005359.1|', hsp.query_id)
        self.assertEqual(1, hsp.domain_index)
        self.assertEqual(6.7e-25, hsp.evalue_cond)
        self.assertEqual(9.2e-21, hsp.evalue)
        self.assertEqual(74.0, hsp.bitscore)
        self.assertEqual(0.2, hsp.bias)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(107, hsp.hit_end)
        self.assertEqual(6, hsp.query_start)
        self.assertEqual(112, hsp.query_end)
        self.assertEqual(6, hsp.env_start)
        self.assertEqual(113, hsp.env_end)
        self.assertEqual(0.97, hsp.acc_avg)

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(1, counter)

    def test_domtab_30_hmmscan_004(self):
        "Test parsing hmmscan-domtab, hmmscan 3.0, multiple queries (domtab_30_hmmscan_004)"

        tab_file = get_file('domtab_30_hmmscan_004.out')
        qresults = parse(tab_file, self.fmt)
        counter = 0

        qresult = next(qresults)
        counter += 1
        self.assertEqual(2, len(qresult))
        self.assertEqual('gi|126362951:116-221', qresult.id)
        self.assertEqual('-', qresult.accession)
        self.assertEqual(106, qresult.seq_len)
        hit = qresult[0]
        self.assertEqual(1, len(hit))
        self.assertEqual('Ig_3', hit.id)
        self.assertEqual('gi|126362951:116-221', hit.query_id)
        self.assertEqual('PF13927.1', hit.accession)
        self.assertEqual(75, hit.seq_len)
        self.assertEqual(1.4e-09, hit.evalue)
        self.assertEqual(38.2, hit.bitscore)
        self.assertEqual(0.4, hit.bias)
        self.assertEqual('Immunoglobulin domain', hit.description)
        hsp = hit.hsps[0]
        self.assertEqual('Ig_3', hsp.hit_id)
        self.assertEqual('gi|126362951:116-221', hsp.query_id)
        self.assertEqual(1, hsp.domain_index)
        self.assertEqual(3e-13, hsp.evalue_cond)
        self.assertEqual(2.1e-09, hsp.evalue)
        self.assertEqual(37.6, hsp.bitscore)
        self.assertEqual(0.3, hsp.bias)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(73, hsp.hit_end)
        self.assertEqual(8, hsp.query_start)
        self.assertEqual(84, hsp.query_end)
        self.assertEqual(8, hsp.env_start)
        self.assertEqual(88, hsp.env_end)
        self.assertEqual(0.94, hsp.acc_avg)
        hit = qresult[1]
        self.assertEqual(1, len(hit))
        self.assertEqual('Ig_2', hit.id)
        self.assertEqual('gi|126362951:116-221', hit.query_id)
        self.assertEqual('PF13895.1', hit.accession)
        self.assertEqual(80, hit.seq_len)
        self.assertEqual(3.5e-05, hit.evalue)
        self.assertEqual(23.7, hit.bitscore)
        self.assertEqual(0.1, hit.bias)
        self.assertEqual('Immunoglobulin domain', hit.description)
        hsp = hit.hsps[0]
        self.assertEqual('Ig_2', hsp.hit_id)
        self.assertEqual('gi|126362951:116-221', hsp.query_id)
        self.assertEqual(1, hsp.domain_index)
        self.assertEqual(6.2e-09, hsp.evalue_cond)
        self.assertEqual(4.3e-05, hsp.evalue)
        self.assertEqual(23.4, hsp.bitscore)
        self.assertEqual(0.1, hsp.bias)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(80, hsp.hit_end)
        self.assertEqual(8, hsp.query_start)
        self.assertEqual(104, hsp.query_end)
        self.assertEqual(8, hsp.env_start)
        self.assertEqual(104, hsp.env_end)
        self.assertEqual(0.71, hsp.acc_avg)

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(1, counter)


class HmmersearchCases(unittest.TestCase):

    fmt = 'hmmsearch3-domtab'

    def test_domtab_30_hmmsearch_001(self):
        "Test parsing hmmsearch-domtab, hmmsearch 3.0, multiple queries (domtab_30_hmmsearch_001)"

        tab_file = get_file('domtab_30_hmmsearch_001.out')
        qresults = parse(tab_file, self.fmt)

        # first qresult
        # we only want to check the coordinate switch actually
        # so checking the first hsp of the first hit of the qresult is enough
        qresult = next(qresults)
        self.assertEqual(7, len(qresult))
        self.assertEqual('Pkinase', qresult.id)
        self.assertEqual('PF00069.17', qresult.accession)
        self.assertEqual(260, qresult.seq_len)
        hit = qresult[0]
        self.assertEqual(2, len(hit))
        self.assertEqual('sp|Q9WUT3|KS6A2_MOUSE', hit.id)
        self.assertEqual('Pkinase', hit.query_id)
        self.assertEqual('-', hit.accession)
        self.assertEqual(733, hit.seq_len)
        self.assertEqual(8.4e-147, hit.evalue)
        self.assertEqual(492.3, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual('Ribosomal protein S6 kinase alpha-2 OS=Mus musculus GN=Rps6ka2 PE=2 SV=1', hit.description)
        hsp = hit.hsps[0]
        self.assertEqual('sp|Q9WUT3|KS6A2_MOUSE', hsp.hit_id)
        self.assertEqual('Pkinase', hsp.query_id)
        self.assertEqual(1, hsp.domain_index)
        self.assertEqual(4.6e-75, hsp.evalue_cond)
        self.assertEqual(3.5e-70, hsp.evalue)
        self.assertEqual(241.2, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(58, hsp.hit_start)
        self.assertEqual(318, hsp.hit_end)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(260, hsp.query_end)
        self.assertEqual(58, hsp.env_start)
        self.assertEqual(318, hsp.env_end)
        self.assertEqual(0.95, hsp.acc_avg)

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
