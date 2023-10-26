# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO HmmerIO parsers."""


import os
import unittest

from Bio.SearchIO import parse


# test case files are in the Blast directory
TEST_DIR = "Hmmer"
FMT = "hmmer3-text"


def get_file(filename):
    """Return the path of a test file."""
    return os.path.join(TEST_DIR, filename)


class HmmscanCases(unittest.TestCase):
    """Testing hmmscan output."""

    def test_31b1_hmmscan_001(self):
        """Parsing hmmscan 3.1b1 (text_31b1_hmmscan_001)."""
        txt_file = get_file("text_31b1_hmmscan_001.out")
        qresults = parse(txt_file, FMT)
        counter = 0

        # test first qresult
        qresult = next(qresults)
        counter += 1

        self.assertEqual("hmmscan", qresult.program)
        self.assertEqual("/home/bow/db/hmmer/protdb/Pfam-A.hmm", qresult.target)
        self.assertEqual("3.1b1", qresult.version)
        self.assertEqual("random_s00", qresult.id)
        self.assertEqual(22, qresult.seq_len)
        self.assertEqual(0, len(qresult))

        # test fifth result
        for qresult in qresults:
            counter += 1

        self.assertEqual("hmmscan", qresult.program)
        self.assertEqual("/home/bow/db/hmmer/protdb/Pfam-A.hmm", qresult.target)
        self.assertEqual("3.1b1", qresult.version)
        self.assertEqual("gi|125490392|ref|NP_038661.2|", qresult.id)
        self.assertEqual(
            "POU domain, class 5, transcription factor 1 isoform 1 [Mus musculus]",
            qresult.description,
        )
        self.assertEqual(352, qresult.seq_len)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual("Pou", hit.id)
        self.assertEqual("Pou domain - N-terminal to homeobox domain", hit.description)
        self.assertTrue(hit.is_included)
        self.assertEqual(7.6e-37, hit.evalue)
        self.assertEqual(124.8, hit.bitscore)
        self.assertEqual(0.5, hit.bias)
        self.assertEqual(1.5, hit.domain_exp_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        hsp = hit.hsps[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_included)
        self.assertEqual(123.9, hsp.bitscore)
        self.assertEqual(0.5, hsp.bias)
        self.assertEqual(5e-40, hsp.evalue_cond)
        self.assertEqual(1.5e-36, hsp.evalue)
        self.assertEqual(2, hsp.hit_start)
        self.assertEqual(75, hsp.hit_end)
        self.assertEqual(".]", hsp.hit_endtype)
        self.assertEqual(132, hsp.query_start)
        self.assertEqual(205, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual(130, hsp.env_start)
        self.assertEqual(205, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.97, hsp.acc_avg)
        self.assertEqual(
            "eldleeleefakefkqrrikLgltqadvgsalgalyGkefsqttIcrFEalqLslknmckLkpllekWLeeae",
            hsp.hit.seq,
        )
        self.assertEqual(
            "KALQKELEQFAKLLKQKRITLGYTQADVGLTLGVLFGKVFSQTTICRFEALQLSLKNMCKLRPLLEKWVEEAD",
            hsp.query.seq,
        )
        self.assertEqual(
            "67899******************************************************************96",
            hsp.aln_annotation["PP"],
        )

        # last hit
        hit = qresult[4]
        self.assertEqual("DUF521", hit.id)
        self.assertEqual("Protein of unknown function (DUF521)", hit.description)
        self.assertFalse(hit.is_included)
        self.assertEqual(0.15, hit.evalue)
        self.assertEqual(10.5, hit.bitscore)
        self.assertEqual(0.1, hit.bias)
        self.assertEqual(1.4, hit.domain_exp_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        hsp = hit.hsps[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertFalse(hsp.is_included)
        self.assertEqual(9.6, hsp.bitscore)
        self.assertEqual(0.1, hsp.bias)
        self.assertEqual(9.4e-05, hsp.evalue_cond)
        self.assertEqual(0.28, hsp.evalue)
        self.assertEqual(272, hsp.hit_start)
        self.assertEqual(334, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)
        self.assertEqual(220, hsp.query_start)
        self.assertEqual(280, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual(196, hsp.env_start)
        self.assertEqual(294, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.77, hsp.acc_avg)
        self.assertEqual(
            "adlaavleelnkakkeevdlvvlGcPhlsleeleelaellkgrkkkvsvelvvttsravlsk",
            hsp.hit.seq,
        )
        self.assertEqual(
            "QARKRKRTSIENRVRWSLETMFLKCPKPSLQQITHIANQLGLEK--DVVRVWFCNRRQKGKR",
            hsp.query.seq,
        )
        self.assertEqual(
            "345666667778888899************************99..9999999988876554",
            hsp.aln_annotation["PP"],
        )

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(5, counter)

    def test_30_hmmscan_001(self):
        """Parsing hmmscan 3.0 (text_30_hmmscan_001)."""
        txt_file = get_file("text_30_hmmscan_001.out")
        qresults = parse(txt_file, FMT)
        counter = 0

        # test first qresult
        qresult = next(qresults)
        counter += 1

        self.assertEqual("hmmscan", qresult.program)
        self.assertEqual("/home/bow/db/hmmer/Pfam-A.hmm", qresult.target)
        self.assertEqual("3.0", qresult.version)
        self.assertEqual("random_s00", qresult.id)
        self.assertEqual(32, qresult.seq_len)
        self.assertEqual(0, len(qresult))

        # test second result
        qresult = next(qresults)
        counter += 1

        self.assertEqual("hmmscan", qresult.program)
        self.assertEqual("/home/bow/db/hmmer/Pfam-A.hmm", qresult.target)
        self.assertEqual("3.0", qresult.version)
        self.assertEqual("gi|4885477|ref|NP_005359.1|", qresult.id)
        self.assertEqual("myoglobin [Homo sapiens]", qresult.description)
        self.assertEqual(154, qresult.seq_len)
        self.assertEqual(1, len(qresult))

        hit = qresult[0]
        self.assertEqual("Globin", hit.id)
        self.assertEqual("Globin", hit.description)
        self.assertTrue(hit.is_included)
        self.assertEqual(6e-21, hit.evalue)
        self.assertEqual(74.6, hit.bitscore)
        self.assertEqual(0.3, hit.bias)
        self.assertEqual(1.3, hit.domain_exp_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        hsp = hit.hsps[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_included)
        self.assertEqual(74.0, hsp.bitscore)
        self.assertEqual(0.2, hsp.bias)
        self.assertEqual(6.7e-25, hsp.evalue_cond)
        self.assertEqual(9.2e-21, hsp.evalue)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(107, hsp.hit_end)
        self.assertEqual("[.", hsp.hit_endtype)
        self.assertEqual(6, hsp.query_start)
        self.assertEqual(112, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual(6, hsp.env_start)
        self.assertEqual(113, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.97, hsp.acc_avg)
        self.assertEqual(
            "HHHHHHHHHHHHCHHHHHHHHHHHHHHHHHSGGGGGGGCCCTTTT.HHHHHTSCHHHHHHHHHHHHHHHHHHCTTSHHHHHHHHHHHHHHHHTT-.--HHHHCCHHHHH",
            hsp.aln_annotation["CS"],
        )
        self.assertEqual(
            "qkalvkaswekvkanaeeigaeilkrlfkaypdtkklFkkfgdls.aedlksspkfkahakkvlaaldeavknldnddnlkaalkklgarHakrg.vdpanfklfgeal",
            hsp.hit.seq,
        )
        self.assertEqual(
            "EWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKsEDEMKASEDLKKHGATVLTALGGILKK---KGHHEAEIKPLAQSHATKHkIPVKYLEFISECI",
            hsp.query.seq,
        )
        self.assertEqual(
            "5789*********************************************************************...6899***********************999998",
            hsp.aln_annotation["PP"],
        )

        # test third result
        qresult = next(qresults)
        counter += 1

        self.assertEqual("hmmscan", qresult.program)
        self.assertEqual("/home/bow/db/hmmer/Pfam-A.hmm", qresult.target)
        self.assertEqual("3.0", qresult.version)
        self.assertEqual("gi|126362951:116-221", qresult.id)
        self.assertEqual(
            "leukocyte immunoglobulin-like receptor subfamily B member 1 isoform 2 precursor [Homo sapiens]",
            qresult.description,
        )
        self.assertEqual(106, qresult.seq_len)
        self.assertEqual(2, len(qresult))

        hit = qresult[0]
        self.assertEqual("Ig_3", hit.id)
        self.assertEqual("Immunoglobulin domain", hit.description)
        self.assertTrue(hit.is_included)
        self.assertEqual(1.4e-09, hit.evalue)
        self.assertEqual(38.2, hit.bitscore)
        self.assertEqual(0.4, hit.bias)
        self.assertEqual(1.3, hit.domain_exp_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        hsp = hit.hsps[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_included)
        self.assertEqual(37.6, hsp.bitscore)
        self.assertEqual(0.3, hsp.bias)
        self.assertEqual(3e-13, hsp.evalue_cond)
        self.assertEqual(2.1e-09, hsp.evalue)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(73, hsp.hit_end)
        self.assertEqual("[.", hsp.hit_endtype)
        self.assertEqual(8, hsp.query_start)
        self.assertEqual(84, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual(8, hsp.env_start)
        self.assertEqual(88, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.94, hsp.acc_avg)
        self.assertEqual(
            "kPvisvspsptvtsggnvtLtCsaeggpppptisWy.....ietppelqgsegssssestLtissvtsedsgtYtCva",
            hsp.hit.seq,
        )
        self.assertEqual(
            "KPTLSAQPSPVVNSGGNVILQCDSQVA--FDGFSLCkegedEHPQCLNSQPHARGSSRAIFSVGPVSPSRRWWYRCYA",
            hsp.query.seq,
        )
        self.assertEqual(
            "8************************99..78888888****************************************9",
            hsp.aln_annotation["PP"],
        )

        hit = qresult[-1]
        self.assertEqual("Ig_2", hit.id)
        self.assertEqual("Immunoglobulin domain", hit.description)
        self.assertEqual(3.5e-05, hit.evalue)
        self.assertEqual(23.7, hit.bitscore)
        self.assertEqual(0.1, hit.bias)
        self.assertEqual(1.1, hit.domain_exp_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        hsp = hit.hsps[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_included)
        self.assertEqual(23.4, hsp.bitscore)
        self.assertEqual(0.1, hsp.bias)
        self.assertEqual(6.2e-09, hsp.evalue_cond)
        self.assertEqual(4.3e-05, hsp.evalue)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(80, hsp.hit_end)
        self.assertEqual("[]", hsp.hit_endtype)
        self.assertEqual(8, hsp.query_start)
        self.assertEqual(104, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual(8, hsp.env_start)
        self.assertEqual(104, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.71, hsp.acc_avg)
        self.assertEqual(
            "kpvlvapp.svvtegenvtLtCsapgnptprvqwykdg.vels......qsqnq........lfipnvsaedsgtYtCra....rnseggktstsveltv",
            hsp.hit.seq,
        )
        self.assertEqual(
            "KPTLSAQPsPVVNSGGNVILQCDSQVA-FDGFSLCKEGeDEHPqclnsqP---HargssraiFSVGPVSPSRRWWYRCYAydsnSPYEWSLPSDLLELLV",
            hsp.query.seq,
        )
        self.assertEqual(
            "799998885779*************85.899***9988655554443320...134455543444669************88443344588888888766",
            hsp.aln_annotation["PP"],
        )

        # test fourth result
        qresult = next(qresults)
        counter += 1

        self.assertEqual("hmmscan", qresult.program)
        self.assertEqual("/home/bow/db/hmmer/Pfam-A.hmm", qresult.target)
        self.assertEqual("3.0", qresult.version)
        self.assertEqual("gi|22748937|ref|NP_065801.1|", qresult.id)
        self.assertEqual("exportin-5 [Homo sapiens]", qresult.description)
        self.assertEqual(1204, qresult.seq_len)
        self.assertEqual(2, len(qresult))

        hit = qresult[0]
        self.assertEqual("Xpo1", hit.id)
        self.assertEqual("Exportin 1-like protein", hit.description)
        self.assertTrue(hit.is_included)
        self.assertEqual(7.8e-34, hit.evalue)
        self.assertEqual(116.6, hit.bitscore)
        self.assertEqual(7.8, hit.bias)
        self.assertEqual(2.8, hit.domain_exp_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))

        hsp = hit.hsps[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_included)
        self.assertEqual(116.1, hsp.bitscore)
        self.assertEqual(3.4, hsp.bias)
        self.assertEqual(1.6e-37, hsp.evalue_cond)
        self.assertEqual(1.1e-33, hsp.evalue)
        self.assertEqual(1, hsp.hit_start)
        self.assertEqual(148, hsp.hit_end)
        self.assertEqual(".]", hsp.hit_endtype)
        self.assertEqual(109, hsp.query_start)
        self.assertEqual(271, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual(108, hsp.env_start)
        self.assertEqual(271, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.98, hsp.acc_avg)
        self.assertEqual(
            "HHHHHHHHHHHHHHHHHHTTTTSTTHHHHHHHHHHG-HHHHHHHHHHHHHHHHHHCCS-TTTS-CCCHHHHHHHCHHHHHHHHHHHHHHHC-TT-..................HHHHHHHHHHHHHHCTTS-CHHCHCS...HHHHHCHHCCSCCCHHHHHHHH",
            hsp.aln_annotation["CS"],
        )
        self.assertEqual(
            "kflrnklaealaelflqeypnqWpsffddllsllssspsglelllriLkvlpeEiadfsrskleqerrnelkdllrsqvqkilelllqileqsvskk...............sselveatLkclsswvswidiglivnsp..llsllfqlLndpelreaAvecL",
            hsp.hit.seq,
        )
        self.assertEqual(
            "NHIKDALSRIVVEMIKREWPQHWPDMLIELDTLSKQGETQTELVMFILLRLAEDVVTF--QTLPPQRRRDIQQTLTQNMERIFSFLLNTLQENVNKYqqvktdtsqeskaqaNCRVGVAALNTLAGYIDWVSMSHITAENckLLEILCLLLNEQELQLGAAECL",
            hsp.query.seq,
        )
        self.assertEqual(
            "89******************************************************99..79*********************************99*****************************************8889*********************8",
            hsp.aln_annotation["PP"],
        )

        hsp = hit.hsps[-1]
        self.assertEqual(2, hsp.domain_index)
        self.assertFalse(hsp.is_included)
        self.assertEqual(-1.8, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(0.35, hsp.evalue_cond)
        self.assertEqual(2.4e03, hsp.evalue)
        self.assertEqual(111, hsp.hit_start)
        self.assertEqual(139, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)
        self.assertEqual(498, hsp.query_start)
        self.assertEqual(525, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual(495, hsp.env_start)
        self.assertEqual(529, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.86, hsp.acc_avg)
        self.assertEqual("HHCTTS-CHHCHCS.HHHHHCHHCCSCC", hsp.aln_annotation["CS"])
        self.assertEqual("swvswidiglivnspllsllfqlLndpe", hsp.hit.seq)
        self.assertEqual("SFVQWEAMTLFLES-VITQMFRTLNREE", hsp.query.seq)
        self.assertEqual("899*********98.8888899998776", hsp.aln_annotation["PP"])

        hit = qresult[-1]
        self.assertEqual("IBN_N", hit.id)
        self.assertEqual("Importin-beta N-terminal domain", hit.description)
        self.assertTrue(hit.is_included)
        self.assertEqual(0.0039, hit.evalue)
        self.assertEqual(16.9, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual(2.7, hit.domain_exp_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))

        hsp = hit.hsps[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_included)
        self.assertEqual(14.0, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(4.8e-06, hsp.evalue_cond)
        self.assertEqual(0.033, hsp.evalue)
        self.assertEqual(3, hsp.hit_start)
        self.assertEqual(75, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)
        self.assertEqual(35, hsp.query_start)
        self.assertEqual(98, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual(32, hsp.env_start)
        self.assertEqual(100, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.87, hsp.acc_avg)
        self.assertEqual(
            "HHHHHHHSCTHHHHHHHHHHHTTTSTHHHHHHHHHHHHHHHHHSCCHHHHHHHHCS-HHHHHHHHHHHHHHH",
            hsp.aln_annotation["CS"],
        )
        self.assertEqual(
            "qLnqlekqkPgflsallqilanksldlevRqlAalyLknlItkhWkseeaqrqqqlpeeekelIrnnllnll",
            hsp.hit.seq,
        )
        self.assertEqual(
            "FCEEFKEKCPICVPCGLRLA-EKTQVAIVRHFGLQILEHVVKFRWN--------GMSRLEKVYLKNSVMELI",
            hsp.query.seq,
        )
        self.assertEqual(
            "56788886699*********.6555899******************........999999****99999887",
            hsp.aln_annotation["PP"],
        )

        hsp = hit.hsps[-1]
        self.assertEqual(2, hsp.domain_index)
        self.assertFalse(hsp.is_included)
        self.assertEqual(-3.3, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(1.2, hsp.evalue_cond)
        self.assertEqual(8e03, hsp.evalue)
        self.assertEqual(56, hsp.hit_start)
        self.assertEqual(75, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)
        self.assertEqual(167, hsp.query_start)
        self.assertEqual(186, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual(164, hsp.env_start)
        self.assertEqual(187, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.85, hsp.acc_avg)
        self.assertEqual("HCS-HHHHHHHHHHHHHHH", hsp.aln_annotation["CS"])
        self.assertEqual("qqlpeeekelIrnnllnll", hsp.hit.seq)
        self.assertEqual("QTLPPQRRRDIQQTLTQNM", hsp.query.seq)
        self.assertEqual("6899*******99998865", hsp.aln_annotation["PP"])

        # test fifth result
        qresult = next(qresults)
        counter += 1

        self.assertEqual("hmmscan", qresult.program)
        self.assertEqual("/home/bow/db/hmmer/Pfam-A.hmm", qresult.target)
        self.assertEqual("3.0", qresult.version)
        self.assertEqual("gi|125490392|ref|NP_038661.2|", qresult.id)
        self.assertEqual(
            "POU domain, class 5, transcription factor 1 isoform 1 [Mus musculus]",
            qresult.description,
        )
        self.assertEqual(352, qresult.seq_len)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual("Pou", hit.id)
        self.assertEqual("Pou domain - N-terminal to homeobox domain", hit.description)
        self.assertTrue(hit.is_included)
        self.assertEqual(7e-37, hit.evalue)
        self.assertEqual(124.8, hit.bitscore)
        self.assertEqual(0.5, hit.bias)
        self.assertEqual(1.5, hit.domain_exp_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        hsp = hit.hsps[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_included)
        self.assertEqual(123.9, hsp.bitscore)
        self.assertEqual(0.3, hsp.bias)
        self.assertEqual(5e-40, hsp.evalue_cond)
        self.assertEqual(1.4e-36, hsp.evalue)
        self.assertEqual(2, hsp.hit_start)
        self.assertEqual(75, hsp.hit_end)
        self.assertEqual(".]", hsp.hit_endtype)
        self.assertEqual(132, hsp.query_start)
        self.assertEqual(205, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual(130, hsp.env_start)
        self.assertEqual(205, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.97, hsp.acc_avg)
        self.assertEqual(
            "eldleeleefakefkqrrikLgltqadvgsalgalyGkefsqttIcrFEalqLslknmckLkpllekWLeeae",
            hsp.hit.seq,
        )
        self.assertEqual(
            "KALQKELEQFAKLLKQKRITLGYTQADVGLTLGVLFGKVFSQTTICRFEALQLSLKNMCKLRPLLEKWVEEAD",
            hsp.query.seq,
        )
        self.assertEqual(
            "67899******************************************************************96",
            hsp.aln_annotation["PP"],
        )

        hit = qresult[1]
        self.assertEqual("Homeobox", hit.id)
        self.assertEqual("Homeobox domain", hit.description)
        self.assertTrue(hit.is_included)
        self.assertEqual(2.1e-18, hit.evalue)
        self.assertEqual(65.5, hit.bitscore)
        self.assertEqual(1.1, hit.bias)
        self.assertEqual(1.5, hit.domain_exp_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        hsp = hit.hsps[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_included)
        self.assertEqual(64.6, hsp.bitscore)
        self.assertEqual(0.7, hsp.bias)
        self.assertEqual(1.5e-21, hsp.evalue_cond)
        self.assertEqual(4.1e-18, hsp.evalue)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(57, hsp.hit_end)
        self.assertEqual("[]", hsp.hit_endtype)
        self.assertEqual(223, hsp.query_start)
        self.assertEqual(280, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual(223, hsp.env_start)
        self.assertEqual(280, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.98, hsp.acc_avg)
        self.assertEqual(
            "SS--SS--HHHHHHHHHHCCTSSS--HHHHHHHHHH----HHHHHHHHHHHHHHHHH",
            hsp.aln_annotation["CS"],
        )
        self.assertEqual(
            "rrkRttftkeqleeLeelFeknrypsaeereeLAkklgLterqVkvWFqNrRakekk", hsp.hit.seq
        )
        self.assertEqual(
            "KRKRTSIENRVRWSLETMFLKCPKPSLQQITHIANQLGLEKDVVRVWFCNRRQKGKR", hsp.query.seq
        )
        self.assertEqual(
            "79****************************************************997",
            hsp.aln_annotation["PP"],
        )

        hit = qresult[2]
        self.assertEqual("HTH_31", hit.id)
        self.assertEqual("Helix-turn-helix domain", hit.description)
        self.assertFalse(hit.is_included)
        self.assertEqual(0.012, hit.evalue)
        self.assertEqual(15.6, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual(2.2, hit.domain_exp_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))

        hsp = hit.hsps[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertFalse(hsp.is_included)
        self.assertEqual(12.0, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(5.7e-05, hsp.evalue_cond)
        self.assertEqual(0.16, hsp.evalue)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(35, hsp.hit_end)
        self.assertEqual("[.", hsp.hit_endtype)
        self.assertEqual(140, hsp.query_start)
        self.assertEqual(181, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual(140, hsp.env_start)
        self.assertEqual(184, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.96, hsp.acc_avg)
        self.assertEqual("aLGarLralReraGLtqeevAerlg......vSastlsrlE", hsp.hit.seq)
        self.assertEqual("QFAKLLKQKRITLGYTQADVGLTLGvlfgkvFSQTTICRFE", hsp.query.seq)
        self.assertEqual(
            "6999***********************************99", hsp.aln_annotation["PP"]
        )

        hsp = hit.hsps[1]
        self.assertEqual(2, hsp.domain_index)
        self.assertFalse(hsp.is_included)
        self.assertEqual(0.8, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(0.19, hsp.evalue_cond)
        self.assertEqual(5.2e02, hsp.evalue)
        self.assertEqual(38, hsp.hit_start)
        self.assertEqual(62, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)
        self.assertEqual(244, hsp.query_start)
        self.assertEqual(268, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual(242, hsp.env_start)
        self.assertEqual(270, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.86, hsp.acc_avg)
        self.assertEqual("rgrpsaavlaalaralgldpaera", hsp.hit.seq)
        self.assertEqual("CPKPSLQQITHIANQLGLEKDVVR", hsp.query.seq)
        self.assertEqual("678**************9988765", hsp.aln_annotation["PP"])

        hit = qresult[3]
        self.assertEqual("Homeobox_KN", hit.id)
        self.assertEqual("Homeobox KN domain", hit.description)
        self.assertFalse(hit.is_included)
        self.assertEqual(0.039, hit.evalue)
        self.assertEqual(13.5, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual(1.6, hit.domain_exp_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        hsp = hit.hsps[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertFalse(hsp.is_included)
        self.assertEqual(12.3, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(3.5e-05, hsp.evalue_cond)
        self.assertEqual(0.095, hsp.evalue)
        self.assertEqual(6, hsp.hit_start)
        self.assertEqual(39, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)
        self.assertEqual(243, hsp.query_start)
        self.assertEqual(276, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual(240, hsp.env_start)
        self.assertEqual(277, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.91, hsp.acc_avg)
        self.assertEqual("hnPYPskevkeelakqTglsrkqidnWFiNaRr", hsp.hit.seq)
        self.assertEqual("KCPKPSLQQITHIANQLGLEKDVVRVWFCNRRQ", hsp.query.seq)
        self.assertEqual("56779*************************996", hsp.aln_annotation["PP"])

        hit = qresult[4]
        self.assertEqual("DUF521", hit.id)
        self.assertEqual("Protein of unknown function (DUF521)", hit.description)
        self.assertFalse(hit.is_included)
        self.assertEqual(0.14, hit.evalue)
        self.assertEqual(10.5, hit.bitscore)
        self.assertEqual(0.1, hit.bias)
        self.assertEqual(1.4, hit.domain_exp_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        hsp = hit.hsps[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertFalse(hsp.is_included)
        self.assertEqual(9.6, hsp.bitscore)
        self.assertEqual(0.1, hsp.bias)
        self.assertEqual(9.4e-05, hsp.evalue_cond)
        self.assertEqual(0.26, hsp.evalue)
        self.assertEqual(272, hsp.hit_start)
        self.assertEqual(334, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)
        self.assertEqual(220, hsp.query_start)
        self.assertEqual(280, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual(196, hsp.env_start)
        self.assertEqual(294, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.77, hsp.acc_avg)
        self.assertEqual(
            "adlaavleelnkakkeevdlvvlGcPhlsleeleelaellkgrkkkvsvelvvttsravlsk",
            hsp.hit.seq,
        )
        self.assertEqual(
            "QARKRKRTSIENRVRWSLETMFLKCPKPSLQQITHIANQLGLEK--DVVRVWFCNRRQKGKR",
            hsp.query.seq,
        )
        self.assertEqual(
            "345666667778888899************************99..9999999988876554",
            hsp.aln_annotation["PP"],
        )

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(5, counter)

    def test_30_hmmscan_002(self):
        """Parsing hmmscan 3.0 (text_30_hmmscan_002)."""
        txt_file = get_file("text_30_hmmscan_002.out")
        qresults = parse(txt_file, FMT)
        counter = 0

        # test first qresult
        qresult = next(qresults)
        counter += 1

        self.assertEqual("hmmscan", qresult.program)
        self.assertEqual("/home/bow/db/hmmer/Pfam-A.hmm", qresult.target)
        self.assertEqual("3.0", qresult.version)
        self.assertEqual("random_s00", qresult.id)
        self.assertEqual(32, qresult.seq_len)
        self.assertEqual(0, len(qresult))

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(1, counter)

    def test_30_hmmscan_003(self):
        """Parsing hmmscan 3.0 (text_30_hmmscan_003)."""
        txt_file = get_file("text_30_hmmscan_003.out")
        qresults = parse(txt_file, FMT)
        counter = 0

        # test first qresult
        qresult = next(qresults)
        counter += 1

        self.assertEqual("hmmscan", qresult.program)
        self.assertEqual("/home/bow/db/hmmer/Pfam-A.hmm", qresult.target)
        self.assertEqual("3.0", qresult.version)
        self.assertEqual("gi|4885477|ref|NP_005359.1|", qresult.id)
        self.assertEqual("myoglobin [Homo sapiens]", qresult.description)
        self.assertEqual(154, qresult.seq_len)
        self.assertEqual(1, len(qresult))

        hit = qresult[0]
        self.assertEqual("Globin", hit.id)
        self.assertEqual("Globin", hit.description)
        self.assertTrue(hit.is_included)
        self.assertEqual(6e-21, hit.evalue)
        self.assertEqual(74.6, hit.bitscore)
        self.assertEqual(0.3, hit.bias)
        self.assertEqual(1.3, hit.domain_exp_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        hsp = hit.hsps[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_included)
        self.assertEqual(74.0, hsp.bitscore)
        self.assertEqual(0.2, hsp.bias)
        self.assertEqual(6.7e-25, hsp.evalue_cond)
        self.assertEqual(9.2e-21, hsp.evalue)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(107, hsp.hit_end)
        self.assertEqual("[.", hsp.hit_endtype)
        self.assertEqual(6, hsp.query_start)
        self.assertEqual(112, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual(6, hsp.env_start)
        self.assertEqual(113, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.97, hsp.acc_avg)
        self.assertEqual(
            "HHHHHHHHHHHHCHHHHHHHHHHHHHHHHHSGGGGGGGCCCTTTT.HHHHHTSCHHHHHHHHHHHHHHHHHHCTTSHHHHHHHHHHHHHHHHTT-.--HHHHCCHHHHH",
            hsp.aln_annotation["CS"],
        )
        self.assertEqual(
            "qkalvkaswekvkanaeeigaeilkrlfkaypdtkklFkkfgdls.aedlksspkfkahakkvlaaldeavknldnddnlkaalkklgarHakrg.vdpanfklfgeal",
            hsp.hit.seq,
        )
        self.assertEqual(
            "EWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKsEDEMKASEDLKKHGATVLTALGGILKK---KGHHEAEIKPLAQSHATKHkIPVKYLEFISECI",
            hsp.query.seq,
        )
        self.assertEqual(
            "5789*********************************************************************...6899***********************999998",
            hsp.aln_annotation["PP"],
        )

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(1, counter)

    def test_30_hmmscan_004(self):
        """Parsing hmmscan 3.0 (text_30_hmmscan_004)."""
        txt_file = get_file("text_30_hmmscan_004.out")
        qresults = parse(txt_file, FMT)
        counter = 0

        # test first qresult
        qresult = next(qresults)
        counter += 1

        self.assertEqual("hmmscan", qresult.program)
        self.assertEqual("/home/bow/db/hmmer/Pfam-A.hmm", qresult.target)
        self.assertEqual("3.0", qresult.version)
        self.assertEqual("gi|126362951:116-221", qresult.id)
        self.assertEqual(
            "leukocyte immunoglobulin-like receptor subfamily B member 1 isoform 2 precursor [Homo sapiens]",
            qresult.description,
        )
        self.assertEqual(106, qresult.seq_len)
        self.assertEqual(2, len(qresult))

        hit = qresult[0]
        self.assertEqual("Ig_3", hit.id)
        self.assertEqual("Immunoglobulin domain", hit.description)
        self.assertTrue(hit.is_included)
        self.assertEqual(1.4e-09, hit.evalue)
        self.assertEqual(38.2, hit.bitscore)
        self.assertEqual(0.4, hit.bias)
        self.assertEqual(1.3, hit.domain_exp_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        hsp = hit.hsps[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_included)
        self.assertEqual(37.6, hsp.bitscore)
        self.assertEqual(0.3, hsp.bias)
        self.assertEqual(3e-13, hsp.evalue_cond)
        self.assertEqual(2.1e-09, hsp.evalue)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(73, hsp.hit_end)
        self.assertEqual("[.", hsp.hit_endtype)
        self.assertEqual(8, hsp.query_start)
        self.assertEqual(84, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual(8, hsp.env_start)
        self.assertEqual(88, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.94, hsp.acc_avg)
        self.assertEqual(
            "kPvisvspsptvtsggnvtLtCsaeggpppptisWy.....ietppelqgsegssssestLtissvtsedsgtYtCva",
            hsp.hit.seq,
        )
        self.assertEqual(
            "KPTLSAQPSPVVNSGGNVILQCDSQVA--FDGFSLCkegedEHPQCLNSQPHARGSSRAIFSVGPVSPSRRWWYRCYA",
            hsp.query.seq,
        )
        self.assertEqual(
            "8************************99..78888888****************************************9",
            hsp.aln_annotation["PP"],
        )

        hit = qresult[-1]
        self.assertEqual("Ig_2", hit.id)
        self.assertEqual("Immunoglobulin domain", hit.description)
        self.assertEqual(3.5e-05, hit.evalue)
        self.assertEqual(23.7, hit.bitscore)
        self.assertEqual(0.1, hit.bias)
        self.assertEqual(1.1, hit.domain_exp_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        hsp = hit.hsps[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_included)
        self.assertEqual(23.4, hsp.bitscore)
        self.assertEqual(0.1, hsp.bias)
        self.assertEqual(6.2e-09, hsp.evalue_cond)
        self.assertEqual(4.3e-05, hsp.evalue)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(80, hsp.hit_end)
        self.assertEqual("[]", hsp.hit_endtype)
        self.assertEqual(8, hsp.query_start)
        self.assertEqual(104, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual(8, hsp.env_start)
        self.assertEqual(104, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.71, hsp.acc_avg)
        self.assertEqual(
            "kpvlvapp.svvtegenvtLtCsapgnptprvqwykdg.vels......qsqnq........lfipnvsaedsgtYtCra....rnseggktstsveltv",
            hsp.hit.seq,
        )
        self.assertEqual(
            "KPTLSAQPsPVVNSGGNVILQCDSQVA-FDGFSLCKEGeDEHPqclnsqP---HargssraiFSVGPVSPSRRWWYRCYAydsnSPYEWSLPSDLLELLV",
            hsp.query.seq,
        )
        self.assertEqual(
            "799998885779*************85.899***9988655554443320...134455543444669************88443344588888888766",
            hsp.aln_annotation["PP"],
        )

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(1, counter)

    def test_30_hmmscan_005(self):
        """Parsing hmmscan 3.0 (text_30_hmmscan_005)."""
        txt_file = get_file("text_30_hmmscan_005.out")
        qresults = parse(txt_file, FMT)
        counter = 0

        # test first result
        qresult = next(qresults)
        counter += 1

        self.assertEqual("hmmscan", qresult.program)
        self.assertEqual("/home/bow/db/hmmer/Pfam-A.hmm", qresult.target)
        self.assertEqual("3.0", qresult.version)
        self.assertEqual("gi|22748937|ref|NP_065801.1|", qresult.id)
        self.assertEqual("exportin-5 [Homo sapiens]", qresult.description)
        self.assertEqual(1204, qresult.seq_len)
        self.assertEqual(2, len(qresult))

        hit = qresult[0]
        self.assertEqual("Xpo1", hit.id)
        self.assertEqual("Exportin 1-like protein", hit.description)
        self.assertTrue(hit.is_included)
        self.assertEqual(7.8e-34, hit.evalue)
        self.assertEqual(116.6, hit.bitscore)
        self.assertEqual(7.8, hit.bias)
        self.assertEqual(2.8, hit.domain_exp_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))

        hsp = hit.hsps[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_included)
        self.assertEqual(116.1, hsp.bitscore)
        self.assertEqual(3.4, hsp.bias)
        self.assertEqual(1.6e-37, hsp.evalue_cond)
        self.assertEqual(1.1e-33, hsp.evalue)
        self.assertEqual(1, hsp.hit_start)
        self.assertEqual(148, hsp.hit_end)
        self.assertEqual(".]", hsp.hit_endtype)
        self.assertEqual(109, hsp.query_start)
        self.assertEqual(271, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual(108, hsp.env_start)
        self.assertEqual(271, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.98, hsp.acc_avg)
        self.assertEqual(
            "HHHHHHHHHHHHHHHHHHTTTTSTTHHHHHHHHHHG-HHHHHHHHHHHHHHHHHHCCS-TTTS-CCCHHHHHHHCHHHHHHHHHHHHHHHC-TT-..................HHHHHHHHHHHHHHCTTS-CHHCHCS...HHHHHCHHCCSCCCHHHHHHHH",
            hsp.aln_annotation["CS"],
        )
        self.assertEqual(
            "kflrnklaealaelflqeypnqWpsffddllsllssspsglelllriLkvlpeEiadfsrskleqerrnelkdllrsqvqkilelllqileqsvskk...............sselveatLkclsswvswidiglivnsp..llsllfqlLndpelreaAvecL",
            hsp.hit.seq,
        )
        self.assertEqual(
            "NHIKDALSRIVVEMIKREWPQHWPDMLIELDTLSKQGETQTELVMFILLRLAEDVVTF--QTLPPQRRRDIQQTLTQNMERIFSFLLNTLQENVNKYqqvktdtsqeskaqaNCRVGVAALNTLAGYIDWVSMSHITAENckLLEILCLLLNEQELQLGAAECL",
            hsp.query.seq,
        )
        self.assertEqual(
            "89******************************************************99..79*********************************99*****************************************8889*********************8",
            hsp.aln_annotation["PP"],
        )

        hsp = hit.hsps[-1]
        self.assertEqual(2, hsp.domain_index)
        self.assertFalse(hsp.is_included)
        self.assertEqual(-1.8, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(0.35, hsp.evalue_cond)
        self.assertEqual(2.4e03, hsp.evalue)
        self.assertEqual(111, hsp.hit_start)
        self.assertEqual(139, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)
        self.assertEqual(498, hsp.query_start)
        self.assertEqual(525, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual(495, hsp.env_start)
        self.assertEqual(529, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.86, hsp.acc_avg)
        self.assertEqual("HHCTTS-CHHCHCS.HHHHHCHHCCSCC", hsp.aln_annotation["CS"])
        self.assertEqual("swvswidiglivnspllsllfqlLndpe", hsp.hit.seq)
        self.assertEqual("SFVQWEAMTLFLES-VITQMFRTLNREE", hsp.query.seq)
        self.assertEqual("899*********98.8888899998776", hsp.aln_annotation["PP"])

        hit = qresult[-1]
        self.assertEqual("IBN_N", hit.id)
        self.assertEqual("Importin-beta N-terminal domain", hit.description)
        self.assertTrue(hit.is_included)
        self.assertEqual(0.0039, hit.evalue)
        self.assertEqual(16.9, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual(2.7, hit.domain_exp_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))

        hsp = hit.hsps[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_included)
        self.assertEqual(14.0, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(4.8e-06, hsp.evalue_cond)
        self.assertEqual(0.033, hsp.evalue)
        self.assertEqual(3, hsp.hit_start)
        self.assertEqual(75, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)
        self.assertEqual(35, hsp.query_start)
        self.assertEqual(98, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual(32, hsp.env_start)
        self.assertEqual(100, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.87, hsp.acc_avg)
        self.assertEqual(
            "HHHHHHHSCTHHHHHHHHHHHTTTSTHHHHHHHHHHHHHHHHHSCCHHHHHHHHCS-HHHHHHHHHHHHHHH",
            hsp.aln_annotation["CS"],
        )
        self.assertEqual(
            "qLnqlekqkPgflsallqilanksldlevRqlAalyLknlItkhWkseeaqrqqqlpeeekelIrnnllnll",
            hsp.hit.seq,
        )
        self.assertEqual(
            "FCEEFKEKCPICVPCGLRLA-EKTQVAIVRHFGLQILEHVVKFRWN--------GMSRLEKVYLKNSVMELI",
            hsp.query.seq,
        )
        self.assertEqual(
            "56788886699*********.6555899******************........999999****99999887",
            hsp.aln_annotation["PP"],
        )

        hsp = hit.hsps[-1]
        self.assertEqual(2, hsp.domain_index)
        self.assertFalse(hsp.is_included)
        self.assertEqual(-3.3, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(1.2, hsp.evalue_cond)
        self.assertEqual(8e03, hsp.evalue)
        self.assertEqual(56, hsp.hit_start)
        self.assertEqual(75, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)
        self.assertEqual(167, hsp.query_start)
        self.assertEqual(186, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual(164, hsp.env_start)
        self.assertEqual(187, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.85, hsp.acc_avg)
        self.assertEqual("HCS-HHHHHHHHHHHHHHH", hsp.aln_annotation["CS"])
        self.assertEqual("qqlpeeekelIrnnllnll", hsp.hit.seq)
        self.assertEqual("QTLPPQRRRDIQQTLTQNM", hsp.query.seq)
        self.assertEqual("6899*******99998865", hsp.aln_annotation["PP"])

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(1, counter)

    def test_30_hmmscan_006(self):
        """Parsing hmmscan 3.0 (text_30_hmmscan_006)."""
        txt_file = get_file("text_30_hmmscan_006.out")
        qresults = parse(txt_file, FMT)
        counter = 0

        # test first qresult
        qresult = next(qresults)
        counter += 1

        self.assertEqual("hmmscan", qresult.program)
        self.assertEqual("/home/bow/db/hmmer/Pfam-A.hmm", qresult.target)
        self.assertEqual("3.0", qresult.version)
        self.assertEqual("gi|125490392|ref|NP_038661.2|", qresult.id)
        self.assertEqual(
            "POU domain, class 5, transcription factor 1 isoform 1 [Mus musculus]",
            qresult.description,
        )
        self.assertEqual(352, qresult.seq_len)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual("Pou", hit.id)
        self.assertEqual("Pou domain - N-terminal to homeobox domain", hit.description)
        self.assertTrue(hit.is_included)
        self.assertEqual(7e-37, hit.evalue)
        self.assertEqual(124.8, hit.bitscore)
        self.assertEqual(0.5, hit.bias)
        self.assertEqual(1.5, hit.domain_exp_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        hsp = hit.hsps[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_included)
        self.assertEqual(123.9, hsp.bitscore)
        self.assertEqual(0.3, hsp.bias)
        self.assertEqual(5e-40, hsp.evalue_cond)
        self.assertEqual(1.4e-36, hsp.evalue)
        self.assertEqual(2, hsp.hit_start)
        self.assertEqual(75, hsp.hit_end)
        self.assertEqual(".]", hsp.hit_endtype)
        self.assertEqual(132, hsp.query_start)
        self.assertEqual(205, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual(130, hsp.env_start)
        self.assertEqual(205, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.97, hsp.acc_avg)
        self.assertEqual(
            "eldleeleefakefkqrrikLgltqadvgsalgalyGkefsqttIcrFEalqLslknmckLkpllekWLeeae",
            hsp.hit.seq,
        )
        self.assertEqual(
            "KALQKELEQFAKLLKQKRITLGYTQADVGLTLGVLFGKVFSQTTICRFEALQLSLKNMCKLRPLLEKWVEEAD",
            hsp.query.seq,
        )
        self.assertEqual(
            "67899******************************************************************96",
            hsp.aln_annotation["PP"],
        )

        hit = qresult[1]
        self.assertEqual("Homeobox", hit.id)
        self.assertEqual("Homeobox domain", hit.description)
        self.assertTrue(hit.is_included)
        self.assertEqual(2.1e-18, hit.evalue)
        self.assertEqual(65.5, hit.bitscore)
        self.assertEqual(1.1, hit.bias)
        self.assertEqual(1.5, hit.domain_exp_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        hsp = hit.hsps[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_included)
        self.assertEqual(64.6, hsp.bitscore)
        self.assertEqual(0.7, hsp.bias)
        self.assertEqual(1.5e-21, hsp.evalue_cond)
        self.assertEqual(4.1e-18, hsp.evalue)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(57, hsp.hit_end)
        self.assertEqual("[]", hsp.hit_endtype)
        self.assertEqual(223, hsp.query_start)
        self.assertEqual(280, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual(223, hsp.env_start)
        self.assertEqual(280, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.98, hsp.acc_avg)
        self.assertEqual(
            "SS--SS--HHHHHHHHHHCCTSSS--HHHHHHHHHH----HHHHHHHHHHHHHHHHH",
            hsp.aln_annotation["CS"],
        )
        self.assertEqual(
            "rrkRttftkeqleeLeelFeknrypsaeereeLAkklgLterqVkvWFqNrRakekk", hsp.hit.seq
        )
        self.assertEqual(
            "KRKRTSIENRVRWSLETMFLKCPKPSLQQITHIANQLGLEKDVVRVWFCNRRQKGKR", hsp.query.seq
        )
        self.assertEqual(
            "79****************************************************997",
            hsp.aln_annotation["PP"],
        )

        hit = qresult[2]
        self.assertEqual("HTH_31", hit.id)
        self.assertEqual("Helix-turn-helix domain", hit.description)
        self.assertFalse(hit.is_included)
        self.assertEqual(0.012, hit.evalue)
        self.assertEqual(15.6, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual(2.2, hit.domain_exp_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))

        hsp = hit.hsps[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertFalse(hsp.is_included)
        self.assertEqual(12.0, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(5.7e-05, hsp.evalue_cond)
        self.assertEqual(0.16, hsp.evalue)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(35, hsp.hit_end)
        self.assertEqual("[.", hsp.hit_endtype)
        self.assertEqual(140, hsp.query_start)
        self.assertEqual(181, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual(140, hsp.env_start)
        self.assertEqual(184, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.96, hsp.acc_avg)
        self.assertEqual("aLGarLralReraGLtqeevAerlg......vSastlsrlE", hsp.hit.seq)
        self.assertEqual("QFAKLLKQKRITLGYTQADVGLTLGvlfgkvFSQTTICRFE", hsp.query.seq)
        self.assertEqual(
            "6999***********************************99", hsp.aln_annotation["PP"]
        )

        hsp = hit.hsps[1]
        self.assertEqual(2, hsp.domain_index)
        self.assertFalse(hsp.is_included)
        self.assertEqual(0.8, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(0.19, hsp.evalue_cond)
        self.assertEqual(5.2e02, hsp.evalue)
        self.assertEqual(38, hsp.hit_start)
        self.assertEqual(62, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)
        self.assertEqual(244, hsp.query_start)
        self.assertEqual(268, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual(242, hsp.env_start)
        self.assertEqual(270, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.86, hsp.acc_avg)
        self.assertEqual("rgrpsaavlaalaralgldpaera", hsp.hit.seq)
        self.assertEqual("CPKPSLQQITHIANQLGLEKDVVR", hsp.query.seq)
        self.assertEqual("678**************9988765", hsp.aln_annotation["PP"])

        hit = qresult[3]
        self.assertEqual("Homeobox_KN", hit.id)
        self.assertEqual("Homeobox KN domain", hit.description)
        self.assertFalse(hit.is_included)
        self.assertEqual(0.039, hit.evalue)
        self.assertEqual(13.5, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual(1.6, hit.domain_exp_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        hsp = hit.hsps[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertFalse(hsp.is_included)
        self.assertEqual(12.3, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(3.5e-05, hsp.evalue_cond)
        self.assertEqual(0.095, hsp.evalue)
        self.assertEqual(6, hsp.hit_start)
        self.assertEqual(39, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)
        self.assertEqual(243, hsp.query_start)
        self.assertEqual(276, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual(240, hsp.env_start)
        self.assertEqual(277, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.91, hsp.acc_avg)
        self.assertEqual("hnPYPskevkeelakqTglsrkqidnWFiNaRr", hsp.hit.seq)
        self.assertEqual("KCPKPSLQQITHIANQLGLEKDVVRVWFCNRRQ", hsp.query.seq)
        self.assertEqual("56779*************************996", hsp.aln_annotation["PP"])

        hit = qresult[4]
        self.assertEqual("DUF521", hit.id)
        self.assertEqual("Protein of unknown function (DUF521)", hit.description)
        self.assertFalse(hit.is_included)
        self.assertEqual(0.14, hit.evalue)
        self.assertEqual(10.5, hit.bitscore)
        self.assertEqual(0.1, hit.bias)
        self.assertEqual(1.4, hit.domain_exp_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        hsp = hit.hsps[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertFalse(hsp.is_included)
        self.assertEqual(9.6, hsp.bitscore)
        self.assertEqual(0.1, hsp.bias)
        self.assertEqual(9.4e-05, hsp.evalue_cond)
        self.assertEqual(0.26, hsp.evalue)
        self.assertEqual(272, hsp.hit_start)
        self.assertEqual(334, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)
        self.assertEqual(220, hsp.query_start)
        self.assertEqual(280, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual(196, hsp.env_start)
        self.assertEqual(294, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.77, hsp.acc_avg)
        self.assertEqual(
            "adlaavleelnkakkeevdlvvlGcPhlsleeleelaellkgrkkkvsvelvvttsravlsk",
            hsp.hit.seq,
        )
        self.assertEqual(
            "QARKRKRTSIENRVRWSLETMFLKCPKPSLQQITHIANQLGLEK--DVVRVWFCNRRQKGKR",
            hsp.query.seq,
        )
        self.assertEqual(
            "345666667778888899************************99..9999999988876554",
            hsp.aln_annotation["PP"],
        )

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(1, counter)

    def test_30_hmmscan_007(self):
        """Parsing hmmscan 3.0 (text_30_hmmscan_007)."""
        txt_file = get_file("text_30_hmmscan_007.out")
        qresults = parse(txt_file, FMT)
        counter = 0

        # test first qresult
        qresult = next(qresults)
        counter += 1

        self.assertEqual("hmmscan", qresult.program)
        self.assertEqual("/home/bow/db/hmmer/Pfam-A.hmm", qresult.target)
        self.assertEqual("3.0", qresult.version)
        self.assertEqual("gi|125490392|ref|NP_038661.2|", qresult.id)
        self.assertEqual(
            "POU domain, class 5, transcription factor 1 isoform 1 [Mus musculus]",
            qresult.description,
        )
        self.assertEqual(352, qresult.seq_len)
        self.assertEqual(5, len(qresult))

        hit = qresult[0]
        self.assertEqual("Pou", hit.id)
        self.assertEqual("Pou domain - N-terminal to homeobox domain", hit.description)
        self.assertTrue(hit.is_included)
        self.assertEqual(7e-37, hit.evalue)
        self.assertEqual(124.8, hit.bitscore)
        self.assertEqual(0.5, hit.bias)
        self.assertEqual(1.5, hit.domain_exp_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        hsp = hit.hsps[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_included)
        self.assertEqual(123.9, hsp.bitscore)
        self.assertEqual(0.3, hsp.bias)
        self.assertEqual(5e-40, hsp.evalue_cond)
        self.assertEqual(1.4e-36, hsp.evalue)
        self.assertEqual(2, hsp.hit_start)
        self.assertEqual(75, hsp.hit_end)
        self.assertEqual(".]", hsp.hit_endtype)
        self.assertEqual(132, hsp.query_start)
        self.assertEqual(205, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual(130, hsp.env_start)
        self.assertEqual(205, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.97, hsp.acc_avg)

        hit = qresult[1]
        self.assertEqual("Homeobox", hit.id)
        self.assertEqual("Homeobox domain", hit.description)
        self.assertTrue(hit.is_included)
        self.assertEqual(2.1e-18, hit.evalue)
        self.assertEqual(65.5, hit.bitscore)
        self.assertEqual(1.1, hit.bias)
        self.assertEqual(1.5, hit.domain_exp_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        hsp = hit.hsps[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_included)
        self.assertEqual(64.6, hsp.bitscore)
        self.assertEqual(0.7, hsp.bias)
        self.assertEqual(1.5e-21, hsp.evalue_cond)
        self.assertEqual(4.1e-18, hsp.evalue)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(57, hsp.hit_end)
        self.assertEqual("[]", hsp.hit_endtype)
        self.assertEqual(223, hsp.query_start)
        self.assertEqual(280, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual(223, hsp.env_start)
        self.assertEqual(280, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.98, hsp.acc_avg)

        hit = qresult[2]
        self.assertEqual("HTH_31", hit.id)
        self.assertEqual("Helix-turn-helix domain", hit.description)
        self.assertFalse(hit.is_included)
        self.assertEqual(0.012, hit.evalue)
        self.assertEqual(15.6, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual(2.2, hit.domain_exp_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))

        hsp = hit.hsps[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertFalse(hsp.is_included)
        self.assertEqual(12.0, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(5.7e-05, hsp.evalue_cond)
        self.assertEqual(0.16, hsp.evalue)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(35, hsp.hit_end)
        self.assertEqual("[.", hsp.hit_endtype)
        self.assertEqual(140, hsp.query_start)
        self.assertEqual(181, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual(140, hsp.env_start)
        self.assertEqual(184, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.96, hsp.acc_avg)

        hsp = hit.hsps[1]
        self.assertEqual(2, hsp.domain_index)
        self.assertFalse(hsp.is_included)
        self.assertEqual(0.8, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(0.19, hsp.evalue_cond)
        self.assertEqual(5.2e02, hsp.evalue)
        self.assertEqual(38, hsp.hit_start)
        self.assertEqual(62, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)
        self.assertEqual(244, hsp.query_start)
        self.assertEqual(268, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual(242, hsp.env_start)
        self.assertEqual(270, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.86, hsp.acc_avg)

        hit = qresult[3]
        self.assertEqual("Homeobox_KN", hit.id)
        self.assertEqual("Homeobox KN domain", hit.description)
        self.assertFalse(hit.is_included)
        self.assertEqual(0.039, hit.evalue)
        self.assertEqual(13.5, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual(1.6, hit.domain_exp_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        hsp = hit.hsps[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertFalse(hsp.is_included)
        self.assertEqual(12.3, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(3.5e-05, hsp.evalue_cond)
        self.assertEqual(0.095, hsp.evalue)
        self.assertEqual(6, hsp.hit_start)
        self.assertEqual(39, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)
        self.assertEqual(243, hsp.query_start)
        self.assertEqual(276, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual(240, hsp.env_start)
        self.assertEqual(277, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.91, hsp.acc_avg)

        hit = qresult[4]
        self.assertEqual("DUF521", hit.id)
        self.assertEqual("Protein of unknown function (DUF521)", hit.description)
        self.assertFalse(hit.is_included)
        self.assertEqual(0.14, hit.evalue)
        self.assertEqual(10.5, hit.bitscore)
        self.assertEqual(0.1, hit.bias)
        self.assertEqual(1.4, hit.domain_exp_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        hsp = hit.hsps[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertFalse(hsp.is_included)
        self.assertEqual(9.6, hsp.bitscore)
        self.assertEqual(0.1, hsp.bias)
        self.assertEqual(9.4e-05, hsp.evalue_cond)
        self.assertEqual(0.26, hsp.evalue)
        self.assertEqual(272, hsp.hit_start)
        self.assertEqual(334, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)
        self.assertEqual(220, hsp.query_start)
        self.assertEqual(280, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual(196, hsp.env_start)
        self.assertEqual(294, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.77, hsp.acc_avg)

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(1, counter)

    def test_30_hmmscan_008(self):
        """Parsing hmmscan 3.0 (text_30_hmmscan_008)."""
        txt_file = get_file("text_30_hmmscan_008.out")
        qresults = parse(txt_file, FMT)
        counter = 0

        # test first result
        qresult = next(qresults)
        counter += 1

        self.assertEqual("hmmscan", qresult.program)
        self.assertEqual("/home/bow/db/hmmer/Pfam-A.hmm", qresult.target)
        self.assertEqual("3.0", qresult.version)
        self.assertEqual("gi|22748937|ref|NP_065801.1|", qresult.id)
        self.assertEqual("exportin-5 [Homo sapiens]", qresult.description)
        self.assertEqual(1204, qresult.seq_len)
        self.assertEqual(2, len(qresult))

        hit = qresult[0]
        self.assertEqual("Xpo1", hit.id)
        self.assertEqual("Exportin 1-like protein", hit.description)
        self.assertTrue(hit.is_included)
        self.assertEqual(7.8e-34, hit.evalue)
        self.assertEqual(116.6, hit.bitscore)
        self.assertEqual(7.8, hit.bias)
        self.assertEqual(2.8, hit.domain_exp_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))

        hsp = hit.hsps[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_included)
        self.assertEqual(116.1, hsp.bitscore)
        self.assertEqual(3.4, hsp.bias)
        self.assertEqual(1.6e-37, hsp.evalue_cond)
        self.assertEqual(1.1e-33, hsp.evalue)
        self.assertEqual(1, hsp.hit_start)
        self.assertEqual(148, hsp.hit_end)
        self.assertEqual(".]", hsp.hit_endtype)
        self.assertEqual(109, hsp.query_start)
        self.assertEqual(271, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual(108, hsp.env_start)
        self.assertEqual(271, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.98, hsp.acc_avg)
        self.assertEqual(
            "HHHHHHHHHHHHHHHHHHTTTTSTTHHHHHHHHHHG-HHHHHHHHHHHHHHHHHHCCS-TTTS-CCCHHHHHHHCHHHHHHHHHHHHHHHC-TT-..................HHHHHHHHHHHHHHCTTS-CHHCHCS...HHHHHCHHCCSCCCHHHHHHHH",
            hsp.aln_annotation["CS"],
        )
        self.assertEqual(
            "kflrnklaealaelflqeypnqWpsffddllsllssspsglelllriLkvlpeEiadfsrskleqerrnelkdllrsqvqkilelllqileqsvskk...............sselveatLkclsswvswidiglivnsp..llsllfqlLndpelreaAvecL",
            hsp.hit.seq,
        )
        self.assertEqual(
            "NHIKDALSRIVVEMIKREWPQHWPDMLIELDTLSKQGETQTELVMFILLRLAEDVVTF--QTLPPQRRRDIQQTLTQNMERIFSFLLNTLQENVNKYqqvktdtsqeskaqaNCRVGVAALNTLAGYIDWVSMSHITAENckLLEILCLLLNEQELQLGAAECL",
            hsp.query.seq,
        )
        self.assertEqual(
            "89******************************************************99..79*********************************99*****************************************8889*********************8",
            hsp.aln_annotation["PP"],
        )

        hsp = hit.hsps[-1]
        self.assertEqual(2, hsp.domain_index)
        self.assertFalse(hsp.is_included)
        self.assertEqual(-1.8, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(0.35, hsp.evalue_cond)
        self.assertEqual(2.4e03, hsp.evalue)
        self.assertEqual(111, hsp.hit_start)
        self.assertEqual(139, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)
        self.assertEqual(498, hsp.query_start)
        self.assertEqual(525, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual(495, hsp.env_start)
        self.assertEqual(529, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.86, hsp.acc_avg)
        self.assertEqual("HHCTTS-CHHCHCS.HHHHHCHHCCSCC", hsp.aln_annotation["CS"])
        self.assertEqual("swvswidiglivnspllsllfqlLndpe", hsp.hit.seq)
        self.assertEqual("SFVQWEAMTLFLES-VITQMFRTLNREE", hsp.query.seq)
        self.assertEqual("899*********98.8888899998776", hsp.aln_annotation["PP"])

        hit = qresult[-1]
        self.assertEqual("IBN_N", hit.id)
        self.assertEqual("Importin-beta N-terminal domain", hit.description)
        self.assertTrue(hit.is_included)
        self.assertEqual(0.0039, hit.evalue)
        self.assertEqual(16.9, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual(2.7, hit.domain_exp_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))

        hsp = hit.hsps[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_included)
        self.assertEqual(14.0, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(4.8e-06, hsp.evalue_cond)
        self.assertEqual(0.033, hsp.evalue)
        self.assertEqual(3, hsp.hit_start)
        self.assertEqual(75, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)
        self.assertEqual(35, hsp.query_start)
        self.assertEqual(98, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual(32, hsp.env_start)
        self.assertEqual(100, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.87, hsp.acc_avg)
        self.assertEqual(
            "HHHHHHHSCTHHHHHHHHHHHTTTSTHHHHHHHHHHHHHHHHHSCCHHHHHHHHCS-HHHHHHHHHHHHHHH",
            hsp.aln_annotation["CS"],
        )
        self.assertEqual(
            "qLnqlekqkPgflsallqilanksldlevRqlAalyLknlItkhWkseeaqrqqqlpeeekelIrnnllnll",
            hsp.hit.seq,
        )
        self.assertEqual(
            "FCEEFKEKCPICVPCGLRLA-EKTQVAIVRHFGLQILEHVVKFRWN--------GMSRLEKVYLKNSVMELI",
            hsp.query.seq,
        )
        self.assertEqual(
            "56788886699*********.6555899******************........999999****99999887",
            hsp.aln_annotation["PP"],
        )

        hsp = hit.hsps[-1]
        self.assertEqual(2, hsp.domain_index)
        self.assertFalse(hsp.is_included)
        self.assertEqual(-3.3, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(1.2, hsp.evalue_cond)
        self.assertEqual(8e03, hsp.evalue)
        self.assertEqual(56, hsp.hit_start)
        self.assertEqual(75, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)
        self.assertEqual(167, hsp.query_start)
        self.assertEqual(186, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual(164, hsp.env_start)
        self.assertEqual(187, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.85, hsp.acc_avg)
        self.assertEqual("HCS-HHHHHHHHHHHHHHH", hsp.aln_annotation["CS"])
        self.assertEqual("qqlpeeekelIrnnllnll", hsp.hit.seq)
        self.assertEqual("QTLPPQRRRDIQQTLTQNM", hsp.query.seq)
        self.assertEqual("6899*******99998865", hsp.aln_annotation["PP"])

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(1, counter)

    def test_30_hmmscan_009(self):
        """Parsing hmmscan 3.0 (text_30_hmmscan_009)."""
        hmmer_file = get_file("text_30_hmmscan_009.out")
        qresults = parse(hmmer_file, FMT)
        counter = 0

        # test first qresult
        qresult = next(qresults)
        counter += 1
        self.assertEqual("SCO3574", qresult.id)
        self.assertEqual(5, len(qresult.hits))
        self.assertEqual("Esterase", qresult.hits[3].id)

    def test_30_hmmscan_010(self):
        """Parsing hmmscan 3.0 (text_30_hmmscan_010)."""
        hmmer_file = get_file("text_30_hmmscan_010.out")
        qresults = list(parse(hmmer_file, FMT))

        # test the Hit object without HSPs
        hit = qresults[0][-1]
        self.assertFalse(hit)
        self.assertEqual("NRPS-COM_Cterm", hit.id)
        self.assertEqual("", hit.description)
        self.assertEqual("bpsA", hit.query_id)
        self.assertEqual("<unknown description>", hit.query_description)
        self.assertEqual(4.4e-11, hit.evalue)
        self.assertEqual(33.6, hit.bitscore)
        self.assertEqual(10.2, hit.bias)
        self.assertEqual(2.9, hit.domain_exp_num)
        self.assertEqual(0, hit.domain_obs_num)
        self.assertEqual(0, len(hit))


class HmmersearchCases(unittest.TestCase):
    """Test hmmsearch output."""

    def test_31b2_hmmsearch_001(self):
        """Test parsing hmmsearch 3.1b2 (text_31b2_hmmsearch_001)."""
        txt_file = get_file("text_31b2_hmmsearch_001.out")
        qresults = list(parse(txt_file, FMT))

        # Assert we deal with only 1 query
        self.assertEqual(len(qresults), 1)

        # Test whether proper query id is read
        self.assertEqual(qresults[0].id, "infile_sto")

        # Test if proper number of hits is read
        self.assertEqual(len(qresults[0].hits), 10)

    def test_31b2_hmmsearch_002(self):
        """Test parsing hmmsearch 3.1b2 (text_31b2_hmmsearch_002)."""
        txt_file = get_file("text_31b2_hmmsearch_002.out")
        qresults = list(parse(txt_file, FMT))

        # Assert we deal with 2 queries
        self.assertEqual(len(qresults), 2)

        # Test whether proper query id is read
        self.assertEqual([x.id for x in qresults], ["infile_sto", "infile_sto2"])

        # Test if proper number of hits is read
        self.assertEqual([len(x.hits) for x in qresults], [10, 10])

    def test_31b2_hmmscan_001(self):
        """Test parsing hmmscan 3.1b2 (text_31b2_hmmscan_001)."""
        txt_file = get_file("text_31b2_hmmscan_001.out")
        qresults = list(parse(txt_file, FMT))

        # expects only 1 query result
        self.assertEqual(1, len(qresults))

        # in that query result, there are 108 hits
        qresult = qresults[0]
        self.assertEqual(108, len(qresult))

        # making sure all query IDs are equal
        self.assertEqual(
            "Protein-arginine kinase activator protein "
            "OS=Bacillus subtilis (strain 168) GN=mcsA PE=1 SV=1",
            qresult.description,
        )

        # and all hits have no descriptions
        for hit in qresult:
            self.assertEqual("", hit.description)

    def test_31b1_hmmsearch_001(self):
        """Test parsing hmmsearch 3.1b1 (text_31b1_hmmsearch_001)."""
        txt_file = get_file("text_31b1_hmmsearch_001.out")
        qresults = list(parse(txt_file, FMT))

        self.assertEqual(2, len(qresults))

        # first qresult is empty
        qresult = qresults[0]
        self.assertEqual("Globins", qresult.id)
        self.assertEqual(149, qresult.seq_len)
        self.assertEqual(0, len(qresult))

        # second qresult
        qresult = qresults[1]
        self.assertEqual("Pkinase", qresult.id)
        self.assertEqual(260, qresult.seq_len)
        self.assertEqual("PF00069.17", qresult.accession)
        self.assertEqual("Protein kinase domain", qresult.description)
        self.assertEqual(4, len(qresult))

        # first hit, first hsp
        hit = qresult[0]
        self.assertEqual("sp|Q9WUT3|KS6A2_MOUSE", hit.id)
        self.assertEqual(
            "Ribosomal protein S6 kinase alpha-2 OS=Mus musculus GN=Rps6ka2 PE=1 SV=1",
            hit.description,
        )
        self.assertTrue(hit.is_included)
        self.assertEqual(8.5e-147, hit.evalue)
        self.assertEqual(492.3, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual(2.1, hit.domain_exp_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))

        hsp = hit.hsps[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_included)
        self.assertEqual(241.2, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(2.6e-75, hsp.evalue_cond)
        self.assertEqual(3.6e-70, hsp.evalue)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(260, hsp.query_end)
        self.assertEqual("[]", hsp.query_endtype)
        self.assertEqual(58, hsp.hit_start)
        self.assertEqual(318, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)
        self.assertEqual(58, hsp.env_start)
        self.assertEqual(318, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.95, hsp.acc_avg)
        self.assertEqual(
            "EEEEEEEEEETTEEEEEEEE...TTTTEEEEEEEEEHHHCCCCCCHHHHHHHHHHHHHSSSSB--EEEEEEETTEEEEEEE--TS-BHHHHHHHHHST-HHHHHHHHHHHHHHHHHHHHTTEE-S--SGGGEEEETTTEEEE--GTT.E..EECSS-C-S--S-GGGS-HHHHCCS-CTHHHHHHHHHHHHHHHHHHSS-TTSSSHHCCTHHHHSSHHH......TTS.....HHHHHHHHHHT-SSGGGSTT.....HHHHHTSGGG",
            hsp.aln_annotation["CS"],
        )
        self.assertEqual(
            "yelleklGsGsfGkVykakk...kktgkkvAvKilkkeeekskkektavrElkilkklsHpnivkllevfetkdelylvleyveggdlfdllkkegklseeeikkialqilegleylHsngiiHrDLKpeNiLldkkgevkiaDFGlakkleksseklttlvgtreYmAPEvllkakeytkkvDvWslGvilyelltgklpfsgeseedqleliekilkkkleedepkssskseelkdlikkllekdpakRlt.....aeeilkhpwl",
            hsp.query.seq,
        )
        self.assertEqual(
            "FELLKVLGQGSYGKVFLVRKvtgSDAGQLYAMKVLKKATLKVRDRVRSKMERDILAEVNHPFIVKLHYAFQTEGKLYLILDFLRGGDLFTRLSKEVMFTEEDVKFYLAELALALDHLHGLGIIYRDLKPENILLDEEGHIKITDFGLSKEATDHDKRAYSFCGTIEYMAPEVVN-RRGHTQSADWWSFGVLMFEMLTGSLPFQGK---DRKETMALILKAKLGMPQFLS----AEAQSLLRALFKRNPCNRLGagvdgVEEIKRHPFF",
            hsp.hit.seq,
        )
        self.assertEqual(
            "67899**********7666611155667*****************99999****************************************************************************************************************************.******************************...999999999999999998866....99******************9999999*****997",
            hsp.aln_annotation["PP"],
        )

    def test_30_hmmsearch_001(self):
        """Parsing hmmersearch 3.0 (text_30_hmmsearch_001)."""
        txt_file = get_file("text_30_hmmsearch_001.out")
        qresults = parse(txt_file, FMT)
        counter = 0

        qresult = next(qresults)
        counter += 1

        self.assertEqual("hmmsearch", qresult.program)
        self.assertEqual("/home/bow/db/hmmer/uniprot_sprot.fasta", qresult.target)
        self.assertEqual("3.0", qresult.version)
        self.assertEqual("globins4", qresult.id)
        self.assertEqual(149, qresult.seq_len)
        self.assertEqual(0, len(qresult))

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(1, counter)

    def test_30_hmmsearch_002(self):
        """Parsing hmmersearch 3.0 (text_30_hmmsearch_002)."""
        txt_file = get_file("text_30_hmmsearch_002.out")
        qresults = parse(txt_file, FMT)
        counter = 0

        qresult = next(qresults)
        counter += 1

        self.assertEqual("hmmsearch", qresult.program)
        self.assertEqual("/home/bow/db/hmmer/uniprot_sprot.fasta", qresult.target)
        self.assertEqual("3.0", qresult.version)
        self.assertEqual("Pkinase", qresult.id)
        self.assertEqual("PF00069.17", qresult.accession)
        self.assertEqual("Protein kinase domain", qresult.description)
        self.assertEqual(260, qresult.seq_len)
        self.assertEqual(7, len(qresult))

        hit = qresult[0]
        self.assertEqual("sp|Q9WUT3|KS6A2_MOUSE", hit.id)
        self.assertEqual(
            "Ribosomal protein S6 kinase alpha-2 OS=Mus musculus GN=Rps6ka2 PE=2 SV=1",
            hit.description,
        )
        self.assertTrue(hit.is_included)
        self.assertEqual(8.4e-147, hit.evalue)
        self.assertEqual(492.3, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual(2.1, hit.domain_exp_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))

        hsp = hit.hsps[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_included)
        self.assertEqual(241.2, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(4.6e-75, hsp.evalue_cond)
        self.assertEqual(3.5e-70, hsp.evalue)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(260, hsp.query_end)
        self.assertEqual("[]", hsp.query_endtype)
        self.assertEqual(58, hsp.hit_start)
        self.assertEqual(318, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)
        self.assertEqual(58, hsp.env_start)
        self.assertEqual(318, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.95, hsp.acc_avg)
        self.assertEqual(
            "EEEEEEEEEETTEEEEEEEE...TTTTEEEEEEEEEHHHCCCCCCHHHHHHHHHHHHHSSSSB--EEEEEEETTEEEEEEE--TS-BHHHHHHHHHST-HHHHHHHHHHHHHHHHHHHHTTEE-S--SGGGEEEETTTEEEE--GTT.E..EECSS-C-S--S-GGGS-HHHHCCS-CTHHHHHHHHHHHHHHHHHHSS-TTSSSHHCCTHHHHSSHHH......TTS.....HHHHHHHHHHT-SSGGGSTT.....HHHHHTSGGG",
            hsp.aln_annotation["CS"],
        )
        self.assertEqual(
            "yelleklGsGsfGkVykakk...kktgkkvAvKilkkeeekskkektavrElkilkklsHpnivkllevfetkdelylvleyveggdlfdllkkegklseeeikkialqilegleylHsngiiHrDLKpeNiLldkkgevkiaDFGlakkleksseklttlvgtreYmAPEvllkakeytkkvDvWslGvilyelltgklpfsgeseedqleliekilkkkleedepkssskseelkdlikkllekdpakRlt.....aeeilkhpwl",
            hsp.query.seq,
        )
        self.assertEqual(
            "FELLKVLGQGSYGKVFLVRKvtgSDAGQLYAMKVLKKATLKVRDRVRSKMERDILAEVNHPFIVKLHYAFQTEGKLYLILDFLRGGDLFTRLSKEVMFTEEDVKFYLAELALALDHLHGLGIIYRDLKPENILLDEEGHIKITDFGLSKEATDHDKRAYSFCGTIEYMAPEVVN-RRGHTQSADWWSFGVLMFEMLTGSLPFQGK---DRKETMALILKAKLGMPQFLS----AEAQSLLRALFKRNPCNRLGagvdgVEEIKRHPFF",
            hsp.hit.seq,
        )
        self.assertEqual(
            "67899**********7666611155667*****************99999****************************************************************************************************************************.******************************...999999999999999998866....99******************9999999*****997",
            hsp.aln_annotation["PP"],
        )

        hsp = hit.hsps[1]
        self.assertEqual(2, hsp.domain_index)
        self.assertTrue(hsp.is_included)
        self.assertEqual(249.3, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(1.5e-77, hsp.evalue_cond)
        self.assertEqual(1.1e-72, hsp.evalue)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(260, hsp.query_end)
        self.assertEqual("[]", hsp.query_endtype)
        self.assertEqual(414, hsp.hit_start)
        self.assertEqual(672, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)
        self.assertEqual(414, hsp.env_start)
        self.assertEqual(672, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.97, hsp.acc_avg)
        self.assertEqual(
            "EEEEEEEEEETTEEEEEEEETTTTEEEEEEEEEHHHCCCCCCHHHHHHHHHHHHH.SSSSB--EEEEEEETTEEEEEEE--TS-BHHHHHHHHHST-HHHHHHHHHHHHHHHHHHHHTTEE-S--SGGGEEEETTTEE....EE--GTT.E..EECSS-C-S--S-GGGS-HHHHCCS-CTHHHHHHHHHHHHHHHHHHSS-TTSSSHHCCTHHHHSSHHH......TTS.....HHHHHHHHHHT-SSGGGSTTHHHHHTSGGG",
            hsp.aln_annotation["CS"],
        )
        self.assertEqual(
            "yelleklGsGsfGkVykakkkktgkkvAvKilkkeeekskkektavrElkilkkl.sHpnivkllevfetkdelylvleyveggdlfdllkkegklseeeikkialqilegleylHsngiiHrDLKpeNiLldkkgev....kiaDFGlakkleksseklttlvgtreYmAPEvllkakeytkkvDvWslGvilyelltgklpfsgeseedqleliekilkkkleedepkssskseelkdlikkllekdpakRltaeeilkhpwl",
            hsp.query.seq,
        )
        self.assertEqual(
            "YEIKEDIGVGSYSVCKRCVHKATDAEYAVKIIDKSKRDPSE------EIEILLRYgQHPNIITLKDVYDDGKYVYLVMELMRGGELLDRILRQRCFSEREASDVLYTIARTMDYLHSQGVVHRDLKPSNILYMDESGNpesiRICDFGFAKQLRAENGLLMTPCYTANFVAPEVLK-RQGYDAACDVWSLGILLYTMLAGFTPFANGPDDTPEEILARIGSGKYALSGGNWDSISDAAKDVVSKMLHVDPQQRLTAVQVLKHPWI",
            hsp.hit.seq,
        )
        self.assertEqual(
            "7899***********************************98......9*******99**************************************************************************98544444888**********************************.***************************************************************************************7",
            hsp.aln_annotation["PP"],
        )

        hit = qresult[-1]
        self.assertEqual("sp|P18654|KS6A3_MOUSE", hit.id)
        self.assertEqual(
            "Ribosomal protein S6 kinase alpha-3 OS=Mus musculus GN=Rps6ka3 PE=1 SV=2",
            hit.description,
        )
        self.assertTrue(hit.is_included)
        self.assertEqual(5e-144, hit.evalue)
        self.assertEqual(483.2, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual(2.2, hit.domain_exp_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))

        hsp = hit.hsps[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_included)
        self.assertEqual(240.3, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(8.7e-75, hsp.evalue_cond)
        self.assertEqual(6.6e-70, hsp.evalue)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(260, hsp.query_end)
        self.assertEqual("[]", hsp.query_endtype)
        self.assertEqual(67, hsp.hit_start)
        self.assertEqual(327, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)
        self.assertEqual(67, hsp.env_start)
        self.assertEqual(327, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.95, hsp.acc_avg)
        self.assertEqual(
            "EEEEEEEEEETTEEEEEEEETTTTE...EEEEEEEEHHHCCCCCCHHHHHHHHHHHHHSSSSB--EEEEEEETTEEEEEEE--TS-BHHHHHHHHHST-HHHHHHHHHHHHHHHHHHHHTTEE-S--SGGGEEEETTTEEEE--GTT.E..EECSS-C-S--S-GGGS-HHHHCCS-CTHHHHHHHHHHHHHHHHHHSS-TTSSSHHCCTHHHHSSHHH......TTS.....HHHHHHHHHHT-SSGGGSTT.....HHHHHTSGGG",
            hsp.aln_annotation["CS"],
        )
        self.assertEqual(
            "yelleklGsGsfGkVykakkkktgk...kvAvKilkkeeekskkektavrElkilkklsHpnivkllevfetkdelylvleyveggdlfdllkkegklseeeikkialqilegleylHsngiiHrDLKpeNiLldkkgevkiaDFGlakkleksseklttlvgtreYmAPEvllkakeytkkvDvWslGvilyelltgklpfsgeseedqleliekilkkkleedepkssskseelkdlikkllekdpakRlt.....aeeilkhpwl",
            hsp.query.seq,
        )
        self.assertEqual(
            "FELLKVLGQGSFGKVFLVKKISGSDarqLYAMKVLKKATLKVRDRVRTKMERDILVEVNHPFIVKLHYAFQTEGKLYLILDFLRGGDLFTRLSKEVMFTEEDVKFYLAELALALDHLHSLGIIYRDLKPENILLDEEGHIKLTDFGLSKESIDHEKKAYSFCGTVEYMAPEVVN-RRGHTQSADWWSFGVLMFEMLTGTLPFQGK---DRKETMTMILKAKLGMPQFLS----PEAQSLLRMLFKRNPANRLGagpdgVEEIKRHSFF",
            hsp.hit.seq,
        )
        self.assertEqual(
            "67899**********8888876655455****************999999****************************************************************************************************9***********************.******************************...999999999999988888866....9******************9888888999999886",
            hsp.aln_annotation["PP"],
        )

        hsp = hit.hsps[1]
        self.assertEqual(2, hsp.domain_index)
        self.assertTrue(hsp.is_included)
        self.assertEqual(241.0, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(5.1e-75, hsp.evalue_cond)
        self.assertEqual(3.9e-70, hsp.evalue)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(260, hsp.query_end)
        self.assertEqual("[]", hsp.query_endtype)
        self.assertEqual(421, hsp.hit_start)
        self.assertEqual(679, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)
        self.assertEqual(421, hsp.env_start)
        self.assertEqual(679, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.97, hsp.acc_avg)
        self.assertEqual(
            "EEEEEEEEEETTEEEEEEEETTTTEEEEEEEEEHHHCCCCCCHHHHHHHHHHHHH.SSSSB--EEEEEEETTEEEEEEE--TS-BHHHHHHHHHST-HHHHHHHHHHHHHHHHHHHHTTEE-S--SGGGEEEETTTEE....EE--GTT.E..EECSS-C-S--S-GGGS-HHHHCCS-CTHHHHHHHHHHHHHHHHHHSS-TTSSSHHCCTHHHHSSHHH......TTS.....HHHHHHHHHHT-SSGGGSTTHHHHHTSGGG",
            hsp.aln_annotation["CS"],
        )
        self.assertEqual(
            "yelleklGsGsfGkVykakkkktgkkvAvKilkkeeekskkektavrElkilkkl.sHpnivkllevfetkdelylvleyveggdlfdllkkegklseeeikkialqilegleylHsngiiHrDLKpeNiLldkkgev....kiaDFGlakkleksseklttlvgtreYmAPEvllkakeytkkvDvWslGvilyelltgklpfsgeseedqleliekilkkkleedepkssskseelkdlikkllekdpakRltaeeilkhpwl",
            hsp.query.seq,
        )
        self.assertEqual(
            "YEVKEDIGVGSYSVCKRCIHKATNMEFAVKIIDKSKRDPTE------EIEILLRYgQHPNIITLKDVYDDGKYVYVVTELMKGGELLDKILRQKFFSEREASAVLFTITKTVEYLHAQGVVHRDLKPSNILYVDESGNpesiRICDFGFAKQLRAENGLLMTPCYTANFVAPEVLK-RQGYDAACDIWSLGVLLYTMLTGYTPFANGPDDTPEEILARIGSGKFSLSGGYWNSVSDTAKDLVSKMLHVDPHQRLTAALVLRHPWI",
            hsp.hit.seq,
        )
        self.assertEqual(
            "7899***********************************88......9*******99**********************************9***************************************98554444888**********************************.***************************************************************************************7",
            hsp.aln_annotation["PP"],
        )

        # test if we've properly finished iteration
        # self.assertRaises(StopIteration, next, qresults)
        # self.assertEqual(1, counter)

    def test_30_hmmsearch_003(self):
        """Parsing hmmersearch 3.0 (text_30_hmmsearch_003)."""
        txt_file = get_file("text_30_hmmsearch_003.out")
        qresults = parse(txt_file, FMT)
        counter = 0

        qresult = next(qresults)
        counter += 1

        self.assertEqual("hmmsearch", qresult.program)
        self.assertEqual("/home/bow/db/hmmer/uniprot_sprot.fasta", qresult.target)
        self.assertEqual("3.0", qresult.version)
        self.assertEqual("Pkinase", qresult.id)
        self.assertEqual("PF00069.17", qresult.accession)
        self.assertEqual("Protein kinase domain", qresult.description)
        self.assertEqual(260, qresult.seq_len)
        self.assertEqual(7, len(qresult))

        hit = qresult[0]
        self.assertEqual("sp|Q9WUT3|KS6A2_MOUSE", hit.id)
        self.assertEqual(
            "Ribosomal protein S6 kinase alpha-2 OS=Mus musculus GN=Rps6ka2 PE=2 SV=1",
            hit.description,
        )
        self.assertTrue(hit.is_included)
        self.assertEqual(8.4e-147, hit.evalue)
        self.assertEqual(492.3, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual(2.1, hit.domain_exp_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))

        hsp = hit.hsps[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_included)
        self.assertEqual(241.2, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(4.6e-75, hsp.evalue_cond)
        self.assertEqual(3.5e-70, hsp.evalue)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(260, hsp.query_end)
        self.assertEqual("[]", hsp.query_endtype)
        self.assertEqual(58, hsp.hit_start)
        self.assertEqual(318, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)
        self.assertEqual(58, hsp.env_start)
        self.assertEqual(318, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.95, hsp.acc_avg)

        hsp = hit.hsps[1]
        self.assertEqual(2, hsp.domain_index)
        self.assertTrue(hsp.is_included)
        self.assertEqual(249.3, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(1.5e-77, hsp.evalue_cond)
        self.assertEqual(1.1e-72, hsp.evalue)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(260, hsp.query_end)
        self.assertEqual("[]", hsp.query_endtype)
        self.assertEqual(414, hsp.hit_start)
        self.assertEqual(672, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)
        self.assertEqual(414, hsp.env_start)
        self.assertEqual(672, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.97, hsp.acc_avg)

        hit = qresult[-1]
        self.assertEqual("sp|P18654|KS6A3_MOUSE", hit.id)
        self.assertEqual(
            "Ribosomal protein S6 kinase alpha-3 OS=Mus musculus GN=Rps6ka3 PE=1 SV=2",
            hit.description,
        )
        self.assertTrue(hit.is_included)
        self.assertEqual(5e-144, hit.evalue)
        self.assertEqual(483.2, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual(2.2, hit.domain_exp_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))

        hsp = hit.hsps[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_included)
        self.assertEqual(240.3, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(8.7e-75, hsp.evalue_cond)
        self.assertEqual(6.6e-70, hsp.evalue)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(260, hsp.query_end)
        self.assertEqual("[]", hsp.query_endtype)
        self.assertEqual(67, hsp.hit_start)
        self.assertEqual(327, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)
        self.assertEqual(67, hsp.env_start)
        self.assertEqual(327, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.95, hsp.acc_avg)

        hsp = hit.hsps[1]
        self.assertEqual(2, hsp.domain_index)
        self.assertTrue(hsp.is_included)
        self.assertEqual(241.0, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(5.1e-75, hsp.evalue_cond)
        self.assertEqual(3.9e-70, hsp.evalue)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(260, hsp.query_end)
        self.assertEqual("[]", hsp.query_endtype)
        self.assertEqual(421, hsp.hit_start)
        self.assertEqual(679, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)
        self.assertEqual(421, hsp.env_start)
        self.assertEqual(679, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.97, hsp.acc_avg)

        # test if we've properly finished iteration
        # self.assertRaises(StopIteration, next, qresults)
        # self.assertEqual(1, counter)

    def test_30_hmmsearch_004(self):
        """Parsing hmmersearch 3.0 (text_30_hmmsearch_004)."""
        txt_file = get_file("text_30_hmmsearch_004.out")
        qresults = parse(txt_file, FMT)
        counter = 0

        qresult = next(qresults)
        counter += 1

        self.assertEqual("hmmsearch", qresult.program)
        self.assertEqual("/home/bow/db/hmmer/uniprot_sprot.fasta", qresult.target)
        self.assertEqual("3.0", qresult.version)
        self.assertEqual("Pkinase", qresult.id)
        self.assertEqual("PF00069.17", qresult.accession)
        self.assertEqual("Protein kinase domain", qresult.description)
        self.assertEqual(260, qresult.seq_len)
        self.assertEqual(7, len(qresult))

        hit = qresult[0]
        self.assertEqual("sp|Q9WUT3|KS6A2_MOUSE", hit.id)
        self.assertEqual(
            "Ribosomal protein S6 kinase alpha-2 OS=Mus musculus GN=Rps6ka2 PE=2 SV=1",
            hit.description,
        )
        self.assertTrue(hit.is_included)
        self.assertEqual(8.4e-147, hit.evalue)
        self.assertEqual(492.3, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual(2.1, hit.domain_exp_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))

        hsp = hit.hsps[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_included)
        self.assertEqual(241.2, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(4.6e-75, hsp.evalue_cond)
        self.assertEqual(3.5e-70, hsp.evalue)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(260, hsp.query_end)
        self.assertEqual("[]", hsp.query_endtype)
        self.assertEqual(58, hsp.hit_start)
        self.assertEqual(318, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)
        self.assertEqual(58, hsp.env_start)
        self.assertEqual(318, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.95, hsp.acc_avg)
        self.assertEqual(
            "EEEEEEEEEETTEEEEEEEE...TTTTEEEEEEEEEHHHCCCCCCHHHHHHHHHHHHHSSSSB--EEEEEEETTEEEEEEE--TS-BHHHHHHHHHST-HHHHHHHHHHHHHHHHHHHHTTEE-S--SGGGEEEETTTEEEE--GTT.E..EECSS-C-S--S-GGGS-HHHHCCS-CTHHHHHHHHHHHHHHHHHHSS-TTSSSHHCCTHHHHSSHHH......TTS.....HHHHHHHHHHT-SSGGGSTT.....HHHHHTSGGG",
            hsp.aln_annotation["CS"],
        )
        self.assertEqual(
            "yelleklGsGsfGkVykakk...kktgkkvAvKilkkeeekskkektavrElkilkklsHpnivkllevfetkdelylvleyveggdlfdllkkegklseeeikkialqilegleylHsngiiHrDLKpeNiLldkkgevkiaDFGlakkleksseklttlvgtreYmAPEvllkakeytkkvDvWslGvilyelltgklpfsgeseedqleliekilkkkleedepkssskseelkdlikkllekdpakRlt.....aeeilkhpwl",
            hsp.query.seq,
        )
        self.assertEqual(
            "FELLKVLGQGSYGKVFLVRKvtgSDAGQLYAMKVLKKATLKVRDRVRSKMERDILAEVNHPFIVKLHYAFQTEGKLYLILDFLRGGDLFTRLSKEVMFTEEDVKFYLAELALALDHLHGLGIIYRDLKPENILLDEEGHIKITDFGLSKEATDHDKRAYSFCGTIEYMAPEVVN-RRGHTQSADWWSFGVLMFEMLTGSLPFQGK---DRKETMALILKAKLGMPQFLS----AEAQSLLRALFKRNPCNRLGagvdgVEEIKRHPFF",
            hsp.hit.seq,
        )
        self.assertEqual(
            "67899**********7666611155667*****************99999****************************************************************************************************************************.******************************...999999999999999998866....99******************9999999*****997",
            hsp.aln_annotation["PP"],
        )

        hsp = hit.hsps[1]
        self.assertEqual(2, hsp.domain_index)
        self.assertTrue(hsp.is_included)
        self.assertEqual(249.3, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(1.5e-77, hsp.evalue_cond)
        self.assertEqual(1.1e-72, hsp.evalue)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(260, hsp.query_end)
        self.assertEqual("[]", hsp.query_endtype)
        self.assertEqual(414, hsp.hit_start)
        self.assertEqual(672, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)
        self.assertEqual(414, hsp.env_start)
        self.assertEqual(672, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.97, hsp.acc_avg)
        self.assertEqual(
            "EEEEEEEEEETTEEEEEEEETTTTEEEEEEEEEHHHCCCCCCHHHHHHHHHHHHH.SSSSB--EEEEEEETTEEEEEEE--TS-BHHHHHHHHHST-HHHHHHHHHHHHHHHHHHHHTTEE-S--SGGGEEEETTTEE....EE--GTT.E..EECSS-C-S--S-GGGS-HHHHCCS-CTHHHHHHHHHHHHHHHHHHSS-TTSSSHHCCTHHHHSSHHH......TTS.....HHHHHHHHHHT-SSGGGSTTHHHHHTSGGG",
            hsp.aln_annotation["CS"],
        )
        self.assertEqual(
            "yelleklGsGsfGkVykakkkktgkkvAvKilkkeeekskkektavrElkilkkl.sHpnivkllevfetkdelylvleyveggdlfdllkkegklseeeikkialqilegleylHsngiiHrDLKpeNiLldkkgev....kiaDFGlakkleksseklttlvgtreYmAPEvllkakeytkkvDvWslGvilyelltgklpfsgeseedqleliekilkkkleedepkssskseelkdlikkllekdpakRltaeeilkhpwl",
            hsp.query.seq,
        )
        self.assertEqual(
            "YEIKEDIGVGSYSVCKRCVHKATDAEYAVKIIDKSKRDPSE------EIEILLRYgQHPNIITLKDVYDDGKYVYLVMELMRGGELLDRILRQRCFSEREASDVLYTIARTMDYLHSQGVVHRDLKPSNILYMDESGNpesiRICDFGFAKQLRAENGLLMTPCYTANFVAPEVLK-RQGYDAACDVWSLGILLYTMLAGFTPFANGPDDTPEEILARIGSGKYALSGGNWDSISDAAKDVVSKMLHVDPQQRLTAVQVLKHPWI",
            hsp.hit.seq,
        )
        self.assertEqual(
            "7899***********************************98......9*******99**************************************************************************98544444888**********************************.***************************************************************************************7",
            hsp.aln_annotation["PP"],
        )

        hit = qresult[-1]
        self.assertEqual("sp|P18654|KS6A3_MOUSE", hit.id)
        self.assertEqual(
            "Ribosomal protein S6 kinase alpha-3 OS=Mus musculus GN=Rps6ka3 PE=1 SV=2",
            hit.description,
        )
        self.assertTrue(hit.is_included)
        self.assertEqual(5e-144, hit.evalue)
        self.assertEqual(483.2, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual(2.2, hit.domain_exp_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))

        hsp = hit.hsps[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_included)
        self.assertEqual(240.3, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(8.7e-75, hsp.evalue_cond)
        self.assertEqual(6.6e-70, hsp.evalue)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(260, hsp.query_end)
        self.assertEqual("[]", hsp.query_endtype)
        self.assertEqual(67, hsp.hit_start)
        self.assertEqual(327, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)
        self.assertEqual(67, hsp.env_start)
        self.assertEqual(327, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.95, hsp.acc_avg)
        self.assertEqual(
            "EEEEEEEEEETTEEEEEEEETTTTE...EEEEEEEEHHHCCCCCCHHHHHHHHHHHHHSSSSB--EEEEEEETTEEEEEEE--TS-BHHHHHHHHHST-HHHHHHHHHHHHHHHHHHHHTTEE-S--SGGGEEEETTTEEEE--GTT.E..EECSS-C-S--S-GGGS-HHHHCCS-CTHHHHHHHHHHHHHHHHHHSS-TTSSSHHCCTHHHHSSHHH......TTS.....HHHHHHHHHHT-SSGGGSTT.....HHHHHTSGGG",
            hsp.aln_annotation["CS"],
        )
        self.assertEqual(
            "yelleklGsGsfGkVykakkkktgk...kvAvKilkkeeekskkektavrElkilkklsHpnivkllevfetkdelylvleyveggdlfdllkkegklseeeikkialqilegleylHsngiiHrDLKpeNiLldkkgevkiaDFGlakkleksseklttlvgtreYmAPEvllkakeytkkvDvWslGvilyelltgklpfsgeseedqleliekilkkkleedepkssskseelkdlikkllekdpakRlt.....aeeilkhpwl",
            hsp.query.seq,
        )
        self.assertEqual(
            "FELLKVLGQGSFGKVFLVKKISGSDarqLYAMKVLKKATLKVRDRVRTKMERDILVEVNHPFIVKLHYAFQTEGKLYLILDFLRGGDLFTRLSKEVMFTEEDVKFYLAELALALDHLHSLGIIYRDLKPENILLDEEGHIKLTDFGLSKESIDHEKKAYSFCGTVEYMAPEVVN-RRGHTQSADWWSFGVLMFEMLTGTLPFQGK---DRKETMTMILKAKLGMPQFLS----PEAQSLLRMLFKRNPANRLGagpdgVEEIKRHSFF",
            hsp.hit.seq,
        )
        self.assertEqual(
            "67899**********8888876655455****************999999****************************************************************************************************9***********************.******************************...999999999999988888866....9******************9888888999999886",
            hsp.aln_annotation["PP"],
        )

        hsp = hit.hsps[1]
        self.assertEqual(2, hsp.domain_index)
        self.assertTrue(hsp.is_included)
        self.assertEqual(241.0, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(5.1e-75, hsp.evalue_cond)
        self.assertEqual(3.9e-70, hsp.evalue)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(260, hsp.query_end)
        self.assertEqual("[]", hsp.query_endtype)
        self.assertEqual(421, hsp.hit_start)
        self.assertEqual(679, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)
        self.assertEqual(421, hsp.env_start)
        self.assertEqual(679, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.97, hsp.acc_avg)
        self.assertEqual(
            "EEEEEEEEEETTEEEEEEEETTTTEEEEEEEEEHHHCCCCCCHHHHHHHHHHHHH.SSSSB--EEEEEEETTEEEEEEE--TS-BHHHHHHHHHST-HHHHHHHHHHHHHHHHHHHHTTEE-S--SGGGEEEETTTEE....EE--GTT.E..EECSS-C-S--S-GGGS-HHHHCCS-CTHHHHHHHHHHHHHHHHHHSS-TTSSSHHCCTHHHHSSHHH......TTS.....HHHHHHHHHHT-SSGGGSTTHHHHHTSGGG",
            hsp.aln_annotation["CS"],
        )
        self.assertEqual(
            "yelleklGsGsfGkVykakkkktgkkvAvKilkkeeekskkektavrElkilkkl.sHpnivkllevfetkdelylvleyveggdlfdllkkegklseeeikkialqilegleylHsngiiHrDLKpeNiLldkkgev....kiaDFGlakkleksseklttlvgtreYmAPEvllkakeytkkvDvWslGvilyelltgklpfsgeseedqleliekilkkkleedepkssskseelkdlikkllekdpakRltaeeilkhpwl",
            hsp.query.seq,
        )
        self.assertEqual(
            "YEVKEDIGVGSYSVCKRCIHKATNMEFAVKIIDKSKRDPTE------EIEILLRYgQHPNIITLKDVYDDGKYVYVVTELMKGGELLDKILRQKFFSEREASAVLFTITKTVEYLHAQGVVHRDLKPSNILYVDESGNpesiRICDFGFAKQLRAENGLLMTPCYTANFVAPEVLK-RQGYDAACDIWSLGVLLYTMLTGYTPFANGPDDTPEEILARIGSGKFSLSGGYWNSVSDTAKDLVSKMLHVDPHQRLTAALVLRHPWI",
            hsp.hit.seq,
        )
        self.assertEqual(
            "7899***********************************88......9*******99**********************************9***************************************98554444888**********************************.***************************************************************************************7",
            hsp.aln_annotation["PP"],
        )

        # test if we've properly finished iteration
        # self.assertRaises(StopIteration, next, qresults)
        # self.assertEqual(1, counter)

    def test_30_hmmsearch_005(self):
        """Parsing hmmersearch 3.0 (text_30_hmmsearch_005)."""
        txt_file = get_file("text_30_hmmsearch_005.out")
        qresults = parse(txt_file, FMT)
        counter = 0

        # test first qresult
        qresult = next(qresults)
        counter += 1

        self.assertEqual("hmmsearch", qresult.program)
        self.assertEqual("/home/bow/db/hmmer/uniprot_sprot.fasta", qresult.target)
        self.assertEqual("3.0", qresult.version)
        self.assertEqual("globins4", qresult.id)
        self.assertEqual(149, qresult.seq_len)
        self.assertEqual(0, len(qresult))

        # test second qresult
        qresult = next(qresults)
        counter += 1

        self.assertEqual("hmmsearch", qresult.program)
        self.assertEqual("/home/bow/db/hmmer/uniprot_sprot.fasta", qresult.target)
        self.assertEqual("3.0", qresult.version)
        self.assertEqual("Pkinase", qresult.id)
        self.assertEqual("PF00069.17", qresult.accession)
        self.assertEqual("Protein kinase domain", qresult.description)
        self.assertEqual(260, qresult.seq_len)
        self.assertEqual(7, len(qresult))

        hit = qresult[0]
        self.assertEqual("sp|Q9WUT3|KS6A2_MOUSE", hit.id)
        self.assertEqual(
            "Ribosomal protein S6 kinase alpha-2 OS=Mus musculus GN=Rps6ka2 PE=2 SV=1",
            hit.description,
        )
        self.assertTrue(hit.is_included)
        self.assertEqual(8.4e-147, hit.evalue)
        self.assertEqual(492.3, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual(2.1, hit.domain_exp_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))

        hsp = hit.hsps[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_included)
        self.assertEqual(241.2, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(4.6e-75, hsp.evalue_cond)
        self.assertEqual(3.5e-70, hsp.evalue)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(260, hsp.query_end)
        self.assertEqual("[]", hsp.query_endtype)
        self.assertEqual(58, hsp.hit_start)
        self.assertEqual(318, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)
        self.assertEqual(58, hsp.env_start)
        self.assertEqual(318, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.95, hsp.acc_avg)
        self.assertEqual(
            "EEEEEEEEEETTEEEEEEEE...TTTTEEEEEEEEEHHHCCCCCCHHHHHHHHHHHHHSSSSB--EEEEEEETTEEEEEEE--TS-BHHHHHHHHHST-HHHHHHHHHHHHHHHHHHHHTTEE-S--SGGGEEEETTTEEEE--GTT.E..EECSS-C-S--S-GGGS-HHHHCCS-CTHHHHHHHHHHHHHHHHHHSS-TTSSSHHCCTHHHHSSHHH......TTS.....HHHHHHHHHHT-SSGGGSTT.....HHHHHTSGGG",
            hsp.aln_annotation["CS"],
        )
        self.assertEqual(
            "yelleklGsGsfGkVykakk...kktgkkvAvKilkkeeekskkektavrElkilkklsHpnivkllevfetkdelylvleyveggdlfdllkkegklseeeikkialqilegleylHsngiiHrDLKpeNiLldkkgevkiaDFGlakkleksseklttlvgtreYmAPEvllkakeytkkvDvWslGvilyelltgklpfsgeseedqleliekilkkkleedepkssskseelkdlikkllekdpakRlt.....aeeilkhpwl",
            hsp.query.seq,
        )
        self.assertEqual(
            "FELLKVLGQGSYGKVFLVRKvtgSDAGQLYAMKVLKKATLKVRDRVRSKMERDILAEVNHPFIVKLHYAFQTEGKLYLILDFLRGGDLFTRLSKEVMFTEEDVKFYLAELALALDHLHGLGIIYRDLKPENILLDEEGHIKITDFGLSKEATDHDKRAYSFCGTIEYMAPEVVN-RRGHTQSADWWSFGVLMFEMLTGSLPFQGK---DRKETMALILKAKLGMPQFLS----AEAQSLLRALFKRNPCNRLGagvdgVEEIKRHPFF",
            hsp.hit.seq,
        )
        self.assertEqual(
            "67899**********7666611155667*****************99999****************************************************************************************************************************.******************************...999999999999999998866....99******************9999999*****997",
            hsp.aln_annotation["PP"],
        )

        hsp = hit.hsps[1]
        self.assertEqual(2, hsp.domain_index)
        self.assertTrue(hsp.is_included)
        self.assertEqual(249.3, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(1.5e-77, hsp.evalue_cond)
        self.assertEqual(1.1e-72, hsp.evalue)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(260, hsp.query_end)
        self.assertEqual("[]", hsp.query_endtype)
        self.assertEqual(414, hsp.hit_start)
        self.assertEqual(672, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)
        self.assertEqual(414, hsp.env_start)
        self.assertEqual(672, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.97, hsp.acc_avg)
        self.assertEqual(
            "EEEEEEEEEETTEEEEEEEETTTTEEEEEEEEEHHHCCCCCCHHHHHHHHHHHHH.SSSSB--EEEEEEETTEEEEEEE--TS-BHHHHHHHHHST-HHHHHHHHHHHHHHHHHHHHTTEE-S--SGGGEEEETTTEE....EE--GTT.E..EECSS-C-S--S-GGGS-HHHHCCS-CTHHHHHHHHHHHHHHHHHHSS-TTSSSHHCCTHHHHSSHHH......TTS.....HHHHHHHHHHT-SSGGGSTTHHHHHTSGGG",
            hsp.aln_annotation["CS"],
        )
        self.assertEqual(
            "yelleklGsGsfGkVykakkkktgkkvAvKilkkeeekskkektavrElkilkkl.sHpnivkllevfetkdelylvleyveggdlfdllkkegklseeeikkialqilegleylHsngiiHrDLKpeNiLldkkgev....kiaDFGlakkleksseklttlvgtreYmAPEvllkakeytkkvDvWslGvilyelltgklpfsgeseedqleliekilkkkleedepkssskseelkdlikkllekdpakRltaeeilkhpwl",
            hsp.query.seq,
        )
        self.assertEqual(
            "YEIKEDIGVGSYSVCKRCVHKATDAEYAVKIIDKSKRDPSE------EIEILLRYgQHPNIITLKDVYDDGKYVYLVMELMRGGELLDRILRQRCFSEREASDVLYTIARTMDYLHSQGVVHRDLKPSNILYMDESGNpesiRICDFGFAKQLRAENGLLMTPCYTANFVAPEVLK-RQGYDAACDVWSLGILLYTMLAGFTPFANGPDDTPEEILARIGSGKYALSGGNWDSISDAAKDVVSKMLHVDPQQRLTAVQVLKHPWI",
            hsp.hit.seq,
        )
        self.assertEqual(
            "7899***********************************98......9*******99**************************************************************************98544444888**********************************.***************************************************************************************7",
            hsp.aln_annotation["PP"],
        )

        hit = qresult[-1]
        self.assertEqual("sp|P18654|KS6A3_MOUSE", hit.id)
        self.assertEqual(
            "Ribosomal protein S6 kinase alpha-3 OS=Mus musculus GN=Rps6ka3 PE=1 SV=2",
            hit.description,
        )
        self.assertTrue(hit.is_included)
        self.assertEqual(5e-144, hit.evalue)
        self.assertEqual(483.2, hit.bitscore)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual(2.2, hit.domain_exp_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))

        hsp = hit.hsps[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_included)
        self.assertEqual(240.3, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(8.7e-75, hsp.evalue_cond)
        self.assertEqual(6.6e-70, hsp.evalue)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(260, hsp.query_end)
        self.assertEqual("[]", hsp.query_endtype)
        self.assertEqual(67, hsp.hit_start)
        self.assertEqual(327, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)
        self.assertEqual(67, hsp.env_start)
        self.assertEqual(327, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.95, hsp.acc_avg)
        self.assertEqual(
            "EEEEEEEEEETTEEEEEEEETTTTE...EEEEEEEEHHHCCCCCCHHHHHHHHHHHHHSSSSB--EEEEEEETTEEEEEEE--TS-BHHHHHHHHHST-HHHHHHHHHHHHHHHHHHHHTTEE-S--SGGGEEEETTTEEEE--GTT.E..EECSS-C-S--S-GGGS-HHHHCCS-CTHHHHHHHHHHHHHHHHHHSS-TTSSSHHCCTHHHHSSHHH......TTS.....HHHHHHHHHHT-SSGGGSTT.....HHHHHTSGGG",
            hsp.aln_annotation["CS"],
        )
        self.assertEqual(
            "yelleklGsGsfGkVykakkkktgk...kvAvKilkkeeekskkektavrElkilkklsHpnivkllevfetkdelylvleyveggdlfdllkkegklseeeikkialqilegleylHsngiiHrDLKpeNiLldkkgevkiaDFGlakkleksseklttlvgtreYmAPEvllkakeytkkvDvWslGvilyelltgklpfsgeseedqleliekilkkkleedepkssskseelkdlikkllekdpakRlt.....aeeilkhpwl",
            hsp.query.seq,
        )
        self.assertEqual(
            "FELLKVLGQGSFGKVFLVKKISGSDarqLYAMKVLKKATLKVRDRVRTKMERDILVEVNHPFIVKLHYAFQTEGKLYLILDFLRGGDLFTRLSKEVMFTEEDVKFYLAELALALDHLHSLGIIYRDLKPENILLDEEGHIKLTDFGLSKESIDHEKKAYSFCGTVEYMAPEVVN-RRGHTQSADWWSFGVLMFEMLTGTLPFQGK---DRKETMTMILKAKLGMPQFLS----PEAQSLLRMLFKRNPANRLGagpdgVEEIKRHSFF",
            hsp.hit.seq,
        )
        self.assertEqual(
            "67899**********8888876655455****************999999****************************************************************************************************9***********************.******************************...999999999999988888866....9******************9888888999999886",
            hsp.aln_annotation["PP"],
        )

        hsp = hit.hsps[1]
        self.assertEqual(2, hsp.domain_index)
        self.assertTrue(hsp.is_included)
        self.assertEqual(241.0, hsp.bitscore)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(5.1e-75, hsp.evalue_cond)
        self.assertEqual(3.9e-70, hsp.evalue)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(260, hsp.query_end)
        self.assertEqual("[]", hsp.query_endtype)
        self.assertEqual(421, hsp.hit_start)
        self.assertEqual(679, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)
        self.assertEqual(421, hsp.env_start)
        self.assertEqual(679, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.97, hsp.acc_avg)
        self.assertEqual(
            "EEEEEEEEEETTEEEEEEEETTTTEEEEEEEEEHHHCCCCCCHHHHHHHHHHHHH.SSSSB--EEEEEEETTEEEEEEE--TS-BHHHHHHHHHST-HHHHHHHHHHHHHHHHHHHHTTEE-S--SGGGEEEETTTEE....EE--GTT.E..EECSS-C-S--S-GGGS-HHHHCCS-CTHHHHHHHHHHHHHHHHHHSS-TTSSSHHCCTHHHHSSHHH......TTS.....HHHHHHHHHHT-SSGGGSTTHHHHHTSGGG",
            hsp.aln_annotation["CS"],
        )
        self.assertEqual(
            "yelleklGsGsfGkVykakkkktgkkvAvKilkkeeekskkektavrElkilkkl.sHpnivkllevfetkdelylvleyveggdlfdllkkegklseeeikkialqilegleylHsngiiHrDLKpeNiLldkkgev....kiaDFGlakkleksseklttlvgtreYmAPEvllkakeytkkvDvWslGvilyelltgklpfsgeseedqleliekilkkkleedepkssskseelkdlikkllekdpakRltaeeilkhpwl",
            hsp.query.seq,
        )
        self.assertEqual(
            "YEVKEDIGVGSYSVCKRCIHKATNMEFAVKIIDKSKRDPTE------EIEILLRYgQHPNIITLKDVYDDGKYVYVVTELMKGGELLDKILRQKFFSEREASAVLFTITKTVEYLHAQGVVHRDLKPSNILYVDESGNpesiRICDFGFAKQLRAENGLLMTPCYTANFVAPEVLK-RQGYDAACDIWSLGVLLYTMLTGYTPFANGPDDTPEEILARIGSGKFSLSGGYWNSVSDTAKDLVSKMLHVDPHQRLTAALVLRHPWI",
            hsp.hit.seq,
        )
        self.assertEqual(
            "7899***********************************88......9*******99**********************************9***************************************98554444888**********************************.***************************************************************************************7",
            hsp.aln_annotation["PP"],
        )

        # test if we've properly finished iteration
        self.assertRaises(StopIteration, next, qresults)
        self.assertEqual(2, counter)


class PhmmerCases(unittest.TestCase):
    """Testing phmmer output."""

    def test_31b2_phmmer_001(self):
        """Parsing phmmer 3.1b2 (text_31b2_phmmer_001)."""
        txt_file = get_file("text_31b2_phmmer_001.out")
        qresults = list(parse(txt_file, FMT))

        # first qresult
        qresult = qresults[0]

        self.assertEqual("phmmer", qresult.program)
        self.assertEqual(
            "/home/bow/devel/sandbox/biopython-sandbox/db/hmmer/protdb/uniprot_sprot.fasta",
            qresult.target,
        )
        self.assertEqual("3.1b2", qresult.version)
        self.assertEqual("sp|Q6GZX4|001R_FRG3G", qresult.id)
        self.assertEqual(256, qresult.seq_len)
        self.assertEqual(13, len(qresult))

        hit = qresult[0]
        self.assertEqual("sp|Q6GZX4|001R_FRG3G", hit.id)
        self.assertEqual(
            "Putative transcription factor 001R OS=Frog virus 3 (isolate Goorha) GN=FV3-001R PE=4 SV=1",
            hit.description,
        )
        self.assertTrue(hit.is_included)
        self.assertEqual(1.1e-176, hit.evalue)
        self.assertEqual(590.0, hit.bitscore)
        self.assertEqual(1.4, hit.bias)
        self.assertEqual(1.0, hit.domain_exp_num)
        self.assertEqual(1, hit.domain_obs_num)
        self.assertEqual(1, len(hit))

        hsp = hit.hsps[0]
        self.assertEqual(1, hsp.domain_index)
        self.assertTrue(hsp.is_included)
        self.assertEqual(589.9, hsp.bitscore)
        self.assertEqual(1.4, hsp.bias)
        self.assertEqual(3e-181, hsp.evalue_cond)
        self.assertEqual(1.3e-176, hsp.evalue)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(256, hsp.hit_end)
        self.assertEqual("[]", hsp.hit_endtype)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(256, hsp.query_end)
        self.assertEqual("[]", hsp.query_endtype)
        self.assertEqual(0, hsp.env_start)
        self.assertEqual(256, hsp.env_end)
        self.assertEqual("[]", hsp.env_endtype)
        self.assertEqual(1.00, hsp.acc_avg)
        self.assertEqual(
            "mafsaedvlkeydrrrrmealllslyypndrklldykewspprvqvecpkapvewnnppsekglivghfsgikykgekaqasevdvnkm",
            hsp.query.seq[:89],
        )
        self.assertEqual(
            "MAFSAEDVLKEYDRRRRMEALLLSLYYPNDRKLLDYKEWSPPRVQVECPKAPVEWNNPPSEKGLIVGHFSGIKYKGEKAQASEVDVNKM",
            hsp.hit.seq[:89],
        )
        self.assertEqual(
            "89***************************************************************************************",
            hsp.aln_annotation["PP"][:89],
        )

        # last query, last hit
        qresult = qresults[-1]

        self.assertEqual("phmmer", qresult.program)
        self.assertEqual(
            "/home/bow/devel/sandbox/biopython-sandbox/db/hmmer/protdb/uniprot_sprot.fasta",
            qresult.target,
        )
        self.assertEqual("3.1b2", qresult.version)
        self.assertEqual("sp|Q197F7|003L_IIV3", qresult.id)
        self.assertEqual(156, qresult.seq_len)
        self.assertEqual(6, len(qresult))

        hit = qresult[-1]
        self.assertEqual("sp|P04060|RNAS1_HYSCR", hit.id)
        self.assertEqual(
            "Ribonuclease pancreatic OS=Hystrix cristata GN=RNASE1 PE=1 SV=1",
            hit.description,
        )
        self.assertFalse(hit.is_included)
        self.assertEqual(2.1, hit.evalue)
        self.assertEqual(13.1, hit.bitscore)
        self.assertEqual(0.7, hit.bias)
        self.assertEqual(1.8, hit.domain_exp_num)
        self.assertEqual(2, hit.domain_obs_num)
        self.assertEqual(2, len(hit))

        hsp = hit.hsps[-1]
        self.assertEqual(2, hsp.domain_index)
        self.assertFalse(hsp.is_included)
        self.assertEqual(12.1, hsp.bitscore)
        self.assertEqual(0.1, hsp.bias)
        self.assertEqual(4.7e-05, hsp.evalue_cond)
        self.assertEqual(4.2, hsp.evalue)
        self.assertEqual(91, hsp.hit_start)
        self.assertEqual(121, hsp.hit_end)
        self.assertEqual("..", hsp.hit_endtype)
        self.assertEqual(7, hsp.query_start)
        self.assertEqual(37, hsp.query_end)
        self.assertEqual("..", hsp.query_endtype)
        self.assertEqual(84, hsp.env_start)
        self.assertEqual(127, hsp.env_end)
        self.assertEqual("..", hsp.env_endtype)
        self.assertEqual(0.84, hsp.acc_avg)
        self.assertEqual("cpqswygspqlereivckmsgaphypnyyp", hsp.query.seq)
        self.assertEqual("YPDCSYGMSQLERSIVVACEGSPYVPVHFD", hsp.hit.seq)
        self.assertEqual("68889******************9887764", hsp.aln_annotation["PP"])


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
