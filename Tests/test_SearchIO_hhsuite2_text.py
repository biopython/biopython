# Copyright 2019 by Jens Thomas.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO HhsuiteIO parsers."""


import os
import unittest

from Bio.SearchIO import parse


# test case files are in the Blast directory
TEST_DIR = "HHsuite"
FMT = "hhsuite2-text"


def get_file(filename):
    """Return the path of a test file."""
    return os.path.join(TEST_DIR, filename)


class HhsuiteCases(unittest.TestCase):
    """Test hhsuite2 output."""

    def test_2uvo(self):
        """Parsing 2uvo."""
        txt_file = get_file("2uvo_hhblits.hhr")
        qresults = parse(txt_file, FMT)

        # test first and only qresult
        qresult = next(qresults)

        num_hits = 16
        self.assertEqual("HHSUITE", qresult.program)
        self.assertEqual("2UVO:A|PDBID|CHAIN|SEQUENCE", qresult.id)
        self.assertEqual(171, qresult.seq_len)
        self.assertEqual(num_hits, len(qresult))

        hit = qresult[0]
        self.assertEqual("2uvo_A", hit.id)
        self.assertEqual(
            "Agglutinin isolectin 1; carbohydrate-binding protein, hevein domain, chitin-binding,"
            " GERM agglutinin, chitin-binding protein; HET: NDG NAG GOL; 1.40A {Triticum aestivum}"
            " PDB: 1wgc_A* 2cwg_A* 2x3t_A* 4aml_A* 7wga_A 9wga_A 2wgc_A 1wgt_A 1k7t_A* 1k7v_A* 1k7u_A"
            " 2x52_A* 1t0w_A*",
            hit.description,
        )
        self.assertTrue(hit.is_included)
        self.assertEqual(3.7e-34, hit.evalue)
        self.assertEqual(210.31, hit.score)
        self.assertEqual(2, len(hit))

        hsp = hit.hsps[0]
        self.assertTrue(hsp.is_included)
        self.assertEqual(0, hsp.output_index)
        self.assertEqual(99.95, hsp.prob)
        self.assertEqual(210.31, hsp.score)
        self.assertEqual(3.7e-34, hsp.evalue)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(171, hsp.hit_end)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(171, hsp.query_end)
        self.assertEqual(
            "ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYC"
            "GAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCS"
            "KWGSCGIGPGYCGAGCQSGGCDG",
            hsp.hit.seq,
        )
        self.assertEqual(
            "ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYC"
            "GAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCS"
            "KWGSCGIGPGYCGAGCQSGGCDG",
            hsp.query.seq,
        )

        # Check last hit
        hit = qresult[num_hits - 1]
        self.assertEqual("4z8i_A", hit.id)
        self.assertEqual(
            "BBTPGRP3, peptidoglycan recognition protein 3; chitin-binding domain, "
            "AM hydrolase; 2.70A {Branchiostoma belcheri tsingtauense}",
            hit.description,
        )
        self.assertTrue(hit.is_included)
        self.assertEqual(0.11, hit.evalue)
        self.assertEqual(36.29, hit.score)
        self.assertEqual(2, len(hit))

        # Check we can get the original last HSP from the file.
        num_hsps = 32
        self.assertEqual(num_hsps, len(qresult.hsps))

        hsp = qresult.hsps[-1]
        self.assertTrue(hsp.is_included)
        self.assertEqual(num_hsps - 1, hsp.output_index)
        self.assertEqual(2.6, hsp.evalue)
        self.assertEqual(25.90, hsp.score)
        self.assertEqual(40.43, hsp.prob)
        self.assertEqual(10, hsp.hit_start)
        self.assertEqual(116, hsp.hit_end)
        self.assertEqual(53, hsp.query_start)
        self.assertEqual(163, hsp.query_end)
        self.assertEqual(
            "XCXXXXCCXXXXXCXXXXXXCXXXCXXXXCXXXXXCXXX--XXXCXXXXCCXXXXXCXXXXXXCXXXCXXXXCXXXXXCX"
            "XX--XXXCXXXXCCXXXXXCXXXXXXCXXX",
            hsp.hit.seq,
        )
        self.assertEqual(
            "TCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCG"
            "KDAGGRVCTNNYCCSKWGSCGIGPGYCGAG",
            hsp.query.seq,
        )

    def test_2uvo_onlyheader(self):
        """Parsing 4uvo with only header present."""
        txt_file = get_file("2uvo_hhblits_onlyheader.hhr")
        qresults = parse(txt_file, FMT)

        with self.assertRaises(RuntimeError):
            next(qresults)

    def test_2uvo_emptytable(self):
        """Parsing 4uvo with empty results table."""
        txt_file = get_file("2uvo_hhblits_emptytable.hhr")
        qresults = parse(txt_file, FMT)

        with self.assertRaises(RuntimeError):
            next(qresults)

    def test_allx(self):
        """Parsing allx.hhr file."""
        txt_file = get_file("allx.hhr")
        qresults = parse(txt_file, FMT)

        # test first and only qresult
        qresult = next(qresults)

        num_hits = 10
        self.assertEqual("HHSUITE", qresult.program)
        self.assertEqual("Only X amino acids", qresult.id)
        self.assertEqual(39, qresult.seq_len)
        self.assertEqual(num_hits, len(qresult))

        hit = qresult[0]
        self.assertEqual("1klr_A", hit.id)
        self.assertEqual(
            "Zinc finger Y-chromosomal protein; transcription; NMR {Synthetic} SCOP: g.37.1.1 PDB: "
            "5znf_A 1kls_A 1xrz_A* 7znf_A",
            hit.description,
        )
        self.assertTrue(hit.is_included)
        self.assertEqual(3.4e04, hit.evalue)
        self.assertEqual(-0.01, hit.score)
        self.assertEqual(1, len(hit))

        hsp = hit.hsps[0]
        self.assertTrue(hsp.is_included)
        self.assertEqual(0, hsp.output_index)
        self.assertEqual(3.4e04, hsp.evalue)
        self.assertEqual(-0.01, hsp.score)
        self.assertEqual(0.04, hsp.prob)
        self.assertEqual(23, hsp.hit_start)
        self.assertEqual(24, hsp.hit_end)
        self.assertEqual(38, hsp.query_start)
        self.assertEqual(39, hsp.query_end)
        self.assertEqual("T", hsp.hit.seq)
        self.assertEqual("X", hsp.query.seq)

        # Check last hit
        hit = qresult[num_hits - 1]
        self.assertEqual("1zfd_A", hit.id)
        self.assertEqual(
            "SWI5; DNA binding motif, zinc finger DNA binding domain; NMR {Saccharomyces cerevisiae}"
            " SCOP: g.37.1.1",
            hit.description,
        )
        self.assertTrue(hit.is_included)
        self.assertEqual(3.6e04, hit.evalue)
        self.assertEqual(0.03, hit.score)
        self.assertEqual(1, len(hit))

        # Check we can get the original last HSP from the file.
        num_hsps = num_hits
        self.assertEqual(num_hsps, len(qresult.hsps))
        hsp = qresult.hsps[-1]

        self.assertTrue(hsp.is_included)
        self.assertEqual(num_hsps - 1, hsp.output_index)
        self.assertEqual(3.6e04, hsp.evalue)
        self.assertEqual(0.03, hsp.score)
        self.assertEqual(0.03, hsp.prob)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(1, hsp.hit_end)
        self.assertEqual(3, hsp.query_start)
        self.assertEqual(4, hsp.query_end)
        self.assertEqual("D", hsp.hit.seq)
        self.assertEqual("X", hsp.query.seq)

    def test_4y9h_nossm(self):
        """Parsing 4y9h_hhsearch_server_NOssm.hhr file."""
        txt_file = get_file("4y9h_hhsearch_server_NOssm.hhr")
        qresults = parse(txt_file, FMT)

        # test first and only qresult
        qresult = next(qresults)

        num_hits = 29
        self.assertEqual("HHSUITE", qresult.program)
        self.assertEqual("4Y9H:A|PDBID|CHAIN|SEQUENCE", qresult.id)
        self.assertEqual(226, qresult.seq_len)
        self.assertEqual(num_hits, len(qresult))

        hit = qresult[0]
        self.assertEqual("5ZIM_A", hit.id)
        self.assertEqual(
            "Bacteriorhodopsin; proton pump, membrane protein, PROTON; HET: L2P, RET; 1.25A {Halobacterium"
            " salinarum}; Related PDB entries: 1R84_A 1KG8_A 1KME_B 1KGB_A 1KG9_A 1KME_A 4X31_A 5ZIL_A 1E0P_A "
            "4X32_A 5ZIN_A 1S53_B 1S51_B 1S53_A 1S54_A 1F50_A 1S54_B 1S51_A 1F4Z_A 5J7A_A 1S52_B 1S52_A 4Y9H_A "
            "3T45_A 3T45_C 3T45_B 1C3W_A 1L0M_A",
            hit.description,
        )
        self.assertTrue(hit.is_included)
        self.assertEqual(2.1e-48, hit.evalue)
        self.assertEqual(320.44, hit.score)
        self.assertEqual(1, len(hit))

        hsp = hit.hsps[0]
        self.assertTrue(hsp.is_included)
        self.assertEqual(0, hsp.output_index)
        self.assertEqual(2.1e-48, hsp.evalue)
        self.assertEqual(320.44, hsp.score)
        self.assertEqual(100.00, hsp.prob)
        self.assertEqual(1, hsp.hit_start)
        self.assertEqual(227, hsp.hit_end)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(226, hsp.query_end)
        self.assertEqual(
            "GRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYAD"
            "WLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRP"
            "EVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFG",
            hsp.hit.seq,
        )
        self.assertEqual(
            "GRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYAD"
            "WLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRP"
            "EVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFG",
            hsp.query.seq,
        )

        # Check last hit
        hit = qresult[num_hits - 1]
        self.assertEqual("5ABB_Z", hit.id)
        self.assertEqual(
            "PROTEIN TRANSLOCASE SUBUNIT SECY, PROTEIN; TRANSLATION, RIBOSOME, MEMBRANE PROTEIN, "
            "TRANSLOCON; 8.0A {ESCHERICHIA COLI}",
            hit.description,
        )
        self.assertTrue(hit.is_included)
        self.assertEqual(3.3e-05, hit.evalue)
        self.assertEqual(51.24, hit.score)
        self.assertEqual(1, len(hit))

        # Check we can get the original last HSP from the file.
        num_hsps = num_hits
        self.assertEqual(num_hsps, len(qresult.hsps))
        hsp = qresult.hsps[-1]

        self.assertTrue(hsp.is_included)
        self.assertEqual(num_hsps - 1, hsp.output_index)
        self.assertEqual(3.3e-05, hsp.evalue)
        self.assertEqual(51.24, hsp.score)
        self.assertEqual(96.55, hsp.prob)
        self.assertEqual(14, hsp.hit_start)
        self.assertEqual(65, hsp.hit_end)
        self.assertEqual(7, hsp.query_start)
        self.assertEqual(59, hsp.query_end)
        self.assertEqual(
            "FWLVTAALLASTVFFFVERDRVS-AKWKTSLTVSGLVTGIAFWHYMYMRGVW", hsp.hit.seq
        )
        self.assertEqual(
            "LALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGY", hsp.query.seq
        )

    def test_q9bsu1(self):
        """Parsing hhsearch_q9bsu1_uniclust_w_ss_pfamA_30.hhr file."""
        txt_file = get_file("hhsearch_q9bsu1_uniclust_w_ss_pfamA_30.hhr")
        qresults = parse(txt_file, FMT)

        # test first and only qresult
        qresult = next(qresults)

        num_hits = 12
        self.assertEqual("HHSUITE", qresult.program)
        self.assertEqual(
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70"
            " PE=1 SV=1",
            qresult.id,
        )
        self.assertEqual(422, qresult.seq_len)
        self.assertEqual(num_hits, len(qresult))

        hit = qresult[0]
        self.assertEqual("PF03676.13", hit.id)
        self.assertEqual(
            "UPF0183 ; Uncharacterised protein family (UPF0183)", hit.description
        )
        self.assertTrue(hit.is_included)
        self.assertEqual(2e-106, hit.evalue)
        self.assertEqual(822.75, hit.score)
        self.assertEqual(1, len(hit))

        hsp = hit.hsps[0]
        self.assertTrue(hsp.is_included)
        self.assertEqual(0, hsp.output_index)
        self.assertEqual(2e-106, hsp.evalue)
        self.assertEqual(822.75, hsp.score)
        self.assertEqual(100.00, hsp.prob)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(395, hsp.hit_end)
        self.assertEqual(10, hsp.query_start)
        self.assertEqual(407, hsp.query_end)
        self.assertEqual(
            "SLGNEQWEFTLGMPLAQAVAILQKHCRIIKNVQVLYSEQSPLSHDLILNLTQDGIKLMFDAFNQRLKVIEVCDLTKVKLK"
            "YCGVHFNSQAIAPTIEQIDQSFGATHPGVYNSAEQLFHLNFRGLSFSFQLDSWTEAPKYEPNFAHGLASLQIPHGATVKR"
            "MYIYSGNSLQDTKAPMMPLSCFLGNVYAESVDVLRDGTGPAGLRLRLLAAGCGPGLLADAKMRVFERSVYFGDSCQDVLS"
            "MLGSPHKVFYKSEDKMKIHSPSPHKQVPSKCNDYFFNYFTLGVDILFDANTHKVKKFVLHTNYPGHYNFNIYHRCEFKIP"
            "LAIKKENADGQTE--TCTTYSKWDNIQELLGHPVEKPVVLHRSSSPNNTNPFGSTFCFGLQRMIFEVMQNNHIASVTLY",
            hsp.query.seq,
        )
        self.assertEqual(
            "EQWE----FALGMPLAQAISILQKHCRIIKNVQVLYSEQMPLSHDLILNLTQDGIKLLFDACNQRLKVIEVYDLTKVKLK"
            "YCGVHFNSQAIAPTIEQIDQSFGATHPGVYNAAEQLFHLNFRGLSFSFQLDSWSEAPKYEPNFAHGLASLQIPHGATVKR"
            "MYIYSGNNLQETKAPAMPLACFLGNVYAECVEVLRDGAGPLGLKLRLLTAGCGPGVLADTKVRAVERSIYFGDSCQDVLS"
            "ALGSPHKVFYKSEDKMKIHSPSPHKQVPSKCNDYFFNYYILGVDILFDSTTHLVKKFVLHTNFPGHYNFNIYHRCDFKIP"
            "LIIKKDGADAHSEDCILTTYSKWDQIQELLGHPMEKPVVLHRSSSANNTNPFGSTFCFGLQRMIFEVMQNNHIASVTLY",
            hsp.hit.seq,
        )

        # Check last hit
        hit = qresult[num_hits - 1]
        self.assertEqual("PF10049.8", hit.id)
        self.assertEqual(
            "DUF2283 ; Protein of unknown function (DUF2283)", hit.description
        )
        self.assertTrue(hit.is_included)
        self.assertEqual(78, hit.evalue)
        self.assertEqual(19.81, hit.score)
        self.assertEqual(1, len(hit))

        # Check we can get the original last HSP from the file.
        num_hsps = 16
        self.assertEqual(num_hsps, len(qresult.hsps))

        hsp = qresult.hsps[-1]
        self.assertTrue(hsp.is_included)
        self.assertEqual(num_hsps - 1, hsp.output_index)
        self.assertEqual(78, hsp.evalue)
        self.assertEqual(19.81, hsp.score)
        self.assertEqual(20.88, hsp.prob)
        self.assertEqual(25, hsp.hit_start)
        self.assertEqual(48, hsp.hit_end)
        self.assertEqual(61, hsp.query_start)
        self.assertEqual(85, hsp.query_end)
        self.assertEqual("APNVIFDYDA-EGRIVGIELLDAR", hsp.hit.seq)
        self.assertEqual("QDGIKLMFDAFNQRLKVIEVCDLT", hsp.query.seq)

    def test_4p79(self):
        """Parsing 4p79_hhsearch_server_NOssm.hhr file."""
        txt_file = get_file("4p79_hhsearch_server_NOssm.hhr")
        qresults = parse(txt_file, FMT)

        # test first and only qresult
        qresult = next(qresults)

        num_hits = 8
        self.assertEqual("HHSUITE", qresult.program)
        self.assertEqual("4P79:A|PDBID|CHAIN|SEQUENCE", qresult.id)
        self.assertEqual(198, qresult.seq_len)
        self.assertEqual(num_hits, len(qresult))

        hit = qresult[0]
        self.assertEqual("4P79_A", hit.id)
        self.assertEqual(
            "cell adhesion protein; cell adhesion, tight junction, membrane; HET: OLC"
            ", MSE; 2.4A {Mus musculus}",
            hit.description,
        )
        self.assertTrue(hit.is_included)
        self.assertEqual(6.8e-32, hit.evalue)
        self.assertEqual(194.63, hit.score)
        self.assertEqual(1, len(hit))

        hsp = hit.hsps[0]
        self.assertTrue(hsp.is_included)
        self.assertEqual(0, hsp.output_index)
        self.assertEqual(6.8e-32, hsp.evalue)
        self.assertEqual(194.63, hsp.score)
        self.assertEqual(99.94, hsp.prob)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(198, hsp.hit_end)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(198, hsp.query_end)
        self.assertEqual(
            "GSEFMSVAVETFGFFMSALGLLMLGLTLSNSYWRVSTVHGNVITTNTIFENLWYSCATDSLGVSNCWDFPSMLALSGYVQ"
            "GCRALMITAILLGFLGLFLGMVGLRATNVGNMDLSKKAKLLAIAGTLHILAGACGMVAISWYAVNITTDFFNPLYAGTKY"
            "ELGPALYLGWSASLLSILGGICVFSTAAASSKEEPATR",
            hsp.query.seq,
        )
        self.assertEqual(
            "GSEFMSVAVETFGFFMSALGLLMLGLTLSNSYWRVSTVHGNVITTNTIFENLWYSCATDSLGVSNCWDFPSMLALSGYVQ"
            "GCRALMITAILLGFLGLFLGMVGLRATNVGNMDLSKKAKLLAIAGTLHILAGACGMVAISWYAVNITTDFFNPLYAGTKY"
            "ELGPALYLGWSASLLSILGGICVFSTAAASSKEEPATR",
            hsp.hit.seq,
        )

        # Check last hit
        hit = qresult[num_hits - 1]
        self.assertEqual("5YQ7_F", hit.id)
        self.assertEqual(
            "Beta subunit of light-harvesting 1; Photosynthetic core complex, PHOTOSYNTHESIS; "
            "HET: MQE, BCL, HEM, KGD, BPH;{Roseiflexus castenholzii}; Related PDB entries: 5YQ7_V"
            " 5YQ7_3 5YQ7_T 5YQ7_J 5YQ7_9 5YQ7_N 5YQ7_A 5YQ7_P 5YQ7_H 5YQ7_D 5YQ7_5 5YQ7_7 5YQ7_1 "
            "5YQ7_R",
            hit.description,
        )
        self.assertTrue(hit.is_included)
        self.assertEqual(6.7, hit.evalue)
        self.assertEqual(20.51, hit.score)
        self.assertEqual(1, len(hit))

        # Check we can get the original last HSP from the file.
        num_hsps = num_hits
        self.assertEqual(num_hsps, len(qresult.hsps))

        hsp = qresult.hsps[-1]
        self.assertTrue(hsp.is_included)
        self.assertEqual(num_hsps - 1, hsp.output_index)
        self.assertEqual(6.7, hsp.evalue)
        self.assertEqual(20.51, hsp.score)
        self.assertEqual(52.07, hsp.prob)
        self.assertEqual(8, hsp.hit_start)
        self.assertEqual(42, hsp.hit_end)
        self.assertEqual(5, hsp.query_start)
        self.assertEqual(37, hsp.query_end)
        self.assertEqual("RTSVVVSTLLGLVMALLIHFVVLSSGAFNWLRAP", hsp.hit.seq)
        self.assertEqual("SVAVETFGFFMSALGLLMLGLTLSNS--YWRVST", hsp.query.seq)

    def test_9590198(self):
        """Parsing hhpred_9590198.hhr file."""
        txt_file = get_file("hhpred_9590198.hhr")
        qresults = parse(txt_file, FMT)

        # test first and only qresult
        qresult = next(qresults)

        num_hits = 22
        self.assertEqual("HHSUITE", qresult.program)
        self.assertEqual(
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70"
            " PE=1 SV=1",
            qresult.id,
        )
        self.assertEqual(422, qresult.seq_len)
        self.assertEqual(num_hits, len(qresult))

        hit = qresult[0]
        self.assertEqual("PF03676.14", hit.id)
        self.assertEqual(
            "UPF0183 ; Uncharacterised protein family (UPF0183)", hit.description
        )
        self.assertTrue(hit.is_included)
        self.assertEqual(9.9e-102, hit.evalue)
        self.assertEqual(792.76, hit.score)
        self.assertEqual(1, len(hit))

        hsp = hit.hsps[0]
        self.assertTrue(hsp.is_included)
        self.assertEqual(0, hsp.output_index)
        self.assertEqual(9.9e-102, hsp.evalue)
        self.assertEqual(792.76, hsp.score)
        self.assertEqual(100.00, hsp.prob)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(394, hsp.hit_end)
        self.assertEqual(21, hsp.query_start)
        self.assertEqual(407, hsp.query_end)
        self.assertEqual(
            "GMHFSQSVAIIQSQVGTIRGVQVLYSDQNPLSVDLVINMPQDGMRLIFDPVAQRLKIIEIYNMKLVKLRYSGMCFNSPEI"
            "TPSIEQVEHCFGATHPGLYDSQRHLFALNFRGLSFYFPVDS-----KFEPGYAHGLGSLQFPNGGSPVVSRTTIYYGSQH"
            "QLSSNTSSRVSGVPLPDLPLSCYRQQLHLRRCDVLRNTTSTMGLRLHMFTEGT--SRALEPSQVALVRVVRFGDSCQGVA"
            "RALGAPARLYYKADDKMRIHRPTARRR-PPPASDYLFNYFTLGLDVLFDARTNQVKKFVLHTNYPGHYNFNMYHRCEFEL"
            "TVQPD-KSEAHSLVESGGGVAVTAYSKWEVVSRAL-RVCERPVVLNRASSTNTTNPFGSTFCYGYQDIIFEVMSNNYIAS"
            "ITLY",
            hsp.hit.seq,
        )
        self.assertEqual(
            "GMPLAQAVAILQKHCRIIKNVQVLYSEQSPLSHDLILNLTQDGIKLMFDAFNQRLKVIEVCDLTKVKLKYCGVHFNSQAI"
            "APTIEQIDQSFGATHPGVYNSAEQLFHLNFRGLSFSFQLDSWTEAPKYEPNFAHGLASLQIPHGA--TVKRMYIYSGNSL"
            "Q---------DTKA-PMMPLSCFLGNVYAESVDVLRDGTGPAGLRLRLLAAGCGPGLLADAKMRVFERSVYFGDSCQDVL"
            "SMLGSPHKVFYKSEDKMKIHSPSPHKQVPSKCNDYFFNYFTLGVDILFDANTHKVKKFVLHTNYPGHYNFNIYHRCEFKI"
            "PLAIKKENADG------QTETCTTYSKWDNIQELLGHPVEKPVVLHRSSSPNNTNPFGSTFCFGLQRMIFEVMQNNHIAS"
            "VTLY",
            hsp.query.seq,
        )

        # Check last hit
        hit = qresult[num_hits - 1]
        self.assertEqual("4IL7_A", hit.id)
        self.assertEqual(
            "Putative uncharacterized protein; partial jelly roll fold, hypothetical; 1.4A "
            "{Sulfolobus turreted icosahedral virus}",
            hit.description,
        )
        self.assertTrue(hit.is_included)
        self.assertEqual(6.8e02, hit.evalue)
        self.assertEqual(22.72, hit.score)
        self.assertEqual(1, len(hit))

        # Check we can get the original last HSP from the file.
        num_hsps = 34
        self.assertEqual(num_hsps, len(qresult.hsps))

        hsp = qresult.hsps[-1]
        self.assertTrue(hsp.is_included)
        self.assertEqual(num_hsps - 1, hsp.output_index)
        self.assertEqual(3.9e02, hsp.evalue)
        self.assertEqual(22.84, hsp.score)
        self.assertEqual(21.56, hsp.prob)
        self.assertEqual(7, hsp.hit_start)
        self.assertEqual(96, hsp.hit_end)
        self.assertEqual(18, hsp.query_start)
        self.assertEqual(114, hsp.query_end)
        self.assertEqual(
            "FTLGMPLAQAVAILQKHCRIIKNVQVLYSEQSPLSHDLILNLTQDGIKLMFDAFNQRLKVIEVCDLTKVKLKYCGVH-FN"
            "SQAIAPTIEQIDQSFGA",
            hsp.query.seq,
        )
        self.assertEqual(
            "IQFGMDRTLVWQLAGADQSCSDQVERIICYNNPDH-------YGPQGHFFFNA-ADKLIHKRQMELFPAPKPTMRLATYN"
            "KTQTGMTEAQFWAAVPS",
            hsp.hit.seq,
        )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
