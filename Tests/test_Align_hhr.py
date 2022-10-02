# Copyright 2020 by Michiel de Hoon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Bio.Align.hhr module."""
import os
import unittest


from Bio import Align

import numpy


class Align_hhr_2uvo_hhblits(unittest.TestCase):
    path = os.path.join("HHsuite", "2uvo_hhblits.hhr")

    def test_reading(self):
        alignments = Align.parse(self.path, "hhr")
        self.assertEqual(alignments.metadata["No_of_seqs"], (1560, 4005))
        self.assertAlmostEqual(alignments.metadata["Neff"], 8.3)
        self.assertEqual(alignments.metadata["Searched_HMMs"], 34)
        self.assertEqual(alignments.metadata["Rundate"], "Fri Feb 15 16:34:13 2019")
        self.assertEqual(
            alignments.metadata["Command line"], "hhblits -i 2uvoAh.fasta -d /pdb70"
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 99.95)
        self.assertAlmostEqual(alignment.annotations["E-value"], 3.7e-34)
        self.assertAlmostEqual(alignment.annotations["Score"], 210.31)
        self.assertAlmostEqual(alignment.annotations["Identities"], 100)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 2.050)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 166.9)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[0:171],
            "ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYCGAGCQSGGCDG",
        )
        self.assertEqual(
            alignment[1],
            "ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYCGAGCQSGGCDG",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "~~cg~~~~~~~c~~~~CCs~~g~CG~~~~~c~~~c~~~~c~~~~~Cg~~~~~~~c~~~~CCs~~g~CG~~~~~c~~~c~~~~~~~~~~Cg~~~~~~~c~~~~CCS~~G~CG~~~~~C~~~Cq~~~c~~~~~Cg~~~~~~~c~~~~CCS~~G~CG~~~~~C~~gCq~~~c~~",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "~~cg~~~~~~~c~~~~CCS~~g~Cg~~~~~Cg~gC~~~~c~~~~~cg~~~~~~~c~~~~CCs~~g~Cg~~~~~c~~~c~~~~~~~~~~cg~~~~~~~c~~~~CCs~~g~CG~~~~~C~~gCq~~~c~~~~~cg~~~~~~~c~~~~ccs~~g~Cg~~~~~C~~~cq~~~~~~",
        )
        self.assertEqual(alignment.target.id, "2uvo_A")
        self.assertEqual(
            alignment.target.seq[0:171],
            "ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYCGAGCQSGGCDG",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "2uvo_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Agglutinin isolectin 1; carbohydrate-binding protein, hevein domain, chitin-binding, GERM agglutinin, chitin-binding protein; HET: NDG NAG GOL; 1.40A {Triticum aestivum} PDB: 1wgc_A* 2cwg_A* 2x3t_A* 4aml_A* 7wga_A 9wga_A 2wgc_A 1wgt_A 1k7t_A* 1k7v_A* 1k7u_A 2x52_A* 1t0w_A*",
        )
        self.assertEqual(
            alignment[0],
            "ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYCGAGCQSGGCDG",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "CBCBGGGTTBBCGGGCEECTTSBEEBSHHHHSTTCCBSSCSSCCBCBGGGTTBCCSTTCEECTTSBEEBSHHHHSTTCCBSSCSSCCBCBGGGTTBCCGGGCEECTTSBEEBSHHHHSTTCCBSSCSSCCCCBTTTTTBCCSTTCEECTTSCEEBSHHHHSTTCCBSSCC ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "CCCCCCCCCcCCCCCCeeCCCCeECCCcccccCCccccccccccccCcccCCcccCCccccCCCceeCCCccccCCCcccccccccccccccccCCCCCCCcccCCCCccCCCcccccCCCcCCccccccccccccccccCCCCCCcCCCCEecCchhhcccccccCCCCC",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            "799999999999999999999999999999999999999999999999999999999999999999999999999999999999999999998899999999999999999999999999999999999999999999999999999999999999999999999999986",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  0, 171],
                             [  0, 171]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['E', 'R', 'C', 'G', 'E', 'Q', 'G', 'S', 'N', 'M', 'E', 'C', 'P',
              'N', 'N', 'L', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'M',
              'G', 'G', 'D', 'Y', 'C', 'G', 'K', 'G', 'C', 'Q', 'N', 'G', 'A',
              'C', 'W', 'T', 'S', 'K', 'R', 'C', 'G', 'S', 'Q', 'A', 'G', 'G',
              'A', 'T', 'C', 'T', 'N', 'N', 'Q', 'C', 'C', 'S', 'Q', 'Y', 'G',
              'Y', 'C', 'G', 'F', 'G', 'A', 'E', 'Y', 'C', 'G', 'A', 'G', 'C',
              'Q', 'G', 'G', 'P', 'C', 'R', 'A', 'D', 'I', 'K', 'C', 'G', 'S',
              'Q', 'A', 'G', 'G', 'K', 'L', 'C', 'P', 'N', 'N', 'L', 'C', 'C',
              'S', 'Q', 'W', 'G', 'F', 'C', 'G', 'L', 'G', 'S', 'E', 'F', 'C',
              'G', 'G', 'G', 'C', 'Q', 'S', 'G', 'A', 'C', 'S', 'T', 'D', 'K',
              'P', 'C', 'G', 'K', 'D', 'A', 'G', 'G', 'R', 'V', 'C', 'T', 'N',
              'N', 'Y', 'C', 'C', 'S', 'K', 'W', 'G', 'S', 'C', 'G', 'I', 'G',
              'P', 'G', 'Y', 'C', 'G', 'A', 'G', 'C', 'Q', 'S', 'G', 'G', 'C',
              'D', 'G'],
             ['E', 'R', 'C', 'G', 'E', 'Q', 'G', 'S', 'N', 'M', 'E', 'C', 'P',
              'N', 'N', 'L', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'M',
              'G', 'G', 'D', 'Y', 'C', 'G', 'K', 'G', 'C', 'Q', 'N', 'G', 'A',
              'C', 'W', 'T', 'S', 'K', 'R', 'C', 'G', 'S', 'Q', 'A', 'G', 'G',
              'A', 'T', 'C', 'T', 'N', 'N', 'Q', 'C', 'C', 'S', 'Q', 'Y', 'G',
              'Y', 'C', 'G', 'F', 'G', 'A', 'E', 'Y', 'C', 'G', 'A', 'G', 'C',
              'Q', 'G', 'G', 'P', 'C', 'R', 'A', 'D', 'I', 'K', 'C', 'G', 'S',
              'Q', 'A', 'G', 'G', 'K', 'L', 'C', 'P', 'N', 'N', 'L', 'C', 'C',
              'S', 'Q', 'W', 'G', 'F', 'C', 'G', 'L', 'G', 'S', 'E', 'F', 'C',
              'G', 'G', 'G', 'C', 'Q', 'S', 'G', 'A', 'C', 'S', 'T', 'D', 'K',
              'P', 'C', 'G', 'K', 'D', 'A', 'G', 'G', 'R', 'V', 'C', 'T', 'N',
              'N', 'Y', 'C', 'C', 'S', 'K', 'W', 'G', 'S', 'C', 'G', 'I', 'G',
              'P', 'G', 'Y', 'C', 'G', 'A', 'G', 'C', 'Q', 'S', 'G', 'G', 'C',
              'D', 'G']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 99.92)
        self.assertAlmostEqual(alignment.annotations["E-value"], 1.1e-30)
        self.assertAlmostEqual(alignment.annotations["Score"], 190.44)
        self.assertAlmostEqual(alignment.annotations["Identities"], 49)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.254)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 148.8)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[1:169],
            "RCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYCGAGCQSGGC",
        )
        self.assertEqual(
            alignment[1],
            "RCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYCGAGCQSGGC",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            " ~cg~~~~~~~c~~~~CCs~~g~CG~~~~~c~~~c~~~~c~~~~~Cg~~~~~~~c~~~~CCs~~g~CG~~~~~c~~~c~~~~~~~~~~Cg~~~~~~~c~~~~CCS~~G~CG~~~~~C~~~Cq~~~c~~~~~Cg~~~~~~~c~~~~CCS~~G~CG~~~~~C~~gCq~~~c  ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            " ~~g~~~~~~~c~~~~CCs~~g~Cg~~~~~c~~~C~~~~c~~~~~cg~~~~~c~~~~CCs~~G~CG~~~~~C~~~C~~~~~~~~~~Cg~~~~~c~~~~Ccs~~G~CGt~~~~C~~~cq~~~c~~~~~cg~~~~~c~~~~Ccs~~g~Cg~~~~~C~~~cq~~~~ ",
        )
        self.assertEqual(alignment.target.id, "2wga")
        self.assertEqual(
            alignment.target.seq[1:163],
            "GXGCXGXXMYCSTNNCCXXWESCGSGGYXCGEGCNLGACQXGXPCXXPGSVCTNLHCCARGGHCGMGSGYCGXGCXGGACXADIXCGXGXXXCPTDSCCGGWGXCGNGXEFCGXGCXVGGCAAXSPCGXPGSXCTLDKCCSGXGACXSGSGGCGXGCXAGGC",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "2wga")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "; lectin (agglutinin); NMR {}",
        )
        self.assertEqual(
            alignment[0],
            "GXGCXGXXMYCSTNNCCXXWESCGSGGYXCGEGCNLGACQXGXPCXX--PGSVCTNLHCCARGGHCGMGSGYCGXGCXGGACXADIXCGXG--XXXCPTDSCCGGWGXCGNGXEFCGXGCXVGGCAAXSPCGXP--GSXCTLDKCCSGXGACXSGSGGCGXGCXAGGC",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            " CSGGGSSCCCCSTTCEECTTSCEECSTTTTSTTCCSSSCSSCCCSSSSSCCCSTTCEECTTSCEESSHHHHSSCCSSSSCSSCCCCTTSSSCCSTTCBCCSSSCCBCSHHHHSTTCCSSSCSSCCCCSSSCCCCSTTCEECSSSSEECSTTTTSSCCSSSSC ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            " ccCCCcccccCCCCceECCcceECCCCccccCccccCcccccceeccCCCcCCCCcccCCCceeCCCCcccCCCccccccccccccCcccccCCCCCccCCCCCccCccccccCCccccccccccccCCCcccCCcccccCCCCceeCCccccCCCCcCCCC ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            " 699999999999999999999999999999999999888888899873567889999999999999999999999988888888999744568899999999999999999999999999998899987356788999999999999999999999999876 ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  1,  48,  48,  90,  90, 131, 131, 163],
                             [  1,  48,  50,  92,  94, 135, 137, 169]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['G', 'X', 'G', 'C', 'X', 'G', 'X', 'X', 'M', 'Y', 'C', 'S', 'T',
              'N', 'N', 'C', 'C', 'X', 'X', 'W', 'E', 'S', 'C', 'G', 'S', 'G',
              'G', 'Y', 'X', 'C', 'G', 'E', 'G', 'C', 'N', 'L', 'G', 'A', 'C',
              'Q', 'X', 'G', 'X', 'P', 'C', 'X', 'X', '-', '-', 'P', 'G', 'S',
              'V', 'C', 'T', 'N', 'L', 'H', 'C', 'C', 'A', 'R', 'G', 'G', 'H',
              'C', 'G', 'M', 'G', 'S', 'G', 'Y', 'C', 'G', 'X', 'G', 'C', 'X',
              'G', 'G', 'A', 'C', 'X', 'A', 'D', 'I', 'X', 'C', 'G', 'X', 'G',
              '-', '-', 'X', 'X', 'X', 'C', 'P', 'T', 'D', 'S', 'C', 'C', 'G',
              'G', 'W', 'G', 'X', 'C', 'G', 'N', 'G', 'X', 'E', 'F', 'C', 'G',
              'X', 'G', 'C', 'X', 'V', 'G', 'G', 'C', 'A', 'A', 'X', 'S', 'P',
              'C', 'G', 'X', 'P', '-', '-', 'G', 'S', 'X', 'C', 'T', 'L', 'D',
              'K', 'C', 'C', 'S', 'G', 'X', 'G', 'A', 'C', 'X', 'S', 'G', 'S',
              'G', 'G', 'C', 'G', 'X', 'G', 'C', 'X', 'A', 'G', 'G', 'C'],
             ['R', 'C', 'G', 'E', 'Q', 'G', 'S', 'N', 'M', 'E', 'C', 'P', 'N',
              'N', 'L', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'M', 'G',
              'G', 'D', 'Y', 'C', 'G', 'K', 'G', 'C', 'Q', 'N', 'G', 'A', 'C',
              'W', 'T', 'S', 'K', 'R', 'C', 'G', 'S', 'Q', 'A', 'G', 'G', 'A',
              'T', 'C', 'T', 'N', 'N', 'Q', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y',
              'C', 'G', 'F', 'G', 'A', 'E', 'Y', 'C', 'G', 'A', 'G', 'C', 'Q',
              'G', 'G', 'P', 'C', 'R', 'A', 'D', 'I', 'K', 'C', 'G', 'S', 'Q',
              'A', 'G', 'G', 'K', 'L', 'C', 'P', 'N', 'N', 'L', 'C', 'C', 'S',
              'Q', 'W', 'G', 'F', 'C', 'G', 'L', 'G', 'S', 'E', 'F', 'C', 'G',
              'G', 'G', 'C', 'Q', 'S', 'G', 'A', 'C', 'S', 'T', 'D', 'K', 'P',
              'C', 'G', 'K', 'D', 'A', 'G', 'G', 'R', 'V', 'C', 'T', 'N', 'N',
              'Y', 'C', 'C', 'S', 'K', 'W', 'G', 'S', 'C', 'G', 'I', 'G', 'P',
              'G', 'Y', 'C', 'G', 'A', 'G', 'C', 'Q', 'S', 'G', 'G', 'C']],
            dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 99.77)
        self.assertAlmostEqual(alignment.annotations["E-value"], 5.2e-24)
        self.assertAlmostEqual(alignment.annotations["Score"], 148.19)
        self.assertAlmostEqual(alignment.annotations["Identities"], 50)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.303)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 103.7)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[0:124],
            "ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSG",
        )
        self.assertEqual(
            alignment[1],
            "ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSG",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "~~cg~~~~~~~c~~~~CCs~~g~CG~~~~~c~~~c~~~~c~~~~~Cg~~~~~~~c~~~~CCs~~g~CG~~~~~c~~~c~~~~~~~~~~Cg~~~~~~~c~~~~CCS~~G~CG~~~~~C~~~Cq~~                                               ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            " ~~cg~~~~~~~c~~g~CCs~~g~CG~~~~~Cg~gCq~~c~~~~Cg~~~~~~~c~~~~CCs~~G~CG~~~~~C~~~Cq~~c~~~~cg~~~~~~~c~~~~CCs~~g~CG~~~~~C~~gCq~~     ",
        )
        self.assertEqual(alignment.target.id, "1ulk_A")
        self.assertEqual(
            alignment.target.seq[1:121],
            "PVCGVRASGRVCPDGYCCSQWGYCGTTEEYCGKGCQSQCDYNRCGKEFGGKECHDELCCSQYGWCGNSDGHCGEGCQSQCSYWRCGKDFGGRLCTEDMCCSQYGWCGLTDDHCEDGCQSQ",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "1ulk_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Lectin-C; chitin-binding protein, hevein domain, PL-C, sugar binding protein; 1.80A {Phytolacca americana} SCOP: g.3.1.1 g.3.1.1 g.3.1.1",
        )
        self.assertEqual(
            alignment[0],
            "PVCGVRASGRVCPDGYCCSQWGYCGTTEEYCGKGCQSQ-C-DYNRCGKEFGGKECHDELCCSQYGWCGNSDGHCGEGCQSQ-C-SYWRCGKDFGGRLCTEDMCCSQYGWCGLTDDHCEDGCQSQ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            " CCSBGGGTTBCCGGGCEECTTSCEESSHHHHSTTCCBCTTTTBCBGGGTTBCCGGGCEECTTSBEECSHHHHSTTCCBCTTTTBCBGGGTTBCCSTTCEECTTSBEECSHHHHSTTCCBC     ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            " CcCcCCCCCCCCCCCCeECCCCeeCCCccccCCCccccceeeeccccccCCCCCCccccCCCcccccCcccccCCcccccCcccccccCCCccCCCCcccccCceecCcccccCcccccc     ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            " 389999999999999999999999999999999999752346788765567888999999999999999999999986422467877554567899999999999999999999999973     ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  1,  39,  39,  40,  40,  80,  80,  81,  81, 121],
                             [  0,  38,  39,  40,  41,  81,  82,  83,  84, 124]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['P', 'V', 'C', 'G', 'V', 'R', 'A', 'S', 'G', 'R', 'V', 'C', 'P',
              'D', 'G', 'Y', 'C', 'C', 'S', 'Q', 'W', 'G', 'Y', 'C', 'G', 'T',
              'T', 'E', 'E', 'Y', 'C', 'G', 'K', 'G', 'C', 'Q', 'S', 'Q', '-',
              'C', '-', 'D', 'Y', 'N', 'R', 'C', 'G', 'K', 'E', 'F', 'G', 'G',
              'K', 'E', 'C', 'H', 'D', 'E', 'L', 'C', 'C', 'S', 'Q', 'Y', 'G',
              'W', 'C', 'G', 'N', 'S', 'D', 'G', 'H', 'C', 'G', 'E', 'G', 'C',
              'Q', 'S', 'Q', '-', 'C', '-', 'S', 'Y', 'W', 'R', 'C', 'G', 'K',
              'D', 'F', 'G', 'G', 'R', 'L', 'C', 'T', 'E', 'D', 'M', 'C', 'C',
              'S', 'Q', 'Y', 'G', 'W', 'C', 'G', 'L', 'T', 'D', 'D', 'H', 'C',
              'E', 'D', 'G', 'C', 'Q', 'S', 'Q'],
             ['E', 'R', 'C', 'G', 'E', 'Q', 'G', 'S', 'N', 'M', 'E', 'C', 'P',
              'N', 'N', 'L', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'M',
              'G', 'G', 'D', 'Y', 'C', 'G', 'K', 'G', 'C', 'Q', 'N', 'G', 'A',
              'C', 'W', 'T', 'S', 'K', 'R', 'C', 'G', 'S', 'Q', 'A', 'G', 'G',
              'A', 'T', 'C', 'T', 'N', 'N', 'Q', 'C', 'C', 'S', 'Q', 'Y', 'G',
              'Y', 'C', 'G', 'F', 'G', 'A', 'E', 'Y', 'C', 'G', 'A', 'G', 'C',
              'Q', 'G', 'G', 'P', 'C', 'R', 'A', 'D', 'I', 'K', 'C', 'G', 'S',
              'Q', 'A', 'G', 'G', 'K', 'L', 'C', 'P', 'N', 'N', 'L', 'C', 'C',
              'S', 'Q', 'W', 'G', 'F', 'C', 'G', 'L', 'G', 'S', 'E', 'F', 'C',
              'G', 'G', 'G', 'C', 'Q', 'S', 'G']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 99.77)
        self.assertAlmostEqual(alignment.annotations["E-value"], 6.8e-24)
        self.assertAlmostEqual(alignment.annotations["Score"], 147.59)
        self.assertAlmostEqual(alignment.annotations["Identities"], 52)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.251)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 102.2)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[43:167],
            "KRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYCGAGCQSG",
        )
        self.assertEqual(
            alignment[1],
            "KRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYCGAGCQSG",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                           ~~Cg~~~~~~~c~~~~CCs~~g~CG~~~~~c~~~c~~~~~~~~~~Cg~~~~~~~c~~~~CCS~~G~CG~~~~~C~~~Cq~~~c~~~~~Cg~~~~~~~c~~~~CCS~~G~CG~~~~~C~~gCq~~    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            " ~~cg~~~~~~~c~~g~CCs~~g~CG~~~~~Cg~gCq~~c~~~~Cg~~~~~~~c~~~~CCs~~G~CG~~~~~C~~~Cq~~c~~~~cg~~~~~~~c~~~~CCs~~g~CG~~~~~C~~gCq~~     ",
        )
        self.assertEqual(alignment.target.id, "1ulk_A")
        self.assertEqual(
            alignment.target.seq[1:121],
            "PVCGVRASGRVCPDGYCCSQWGYCGTTEEYCGKGCQSQCDYNRCGKEFGGKECHDELCCSQYGWCGNSDGHCGEGCQSQCSYWRCGKDFGGRLCTEDMCCSQYGWCGLTDDHCEDGCQSQ",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "1ulk_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Lectin-C; chitin-binding protein, hevein domain, PL-C, sugar binding protein; 1.80A {Phytolacca americana} SCOP: g.3.1.1 g.3.1.1 g.3.1.1",
        )
        self.assertEqual(
            alignment[0],
            "PVCGVRASGRVCPDGYCCSQWGYCGTTEEYCGKGCQSQ-C-DYNRCGKEFGGKECHDELCCSQYGWCGNSDGHCGEGCQSQ-C-SYWRCGKDFGGRLCTEDMCCSQYGWCGLTDDHCEDGCQSQ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            " CCSBGGGTTBCCGGGCEECTTSCEESSHHHHSTTCCBCTTTTBCBGGGTTBCCGGGCEECTTSBEECSHHHHSTTCCBCTTTTBCBGGGTTBCCSTTCEECTTSBEECSHHHHSTTCCBC     ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            " CcCcCCCCCCCCCCCCeECCCCeeCCCccccCCCccccceeeeccccccCCCCCCccccCCCcccccCcccccCCcccccCcccccccCCCccCCCCcccccCceecCcccccCcccccc     ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            " 468888778999999999999999999999999998642345777654456788899999999999999999999997422467787555668999999999999999999999999984     ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  1,  39,  39,  40,  40,  80,  80,  81,  81, 121],
                             [ 43,  81,  82,  83,  84, 124, 125, 126, 127, 167]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['P', 'V', 'C', 'G', 'V', 'R', 'A', 'S', 'G', 'R', 'V', 'C', 'P',
              'D', 'G', 'Y', 'C', 'C', 'S', 'Q', 'W', 'G', 'Y', 'C', 'G', 'T',
              'T', 'E', 'E', 'Y', 'C', 'G', 'K', 'G', 'C', 'Q', 'S', 'Q', '-',
              'C', '-', 'D', 'Y', 'N', 'R', 'C', 'G', 'K', 'E', 'F', 'G', 'G',
              'K', 'E', 'C', 'H', 'D', 'E', 'L', 'C', 'C', 'S', 'Q', 'Y', 'G',
              'W', 'C', 'G', 'N', 'S', 'D', 'G', 'H', 'C', 'G', 'E', 'G', 'C',
              'Q', 'S', 'Q', '-', 'C', '-', 'S', 'Y', 'W', 'R', 'C', 'G', 'K',
              'D', 'F', 'G', 'G', 'R', 'L', 'C', 'T', 'E', 'D', 'M', 'C', 'C',
              'S', 'Q', 'Y', 'G', 'W', 'C', 'G', 'L', 'T', 'D', 'D', 'H', 'C',
              'E', 'D', 'G', 'C', 'Q', 'S', 'Q'],
             ['K', 'R', 'C', 'G', 'S', 'Q', 'A', 'G', 'G', 'A', 'T', 'C', 'T',
              'N', 'N', 'Q', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'F',
              'G', 'A', 'E', 'Y', 'C', 'G', 'A', 'G', 'C', 'Q', 'G', 'G', 'P',
              'C', 'R', 'A', 'D', 'I', 'K', 'C', 'G', 'S', 'Q', 'A', 'G', 'G',
              'K', 'L', 'C', 'P', 'N', 'N', 'L', 'C', 'C', 'S', 'Q', 'W', 'G',
              'F', 'C', 'G', 'L', 'G', 'S', 'E', 'F', 'C', 'G', 'G', 'G', 'C',
              'Q', 'S', 'G', 'A', 'C', 'S', 'T', 'D', 'K', 'P', 'C', 'G', 'K',
              'D', 'A', 'G', 'G', 'R', 'V', 'C', 'T', 'N', 'N', 'Y', 'C', 'C',
              'S', 'K', 'W', 'G', 'S', 'C', 'G', 'I', 'G', 'P', 'G', 'Y', 'C',
              'G', 'A', 'G', 'C', 'Q', 'S', 'G']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 99.72)
        self.assertAlmostEqual(alignment.annotations["E-value"], 9.7e-23)
        self.assertAlmostEqual(alignment.annotations["Score"], 148.46)
        self.assertAlmostEqual(alignment.annotations["Identities"], 61)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.542)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 111.6)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[44:169],
            "RCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYCGAGCQSGGC",
        )
        self.assertEqual(
            alignment[1],
            "RCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYCGAGCQSGGC",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                            ~Cg~~~~~~~c~~~~CCs~~g~CG~~~~~c~~~c~~~~~~~~~~Cg~~~~~~~c~~~~CCS~~G~CG~~~~~C~~~Cq~~~c~~~~~Cg~~~~~~~c~~~~CCS~~G~CG~~~~~C~~gCq~~~c  ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            " ~cg~~~~~~~c~~~~CCS~~g~Cg~~~~~Cg~gC~~~~c~~~~~cg~~~~~~~c~~~~CCs~~g~Cg~~~~~c~~~c~~~~~~~~~~cg~~~~~~~c~~~~CCs~~g~CG~~~~~C~~gCq~~~c                                             ",
        )
        self.assertEqual(alignment.target.id, "2uvo_A")
        self.assertEqual(
            alignment.target.seq[1:126],
            "RCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGAC",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "2uvo_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Agglutinin isolectin 1; carbohydrate-binding protein, hevein domain, chitin-binding, GERM agglutinin, chitin-binding protein; HET: NDG NAG GOL; 1.40A {Triticum aestivum} PDB: 1wgc_A* 2cwg_A* 2x3t_A* 4aml_A* 7wga_A 9wga_A 2wgc_A 1wgt_A 1k7t_A* 1k7v_A* 1k7u_A 2x52_A* 1t0w_A*",
        )
        self.assertEqual(
            alignment[0],
            "RCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGAC",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            " BCBGGGTTBBCGGGCEECTTSBEEBSHHHHSTTCCBSSCSSCCBCBGGGTTBCCSTTCEECTTSBEEBSHHHHSTTCCBSSCSSCCBCBGGGTTBCCGGGCEECTTSBEEBSHHHHSTTCCBSSC                                             ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            " CCCCCCCCcCCCCCCeeCCCCeECCCcccccCCccccccccccccCcccCCcccCCccccCCCceeCCCccccCCCcccccccccccccccccCCCCCCCcccCCCCccCCCcccccCCCcCCcc                                             ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            " 57776556789999999999999999999999999888888899998654567889999999999999999999999988888889999754456789999999999999999999999998755                                             ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  1, 126],
                             [ 44, 169]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['R', 'C', 'G', 'E', 'Q', 'G', 'S', 'N', 'M', 'E', 'C', 'P', 'N',
              'N', 'L', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'M', 'G',
              'G', 'D', 'Y', 'C', 'G', 'K', 'G', 'C', 'Q', 'N', 'G', 'A', 'C',
              'W', 'T', 'S', 'K', 'R', 'C', 'G', 'S', 'Q', 'A', 'G', 'G', 'A',
              'T', 'C', 'T', 'N', 'N', 'Q', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y',
              'C', 'G', 'F', 'G', 'A', 'E', 'Y', 'C', 'G', 'A', 'G', 'C', 'Q',
              'G', 'G', 'P', 'C', 'R', 'A', 'D', 'I', 'K', 'C', 'G', 'S', 'Q',
              'A', 'G', 'G', 'K', 'L', 'C', 'P', 'N', 'N', 'L', 'C', 'C', 'S',
              'Q', 'W', 'G', 'F', 'C', 'G', 'L', 'G', 'S', 'E', 'F', 'C', 'G',
              'G', 'G', 'C', 'Q', 'S', 'G', 'A', 'C'],
             ['R', 'C', 'G', 'S', 'Q', 'A', 'G', 'G', 'A', 'T', 'C', 'T', 'N',
              'N', 'Q', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'F', 'G',
              'A', 'E', 'Y', 'C', 'G', 'A', 'G', 'C', 'Q', 'G', 'G', 'P', 'C',
              'R', 'A', 'D', 'I', 'K', 'C', 'G', 'S', 'Q', 'A', 'G', 'G', 'K',
              'L', 'C', 'P', 'N', 'N', 'L', 'C', 'C', 'S', 'Q', 'W', 'G', 'F',
              'C', 'G', 'L', 'G', 'S', 'E', 'F', 'C', 'G', 'G', 'G', 'C', 'Q',
              'S', 'G', 'A', 'C', 'S', 'T', 'D', 'K', 'P', 'C', 'G', 'K', 'D',
              'A', 'G', 'G', 'R', 'V', 'C', 'T', 'N', 'N', 'Y', 'C', 'C', 'S',
              'K', 'W', 'G', 'S', 'C', 'G', 'I', 'G', 'P', 'G', 'Y', 'C', 'G',
              'A', 'G', 'C', 'Q', 'S', 'G', 'G', 'C']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 99.66)
        self.assertAlmostEqual(alignment.annotations["E-value"], 2.4e-21)
        self.assertAlmostEqual(alignment.annotations["Score"], 140.12)
        self.assertAlmostEqual(alignment.annotations["Identities"], 41)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.182)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[45:170],
            "CGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYCGAGCQSGGCD",
        )
        self.assertEqual(
            alignment[1],
            "CGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYCGAGCQSGGCD",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                             Cg~~~~~~~c~~~~CCs~~g~CG~~~~~c~~~c~~~~~~~~~~Cg~~~~~~~c~~~~CCS~~G~CG~~~~~C~~~Cq~~~c~~~~~Cg~~~~~~~c~~~~CCS~~G~CG~~~~~C~~gCq~~~c~ ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "  ~g~~~~~~~c~~~~CCs~~g~Cg~~~~~c~~~C~~~~c~~~~~cg~~~~~c~~~~CCs~~G~CG~~~~~C~~~C~~~~~~~~~~Cg~~~~~c~~~~Ccs~~G~CGt~~~~C~~~cq~~~c~                                         ",
        )
        self.assertEqual(alignment.target.id, "2wga")
        self.assertEqual(
            alignment.target.seq[2:123],
            "XGCXGXXMYCSTNNCCXXWESCGSGGYXCGEGCNLGACQXGXPCXXPGSVCTNLHCCARGGHCGMGSGYCGXGCXGGACXADIXCGXGXXXCPTDSCCGGWGXCGNGXEFCGXGCXVGGCA",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "2wga")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "; lectin (agglutinin); NMR {}",
        )
        self.assertEqual(
            alignment[0],
            "XGCXGXXMYCSTNNCCXXWESCGSGGYXCGEGCNLGACQXGXPCXX--PGSVCTNLHCCARGGHCGMGSGYCGXGCXGGACXADIXCGXG--XXXCPTDSCCGGWGXCGNGXEFCGXGCXVGGCA",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "  SGGGSSCCCCSTTCEECTTSCEECSTTTTSTTCCSSSCSSCCCSSSSSCCCSTTCEECTTSCEESSHHHHSSCCSSSSCSSCCCCTTSSSCCSTTCBCCSSSCCBCSHHHHSTTCCSSSCS                                         ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "  cCCCcccccCCCCceECCcceECCCCccccCccccCcccccceeccCCCcCCCCcccCCCceeCCCCcccCCCccccccccccccCcccccCCCCCccCCCCCccCccccccCCccccccc                                         ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  2,  48,  48,  90,  90, 123],
                             [ 45,  91,  93, 135, 137, 170]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['X', 'G', 'C', 'X', 'G', 'X', 'X', 'M', 'Y', 'C', 'S', 'T', 'N',
              'N', 'C', 'C', 'X', 'X', 'W', 'E', 'S', 'C', 'G', 'S', 'G', 'G',
              'Y', 'X', 'C', 'G', 'E', 'G', 'C', 'N', 'L', 'G', 'A', 'C', 'Q',
              'X', 'G', 'X', 'P', 'C', 'X', 'X', '-', '-', 'P', 'G', 'S', 'V',
              'C', 'T', 'N', 'L', 'H', 'C', 'C', 'A', 'R', 'G', 'G', 'H', 'C',
              'G', 'M', 'G', 'S', 'G', 'Y', 'C', 'G', 'X', 'G', 'C', 'X', 'G',
              'G', 'A', 'C', 'X', 'A', 'D', 'I', 'X', 'C', 'G', 'X', 'G', '-',
              '-', 'X', 'X', 'X', 'C', 'P', 'T', 'D', 'S', 'C', 'C', 'G', 'G',
              'W', 'G', 'X', 'C', 'G', 'N', 'G', 'X', 'E', 'F', 'C', 'G', 'X',
              'G', 'C', 'X', 'V', 'G', 'G', 'C', 'A'],
             ['C', 'G', 'S', 'Q', 'A', 'G', 'G', 'A', 'T', 'C', 'T', 'N', 'N',
              'Q', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'F', 'G', 'A',
              'E', 'Y', 'C', 'G', 'A', 'G', 'C', 'Q', 'G', 'G', 'P', 'C', 'R',
              'A', 'D', 'I', 'K', 'C', 'G', 'S', 'Q', 'A', 'G', 'G', 'K', 'L',
              'C', 'P', 'N', 'N', 'L', 'C', 'C', 'S', 'Q', 'W', 'G', 'F', 'C',
              'G', 'L', 'G', 'S', 'E', 'F', 'C', 'G', 'G', 'G', 'C', 'Q', 'S',
              'G', 'A', 'C', 'S', 'T', 'D', 'K', 'P', 'C', 'G', 'K', 'D', 'A',
              'G', 'G', 'R', 'V', 'C', 'T', 'N', 'N', 'Y', 'C', 'C', 'S', 'K',
              'W', 'G', 'S', 'C', 'G', 'I', 'G', 'P', 'G', 'Y', 'C', 'G', 'A',
              'G', 'C', 'Q', 'S', 'G', 'G', 'C', 'D']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 99.31)
        self.assertAlmostEqual(alignment.annotations["E-value"], 1.2e-16)
        self.assertAlmostEqual(alignment.annotations["Score"], 102.02)
        self.assertAlmostEqual(alignment.annotations["Identities"], 47)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.275)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 66.8)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[43:123],
            "KRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQS",
        )
        self.assertEqual(
            alignment[1],
            "KRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQS",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                           ~~Cg~~~~~~~c~~~~CCs~~g~CG~~~~~c~~~c~~~~~~~~~~Cg~~~~~~~c~~~~CCS~~G~CG~~~~~C~~~Cq~                                                ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            " ~~cg~~~~~~~C~~g~CCs~~G~Cg~~~~~c~~~c~~~~~~g~Cg~~~~~~~c~~~~CCs~~g~CG~~~~~C~~~Cqs   ",
        )
        self.assertEqual(alignment.target.id, "1uha_A")
        self.assertEqual(
            alignment.target.seq[1:79],
            "PECGERASGKRCPNGKCCSQWGYCGTTDNYCGQGCQSQCDYWRCGRDFGGRLCEEDMCCSKYGWCGYSDDHCEDGCQS",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "1uha_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Lectin-D2; chitin-binding domain, sugar binding protein; 1.50A {Phytolacca americana} SCOP: g.3.1.1 g.3.1.1 PDB: 1ulm_A* 1uln_A",
        )
        self.assertEqual(
            alignment[0],
            "PECGERASGKRCPNGKCCSQWGYCGTTDNYCGQGCQSQ--CDYWRCGRDFGGRLCEEDMCCSKYGWCGYSDDHCEDGCQS",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            " CCSBGGGTTBCCGGGCEECTTSCEESSHHHHSTTCCBCTTTTBCBGGGTTBCCSTTCEECTTSBEECSHHHHSTTCCB   ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            " CCCCCCCCCCcCCCCCccCCCccccCccccccCCccccccccccccccceecCCCCCccCCCccccCCcccccccccc   ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            " 467776667789999999999999999999999998642356788765567788999999999999999999999986   ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  1,  39,  39,  79],
                             [ 43,  81,  83, 123]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['P', 'E', 'C', 'G', 'E', 'R', 'A', 'S', 'G', 'K', 'R', 'C', 'P',
              'N', 'G', 'K', 'C', 'C', 'S', 'Q', 'W', 'G', 'Y', 'C', 'G', 'T',
              'T', 'D', 'N', 'Y', 'C', 'G', 'Q', 'G', 'C', 'Q', 'S', 'Q', '-',
              '-', 'C', 'D', 'Y', 'W', 'R', 'C', 'G', 'R', 'D', 'F', 'G', 'G',
              'R', 'L', 'C', 'E', 'E', 'D', 'M', 'C', 'C', 'S', 'K', 'Y', 'G',
              'W', 'C', 'G', 'Y', 'S', 'D', 'D', 'H', 'C', 'E', 'D', 'G', 'C',
              'Q', 'S'],
             ['K', 'R', 'C', 'G', 'S', 'Q', 'A', 'G', 'G', 'A', 'T', 'C', 'T',
              'N', 'N', 'Q', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'F',
              'G', 'A', 'E', 'Y', 'C', 'G', 'A', 'G', 'C', 'Q', 'G', 'G', 'P',
              'C', 'R', 'A', 'D', 'I', 'K', 'C', 'G', 'S', 'Q', 'A', 'G', 'G',
              'K', 'L', 'C', 'P', 'N', 'N', 'L', 'C', 'C', 'S', 'Q', 'W', 'G',
              'F', 'C', 'G', 'L', 'G', 'S', 'E', 'F', 'C', 'G', 'G', 'G', 'C',
              'Q', 'S']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 99.26)
        self.assertAlmostEqual(alignment.annotations["E-value"], 3.7e-16)
        self.assertAlmostEqual(alignment.annotations["Score"], 101.24)
        self.assertAlmostEqual(alignment.annotations["Identities"], 47)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.259)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 71.9)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[0:82],
            "ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGP",
        )
        self.assertEqual(
            alignment[1],
            "ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNGACW----TSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAG-CQGGP",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "~~cg~~~~~~~c~~~~CCs~~g~CG~~~~~c~~~c~~~~c~~~~~Cg~~~~~~~c~~~~CCs~~g~CG~~~~~c~~~c~~~~                                                                                         ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "~~cg~~~~~~~C~~~~CCS~~G~CG~~~~~C~~~Cq~~c~~~~~~~~~Cg~~~~~~~c~~~~CCS~~G~CG~~~~~C~~~~Cq~~~   ",
        )
        self.assertEqual(alignment.target.id, "1en2_A")
        self.assertEqual(
            alignment.target.seq[0:86],
            "ERCGSQGGGSTCPGLRCCSIWGWCGDSEPYCGRTCENKCWSGERSDHRCGAAVGNPPCGQDRCCSVHGWCGGGNDYCSGGNCQYRC",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "1en2_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "UDA, agglutinin isolectin I/agglutinin isolectin V/ AG isolectin VI; hevein domain, superantigen, saccharide binding binding protein; HET: NAG; 1.40A {Urtica dioica} SCOP: g.3.1.1 g.3.1.1 PDB: 1eis_A* 1enm_A* 1ehd_A 1ehh_A* 1iqb_A",
        )
        self.assertEqual(
            alignment[0],
            "ERCGSQGGGSTCPGLRCCSIWGWCGDSEPYCGRTCENK-CWSGERSDHRCGAAVGNPPCGQDRCCSVHGWCGGGNDYCSGGNCQYRC",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "CBCTTTTTSCCCGGGCEEETTSBEESSHHHHSTTEEESCGGGCCTTCBCSGGGTCCCCCTTCEEETTSBEESSHHHHSGGGEEECC   ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "CCCCCCCCCccCCCCcccCCCceecccccccCCCCcCCCcccccCCcccCCcccccccCCCCeECCCceECCCccccCCCCcccCC   ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            "689999988999999999999999999999999999764347889986556788999999999999999999998698753        ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 0, 38, 38, 40, 44, 80, 81, 86],
                             [ 0, 38, 39, 41, 41, 77, 77, 82]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['E', 'R', 'C', 'G', 'S', 'Q', 'G', 'G', 'G', 'S', 'T', 'C', 'P',
              'G', 'L', 'R', 'C', 'C', 'S', 'I', 'W', 'G', 'W', 'C', 'G', 'D',
              'S', 'E', 'P', 'Y', 'C', 'G', 'R', 'T', 'C', 'E', 'N', 'K', '-',
              'C', 'W', 'S', 'G', 'E', 'R', 'S', 'D', 'H', 'R', 'C', 'G', 'A',
              'A', 'V', 'G', 'N', 'P', 'P', 'C', 'G', 'Q', 'D', 'R', 'C', 'C',
              'S', 'V', 'H', 'G', 'W', 'C', 'G', 'G', 'G', 'N', 'D', 'Y', 'C',
              'S', 'G', 'G', 'N', 'C', 'Q', 'Y', 'R', 'C'],
             ['E', 'R', 'C', 'G', 'E', 'Q', 'G', 'S', 'N', 'M', 'E', 'C', 'P',
              'N', 'N', 'L', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'M',
              'G', 'G', 'D', 'Y', 'C', 'G', 'K', 'G', 'C', 'Q', 'N', 'G', 'A',
              'C', 'W', '-', '-', '-', '-', 'T', 'S', 'K', 'R', 'C', 'G', 'S',
              'Q', 'A', 'G', 'G', 'A', 'T', 'C', 'T', 'N', 'N', 'Q', 'C', 'C',
              'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'F', 'G', 'A', 'E', 'Y', 'C',
              'G', 'A', 'G', '-', 'C', 'Q', 'G', 'G', 'P']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 99.22)
        self.assertAlmostEqual(alignment.annotations["E-value"], 7.4e-16)
        self.assertAlmostEqual(alignment.annotations["Score"], 99.76)
        self.assertAlmostEqual(alignment.annotations["Identities"], 47)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.232)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 67.7)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[44:124],
            "RCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSG",
        )
        self.assertEqual(
            alignment[1],
            "RCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCR----ADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGG-CQSG",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                            ~Cg~~~~~~~c~~~~CCs~~g~CG~~~~~c~~~c~~~~~~~~~~Cg~~~~~~~c~~~~CCS~~G~CG~~~~~C~~~Cq~~                                               ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            " ~cg~~~~~~~C~~~~CCS~~G~CG~~~~~C~~~Cq~~c~~~~~~~~~Cg~~~~~~~c~~~~CCS~~G~CG~~~~~C~~~~Cq~~    ",
        )
        self.assertEqual(alignment.target.id, "1en2_A")
        self.assertEqual(
            alignment.target.seq[1:85],
            "RCGSQGGGSTCPGLRCCSIWGWCGDSEPYCGRTCENKCWSGERSDHRCGAAVGNPPCGQDRCCSVHGWCGGGNDYCSGGNCQYR",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "1en2_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "UDA, agglutinin isolectin I/agglutinin isolectin V/ AG isolectin VI; hevein domain, superantigen, saccharide binding binding protein; HET: NAG; 1.40A {Urtica dioica} SCOP: g.3.1.1 g.3.1.1 PDB: 1eis_A* 1enm_A* 1ehd_A 1ehh_A* 1iqb_A",
        )
        self.assertEqual(
            alignment[0],
            "RCGSQGGGSTCPGLRCCSIWGWCGDSEPYCGRTCENK-CWSGERSDHRCGAAVGNPPCGQDRCCSVHGWCGGGNDYCSGGNCQYR",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            " BCTTTTTSCCCGGGCEEETTSBEESSHHHHSTTEEESCGGGCCTTCBCSGGGTCCCCCTTCEEETTSBEESSHHHHSGGGEEEC    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            " CCCCCCCCccCCCCcccCCCceecccccccCCCCcCCCcccccCCcccCCcccccccCCCCeECCCceECCCccccCCCCcccC    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            " 5777666678999999999999999999999999875434788887654567889999999999999999999979974         ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  1,  38,  38,  40,  44,  80,  81,  85],
                             [ 44,  81,  82,  84,  84, 120, 120, 124]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['R', 'C', 'G', 'S', 'Q', 'G', 'G', 'G', 'S', 'T', 'C', 'P', 'G',
              'L', 'R', 'C', 'C', 'S', 'I', 'W', 'G', 'W', 'C', 'G', 'D', 'S',
              'E', 'P', 'Y', 'C', 'G', 'R', 'T', 'C', 'E', 'N', 'K', '-', 'C',
              'W', 'S', 'G', 'E', 'R', 'S', 'D', 'H', 'R', 'C', 'G', 'A', 'A',
              'V', 'G', 'N', 'P', 'P', 'C', 'G', 'Q', 'D', 'R', 'C', 'C', 'S',
              'V', 'H', 'G', 'W', 'C', 'G', 'G', 'G', 'N', 'D', 'Y', 'C', 'S',
              'G', 'G', 'N', 'C', 'Q', 'Y', 'R'],
             ['R', 'C', 'G', 'S', 'Q', 'A', 'G', 'G', 'A', 'T', 'C', 'T', 'N',
              'N', 'Q', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'F', 'G',
              'A', 'E', 'Y', 'C', 'G', 'A', 'G', 'C', 'Q', 'G', 'G', 'P', 'C',
              'R', '-', '-', '-', '-', 'A', 'D', 'I', 'K', 'C', 'G', 'S', 'Q',
              'A', 'G', 'G', 'K', 'L', 'C', 'P', 'N', 'N', 'L', 'C', 'C', 'S',
              'Q', 'W', 'G', 'F', 'C', 'G', 'L', 'G', 'S', 'E', 'F', 'C', 'G',
              'G', 'G', '-', 'C', 'Q', 'S', 'G']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 99.21)
        self.assertAlmostEqual(alignment.annotations["E-value"], 9.6e-16)
        self.assertAlmostEqual(alignment.annotations["Score"], 97.65)
        self.assertAlmostEqual(alignment.annotations["Identities"], 54)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.333)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 65.9)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[87:167],
            "KCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYCGAGCQSG",
        )
        self.assertEqual(
            alignment[1],
            "KCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYCGAGCQSG",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                       ~Cg~~~~~~~c~~~~CCS~~G~CG~~~~~C~~~Cq~~~c~~~~~Cg~~~~~~~c~~~~CCS~~G~CG~~~~~C~~gCq~~    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "  ~cg~~~~~~~C~~g~CCs~~G~Cg~~~~~c~~~c~~~~~~g~Cg~~~~~~~c~~~~CCs~~g~CG~~~~~C~~~Cqs~  ",
        )
        self.assertEqual(alignment.target.id, "1uha_A")
        self.assertEqual(
            alignment.target.seq[2:80],
            "ECGERASGKRCPNGKCCSQWGYCGTTDNYCGQGCQSQCDYWRCGRDFGGRLCEEDMCCSKYGWCGYSDDHCEDGCQSQ",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "1uha_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Lectin-D2; chitin-binding domain, sugar binding protein; 1.50A {Phytolacca americana} SCOP: g.3.1.1 g.3.1.1 PDB: 1ulm_A* 1uln_A",
        )
        self.assertEqual(
            alignment[0],
            "ECGERASGKRCPNGKCCSQWGYCGTTDNYCGQGCQSQ-C-DYWRCGRDFGGRLCEEDMCCSKYGWCGYSDDHCEDGCQSQ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "  CSBGGGTTBCCGGGCEECTTSCEESSHHHHSTTCCBCTTTTBCBGGGTTBCCSTTCEECTTSBEECSHHHHSTTCCBC  ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "  CCCCCCCCCcCCCCCccCCCccccCccccccCCccccccccccccccceecCCCCCccCCCccccCCccccccccccc  ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            "  466654457789999999999999999999999997423467887655678899999999999999999999999974  ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  2,  39,  39,  40,  40,  80],
                             [ 87, 124, 125, 126, 127, 167]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['E', 'C', 'G', 'E', 'R', 'A', 'S', 'G', 'K', 'R', 'C', 'P', 'N',
              'G', 'K', 'C', 'C', 'S', 'Q', 'W', 'G', 'Y', 'C', 'G', 'T', 'T',
              'D', 'N', 'Y', 'C', 'G', 'Q', 'G', 'C', 'Q', 'S', 'Q', '-', 'C',
              '-', 'D', 'Y', 'W', 'R', 'C', 'G', 'R', 'D', 'F', 'G', 'G', 'R',
              'L', 'C', 'E', 'E', 'D', 'M', 'C', 'C', 'S', 'K', 'Y', 'G', 'W',
              'C', 'G', 'Y', 'S', 'D', 'D', 'H', 'C', 'E', 'D', 'G', 'C', 'Q',
              'S', 'Q'],
             ['K', 'C', 'G', 'S', 'Q', 'A', 'G', 'G', 'K', 'L', 'C', 'P', 'N',
              'N', 'L', 'C', 'C', 'S', 'Q', 'W', 'G', 'F', 'C', 'G', 'L', 'G',
              'S', 'E', 'F', 'C', 'G', 'G', 'G', 'C', 'Q', 'S', 'G', 'A', 'C',
              'S', 'T', 'D', 'K', 'P', 'C', 'G', 'K', 'D', 'A', 'G', 'G', 'R',
              'V', 'C', 'T', 'N', 'N', 'Y', 'C', 'C', 'S', 'K', 'W', 'G', 'S',
              'C', 'G', 'I', 'G', 'P', 'G', 'Y', 'C', 'G', 'A', 'G', 'C', 'Q',
              'S', 'G']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 97.79)
        self.assertAlmostEqual(alignment.annotations["E-value"], 7.7e-09)
        self.assertAlmostEqual(alignment.annotations["Score"], 56.76)
        self.assertAlmostEqual(alignment.annotations["Identities"], 68)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.600)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 30.2)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[87:124], "KCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSG"
        )
        self.assertEqual(alignment[1], "KCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGG--GCQSG")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                       ~Cg~~~~~~~c~~~~CCS~~G~CG~~~~~C~~~Cq~~                                               ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            " ~CG~~~~~~~C~~~~CCS~~G~CG~t~~~C~~~~~Cq~~   ",
        )
        self.assertEqual(alignment.target.id, "1wkx_A")
        self.assertEqual(
            alignment.target.seq[1:40], "QCGRQAGGKLCPDNLCCSQWGWCGSTDEYCSPDHNCQSN"
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "1wkx_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Hevein isoform 2; allergen, lectin, agglutinin-toxin motif; 1.70A {Hevea brasiliensis} PDB: 1hev_A 1q9b_A* 4wp4_A",
        )
        self.assertEqual(alignment[0], "QCGRQAGGKLCPDNLCCSQWGWCGSTDEYCSPDHNCQSN")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            " CSBGGGTTBCCSTTCEECTTSCEESSHHHHCGGGTCCBS   ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            " CCCCcCCCcccCCCCeEeecCcccCCcccccCCCCccCC   ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            " 3555444567889999999999999999998679874     ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  1,  33,  35,  40],
                             [ 87, 119, 119, 124]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['Q', 'C', 'G', 'R', 'Q', 'A', 'G', 'G', 'K', 'L', 'C', 'P', 'D',
              'N', 'L', 'C', 'C', 'S', 'Q', 'W', 'G', 'W', 'C', 'G', 'S', 'T',
              'D', 'E', 'Y', 'C', 'S', 'P', 'D', 'H', 'N', 'C', 'Q', 'S', 'N'],
             ['K', 'C', 'G', 'S', 'Q', 'A', 'G', 'G', 'K', 'L', 'C', 'P', 'N',
              'N', 'L', 'C', 'C', 'S', 'Q', 'W', 'G', 'F', 'C', 'G', 'L', 'G',
              'S', 'E', 'F', 'C', 'G', 'G', '-', '-', 'G', 'C', 'Q', 'S', 'G']],
            dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 97.67)
        self.assertAlmostEqual(alignment.annotations["E-value"], 1.9e-08)
        self.assertAlmostEqual(alignment.annotations["Score"], 54.61)
        self.assertAlmostEqual(alignment.annotations["Identities"], 46)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.182)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 32.7)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[0:38], "ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNG"
        )
        self.assertEqual(alignment[1], "ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGK-GCQNG")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "~~cg~~~~~~~c~~~~CCs~~g~CG~~~~~c~~~c~~~                                                                                                                                     ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "~~CG~~~~~~C~~~~CCS~~G~CG~t~~~C~~~~Cq~~   ",
        )
        self.assertEqual(alignment.target.id, "1p9g_A")
        self.assertEqual(
            alignment.target.seq[0:38], "ETCASRCPRPCNAGLCCSIYGYCGSGAAYCGAGNCRCQ"
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "1p9g_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "EAFP 2; antifungal peptide, atomic resolution, antifungal protein; HET: PCA; 0.84A {Eucommia ulmoides} SCOP: g.3.1.1 PDB: 1p9z_A*",
        )
        self.assertEqual(alignment[0], "ETCAS-RCPRPCNAGLCCSIYGYCGSGAAYCGAGNCRCQ")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "CCGGGGTTCCSCTTCEEETTSCEECSHHHHSTTTEEEC   ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "CCcCCcCCcccCCCCeECccceeCCCccccCCCccccC   ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            "6899567789999999999999999999999849864    ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 0,  5,  5, 32, 33, 38],
                             [ 0,  5,  6, 33, 33, 38]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['E', 'T', 'C', 'A', 'S', '-', 'R', 'C', 'P', 'R', 'P', 'C', 'N',
              'A', 'G', 'L', 'C', 'C', 'S', 'I', 'Y', 'G', 'Y', 'C', 'G', 'S',
              'G', 'A', 'A', 'Y', 'C', 'G', 'A', 'G', 'N', 'C', 'R', 'C', 'Q'],
             ['E', 'R', 'C', 'G', 'E', 'Q', 'G', 'S', 'N', 'M', 'E', 'C', 'P',
              'N', 'N', 'L', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'M',
              'G', 'G', 'D', 'Y', 'C', 'G', 'K', '-', 'G', 'C', 'Q', 'N', 'G']],
            dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 97.66)
        self.assertAlmostEqual(alignment.annotations["E-value"], 1.9e-08)
        self.assertAlmostEqual(alignment.annotations["Score"], 54.54)
        self.assertAlmostEqual(alignment.annotations["Identities"], 42)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.206)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 29.2)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[87:124], "KCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSG"
        )
        self.assertEqual(alignment[1], "KCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGG-CQSG")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                       ~Cg~~~~~~~c~~~~CCS~~G~CG~~~~~C~~~Cq~~                                               ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            " ~CG~~~~~~C~~~~CCS~~G~CG~t~~~C~~~~Cq~~   ",
        )
        self.assertEqual(alignment.target.id, "1p9g_A")
        self.assertEqual(
            alignment.target.seq[1:38], "TCASRCPRPCNAGLCCSIYGYCGSGAAYCGAGNCRCQ"
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "1p9g_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "EAFP 2; antifungal peptide, atomic resolution, antifungal protein; HET: PCA; 0.84A {Eucommia ulmoides} SCOP: g.3.1.1 PDB: 1p9z_A*",
        )
        self.assertEqual(alignment[0], "TCAS-RCPRPCNAGLCCSIYGYCGSGAAYCGAGNCRCQ")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            " CGGGGTTCCSCTTCEEETTSCEECSHHHHSTTTEEEC   ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            " CcCCcCCcccCCCCeECccceeCCCccccCCCccccC   ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            " 355333567889999999999999999999849874    ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  1,   5,   5,  33,  34,  38],
                             [ 87,  91,  92, 120, 120, 124]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['T', 'C', 'A', 'S', '-', 'R', 'C', 'P', 'R', 'P', 'C', 'N', 'A',
              'G', 'L', 'C', 'C', 'S', 'I', 'Y', 'G', 'Y', 'C', 'G', 'S', 'G',
              'A', 'A', 'Y', 'C', 'G', 'A', 'G', 'N', 'C', 'R', 'C', 'Q'],
             ['K', 'C', 'G', 'S', 'Q', 'A', 'G', 'G', 'K', 'L', 'C', 'P', 'N',
              'N', 'L', 'C', 'C', 'S', 'Q', 'W', 'G', 'F', 'C', 'G', 'L', 'G',
              'S', 'E', 'F', 'C', 'G', 'G', 'G', '-', 'C', 'Q', 'S', 'G']],
            dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 97.65)
        self.assertAlmostEqual(alignment.annotations["E-value"], 2e-08)
        self.assertAlmostEqual(alignment.annotations["Score"], 54.99)
        self.assertAlmostEqual(alignment.annotations["Identities"], 50)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.302)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 34.5)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[0:38], "ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNG"
        )
        self.assertEqual(alignment[1], "ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGK--GCQNG")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "~~cg~~~~~~~c~~~~CCs~~g~CG~~~~~c~~~c~~~                                                                                                                                     ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "~~CG~~~~~~~C~~~~CCS~~G~CG~t~~~C~~~~~Cq~~   ",
        )
        self.assertEqual(alignment.target.id, "1wkx_A")
        self.assertEqual(
            alignment.target.seq[0:40], "EQCGRQAGGKLCPDNLCCSQWGWCGSTDEYCSPDHNCQSN"
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "1wkx_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Hevein isoform 2; allergen, lectin, agglutinin-toxin motif; 1.70A {Hevea brasiliensis} PDB: 1hev_A 1q9b_A* 4wp4_A",
        )
        self.assertEqual(alignment[0], "EQCGRQAGGKLCPDNLCCSQWGWCGSTDEYCSPDHNCQSN")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "CCSBGGGTTBCCSTTCEECTTSCEESSHHHHCGGGTCCBS   ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "CCCCCcCCCcccCCCCeEeecCcccCCcccccCCCCccCC   ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            "57999988889999999999999999999999679864     ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 0, 33, 35, 40],
                             [ 0, 33, 33, 38]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['E', 'Q', 'C', 'G', 'R', 'Q', 'A', 'G', 'G', 'K', 'L', 'C', 'P',
              'D', 'N', 'L', 'C', 'C', 'S', 'Q', 'W', 'G', 'W', 'C', 'G', 'S',
              'T', 'D', 'E', 'Y', 'C', 'S', 'P', 'D', 'H', 'N', 'C', 'Q', 'S',
              'N'],
             ['E', 'R', 'C', 'G', 'E', 'Q', 'G', 'S', 'N', 'M', 'E', 'C', 'P',
              'N', 'N', 'L', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'M',
              'G', 'G', 'D', 'Y', 'C', 'G', 'K', '-', '-', 'G', 'C', 'Q', 'N',
              'G']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 97.62)
        self.assertAlmostEqual(alignment.annotations["E-value"], 2.4e-08)
        self.assertAlmostEqual(alignment.annotations["Score"], 55.16)
        self.assertAlmostEqual(alignment.annotations["Identities"], 51)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.325)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 33.8)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[1:38], "RCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNG"
        )
        self.assertEqual(alignment[1], "RCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNG")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            " ~cg~~~~~~~c~~~~CCs~~g~CG~~~~~c~~~c~~~                                                                                                                                     ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "   ~CG~~~~~~~C~~~~CCs~~G~CG~t~~~C~~gCq~~     ",
        )
        self.assertEqual(alignment.target.id, "4mpi_A")
        self.assertEqual(
            alignment.target.seq[3:40], "QCGRQAGGALCPGGLCCSQYGWCANTPEYCGSGCQSQ"
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "4mpi_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Class I chitinase; hevein-like domain, chitin oligomers, sugar binding protein; HET: MES; 1.60A {Hevea brasiliensis subsp}",
        )
        self.assertEqual(alignment[0], "QCGRQAGGALCPGGLCCSQYGWCANTPEYCGSGCQSQ")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "   BCBGGGTTBCCGGGCEECTTSBEECSHHHHSTTCCBC     ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "   ccCCcCCCcccCCCCcCcccceecCCccccccccccc     ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            "   6888887789999999999999999999999999864     ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 3, 40],
                             [ 1, 38]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['Q', 'C', 'G', 'R', 'Q', 'A', 'G', 'G', 'A', 'L', 'C', 'P', 'G',
              'G', 'L', 'C', 'C', 'S', 'Q', 'Y', 'G', 'W', 'C', 'A', 'N', 'T',
              'P', 'E', 'Y', 'C', 'G', 'S', 'G', 'C', 'Q', 'S', 'Q'],
             ['R', 'C', 'G', 'E', 'Q', 'G', 'S', 'N', 'M', 'E', 'C', 'P', 'N',
              'N', 'L', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'M', 'G',
              'G', 'D', 'Y', 'C', 'G', 'K', 'G', 'C', 'Q', 'N', 'G']],
            dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 97.59)
        self.assertAlmostEqual(alignment.annotations["E-value"], 3e-08)
        self.assertAlmostEqual(alignment.annotations["Score"], 54.44)
        self.assertAlmostEqual(alignment.annotations["Identities"], 59)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.493)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 30.2)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[87:124], "KCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSG"
        )
        self.assertEqual(alignment[1], "KCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGG-CQSG")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                       ~Cg~~~~~~~c~~~~CCS~~G~CG~~~~~C~~~Cq~~                                               ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "  ~CG~~~~~~~C~~~~CCS~~G~CG~t~~~C~~~~Cq~~    ",
        )
        self.assertEqual(alignment.target.id, "2lb7_A")
        self.assertEqual(
            alignment.target.seq[2:40], "RCGDQARGAKCPNCLCCGKYGFCGSGDAYCGAGSCQSQ"
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "2lb7_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "WAMP-1A, antimicrobial peptide 1A; antimicrobial protein; NMR {Triticum kiharae}",
        )
        self.assertEqual(alignment[0], "RCGDQARGAKCPNCLCCGKYGFCGSGDAYCGAGSCQSQ")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "  ECBGGGTTBCCCTTCEEETTTEEECSHHHHSTTSEEEC    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "  CCcCCCCCcccCCCCcCCcceeecCCccccCCCCccCC    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            "  4555444567889999999999999999999868864     ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  2,  35,  36,  40],
                             [ 87, 120, 120, 124]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['R', 'C', 'G', 'D', 'Q', 'A', 'R', 'G', 'A', 'K', 'C', 'P', 'N',
              'C', 'L', 'C', 'C', 'G', 'K', 'Y', 'G', 'F', 'C', 'G', 'S', 'G',
              'D', 'A', 'Y', 'C', 'G', 'A', 'G', 'S', 'C', 'Q', 'S', 'Q'],
             ['K', 'C', 'G', 'S', 'Q', 'A', 'G', 'G', 'K', 'L', 'C', 'P', 'N',
              'N', 'L', 'C', 'C', 'S', 'Q', 'W', 'G', 'F', 'C', 'G', 'L', 'G',
              'S', 'E', 'F', 'C', 'G', 'G', 'G', '-', 'C', 'Q', 'S', 'G']],
            dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 97.57)
        self.assertAlmostEqual(alignment.annotations["E-value"], 3.3e-08)
        self.assertAlmostEqual(alignment.annotations["Score"], 54.59)
        self.assertAlmostEqual(alignment.annotations["Identities"], 61)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.420)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 31.3)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[86:124], "IKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSG"
        )
        self.assertEqual(alignment[1], "IKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSG")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                      ~~Cg~~~~~~~c~~~~CCS~~G~CG~~~~~C~~~Cq~~                                               ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "  ~~CG~~~~~~~C~~~~CCs~~G~CG~t~~~C~~gCq~~     ",
        )
        self.assertEqual(alignment.target.id, "4mpi_A")
        self.assertEqual(
            alignment.target.seq[2:40], "EQCGRQAGGALCPGGLCCSQYGWCANTPEYCGSGCQSQ"
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "4mpi_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Class I chitinase; hevein-like domain, chitin oligomers, sugar binding protein; HET: MES; 1.60A {Hevea brasiliensis subsp}",
        )
        self.assertEqual(alignment[0], "EQCGRQAGGALCPGGLCCSQYGWCANTPEYCGSGCQSQ")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "  CBCBGGGTTBCCGGGCEECTTSBEECSHHHHSTTCCBC     ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "  cccCCcCCCcccCCCCcCcccceecCCccccccccccc     ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            "  45665544567889999999999999999999999873     ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  2,  40],
                             [ 86, 124]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['E', 'Q', 'C', 'G', 'R', 'Q', 'A', 'G', 'G', 'A', 'L', 'C', 'P',
              'G', 'G', 'L', 'C', 'C', 'S', 'Q', 'Y', 'G', 'W', 'C', 'A', 'N',
              'T', 'P', 'E', 'Y', 'C', 'G', 'S', 'G', 'C', 'Q', 'S', 'Q'],
             ['I', 'K', 'C', 'G', 'S', 'Q', 'A', 'G', 'G', 'K', 'L', 'C', 'P',
              'N', 'N', 'L', 'C', 'C', 'S', 'Q', 'W', 'G', 'F', 'C', 'G', 'L',
              'G', 'S', 'E', 'F', 'C', 'G', 'G', 'G', 'C', 'Q', 'S', 'G']],
            dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 97.46)
        self.assertAlmostEqual(alignment.annotations["E-value"], 6.6e-08)
        self.assertAlmostEqual(alignment.annotations["Score"], 52.98)
        self.assertAlmostEqual(alignment.annotations["Identities"], 55)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.488)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 34.0)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[0:38], "ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNG"
        )
        self.assertEqual(alignment[1], "ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKG-CQNG")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "~~cg~~~~~~~c~~~~CCs~~g~CG~~~~~c~~~c~~~                                                                                                                                     ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            " ~~CG~~~~~~~C~~~~CCS~~G~CG~t~~~C~~~~Cq~~    ",
        )
        self.assertEqual(alignment.target.id, "2lb7_A")
        self.assertEqual(
            alignment.target.seq[1:40], "QRCGDQARGAKCPNCLCCGKYGFCGSGDAYCGAGSCQSQ"
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "2lb7_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "WAMP-1A, antimicrobial peptide 1A; antimicrobial protein; NMR {Triticum kiharae}",
        )
        self.assertEqual(alignment[0], "QRCGDQARGAKCPNCLCCGKYGFCGSGDAYCGAGSCQSQ")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            " EECBGGGTTBCCCTTCEEETTTEEECSHHHHSTTSEEEC    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            " CCCcCCCCCcccCCCCcCCcceeecCCccccCCCCccCC    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            " 47899987789999999999999999999999868864     ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 1, 35, 36, 40],
                             [ 0, 34, 34, 38]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['Q', 'R', 'C', 'G', 'D', 'Q', 'A', 'R', 'G', 'A', 'K', 'C', 'P',
              'N', 'C', 'L', 'C', 'C', 'G', 'K', 'Y', 'G', 'F', 'C', 'G', 'S',
              'G', 'D', 'A', 'Y', 'C', 'G', 'A', 'G', 'S', 'C', 'Q', 'S', 'Q'],
             ['E', 'R', 'C', 'G', 'E', 'Q', 'G', 'S', 'N', 'M', 'E', 'C', 'P',
              'N', 'N', 'L', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'M',
              'G', 'G', 'D', 'Y', 'C', 'G', 'K', 'G', '-', 'C', 'Q', 'N', 'G']],
            dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 97.25)
        self.assertAlmostEqual(alignment.annotations["E-value"], 2.2e-07)
        self.assertAlmostEqual(alignment.annotations["Score"], 48.31)
        self.assertAlmostEqual(alignment.annotations["Identities"], 50)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.296)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 27.1)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[41:75], "TSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCG"
        )
        self.assertEqual(alignment[1], "TSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCG")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                         ~~~~Cg~~~~~~~c~~~~CCs~~g~CG~~~~~c~                                                                                                ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "  ~~~~CG~~~g~C~~g~CCS~~G~CG~~~~~C~ ",
        )
        self.assertEqual(alignment.target.id, "2kus_A")
        self.assertEqual(alignment.target.seq[2:34], "PNGQCGPGWGGCRGGLCCSQYGYCGSGPKYCA")
        self.assertEqual(alignment.target.annotations["hmm_name"], "2kus_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "SM-AMP-1.1A; plant antimicrobial peptide, chitin-binding peptide, antimic protein; NMR {Stellaria media}",
        )
        self.assertEqual(alignment[0], "PNGQCGPGWG--GCRGGLCCSQYGYCGSGPKYCA")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "  TTCBCBTTTBCCCTTCEECTTSBEECSHHHHC ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "  CCcccCCCCCcCCCCcEECCCceecCChhhhC ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            "  45678875436899999999999999999986 ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 2, 12, 12, 34],
                             [41, 51, 53, 75]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['P', 'N', 'G', 'Q', 'C', 'G', 'P', 'G', 'W', 'G', '-', '-', 'G',
              'C', 'R', 'G', 'G', 'L', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y', 'C',
              'G', 'S', 'G', 'P', 'K', 'Y', 'C', 'A'],
             ['T', 'S', 'K', 'R', 'C', 'G', 'S', 'Q', 'A', 'G', 'G', 'A', 'T',
              'C', 'T', 'N', 'N', 'Q', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y', 'C',
              'G', 'F', 'G', 'A', 'E', 'Y', 'C', 'G']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 97.13)
        self.assertAlmostEqual(alignment.annotations["E-value"], 4.2e-07)
        self.assertAlmostEqual(alignment.annotations["Score"], 47.19)
        self.assertAlmostEqual(alignment.annotations["Identities"], 44)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.207)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 26.5)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[84:118], "ADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCG"
        )
        self.assertEqual(alignment[1], "ADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCG")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                    ~~~~Cg~~~~~~~c~~~~CCS~~G~CG~~~~~C~                                                     ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "  ~~~~CG~~~g~C~~g~CCS~~G~CG~~~~~C~ ",
        )
        self.assertEqual(alignment.target.id, "2kus_A")
        self.assertEqual(alignment.target.seq[2:34], "PNGQCGPGWGGCRGGLCCSQYGYCGSGPKYCA")
        self.assertEqual(alignment.target.annotations["hmm_name"], "2kus_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "SM-AMP-1.1A; plant antimicrobial peptide, chitin-binding peptide, antimic protein; NMR {Stellaria media}",
        )
        self.assertEqual(alignment[0], "PNGQCGPGWG--GCRGGLCCSQYGYCGSGPKYCA")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "  TTCBCBTTTBCCCTTCEECTTSBEECSHHHHC ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "  CCcccCCCCCcCCCCcEECCCceecCChhhhC ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            "  45677775436889999999999999999986 ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  2,  12,  12,  34],
                             [ 84,  94,  96, 118]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['P', 'N', 'G', 'Q', 'C', 'G', 'P', 'G', 'W', 'G', '-', '-', 'G',
              'C', 'R', 'G', 'G', 'L', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y', 'C',
              'G', 'S', 'G', 'P', 'K', 'Y', 'C', 'A'],
             ['A', 'D', 'I', 'K', 'C', 'G', 'S', 'Q', 'A', 'G', 'G', 'K', 'L',
              'C', 'P', 'N', 'N', 'L', 'C', 'C', 'S', 'Q', 'W', 'G', 'F', 'C',
              'G', 'L', 'G', 'S', 'E', 'F', 'C', 'G']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 97.04)
        self.assertAlmostEqual(alignment.annotations["E-value"], 6.8e-07)
        self.assertAlmostEqual(alignment.annotations["Score"], 44.77)
        self.assertAlmostEqual(alignment.annotations["Identities"], 48)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.444)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 22.0)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(alignment.query.seq[94:119], "GKLCPNNLCCSQWGFCGLGSEFCGG")
        self.assertEqual(alignment[1], "GKLCPNNLCCSQWGFCGLGSEFCGG")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                              ~~~c~~~~CCS~~G~CG~~~~~C~~                                                    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "     ~~~C~~~~CCS~~G~CG~t~~~C~~",
        )
        self.assertEqual(alignment.target.id, "1mmc_A")
        self.assertEqual(alignment.target.seq[5:30], "RGRCPSGMCCSQFGYCGKGPKYCGR")
        self.assertEqual(alignment.target.annotations["hmm_name"], "1mmc_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "AC-AMP2, antimicrobial peptide 2; antifungal antimicrobial, chitin-binding; NMR {Amaranthus caudatus} SCOP: g.3.1.2 PDB: 1zuv_A 1zwu_A* 1znt_A*",
        )
        self.assertEqual(alignment[0], "RGRCPSGMCCSQFGYCGKGPKYCGR")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "     SSCCSTTCEECTTSCEESSHHHHCC",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "     cCCCCCCCcccccceeCCchHhhCc",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            "     3468899999999999999999963",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  5,  30],
                             [ 94, 119]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['R', 'G', 'R', 'C', 'P', 'S', 'G', 'M', 'C', 'C', 'S', 'Q', 'F',
              'G', 'Y', 'C', 'G', 'K', 'G', 'P', 'K', 'Y', 'C', 'G', 'R'],
             ['G', 'K', 'L', 'C', 'P', 'N', 'N', 'L', 'C', 'C', 'S', 'Q', 'W',
              'G', 'F', 'C', 'G', 'L', 'G', 'S', 'E', 'F', 'C', 'G', 'G']],
            dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 96.99)
        self.assertAlmostEqual(alignment.annotations["E-value"], 8.4e-07)
        self.assertAlmostEqual(alignment.annotations["Score"], 44.41)
        self.assertAlmostEqual(alignment.annotations["Identities"], 52)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.425)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 23.6)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(alignment.query.seq[6:33], "GSNMECPNNLCCSQYGYCGMGGDYCGK")
        self.assertEqual(alignment[1], "GSNMECPNNLCCSQYGYCGMGGDYCGK")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "      ~~~~~c~~~~CCs~~g~CG~~~~~c~~                                                                                                                                          ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "   ~~~~~C~~~~CCS~~G~CG~t~~~C~~",
        )
        self.assertEqual(alignment.target.id, "1mmc_A")
        self.assertEqual(alignment.target.seq[3:30], "CVRGRCPSGMCCSQFGYCGKGPKYCGR")
        self.assertEqual(alignment.target.annotations["hmm_name"], "1mmc_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "AC-AMP2, antimicrobial peptide 2; antifungal antimicrobial, chitin-binding; NMR {Amaranthus caudatus} SCOP: g.3.1.2 PDB: 1zuv_A 1zwu_A* 1znt_A*",
        )
        self.assertEqual(alignment[0], "CVRGRCPSGMCCSQFGYCGKGPKYCGR")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "   CSSSCCSTTCEECTTSCEESSHHHHCC",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "   CccCCCCCCCcccccceeCCchHhhCc",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            "   345689999999999999999999973",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 3, 30],
                             [ 6, 33]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['C', 'V', 'R', 'G', 'R', 'C', 'P', 'S', 'G', 'M', 'C', 'C', 'S',
              'Q', 'F', 'G', 'Y', 'C', 'G', 'K', 'G', 'P', 'K', 'Y', 'C', 'G',
              'R'],
             ['G', 'S', 'N', 'M', 'E', 'C', 'P', 'N', 'N', 'L', 'C', 'C', 'S',
              'Q', 'Y', 'G', 'Y', 'C', 'G', 'M', 'G', 'G', 'D', 'Y', 'C', 'G',
              'K']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 96.90)
        self.assertAlmostEqual(alignment.annotations["E-value"], 1.3e-06)
        self.assertAlmostEqual(alignment.annotations["Score"], 43.64)
        self.assertAlmostEqual(alignment.annotations["Identities"], 48)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.477)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 20.7)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(alignment.query.seq[95:118], "KLCPNNLCCSQWGFCGLGSEFCG")
        self.assertEqual(alignment[1], "KLCPNNLCCSQWGFCGLGSEFCG")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                               ~~c~~~~CCS~~G~CG~~~~~C~                                                     ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "       ~~C~~~~CCs~~G~CGtt~~~C~",
        )
        self.assertEqual(alignment.target.id, "2n1s_A")
        self.assertEqual(alignment.target.seq[7:30], "GRCSGGLCCSKYGYCGSGPAYCG")
        self.assertEqual(alignment.target.annotations["hmm_name"], "2n1s_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "AMP-2; antimicrobial peptide, ICK, cystine knot inhibitor, cystine antimicrobial protein; NMR {Stellaria media}",
        )
        self.assertEqual(alignment[0], "GRCSGGLCCSKYGYCGSGPAYCG")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "       TBCSTTCEECTTSBEECSHHHHC",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "       CCCCCCCccccccccCcchhhcC",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            "       36888999999999999999985",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  7,  30],
                             [ 95, 118]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['G', 'R', 'C', 'S', 'G', 'G', 'L', 'C', 'C', 'S', 'K', 'Y', 'G',
              'Y', 'C', 'G', 'S', 'G', 'P', 'A', 'Y', 'C', 'G'],
             ['K', 'L', 'C', 'P', 'N', 'N', 'L', 'C', 'C', 'S', 'Q', 'W', 'G',
              'F', 'C', 'G', 'L', 'G', 'S', 'E', 'F', 'C', 'G']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 96.77)
        self.assertAlmostEqual(alignment.annotations["E-value"], 2.3e-06)
        self.assertAlmostEqual(alignment.annotations["Score"], 42.68)
        self.assertAlmostEqual(alignment.annotations["Identities"], 56)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.468)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 23.6)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(alignment.query.seq[1:32], "RCGEQGSNMECPNNLCCSQYGYCGMGGDYCG")
        self.assertEqual(alignment[1], "RCGEQGSNMECPNNLCCSQYGYCGMGGDYCG")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            " ~cg~~~~~~~c~~~~CCs~~g~CG~~~~~c~                                                                                                                                           ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "   ~cg~~~C~~~~CCs~~G~CGtt~~~C~",
        )
        self.assertEqual(alignment.target.id, "2n1s_A")
        self.assertEqual(alignment.target.seq[3:30], "QCYRGRCSGGLCCSKYGYCGSGPAYCG")
        self.assertEqual(alignment.target.annotations["hmm_name"], "2n1s_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "AMP-2; antimicrobial peptide, ICK, cystine knot inhibitor, cystine antimicrobial protein; NMR {Stellaria media}",
        )
        self.assertEqual(alignment[0], "QCY----RGRCSGGLCCSKYGYCGSGPAYCG")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "   BCBTTBCSTTCEECTTSBEECSHHHHC",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "   hcCCCCCCCCCccccccccCcchhhcC",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            "   566258999999999999999999985",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 3,  6,  6, 30],
                             [ 1,  4,  8, 32]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['Q', 'C', 'Y', '-', '-', '-', '-', 'R', 'G', 'R', 'C', 'S', 'G',
              'G', 'L', 'C', 'C', 'S', 'K', 'Y', 'G', 'Y', 'C', 'G', 'S', 'G',
              'P', 'A', 'Y', 'C', 'G'],
             ['R', 'C', 'G', 'E', 'Q', 'G', 'S', 'N', 'M', 'E', 'C', 'P', 'N',
              'N', 'L', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'M', 'G',
              'G', 'D', 'Y', 'C', 'G']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 96.11)
        self.assertAlmostEqual(alignment.annotations["E-value"], 2.5e-05)
        self.assertAlmostEqual(alignment.annotations["Score"], 59.25)
        self.assertAlmostEqual(alignment.annotations["Identities"], 66)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.556)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 32.2)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[86:124], "IKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSG"
        )
        self.assertEqual(alignment[1], "IKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSG")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                      ~~Cg~~~~~~~c~~~~CCS~~G~CG~~~~~C~~~Cq~~                                               ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            " ~~cG~~~~~~~c~~~~ccs~~g~cg~~~~~C~~~cq~~                                                                                                                                                                                                                                                                              ",
        )
        self.assertEqual(alignment.target.id, "2dkv_A")
        self.assertEqual(
            alignment.target.seq[1:39], "EQCGAQAGGARCPNCLCCSRWGWCGTTSDFCGDGCQSQ"
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "2dkv_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Chitinase; whole structure, oryza sativa L. japonica, hydrolase; HET: MES; 2.00A {Oryza sativa japonica group} PDB: 3iwr_A*",
        )
        self.assertEqual(alignment[0], "EQCGAQAGGARCPNCLCCSRWGWCGTTSDFCGDGCQSQ")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            " CBCSTTTTTCCCGGGCEECTTSBEESSHHHHSTTCCBC                                                                                                                                                                                                                                                                              ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            " CCCCCCCCCCcCCCCCeeCcCCcccCCccccCccccCC                                                                                                                                                                                                                                                                              ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            " 35666555578999999999999999999999999975                                                                                                                                                                                                                                                                              ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  1,  39],
                             [ 86, 124]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['E', 'Q', 'C', 'G', 'A', 'Q', 'A', 'G', 'G', 'A', 'R', 'C', 'P',
              'N', 'C', 'L', 'C', 'C', 'S', 'R', 'W', 'G', 'W', 'C', 'G', 'T',
              'T', 'S', 'D', 'F', 'C', 'G', 'D', 'G', 'C', 'Q', 'S', 'Q'],
             ['I', 'K', 'C', 'G', 'S', 'Q', 'A', 'G', 'G', 'K', 'L', 'C', 'P',
              'N', 'N', 'L', 'C', 'C', 'S', 'Q', 'W', 'G', 'F', 'C', 'G', 'L',
              'G', 'S', 'E', 'F', 'C', 'G', 'G', 'G', 'C', 'Q', 'S', 'G']],
            dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 95.96)
        self.assertAlmostEqual(alignment.annotations["E-value"], 3.8e-05)
        self.assertAlmostEqual(alignment.annotations["Score"], 58.17)
        self.assertAlmostEqual(alignment.annotations["Identities"], 53)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.380)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 33.1)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[43:81], "KRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGG"
        )
        self.assertEqual(alignment[1], "KRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGG")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                           ~~Cg~~~~~~~c~~~~CCs~~g~CG~~~~~c~~~c~~~                                                                                          ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            " ~~cG~~~~~~~c~~~~ccs~~g~cg~~~~~C~~~cq~~                                                                                                                                                                                                                                                                              ",
        )
        self.assertEqual(alignment.target.id, "2dkv_A")
        self.assertEqual(
            alignment.target.seq[1:39], "EQCGAQAGGARCPNCLCCSRWGWCGTTSDFCGDGCQSQ"
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "2dkv_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Chitinase; whole structure, oryza sativa L. japonica, hydrolase; HET: MES; 2.00A {Oryza sativa japonica group} PDB: 3iwr_A*",
        )
        self.assertEqual(alignment[0], "EQCGAQAGGARCPNCLCCSRWGWCGTTSDFCGDGCQSQ")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            " CBCSTTTTTCCCGGGCEECTTSBEESSHHHHSTTCCBC                                                                                                                                                                                                                                                                              ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            " CCCCCCCCCCcCCCCCeeCcCCcccCCccccCccccCC                                                                                                                                                                                                                                                                              ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            " 46777666678999999999999999999999999865                                                                                                                                                                                                                                                                              ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 1, 39],
                             [43, 81]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['E', 'Q', 'C', 'G', 'A', 'Q', 'A', 'G', 'G', 'A', 'R', 'C', 'P',
              'N', 'C', 'L', 'C', 'C', 'S', 'R', 'W', 'G', 'W', 'C', 'G', 'T',
              'T', 'S', 'D', 'F', 'C', 'G', 'D', 'G', 'C', 'Q', 'S', 'Q'],
             ['K', 'R', 'C', 'G', 'S', 'Q', 'A', 'G', 'G', 'A', 'T', 'C', 'T',
              'N', 'N', 'Q', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'F',
              'G', 'A', 'E', 'Y', 'C', 'G', 'A', 'G', 'C', 'Q', 'G', 'G']],
            dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 85.64)
        self.assertAlmostEqual(alignment.annotations["E-value"], 0.043)
        self.assertAlmostEqual(alignment.annotations["Score"], 38.99)
        self.assertAlmostEqual(alignment.annotations["Identities"], 41)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.195)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(alignment.query.seq[143:165], "NYCCSKWGSCGIGPGYCGAGCQ")
        self.assertEqual(alignment[1], "NYCCSKWGSCGIGPGYCG-AGCQ")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                                                                               ~~CCS~~G~CG~~~~~C~~gCq      ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                                               ~~CCs~~~~Cg~~~~~C~~~~c~                                                                                                                                                                                          ",
        )
        self.assertEqual(alignment.target.id, "4zxm_A")
        self.assertEqual(alignment.target.seq[47:70], "DHCCSEWGWCGRETSHCTCSSCV")
        self.assertEqual(alignment.target.annotations["hmm_name"], "4zxm_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "PGRP domain of peptidoglycan recognition protein; amidase, hydrolase; 2.80A {Branchiostoma belcheri tsingtauense}",
        )
        self.assertEqual(alignment[0], "DHCCSEWGWCGRETSHCTCSSCV")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "                                                                                                                                                                                                                                                                ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                                               CCCCCCCCeEeCCCCCcCCcccc                                                                                                                                                                                          ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            "                                               6899999999999999973565                                                                                                                                                                                           ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 47,  65,  66,  70],
                             [143, 161, 161, 165]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['D', 'H', 'C', 'C', 'S', 'E', 'W', 'G', 'W', 'C', 'G', 'R', 'E',
              'T', 'S', 'H', 'C', 'T', 'C', 'S', 'S', 'C', 'V'],
             ['N', 'Y', 'C', 'C', 'S', 'K', 'W', 'G', 'S', 'C', 'G', 'I', 'G',
              'P', 'G', 'Y', 'C', 'G', '-', 'A', 'G', 'C', 'Q']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 80.97)
        self.assertAlmostEqual(alignment.annotations["E-value"], 0.11)
        self.assertAlmostEqual(alignment.annotations["Score"], 32.66)
        self.assertAlmostEqual(alignment.annotations["Identities"], 20)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.658)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 91.1)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[10:166],
            "ECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYCGAGCQS",
        )
        self.assertEqual(
            alignment[1],
            "ECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYCGAGCQS",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "          ~c~~~~CCs~~g~CG~~~~~c~~~c~~~~c~~~~~Cg~~~~~~~c~~~~CCs~~g~CG~~~~~c~~~c~~~~~~~~~~Cg~~~~~~~c~~~~CCS~~G~CG~~~~~C~~~Cq~~~c~~~~~Cg~~~~~~~c~~~~CCS~~G~CG~~~~~C~~gCq~     ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "          ~c~~~~cc~~~~~c~~~~~~c~~~c~~~~c~~~~~c~~~~~~c~~~~cc~~~~~c~~~~~~c~~~c~~~~c~~~~~c~~~~~~c~~~~cc~~~~~c~~~~~~c~~~c~~~~c~~~~~c~~~~~~c~~~~cc~~~~~c~~~~~~c~~~c~~    ",
        )
        self.assertEqual(alignment.target.id, "1wga")
        self.assertEqual(
            alignment.target.seq[10:160],
            "XCXXXXCCXXXXXCXXXXXXCXXXCXXXXCXXXXXCXXXXXXCXXXXCCXXXXXCXXXXXXCXXXCXXXXCXXXXXCXXXXXXCXXXXCCXXXXXCXXXXXXCXXXCXXXXCXXXXXCXXXXXXCXXXXCCXXXXXCXXXXXXCXXXCXX",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "1wga")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "; lectin (agglutinin); NMR {}",
        )
        self.assertEqual(
            alignment[0],
            "XCXXXXCCXXXXXCXXXXXXCXXXCXXXXCXXXXXCXXX--XXXCXXXXCCXXXXXCXXXXXXCXXXCXXXXCXXXXXCXXX--XXXCXXXXCCXXXXXCXXXXXXCXXXCXXXXCXXXXXCXXX--XXXCXXXXCCXXXXXCXXXXXXCXXXCXX",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "          cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            "          456667777777777777777665554434333344332235666778877777777777776665544333333343323456667788888888777777776665544433334443234566677888777787777777665543    ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 10,  49,  49,  90,  90, 131, 131, 160],
                             [ 10,  49,  51,  92,  94, 135, 137, 166]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['X', 'C', 'X', 'X', 'X', 'X', 'C', 'C', 'X', 'X', 'X', 'X', 'X',
              'C', 'X', 'X', 'X', 'X', 'X', 'X', 'C', 'X', 'X', 'X', 'C', 'X',
              'X', 'X', 'X', 'C', 'X', 'X', 'X', 'X', 'X', 'C', 'X', 'X', 'X',
              '-', '-', 'X', 'X', 'X', 'C', 'X', 'X', 'X', 'X', 'C', 'C', 'X',
              'X', 'X', 'X', 'X', 'C', 'X', 'X', 'X', 'X', 'X', 'X', 'C', 'X',
              'X', 'X', 'C', 'X', 'X', 'X', 'X', 'C', 'X', 'X', 'X', 'X', 'X',
              'C', 'X', 'X', 'X', '-', '-', 'X', 'X', 'X', 'C', 'X', 'X', 'X',
              'X', 'C', 'C', 'X', 'X', 'X', 'X', 'X', 'C', 'X', 'X', 'X', 'X',
              'X', 'X', 'C', 'X', 'X', 'X', 'C', 'X', 'X', 'X', 'X', 'C', 'X',
              'X', 'X', 'X', 'X', 'C', 'X', 'X', 'X', '-', '-', 'X', 'X', 'X',
              'C', 'X', 'X', 'X', 'X', 'C', 'C', 'X', 'X', 'X', 'X', 'X', 'C',
              'X', 'X', 'X', 'X', 'X', 'X', 'C', 'X', 'X', 'X', 'C', 'X', 'X'],
             ['E', 'C', 'P', 'N', 'N', 'L', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y',
              'C', 'G', 'M', 'G', 'G', 'D', 'Y', 'C', 'G', 'K', 'G', 'C', 'Q',
              'N', 'G', 'A', 'C', 'W', 'T', 'S', 'K', 'R', 'C', 'G', 'S', 'Q',
              'A', 'G', 'G', 'A', 'T', 'C', 'T', 'N', 'N', 'Q', 'C', 'C', 'S',
              'Q', 'Y', 'G', 'Y', 'C', 'G', 'F', 'G', 'A', 'E', 'Y', 'C', 'G',
              'A', 'G', 'C', 'Q', 'G', 'G', 'P', 'C', 'R', 'A', 'D', 'I', 'K',
              'C', 'G', 'S', 'Q', 'A', 'G', 'G', 'K', 'L', 'C', 'P', 'N', 'N',
              'L', 'C', 'C', 'S', 'Q', 'W', 'G', 'F', 'C', 'G', 'L', 'G', 'S',
              'E', 'F', 'C', 'G', 'G', 'G', 'C', 'Q', 'S', 'G', 'A', 'C', 'S',
              'T', 'D', 'K', 'P', 'C', 'G', 'K', 'D', 'A', 'G', 'G', 'R', 'V',
              'C', 'T', 'N', 'N', 'Y', 'C', 'C', 'S', 'K', 'W', 'G', 'S', 'C',
              'G', 'I', 'G', 'P', 'G', 'Y', 'C', 'G', 'A', 'G', 'C', 'Q', 'S']],
            dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 80.29)
        self.assertAlmostEqual(alignment.annotations["E-value"], 0.11)
        self.assertAlmostEqual(alignment.annotations["Score"], 36.76)
        self.assertAlmostEqual(alignment.annotations["Identities"], 37)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.967)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[1:36], "RCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQ"
        )
        self.assertEqual(alignment[1], "RCGEQG-----SNMECPN---NLCCSQYGYCGMGGDYCGK-GCQ")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            " ~cg~~~~~~~c~~~~CCs~~g~CG~~~~~c~~~c~                                                                                                                                       ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                          ~Cg~~~~~~~~~~~~C~~~~~~~CCs~~~~Cg~~~~~C~~~~c~                                                                                                                                                                                          ",
        )
        self.assertEqual(alignment.target.id, "4zxm_A")
        self.assertEqual(
            alignment.target.seq[26:70], "RCGPNYPAPDANPGECNPHAVDHCCSEWGWCGRETSHCTCSSCV"
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "4zxm_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "PGRP domain of peptidoglycan recognition protein; amidase, hydrolase; 2.80A {Branchiostoma belcheri tsingtauense}",
        )
        self.assertEqual(alignment[0], "RCGPNYPAPDANPGECNPHAVDHCCSEWGWCGRETSHCTCSSCV")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "                                                                                                                                                                                                                                                                ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                          CCCCCCCCCCCCCcccCCCCCCCCCCCCCeEeCCCCCcCCcccc                                                                                                                                                                                          ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            "                          57765423467657899999999999999974576                                                                                                                                                                                                   ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[26, 32, 37, 44, 47, 66, 67, 70],
                             [ 1,  7,  7, 14, 14, 33, 33, 36]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['R', 'C', 'G', 'P', 'N', 'Y', 'P', 'A', 'P', 'D', 'A', 'N', 'P',
              'G', 'E', 'C', 'N', 'P', 'H', 'A', 'V', 'D', 'H', 'C', 'C', 'S',
              'E', 'W', 'G', 'W', 'C', 'G', 'R', 'E', 'T', 'S', 'H', 'C', 'T',
              'C', 'S', 'S', 'C', 'V'],
             ['R', 'C', 'G', 'E', 'Q', 'G', '-', '-', '-', '-', '-', 'S', 'N',
              'M', 'E', 'C', 'P', 'N', '-', '-', '-', 'N', 'L', 'C', 'C', 'S',
              'Q', 'Y', 'G', 'Y', 'C', 'G', 'M', 'G', 'G', 'D', 'Y', 'C', 'G',
              'K', '-', 'G', 'C', 'Q']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 80.27)
        self.assertAlmostEqual(alignment.annotations["E-value"], 0.11)
        self.assertAlmostEqual(alignment.annotations["Score"], 36.29)
        self.assertAlmostEqual(alignment.annotations["Identities"], 38)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.185)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 19.7)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(alignment.query.seq[99:123], "NNLCCSQWGFCGLGSEFCGGGCQS")
        self.assertEqual(alignment[1], "NNLCCSQWGFCGLGSEFCGG-GCQS")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                                   ~~~CCS~~G~CG~~~~~C~~~Cq~                                                ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                             ~~~CCs~~~~cg~~~~~c~~~~c~d                                                                                                                                                                                      ",
        )
        self.assertEqual(alignment.target.id, "4z8i_A")
        self.assertEqual(alignment.target.seq[29:54], "VDHCCSEWGWCGRETSHCTCSSCVD")
        self.assertEqual(alignment.target.annotations["hmm_name"], "4z8i_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "BBTPGRP3, peptidoglycan recognition protein 3; chitin-binding domain, AM hydrolase; 2.70A {Branchiostoma belcheri tsingtauense}",
        )
        self.assertEqual(alignment[0], "VDHCCSEWGWCGRETSHCTCSSCVD")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "                             SCCEECTTSBEECSHHHHHSTTCEE                                                                                                                                                                                      ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                             CCCCCCCCCEEeCCcccccCCcccc                                                                                                                                                                                      ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            "                             368999999999999999835553                                                                                                                                                                                       ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 29,  49,  50,  54],
                             [ 99, 119, 119, 123]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['V', 'D', 'H', 'C', 'C', 'S', 'E', 'W', 'G', 'W', 'C', 'G', 'R',
              'E', 'T', 'S', 'H', 'C', 'T', 'C', 'S', 'S', 'C', 'V', 'D'],
             ['N', 'N', 'L', 'C', 'C', 'S', 'Q', 'W', 'G', 'F', 'C', 'G', 'L',
              'G', 'S', 'E', 'F', 'C', 'G', 'G', '-', 'G', 'C', 'Q', 'S']],
            dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 79.58)
        self.assertAlmostEqual(alignment.annotations["E-value"], 0.12)
        self.assertAlmostEqual(alignment.annotations["Score"], 36.06)
        self.assertAlmostEqual(alignment.annotations["Identities"], 35)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.927)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 27.4)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[0:37], "ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQN"
        )
        self.assertEqual(alignment[1], "ERCGEQG-----SNMECPN---NLCCSQYGYCGMGGDYCGK-GCQN")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "~~cg~~~~~~~c~~~~CCs~~g~CG~~~~~c~~~c~~                                                                                                                                      ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "        ~rCg~~~~~~~~~~~~C~~~~~~~CCs~~~~cg~~~~~c~~~~c~d                                                                                                                                                                                      ",
        )
        self.assertEqual(alignment.target.id, "4z8i_A")
        self.assertEqual(
            alignment.target.seq[8:54], "GRCGPNYPAPDANPGECNPHAVDHCCSEWGWCGRETSHCTCSSCVD"
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "4z8i_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "BBTPGRP3, peptidoglycan recognition protein 3; chitin-binding domain, AM hydrolase; 2.70A {Branchiostoma belcheri tsingtauense}",
        )
        self.assertEqual(alignment[0], "GRCGPNYPAPDANPGECNPHAVDHCCSEWGWCGRETSHCTCSSCVD")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "        SBCBSSSCBTTBSSBBCCTTSSCCEECTTSBEECSHHHHHSTTCEE                                                                                                                                                                                      ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "        CCCCCCCCCCCCCCcccCCCCCCCCCCCCCEEeCCcccccCCcccc                                                                                                                                                                                      ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            "        3677664334565378999999999999999845654                                                                                                                                                                                               ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 8, 15, 20, 27, 30, 49, 50, 54],
                             [ 0,  7,  7, 14, 14, 33, 33, 37]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['G', 'R', 'C', 'G', 'P', 'N', 'Y', 'P', 'A', 'P', 'D', 'A', 'N',
              'P', 'G', 'E', 'C', 'N', 'P', 'H', 'A', 'V', 'D', 'H', 'C', 'C',
              'S', 'E', 'W', 'G', 'W', 'C', 'G', 'R', 'E', 'T', 'S', 'H', 'C',
              'T', 'C', 'S', 'S', 'C', 'V', 'D'],
             ['E', 'R', 'C', 'G', 'E', 'Q', 'G', '-', '-', '-', '-', '-', 'S',
              'N', 'M', 'E', 'C', 'P', 'N', '-', '-', '-', 'N', 'L', 'C', 'C',
              'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'M', 'G', 'G', 'D', 'Y', 'C',
              'G', 'K', '-', 'G', 'C', 'Q', 'N']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 40.43)
        self.assertAlmostEqual(alignment.annotations["E-value"], 2.6)
        self.assertAlmostEqual(alignment.annotations["Score"], 25.90)
        self.assertAlmostEqual(alignment.annotations["Identities"], 20)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.652)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 54.7)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[53:163],
            "TCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYCGAG",
        )
        self.assertEqual(
            alignment[1],
            "TCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYCGAG",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                     ~c~~~~CCs~~g~CG~~~~~c~~~c~~~~~~~~~~Cg~~~~~~~c~~~~CCS~~G~CG~~~~~C~~~Cq~~~c~~~~~Cg~~~~~~~c~~~~CCS~~G~CG~~~~~C~~g        ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "          ~c~~~~cc~~~~~c~~~~~~c~~~c~~~~c~~~~~c~~~~~~c~~~~cc~~~~~c~~~~~~c~~~c~~~~c~~~~~c~~~~~~c~~~~cc~~~~~c~~~~~~c~~~                                                ",
        )
        self.assertEqual(alignment.target.id, "1wga")
        self.assertEqual(
            alignment.target.seq[10:116],
            "XCXXXXCCXXXXXCXXXXXXCXXXCXXXXCXXXXXCXXXXXXCXXXXCCXXXXXCXXXXXXCXXXCXXXXCXXXXXCXXXXXXCXXXXCCXXXXXCXXXXXXCXXX",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "1wga")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "; lectin (agglutinin); NMR {}",
        )
        self.assertEqual(
            alignment[0],
            "XCXXXXCCXXXXXCXXXXXXCXXXCXXXXCXXXXXCXXX--XXXCXXXXCCXXXXXCXXXXXXCXXXCXXXXCXXXXXCXXX--XXXCXXXXCCXXXXXCXXXXXXCXXX",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "          cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                                                ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            "          3445566666666666665665555433332233333212346666777777777777666666555443322233332223344455555555555555544433                                                ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 10,  49,  49,  90,  90, 116],
                             [ 53,  92,  94, 135, 137, 163]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['X', 'C', 'X', 'X', 'X', 'X', 'C', 'C', 'X', 'X', 'X', 'X', 'X',
              'C', 'X', 'X', 'X', 'X', 'X', 'X', 'C', 'X', 'X', 'X', 'C', 'X',
              'X', 'X', 'X', 'C', 'X', 'X', 'X', 'X', 'X', 'C', 'X', 'X', 'X',
              '-', '-', 'X', 'X', 'X', 'C', 'X', 'X', 'X', 'X', 'C', 'C', 'X',
              'X', 'X', 'X', 'X', 'C', 'X', 'X', 'X', 'X', 'X', 'X', 'C', 'X',
              'X', 'X', 'C', 'X', 'X', 'X', 'X', 'C', 'X', 'X', 'X', 'X', 'X',
              'C', 'X', 'X', 'X', '-', '-', 'X', 'X', 'X', 'C', 'X', 'X', 'X',
              'X', 'C', 'C', 'X', 'X', 'X', 'X', 'X', 'C', 'X', 'X', 'X', 'X',
              'X', 'X', 'C', 'X', 'X', 'X'],
             ['T', 'C', 'T', 'N', 'N', 'Q', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y',
              'C', 'G', 'F', 'G', 'A', 'E', 'Y', 'C', 'G', 'A', 'G', 'C', 'Q',
              'G', 'G', 'P', 'C', 'R', 'A', 'D', 'I', 'K', 'C', 'G', 'S', 'Q',
              'A', 'G', 'G', 'K', 'L', 'C', 'P', 'N', 'N', 'L', 'C', 'C', 'S',
              'Q', 'W', 'G', 'F', 'C', 'G', 'L', 'G', 'S', 'E', 'F', 'C', 'G',
              'G', 'G', 'C', 'Q', 'S', 'G', 'A', 'C', 'S', 'T', 'D', 'K', 'P',
              'C', 'G', 'K', 'D', 'A', 'G', 'G', 'R', 'V', 'C', 'T', 'N', 'N',
              'Y', 'C', 'C', 'S', 'K', 'W', 'G', 'S', 'C', 'G', 'I', 'G', 'P',
              'G', 'Y', 'C', 'G', 'A', 'G']], dtype='U')
                # fmt: on
            )
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_length(self):
        """Test getting the number of alignments without parsing the file."""
        stream = open(self.path)
        alignments = Align.parse(stream, "hhr")
        stream.close()
        self.assertEqual(len(alignments), 32)


class Align_hhr_2uvo_hhsearch(unittest.TestCase):
    path = os.path.join("HHsuite", "2uvo_hhsearch.hhr")

    def test_reading(self):
        alignments = Align.parse(self.path, "hhr")
        self.assertEqual(alignments.metadata["No_of_seqs"], (1, 4))
        self.assertAlmostEqual(alignments.metadata["Neff"], 1.0)
        self.assertEqual(alignments.metadata["Searched_HMMs"], 38388)
        self.assertEqual(alignments.metadata["Rundate"], "Fri Feb  1 13:49:32 2019")
        self.assertEqual(
            alignments.metadata["Command line"],
            "hhsearch -i 2uvo.fasta -d /pdb70_hhm_db",
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 100.00)
        self.assertAlmostEqual(alignment.annotations["E-value"], 4.6e-42)
        self.assertAlmostEqual(alignment.annotations["Score"], 249.39)
        self.assertAlmostEqual(alignment.annotations["Identities"], 100)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 2.050)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 166.9)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[0:171],
            "ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYCGAGCQSGGCDG",
        )
        self.assertEqual(
            alignment[1],
            "ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYCGAGCQSGGCDG",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "ercgeqgsnmecpnnlccsqygycgmggdycgkgcqngacwtskrcgsqaggatctnnqccsqygycgfgaeycgagcqggpcradikcgsqaggklcpnnlccsqwgfcglgsefcgggcqsgacstdkpcgkdaggrvctnnyccskwgscgigpgycgagcqsggcdg",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "~~cg~~~~~~~c~~~~CCS~~g~Cg~~~~~Cg~gC~~~~c~~~~~cg~~~~~~~c~~~~CCs~~g~Cg~~~~~c~~~c~~~~~~~~~~cg~~~~~~~c~~~~CCs~~g~CG~~~~~C~~gCq~~~c~~~~~cg~~~~~~~c~~~~ccs~~g~Cg~~~~~C~~~cq~~~~~~",
        )
        self.assertEqual(alignment.target.id, "2uvo_A")
        self.assertEqual(
            alignment.target.seq[0:171],
            "ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYCGAGCQSGGCDG",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "2uvo_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Agglutinin isolectin 1; carbohydrate-binding protein, hevein domain, chitin-binding, GERM agglutinin, chitin-binding protein; HET: NDG NAG GOL; 1.40A {Triticum aestivum} PDB: 1wgc_A* 2cwg_A* 2x3t_A* 4aml_A* 7wga_A 9wga_A 2wgc_A 1wgt_A 1k7t_A* 1k7v_A* 1k7u_A 2x52_A* 1t0w_A*",
        )
        self.assertEqual(
            alignment[0],
            "ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYCGAGCQSGGCDG",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "CBCBGGGTTBBCGGGCEECTTSBEEBSHHHHSTTCCBSSCSSCCBCBGGGTTBCCSTTCEECTTSBEEBSHHHHSTTCCBSSCSSCCBCBGGGTTBCCGGGCEECTTSBEEBSHHHHSTTCCBSSCSSCCCCBTTTTTBCCSTTCEECTTSCEEBSHHHHSTTCCBSSCC ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "CCCCCCCCCcCCCCCCeeCCCCeECCCcccccCCccccccccccccCcccCCcccCCccccCCCceeCCCccccCCCcccccccccccccccccCCCCCCCcccCCCCccCCCcccccCCCcCCccccccccccccccccCCCCCCcCCCCEecCchhhcccccccCCCCC",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            "799999999999999999999999999999999999999999999999999999999999999999999999999999999999999999998899999999999999999999999999999999999999999999999999999999999999999999999999986",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  0, 171],
                             [  0, 171]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['E', 'R', 'C', 'G', 'E', 'Q', 'G', 'S', 'N', 'M', 'E', 'C', 'P',
              'N', 'N', 'L', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'M',
              'G', 'G', 'D', 'Y', 'C', 'G', 'K', 'G', 'C', 'Q', 'N', 'G', 'A',
              'C', 'W', 'T', 'S', 'K', 'R', 'C', 'G', 'S', 'Q', 'A', 'G', 'G',
              'A', 'T', 'C', 'T', 'N', 'N', 'Q', 'C', 'C', 'S', 'Q', 'Y', 'G',
              'Y', 'C', 'G', 'F', 'G', 'A', 'E', 'Y', 'C', 'G', 'A', 'G', 'C',
              'Q', 'G', 'G', 'P', 'C', 'R', 'A', 'D', 'I', 'K', 'C', 'G', 'S',
              'Q', 'A', 'G', 'G', 'K', 'L', 'C', 'P', 'N', 'N', 'L', 'C', 'C',
              'S', 'Q', 'W', 'G', 'F', 'C', 'G', 'L', 'G', 'S', 'E', 'F', 'C',
              'G', 'G', 'G', 'C', 'Q', 'S', 'G', 'A', 'C', 'S', 'T', 'D', 'K',
              'P', 'C', 'G', 'K', 'D', 'A', 'G', 'G', 'R', 'V', 'C', 'T', 'N',
              'N', 'Y', 'C', 'C', 'S', 'K', 'W', 'G', 'S', 'C', 'G', 'I', 'G',
              'P', 'G', 'Y', 'C', 'G', 'A', 'G', 'C', 'Q', 'S', 'G', 'G', 'C',
              'D', 'G'],
             ['E', 'R', 'C', 'G', 'E', 'Q', 'G', 'S', 'N', 'M', 'E', 'C', 'P',
              'N', 'N', 'L', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'M',
              'G', 'G', 'D', 'Y', 'C', 'G', 'K', 'G', 'C', 'Q', 'N', 'G', 'A',
              'C', 'W', 'T', 'S', 'K', 'R', 'C', 'G', 'S', 'Q', 'A', 'G', 'G',
              'A', 'T', 'C', 'T', 'N', 'N', 'Q', 'C', 'C', 'S', 'Q', 'Y', 'G',
              'Y', 'C', 'G', 'F', 'G', 'A', 'E', 'Y', 'C', 'G', 'A', 'G', 'C',
              'Q', 'G', 'G', 'P', 'C', 'R', 'A', 'D', 'I', 'K', 'C', 'G', 'S',
              'Q', 'A', 'G', 'G', 'K', 'L', 'C', 'P', 'N', 'N', 'L', 'C', 'C',
              'S', 'Q', 'W', 'G', 'F', 'C', 'G', 'L', 'G', 'S', 'E', 'F', 'C',
              'G', 'G', 'G', 'C', 'Q', 'S', 'G', 'A', 'C', 'S', 'T', 'D', 'K',
              'P', 'C', 'G', 'K', 'D', 'A', 'G', 'G', 'R', 'V', 'C', 'T', 'N',
              'N', 'Y', 'C', 'C', 'S', 'K', 'W', 'G', 'S', 'C', 'G', 'I', 'G',
              'P', 'G', 'Y', 'C', 'G', 'A', 'G', 'C', 'Q', 'S', 'G', 'G', 'C',
              'D', 'G']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 99.95)
        self.assertAlmostEqual(alignment.annotations["E-value"], 2.8e-33)
        self.assertAlmostEqual(alignment.annotations["Score"], 204.56)
        self.assertAlmostEqual(alignment.annotations["Identities"], 49)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.252)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 153.2)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[1:169],
            "RCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYCGAGCQSGGC",
        )
        self.assertEqual(
            alignment[1],
            "RCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYCGAGCQSGGC",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            " rcgeqgsnmecpnnlccsqygycgmggdycgkgcqngacwtskrcgsqaggatctnnqccsqygycgfgaeycgagcqggpcradikcgsqaggklcpnnlccsqwgfcglgsefcgggcqsgacstdkpcgkdaggrvctnnyccskwgscgigpgycgagcqsggc  ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            " ~~g~~~~~~~c~~~~CCs~~g~Cg~~~~~c~~~C~~~~c~~~~~cg~~~~~c~~~~CCs~~G~CG~~~~~C~~~C~~~~~~~~~~Cg~~~~~c~~~~Ccs~~G~CGt~~~~C~~~cq~~~c~~~~~cg~~~~~c~~~~Ccs~~g~Cg~~~~~C~~~cq~~~~ ",
        )
        self.assertEqual(alignment.target.id, "2wga")
        self.assertEqual(
            alignment.target.seq[1:163],
            "GXGCXGXXMYCSTNNCCXXWESCGSGGYXCGEGCNLGACQXGXPCXXPGSVCTNLHCCARGGHCGMGSGYCGXGCXGGACXADIXCGXGXXXCPTDSCCGGWGXCGNGXEFCGXGCXVGGCAAXSPCGXPGSXCTLDKCCSGXGACXSGSGGCGXGCXAGGC",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "2wga")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "; lectin (agglutinin); NMR {}",
        )
        self.assertEqual(
            alignment[0],
            "GXGCXGXXMYCSTNNCCXXWESCGSGGYXCGEGCNLGACQXGXPCXX--PGSVCTNLHCCARGGHCGMGSGYCGXGCXGGACXADIXCGXG--XXXCPTDSCCGGWGXCGNGXEFCGXGCXVGGCAAXSPCGX--PGSXCTLDKCCSGXGACXSGSGGCGXGCXAGGC",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            " CSGGGSSCCCCSTTCEECTTSCEECSTTTTSTTCCSSSCSSCCCSSSSSCCCSTTCEECTTSCEESSHHHHSSCCSSSSCSSCCCCTTSSSCCSTTCBCCSSSCCBCSHHHHSTTCCSSSCSSCCCCSSSCCCCSTTCEECSSSSEECSTTTTSSCCSSSSC ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            " ccCCCcccccCCCCceECCcceECCCCccccCccccCcccccceeccCCCcCCCCcccCCCceeCCCCcccCCCccccccccccccCcccccCCCCCccCCCCCccCccccccCCccccccccccccCCCcccCCcccccCCCCceeCCccccCCCCcCCCC ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            " 588899999999999999999999999999999999999999999985789999999999999999999999999999999999999855679999999999999999999999999999999999988367899999999999999999999999999876 ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  1,  48,  48,  90,  90, 130, 130, 163],
                             [  1,  48,  50,  92,  94, 134, 136, 169]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['G', 'X', 'G', 'C', 'X', 'G', 'X', 'X', 'M', 'Y', 'C', 'S', 'T',
              'N', 'N', 'C', 'C', 'X', 'X', 'W', 'E', 'S', 'C', 'G', 'S', 'G',
              'G', 'Y', 'X', 'C', 'G', 'E', 'G', 'C', 'N', 'L', 'G', 'A', 'C',
              'Q', 'X', 'G', 'X', 'P', 'C', 'X', 'X', '-', '-', 'P', 'G', 'S',
              'V', 'C', 'T', 'N', 'L', 'H', 'C', 'C', 'A', 'R', 'G', 'G', 'H',
              'C', 'G', 'M', 'G', 'S', 'G', 'Y', 'C', 'G', 'X', 'G', 'C', 'X',
              'G', 'G', 'A', 'C', 'X', 'A', 'D', 'I', 'X', 'C', 'G', 'X', 'G',
              '-', '-', 'X', 'X', 'X', 'C', 'P', 'T', 'D', 'S', 'C', 'C', 'G',
              'G', 'W', 'G', 'X', 'C', 'G', 'N', 'G', 'X', 'E', 'F', 'C', 'G',
              'X', 'G', 'C', 'X', 'V', 'G', 'G', 'C', 'A', 'A', 'X', 'S', 'P',
              'C', 'G', 'X', '-', '-', 'P', 'G', 'S', 'X', 'C', 'T', 'L', 'D',
              'K', 'C', 'C', 'S', 'G', 'X', 'G', 'A', 'C', 'X', 'S', 'G', 'S',
              'G', 'G', 'C', 'G', 'X', 'G', 'C', 'X', 'A', 'G', 'G', 'C'],
             ['R', 'C', 'G', 'E', 'Q', 'G', 'S', 'N', 'M', 'E', 'C', 'P', 'N',
              'N', 'L', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'M', 'G',
              'G', 'D', 'Y', 'C', 'G', 'K', 'G', 'C', 'Q', 'N', 'G', 'A', 'C',
              'W', 'T', 'S', 'K', 'R', 'C', 'G', 'S', 'Q', 'A', 'G', 'G', 'A',
              'T', 'C', 'T', 'N', 'N', 'Q', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y',
              'C', 'G', 'F', 'G', 'A', 'E', 'Y', 'C', 'G', 'A', 'G', 'C', 'Q',
              'G', 'G', 'P', 'C', 'R', 'A', 'D', 'I', 'K', 'C', 'G', 'S', 'Q',
              'A', 'G', 'G', 'K', 'L', 'C', 'P', 'N', 'N', 'L', 'C', 'C', 'S',
              'Q', 'W', 'G', 'F', 'C', 'G', 'L', 'G', 'S', 'E', 'F', 'C', 'G',
              'G', 'G', 'C', 'Q', 'S', 'G', 'A', 'C', 'S', 'T', 'D', 'K', 'P',
              'C', 'G', 'K', 'D', 'A', 'G', 'G', 'R', 'V', 'C', 'T', 'N', 'N',
              'Y', 'C', 'C', 'S', 'K', 'W', 'G', 'S', 'C', 'G', 'I', 'G', 'P',
              'G', 'Y', 'C', 'G', 'A', 'G', 'C', 'Q', 'S', 'G', 'G', 'C']],
            dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 99.84)
        self.assertAlmostEqual(alignment.annotations["E-value"], 1.1e-25)
        self.assertAlmostEqual(alignment.annotations["Score"], 163.39)
        self.assertAlmostEqual(alignment.annotations["Identities"], 60)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.533)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 121.1)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[43:169],
            "KRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYCGAGCQSGGC",
        )
        self.assertEqual(
            alignment[1],
            "KRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYCGAGCQSGGC",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                           krcgsqaggatctnnqccsqygycgfgaeycgagcqggpcradikcgsqaggklcpnnlccsqwgfcglgsefcgggcqsgacstdkpcgkdaggrvctnnyccskwgscgigpgycgagcqsggc  ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "~~cg~~~~~~~c~~~~CCS~~g~Cg~~~~~Cg~gC~~~~c~~~~~cg~~~~~~~c~~~~CCs~~g~Cg~~~~~c~~~c~~~~~~~~~~cg~~~~~~~c~~~~CCs~~g~CG~~~~~C~~gCq~~~c                                             ",
        )
        self.assertEqual(alignment.target.id, "2uvo_A")
        self.assertEqual(
            alignment.target.seq[0:126],
            "ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGAC",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "2uvo_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Agglutinin isolectin 1; carbohydrate-binding protein, hevein domain, chitin-binding, GERM agglutinin, chitin-binding protein; HET: NDG NAG GOL; 1.40A {Triticum aestivum} PDB: 1wgc_A* 2cwg_A* 2x3t_A* 4aml_A* 7wga_A 9wga_A 2wgc_A 1wgt_A 1k7t_A* 1k7v_A* 1k7u_A 2x52_A* 1t0w_A*",
        )
        self.assertEqual(
            alignment[0],
            "ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGAC",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "CBCBGGGTTBBCGGGCEECTTSBEEBSHHHHSTTCCBSSCSSCCBCBGGGTTBCCSTTCEECTTSBEEBSHHHHSTTCCBSSCSSCCBCBGGGTTBCCGGGCEECTTSBEEBSHHHHSTTCCBSSC                                             ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "CCCCCCCCCcCCCCCCeeCCCCeECCCcccccCCccccccccccccCcccCCcccCCccccCCCceeCCCccccCCCcccccccccccccccccCCCCCCCcccCCCCccCCCcccccCCCcCCcc                                             ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            "489999999999999999999999999999999999999999999999998899999999999999999999999999999999999999988889999999999999999999999999998644                                             ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  0, 126],
                             [ 43, 169]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['E', 'R', 'C', 'G', 'E', 'Q', 'G', 'S', 'N', 'M', 'E', 'C', 'P',
              'N', 'N', 'L', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'M',
              'G', 'G', 'D', 'Y', 'C', 'G', 'K', 'G', 'C', 'Q', 'N', 'G', 'A',
              'C', 'W', 'T', 'S', 'K', 'R', 'C', 'G', 'S', 'Q', 'A', 'G', 'G',
              'A', 'T', 'C', 'T', 'N', 'N', 'Q', 'C', 'C', 'S', 'Q', 'Y', 'G',
              'Y', 'C', 'G', 'F', 'G', 'A', 'E', 'Y', 'C', 'G', 'A', 'G', 'C',
              'Q', 'G', 'G', 'P', 'C', 'R', 'A', 'D', 'I', 'K', 'C', 'G', 'S',
              'Q', 'A', 'G', 'G', 'K', 'L', 'C', 'P', 'N', 'N', 'L', 'C', 'C',
              'S', 'Q', 'W', 'G', 'F', 'C', 'G', 'L', 'G', 'S', 'E', 'F', 'C',
              'G', 'G', 'G', 'C', 'Q', 'S', 'G', 'A', 'C'],
             ['K', 'R', 'C', 'G', 'S', 'Q', 'A', 'G', 'G', 'A', 'T', 'C', 'T',
              'N', 'N', 'Q', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'F',
              'G', 'A', 'E', 'Y', 'C', 'G', 'A', 'G', 'C', 'Q', 'G', 'G', 'P',
              'C', 'R', 'A', 'D', 'I', 'K', 'C', 'G', 'S', 'Q', 'A', 'G', 'G',
              'K', 'L', 'C', 'P', 'N', 'N', 'L', 'C', 'C', 'S', 'Q', 'W', 'G',
              'F', 'C', 'G', 'L', 'G', 'S', 'E', 'F', 'C', 'G', 'G', 'G', 'C',
              'Q', 'S', 'G', 'A', 'C', 'S', 'T', 'D', 'K', 'P', 'C', 'G', 'K',
              'D', 'A', 'G', 'G', 'R', 'V', 'C', 'T', 'N', 'N', 'Y', 'C', 'C',
              'S', 'K', 'W', 'G', 'S', 'C', 'G', 'I', 'G', 'P', 'G', 'Y', 'C',
              'G', 'A', 'G', 'C', 'Q', 'S', 'G', 'G', 'C']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 99.84)
        self.assertAlmostEqual(alignment.annotations["E-value"], 1.5e-25)
        self.assertAlmostEqual(alignment.annotations["Score"], 157.96)
        self.assertAlmostEqual(alignment.annotations["Identities"], 52)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.251)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 110.7)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[43:167],
            "KRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYCGAGCQSG",
        )
        self.assertEqual(
            alignment[1],
            "KRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYCGAGCQSG",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                           krcgsqaggatctnnqccsqygycgfgaeycgagcqggpcradikcgsqaggklcpnnlccsqwgfcglgsefcgggcqsgacstdkpcgkdaggrvctnnyccskwgscgigpgycgagcqsg    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            " ~~cg~~~~~~~c~~g~CCs~~g~CG~~~~~Cg~gCq~~c~~~~Cg~~~~~~~c~~~~CCs~~G~CG~~~~~C~~~Cq~~c~~~~cg~~~~~~~c~~~~CCs~~g~CG~~~~~C~~gCq~~     ",
        )
        self.assertEqual(alignment.target.id, "1ulk_A")
        self.assertEqual(
            alignment.target.seq[1:121],
            "PVCGVRASGRVCPDGYCCSQWGYCGTTEEYCGKGCQSQCDYNRCGKEFGGKECHDELCCSQYGWCGNSDGHCGEGCQSQCSYWRCGKDFGGRLCTEDMCCSQYGWCGLTDDHCEDGCQSQ",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "1ulk_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Lectin-C; chitin-binding protein, hevein domain, PL-C, sugar binding protein; 1.80A {Phytolacca americana} SCOP: g.3.1.1 g.3.1.1 g.3.1.1",
        )
        self.assertEqual(
            alignment[0],
            "PVCGVRASGRVCPDGYCCSQWGYCGTTEEYCGKGCQSQ-CD-YNRCGKEFGGKECHDELCCSQYGWCGNSDGHCGEGCQSQ-C-SYWRCGKDFGGRLCTEDMCCSQYGWCGLTDDHCEDGCQSQ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            " CCSBGGGTTBCCGGGCEECTTSCEESSHHHHSTTCCBCTTTTBCBGGGTTBCCGGGCEECTTSBEECSHHHHSTTCCBCTTTTBCBGGGTTBCCSTTCEECTTSBEECSHHHHSTTCCBC     ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            " CcCcCCCCCCCCCCCCeECCCCeeCCCccccCCCccccceeeeccccccCCCCCCccccCCCcccccCcccccCCcccccCcccccccCCCccCCCCcccccCceecCcccccCcccccc     ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            " 699999999999999999999999999999999999764355788888889999999999999999999999999998633578999888999999999999999999999999999984     ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  1,  39,  39,  41,  41,  80,  80,  81,  81, 121],
                             [ 43,  81,  82,  84,  85, 124, 125, 126, 127, 167]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['P', 'V', 'C', 'G', 'V', 'R', 'A', 'S', 'G', 'R', 'V', 'C', 'P',
              'D', 'G', 'Y', 'C', 'C', 'S', 'Q', 'W', 'G', 'Y', 'C', 'G', 'T',
              'T', 'E', 'E', 'Y', 'C', 'G', 'K', 'G', 'C', 'Q', 'S', 'Q', '-',
              'C', 'D', '-', 'Y', 'N', 'R', 'C', 'G', 'K', 'E', 'F', 'G', 'G',
              'K', 'E', 'C', 'H', 'D', 'E', 'L', 'C', 'C', 'S', 'Q', 'Y', 'G',
              'W', 'C', 'G', 'N', 'S', 'D', 'G', 'H', 'C', 'G', 'E', 'G', 'C',
              'Q', 'S', 'Q', '-', 'C', '-', 'S', 'Y', 'W', 'R', 'C', 'G', 'K',
              'D', 'F', 'G', 'G', 'R', 'L', 'C', 'T', 'E', 'D', 'M', 'C', 'C',
              'S', 'Q', 'Y', 'G', 'W', 'C', 'G', 'L', 'T', 'D', 'D', 'H', 'C',
              'E', 'D', 'G', 'C', 'Q', 'S', 'Q'],
             ['K', 'R', 'C', 'G', 'S', 'Q', 'A', 'G', 'G', 'A', 'T', 'C', 'T',
              'N', 'N', 'Q', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'F',
              'G', 'A', 'E', 'Y', 'C', 'G', 'A', 'G', 'C', 'Q', 'G', 'G', 'P',
              'C', 'R', 'A', 'D', 'I', 'K', 'C', 'G', 'S', 'Q', 'A', 'G', 'G',
              'K', 'L', 'C', 'P', 'N', 'N', 'L', 'C', 'C', 'S', 'Q', 'W', 'G',
              'F', 'C', 'G', 'L', 'G', 'S', 'E', 'F', 'C', 'G', 'G', 'G', 'C',
              'Q', 'S', 'G', 'A', 'C', 'S', 'T', 'D', 'K', 'P', 'C', 'G', 'K',
              'D', 'A', 'G', 'G', 'R', 'V', 'C', 'T', 'N', 'N', 'Y', 'C', 'C',
              'S', 'K', 'W', 'G', 'S', 'C', 'G', 'I', 'G', 'P', 'G', 'Y', 'C',
              'G', 'A', 'G', 'C', 'Q', 'S', 'G']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 99.82)
        self.assertAlmostEqual(alignment.annotations["E-value"], 6.5e-25)
        self.assertAlmostEqual(alignment.annotations["Score"], 154.69)
        self.assertAlmostEqual(alignment.annotations["Identities"], 50)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.299)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 109.2)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[0:123],
            "ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQS",
        )
        self.assertEqual(
            alignment[1],
            "ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQS",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "ercgeqgsnmecpnnlccsqygycgmggdycgkgcqngacwtskrcgsqaggatctnnqccsqygycgfgaeycgagcqggpcradikcgsqaggklcpnnlccsqwgfcglgsefcgggcqs                                                ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            " ~~cg~~~~~~~c~~g~CCs~~g~CG~~~~~Cg~gCq~~c~~~~Cg~~~~~~~c~~~~CCs~~G~CG~~~~~C~~~Cq~~c~~~~cg~~~~~~~c~~~~CCs~~g~CG~~~~~C~~gCq~      ",
        )
        self.assertEqual(alignment.target.id, "1ulk_A")
        self.assertEqual(
            alignment.target.seq[1:120],
            "PVCGVRASGRVCPDGYCCSQWGYCGTTEEYCGKGCQSQCDYNRCGKEFGGKECHDELCCSQYGWCGNSDGHCGEGCQSQCSYWRCGKDFGGRLCTEDMCCSQYGWCGLTDDHCEDGCQS",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "1ulk_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Lectin-C; chitin-binding protein, hevein domain, PL-C, sugar binding protein; 1.80A {Phytolacca americana} SCOP: g.3.1.1 g.3.1.1 g.3.1.1",
        )
        self.assertEqual(
            alignment[0],
            "PVCGVRASGRVCPDGYCCSQWGYCGTTEEYCGKGCQSQ-CD-YNRCGKEFGGKECHDELCCSQYGWCGNSDGHCGEGCQSQ-CS-YWRCGKDFGGRLCTEDMCCSQYGWCGLTDDHCEDGCQS",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            " CCSBGGGTTBCCGGGCEECTTSCEESSHHHHSTTCCBCTTTTBCBGGGTTBCCGGGCEECTTSBEECSHHHHSTTCCBCTTTTBCBGGGTTBCCSTTCEECTTSBEECSHHHHSTTCCB      ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            " CcCcCCCCCCCCCCCCeECCCCeeCCCccccCCCccccceeeeccccccCCCCCCccccCCCcccccCcccccCCcccccCcccccccCCCccCCCCcccccCceecCcccccCccccc      ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            " 48999999999999999999999999999999999976444589988889999999999999999999999999999753356888888888999999999999999999999999997      ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  1,  39,  39,  41,  41,  80,  80,  82,  82, 120],
                             [  0,  38,  39,  41,  42,  81,  82,  84,  85, 123]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['P', 'V', 'C', 'G', 'V', 'R', 'A', 'S', 'G', 'R', 'V', 'C', 'P',
              'D', 'G', 'Y', 'C', 'C', 'S', 'Q', 'W', 'G', 'Y', 'C', 'G', 'T',
              'T', 'E', 'E', 'Y', 'C', 'G', 'K', 'G', 'C', 'Q', 'S', 'Q', '-',
              'C', 'D', '-', 'Y', 'N', 'R', 'C', 'G', 'K', 'E', 'F', 'G', 'G',
              'K', 'E', 'C', 'H', 'D', 'E', 'L', 'C', 'C', 'S', 'Q', 'Y', 'G',
              'W', 'C', 'G', 'N', 'S', 'D', 'G', 'H', 'C', 'G', 'E', 'G', 'C',
              'Q', 'S', 'Q', '-', 'C', 'S', '-', 'Y', 'W', 'R', 'C', 'G', 'K',
              'D', 'F', 'G', 'G', 'R', 'L', 'C', 'T', 'E', 'D', 'M', 'C', 'C',
              'S', 'Q', 'Y', 'G', 'W', 'C', 'G', 'L', 'T', 'D', 'D', 'H', 'C',
              'E', 'D', 'G', 'C', 'Q', 'S'],
             ['E', 'R', 'C', 'G', 'E', 'Q', 'G', 'S', 'N', 'M', 'E', 'C', 'P',
              'N', 'N', 'L', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'M',
              'G', 'G', 'D', 'Y', 'C', 'G', 'K', 'G', 'C', 'Q', 'N', 'G', 'A',
              'C', 'W', 'T', 'S', 'K', 'R', 'C', 'G', 'S', 'Q', 'A', 'G', 'G',
              'A', 'T', 'C', 'T', 'N', 'N', 'Q', 'C', 'C', 'S', 'Q', 'Y', 'G',
              'Y', 'C', 'G', 'F', 'G', 'A', 'E', 'Y', 'C', 'G', 'A', 'G', 'C',
              'Q', 'G', 'G', 'P', 'C', 'R', 'A', 'D', 'I', 'K', 'C', 'G', 'S',
              'Q', 'A', 'G', 'G', 'K', 'L', 'C', 'P', 'N', 'N', 'L', 'C', 'C',
              'S', 'Q', 'W', 'G', 'F', 'C', 'G', 'L', 'G', 'S', 'E', 'F', 'C',
              'G', 'G', 'G', 'C', 'Q', 'S']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 99.78)
        self.assertAlmostEqual(alignment.annotations["E-value"], 1.7e-23)
        self.assertAlmostEqual(alignment.annotations["Score"], 152.97)
        self.assertAlmostEqual(alignment.annotations["Identities"], 41)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.180)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 114.3)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[44:169],
            "RCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYCGAGCQSGGC",
        )
        self.assertEqual(
            alignment[1],
            "RCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYCGAGCQSGGC",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                            rcgsqaggatctnnqccsqygycgfgaeycgagcqggpcradikcgsqaggklcpnnlccsqwgfcglgsefcgggcqsgacstdkpcgkdaggrvctnnyccskwgscgigpgycgagcqsggc  ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            " ~~g~~~~~~~c~~~~CCs~~g~Cg~~~~~c~~~C~~~~c~~~~~cg~~~~~c~~~~CCs~~G~CG~~~~~C~~~C~~~~~~~~~~Cg~~~~~c~~~~Ccs~~G~CGt~~~~C~~~cq~~~c                                          ",
        )
        self.assertEqual(alignment.target.id, "2wga")
        self.assertEqual(
            alignment.target.seq[1:122],
            "GXGCXGXXMYCSTNNCCXXWESCGSGGYXCGEGCNLGACQXGXPCXXPGSVCTNLHCCARGGHCGMGSGYCGXGCXGGACXADIXCGXGXXXCPTDSCCGGWGXCGNGXEFCGXGCXVGGC",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "2wga")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "; lectin (agglutinin); NMR {}",
        )
        self.assertEqual(
            alignment[0],
            "GXGCXGXXMYCSTNNCCXXWESCGSGGYXCGEGCNLGACQXGXPCXX--PGSVCTNLHCCARGGHCGMGSGYCGXGCXGGACXADIXCGXG--XXXCPTDSCCGGWGXCGNGXEFCGXGCXVGGC",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            " CSGGGSSCCCCSTTCEECTTSCEECSTTTTSTTCCSSSCSSCCCSSSSSCCCSTTCEECTTSCEESSHHHHSSCCSSSSCSSCCCCTTSSSCCSTTCBCCSSSCCBCSHHHHSTTCCSSSC                                          ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            " ccCCCcccccCCCCceECCcceECCCCccccCccccCcccccceeccCCCcCCCCcccCCCceeCCCCcccCCCccccccccccccCcccccCCCCCccCCCCCccCccccccCCcccccc                                          ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            " 4889999999999999999999999999999999999999999999857889999999999999999999999999999999999998667899999999999999999999999998654                                          ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  1,  48,  48,  90,  90, 122],
                             [ 44,  91,  93, 135, 137, 169]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['G', 'X', 'G', 'C', 'X', 'G', 'X', 'X', 'M', 'Y', 'C', 'S', 'T',
              'N', 'N', 'C', 'C', 'X', 'X', 'W', 'E', 'S', 'C', 'G', 'S', 'G',
              'G', 'Y', 'X', 'C', 'G', 'E', 'G', 'C', 'N', 'L', 'G', 'A', 'C',
              'Q', 'X', 'G', 'X', 'P', 'C', 'X', 'X', '-', '-', 'P', 'G', 'S',
              'V', 'C', 'T', 'N', 'L', 'H', 'C', 'C', 'A', 'R', 'G', 'G', 'H',
              'C', 'G', 'M', 'G', 'S', 'G', 'Y', 'C', 'G', 'X', 'G', 'C', 'X',
              'G', 'G', 'A', 'C', 'X', 'A', 'D', 'I', 'X', 'C', 'G', 'X', 'G',
              '-', '-', 'X', 'X', 'X', 'C', 'P', 'T', 'D', 'S', 'C', 'C', 'G',
              'G', 'W', 'G', 'X', 'C', 'G', 'N', 'G', 'X', 'E', 'F', 'C', 'G',
              'X', 'G', 'C', 'X', 'V', 'G', 'G', 'C'],
             ['R', 'C', 'G', 'S', 'Q', 'A', 'G', 'G', 'A', 'T', 'C', 'T', 'N',
              'N', 'Q', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'F', 'G',
              'A', 'E', 'Y', 'C', 'G', 'A', 'G', 'C', 'Q', 'G', 'G', 'P', 'C',
              'R', 'A', 'D', 'I', 'K', 'C', 'G', 'S', 'Q', 'A', 'G', 'G', 'K',
              'L', 'C', 'P', 'N', 'N', 'L', 'C', 'C', 'S', 'Q', 'W', 'G', 'F',
              'C', 'G', 'L', 'G', 'S', 'E', 'F', 'C', 'G', 'G', 'G', 'C', 'Q',
              'S', 'G', 'A', 'C', 'S', 'T', 'D', 'K', 'P', 'C', 'G', 'K', 'D',
              'A', 'G', 'G', 'R', 'V', 'C', 'T', 'N', 'N', 'Y', 'C', 'C', 'S',
              'K', 'W', 'G', 'S', 'C', 'G', 'I', 'G', 'P', 'G', 'Y', 'C', 'G',
              'A', 'G', 'C', 'Q', 'S', 'G', 'G', 'C']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 99.54)
        self.assertAlmostEqual(alignment.annotations["E-value"], 1.1e-18)
        self.assertAlmostEqual(alignment.annotations["Score"], 113.58)
        self.assertAlmostEqual(alignment.annotations["Identities"], 48)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.312)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 73.1)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[43:124],
            "KRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSG",
        )
        self.assertEqual(
            alignment[1],
            "KRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSG",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                           krcgsqaggatctnnqccsqygycgfgaeycgagcqggpcradikcgsqaggklcpnnlccsqwgfcglgsefcgggcqsg                                               ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            " ~~cg~~~~~~~C~~g~CCs~~G~Cg~~~~~c~~~c~~~~~~g~Cg~~~~~~~c~~~~CCs~~g~CG~~~~~C~~~Cqs~  ",
        )
        self.assertEqual(alignment.target.id, "1uha_A")
        self.assertEqual(
            alignment.target.seq[1:80],
            "PECGERASGKRCPNGKCCSQWGYCGTTDNYCGQGCQSQCDYWRCGRDFGGRLCEEDMCCSKYGWCGYSDDHCEDGCQSQ",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "1uha_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Lectin-D2; chitin-binding domain, sugar binding protein; 1.50A {Phytolacca americana} SCOP: g.3.1.1 g.3.1.1 PDB: 1ulm_A* 1uln_A",
        )
        self.assertEqual(
            alignment[0],
            "PECGERASGKRCPNGKCCSQWGYCGTTDNYCGQGCQSQ-C-DYWRCGRDFGGRLCEEDMCCSKYGWCGYSDDHCEDGCQSQ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            " CCSBGGGTTBCCGGGCEECTTSCEESSHHHHSTTCCBCTTTTBCBGGGTTBCCSTTCEECTTSBEECSHHHHSTTCCBC  ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            " CCCCCCCCCCcCCCCCccCCCccccCccccccCCccccccccccccccceecCCCCCccCCCccccCCccccccccccc  ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            " 5899999999999999999999999999999999996534678998888999999999999999999999999999984  ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  1,  39,  39,  40,  40,  80],
                             [ 43,  81,  82,  83,  84, 124]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['P', 'E', 'C', 'G', 'E', 'R', 'A', 'S', 'G', 'K', 'R', 'C', 'P',
              'N', 'G', 'K', 'C', 'C', 'S', 'Q', 'W', 'G', 'Y', 'C', 'G', 'T',
              'T', 'D', 'N', 'Y', 'C', 'G', 'Q', 'G', 'C', 'Q', 'S', 'Q', '-',
              'C', '-', 'D', 'Y', 'W', 'R', 'C', 'G', 'R', 'D', 'F', 'G', 'G',
              'R', 'L', 'C', 'E', 'E', 'D', 'M', 'C', 'C', 'S', 'K', 'Y', 'G',
              'W', 'C', 'G', 'Y', 'S', 'D', 'D', 'H', 'C', 'E', 'D', 'G', 'C',
              'Q', 'S', 'Q'],
             ['K', 'R', 'C', 'G', 'S', 'Q', 'A', 'G', 'G', 'A', 'T', 'C', 'T',
              'N', 'N', 'Q', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'F',
              'G', 'A', 'E', 'Y', 'C', 'G', 'A', 'G', 'C', 'Q', 'G', 'G', 'P',
              'C', 'R', 'A', 'D', 'I', 'K', 'C', 'G', 'S', 'Q', 'A', 'G', 'G',
              'K', 'L', 'C', 'P', 'N', 'N', 'L', 'C', 'C', 'S', 'Q', 'W', 'G',
              'F', 'C', 'G', 'L', 'G', 'S', 'E', 'F', 'C', 'G', 'G', 'G', 'C',
              'Q', 'S', 'G']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 99.54)
        self.assertAlmostEqual(alignment.annotations["E-value"], 1.2e-18)
        self.assertAlmostEqual(alignment.annotations["Score"], 115.94)
        self.assertAlmostEqual(alignment.annotations["Identities"], 47)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.232)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 74.1)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[44:124],
            "RCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSG",
        )
        self.assertEqual(
            alignment[1],
            "RCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRA----DIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGG-CQSG",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                            rcgsqaggatctnnqccsqygycgfgaeycgagcqggpcradikcgsqaggklcpnnlccsqwgfcglgsefcgggcqsg                                               ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            " ~cg~~~~~~~C~~~~CCS~~G~CG~~~~~C~~~Cq~~c~~~~~~~~~Cg~~~~~~~c~~~~CCS~~G~CG~~~~~C~~~~Cq~~    ",
        )
        self.assertEqual(alignment.target.id, "1en2_A")
        self.assertEqual(
            alignment.target.seq[1:85],
            "RCGSQGGGSTCPGLRCCSIWGWCGDSEPYCGRTCENKCWSGERSDHRCGAAVGNPPCGQDRCCSVHGWCGGGNDYCSGGNCQYR",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "1en2_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "UDA, agglutinin isolectin I/agglutinin isolectin V/ AG isolectin VI; hevein domain, superantigen, saccharide binding binding protein; HET: NAG; 1.40A {Urtica dioica} SCOP: g.3.1.1 g.3.1.1 PDB: 1eis_A* 1enm_A* 1ehd_A 1ehh_A* 1iqb_A",
        )
        self.assertEqual(
            alignment[0],
            "RCGSQGGGSTCPGLRCCSIWGWCGDSEPYCGRTCENK-CWSGERSDHRCGAAVGNPPCGQDRCCSVHGWCGGGNDYCSGGNCQYR",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            " BCTTTTTSCCCGGGCEEETTSBEESSHHHHSTTEEESCGGGCCTTCBCSGGGTCCCCCTTCEEETTSBEESSHHHHSGGGEEEC    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            " CCCCCCCCccCCCCcccCCCceecccccccCCCCcCCCcccccCCcccCCcccccccCCCCeECCCceECCCccccCCCCcccC    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            " 7999999999999999999999999999999999976665889999888899999999999999999999999989985         ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  1,  38,  38,  41,  45,  80,  81,  85],
                             [ 44,  81,  82,  85,  85, 120, 120, 124]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['R', 'C', 'G', 'S', 'Q', 'G', 'G', 'G', 'S', 'T', 'C', 'P', 'G',
              'L', 'R', 'C', 'C', 'S', 'I', 'W', 'G', 'W', 'C', 'G', 'D', 'S',
              'E', 'P', 'Y', 'C', 'G', 'R', 'T', 'C', 'E', 'N', 'K', '-', 'C',
              'W', 'S', 'G', 'E', 'R', 'S', 'D', 'H', 'R', 'C', 'G', 'A', 'A',
              'V', 'G', 'N', 'P', 'P', 'C', 'G', 'Q', 'D', 'R', 'C', 'C', 'S',
              'V', 'H', 'G', 'W', 'C', 'G', 'G', 'G', 'N', 'D', 'Y', 'C', 'S',
              'G', 'G', 'N', 'C', 'Q', 'Y', 'R'],
             ['R', 'C', 'G', 'S', 'Q', 'A', 'G', 'G', 'A', 'T', 'C', 'T', 'N',
              'N', 'Q', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'F', 'G',
              'A', 'E', 'Y', 'C', 'G', 'A', 'G', 'C', 'Q', 'G', 'G', 'P', 'C',
              'R', 'A', '-', '-', '-', '-', 'D', 'I', 'K', 'C', 'G', 'S', 'Q',
              'A', 'G', 'G', 'K', 'L', 'C', 'P', 'N', 'N', 'L', 'C', 'C', 'S',
              'Q', 'W', 'G', 'F', 'C', 'G', 'L', 'G', 'S', 'E', 'F', 'C', 'G',
              'G', 'G', '-', 'C', 'Q', 'S', 'G']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 99.41)
        self.assertAlmostEqual(alignment.annotations["E-value"], 5.1e-17)
        self.assertAlmostEqual(alignment.annotations["Score"], 108.07)
        self.assertAlmostEqual(alignment.annotations["Identities"], 48)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.287)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 74.0)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[0:81],
            "ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGG",
        )
        self.assertEqual(
            alignment[1],
            "ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWT----SKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAG-CQGG",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "ercgeqgsnmecpnnlccsqygycgmggdycgkgcqngacwtskrcgsqaggatctnnqccsqygycgfgaeycgagcqgg                                                                                          ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "~~cg~~~~~~~C~~~~CCS~~G~CG~~~~~C~~~Cq~~c~~~~~~~~~Cg~~~~~~~c~~~~CCS~~G~CG~~~~~C~~~~Cq~~    ",
        )
        self.assertEqual(alignment.target.id, "1en2_A")
        self.assertEqual(
            alignment.target.seq[0:85],
            "ERCGSQGGGSTCPGLRCCSIWGWCGDSEPYCGRTCENKCWSGERSDHRCGAAVGNPPCGQDRCCSVHGWCGGGNDYCSGGNCQYR",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "1en2_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "UDA, agglutinin isolectin I/agglutinin isolectin V/ AG isolectin VI; hevein domain, superantigen, saccharide binding binding protein; HET: NAG; 1.40A {Urtica dioica} SCOP: g.3.1.1 g.3.1.1 PDB: 1eis_A* 1enm_A* 1ehd_A 1ehh_A* 1iqb_A",
        )
        self.assertEqual(
            alignment[0],
            "ERCGSQGGGSTCPGLRCCSIWGWCGDSEPYCGRTCENK-CWSGERSDHRCGAAVGNPPCGQDRCCSVHGWCGGGNDYCSGGNCQYR",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "CBCTTTTTSCCCGGGCEEETTSBEESSHHHHSTTEEESCGGGCCTTCBCSGGGTCCCCCTTCEEETTSBEESSHHHHSGGGEEEC    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "CCCCCCCCCccCCCCcccCCCceecccccccCCCCcCCCcccccCCcccCCcccccccCCCCeECCCceECCCccccCCCCcccC    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            "68999999999999999999999999999999999986764678998888889999999999999999999999988863         ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 0, 38, 38, 41, 45, 80, 81, 85],
                             [ 0, 38, 39, 42, 42, 77, 77, 81]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['E', 'R', 'C', 'G', 'S', 'Q', 'G', 'G', 'G', 'S', 'T', 'C', 'P',
              'G', 'L', 'R', 'C', 'C', 'S', 'I', 'W', 'G', 'W', 'C', 'G', 'D',
              'S', 'E', 'P', 'Y', 'C', 'G', 'R', 'T', 'C', 'E', 'N', 'K', '-',
              'C', 'W', 'S', 'G', 'E', 'R', 'S', 'D', 'H', 'R', 'C', 'G', 'A',
              'A', 'V', 'G', 'N', 'P', 'P', 'C', 'G', 'Q', 'D', 'R', 'C', 'C',
              'S', 'V', 'H', 'G', 'W', 'C', 'G', 'G', 'G', 'N', 'D', 'Y', 'C',
              'S', 'G', 'G', 'N', 'C', 'Q', 'Y', 'R'],
             ['E', 'R', 'C', 'G', 'E', 'Q', 'G', 'S', 'N', 'M', 'E', 'C', 'P',
              'N', 'N', 'L', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'M',
              'G', 'G', 'D', 'Y', 'C', 'G', 'K', 'G', 'C', 'Q', 'N', 'G', 'A',
              'C', 'W', 'T', '-', '-', '-', '-', 'S', 'K', 'R', 'C', 'G', 'S',
              'Q', 'A', 'G', 'G', 'A', 'T', 'C', 'T', 'N', 'N', 'Q', 'C', 'C',
              'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'F', 'G', 'A', 'E', 'Y', 'C',
              'G', 'A', 'G', '-', 'C', 'Q', 'G', 'G']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 99.38)
        self.assertAlmostEqual(alignment.annotations["E-value"], 1e-16)
        self.assertAlmostEqual(alignment.annotations["Score"], 104.25)
        self.assertAlmostEqual(alignment.annotations["Identities"], 53)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.284)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 72.3)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[87:167],
            "KCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYCGAGCQSG",
        )
        self.assertEqual(
            alignment[1],
            "KCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYCGAGCQSG",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                       kcgsqaggklcpnnlccsqwgfcglgsefcgggcqsgacstdkpcgkdaggrvctnnyccskwgscgigpgycgagcqsg    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "  ~cg~~~~~~~C~~g~CCs~~G~Cg~~~~~c~~~c~~~~~~g~Cg~~~~~~~c~~~~CCs~~g~CG~~~~~C~~~Cqs~  ",
        )
        self.assertEqual(alignment.target.id, "1uha_A")
        self.assertEqual(
            alignment.target.seq[2:80],
            "ECGERASGKRCPNGKCCSQWGYCGTTDNYCGQGCQSQCDYWRCGRDFGGRLCEEDMCCSKYGWCGYSDDHCEDGCQSQ",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "1uha_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Lectin-D2; chitin-binding domain, sugar binding protein; 1.50A {Phytolacca americana} SCOP: g.3.1.1 g.3.1.1 PDB: 1ulm_A* 1uln_A",
        )
        self.assertEqual(
            alignment[0],
            "ECGERASGKRCPNGKCCSQWGYCGTTDNYCGQGCQSQ--CDYWRCGRDFGGRLCEEDMCCSKYGWCGYSDDHCEDGCQSQ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "  CSBGGGTTBCCGGGCEECTTSCEESSHHHHSTTCCBCTTTTBCBGGGTTBCCSTTCEECTTSBEECSHHHHSTTCCBC  ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "  CCCCCCCCCcCCCCCccCCCccccCccccccCCccccccccccccccceecCCCCCccCCCccccCCccccccccccc  ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            "  588999999999999999999999999999999998624678999889999999999999999999999999999984  ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  2,  39,  39,  80],
                             [ 87, 124, 126, 167]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['E', 'C', 'G', 'E', 'R', 'A', 'S', 'G', 'K', 'R', 'C', 'P', 'N',
              'G', 'K', 'C', 'C', 'S', 'Q', 'W', 'G', 'Y', 'C', 'G', 'T', 'T',
              'D', 'N', 'Y', 'C', 'G', 'Q', 'G', 'C', 'Q', 'S', 'Q', '-', '-',
              'C', 'D', 'Y', 'W', 'R', 'C', 'G', 'R', 'D', 'F', 'G', 'G', 'R',
              'L', 'C', 'E', 'E', 'D', 'M', 'C', 'C', 'S', 'K', 'Y', 'G', 'W',
              'C', 'G', 'Y', 'S', 'D', 'D', 'H', 'C', 'E', 'D', 'G', 'C', 'Q',
              'S', 'Q'],
             ['K', 'C', 'G', 'S', 'Q', 'A', 'G', 'G', 'K', 'L', 'C', 'P', 'N',
              'N', 'L', 'C', 'C', 'S', 'Q', 'W', 'G', 'F', 'C', 'G', 'L', 'G',
              'S', 'E', 'F', 'C', 'G', 'G', 'G', 'C', 'Q', 'S', 'G', 'A', 'C',
              'S', 'T', 'D', 'K', 'P', 'C', 'G', 'K', 'D', 'A', 'G', 'G', 'R',
              'V', 'C', 'T', 'N', 'N', 'Y', 'C', 'C', 'S', 'K', 'W', 'G', 'S',
              'C', 'G', 'I', 'G', 'P', 'G', 'Y', 'C', 'G', 'A', 'G', 'C', 'Q',
              'S', 'G']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 98.20)
        self.assertAlmostEqual(alignment.annotations["E-value"], 1.6e-09)
        self.assertAlmostEqual(alignment.annotations["Score"], 66.41)
        self.assertAlmostEqual(alignment.annotations["Identities"], 59)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.390)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 36.7)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[85:127], "DIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACS"
        )
        self.assertEqual(alignment[1], "DIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACS")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                     dikcgsqaggklcpnnlccsqwgfcglgsefcgggcqsgacs                                            ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            " ~~~CG~~~~~~~C~~~~CCs~~G~CG~t~~~C~~gCq~~c~   ",
        )
        self.assertEqual(alignment.target.id, "4mpi_A")
        self.assertEqual(
            alignment.target.seq[1:42], "MEQCGRQAGGALCPGGLCCSQYGWCANTPEYCGSGCQSQCD"
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "4mpi_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Class I chitinase; hevein-like domain, chitin oligomers, sugar binding protein; HET: MES; 1.60A {Hevea brasiliensis subsp}",
        )
        self.assertEqual(alignment[0], "MEQCGRQAGGALCPGGLCCSQYGWCANTPEYCGSGCQSQ-CD")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            " CCBCBGGGTTBCCGGGCEECTTSBEECSHHHHSTTCCBCTT   ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            " ccccCCcCCCcccCCCCcCcccceecCCcccccccccccCC   ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            " 45688888889999999999999999999999999998554   ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  1,  40,  40,  42],
                             [ 85, 124, 125, 127]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['M', 'E', 'Q', 'C', 'G', 'R', 'Q', 'A', 'G', 'G', 'A', 'L', 'C',
              'P', 'G', 'G', 'L', 'C', 'C', 'S', 'Q', 'Y', 'G', 'W', 'C', 'A',
              'N', 'T', 'P', 'E', 'Y', 'C', 'G', 'S', 'G', 'C', 'Q', 'S', 'Q',
              '-', 'C', 'D'],
             ['D', 'I', 'K', 'C', 'G', 'S', 'Q', 'A', 'G', 'G', 'K', 'L', 'C',
              'P', 'N', 'N', 'L', 'C', 'C', 'S', 'Q', 'W', 'G', 'F', 'C', 'G',
              'L', 'G', 'S', 'E', 'F', 'C', 'G', 'G', 'G', 'C', 'Q', 'S', 'G',
              'A', 'C', 'S']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 98.12)
        self.assertAlmostEqual(alignment.annotations["E-value"], 3.1e-09)
        self.assertAlmostEqual(alignment.annotations["Score"], 64.73)
        self.assertAlmostEqual(alignment.annotations["Identities"], 68)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.600)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 34.0)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[87:124], "KCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSG"
        )
        self.assertEqual(alignment[1], "KCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGG--GCQSG")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                       kcgsqaggklcpnnlccsqwgfcglgsefcgggcqsg                                               ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            " ~CG~~~~~~~C~~~~CCS~~G~CG~t~~~C~~~~~Cq~~   ",
        )
        self.assertEqual(alignment.target.id, "1wkx_A")
        self.assertEqual(
            alignment.target.seq[1:40], "QCGRQAGGKLCPDNLCCSQWGWCGSTDEYCSPDHNCQSN"
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "1wkx_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Hevein isoform 2; allergen, lectin, agglutinin-toxin motif; 1.70A {Hevea brasiliensis} PDB: 1hev_A 1q9b_A* 4wp4_A",
        )
        self.assertEqual(alignment[0], "QCGRQAGGKLCPDNLCCSQWGWCGSTDEYCSPDHNCQSN")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            " CSBGGGTTBCCSTTCEECTTSCEESSHHHHCGGGTCCBS   ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            " CCCCcCCCcccCCCCeEeecCcccCCcccccCCCCccCC   ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            " 4788888899999999999999999999999689985     ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  1,  33,  35,  40],
                             [ 87, 119, 119, 124]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['Q', 'C', 'G', 'R', 'Q', 'A', 'G', 'G', 'K', 'L', 'C', 'P', 'D',
              'N', 'L', 'C', 'C', 'S', 'Q', 'W', 'G', 'W', 'C', 'G', 'S', 'T',
              'D', 'E', 'Y', 'C', 'S', 'P', 'D', 'H', 'N', 'C', 'Q', 'S', 'N'],
             ['K', 'C', 'G', 'S', 'Q', 'A', 'G', 'G', 'K', 'L', 'C', 'P', 'N',
              'N', 'L', 'C', 'C', 'S', 'Q', 'W', 'G', 'F', 'C', 'G', 'L', 'G',
              'S', 'E', 'F', 'C', 'G', 'G', '-', '-', 'G', 'C', 'Q', 'S', 'G']],
            dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 98.00)
        self.assertAlmostEqual(alignment.annotations["E-value"], 8e-09)
        self.assertAlmostEqual(alignment.annotations["Score"], 63.24)
        self.assertAlmostEqual(alignment.annotations["Identities"], 49)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.241)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 36.5)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[1:43], "RCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTS"
        )
        self.assertEqual(alignment[1], "RCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTS")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            " rcgeqgsnmecpnnlccsqygycgmggdycgkgcqngacwts                                                                                                                                ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "   ~CG~~~~~~~C~~~~CCs~~G~CG~t~~~C~~gCq~~c~~~ ",
        )
        self.assertEqual(alignment.target.id, "4mpi_A")
        self.assertEqual(
            alignment.target.seq[3:44], "QCGRQAGGALCPGGLCCSQYGWCANTPEYCGSGCQSQCDGG"
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "4mpi_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Class I chitinase; hevein-like domain, chitin oligomers, sugar binding protein; HET: MES; 1.60A {Hevea brasiliensis subsp}",
        )
        self.assertEqual(alignment[0], "QCGRQAGGALCPGGLCCSQYGWCANTPEYCGSGCQSQ-CDGG")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "   BCBGGGTTBCCGGGCEECTTSBEECSHHHHSTTCCBCTTCC ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "   ccCCcCCCcccCCCCcCcccceecCCcccccccccccCCCC ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            "   68888888999999999999999999999999999756543 ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 3, 40, 40, 44],
                             [ 1, 38, 39, 43]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['Q', 'C', 'G', 'R', 'Q', 'A', 'G', 'G', 'A', 'L', 'C', 'P', 'G',
              'G', 'L', 'C', 'C', 'S', 'Q', 'Y', 'G', 'W', 'C', 'A', 'N', 'T',
              'P', 'E', 'Y', 'C', 'G', 'S', 'G', 'C', 'Q', 'S', 'Q', '-', 'C',
              'D', 'G', 'G'],
             ['R', 'C', 'G', 'E', 'Q', 'G', 'S', 'N', 'M', 'E', 'C', 'P', 'N',
              'N', 'L', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'M', 'G',
              'G', 'D', 'Y', 'C', 'G', 'K', 'G', 'C', 'Q', 'N', 'G', 'A', 'C',
              'W', 'T', 'S']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 97.97)
        self.assertAlmostEqual(alignment.annotations["E-value"], 1e-08)
        self.assertAlmostEqual(alignment.annotations["Score"], 61.74)
        self.assertAlmostEqual(alignment.annotations["Identities"], 59)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.493)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 33.9)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[87:124], "KCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSG"
        )
        self.assertEqual(alignment[1], "KCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGG-CQSG")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                       kcgsqaggklcpnnlccsqwgfcglgsefcgggcqsg                                               ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "  ~CG~~~~~~~C~~~~CCS~~G~CG~t~~~C~~~~Cq~~    ",
        )
        self.assertEqual(alignment.target.id, "2lb7_A")
        self.assertEqual(
            alignment.target.seq[2:40], "RCGDQARGAKCPNCLCCGKYGFCGSGDAYCGAGSCQSQ"
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "2lb7_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "WAMP-1A, antimicrobial peptide 1A; antimicrobial protein; NMR {Triticum kiharae}",
        )
        self.assertEqual(alignment[0], "RCGDQARGAKCPNCLCCGKYGFCGSGDAYCGAGSCQSQ")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "  ECBGGGTTBCCCTTCEEETTTEEECSHHHHSTTSEEEC    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "  CCcCCCCCcccCCCCcCCcceeecCCccccCCCCccCC    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            "  5788888899999999999999999999999879976     ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  2,  35,  36,  40],
                             [ 87, 120, 120, 124]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['R', 'C', 'G', 'D', 'Q', 'A', 'R', 'G', 'A', 'K', 'C', 'P', 'N',
              'C', 'L', 'C', 'C', 'G', 'K', 'Y', 'G', 'F', 'C', 'G', 'S', 'G',
              'D', 'A', 'Y', 'C', 'G', 'A', 'G', 'S', 'C', 'Q', 'S', 'Q'],
             ['K', 'C', 'G', 'S', 'Q', 'A', 'G', 'G', 'K', 'L', 'C', 'P', 'N',
              'N', 'L', 'C', 'C', 'S', 'Q', 'W', 'G', 'F', 'C', 'G', 'L', 'G',
              'S', 'E', 'F', 'C', 'G', 'G', 'G', '-', 'C', 'Q', 'S', 'G']],
            dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 97.88)
        self.assertAlmostEqual(alignment.annotations["E-value"], 2e-08)
        self.assertAlmostEqual(alignment.annotations["Score"], 62.00)
        self.assertAlmostEqual(alignment.annotations["Identities"], 49)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.238)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 33.7)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[0:38], "ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNG"
        )
        self.assertEqual(alignment[1], "ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKG-CQNG")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "ercgeqgsnmecpnnlccsqygycgmggdycgkgcqng                                                                                                                                     ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "~~CG~~~~~~C~~~~CCS~~G~CG~t~~~C~~~~Cq~~   ",
        )
        self.assertEqual(alignment.target.id, "1p9g_A")
        self.assertEqual(
            alignment.target.seq[0:38], "ETCASRCPRPCNAGLCCSIYGYCGSGAAYCGAGNCRCQ"
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "1p9g_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "EAFP 2; antifungal peptide, atomic resolution, antifungal protein; HET: PCA; 0.84A {Eucommia ulmoides} SCOP: g.3.1.1 PDB: 1p9z_A*",
        )
        self.assertEqual(alignment[0], "ETCA-SRCPRPCNAGLCCSIYGYCGSGAAYCGAGNCRCQ")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "CCGGGGTTCCSCTTCEEETTSCEECSHHHHSTTTEEEC   ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "CCcCCcCCcccCCCCeECccceeCCCccccCCCccccC   ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            "6899578889999999999999999999999959976    ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 0,  4,  4, 33, 34, 38],
                             [ 0,  4,  5, 34, 34, 38]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['E', 'T', 'C', 'A', '-', 'S', 'R', 'C', 'P', 'R', 'P', 'C', 'N',
              'A', 'G', 'L', 'C', 'C', 'S', 'I', 'Y', 'G', 'Y', 'C', 'G', 'S',
              'G', 'A', 'A', 'Y', 'C', 'G', 'A', 'G', 'N', 'C', 'R', 'C', 'Q'],
             ['E', 'R', 'C', 'G', 'E', 'Q', 'G', 'S', 'N', 'M', 'E', 'C', 'P',
              'N', 'N', 'L', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'M',
              'G', 'G', 'D', 'Y', 'C', 'G', 'K', 'G', '-', 'C', 'Q', 'N', 'G']],
            dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 97.83)
        self.assertAlmostEqual(alignment.annotations["E-value"], 2.8e-08)
        self.assertAlmostEqual(alignment.annotations["Score"], 58.62)
        self.assertAlmostEqual(alignment.annotations["Identities"], 44)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.336)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 25.7)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(alignment.query.seq[91:118], "QAGGKLCPNNLCCSQWGFCGLGSEFCG")
        self.assertEqual(alignment[1], "QAGGKLCPNNLCCSQWGFCGLGSEFCG")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                           qaggklcpnnlccsqwgfcglgsefcg                                                     ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "  ~~~~~~C~~~~CCS~~G~CG~t~~~C~ ",
        )
        self.assertEqual(alignment.target.id, "1mmc_A")
        self.assertEqual(alignment.target.seq[2:29], "ECVRGRCPSGMCCSQFGYCGKGPKYCG")
        self.assertEqual(alignment.target.annotations["hmm_name"], "1mmc_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "AC-AMP2, antimicrobial peptide 2; antifungal antimicrobial, chitin-binding; NMR {Amaranthus caudatus} SCOP: g.3.1.2 PDB: 1zuv_A 1zwu_A* 1znt_A*",
        )
        self.assertEqual(alignment[0], "ECVRGRCPSGMCCSQFGYCGKGPKYCG")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "  CCSSSCCSTTCEECTTSCEESSHHHHC ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "  cCccCCCCCCCcccccceeCCchHhhC ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            "  788889999999999999999999996 ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  2,  29],
                             [ 91, 118]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['E', 'C', 'V', 'R', 'G', 'R', 'C', 'P', 'S', 'G', 'M', 'C', 'C',
              'S', 'Q', 'F', 'G', 'Y', 'C', 'G', 'K', 'G', 'P', 'K', 'Y', 'C',
              'G'],
             ['Q', 'A', 'G', 'G', 'K', 'L', 'C', 'P', 'N', 'N', 'L', 'C', 'C',
              'S', 'Q', 'W', 'G', 'F', 'C', 'G', 'L', 'G', 'S', 'E', 'F', 'C',
              'G']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 97.83)
        self.assertAlmostEqual(alignment.annotations["E-value"], 2.9e-08)
        self.assertAlmostEqual(alignment.annotations["Score"], 61.32)
        self.assertAlmostEqual(alignment.annotations["Identities"], 39)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.187)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 32.1)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[87:124], "KCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSG"
        )
        self.assertEqual(alignment[1], "KCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGG-CQSG")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                       kcgsqaggklcpnnlccsqwgfcglgsefcgggcqsg                                               ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            " ~CG~~~~~~C~~~~CCS~~G~CG~t~~~C~~~~Cq~~   ",
        )
        self.assertEqual(alignment.target.id, "1p9g_A")
        self.assertEqual(
            alignment.target.seq[1:38], "TCASRCPRPCNAGLCCSIYGYCGSGAAYCGAGNCRCQ"
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "1p9g_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "EAFP 2; antifungal peptide, atomic resolution, antifungal protein; HET: PCA; 0.84A {Eucommia ulmoides} SCOP: g.3.1.1 PDB: 1p9z_A*",
        )
        self.assertEqual(alignment[0], "TCA-SRCPRPCNAGLCCSIYGYCGSGAAYCGAGNCRCQ")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            " CGGGGTTCCSCTTCEEETTSCEECSHHHHSTTTEEEC   ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            " CcCCcCCcccCCCCeECccceeCCCccccCCCccccC   ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            " 477477889999999999999999999999859986    ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  1,   4,   4,  33,  34,  38],
                             [ 87,  90,  91, 120, 120, 124]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['T', 'C', 'A', '-', 'S', 'R', 'C', 'P', 'R', 'P', 'C', 'N', 'A',
              'G', 'L', 'C', 'C', 'S', 'I', 'Y', 'G', 'Y', 'C', 'G', 'S', 'G',
              'A', 'A', 'Y', 'C', 'G', 'A', 'G', 'N', 'C', 'R', 'C', 'Q'],
             ['K', 'C', 'G', 'S', 'Q', 'A', 'G', 'G', 'K', 'L', 'C', 'P', 'N',
              'N', 'L', 'C', 'C', 'S', 'Q', 'W', 'G', 'F', 'C', 'G', 'L', 'G',
              'S', 'E', 'F', 'C', 'G', 'G', 'G', '-', 'C', 'Q', 'S', 'G']],
            dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 97.76)
        self.assertAlmostEqual(alignment.annotations["E-value"], 4.7e-08)
        self.assertAlmostEqual(alignment.annotations["Score"], 59.52)
        self.assertAlmostEqual(alignment.annotations["Identities"], 50)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.302)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 34.7)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[0:38], "ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNG"
        )
        self.assertEqual(alignment[1], "ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGK--GCQNG")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "ercgeqgsnmecpnnlccsqygycgmggdycgkgcqng                                                                                                                                     ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "~~CG~~~~~~~C~~~~CCS~~G~CG~t~~~C~~~~~Cq~~   ",
        )
        self.assertEqual(alignment.target.id, "1wkx_A")
        self.assertEqual(
            alignment.target.seq[0:40], "EQCGRQAGGKLCPDNLCCSQWGWCGSTDEYCSPDHNCQSN"
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "1wkx_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Hevein isoform 2; allergen, lectin, agglutinin-toxin motif; 1.70A {Hevea brasiliensis} PDB: 1hev_A 1q9b_A* 4wp4_A",
        )
        self.assertEqual(alignment[0], "EQCGRQAGGKLCPDNLCCSQWGWCGSTDEYCSPDHNCQSN")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "CCSBGGGTTBCCSTTCEECTTSCEESSHHHHCGGGTCCBS   ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "CCCCCcCCCcccCCCCeEeecCcccCCcccccCCCCccCC   ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            "58999988899999999999999999999999778864     ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 0, 33, 35, 40],
                             [ 0, 33, 33, 38]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['E', 'Q', 'C', 'G', 'R', 'Q', 'A', 'G', 'G', 'K', 'L', 'C', 'P',
              'D', 'N', 'L', 'C', 'C', 'S', 'Q', 'W', 'G', 'W', 'C', 'G', 'S',
              'T', 'D', 'E', 'Y', 'C', 'S', 'P', 'D', 'H', 'N', 'C', 'Q', 'S',
              'N'],
             ['E', 'R', 'C', 'G', 'E', 'Q', 'G', 'S', 'N', 'M', 'E', 'C', 'P',
              'N', 'N', 'L', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'M',
              'G', 'G', 'D', 'Y', 'C', 'G', 'K', '-', '-', 'G', 'C', 'Q', 'N',
              'G']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 97.58)
        self.assertAlmostEqual(alignment.annotations["E-value"], 1.5e-07)
        self.assertAlmostEqual(alignment.annotations["Score"], 56.59)
        self.assertAlmostEqual(alignment.annotations["Identities"], 55)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.488)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 34.3)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[0:38], "ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNG"
        )
        self.assertEqual(alignment[1], "ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKG-CQNG")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "ercgeqgsnmecpnnlccsqygycgmggdycgkgcqng                                                                                                                                     ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            " ~~CG~~~~~~~C~~~~CCS~~G~CG~t~~~C~~~~Cq~~    ",
        )
        self.assertEqual(alignment.target.id, "2lb7_A")
        self.assertEqual(
            alignment.target.seq[1:40], "QRCGDQARGAKCPNCLCCGKYGFCGSGDAYCGAGSCQSQ"
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "2lb7_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "WAMP-1A, antimicrobial peptide 1A; antimicrobial protein; NMR {Triticum kiharae}",
        )
        self.assertEqual(alignment[0], "QRCGDQARGAKCPNCLCCGKYGFCGSGDAYCGAGSCQSQ")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            " EECBGGGTTBCCCTTCEEETTTEEECSHHHHSTTSEEEC    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            " CCCcCCCCCcccCCCCcCCcceeecCCccccCCCCccCC    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            " 47888888889999999999999999999999878865     ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 1, 35, 36, 40],
                             [ 0, 34, 34, 38]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['Q', 'R', 'C', 'G', 'D', 'Q', 'A', 'R', 'G', 'A', 'K', 'C', 'P',
              'N', 'C', 'L', 'C', 'C', 'G', 'K', 'Y', 'G', 'F', 'C', 'G', 'S',
              'G', 'D', 'A', 'Y', 'C', 'G', 'A', 'G', 'S', 'C', 'Q', 'S', 'Q'],
             ['E', 'R', 'C', 'G', 'E', 'Q', 'G', 'S', 'N', 'M', 'E', 'C', 'P',
              'N', 'N', 'L', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'M',
              'G', 'G', 'D', 'Y', 'C', 'G', 'K', 'G', '-', 'C', 'Q', 'N', 'G']],
            dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 97.58)
        self.assertAlmostEqual(alignment.annotations["E-value"], 1.6e-07)
        self.assertAlmostEqual(alignment.annotations["Score"], 55.11)
        self.assertAlmostEqual(alignment.annotations["Identities"], 42)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.250)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 23.2)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(alignment.query.seq[92:118], "AGGKLCPNNLCCSQWGFCGLGSEFCG")
        self.assertEqual(alignment[1], "AGGKLCPNNLCCSQWGFCGLGSEFCG")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                            aggklcpnnlccsqwgfcglgsefcg                                                     ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "    cg~~~C~~~~CCs~~G~CGtt~~~C~",
        )
        self.assertEqual(alignment.target.id, "2n1s_A")
        self.assertEqual(alignment.target.seq[4:30], "CYRGRCSGGLCCSKYGYCGSGPAYCG")
        self.assertEqual(alignment.target.annotations["hmm_name"], "2n1s_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "AMP-2; antimicrobial peptide, ICK, cystine knot inhibitor, cystine antimicrobial protein; NMR {Stellaria media}",
        )
        self.assertEqual(alignment[0], "CYRGRCSGGLCCSKYGYCGSGPAYCG")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "    CBTTBCSTTCEECTTSBEECSHHHHC",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "    cCCCCCCCCCccccccccCcchhhcC",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            "    34568999999999999999999985",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  4,  30],
                             [ 92, 118]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['C', 'Y', 'R', 'G', 'R', 'C', 'S', 'G', 'G', 'L', 'C', 'C', 'S',
              'K', 'Y', 'G', 'Y', 'C', 'G', 'S', 'G', 'P', 'A', 'Y', 'C', 'G'],
             ['A', 'G', 'G', 'K', 'L', 'C', 'P', 'N', 'N', 'L', 'C', 'C', 'S',
              'Q', 'W', 'G', 'F', 'C', 'G', 'L', 'G', 'S', 'E', 'F', 'C', 'G']],
            dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 97.56)
        self.assertAlmostEqual(alignment.annotations["E-value"], 1.7e-07)
        self.assertAlmostEqual(alignment.annotations["Score"], 55.43)
        self.assertAlmostEqual(alignment.annotations["Identities"], 42)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.160)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 29.5)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[83:118], "RADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCG"
        )
        self.assertEqual(alignment[1], "RADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCG")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                   radikcgsqaggklcpnnlccsqwgfcglgsefcg                                                     ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            " ~~~~~CG~~~g~C~~g~CCS~~G~CG~~~~~C~ ",
        )
        self.assertEqual(alignment.target.id, "2kus_A")
        self.assertEqual(
            alignment.target.seq[1:34], "GPNGQCGPGWGGCRGGLCCSQYGYCGSGPKYCA"
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "2kus_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "SM-AMP-1.1A; plant antimicrobial peptide, chitin-binding peptide, antimic protein; NMR {Stellaria media}",
        )
        self.assertEqual(alignment[0], "GPNGQCGPGWG--GCRGGLCCSQYGYCGSGPKYCA")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            " CTTCBCBTTTBCCCTTCEECTTSBEECSHHHHC ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            " CCCcccCCCCCcCCCCcEECCCceecCChhhhC ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            " 467889988876999999999999999999986 ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  1,  12,  12,  34],
                             [ 83,  94,  96, 118]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['G', 'P', 'N', 'G', 'Q', 'C', 'G', 'P', 'G', 'W', 'G', '-', '-',
              'G', 'C', 'R', 'G', 'G', 'L', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y',
              'C', 'G', 'S', 'G', 'P', 'K', 'Y', 'C', 'A'],
             ['R', 'A', 'D', 'I', 'K', 'C', 'G', 'S', 'Q', 'A', 'G', 'G', 'K',
              'L', 'C', 'P', 'N', 'N', 'L', 'C', 'C', 'S', 'Q', 'W', 'G', 'F',
              'C', 'G', 'L', 'G', 'S', 'E', 'F', 'C', 'G']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 97.24)
        self.assertAlmostEqual(alignment.annotations["E-value"], 1.1e-06)
        self.assertAlmostEqual(alignment.annotations["Score"], 69.81)
        self.assertAlmostEqual(alignment.annotations["Identities"], 68)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.622)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 35.1)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[87:124], "KCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSG"
        )
        self.assertEqual(alignment[1], "KCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSG")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                       kcgsqaggklcpnnlccsqwgfcglgsefcgggcqsg                                               ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "  ~cG~~~~~~~c~~~~ccs~~g~cg~~~~~C~~~cq~~                                                                                                                                                                                                                                                                              ",
        )
        self.assertEqual(alignment.target.id, "2dkv_A")
        self.assertEqual(
            alignment.target.seq[2:39], "QCGAQAGGARCPNCLCCSRWGWCGTTSDFCGDGCQSQ"
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "2dkv_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Chitinase; whole structure, oryza sativa L. japonica, hydrolase; HET: MES; 2.00A {Oryza sativa japonica group} PDB: 3iwr_A*",
        )
        self.assertEqual(alignment[0], "QCGAQAGGARCPNCLCCSRWGWCGTTSDFCGDGCQSQ")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "  BCSTTTTTCCCGGGCEECTTSBEESSHHHHSTTCCBC                                                                                                                                                                                                                                                                              ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "  CCCCCCCCCcCCCCCeeCcCCcccCCccccCccccCC                                                                                                                                                                                                                                                                              ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            "  4888889999999999999999999999999999987                                                                                                                                                                                                                                                                              ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  2,  39],
                             [ 87, 124]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['Q', 'C', 'G', 'A', 'Q', 'A', 'G', 'G', 'A', 'R', 'C', 'P', 'N',
              'C', 'L', 'C', 'C', 'S', 'R', 'W', 'G', 'W', 'C', 'G', 'T', 'T',
              'S', 'D', 'F', 'C', 'G', 'D', 'G', 'C', 'Q', 'S', 'Q'],
             ['K', 'C', 'G', 'S', 'Q', 'A', 'G', 'G', 'K', 'L', 'C', 'P', 'N',
              'N', 'L', 'C', 'C', 'S', 'Q', 'W', 'G', 'F', 'C', 'G', 'L', 'G',
              'S', 'E', 'F', 'C', 'G', 'G', 'G', 'C', 'Q', 'S', 'G']],
            dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 97.14)
        self.assertAlmostEqual(alignment.annotations["E-value"], 1.9e-06)
        self.assertAlmostEqual(alignment.annotations["Score"], 50.76)
        self.assertAlmostEqual(alignment.annotations["Identities"], 50)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.394)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 24.5)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(alignment.query.seq[5:33], "QGSNMECPNNLCCSQYGYCGMGGDYCGK")
        self.assertEqual(alignment[1], "QGSNMECPNNLCCSQYGYCGMGGDYCGK")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "     qgsnmecpnnlccsqygycgmggdycgk                                                                                                                                          ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "  ~~~~~~C~~~~CCS~~G~CG~t~~~C~~",
        )
        self.assertEqual(alignment.target.id, "1mmc_A")
        self.assertEqual(alignment.target.seq[2:30], "ECVRGRCPSGMCCSQFGYCGKGPKYCGR")
        self.assertEqual(alignment.target.annotations["hmm_name"], "1mmc_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "AC-AMP2, antimicrobial peptide 2; antifungal antimicrobial, chitin-binding; NMR {Amaranthus caudatus} SCOP: g.3.1.2 PDB: 1zuv_A 1zwu_A* 1znt_A*",
        )
        self.assertEqual(alignment[0], "ECVRGRCPSGMCCSQFGYCGKGPKYCGR")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "  CCSSSCCSTTCEECTTSCEESSHHHHCC",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "  cCccCCCCCCCcccccceeCCchHhhCc",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            "  4456689999999999999999999974",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 2, 30],
                             [ 5, 33]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['E', 'C', 'V', 'R', 'G', 'R', 'C', 'P', 'S', 'G', 'M', 'C', 'C',
              'S', 'Q', 'F', 'G', 'Y', 'C', 'G', 'K', 'G', 'P', 'K', 'Y', 'C',
              'G', 'R'],
             ['Q', 'G', 'S', 'N', 'M', 'E', 'C', 'P', 'N', 'N', 'L', 'C', 'C',
              'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'M', 'G', 'G', 'D', 'Y', 'C',
              'G', 'K']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 97.02)
        self.assertAlmostEqual(alignment.annotations["E-value"], 3.5e-06)
        self.assertAlmostEqual(alignment.annotations["Score"], 49.35)
        self.assertAlmostEqual(alignment.annotations["Identities"], 56)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.468)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 23.7)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(alignment.query.seq[1:32], "RCGEQGSNMECPNNLCCSQYGYCGMGGDYCG")
        self.assertEqual(alignment[1], "RCGEQGSNMECPNNLCCSQYGYCGMGGDYCG")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            " rcgeqgsnmecpnnlccsqygycgmggdycg                                                                                                                                           ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "   ~cg~~~C~~~~CCs~~G~CGtt~~~C~",
        )
        self.assertEqual(alignment.target.id, "2n1s_A")
        self.assertEqual(alignment.target.seq[3:30], "QCYRGRCSGGLCCSKYGYCGSGPAYCG")
        self.assertEqual(alignment.target.annotations["hmm_name"], "2n1s_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "AMP-2; antimicrobial peptide, ICK, cystine knot inhibitor, cystine antimicrobial protein; NMR {Stellaria media}",
        )
        self.assertEqual(alignment[0], "QCY----RGRCSGGLCCSKYGYCGSGPAYCG")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "   BCBTTBCSTTCEECTTSBEECSHHHHC",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "   hcCCCCCCCCCccccccccCcchhhcC",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            "   566348999999999999999999985",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 3,  6,  6, 30],
                             [ 1,  4,  8, 32]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['Q', 'C', 'Y', '-', '-', '-', '-', 'R', 'G', 'R', 'C', 'S', 'G',
              'G', 'L', 'C', 'C', 'S', 'K', 'Y', 'G', 'Y', 'C', 'G', 'S', 'G',
              'P', 'A', 'Y', 'C', 'G'],
             ['R', 'C', 'G', 'E', 'Q', 'G', 'S', 'N', 'M', 'E', 'C', 'P', 'N',
              'N', 'L', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'M', 'G',
              'G', 'D', 'Y', 'C', 'G']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 96.79)
        self.assertAlmostEqual(alignment.annotations["E-value"], 9.9e-06)
        self.assertAlmostEqual(alignment.annotations["Score"], 64.37)
        self.assertAlmostEqual(alignment.annotations["Identities"], 53)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.380)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 35.8)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[43:81], "KRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGG"
        )
        self.assertEqual(alignment[1], "KRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGG")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                           krcgsqaggatctnnqccsqygycgfgaeycgagcqgg                                                                                          ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            " ~~cG~~~~~~~c~~~~ccs~~g~cg~~~~~C~~~cq~~                                                                                                                                                                                                                                                                              ",
        )
        self.assertEqual(alignment.target.id, "2dkv_A")
        self.assertEqual(
            alignment.target.seq[1:39], "EQCGAQAGGARCPNCLCCSRWGWCGTTSDFCGDGCQSQ"
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "2dkv_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Chitinase; whole structure, oryza sativa L. japonica, hydrolase; HET: MES; 2.00A {Oryza sativa japonica group} PDB: 3iwr_A*",
        )
        self.assertEqual(alignment[0], "EQCGAQAGGARCPNCLCCSRWGWCGTTSDFCGDGCQSQ")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            " CBCSTTTTTCCCGGGCEECTTSBEESSHHHHSTTCCBC                                                                                                                                                                                                                                                                              ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            " CCCCCCCCCCcCCCCCeeCcCCcccCCccccCccccCC                                                                                                                                                                                                                                                                              ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            " 47999999999999999999999999999999999965                                                                                                                                                                                                                                                                              ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 1, 39],
                             [43, 81]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['E', 'Q', 'C', 'G', 'A', 'Q', 'A', 'G', 'G', 'A', 'R', 'C', 'P',
              'N', 'C', 'L', 'C', 'C', 'S', 'R', 'W', 'G', 'W', 'C', 'G', 'T',
              'T', 'S', 'D', 'F', 'C', 'G', 'D', 'G', 'C', 'Q', 'S', 'Q'],
             ['K', 'R', 'C', 'G', 'S', 'Q', 'A', 'G', 'G', 'A', 'T', 'C', 'T',
              'N', 'N', 'Q', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'F',
              'G', 'A', 'E', 'Y', 'C', 'G', 'A', 'G', 'C', 'Q', 'G', 'G']],
            dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 96.68)
        self.assertAlmostEqual(alignment.annotations["E-value"], 1.6e-05)
        self.assertAlmostEqual(alignment.annotations["Score"], 46.90)
        self.assertAlmostEqual(alignment.annotations["Identities"], 50)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.287)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 28.1)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[42:76], "SKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGA"
        )
        self.assertEqual(alignment[1], "SKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGA")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                          skrcgsqaggatctnnqccsqygycgfgaeycga                                                                                               ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "   ~~~CG~~~g~C~~g~CCS~~G~CG~~~~~C~~",
        )
        self.assertEqual(alignment.target.id, "2kus_A")
        self.assertEqual(alignment.target.seq[3:35], "NGQCGPGWGGCRGGLCCSQYGYCGSGPKYCAH")
        self.assertEqual(alignment.target.annotations["hmm_name"], "2kus_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "SM-AMP-1.1A; plant antimicrobial peptide, chitin-binding peptide, antimic protein; NMR {Stellaria media}",
        )
        self.assertEqual(alignment[0], "NGQCGPGWG--GCRGGLCCSQYGYCGSGPKYCAH")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "   TCBCBTTTBCCCTTCEECTTSBEECSHHHHCC",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "   CcccCCCCCcCCCCcEECCCceecCChhhhCc",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            "   45788888779999999999999999999963",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 3, 12, 12, 35],
                             [42, 51, 53, 76]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['N', 'G', 'Q', 'C', 'G', 'P', 'G', 'W', 'G', '-', '-', 'G', 'C',
              'R', 'G', 'G', 'L', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y', 'C', 'G',
              'S', 'G', 'P', 'K', 'Y', 'C', 'A', 'H'],
             ['S', 'K', 'R', 'C', 'G', 'S', 'Q', 'A', 'G', 'G', 'A', 'T', 'C',
              'T', 'N', 'N', 'Q', 'C', 'C', 'S', 'Q', 'Y', 'G', 'Y', 'C', 'G',
              'F', 'G', 'A', 'E', 'Y', 'C', 'G', 'A']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 92.12)
        self.assertAlmostEqual(alignment.annotations["E-value"], 0.024)
        self.assertAlmostEqual(alignment.annotations["Score"], 45.43)
        self.assertAlmostEqual(alignment.annotations["Identities"], 20)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.652)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 136.2)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[8:169],
            "NMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYCGAGCQSGGC",
        )
        self.assertEqual(
            alignment[1],
            "NMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYCGAGCQSGGC",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "        nmecpnnlccsqygycgmggdycgkgcqngacwtskrcgsqaggatctnnqccsqygycgfgaeycgagcqggpcradikcgsqaggklcpnnlccsqwgfcglgsefcgggcqsgacstdkpcgkdaggrvctnnyccskwgscgigpgycgagcqsggc  ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "        ~~~c~~~~cc~~~~~c~~~~~~c~~~c~~~~c~~~~~c~~~~~~c~~~~cc~~~~~c~~~~~~c~~~c~~~~c~~~~~c~~~~~~c~~~~cc~~~~~c~~~~~~c~~~c~~~~c~~~~~c~~~~~~c~~~~cc~~~~~c~~~~~~c~~~c~~~~c ",
        )
        self.assertEqual(alignment.target.id, "1wga")
        self.assertEqual(
            alignment.target.seq[8:163],
            "XXXCXXXXCCXXXXXCXXXXXXCXXXCXXXXCXXXXXCXXXXXXCXXXXCCXXXXXCXXXXXXCXXXCXXXXCXXXXXCXXXXXXCXXXXCCXXXXXCXXXXXXCXXXCXXXXCXXXXXCXXXXXXCXXXXCCXXXXXCXXXXXXCXXXCXXXXC",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "1wga")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "; lectin (agglutinin); NMR {}",
        )
        self.assertEqual(
            alignment[0],
            "XXXCXXXXCCXXXXXCXXXXXXCXXXCXXXXCXXXXXCXXX--XXXCXXXXCCXXXXXCXXXXXXCXXXCXXXXCXXXXXCXXX--XXXCXXXXCCXXXXXCXXXXXXCXXXCXXXXCXXXXXCXXX--XXXCXXXXCCXXXXXCXXXXXXCXXXCXXXXC",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "        ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            "        45688889999999999999999999999999999988754467999999999999999999999999999999999887544679999999999999999999999999999999998876446799999999999998888888877765444 ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  8,  49,  49,  90,  90, 131, 131, 163],
                             [  8,  49,  51,  92,  94, 135, 137, 169]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['X', 'X', 'X', 'C', 'X', 'X', 'X', 'X', 'C', 'C', 'X', 'X', 'X',
              'X', 'X', 'C', 'X', 'X', 'X', 'X', 'X', 'X', 'C', 'X', 'X', 'X',
              'C', 'X', 'X', 'X', 'X', 'C', 'X', 'X', 'X', 'X', 'X', 'C', 'X',
              'X', 'X', '-', '-', 'X', 'X', 'X', 'C', 'X', 'X', 'X', 'X', 'C',
              'C', 'X', 'X', 'X', 'X', 'X', 'C', 'X', 'X', 'X', 'X', 'X', 'X',
              'C', 'X', 'X', 'X', 'C', 'X', 'X', 'X', 'X', 'C', 'X', 'X', 'X',
              'X', 'X', 'C', 'X', 'X', 'X', '-', '-', 'X', 'X', 'X', 'C', 'X',
              'X', 'X', 'X', 'C', 'C', 'X', 'X', 'X', 'X', 'X', 'C', 'X', 'X',
              'X', 'X', 'X', 'X', 'C', 'X', 'X', 'X', 'C', 'X', 'X', 'X', 'X',
              'C', 'X', 'X', 'X', 'X', 'X', 'C', 'X', 'X', 'X', '-', '-', 'X',
              'X', 'X', 'C', 'X', 'X', 'X', 'X', 'C', 'C', 'X', 'X', 'X', 'X',
              'X', 'C', 'X', 'X', 'X', 'X', 'X', 'X', 'C', 'X', 'X', 'X', 'C',
              'X', 'X', 'X', 'X', 'C'],
             ['N', 'M', 'E', 'C', 'P', 'N', 'N', 'L', 'C', 'C', 'S', 'Q', 'Y',
              'G', 'Y', 'C', 'G', 'M', 'G', 'G', 'D', 'Y', 'C', 'G', 'K', 'G',
              'C', 'Q', 'N', 'G', 'A', 'C', 'W', 'T', 'S', 'K', 'R', 'C', 'G',
              'S', 'Q', 'A', 'G', 'G', 'A', 'T', 'C', 'T', 'N', 'N', 'Q', 'C',
              'C', 'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'F', 'G', 'A', 'E', 'Y',
              'C', 'G', 'A', 'G', 'C', 'Q', 'G', 'G', 'P', 'C', 'R', 'A', 'D',
              'I', 'K', 'C', 'G', 'S', 'Q', 'A', 'G', 'G', 'K', 'L', 'C', 'P',
              'N', 'N', 'L', 'C', 'C', 'S', 'Q', 'W', 'G', 'F', 'C', 'G', 'L',
              'G', 'S', 'E', 'F', 'C', 'G', 'G', 'G', 'C', 'Q', 'S', 'G', 'A',
              'C', 'S', 'T', 'D', 'K', 'P', 'C', 'G', 'K', 'D', 'A', 'G', 'G',
              'R', 'V', 'C', 'T', 'N', 'N', 'Y', 'C', 'C', 'S', 'K', 'W', 'G',
              'S', 'C', 'G', 'I', 'G', 'P', 'G', 'Y', 'C', 'G', 'A', 'G', 'C',
              'Q', 'S', 'G', 'G', 'C']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 88.54)
        self.assertAlmostEqual(alignment.annotations["E-value"], 0.1)
        self.assertAlmostEqual(alignment.annotations["Score"], 39.32)
        self.assertAlmostEqual(alignment.annotations["Identities"], 37)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.037)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 27.7)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[83:118], "RADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCG"
        )
        self.assertEqual(alignment[1], "RADIKCGSQA-----GGKLCPN---NLCCSQWGFCGLGSEFCG")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                   radikcgsqaggklcpnnlccsqwgfcglgsefcg                                                     ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "     r~d~rCg~~~~~~~~~~~~C~~~~~~~CCs~~~~cg~~~~~c~                                                                                                                                                                                            ",
        )
        self.assertEqual(alignment.target.id, "4z8i_A")
        self.assertEqual(
            alignment.target.seq[5:48], "RSDGRCGPNYPAPDANPGECNPHAVDHCCSEWGWCGRETSHCT"
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "4z8i_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "BBTPGRP3, peptidoglycan recognition protein 3; chitin-binding domain, AM hydrolase; 2.70A {Branchiostoma belcheri tsingtauense}",
        )
        self.assertEqual(alignment[0], "RSDGRCGPNYPAPDANPGECNPHAVDHCCSEWGWCGRETSHCT")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "     CSSSBCBSSSCBTTBSSBBCCTTSSCCEECTTSBEECSHHHHH                                                                                                                                                                                            ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "     CCCCCCCCCCCCCCCCCcccCCCCCCCCCCCCCEEeCCccccc                                                                                                                                                                                            ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            "     55666776655567775899999999999999886                                                                                                                                                                                                    ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  5,  15,  20,  27,  30,  48],
                             [ 83,  93,  93, 100, 100, 118]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['R', 'S', 'D', 'G', 'R', 'C', 'G', 'P', 'N', 'Y', 'P', 'A', 'P',
              'D', 'A', 'N', 'P', 'G', 'E', 'C', 'N', 'P', 'H', 'A', 'V', 'D',
              'H', 'C', 'C', 'S', 'E', 'W', 'G', 'W', 'C', 'G', 'R', 'E', 'T',
              'S', 'H', 'C', 'T'],
             ['R', 'A', 'D', 'I', 'K', 'C', 'G', 'S', 'Q', 'A', '-', '-', '-',
              '-', '-', 'G', 'G', 'K', 'L', 'C', 'P', 'N', '-', '-', '-', 'N',
              'L', 'C', 'C', 'S', 'Q', 'W', 'G', 'F', 'C', 'G', 'L', 'G', 'S',
              'E', 'F', 'C', 'G']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 88.13)
        self.assertAlmostEqual(alignment.annotations["E-value"], 0.12)
        self.assertAlmostEqual(alignment.annotations["Score"], 39.05)
        self.assertAlmostEqual(alignment.annotations["Identities"], 33)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.005)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 26.8)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[42:75], "SKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCG"
        )
        self.assertEqual(alignment[1], "SKRCGSQA-----GGATCTN---NQCCSQYGYCGFGAEYCG")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                          skrcgsqaggatctnnqccsqygycgfgaeycg                                                                                                ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "       d~rCg~~~~~~~~~~~~C~~~~~~~CCs~~~~cg~~~~~c~                                                                                                                                                                                            ",
        )
        self.assertEqual(alignment.target.id, "4z8i_A")
        self.assertEqual(
            alignment.target.seq[7:48], "DGRCGPNYPAPDANPGECNPHAVDHCCSEWGWCGRETSHCT"
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "4z8i_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "BBTPGRP3, peptidoglycan recognition protein 3; chitin-binding domain, AM hydrolase; 2.70A {Branchiostoma belcheri tsingtauense}",
        )
        self.assertEqual(alignment[0], "DGRCGPNYPAPDANPGECNPHAVDHCCSEWGWCGRETSHCT")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "       SSBCBSSSCBTTBSSBBCCTTSSCCEECTTSBEECSHHHHH                                                                                                                                                                                            ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "       CCCCCCCCCCCCCCCcccCCCCCCCCCCCCCEEeCCccccc                                                                                                                                                                                            ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            "       346666656678875899999999999999997                                                                                                                                                                                                    ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 7, 15, 20, 27, 30, 48],
                             [42, 50, 50, 57, 57, 75]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['D', 'G', 'R', 'C', 'G', 'P', 'N', 'Y', 'P', 'A', 'P', 'D', 'A',
              'N', 'P', 'G', 'E', 'C', 'N', 'P', 'H', 'A', 'V', 'D', 'H', 'C',
              'C', 'S', 'E', 'W', 'G', 'W', 'C', 'G', 'R', 'E', 'T', 'S', 'H',
              'C', 'T'],
             ['S', 'K', 'R', 'C', 'G', 'S', 'Q', 'A', '-', '-', '-', '-', '-',
              'G', 'G', 'A', 'T', 'C', 'T', 'N', '-', '-', '-', 'N', 'Q', 'C',
              'C', 'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'F', 'G', 'A', 'E', 'Y',
              'C', 'G']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 85.66)
        self.assertAlmostEqual(alignment.annotations["E-value"], 0.21)
        self.assertAlmostEqual(alignment.annotations["Score"], 38.06)
        self.assertAlmostEqual(alignment.annotations["Identities"], 37)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.037)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[83:118], "RADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCG"
        )
        self.assertEqual(alignment[1], "RADIKCGSQA-----GGKLCPN---NLCCSQWGFCGLGSEFCG")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                   radikcgsqaggklcpnnlccsqwgfcglgsefcg                                                     ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                      ~~d~~Cg~~~~~~~~~~~~C~~~~~~~CCs~~~~Cg~~~~~C~                                                                                                                                                                                               ",
        )
        self.assertEqual(alignment.target.id, "4zxm_A")
        self.assertEqual(
            alignment.target.seq[22:65], "RSDGRCGPNYPAPDANPGECNPHAVDHCCSEWGWCGRETSHCT"
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "4zxm_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "PGRP domain of peptidoglycan recognition protein; amidase, hydrolase; 2.80A {Branchiostoma belcheri tsingtauense}",
        )
        self.assertEqual(alignment[0], "RSDGRCGPNYPAPDANPGECNPHAVDHCCSEWGWCGRETSHCT")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "                                                                                                                                                                                                                                                                ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                      CCCCCCCCCCCCCCCCCcccCCCCCCCCCCCCCeEeCCCCCcC                                                                                                                                                                                               ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            "                      34555555544457776899999999999999875                                                                                                                                                                                                       ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 22,  32,  37,  44,  47,  65],
                             [ 83,  93,  93, 100, 100, 118]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['R', 'S', 'D', 'G', 'R', 'C', 'G', 'P', 'N', 'Y', 'P', 'A', 'P',
              'D', 'A', 'N', 'P', 'G', 'E', 'C', 'N', 'P', 'H', 'A', 'V', 'D',
              'H', 'C', 'C', 'S', 'E', 'W', 'G', 'W', 'C', 'G', 'R', 'E', 'T',
              'S', 'H', 'C', 'T'],
             ['R', 'A', 'D', 'I', 'K', 'C', 'G', 'S', 'Q', 'A', '-', '-', '-',
              '-', '-', 'G', 'G', 'K', 'L', 'C', 'P', 'N', '-', '-', '-', 'N',
              'L', 'C', 'C', 'S', 'Q', 'W', 'G', 'F', 'C', 'G', 'L', 'G', 'S',
              'E', 'F', 'C', 'G']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 85.31)
        self.assertAlmostEqual(alignment.annotations["E-value"], 0.23)
        self.assertAlmostEqual(alignment.annotations["Score"], 37.89)
        self.assertAlmostEqual(alignment.annotations["Identities"], 33)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.023)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[43:79], "KRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQ"
        )
        self.assertEqual(alignment[1], "KRCGSQA-----GGATCTN---NQCCSQYGYCGFGAEYCG-AGCQ")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                           krcgsqaggatctnnqccsqygycgfgaeycgagcq                                                                                            ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                         ~~Cg~~~~~~~~~~~~C~~~~~~~CCs~~~~Cg~~~~~C~~~~c~                                                                                                                                                                                          ",
        )
        self.assertEqual(alignment.target.id, "4zxm_A")
        self.assertEqual(
            alignment.target.seq[25:70], "GRCGPNYPAPDANPGECNPHAVDHCCSEWGWCGRETSHCTCSSCV"
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "4zxm_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "PGRP domain of peptidoglycan recognition protein; amidase, hydrolase; 2.80A {Branchiostoma belcheri tsingtauense}",
        )
        self.assertEqual(alignment[0], "GRCGPNYPAPDANPGECNPHAVDHCCSEWGWCGRETSHCTCSSCV")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "                                                                                                                                                                                                                                                                ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                         CCCCCCCCCCCCCCcccCCCCCCCCCCCCCeEeCCCCCcCCcccc                                                                                                                                                                                          ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            "                         456665555788768999999999999999864555                                                                                                                                                                                                   ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[25, 32, 37, 44, 47, 65, 66, 70],
                             [43, 50, 50, 57, 57, 75, 75, 79]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['G', 'R', 'C', 'G', 'P', 'N', 'Y', 'P', 'A', 'P', 'D', 'A', 'N',
              'P', 'G', 'E', 'C', 'N', 'P', 'H', 'A', 'V', 'D', 'H', 'C', 'C',
              'S', 'E', 'W', 'G', 'W', 'C', 'G', 'R', 'E', 'T', 'S', 'H', 'C',
              'T', 'C', 'S', 'S', 'C', 'V'],
             ['K', 'R', 'C', 'G', 'S', 'Q', 'A', '-', '-', '-', '-', '-', 'G',
              'G', 'A', 'T', 'C', 'T', 'N', '-', '-', '-', 'N', 'Q', 'C', 'C',
              'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'F', 'G', 'A', 'E', 'Y', 'C',
              'G', '-', 'A', 'G', 'C', 'Q']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 72.00)
        self.assertAlmostEqual(alignment.annotations["E-value"], 1.3)
        self.assertAlmostEqual(alignment.annotations["Score"], 35.47)
        self.assertAlmostEqual(alignment.annotations["Identities"], 20)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.657)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 98.0)
        self.assertEqual(alignment.query.id, "2UVO:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[8:126],
            "NMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGAC",
        )
        self.assertEqual(
            alignment[1],
            "NMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGAC",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "        nmecpnnlccsqygycgmggdycgkgcqngacwtskrcgsqaggatctnnqccsqygycgfgaeycgagcqggpcradikcgsqaggklcpnnlccsqwgfcglgsefcgggcqsgac                                             ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                                                 ~~~c~~~~cc~~~~~c~~~~~~c~~~c~~~~c~~~~~c~~~~~~c~~~~cc~~~~~c~~~~~~c~~~c~~~~c~~~~~c~~~~~~c~~~~cc~~~~~c~~~~~~c~~~c~~~~c ",
        )
        self.assertEqual(alignment.target.id, "1wga")
        self.assertEqual(
            alignment.target.seq[49:163],
            "XXXCXXXXCCXXXXXCXXXXXXCXXXCXXXXCXXXXXCXXXXXXCXXXXCCXXXXXCXXXXXXCXXXCXXXXCXXXXXCXXXXXXCXXXXCCXXXXXCXXXXXXCXXXCXXXXC",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "1wga")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "; lectin (agglutinin); NMR {}",
        )
        self.assertEqual(
            alignment[0],
            "XXXCXXXXCCXXXXXCXXXXXXCXXXCXXXXCXXXXXCXXX--XXXCXXXXCCXXXXXCXXXXXXCXXXCXXXXCXXXXXCXXX--XXXCXXXXCCXXXXXCXXXXXXCXXXCXXXXC",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                                                 cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Confidence"],
            "                                                 456788899999999999999999999999999998887544679999999999999999999999999999999998875446799999999999998888888776655443 ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 49,  90,  90, 131, 131, 163],
                             [  8,  49,  51,  92,  94, 126]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['X', 'X', 'X', 'C', 'X', 'X', 'X', 'X', 'C', 'C', 'X', 'X', 'X',
              'X', 'X', 'C', 'X', 'X', 'X', 'X', 'X', 'X', 'C', 'X', 'X', 'X',
              'C', 'X', 'X', 'X', 'X', 'C', 'X', 'X', 'X', 'X', 'X', 'C', 'X',
              'X', 'X', '-', '-', 'X', 'X', 'X', 'C', 'X', 'X', 'X', 'X', 'C',
              'C', 'X', 'X', 'X', 'X', 'X', 'C', 'X', 'X', 'X', 'X', 'X', 'X',
              'C', 'X', 'X', 'X', 'C', 'X', 'X', 'X', 'X', 'C', 'X', 'X', 'X',
              'X', 'X', 'C', 'X', 'X', 'X', '-', '-', 'X', 'X', 'X', 'C', 'X',
              'X', 'X', 'X', 'C', 'C', 'X', 'X', 'X', 'X', 'X', 'C', 'X', 'X',
              'X', 'X', 'X', 'X', 'C', 'X', 'X', 'X', 'C', 'X', 'X', 'X', 'X',
              'C'],
             ['N', 'M', 'E', 'C', 'P', 'N', 'N', 'L', 'C', 'C', 'S', 'Q', 'Y',
              'G', 'Y', 'C', 'G', 'M', 'G', 'G', 'D', 'Y', 'C', 'G', 'K', 'G',
              'C', 'Q', 'N', 'G', 'A', 'C', 'W', 'T', 'S', 'K', 'R', 'C', 'G',
              'S', 'Q', 'A', 'G', 'G', 'A', 'T', 'C', 'T', 'N', 'N', 'Q', 'C',
              'C', 'S', 'Q', 'Y', 'G', 'Y', 'C', 'G', 'F', 'G', 'A', 'E', 'Y',
              'C', 'G', 'A', 'G', 'C', 'Q', 'G', 'G', 'P', 'C', 'R', 'A', 'D',
              'I', 'K', 'C', 'G', 'S', 'Q', 'A', 'G', 'G', 'K', 'L', 'C', 'P',
              'N', 'N', 'L', 'C', 'C', 'S', 'Q', 'W', 'G', 'F', 'C', 'G', 'L',
              'G', 'S', 'E', 'F', 'C', 'G', 'G', 'G', 'C', 'Q', 'S', 'G', 'A',
              'C']], dtype='U')
                # fmt: on
            )
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_length(self):
        """Test getting the number of alignments without parsing the file."""
        stream = open(self.path)
        alignments = Align.parse(stream, "hhr")
        stream.close()
        self.assertEqual(len(alignments), 32)


class Align_hhr_allx(unittest.TestCase):
    path = os.path.join("HHsuite", "allx.hhr")

    def test_reading(self):
        alignments = Align.parse(self.path, "hhr")
        self.assertEqual(alignments.metadata["No_of_seqs"], (1, 1))
        self.assertAlmostEqual(alignments.metadata["Neff"], 1.0)
        self.assertEqual(alignments.metadata["Searched_HMMs"], 38388)
        self.assertEqual(alignments.metadata["Rundate"], "Fri Feb 15 16:24:19 2019")
        self.assertEqual(
            alignments.metadata["Command line"],
            "hhsearch -i allx.fasta -d /pdb70_hhm_db",
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 0.04)
        self.assertAlmostEqual(alignment.annotations["E-value"], 3.4e04)
        self.assertAlmostEqual(alignment.annotations["Score"], -0.01)
        self.assertAlmostEqual(alignment.annotations["Identities"], 0)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.427)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertEqual(alignment.query.id, "Only X amino acids")
        self.assertEqual(alignment.query.seq[38:39], "X")
        self.assertEqual(alignment[1], "X")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                      ~",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                       ~      ",
        )
        self.assertEqual(alignment.target.id, "1klr_A")
        self.assertEqual(alignment.target.seq[23:24], "T")
        self.assertEqual(alignment.target.annotations["hmm_name"], "1klr_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Zinc finger Y-chromosomal protein; transcription; NMR {Synthetic} SCOP: g.37.1.1 PDB: 5znf_A 1kls_A 1xrz_A* 7znf_A",
        )
        self.assertEqual(alignment[0], "T")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "                       H      ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                       H      ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[23, 24],
                             [38, 39]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['T'],
             ['X']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 0.04)
        self.assertAlmostEqual(alignment.annotations["E-value"], 3.4e04)
        self.assertAlmostEqual(alignment.annotations["Score"], 0.00)
        self.assertAlmostEqual(alignment.annotations["Identities"], 0)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.158)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertEqual(alignment.query.id, "Only X amino acids")
        self.assertEqual(alignment.query.seq[3:4], "X")
        self.assertEqual(alignment[1], "X")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "   ~                                   ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "         ~                   ",
        )
        self.assertEqual(alignment.target.id, "5ion_A")
        self.assertEqual(alignment.target.seq[9:10], "G")
        self.assertEqual(alignment.target.annotations["hmm_name"], "5ion_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Zinc finger and BTB domain-containing protein 17; C2H2 zinc finger, MIZ-1, ZBTB17, transcription factor, trans; NMR {Homo sapiens}",
        )
        self.assertEqual(alignment[0], "G")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "         T                   ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "         C                   ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 9, 10],
                             [ 3,  4]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['G'],
             ['X']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 0.04)
        self.assertAlmostEqual(alignment.annotations["E-value"], 3.4e04)
        self.assertAlmostEqual(alignment.annotations["Score"], 0.04)
        self.assertAlmostEqual(alignment.annotations["Identities"], 0)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.575)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertEqual(alignment.query.id, "Only X amino acids")
        self.assertEqual(alignment.query.seq[3:4], "X")
        self.assertEqual(alignment[1], "X")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "   ~                                   ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "~                           ",
        )
        self.assertEqual(alignment.target.id, "2jvx_A")
        self.assertEqual(alignment.target.seq[0:1], "S")
        self.assertEqual(alignment.target.annotations["hmm_name"], "2jvx_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "NF-kappa-B essential modulator; CCHC classical zinc finger, NEMO zinc finger, beta-BETA- alpha fold, coiled coil, cytoplasm, disease mutation; NMR {Synthetic} PDB: 2jvy_A",
        )
        self.assertEqual(alignment[0], "S")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "C                           ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "C                           ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[0, 1],
                             [3, 4]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['S'],
             ['X']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 0.04)
        self.assertAlmostEqual(alignment.annotations["E-value"], 3.4e04)
        self.assertAlmostEqual(alignment.annotations["Score"], 0.03)
        self.assertAlmostEqual(alignment.annotations["Identities"], 0)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.158)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertEqual(alignment.query.id, "Only X amino acids")
        self.assertEqual(alignment.query.seq[3:4], "X")
        self.assertEqual(alignment[1], "X")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "   ~                                   ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "          ~                  ",
        )
        self.assertEqual(alignment.target.id, "2ab3_A")
        self.assertEqual(alignment.target.seq[10:11], "G")
        self.assertEqual(alignment.target.annotations["hmm_name"], "2ab3_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "ZNF29; zinc finger protein, beta BETA alpha, RREIIB-TR, RNA binding protein; NMR {Escherichia coli} SCOP: k.12.1.1 PDB: 2ab7_A",
        )
        self.assertEqual(alignment[0], "G")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "          C                  ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "          C                  ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[10, 11],
                             [ 3,  4]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['G'],
             ['X']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 0.04)
        self.assertAlmostEqual(alignment.annotations["E-value"], 3.5e04)
        self.assertAlmostEqual(alignment.annotations["Score"], 0.02)
        self.assertAlmostEqual(alignment.annotations["Identities"], 0)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.427)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertEqual(alignment.query.id, "Only X amino acids")
        self.assertEqual(alignment.query.seq[3:4], "X")
        self.assertEqual(alignment[1], "X")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "   ~                                   ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "     ~                         ",
        )
        self.assertEqual(alignment.target.id, "1sp2_A")
        self.assertEqual(alignment.target.seq[5:6], "T")
        self.assertEqual(alignment.target.annotations["hmm_name"], "1sp2_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "SP1F2; zinc finger, transcription activation; NMR {Homo sapiens} SCOP: g.37.1.1 PDB: 1va2_A",
        )
        self.assertEqual(alignment[0], "T")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "     C                         ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "     C                         ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[5, 6],
                             [3, 4]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['T'],
             ['X']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 0.04)
        self.assertAlmostEqual(alignment.annotations["E-value"], 3.5e04)
        self.assertAlmostEqual(alignment.annotations["Score"], 0.00)
        self.assertAlmostEqual(alignment.annotations["Identities"], 0)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.251)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertEqual(alignment.query.id, "Only X amino acids")
        self.assertEqual(alignment.query.seq[37:38], "X")
        self.assertEqual(alignment[1], "X")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                     ~ ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                            ~ ",
        )
        self.assertEqual(alignment.target.id, "2lvr_A")
        self.assertEqual(alignment.target.seq[28:29], "E")
        self.assertEqual(alignment.target.annotations["hmm_name"], "2lvr_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Zinc finger and BTB domain-containing protein 17; C2H2 zinc finger, classical zinc finger, transcription; NMR {Homo sapiens}",
        )
        self.assertEqual(alignment[0], "E")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "                            C ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                            C ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[28, 29],
                             [37, 38]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['E'],
             ['X']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 0.04)
        self.assertAlmostEqual(alignment.annotations["E-value"], 3.5e04)
        self.assertAlmostEqual(alignment.annotations["Score"], 0.02)
        self.assertAlmostEqual(alignment.annotations["Identities"], 0)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.155)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertEqual(alignment.query.id, "Only X amino acids")
        self.assertEqual(alignment.query.seq[1:2], "X")
        self.assertEqual(alignment[1], "X")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            " ~                                     ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "~                           ",
        )
        self.assertEqual(alignment.target.id, "2kvf_A")
        self.assertEqual(alignment.target.seq[0:1], "M")
        self.assertEqual(alignment.target.annotations["hmm_name"], "2kvf_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Zinc finger and BTB domain-containing protein 32; protein/DNA, metal-binding, transcription; NMR {Mus musculus}",
        )
        self.assertEqual(alignment[0], "M")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "C                           ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "C                           ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[0, 1],
                             [1, 2]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['M'],
             ['X']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 0.03)
        self.assertAlmostEqual(alignment.annotations["E-value"], 3.6e04)
        self.assertAlmostEqual(alignment.annotations["Score"], 0.02)
        self.assertAlmostEqual(alignment.annotations["Identities"], 0)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.170)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertEqual(alignment.query.id, "Only X amino acids")
        self.assertEqual(alignment.query.seq[35:36], "X")
        self.assertEqual(alignment[1], "X")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                   ~   ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                ~         ",
        )
        self.assertEqual(alignment.target.id, "1dsq_A")
        self.assertEqual(alignment.target.seq[16:17], "D")
        self.assertEqual(alignment.target.annotations["hmm_name"], "1dsq_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Nucleic acid binding protein P14; CCHC type zinc finger, virus/viral protein; NMR {Mouse mammary tumor virus} SCOP: g.40.1.1",
        )
        self.assertEqual(alignment[0], "D")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"], "                T         "
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"], "                h         "
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[16, 17],
                             [35, 36]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['D'],
             ['X']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 0.03)
        self.assertAlmostEqual(alignment.annotations["E-value"], 3.6e04)
        self.assertAlmostEqual(alignment.annotations["Score"], 0.03)
        self.assertAlmostEqual(alignment.annotations["Identities"], 0)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.209)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertEqual(alignment.query.id, "Only X amino acids")
        self.assertEqual(alignment.query.seq[35:36], "X")
        self.assertEqual(alignment[1], "X")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                   ~   ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                            ~",
        )
        self.assertEqual(alignment.target.id, "5a7u_A")
        self.assertEqual(alignment.target.seq[28:29], "N")
        self.assertEqual(alignment.target.annotations["hmm_name"], "5a7u_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "ADR1, regulatory protein ADR1; protein folding, translation, ribosome, zinc finger, SECM, translational arrest peptide, cryo-EM; 4.80A {Saccharomyces cerevisiae} PDB: 1paa_A",
        )
        self.assertEqual(alignment[0], "N")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "                             ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                            C",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[28, 29],
                             [35, 36]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['N'],
             ['X']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 0.03)
        self.assertAlmostEqual(alignment.annotations["E-value"], 3.6e04)
        self.assertAlmostEqual(alignment.annotations["Score"], 0.03)
        self.assertAlmostEqual(alignment.annotations["Identities"], 0)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.170)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertEqual(alignment.query.id, "Only X amino acids")
        self.assertEqual(alignment.query.seq[3:4], "X")
        self.assertEqual(alignment[1], "X")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "   ~                                   ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "~                               ",
        )
        self.assertEqual(alignment.target.id, "1zfd_A")
        self.assertEqual(alignment.target.seq[0:1], "D")
        self.assertEqual(alignment.target.annotations["hmm_name"], "1zfd_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "SWI5; DNA binding motif, zinc finger DNA binding domain; NMR {Saccharomyces cerevisiae} SCOP: g.37.1.1",
        )
        self.assertEqual(alignment[0], "D")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "C                               ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "C                               ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[0, 1],
                             [3, 4]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['D'],
             ['X']], dtype='U')
                # fmt: on
            )
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_length(self):
        """Test getting the number of alignments without parsing the file."""
        stream = open(self.path)
        alignments = Align.parse(stream, "hhr")
        stream.close()
        self.assertEqual(len(alignments), 10)


class Align_hhr_4p79_hhsearch_server_NOssm(unittest.TestCase):
    path = os.path.join("HHsuite", "4p79_hhsearch_server_NOssm.hhr")

    def test_reading(self):
        alignments = Align.parse(self.path, "hhr")
        self.assertEqual(alignments.metadata["No_of_seqs"], (110, 1051))
        self.assertAlmostEqual(alignments.metadata["Neff"], 10.453)
        self.assertEqual(alignments.metadata["Searched_HMMs"], 46616)
        self.assertEqual(alignments.metadata["Rundate"], "Thu Nov 29 16:33:45 2018")
        self.assertEqual(
            alignments.metadata["Command line"],
            "hhsearch -cpu 8 -i ../results/full.a3m -d /cluster/toolkit/production/databases/hh-suite/mmcif70/pdb70 -o ../results/4P79_1.hhr -oa3m ../results/4P79_1.a3m -p 50 -Z 250 -loc -z 1 -b 1 -B 250 -ssm 0 -sc 1 -seq 1 -dbstrlen 10000 -norealign -maxres 32000 -contxt /cluster/toolkit/production/bioprogs/tools/hh-suite-build/data/context_data.crf",
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 99.94)
        self.assertAlmostEqual(alignment.annotations["E-value"], 6.8e-32)
        self.assertAlmostEqual(alignment.annotations["Score"], 194.63)
        self.assertAlmostEqual(alignment.annotations["Identities"], 100)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.536)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 10.300)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "CcHHHHHHHHHHHHHHHHHHHHHHHHHHhCcCCcccccCCCccccceEEeccchhhhccCCCCeecccchhHhcCChhHHHHHHHHHHHHHHHHHHHHHHHHHHHhhcccCCChHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHhcCCCcCCCCCceecchHHHHHHHHHHHHHHHHHHHhhhhhcCCCCcccC",
        )
        self.assertEqual(alignment.query.id, "4P79:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[0:198],
            "GSEFMSVAVETFGFFMSALGLLMLGLTLSNSYWRVSTVHGNVITTNTIFENLWYSCATDSLGVSNCWDFPSMLALSGYVQGCRALMITAILLGFLGLFLGMVGLRATNVGNMDLSKKAKLLAIAGTLHILAGACGMVAISWYAVNITTDFFNPLYAGTKYELGPALYLGWSASLLSILGGICVFSTAAASSKEEPATR",
        )
        self.assertEqual(
            alignment[1],
            "GSEFMSVAVETFGFFMSALGLLMLGLTLSNSYWRVSTVHGNVITTNTIFENLWYSCATDSLGVSNCWDFPSMLALSGYVQGCRALMITAILLGFLGLFLGMVGLRATNVGNMDLSKKAKLLAIAGTLHILAGACGMVAISWYAVNITTDFFNPLYAGTKYELGPALYLGWSASLLSILGGICVFSTAAASSKEEPATR",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "~~~m~~~~l~~~~~~l~~~a~~l~~~a~~t~~W~~~~~~~~~~~~~~~~~GLw~~C~~~~~~~~~C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~l~~la~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~g~sf~l~~~~~~l~~~~~~l~~~~~~~~~~~~~~~~",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "~~~m~~~~~~~~~~~l~~~a~~l~~va~~tp~W~~~~~~~~~~~~~~~~~GLw~~C~~~~~~~~~C~~~~~~~~~~~~~~~~~~~~~~s~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~la~~~~~i~~~~~~~~~~~~~~~~~~~~~~~~~g~s~~l~~~a~~l~~~~~i~~~~~~~~~~~~~~~~~",
        )
        self.assertEqual(alignment.target.id, "4P79_A")
        self.assertEqual(
            alignment.target.seq[0:198],
            "GSEFMSVAVETFGFFMSALGLLMLGLTLSNSYWRVSTVHGNVITTNTIFENLWYSCATDSLGVSNCWDFPSMLALSGYVQGCRALMITAILLGFLGLFLGMVGLRATNVGNMDLSKKAKLLAIAGTLHILAGACGMVAISWYAVNITTDFFNPLYAGTKYELGPALYLGWSASLLSILGGICVFSTAAASSKEEPATR",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "4P79_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "cell adhesion protein; cell adhesion, tight junction, membrane; HET: OLC, MSE; 2.4A {Mus musculus}",
        )
        self.assertEqual(
            alignment[0],
            "GSEFMSVAVETFGFFMSALGLLMLGLTLSNSYWRVSTVHGNVITTNTIFENLWYSCATDSLGVSNCWDFPSMLALSGYVQGCRALMITAILLGFLGLFLGMVGLRATNVGNMDLSKKAKLLAIAGTLHILAGACGMVAISWYAVNITTDFFNPLYAGTKYELGPALYLGWSASLLSILGGICVFSTAAASSKEEPATR",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "CHHHHHHHHHHHHHHHHHHHHHHHHHHTCSCSEEECCEEEECSSEEEEECTTSCEEEEECCHHHHHHSSHHHHHHHHHHHHHHHHHHHHHHHHHSSSCCSSCCCHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHTCSSCCSCCEEECHHHHHHHHHHHHHHHHHHHHHHHHHTC                 ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "CChHHHHHHHHHHHHHHHHHHHHHHHHHcCCcccccccCCCccccceEEeccchhhhhcCCCCcccccchHHhcCCHHHHHHHHHHHHHHHHHHHHHHHHHHHHHhccccCCChhHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHhccCcCCCCCCcccchhHHHHHHHHHHHHHHHHHHHhchhccCCCCccCC",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  0, 198],
                             [  0, 198]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['G', 'S', 'E', 'F', 'M', 'S', 'V', 'A', 'V', 'E', 'T', 'F', 'G',
              'F', 'F', 'M', 'S', 'A', 'L', 'G', 'L', 'L', 'M', 'L', 'G', 'L',
              'T', 'L', 'S', 'N', 'S', 'Y', 'W', 'R', 'V', 'S', 'T', 'V', 'H',
              'G', 'N', 'V', 'I', 'T', 'T', 'N', 'T', 'I', 'F', 'E', 'N', 'L',
              'W', 'Y', 'S', 'C', 'A', 'T', 'D', 'S', 'L', 'G', 'V', 'S', 'N',
              'C', 'W', 'D', 'F', 'P', 'S', 'M', 'L', 'A', 'L', 'S', 'G', 'Y',
              'V', 'Q', 'G', 'C', 'R', 'A', 'L', 'M', 'I', 'T', 'A', 'I', 'L',
              'L', 'G', 'F', 'L', 'G', 'L', 'F', 'L', 'G', 'M', 'V', 'G', 'L',
              'R', 'A', 'T', 'N', 'V', 'G', 'N', 'M', 'D', 'L', 'S', 'K', 'K',
              'A', 'K', 'L', 'L', 'A', 'I', 'A', 'G', 'T', 'L', 'H', 'I', 'L',
              'A', 'G', 'A', 'C', 'G', 'M', 'V', 'A', 'I', 'S', 'W', 'Y', 'A',
              'V', 'N', 'I', 'T', 'T', 'D', 'F', 'F', 'N', 'P', 'L', 'Y', 'A',
              'G', 'T', 'K', 'Y', 'E', 'L', 'G', 'P', 'A', 'L', 'Y', 'L', 'G',
              'W', 'S', 'A', 'S', 'L', 'L', 'S', 'I', 'L', 'G', 'G', 'I', 'C',
              'V', 'F', 'S', 'T', 'A', 'A', 'A', 'S', 'S', 'K', 'E', 'E', 'P',
              'A', 'T', 'R'],
             ['G', 'S', 'E', 'F', 'M', 'S', 'V', 'A', 'V', 'E', 'T', 'F', 'G',
              'F', 'F', 'M', 'S', 'A', 'L', 'G', 'L', 'L', 'M', 'L', 'G', 'L',
              'T', 'L', 'S', 'N', 'S', 'Y', 'W', 'R', 'V', 'S', 'T', 'V', 'H',
              'G', 'N', 'V', 'I', 'T', 'T', 'N', 'T', 'I', 'F', 'E', 'N', 'L',
              'W', 'Y', 'S', 'C', 'A', 'T', 'D', 'S', 'L', 'G', 'V', 'S', 'N',
              'C', 'W', 'D', 'F', 'P', 'S', 'M', 'L', 'A', 'L', 'S', 'G', 'Y',
              'V', 'Q', 'G', 'C', 'R', 'A', 'L', 'M', 'I', 'T', 'A', 'I', 'L',
              'L', 'G', 'F', 'L', 'G', 'L', 'F', 'L', 'G', 'M', 'V', 'G', 'L',
              'R', 'A', 'T', 'N', 'V', 'G', 'N', 'M', 'D', 'L', 'S', 'K', 'K',
              'A', 'K', 'L', 'L', 'A', 'I', 'A', 'G', 'T', 'L', 'H', 'I', 'L',
              'A', 'G', 'A', 'C', 'G', 'M', 'V', 'A', 'I', 'S', 'W', 'Y', 'A',
              'V', 'N', 'I', 'T', 'T', 'D', 'F', 'F', 'N', 'P', 'L', 'Y', 'A',
              'G', 'T', 'K', 'Y', 'E', 'L', 'G', 'P', 'A', 'L', 'Y', 'L', 'G',
              'W', 'S', 'A', 'S', 'L', 'L', 'S', 'I', 'L', 'G', 'G', 'I', 'C',
              'V', 'F', 'S', 'T', 'A', 'A', 'A', 'S', 'S', 'K', 'E', 'E', 'P',
              'A', 'T', 'R']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 99.88)
        self.assertAlmostEqual(alignment.annotations["E-value"], 1.8e-27)
        self.assertAlmostEqual(alignment.annotations["Score"], 169.81)
        self.assertAlmostEqual(alignment.annotations["Identities"], 38)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.775)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 10.300)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "   HHHHHHHHHHHHHHHHHHHHHHHHHhCcCCcccccCCCccccceEEeccchhhhccCCCCeecccchhHhcCChhHHHHHHHHHHHHHHHHHHHHHHHHHHHhhcccCCChHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHhcCCCcCCCCCceecchHHHHHHHHHHHHHHHHHHHhh            ",
        )
        self.assertEqual(alignment.query.id, "4P79:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[3:186],
            "FMSVAVETFGFFMSALGLLMLGLTLSNSYWRVSTVHGNVITTNTIFENLWYSCATDSLGVSNCWDFPSMLALSGYVQGCRALMITAILLGFLGLFLGMVGLRATNVGNMDLSKKAKLLAIAGTLHILAGACGMVAISWYAVNITTDFFNPLYAGTKYELGPALYLGWSASLLSILGGICVFST",
        )
        self.assertEqual(
            alignment[1],
            "FMSVAVETFGFFMSALGLLMLGLTLSNSYWRVSTVHGN-VITTNTIFENLWYSCATDSLGVSNCWDFPSMLALSGYVQGCRALMITAILLGFLGLFLGMVGLRATNVGNMDLSKKAKLLAIAGTLHILAGACGMVAISWYAVNITTDFFNPLY-AGTKYELGPALYLGWSASLLSILGGICVFST",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "   m~~~~l~~~~~~l~~~a~~l~~~a~~t~~W~~~~~~~~~~~~~~~~~GLw~~C~~~~~~~~~C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~l~~la~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~g~sf~l~~~~~~l~~~~~~l~~~~            ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "M~~~~~~~~~~~l~~~~~~l~i~a~~tp~W~~~~~~~~~~~~~~~~~~GLw~~C~~~~~~~~~C~~~~~~~~~~~~~~~~~~~~~~s~~~~~~~~~~~i~~~~~~~~~~~~~~~~~~~~~~~~~l~~~a~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~g~s~~l~~~~~~l~~~~~~l~~~~",
        )
        self.assertEqual(alignment.target.id, "3X29_C")
        self.assertEqual(
            alignment.target.seq[0:185],
            "MANSGLQLLGYFLALGGWVGIIASTALPQWKQSSYAGDAIITAVGLYEGLWMSCASQSTGQVQCKLYDSLLALDGHIQSARALMVVAVLLGFVAMVLSVVGMKATRVGDSNPTAKSRVAISGGALFLLAGLCTLTAVSWYATLVTQEFFNPSTPVNARYEFGPALFVGWASAGLAMLGGSFLAAT",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "3X29_C")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Claudin-19, Heat-labile enterotoxin B chain; TOXIN-CELL ADHESION COMPLEX, TIGHT JUNCTION; 3.7A {Mus musculus}; Related PDB entries: 3X29_A",
        )
        self.assertEqual(
            alignment[0],
            "MANSGLQLLGYFLALGGWVGIIASTALPQWKQSSYAGDAIITAVGLYEGLWMSCASQSTGQVQCKLYDSLLALDGHIQSARALMVVAVLLGFVAMVLSVVGMKATRVGDSNPTAKSRVAISGGALFLLAGLCTLTAVSWYATLVTQEFFNPSTPVNARYEFGPALFVGWASAGLAMLGGSFLAAT",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "CHHHHHHHHHHHHHHHHHHHHHHCSCSEEEEECSSSSSSCEEEEECSSEEECCCCCEEECCCCHHHHHHHHHHHHHHHHHHHHHHHHHHHCCHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHCTTSCTTEEEEECTHHHHHHHHHHHHHHHHHHHHTC                      ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "ChhHHHHHHHHHHHHHHHHHHHHHHhCchhhhcccCCCcccccceeeeeccchhhhcCCCcEEceechhHhcCCHHHHHHHHHHHHHHHHHHHHHHHHHHcchhhhcCCCChhHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHhcCCCCCCCCcccccHHHHHHHHHHHHHHHHHHHHhcC",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  0,  38,  39, 153, 154, 185],
                             [  3,  41,  41, 155, 155, 186]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['M', 'A', 'N', 'S', 'G', 'L', 'Q', 'L', 'L', 'G', 'Y', 'F', 'L',
              'A', 'L', 'G', 'G', 'W', 'V', 'G', 'I', 'I', 'A', 'S', 'T', 'A',
              'L', 'P', 'Q', 'W', 'K', 'Q', 'S', 'S', 'Y', 'A', 'G', 'D', 'A',
              'I', 'I', 'T', 'A', 'V', 'G', 'L', 'Y', 'E', 'G', 'L', 'W', 'M',
              'S', 'C', 'A', 'S', 'Q', 'S', 'T', 'G', 'Q', 'V', 'Q', 'C', 'K',
              'L', 'Y', 'D', 'S', 'L', 'L', 'A', 'L', 'D', 'G', 'H', 'I', 'Q',
              'S', 'A', 'R', 'A', 'L', 'M', 'V', 'V', 'A', 'V', 'L', 'L', 'G',
              'F', 'V', 'A', 'M', 'V', 'L', 'S', 'V', 'V', 'G', 'M', 'K', 'A',
              'T', 'R', 'V', 'G', 'D', 'S', 'N', 'P', 'T', 'A', 'K', 'S', 'R',
              'V', 'A', 'I', 'S', 'G', 'G', 'A', 'L', 'F', 'L', 'L', 'A', 'G',
              'L', 'C', 'T', 'L', 'T', 'A', 'V', 'S', 'W', 'Y', 'A', 'T', 'L',
              'V', 'T', 'Q', 'E', 'F', 'F', 'N', 'P', 'S', 'T', 'P', 'V', 'N',
              'A', 'R', 'Y', 'E', 'F', 'G', 'P', 'A', 'L', 'F', 'V', 'G', 'W',
              'A', 'S', 'A', 'G', 'L', 'A', 'M', 'L', 'G', 'G', 'S', 'F', 'L',
              'A', 'A', 'T'],
             ['F', 'M', 'S', 'V', 'A', 'V', 'E', 'T', 'F', 'G', 'F', 'F', 'M',
              'S', 'A', 'L', 'G', 'L', 'L', 'M', 'L', 'G', 'L', 'T', 'L', 'S',
              'N', 'S', 'Y', 'W', 'R', 'V', 'S', 'T', 'V', 'H', 'G', 'N', '-',
              'V', 'I', 'T', 'T', 'N', 'T', 'I', 'F', 'E', 'N', 'L', 'W', 'Y',
              'S', 'C', 'A', 'T', 'D', 'S', 'L', 'G', 'V', 'S', 'N', 'C', 'W',
              'D', 'F', 'P', 'S', 'M', 'L', 'A', 'L', 'S', 'G', 'Y', 'V', 'Q',
              'G', 'C', 'R', 'A', 'L', 'M', 'I', 'T', 'A', 'I', 'L', 'L', 'G',
              'F', 'L', 'G', 'L', 'F', 'L', 'G', 'M', 'V', 'G', 'L', 'R', 'A',
              'T', 'N', 'V', 'G', 'N', 'M', 'D', 'L', 'S', 'K', 'K', 'A', 'K',
              'L', 'L', 'A', 'I', 'A', 'G', 'T', 'L', 'H', 'I', 'L', 'A', 'G',
              'A', 'C', 'G', 'M', 'V', 'A', 'I', 'S', 'W', 'Y', 'A', 'V', 'N',
              'I', 'T', 'T', 'D', 'F', 'F', 'N', 'P', 'L', 'Y', '-', 'A', 'G',
              'T', 'K', 'Y', 'E', 'L', 'G', 'P', 'A', 'L', 'Y', 'L', 'G', 'W',
              'S', 'A', 'S', 'L', 'L', 'S', 'I', 'L', 'G', 'G', 'I', 'C', 'V',
              'F', 'S', 'T']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 99.84)
        self.assertAlmostEqual(alignment.annotations["E-value"], 1.2e-25)
        self.assertAlmostEqual(alignment.annotations["Score"], 176.50)
        self.assertAlmostEqual(alignment.annotations["Identities"], 37)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.639)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 9.700)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            " cHHHHHHHHHHHHHHHHHHHHHHHHHHhCcCCcccccCCCccccceEEeccchhhhccCCCCeecccchhHhcCChhHHHHHHHHHHHHHHHHHHHHHHHHHHHhhcccCCChHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHhcCCCcCCCCCceecchHHHHHHHHHHHHHHHHHHHhhhhhcCCCCccc ",
        )
        self.assertEqual(alignment.query.id, "4P79:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[1:197],
            "SEFMSVAVETFGFFMSALGLLMLGLTLSNSYWRVSTVHGNVITTNTIFENLWYSCATDSLGVSNCWDFPSMLALSGYVQGCRALMITAILLGFLGLFLGMVGLRATNVGNMDLSKKAKLLAIAGTLHILAGACGMVAISWYAVNITTDFFNPLYAGTKYELGPALYLGWSASLLSILGGICVFSTAAASSKEEPAT",
        )
        self.assertEqual(
            alignment[1],
            "SEFMSVAVETFGFFMSALGLLMLGLTLSNSYWRVSTVHGN-VITTNTIFENLWYSCATDSLGVSNCWDFPSMLALSGYVQGCRALMITAILLGFLGLFLGMVGLRATNVGNMDLSKKAKLLAIAGTLHILAGACGMVAISWYAVNITTDFFNPLY-AGTKYELGPALYLGWSASLLSILGGICVFSTAAASSKEEPAT",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            " ~~m~~~~l~~~~~~l~~~a~~l~~~a~~t~~W~~~~~~~~~~~~~~~~~GLw~~C~~~~~~~~~C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~l~~la~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~g~sf~l~~~~~~l~~~~~~l~~~~~~~~~~~~~~~ ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                                                                                                                                                                       ~~m~~~~~~~~~~~l~~~a~il~~va~~tp~W~~~~~~~~~~~~~~~~~~GLw~~C~~~~~~~~~C~~~~~~~~~~~~~~~~~~~~i~s~i~~~~a~i~~i~~~~~~~~~~~~~~~~~~~~~~~~~~~la~l~~~~~v~~~~~~~~~~~~~~~~~~~~~~~~GwS~~l~~~s~~l~lia~~l~~~~~~~~~~~~~~~",
        )
        self.assertEqual(alignment.target.id, "5B2G_A")
        self.assertEqual(
            alignment.target.seq[167:364],
            "KGMASMGLQVMGIALAVLGWLAVMLCCALPMWRVTAFIGSNIVTSQTIWEGLWMNCVVQSTGQMQCKVYDSLLALPQDLQAARALVIISIIVAALGVLLSVVGGKCTNCLEDESAKAKTMIVAGVVFLLAGLMVIVPVSWTAHNIIQDFYNPLVASGQKREMGASLYVGWAASGLLLLGGGLLCCSGPSSGENLYFQ",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "5B2G_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Endolysin,Claudin-4 (E.C.3.2.1.17), Heat-labile enterotoxin B; Membrane protein, Complex, Cell-free protein; HET: MSE; 3.5A {Enterobacteria phage T4}; Related PDB entries: 5B2G_E 5B2G_C 5B2G_G",
        )
        self.assertEqual(
            alignment[0],
            "KGMASMGLQVMGIALAVLGWLAVMLCCALPMWRVTAFIGSNIVTSQTIWEGLWMNCVVQSTGQMQCKVYDSLLALPQDLQAARALVIISIIVAALGVLLSVVGGKCTNCLED-ESAKAKTMIVAGVVFLLAGLMVIVPVSWTAHNIIQDFYNPLVASGQKREMGASLYVGWAASGLLLLGGGLLCCSGPSSGENLYFQ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "                                                                                                                                                                       CCHHHHHHHHHHHHHHHHHHHHSSCSEEEEECSTTCSSCEEEEECSSEEEEEESSSSCEEEECCSSCSCCHHHHHHHHHHHHHHHHHHHHHHTTCCCCCCCCTTHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHCTTSCSSEEEEECTHHHHHHHHHHHHHHHHHHHHCC                         ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                                                                                                                                                                       hhhHHHHHHHHHHHHHHHHHHHHHHHHhCchhhhhcccCCcccCceeEEEeccchheeecCCceEeeechhhhcCCHHHHHHHHHHHHHHHHHHHHHHHHHHhhhcccCCCChhhHHHHHHHHHHHHHHHHHHhHHHHHHHHHHHHHHhcCcccCCCCCcccchHHHHHHHHHHHHHHHHHHHHhcCCCCCCccccC",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[167, 207, 208, 279, 279, 321, 322, 364],
                             [  1,  41,  41, 112, 113, 155, 155, 197]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['K', 'G', 'M', 'A', 'S', 'M', 'G', 'L', 'Q', 'V', 'M', 'G', 'I',
              'A', 'L', 'A', 'V', 'L', 'G', 'W', 'L', 'A', 'V', 'M', 'L', 'C',
              'C', 'A', 'L', 'P', 'M', 'W', 'R', 'V', 'T', 'A', 'F', 'I', 'G',
              'S', 'N', 'I', 'V', 'T', 'S', 'Q', 'T', 'I', 'W', 'E', 'G', 'L',
              'W', 'M', 'N', 'C', 'V', 'V', 'Q', 'S', 'T', 'G', 'Q', 'M', 'Q',
              'C', 'K', 'V', 'Y', 'D', 'S', 'L', 'L', 'A', 'L', 'P', 'Q', 'D',
              'L', 'Q', 'A', 'A', 'R', 'A', 'L', 'V', 'I', 'I', 'S', 'I', 'I',
              'V', 'A', 'A', 'L', 'G', 'V', 'L', 'L', 'S', 'V', 'V', 'G', 'G',
              'K', 'C', 'T', 'N', 'C', 'L', 'E', 'D', '-', 'E', 'S', 'A', 'K',
              'A', 'K', 'T', 'M', 'I', 'V', 'A', 'G', 'V', 'V', 'F', 'L', 'L',
              'A', 'G', 'L', 'M', 'V', 'I', 'V', 'P', 'V', 'S', 'W', 'T', 'A',
              'H', 'N', 'I', 'I', 'Q', 'D', 'F', 'Y', 'N', 'P', 'L', 'V', 'A',
              'S', 'G', 'Q', 'K', 'R', 'E', 'M', 'G', 'A', 'S', 'L', 'Y', 'V',
              'G', 'W', 'A', 'A', 'S', 'G', 'L', 'L', 'L', 'L', 'G', 'G', 'G',
              'L', 'L', 'C', 'C', 'S', 'G', 'P', 'S', 'S', 'G', 'E', 'N', 'L',
              'Y', 'F', 'Q'],
             ['S', 'E', 'F', 'M', 'S', 'V', 'A', 'V', 'E', 'T', 'F', 'G', 'F',
              'F', 'M', 'S', 'A', 'L', 'G', 'L', 'L', 'M', 'L', 'G', 'L', 'T',
              'L', 'S', 'N', 'S', 'Y', 'W', 'R', 'V', 'S', 'T', 'V', 'H', 'G',
              'N', '-', 'V', 'I', 'T', 'T', 'N', 'T', 'I', 'F', 'E', 'N', 'L',
              'W', 'Y', 'S', 'C', 'A', 'T', 'D', 'S', 'L', 'G', 'V', 'S', 'N',
              'C', 'W', 'D', 'F', 'P', 'S', 'M', 'L', 'A', 'L', 'S', 'G', 'Y',
              'V', 'Q', 'G', 'C', 'R', 'A', 'L', 'M', 'I', 'T', 'A', 'I', 'L',
              'L', 'G', 'F', 'L', 'G', 'L', 'F', 'L', 'G', 'M', 'V', 'G', 'L',
              'R', 'A', 'T', 'N', 'V', 'G', 'N', 'M', 'D', 'L', 'S', 'K', 'K',
              'A', 'K', 'L', 'L', 'A', 'I', 'A', 'G', 'T', 'L', 'H', 'I', 'L',
              'A', 'G', 'A', 'C', 'G', 'M', 'V', 'A', 'I', 'S', 'W', 'Y', 'A',
              'V', 'N', 'I', 'T', 'T', 'D', 'F', 'F', 'N', 'P', 'L', 'Y', '-',
              'A', 'G', 'T', 'K', 'Y', 'E', 'L', 'G', 'P', 'A', 'L', 'Y', 'L',
              'G', 'W', 'S', 'A', 'S', 'L', 'L', 'S', 'I', 'L', 'G', 'G', 'I',
              'C', 'V', 'F', 'S', 'T', 'A', 'A', 'A', 'S', 'S', 'K', 'E', 'E',
              'P', 'A', 'T']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 99.82)
        self.assertAlmostEqual(alignment.annotations["E-value"], 1.1e-24)
        self.assertAlmostEqual(alignment.annotations["Score"], 168.32)
        self.assertAlmostEqual(alignment.annotations["Identities"], 15)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.164)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 8.600)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "CcHHHHHHHHHHHHHHHHHHHHHHHHHHhCcCCcccccCCCccccceEEeccchhhhccCCCCeecccchhHhcCChhHHHHHHHHHHHHHHHHHHHHHHHHHHHhhcccCCChHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHhcCCCcCCCCCceecchHHHHHHHHHHHHHHHHHHHhhhhhcCCCCcccC",
        )
        self.assertEqual(alignment.query.id, "4P79:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[0:198],
            "GSEFMSVAVETFGFFMSALGLLMLGLTLSNSYWRVSTVHGNVITTNTIFENLWYSCATDSLGVSNCWDFPSMLALSGYVQGCRALMITAILLGFLGLFLGMVGLRATNVGNMDLSKKAKLLAIAGTLHILAGACGMVAISWYAVNITTDFFNPLYAGTKYELGPALYLGWSASLLSILGGICVFSTAAASSKEEPATR",
        )
        self.assertEqual(
            alignment[1],
            "GSEFMSVAVETFGFFMSALGLLMLGLTLSNSYWRVSTVHGN-----------VITTNTIFENLWYSCATDSLGVSNCWDFPSMLAL-----------SGYVQGCRALMITAILLGFLGLFLGMVGLRATNVGNMDLSKKAKLLAIAGTLHILAGACGMVAISWYAVNITTDFFNPLYAGTKYELGPALYLGWSASLLSILGGICVFSTAAASSKEEPATR",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "~~~m~~~~l~~~~~~l~~~a~~l~~~a~~t~~W~~~~~~~~~~~~~~~~~GLw~~C~~~~~~~~~C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~l~~la~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~g~sf~l~~~~~~l~~~~~~l~~~~~~~~~~~~~~~~",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "m~~~~r~~l~~l~~vl~~~a~vL~~ia~~tp~W~~~~~~~~~~~~~~~~~~~~~~~~~~~~GLw~~C~~~~~~~~~C~~~~~f~~~~~~~~~~~~~~~~~~r~~~~f~ilsl~l~~l~~l~~~~~~~~~~~~~~~~~agil~~lAgl~~~i~vily~~~~~~~~~~~~~~~~~~~~GwSf~~~~~s~~l~~la~vl~v~~~~~~~~~~~~~~                                                                                                               ",
        )
        self.assertEqual(alignment.target.id, "5KK2_G")
        self.assertEqual(
            alignment.target.seq[0:212],
            "MGLFDRGVQMLLTTVGAFAAFSLMTIAVGTDYWLYSRGVCKTKSVSENETSKKNEEVMTHSGLWRTCCLEGNFKGLCKQIDHFPEDADYEADTAEYFLRAVRASSIFPILSVILLFMGGLCIAASEFYKTRHNIILSAGIFFVSAGLSNIIGIIVYISANAGDPSKSDSKKNSYSYGWSFYFGALSFIIAEMVGVLAVHMFIDRHKQLRATA",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "5KK2_G")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Glutamate receptor 2, Voltage-dependent calcium; Tetrameric, ionotropic, glutamate receptors, GluA2; 7.3A {Rattus norvegicus}; Related PDB entries: 5VOT_E 5VOT_H 5VOU_H 5KK2_H 5VOV_E 5VOT_F 5KK2_F 5VOT_G 5VOV_G 5VOU_E 5VOV_H 5VOU_G 5VOV_F 5VOU_F 5KK2_E",
        )
        self.assertEqual(
            alignment[0],
            "MGLFDRGVQMLLTTVGAFAAFSLMTIAVGTDYWLYSRGVCKTKSVSENETSKKNEEVMTHSGLWRTCCLEGNFKGLCKQIDHFPEDADYEADTAEYFLRAVRASSIFPILSVILLFMGGLCIAASEFY--------KTRHNIILSAGIFFVSAGLSNIIGIIVYISANAGDPSKSDSKKNSYSYGWSFYFGALSFIIAEMVGVLAVHMFIDRHKQLRATA",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "CHHHHHHHHHHHHHHHHHHHHHHCCCCCBCCCCCCEECCSSEEEECCCCEEEECCCCCCCCCCSSHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHTTCCCHHHHHHHHHHHHHHHHHHHHHHHHHHHGGGCCCCCHHHHHHHHHHHHHHHHHHHHHHHHCCCCCCCHHHH                                                                                                                                                        ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "CchhhHHHHHHHHHHHHHHHHHHHHHHHhcCCceeeecccCCCCCcccccccccceeeeeccccchhccccccCCceeeccCCCCCccCCcccHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHcccchHHHHHHHHHHHHHHHHHHHHHHHHHHHhcCCcccCccCCCcccccHHHHHHHHHHHHHHHHHHHHHHHHhhhhHHHHHHH                                                                                                               ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  0,  41,  52,  86,  97, 128, 128, 212],
                             [  0,  41,  41,  75,  75, 106, 114, 198]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['M', 'G', 'L', 'F', 'D', 'R', 'G', 'V', 'Q', 'M', 'L', 'L', 'T',
              'T', 'V', 'G', 'A', 'F', 'A', 'A', 'F', 'S', 'L', 'M', 'T', 'I',
              'A', 'V', 'G', 'T', 'D', 'Y', 'W', 'L', 'Y', 'S', 'R', 'G', 'V',
              'C', 'K', 'T', 'K', 'S', 'V', 'S', 'E', 'N', 'E', 'T', 'S', 'K',
              'K', 'N', 'E', 'E', 'V', 'M', 'T', 'H', 'S', 'G', 'L', 'W', 'R',
              'T', 'C', 'C', 'L', 'E', 'G', 'N', 'F', 'K', 'G', 'L', 'C', 'K',
              'Q', 'I', 'D', 'H', 'F', 'P', 'E', 'D', 'A', 'D', 'Y', 'E', 'A',
              'D', 'T', 'A', 'E', 'Y', 'F', 'L', 'R', 'A', 'V', 'R', 'A', 'S',
              'S', 'I', 'F', 'P', 'I', 'L', 'S', 'V', 'I', 'L', 'L', 'F', 'M',
              'G', 'G', 'L', 'C', 'I', 'A', 'A', 'S', 'E', 'F', 'Y', '-', '-',
              '-', '-', '-', '-', '-', '-', 'K', 'T', 'R', 'H', 'N', 'I', 'I',
              'L', 'S', 'A', 'G', 'I', 'F', 'F', 'V', 'S', 'A', 'G', 'L', 'S',
              'N', 'I', 'I', 'G', 'I', 'I', 'V', 'Y', 'I', 'S', 'A', 'N', 'A',
              'G', 'D', 'P', 'S', 'K', 'S', 'D', 'S', 'K', 'K', 'N', 'S', 'Y',
              'S', 'Y', 'G', 'W', 'S', 'F', 'Y', 'F', 'G', 'A', 'L', 'S', 'F',
              'I', 'I', 'A', 'E', 'M', 'V', 'G', 'V', 'L', 'A', 'V', 'H', 'M',
              'F', 'I', 'D', 'R', 'H', 'K', 'Q', 'L', 'R', 'A', 'T', 'A'],
             ['G', 'S', 'E', 'F', 'M', 'S', 'V', 'A', 'V', 'E', 'T', 'F', 'G',
              'F', 'F', 'M', 'S', 'A', 'L', 'G', 'L', 'L', 'M', 'L', 'G', 'L',
              'T', 'L', 'S', 'N', 'S', 'Y', 'W', 'R', 'V', 'S', 'T', 'V', 'H',
              'G', 'N', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',
              'V', 'I', 'T', 'T', 'N', 'T', 'I', 'F', 'E', 'N', 'L', 'W', 'Y',
              'S', 'C', 'A', 'T', 'D', 'S', 'L', 'G', 'V', 'S', 'N', 'C', 'W',
              'D', 'F', 'P', 'S', 'M', 'L', 'A', 'L', '-', '-', '-', '-', '-',
              '-', '-', '-', '-', '-', '-', 'S', 'G', 'Y', 'V', 'Q', 'G', 'C',
              'R', 'A', 'L', 'M', 'I', 'T', 'A', 'I', 'L', 'L', 'G', 'F', 'L',
              'G', 'L', 'F', 'L', 'G', 'M', 'V', 'G', 'L', 'R', 'A', 'T', 'N',
              'V', 'G', 'N', 'M', 'D', 'L', 'S', 'K', 'K', 'A', 'K', 'L', 'L',
              'A', 'I', 'A', 'G', 'T', 'L', 'H', 'I', 'L', 'A', 'G', 'A', 'C',
              'G', 'M', 'V', 'A', 'I', 'S', 'W', 'Y', 'A', 'V', 'N', 'I', 'T',
              'T', 'D', 'F', 'F', 'N', 'P', 'L', 'Y', 'A', 'G', 'T', 'K', 'Y',
              'E', 'L', 'G', 'P', 'A', 'L', 'Y', 'L', 'G', 'W', 'S', 'A', 'S',
              'L', 'L', 'S', 'I', 'L', 'G', 'G', 'I', 'C', 'V', 'F', 'S', 'T',
              'A', 'A', 'A', 'S', 'S', 'K', 'E', 'E', 'P', 'A', 'T', 'R']],
            dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 99.81)
        self.assertAlmostEqual(alignment.annotations["E-value"], 1.2e-24)
        self.assertAlmostEqual(alignment.annotations["Score"], 159.28)
        self.assertAlmostEqual(alignment.annotations["Identities"], 17)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.231)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 10.600)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            " cHHHHHHHHHHHHHHHHHHHHHHHHHHhCcCCcccccCCCccccceEEeccchhhhccCCCCeecccchhHhcCChhHHHHHHHHHHHHHHHHHHHHHHHHHHHhhcccCCChHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHhcCCCcCCCCCceecchHHHHHHHHHHHHHHHHHHHhhhhhcCCCCcccC",
        )
        self.assertEqual(alignment.query.id, "4P79:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[1:198],
            "SEFMSVAVETFGFFMSALGLLMLGLTLSNSYWRVSTVHGNVITTNTIFENLWYSCATDSLGVSNCWDFPSMLALSGYVQGCRALMITAILLGFLGLFLGMVGLRATNVGNMDLSKKAKLLAIAGTLHILAGACGMVAISWYAVNITTDFFNPLYAGTKYELGPALYLGWSASLLSILGGICVFSTAAASSKEEPATR",
        )
        self.assertEqual(
            alignment[1],
            "SEFMSVAVETFGFFMSALGLLMLGLTLSNSYWRVSTVHGNV--ITTNTIFENLWYSCATDSLG-------------VSNCWDFPSML-----------ALSGYVQGCRALMITAILLGFLGLFLGMVGLRATNVGNMDLSKKAKLLAIAGTLHILAGACGMVAISWYAVNITTDFFNPLYAGTKYELGPALYLGWSASLLSILGGICVFSTAAASSKEEPATR",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            " ~~m~~~~l~~~~~~l~~~a~~l~~~a~~t~~W~~~~~~~~~~~~~~~~~GLw~~C~~~~~~~~~C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~l~~la~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~g~sf~l~~~~~~l~~~~~~l~~~~~~~~~~~~~~~~",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "m~~~~~~~~~~~~~l~~~~~~l~~ia~~t~~W~~~~~~~~~~~~~~~~~~~GLw~~C~~~~~~~~~~~c~~~~~~~~~~C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~a~~~~~~~~~~~~~~~~~~~~~~~~~~~~l~~ls~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~g~s~~l~~~a~~l~~l~~~l~~~~~~~~~~~~~~~~        ",
        )
        self.assertEqual(alignment.target.id, "5GJV_E")
        self.assertEqual(
            alignment.target.seq[0:214],
            "MSPTEAPKVRVTLFCILVGIVLAMTAVVSDHWAVLSPHMENHNTTCEAAHFGLWRICTKRIALGEDRSCGPITLPGEKNCSYFRHFNPGESSEIFEFTTQKEYSISAAAISVFSLGFLIMGTICALMAFRKKRDYLLRPASMFYVFAGLCLFVSLEVMRQSVKRMIDSEDTVWIEYYYSWSFACACAAFVLLFLGGISLLLFSLPRMPQNPWES",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "5GJV_E")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Voltage-dependent L-type calcium channel subunit; complex, channel, membrane protein; HET: PC1, BMA, NAG; 3.6A {Oryctolagus cuniculus}; Related PDB entries: 3JBR_E 5GJW_E",
        )
        self.assertEqual(
            alignment[0],
            "MSPTEAPKVRVTLFCILVGIVLAMTAVVSDHWAVLSPHMENHNTTCEAAHFGLWRICTKRIALGEDRSCGPITLPGEKNCSYFRHFNPGESSEIFEFTTQKEYSISAAAISVFSLGFLIMGTICALMAFRK---------KRDYLLRPASMFYVFAGLCLFVSLEVMRQSVKRMIDSEDTVWIEYYYSWSFACACAAFVLLFLGGISLLLFSLPRMPQNPWES",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "CHHHHSHHHHHHHHHHHHHHHHHHHHHHCCCSEEEECCCCEEECCSSCEEEECCCEEEECCCCHHHHHHHHHHHHHHHHHHHHHHHTTTTCCTTHHHHHHHHHHHHHHHHHHHHHHHHHHHHCCEEEECHHHHHHHHHHHHHHHHHHHHHHHHSCCCCSSTTSC                                                          ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "CCccccHHHHHHHHHHHHHHHHHHHHHhCcchhhcCCCCCCCCcchheeecccHHHhchhhccCCCCCCccccCCCcccceeeeccCCCCCCCCCCcccHHHHHHHHHHHHHHHHHHHHHHHHHHHHHccCCChhhHHHHHHHHHHHHHHHHHHHHHHHHHHHHhcCCCCCcceEeecchHHHHHHHHHHHHHHHHHHHHHHcccCCCCCcchh        ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  0,  41,  43,  63,  76,  87,  98, 131, 131, 214],
                             [  1,  42,  42,  62,  62,  73,  73, 106, 115, 198]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['M', 'S', 'P', 'T', 'E', 'A', 'P', 'K', 'V', 'R', 'V', 'T', 'L',
              'F', 'C', 'I', 'L', 'V', 'G', 'I', 'V', 'L', 'A', 'M', 'T', 'A',
              'V', 'V', 'S', 'D', 'H', 'W', 'A', 'V', 'L', 'S', 'P', 'H', 'M',
              'E', 'N', 'H', 'N', 'T', 'T', 'C', 'E', 'A', 'A', 'H', 'F', 'G',
              'L', 'W', 'R', 'I', 'C', 'T', 'K', 'R', 'I', 'A', 'L', 'G', 'E',
              'D', 'R', 'S', 'C', 'G', 'P', 'I', 'T', 'L', 'P', 'G', 'E', 'K',
              'N', 'C', 'S', 'Y', 'F', 'R', 'H', 'F', 'N', 'P', 'G', 'E', 'S',
              'S', 'E', 'I', 'F', 'E', 'F', 'T', 'T', 'Q', 'K', 'E', 'Y', 'S',
              'I', 'S', 'A', 'A', 'A', 'I', 'S', 'V', 'F', 'S', 'L', 'G', 'F',
              'L', 'I', 'M', 'G', 'T', 'I', 'C', 'A', 'L', 'M', 'A', 'F', 'R',
              'K', '-', '-', '-', '-', '-', '-', '-', '-', '-', 'K', 'R', 'D',
              'Y', 'L', 'L', 'R', 'P', 'A', 'S', 'M', 'F', 'Y', 'V', 'F', 'A',
              'G', 'L', 'C', 'L', 'F', 'V', 'S', 'L', 'E', 'V', 'M', 'R', 'Q',
              'S', 'V', 'K', 'R', 'M', 'I', 'D', 'S', 'E', 'D', 'T', 'V', 'W',
              'I', 'E', 'Y', 'Y', 'Y', 'S', 'W', 'S', 'F', 'A', 'C', 'A', 'C',
              'A', 'A', 'F', 'V', 'L', 'L', 'F', 'L', 'G', 'G', 'I', 'S', 'L',
              'L', 'L', 'F', 'S', 'L', 'P', 'R', 'M', 'P', 'Q', 'N', 'P', 'W',
              'E', 'S'],
             ['S', 'E', 'F', 'M', 'S', 'V', 'A', 'V', 'E', 'T', 'F', 'G', 'F',
              'F', 'M', 'S', 'A', 'L', 'G', 'L', 'L', 'M', 'L', 'G', 'L', 'T',
              'L', 'S', 'N', 'S', 'Y', 'W', 'R', 'V', 'S', 'T', 'V', 'H', 'G',
              'N', 'V', '-', '-', 'I', 'T', 'T', 'N', 'T', 'I', 'F', 'E', 'N',
              'L', 'W', 'Y', 'S', 'C', 'A', 'T', 'D', 'S', 'L', 'G', '-', '-',
              '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', 'V', 'S',
              'N', 'C', 'W', 'D', 'F', 'P', 'S', 'M', 'L', '-', '-', '-', '-',
              '-', '-', '-', '-', '-', '-', '-', 'A', 'L', 'S', 'G', 'Y', 'V',
              'Q', 'G', 'C', 'R', 'A', 'L', 'M', 'I', 'T', 'A', 'I', 'L', 'L',
              'G', 'F', 'L', 'G', 'L', 'F', 'L', 'G', 'M', 'V', 'G', 'L', 'R',
              'A', 'T', 'N', 'V', 'G', 'N', 'M', 'D', 'L', 'S', 'K', 'K', 'A',
              'K', 'L', 'L', 'A', 'I', 'A', 'G', 'T', 'L', 'H', 'I', 'L', 'A',
              'G', 'A', 'C', 'G', 'M', 'V', 'A', 'I', 'S', 'W', 'Y', 'A', 'V',
              'N', 'I', 'T', 'T', 'D', 'F', 'F', 'N', 'P', 'L', 'Y', 'A', 'G',
              'T', 'K', 'Y', 'E', 'L', 'G', 'P', 'A', 'L', 'Y', 'L', 'G', 'W',
              'S', 'A', 'S', 'L', 'L', 'S', 'I', 'L', 'G', 'G', 'I', 'C', 'V',
              'F', 'S', 'T', 'A', 'A', 'A', 'S', 'S', 'K', 'E', 'E', 'P', 'A',
              'T', 'R']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 99.79)
        self.assertAlmostEqual(alignment.annotations["E-value"], 5.7e-24)
        self.assertAlmostEqual(alignment.annotations["Score"], 156.78)
        self.assertAlmostEqual(alignment.annotations["Identities"], 14)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.065)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 9.000)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "    HHHHHHHHHHHHHHHHHHHHHHHHhCcCCcccccCCCccccceEEeccchhhhccCCCCeecccchhHhcCChhHHHHHHHHHHHHHHHHHHHHHHHHHHHhhcccCCChHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHhcCCCcCCCCCceecchHHHHHHHHHHHHHHHHHHHhhhhhcCCCCcccC",
        )
        self.assertEqual(alignment.query.id, "4P79:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[4:198],
            "MSVAVETFGFFMSALGLLMLGLTLSNSYWRVSTVHGNVITTNTIFENLWYSCATDSLGVSNCWDFPSMLALSGYVQGCRALMITAILLGFLGLFLGMVGLRATNVGNMDLSKKAKLLAIAGTLHILAGACGMVAISWYAVNITTDFFNPLYAGTKYELGPALYLGWSASLLSILGGICVFSTAAASSKEEPATR",
        )
        self.assertEqual(
            alignment[1],
            "MSVAVETFGFFMSALGLLMLGLTLSNSYWRVSTVHGNVITTNTIFENLWYSCATDSL-GVSNCWDFPSMLAL--SGYVQGCRALMITAILLGFLGLFLGMVGLRATNVGNMDLSKKAKLLAIAGTLHILAGACGMVAISWYAVNITTDFFN-------PLYAGTKYELGPALYLGWSASLLSILGGICVFSTAAASSKEEPATR",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "    ~~~~l~~~~~~l~~~a~~l~~~a~~t~~W~~~~~~~~~~~~~~~~~GLw~~C~~~~~~~~~C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~l~~la~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~g~sf~l~~~~~~l~~~~~~l~~~~~~~~~~~~~~~~",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                           ~~~~~~~~~~~l~~~a~~l~~va~~tp~W~~~~~~~~~~~~~GLw~~C~~~~~~~~~~C~~~~~~~~~~~~~~~~a~~~l~i~~~~l~~~~~i~~~~~~~~~~~~~~~~~~~~~~lagl~~lv~~ivf~~~~~~~~~~~~~~~~~~~~~~~~~~~GwSf~la~~a~~l~~la~il~~~~~~~~~~~~~~~~          ",
        )
        self.assertEqual(alignment.target.id, "6C14_D")
        self.assertEqual(
            alignment.target.seq[27:218],
            "NSRAVGVMWGTLTICFSVLVMALFIQPYWIGDSVSTPQAGYFGLFSYCVGNVLSSELICKGGPLDFSSIPSRAFKTAMFFVALAMFLIIGSIICFSLFFVCNTATVYKICAWMQLAAATGLMIGCLVYPDGWDSSEVRRMCGEQTGKYTLGHCTIRWAFMLAILSIGDALILSFLAFVLGYRQDKLLPDDY",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "6C14_D")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Protocadherin-15, LHFPL tetraspan subfamily member; PCDH15, LHFPL5, protocadherin, tip link;{Mus musculus}; Related PDB entries: 6C14_B",
        )
        self.assertEqual(
            alignment[0],
            "NSRAVGVMWGTLTICFSVLVMALFIQPYWIGDSVSTP----QAGYFGLFSYCVGNVLSSELICKGGPLDFSSIPSRAFKTAMFFVALAMFLIIGSIICFSLFFVC---------NTATVYKICAWMQLAAATGLMIGCLVYPDGWDSSEVRRMCGEQTGKYTLGHCTIRWAFMLAILSIGDALILSFLAFVLGYRQDKLLPDDY",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "                           HHHHHHHHHHHHHHHHHHHHHHHHSCSCSEECCCCEECCCSSSCBEECTTCCCCCBCCSSCCSSSSSCHHHHHHHHHHHHHHHHHHHHHHHHHTTTTSCHHHHHHHHHHHHHHHHHHHHHHHHHHHHHCCCCEECHHHHHHHHHHHHHHHHHHHHHHC                                           ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                           HhHHHHHHHHHHHHHHHHHHHHHHHChHhhCCCCCCCCCceeehhhhhccccCCcceeeecCCCCHhhCCCHHHHHHHHHHHHHHHHHHHHHHHHHHHHHcccchHHHHHHHHHHHHHHHHHHHHHHcccCCCCHHHHHHhcCCCCCCCCCCCcHHHHHHHHHHHHHHHHHHHHHHHHHHHhhhhcCchhh          ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 27,  64,  64,  80,  81,  95,  97, 128, 128, 165, 172, 218],
                             [  4,  41,  45,  61,  61,  75,  75, 106, 115, 152, 152, 198]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['N', 'S', 'R', 'A', 'V', 'G', 'V', 'M', 'W', 'G', 'T', 'L', 'T',
              'I', 'C', 'F', 'S', 'V', 'L', 'V', 'M', 'A', 'L', 'F', 'I', 'Q',
              'P', 'Y', 'W', 'I', 'G', 'D', 'S', 'V', 'S', 'T', 'P', '-', '-',
              '-', '-', 'Q', 'A', 'G', 'Y', 'F', 'G', 'L', 'F', 'S', 'Y', 'C',
              'V', 'G', 'N', 'V', 'L', 'S', 'S', 'E', 'L', 'I', 'C', 'K', 'G',
              'G', 'P', 'L', 'D', 'F', 'S', 'S', 'I', 'P', 'S', 'R', 'A', 'F',
              'K', 'T', 'A', 'M', 'F', 'F', 'V', 'A', 'L', 'A', 'M', 'F', 'L',
              'I', 'I', 'G', 'S', 'I', 'I', 'C', 'F', 'S', 'L', 'F', 'F', 'V',
              'C', '-', '-', '-', '-', '-', '-', '-', '-', '-', 'N', 'T', 'A',
              'T', 'V', 'Y', 'K', 'I', 'C', 'A', 'W', 'M', 'Q', 'L', 'A', 'A',
              'A', 'T', 'G', 'L', 'M', 'I', 'G', 'C', 'L', 'V', 'Y', 'P', 'D',
              'G', 'W', 'D', 'S', 'S', 'E', 'V', 'R', 'R', 'M', 'C', 'G', 'E',
              'Q', 'T', 'G', 'K', 'Y', 'T', 'L', 'G', 'H', 'C', 'T', 'I', 'R',
              'W', 'A', 'F', 'M', 'L', 'A', 'I', 'L', 'S', 'I', 'G', 'D', 'A',
              'L', 'I', 'L', 'S', 'F', 'L', 'A', 'F', 'V', 'L', 'G', 'Y', 'R',
              'Q', 'D', 'K', 'L', 'L', 'P', 'D', 'D', 'Y'],
             ['M', 'S', 'V', 'A', 'V', 'E', 'T', 'F', 'G', 'F', 'F', 'M', 'S',
              'A', 'L', 'G', 'L', 'L', 'M', 'L', 'G', 'L', 'T', 'L', 'S', 'N',
              'S', 'Y', 'W', 'R', 'V', 'S', 'T', 'V', 'H', 'G', 'N', 'V', 'I',
              'T', 'T', 'N', 'T', 'I', 'F', 'E', 'N', 'L', 'W', 'Y', 'S', 'C',
              'A', 'T', 'D', 'S', 'L', '-', 'G', 'V', 'S', 'N', 'C', 'W', 'D',
              'F', 'P', 'S', 'M', 'L', 'A', 'L', '-', '-', 'S', 'G', 'Y', 'V',
              'Q', 'G', 'C', 'R', 'A', 'L', 'M', 'I', 'T', 'A', 'I', 'L', 'L',
              'G', 'F', 'L', 'G', 'L', 'F', 'L', 'G', 'M', 'V', 'G', 'L', 'R',
              'A', 'T', 'N', 'V', 'G', 'N', 'M', 'D', 'L', 'S', 'K', 'K', 'A',
              'K', 'L', 'L', 'A', 'I', 'A', 'G', 'T', 'L', 'H', 'I', 'L', 'A',
              'G', 'A', 'C', 'G', 'M', 'V', 'A', 'I', 'S', 'W', 'Y', 'A', 'V',
              'N', 'I', 'T', 'T', 'D', 'F', 'F', 'N', '-', '-', '-', '-', '-',
              '-', '-', 'P', 'L', 'Y', 'A', 'G', 'T', 'K', 'Y', 'E', 'L', 'G',
              'P', 'A', 'L', 'Y', 'L', 'G', 'W', 'S', 'A', 'S', 'L', 'L', 'S',
              'I', 'L', 'G', 'G', 'I', 'C', 'V', 'F', 'S', 'T', 'A', 'A', 'A',
              'S', 'S', 'K', 'E', 'E', 'P', 'A', 'T', 'R']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 99.53)
        self.assertAlmostEqual(alignment.annotations["E-value"], 1.8e-18)
        self.assertAlmostEqual(alignment.annotations["Score"], 150.79)
        self.assertAlmostEqual(alignment.annotations["Identities"], 13)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.177)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 11.600)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            " cHHHHHHHHHHHHHHHHHHHHHHHHHHhCcCCcccccCCCccccceEEeccchhhhccCCCCeecccchhHhcCChhHHHHHHHHHHHHHHHHHHHHHHHHHHHhhcccCCChHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHhcCCCcCCCCCceecchHHHHHHHHHHHHHHHHHHHhhhhhcCCCCcccC",
        )
        self.assertEqual(alignment.query.id, "4P79:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[1:198],
            "SEFMSVAVETFGFFMSALGLLMLGLTLSNSYWRVSTVHGNVITTNTIFENLWYSCATDSLGVSNCWDFPSMLALSGYVQGCRALMITAILLGFLGLFLGMVGLRATNVGNMDLSKKAKLLAIAGTLHILAGACGMVAISWYAVNITTDFFNPLYAGTKYELGPALYLGWSASLLSILGGICVFSTAAASSKEEPATR",
        )
        self.assertEqual(
            alignment[1],
            "SEFMSVAVETFGFFMSALGLLMLGLTLSNSYW-----------RVSTVHGNVITTNTIFENLWYSCATDSLGVSNCWDFPSML-----------ALSGYVQGCRALMITAILLGFLGLFLGMVGLRATNVGNMDLSKKAKLLAIAGTLHILAGACGMVAISWYAVNITTDFFNPLYAGTKYELGPALYLGWSASLLSILGGICVFSTAAASSKEEPATR",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            " ~~m~~~~l~~~~~~l~~~a~~l~~~a~~t~~W~~~~~~~~~~~~~~~~~GLw~~C~~~~~~~~~C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~l~~la~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~g~sf~l~~~~~~l~~~~~~l~~~~~~~~~~~~~~~~",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   ~~~~~~~~~~~~~~~~~~~~~l~~~~~~~~~W~~~~~~~~~~~~~~~~~~~~~~~~~~~~GLw~~C~~~~~~~~~C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~g~s~~l~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ",
        )
        self.assertEqual(alignment.target.id, "6DLZ_D")
        self.assertEqual(
            alignment.target.seq[819:1030],
            "GLFDRGVQMLLTTVGAFAAFSLMTIAVGTDYWLYSRGVCKTKSVSEDETSKKNEEVMTHSGLWRTCCLEGNFKGLCKQIDHFPEDADYEADTAEYFLRAVRASSIFPILSVILLFMGGLCIAASEFYKTRHNIILSAGIFFVSAGLSNIIGIIVYISANAGDPSKSDSKKNSYSYGWSFYFGALSFIIAEMVGVLAVHMFIDRHKQLTGGA",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "6DLZ_D")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Glutamate receptor 2,Voltage-dependent calcium channel; Ion channel, TRANSPORT PROTEIN; HET: GLU, CYZ; 3.9A {Rattus norvegicus}; Related PDB entries: 5VHY_B 5WEL_C 5WEL_D 5VHZ_E 5WEK_D 5WEN_C 5VHX_E 5VHY_C 5VHZ_F 5VHY_D 5VHZ_B 5WEN_D 5WEM_C 5WEL_B 5WEM_A 5VHX_C 5WEM_D 5WEK_C 5WEK_A 5VHY_A 5VHW_C 5VHX_A 5VHZ_C 5VHW_D 5VHY_E 5WEN_A 5VHX_B 5VHX_D 5VHW_B 5WEL_A 5VHW_A 5VHZ_D 5WEN_B 5WEK_B 5VHY_F 5VHZ_A 5WEM_B 5WEO_B 5KBT_D 5KBS_A 5KBU_B 5WEO_C 5KBU_A 5KBS_C 5KBU_D 5KBT_A 5KBS_B 5KBT_C 5KBT_B 5KBU_C 5WEO_D 5WEO_A 5KBS_D 6DM2_D 6DLZ_A 6DM0_A 6DM0_B 6DM1_D 6DM2_C 6DM0_C 6DM1_B 6DLZ_B 6DM2_B 6DM1_C 6DM2_A 6DLZ_C 6DM1_A 6DM0_D",
        )
        self.assertEqual(
            alignment[0],
            "GLFDRGVQMLLTTVGAFAAFSLMTIAVGTDYWLYSRGVCKTKSVSEDETSKKNEEVMTHSGLWRTCCLEGNFKGLCKQIDHFPEDADYEADTAEYFLRAVRASSIFPILSVILLFMGGLCIAASEFY--------KTRHNIILSAGIFFVSAGLSNIIGIIVYISANAGDPSKSDSKKNSYSYGWSFYFGALSFIIAEMVGVLAVHMFIDRHKQLTGGA",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   CCSTTTHHHHHHHHHHHHHHHHHHHHSCSCSEEEECCSCCCCEEEECSSBCCBCSSSSCSCCCBCCCCCCSSTTTTCHHHHHHHHHHHSCHHHHHHHHHHHHHHHHHHTSSTTCSCCSHHHHHHHHHHHHHHHHHHHHHHHHHHHCCCCSCCEECHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHC                        ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   CchHHHHHHHHHHHHHHHHHHHHHHHhccCcceeeccccccccCCcCccccccccccccccceeeeeecCcccCceeecccCCCCcCCCCCcHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHhcccCCceehhhHHHHHHHHHHHHHHHHHHHHHhcCCcccCCCCcccceecHHHHHHHHHHHHHHHHHHHHHHHHHHHhHhhcCCC ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 819,  851,  862,  902,  913,  946,  946, 1030],
                             [   1,   33,   33,   73,   73,  106,  114,  198]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['G', 'L', 'F', 'D', 'R', 'G', 'V', 'Q', 'M', 'L', 'L', 'T', 'T',
              'V', 'G', 'A', 'F', 'A', 'A', 'F', 'S', 'L', 'M', 'T', 'I', 'A',
              'V', 'G', 'T', 'D', 'Y', 'W', 'L', 'Y', 'S', 'R', 'G', 'V', 'C',
              'K', 'T', 'K', 'S', 'V', 'S', 'E', 'D', 'E', 'T', 'S', 'K', 'K',
              'N', 'E', 'E', 'V', 'M', 'T', 'H', 'S', 'G', 'L', 'W', 'R', 'T',
              'C', 'C', 'L', 'E', 'G', 'N', 'F', 'K', 'G', 'L', 'C', 'K', 'Q',
              'I', 'D', 'H', 'F', 'P', 'E', 'D', 'A', 'D', 'Y', 'E', 'A', 'D',
              'T', 'A', 'E', 'Y', 'F', 'L', 'R', 'A', 'V', 'R', 'A', 'S', 'S',
              'I', 'F', 'P', 'I', 'L', 'S', 'V', 'I', 'L', 'L', 'F', 'M', 'G',
              'G', 'L', 'C', 'I', 'A', 'A', 'S', 'E', 'F', 'Y', '-', '-', '-',
              '-', '-', '-', '-', '-', 'K', 'T', 'R', 'H', 'N', 'I', 'I', 'L',
              'S', 'A', 'G', 'I', 'F', 'F', 'V', 'S', 'A', 'G', 'L', 'S', 'N',
              'I', 'I', 'G', 'I', 'I', 'V', 'Y', 'I', 'S', 'A', 'N', 'A', 'G',
              'D', 'P', 'S', 'K', 'S', 'D', 'S', 'K', 'K', 'N', 'S', 'Y', 'S',
              'Y', 'G', 'W', 'S', 'F', 'Y', 'F', 'G', 'A', 'L', 'S', 'F', 'I',
              'I', 'A', 'E', 'M', 'V', 'G', 'V', 'L', 'A', 'V', 'H', 'M', 'F',
              'I', 'D', 'R', 'H', 'K', 'Q', 'L', 'T', 'G', 'G', 'A'],
             ['S', 'E', 'F', 'M', 'S', 'V', 'A', 'V', 'E', 'T', 'F', 'G', 'F',
              'F', 'M', 'S', 'A', 'L', 'G', 'L', 'L', 'M', 'L', 'G', 'L', 'T',
              'L', 'S', 'N', 'S', 'Y', 'W', '-', '-', '-', '-', '-', '-', '-',
              '-', '-', '-', '-', 'R', 'V', 'S', 'T', 'V', 'H', 'G', 'N', 'V',
              'I', 'T', 'T', 'N', 'T', 'I', 'F', 'E', 'N', 'L', 'W', 'Y', 'S',
              'C', 'A', 'T', 'D', 'S', 'L', 'G', 'V', 'S', 'N', 'C', 'W', 'D',
              'F', 'P', 'S', 'M', 'L', '-', '-', '-', '-', '-', '-', '-', '-',
              '-', '-', '-', 'A', 'L', 'S', 'G', 'Y', 'V', 'Q', 'G', 'C', 'R',
              'A', 'L', 'M', 'I', 'T', 'A', 'I', 'L', 'L', 'G', 'F', 'L', 'G',
              'L', 'F', 'L', 'G', 'M', 'V', 'G', 'L', 'R', 'A', 'T', 'N', 'V',
              'G', 'N', 'M', 'D', 'L', 'S', 'K', 'K', 'A', 'K', 'L', 'L', 'A',
              'I', 'A', 'G', 'T', 'L', 'H', 'I', 'L', 'A', 'G', 'A', 'C', 'G',
              'M', 'V', 'A', 'I', 'S', 'W', 'Y', 'A', 'V', 'N', 'I', 'T', 'T',
              'D', 'F', 'F', 'N', 'P', 'L', 'Y', 'A', 'G', 'T', 'K', 'Y', 'E',
              'L', 'G', 'P', 'A', 'L', 'Y', 'L', 'G', 'W', 'S', 'A', 'S', 'L',
              'L', 'S', 'I', 'L', 'G', 'G', 'I', 'C', 'V', 'F', 'S', 'T', 'A',
              'A', 'A', 'S', 'S', 'K', 'E', 'E', 'P', 'A', 'T', 'R']],
            dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 52.07)
        self.assertAlmostEqual(alignment.annotations["E-value"], 6.7)
        self.assertAlmostEqual(alignment.annotations["Score"], 20.51)
        self.assertAlmostEqual(alignment.annotations["Identities"], 19)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.239)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 3.200)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "     HHHHHHHHHHHHHHHHHHHHHHHhCcCCcccc                                                                                                                                                                 ",
        )
        self.assertEqual(alignment.query.id, "4P79:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(alignment.query.seq[5:37], "SVAVETFGFFMSALGLLMLGLTLSNSYWRVST")
        self.assertEqual(alignment[1], "SVAVETFGFFMSALGLLMLGLTLSNS--YWRVST")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "     ~~~l~~~~~~l~~~a~~l~~~a~~t~~W~~~~                                                                                                                                                                 ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "        Rr~lval~~fl~vLallIHfiLLST~~~NWl~~~",
        )
        self.assertEqual(alignment.target.id, "5YQ7_F")
        self.assertEqual(
            alignment.target.seq[8:42], "RTSVVVSTLLGLVMALLIHFVVLSSGAFNWLRAP"
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "5YQ7_F")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Beta subunit of light-harvesting 1; Photosynthetic core complex, PHOTOSYNTHESIS; HET: MQE, BCL, HEM, KGD, BPH;{Roseiflexus castenholzii}; Related PDB entries: 5YQ7_V 5YQ7_3 5YQ7_T 5YQ7_J 5YQ7_9 5YQ7_N 5YQ7_A 5YQ7_P 5YQ7_H 5YQ7_D 5YQ7_5 5YQ7_7 5YQ7_1 5YQ7_R",
        )
        self.assertEqual(alignment[0], "RTSVVVSTLLGLVMALLIHFVVLSSGAFNWLRAP")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "        HHHHHHHHHHHHHHHHHCCCCCTTSSSCCSCC  ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "        HHHHHHHHHHHHHHHHHHHHHHHcccccccccCC",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 8, 34, 36, 42],
                             [ 5, 31, 31, 37]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['R', 'T', 'S', 'V', 'V', 'V', 'S', 'T', 'L', 'L', 'G', 'L', 'V',
              'M', 'A', 'L', 'L', 'I', 'H', 'F', 'V', 'V', 'L', 'S', 'S', 'G',
              'A', 'F', 'N', 'W', 'L', 'R', 'A', 'P'],
             ['S', 'V', 'A', 'V', 'E', 'T', 'F', 'G', 'F', 'F', 'M', 'S', 'A',
              'L', 'G', 'L', 'L', 'M', 'L', 'G', 'L', 'T', 'L', 'S', 'N', 'S',
              '-', '-', 'Y', 'W', 'R', 'V', 'S', 'T']], dtype='U')
                # fmt: on
            )
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_length(self):
        """Test getting the number of alignments without parsing the file."""
        stream = open(self.path)
        alignments = Align.parse(stream, "hhr")
        stream.close()
        self.assertEqual(len(alignments), 8)


class Align_hhr_4y9h_hhsearch_server_NOssm(unittest.TestCase):
    path = os.path.join("HHsuite", "4y9h_hhsearch_server_NOssm.hhr")

    def test_reading(self):
        alignments = Align.parse(self.path, "hhr")
        self.assertEqual(alignments.metadata["No_of_seqs"], (141, 1242))
        self.assertAlmostEqual(alignments.metadata["Neff"], 8.55177)
        self.assertEqual(alignments.metadata["Searched_HMMs"], 46616)
        self.assertEqual(alignments.metadata["Rundate"], "Thu Nov 29 16:28:55 2018")
        self.assertEqual(
            alignments.metadata["Command line"],
            "hhsearch -cpu 8 -i ../results/full.a3m -d /cluster/toolkit/production/databases/hh-suite/mmcif70/pdb70 -o ../results/4Y9H_1.hhr -oa3m ../results/4Y9H_1.a3m -p 50 -Z 250 -loc -z 1 -b 1 -B 250 -ssm 0 -sc 1 -seq 1 -dbstrlen 10000 -norealign -maxres 32000 -contxt /cluster/toolkit/production/bioprogs/tools/hh-suite-build/data/context_data.crf",
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 100.00)
        self.assertAlmostEqual(alignment.annotations["E-value"], 2.1e-48)
        self.assertAlmostEqual(alignment.annotations["Score"], 320.44)
        self.assertAlmostEqual(alignment.annotations["Identities"], 100)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.574)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 8.500)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "CCcHHHHHHHHHHHHHHHHHHHHHHhcCCCChhhHHHHHHHHHHHHHHHHHHHHHHhCCCceeeccCCCCCcchHHHHHHHHHHHHHHHHHHHHHcCCCHHHHHHHHHHHHHHHHHHHHHHHcccchHHHHHHHHHHHHHHHHHHHHHHhHHHHHHhCCHHHHHHHHHHHHHHHHHHHHHHHHHHHccCCCCCCChhHHHHHHHHHHHHHHHHHHHHHHHhhhhhC",
        )
        self.assertEqual(alignment.query.id, "4Y9H:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[0:226],
            "GRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFG",
        )
        self.assertEqual(
            alignment[1],
            "GRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFG",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "~~~~~~~~~~~~~~m~~~~~~f~~~~~~~~~~~~r~~~~l~~~i~~ia~~~Y~~ma~~~g~~~~~~~~~~~~v~~~RYv~W~ittPlll~~l~~l~~~~~~~~~~~v~~~~~mi~~g~~g~~~~~~~~~~~~~~~s~~~~~~i~~~l~~~~~~~~~~~~~~~~~~~~~l~~~~~i~W~~YPi~w~l~~~g~~~i~~~~~~~~~~ilDi~~K~~f~~~l~~~~~~~~",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            " ~~~~~~~~~~~~~~~~~~~~~f~~~~~~~~~~~~r~~~~l~~~i~~va~~~Y~~ma~~~g~~~~~~~~~~~~v~~~RYv~W~ittPlll~~l~~l~~~~~~~~~~~~~~~~~mi~~g~~g~~~~~~~~~w~~~~~s~~~~~~v~~~l~~~~~~~~~~~~~~~~~~~~~l~~~~~i~W~~YPi~w~l~~~g~~~i~~~~~~i~y~ilDi~~K~~f~~~l~~~~~~~~ ",
        )
        self.assertEqual(alignment.target.id, "5ZIM_A")
        self.assertEqual(
            alignment.target.seq[1:227],
            "GRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFG",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "5ZIM_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Bacteriorhodopsin; proton pump, membrane protein, PROTON; HET: L2P, RET; 1.25A {Halobacterium salinarum}; Related PDB entries: 1R84_A 1KG8_A 1KME_B 1KGB_A 1KG9_A 1KME_A 4X31_A 5ZIL_A 1E0P_A 4X32_A 5ZIN_A 1S53_B 1S51_B 1S53_A 1S54_A 1F50_A 1S54_B 1S51_A 1F4Z_A 5J7A_A 1S52_B 1S52_A 4Y9H_A 3T45_A 3T45_C 3T45_B 1C3W_A 1L0M_A",
        )
        self.assertEqual(
            alignment[0],
            "GRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFG",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            " CSTTHHHHHHHHHHHHHHHHHHHHHHTTCCCHHHHHHHHHHHHHHHHHHHHHHHHHTTTTEEEEEETTEEEEEETHHHHHHHHHHHHHHHHHHHHTTCCHHHHHHHHHHHHHHHHHHHHHHHCSSHHHHHHHHHHHHHHHHHHHHHHHHSCTTSSSSCCHHHHHHHHHHHHHHHHHHHHHHHHHHHSTTTTCCSCHHHHHHHHHHHHCCCCCCHHHHHHTSGGGSC ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            " CCcHHHHHHHHHHHHHHHHHHHHHHhcCCCChhhHHHHHHHHHHHHHHHHHHHHHHcCCCceeeccCCCCCcchHHHHHHHHHHHHHHHHHHHHHcCCCHHHHHHHHHHHHHHHHHHHHHHHcccchHHHHHHHHHHHHHHHHHHHHHHhHHHHHHhCCHHHHHHHHHHHHHHHHHHHHHHHHHHHccCCCCCCChhHHHHHHHHHHHHHHHHHHHHHHHHHHHhC ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  1, 227],
                             [  0, 226]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['G', 'R', 'P', 'E', 'W', 'I', 'W', 'L', 'A', 'L', 'G', 'T', 'A',
              'L', 'M', 'G', 'L', 'G', 'T', 'L', 'Y', 'F', 'L', 'V', 'K', 'G',
              'M', 'G', 'V', 'S', 'D', 'P', 'D', 'A', 'K', 'K', 'F', 'Y', 'A',
              'I', 'T', 'T', 'L', 'V', 'P', 'A', 'I', 'A', 'F', 'T', 'M', 'Y',
              'L', 'S', 'M', 'L', 'L', 'G', 'Y', 'G', 'L', 'T', 'M', 'V', 'P',
              'F', 'G', 'G', 'E', 'Q', 'N', 'P', 'I', 'Y', 'W', 'A', 'R', 'Y',
              'A', 'D', 'W', 'L', 'F', 'T', 'T', 'P', 'L', 'L', 'L', 'L', 'D',
              'L', 'A', 'L', 'L', 'V', 'D', 'A', 'D', 'Q', 'G', 'T', 'I', 'L',
              'A', 'L', 'V', 'G', 'A', 'D', 'G', 'I', 'M', 'I', 'G', 'T', 'G',
              'L', 'V', 'G', 'A', 'L', 'T', 'K', 'V', 'Y', 'S', 'Y', 'R', 'F',
              'V', 'W', 'W', 'A', 'I', 'S', 'T', 'A', 'A', 'M', 'L', 'Y', 'I',
              'L', 'Y', 'V', 'L', 'F', 'F', 'G', 'F', 'T', 'S', 'K', 'A', 'E',
              'S', 'M', 'R', 'P', 'E', 'V', 'A', 'S', 'T', 'F', 'K', 'V', 'L',
              'R', 'N', 'V', 'T', 'V', 'V', 'L', 'W', 'S', 'A', 'Y', 'P', 'V',
              'V', 'W', 'L', 'I', 'G', 'S', 'E', 'G', 'A', 'G', 'I', 'V', 'P',
              'L', 'N', 'I', 'E', 'T', 'L', 'L', 'F', 'M', 'V', 'L', 'D', 'V',
              'S', 'A', 'K', 'V', 'G', 'F', 'G', 'L', 'I', 'L', 'L', 'R', 'S',
              'R', 'A', 'I', 'F', 'G'],
             ['G', 'R', 'P', 'E', 'W', 'I', 'W', 'L', 'A', 'L', 'G', 'T', 'A',
              'L', 'M', 'G', 'L', 'G', 'T', 'L', 'Y', 'F', 'L', 'V', 'K', 'G',
              'M', 'G', 'V', 'S', 'D', 'P', 'D', 'A', 'K', 'K', 'F', 'Y', 'A',
              'I', 'T', 'T', 'L', 'V', 'P', 'A', 'I', 'A', 'F', 'T', 'M', 'Y',
              'L', 'S', 'M', 'L', 'L', 'G', 'Y', 'G', 'L', 'T', 'M', 'V', 'P',
              'F', 'G', 'G', 'E', 'Q', 'N', 'P', 'I', 'Y', 'W', 'A', 'R', 'Y',
              'A', 'D', 'W', 'L', 'F', 'T', 'T', 'P', 'L', 'L', 'L', 'L', 'D',
              'L', 'A', 'L', 'L', 'V', 'D', 'A', 'D', 'Q', 'G', 'T', 'I', 'L',
              'A', 'L', 'V', 'G', 'A', 'D', 'G', 'I', 'M', 'I', 'G', 'T', 'G',
              'L', 'V', 'G', 'A', 'L', 'T', 'K', 'V', 'Y', 'S', 'Y', 'R', 'F',
              'V', 'W', 'W', 'A', 'I', 'S', 'T', 'A', 'A', 'M', 'L', 'Y', 'I',
              'L', 'Y', 'V', 'L', 'F', 'F', 'G', 'F', 'T', 'S', 'K', 'A', 'E',
              'S', 'M', 'R', 'P', 'E', 'V', 'A', 'S', 'T', 'F', 'K', 'V', 'L',
              'R', 'N', 'V', 'T', 'V', 'V', 'L', 'W', 'S', 'A', 'Y', 'P', 'V',
              'V', 'W', 'L', 'I', 'G', 'S', 'E', 'G', 'A', 'G', 'I', 'V', 'P',
              'L', 'N', 'I', 'E', 'T', 'L', 'L', 'F', 'M', 'V', 'L', 'D', 'V',
              'S', 'A', 'K', 'V', 'G', 'F', 'G', 'L', 'I', 'L', 'L', 'R', 'S',
              'R', 'A', 'I', 'F', 'G']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 100.00)
        self.assertAlmostEqual(alignment.annotations["E-value"], 2.9e-45)
        self.assertAlmostEqual(alignment.annotations["Score"], 307.43)
        self.assertAlmostEqual(alignment.annotations["Identities"], 100)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.574)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 7.700)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "CCcHHHHHHHHHHHHHHHHHHHHHHhcCCCChhhHHHHHHHHHHHHHHHHHHHHHHhCCCceeeccCCCCCcchHHHHHHHHHHHHHHHHHHHHHcCCCHHHHHHHHHHHHHHHHHHHHHHHcccchHHHHHHHHHHHHHHHHHHHHHHhHHHHHHhCCHHHHHHHHHHHHHHHHHHHHHHHHHHHccCCCCCCChhHHHHHHHHHHHHHHHHHHHHHHHhhhhhC",
        )
        self.assertEqual(alignment.query.id, "4Y9H:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[0:226],
            "GRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFG",
        )
        self.assertEqual(
            alignment[1],
            "GRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFG",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "~~~~~~~~~~~~~~m~~~~~~f~~~~~~~~~~~~r~~~~l~~~i~~ia~~~Y~~ma~~~g~~~~~~~~~~~~v~~~RYv~W~ittPlll~~l~~l~~~~~~~~~~~v~~~~~mi~~g~~g~~~~~~~~~~~~~~~s~~~~~~i~~~l~~~~~~~~~~~~~~~~~~~~~l~~~~~i~W~~YPi~w~l~~~g~~~i~~~~~~~~~~ilDi~~K~~f~~~l~~~~~~~~",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                  ~~~~~~~~~~~~~~m~~~~l~f~~~~~~~~~~~~r~~~~l~~~i~~iaa~~Y~~ma~g~g~~~~~~~~~~r~v~~~RYvdW~iTtPlll~~l~~la~~~~~~~~~~v~~~~~mi~~G~~g~~~~~~~~kw~~~~~s~~~f~~v~~~l~~~~~~~a~~~~~~~~~~~~~l~~~~~v~W~~YPi~w~l~~eg~~~is~~~e~i~y~ilDv~~K~~f~~~ll~~~~~~~                  ",
        )
        self.assertEqual(alignment.target.id, "1M0K_A")
        self.assertEqual(
            alignment.target.seq[18:244],
            "GRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFG",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "1M0K_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "BACTERIORHODOPSIN; ION PUMP, MEMBRANE PROTEIN, RETINAL; HET: RET, SQU, LI1; 1.43A {Halobacterium salinarum} SCOP: f.13.1.1; Related PDB entries: 4XXJ_B 4XXJ_A 4XXJ_C 6G7L_A 5VN7_B 1M0M_A 4HWL_A 2ZFE_A 4FPD_A 6G7J_A 3VHZ_A 1M0L_A 2ZZL_A 4HYX_A 4HWL_B 4HYX_B 6G7K_A 6G7I_A 3VI0_A 5VN7_A 5VN9_B 5VN9_A 3COC_B 2I21_A 5B34_A 3HAQ_A 1S8L_A 1TN0_B 1O0A_A 1JV6_A 1TN0_A 5BR2_A 2I1X_A 1UCQ_A 1MGY_A 1C8R_A 3HAR_A 1P8U_A 3UTX_A 1Q5I_A 1Q5J_A 1Q5J_B 1JV7_A 1P8I_A 1S8J_A 3COC_A 3HAO_A 1PY6_B 1PXR_B 1Q5I_B 3UTY_A 2NTU_A 1PXS_B 5BR5_A 1PXS_A 3HAS_A 1PXR_A 1TN5_A 2WJL_A 3HAP_A 3COD_A 1TN5_B 3HAO_B 2I20_A 2WJK_A 1R2N_A 3UTV_A 3UTX_B 3HAN_A 1PY6_A 1VJM_A 3UTY_B 5B35_A 3COD_B 3UTW_A 2NTW_A 4OV0_A 1P8H_A 1CWQ_A 6G7H_A 5H2J_A 1BRD_A 5H2H_A 1X0I_1 5A44_A 1QM8_A 1X0S_A 5B6V_A 1BM1_A 1AP9_A 2BRD_A 1IW6_A 3MBV_A 5H2P_A 1QKP_A 5H2O_A 5H2N_A 5B6X_A 1FBB_A 1CWQ_B 4MD2_A 1QKO_A 1IW9_A 1QHJ_A 1AT9_A 1DZE_A 4MD1_A 5H2K_A 1X0K_1 5B6Z_A 5H2L_A 5H2I_A 1FBK_A 1IXF_A 3NS0_A 1BRX_A 3NSB_A 5H2M_A 2AT9_A 5B6Y_A 5B6W_A 5A45_A 1XJI_A 1BRR_B 1BRR_A 1BRR_C",
        )
        self.assertEqual(
            alignment[0],
            "GRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFG",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "                  CSTTHHHHHHHHHHHHHHHHHHHHHTTTCCCHHHHHHHHHHHHHHHHHHHHHHHHHTTTTCEEEEETTEEEEECTHHHHHHHHHHHHHHHHHHHHTTCCHHHHHHHHHHHHHHHHHHHHHHHCSSHHHHHHHHHHHHHHHHHHHHHHHHCCCCCHHHHHHHHHHHHHHHHHHHHHHHHHHHSTTTTCCSCHHHHHHHHHHHHHCCCCCHHHHHHHSGGGCC                       ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                  CchHHHHHHHHHHHHHHHHHHHHHHHcCCCChhhHHHHHHHHHHHHHHHHHHHHHHcCCCceeccCCCCcCccchHHHHHHHhHHHHHHHHHHHHhCCCHHHHHHHHHHHHHHHHHHHHHHHcccchhhHHHHHHHHHHHHHHHHHHHHhHHHHHHhCCHHHHHHHHHHHHHHHHHHHHHHHHHHHccCCCCCCCchHHHHHHHHHHHHHHHHHHHHHHHHHHHhc                  ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 18, 244],
                             [  0, 226]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['G', 'R', 'P', 'E', 'W', 'I', 'W', 'L', 'A', 'L', 'G', 'T', 'A',
              'L', 'M', 'G', 'L', 'G', 'T', 'L', 'Y', 'F', 'L', 'V', 'K', 'G',
              'M', 'G', 'V', 'S', 'D', 'P', 'D', 'A', 'K', 'K', 'F', 'Y', 'A',
              'I', 'T', 'T', 'L', 'V', 'P', 'A', 'I', 'A', 'F', 'T', 'M', 'Y',
              'L', 'S', 'M', 'L', 'L', 'G', 'Y', 'G', 'L', 'T', 'M', 'V', 'P',
              'F', 'G', 'G', 'E', 'Q', 'N', 'P', 'I', 'Y', 'W', 'A', 'R', 'Y',
              'A', 'D', 'W', 'L', 'F', 'T', 'T', 'P', 'L', 'L', 'L', 'L', 'D',
              'L', 'A', 'L', 'L', 'V', 'D', 'A', 'D', 'Q', 'G', 'T', 'I', 'L',
              'A', 'L', 'V', 'G', 'A', 'D', 'G', 'I', 'M', 'I', 'G', 'T', 'G',
              'L', 'V', 'G', 'A', 'L', 'T', 'K', 'V', 'Y', 'S', 'Y', 'R', 'F',
              'V', 'W', 'W', 'A', 'I', 'S', 'T', 'A', 'A', 'M', 'L', 'Y', 'I',
              'L', 'Y', 'V', 'L', 'F', 'F', 'G', 'F', 'T', 'S', 'K', 'A', 'E',
              'S', 'M', 'R', 'P', 'E', 'V', 'A', 'S', 'T', 'F', 'K', 'V', 'L',
              'R', 'N', 'V', 'T', 'V', 'V', 'L', 'W', 'S', 'A', 'Y', 'P', 'V',
              'V', 'W', 'L', 'I', 'G', 'S', 'E', 'G', 'A', 'G', 'I', 'V', 'P',
              'L', 'N', 'I', 'E', 'T', 'L', 'L', 'F', 'M', 'V', 'L', 'D', 'V',
              'S', 'A', 'K', 'V', 'G', 'F', 'G', 'L', 'I', 'L', 'L', 'R', 'S',
              'R', 'A', 'I', 'F', 'G'],
             ['G', 'R', 'P', 'E', 'W', 'I', 'W', 'L', 'A', 'L', 'G', 'T', 'A',
              'L', 'M', 'G', 'L', 'G', 'T', 'L', 'Y', 'F', 'L', 'V', 'K', 'G',
              'M', 'G', 'V', 'S', 'D', 'P', 'D', 'A', 'K', 'K', 'F', 'Y', 'A',
              'I', 'T', 'T', 'L', 'V', 'P', 'A', 'I', 'A', 'F', 'T', 'M', 'Y',
              'L', 'S', 'M', 'L', 'L', 'G', 'Y', 'G', 'L', 'T', 'M', 'V', 'P',
              'F', 'G', 'G', 'E', 'Q', 'N', 'P', 'I', 'Y', 'W', 'A', 'R', 'Y',
              'A', 'D', 'W', 'L', 'F', 'T', 'T', 'P', 'L', 'L', 'L', 'L', 'D',
              'L', 'A', 'L', 'L', 'V', 'D', 'A', 'D', 'Q', 'G', 'T', 'I', 'L',
              'A', 'L', 'V', 'G', 'A', 'D', 'G', 'I', 'M', 'I', 'G', 'T', 'G',
              'L', 'V', 'G', 'A', 'L', 'T', 'K', 'V', 'Y', 'S', 'Y', 'R', 'F',
              'V', 'W', 'W', 'A', 'I', 'S', 'T', 'A', 'A', 'M', 'L', 'Y', 'I',
              'L', 'Y', 'V', 'L', 'F', 'F', 'G', 'F', 'T', 'S', 'K', 'A', 'E',
              'S', 'M', 'R', 'P', 'E', 'V', 'A', 'S', 'T', 'F', 'K', 'V', 'L',
              'R', 'N', 'V', 'T', 'V', 'V', 'L', 'W', 'S', 'A', 'Y', 'P', 'V',
              'V', 'W', 'L', 'I', 'G', 'S', 'E', 'G', 'A', 'G', 'I', 'V', 'P',
              'L', 'N', 'I', 'E', 'T', 'L', 'L', 'F', 'M', 'V', 'L', 'D', 'V',
              'S', 'A', 'K', 'V', 'G', 'F', 'G', 'L', 'I', 'L', 'L', 'R', 'S',
              'R', 'A', 'I', 'F', 'G']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 100.00)
        self.assertAlmostEqual(alignment.annotations["E-value"], 4.4e-43)
        self.assertAlmostEqual(alignment.annotations["Score"], 291.03)
        self.assertAlmostEqual(alignment.annotations["Identities"], 55)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.956)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 8.200)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            " CcHHHHHHHHHHHHHHHHHHHHHHhcCCCChhhHHHHHHHHHHHHHHHHHHHHHHhCCCceeeccCCCCCcchHHHHHHHHHHHHHHHHHHHHHcCCCHHHHHHHHHHHHHHHHHHHHHHHcccchHHHHHHHHHHHHHHHHHHHHHHhHHHHHHhCCHHHHHHHHHHHHHHHHHHHHHHHHHHHccCCCCCCChhHHHHHHHHHHHHHHHHHHHHHHHhhhhhC",
        )
        self.assertEqual(alignment.query.id, "4Y9H:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[1:226],
            "RPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFG",
        )
        self.assertEqual(
            alignment[1],
            "RPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEG-AGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFG",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            " ~~~~~~~~~~~~~m~~~~~~f~~~~~~~~~~~~r~~~~l~~~i~~ia~~~Y~~ma~~~g~~~~~~~~~~~~v~~~RYv~W~ittPlll~~l~~l~~~~~~~~~~~v~~~~~mi~~g~~g~~~~~~~~~~~~~~~s~~~~~~i~~~l~~~~~~~~~~~~~~~~~~~~~l~~~~~i~W~~YPi~w~l~~~g~~~i~~~~~~~~~~ilDi~~K~~f~~~l~~~~~~~~",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "     ~~~~~~lw~~~~~~~~~~~~f~~~~~~~~~~~~R~~~~l~~~i~~va~~~Y~~ma~~~g~~~~~~~~~~~~v~~~RYv~W~iTtPlll~~l~~lag~~~~~~~~~~~~~~~mi~~G~~g~~~~~~~~~w~~~~~s~~~~~~i~~~l~~~~~~~a~~~~~~~~~~~~~l~~~~~v~W~~YPi~w~l~~~g~~~~i~~~~~~i~~~vlDv~~K~~fg~~l~~~~~~~~          ",
        )
        self.assertEqual(alignment.target.id, "4FBZ_A")
        self.assertEqual(
            alignment.target.seq[5:231],
            "GPESIWLWIGTIGMTLGTLYFVGRGRGVRDRKMQEFYIITTFITTIAAAMYFAMATGFGVTEVVVGDEALTIYWARYADWLFTTPLLLLDLGLLAGANRNTIATLIGLDVFMIGTGMIAAFAATPGTRIAWWGISTGALLALLYVLVGTLSKDARGQSPEVASLFGRLRNLVIVLWLLYPVVWILGTEGTFGILPLYWETAAFMVLDLSAKVGFGVVLLRSRSVLR",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "4FBZ_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "deltarhodopsin; 7 transmembrane helices, light-driven proton; HET: RET, L2P, 22B, SQL, BNG, SO4; 2.7A {Haloterrigena thermotolerans} SCOP: f.13.1.0",
        )
        self.assertEqual(
            alignment[0],
            "GPESIWLWIGTIGMTLGTLYFVGRGRGVRDRKMQEFYIITTFITTIAAAMYFAMATGFGVTEVVVGDEALTIYWARYADWLFTTPLLLLDLGLLAGANRNTIATLIGLDVFMIGTGMIAAFAATPGTRIAWWGISTGALLALLYVLVGTLSKDARGQSPEVASLFGRLRNLVIVLWLLYPVVWILGTEGTFGILPLYWETAAFMVLDLSAKVGFGVVLLRSRSVLR",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "     CTTHHHHHHHHHHHHHHHHHHHHHTTTCCSTTHHHHHHHHHHHHHHHHHHHHHHHTTTTCEEEECSSCEEEECTHHHHHHHHHHHHHHHHHHHHHTCCHHHHHHHHHHHHHHHHHHHHHHHCSSHHHHHHHHHHHHHHHHHHHHHCCCCCHHHHTTSCHHHHHHHHHHHHHHHHHHHHHHHHHHHSGGGSSCCSCHHHHHHHHHHHHHCCCCCHHHHHSSCHHHHG          ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "     hhHHHHHHHHHHHHHHHHHHHHHHhcCCCChhHHHHHHHHHHHHHHHHHHHHHHHcCCCceeEecCCcceeeeHHHHHHHHHHHHHHHHHHHHHhCCCHHHHHHHHHHHHHHHHHHHHHHHcCChhhHHHHHHHHHHHHHHHHHHHHHHHHHHHhcCCHHHHHHHHHHHHHHHHHHHHHHHHHHHcccccCCCCChHHHHHHHHHHHHHHHHhHHHHHHHhhHHHH          ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  5, 194, 195, 231],
                             [  1, 190, 190, 226]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['G', 'P', 'E', 'S', 'I', 'W', 'L', 'W', 'I', 'G', 'T', 'I', 'G',
              'M', 'T', 'L', 'G', 'T', 'L', 'Y', 'F', 'V', 'G', 'R', 'G', 'R',
              'G', 'V', 'R', 'D', 'R', 'K', 'M', 'Q', 'E', 'F', 'Y', 'I', 'I',
              'T', 'T', 'F', 'I', 'T', 'T', 'I', 'A', 'A', 'A', 'M', 'Y', 'F',
              'A', 'M', 'A', 'T', 'G', 'F', 'G', 'V', 'T', 'E', 'V', 'V', 'V',
              'G', 'D', 'E', 'A', 'L', 'T', 'I', 'Y', 'W', 'A', 'R', 'Y', 'A',
              'D', 'W', 'L', 'F', 'T', 'T', 'P', 'L', 'L', 'L', 'L', 'D', 'L',
              'G', 'L', 'L', 'A', 'G', 'A', 'N', 'R', 'N', 'T', 'I', 'A', 'T',
              'L', 'I', 'G', 'L', 'D', 'V', 'F', 'M', 'I', 'G', 'T', 'G', 'M',
              'I', 'A', 'A', 'F', 'A', 'A', 'T', 'P', 'G', 'T', 'R', 'I', 'A',
              'W', 'W', 'G', 'I', 'S', 'T', 'G', 'A', 'L', 'L', 'A', 'L', 'L',
              'Y', 'V', 'L', 'V', 'G', 'T', 'L', 'S', 'K', 'D', 'A', 'R', 'G',
              'Q', 'S', 'P', 'E', 'V', 'A', 'S', 'L', 'F', 'G', 'R', 'L', 'R',
              'N', 'L', 'V', 'I', 'V', 'L', 'W', 'L', 'L', 'Y', 'P', 'V', 'V',
              'W', 'I', 'L', 'G', 'T', 'E', 'G', 'T', 'F', 'G', 'I', 'L', 'P',
              'L', 'Y', 'W', 'E', 'T', 'A', 'A', 'F', 'M', 'V', 'L', 'D', 'L',
              'S', 'A', 'K', 'V', 'G', 'F', 'G', 'V', 'V', 'L', 'L', 'R', 'S',
              'R', 'S', 'V', 'L', 'R'],
             ['R', 'P', 'E', 'W', 'I', 'W', 'L', 'A', 'L', 'G', 'T', 'A', 'L',
              'M', 'G', 'L', 'G', 'T', 'L', 'Y', 'F', 'L', 'V', 'K', 'G', 'M',
              'G', 'V', 'S', 'D', 'P', 'D', 'A', 'K', 'K', 'F', 'Y', 'A', 'I',
              'T', 'T', 'L', 'V', 'P', 'A', 'I', 'A', 'F', 'T', 'M', 'Y', 'L',
              'S', 'M', 'L', 'L', 'G', 'Y', 'G', 'L', 'T', 'M', 'V', 'P', 'F',
              'G', 'G', 'E', 'Q', 'N', 'P', 'I', 'Y', 'W', 'A', 'R', 'Y', 'A',
              'D', 'W', 'L', 'F', 'T', 'T', 'P', 'L', 'L', 'L', 'L', 'D', 'L',
              'A', 'L', 'L', 'V', 'D', 'A', 'D', 'Q', 'G', 'T', 'I', 'L', 'A',
              'L', 'V', 'G', 'A', 'D', 'G', 'I', 'M', 'I', 'G', 'T', 'G', 'L',
              'V', 'G', 'A', 'L', 'T', 'K', 'V', 'Y', 'S', 'Y', 'R', 'F', 'V',
              'W', 'W', 'A', 'I', 'S', 'T', 'A', 'A', 'M', 'L', 'Y', 'I', 'L',
              'Y', 'V', 'L', 'F', 'F', 'G', 'F', 'T', 'S', 'K', 'A', 'E', 'S',
              'M', 'R', 'P', 'E', 'V', 'A', 'S', 'T', 'F', 'K', 'V', 'L', 'R',
              'N', 'V', 'T', 'V', 'V', 'L', 'W', 'S', 'A', 'Y', 'P', 'V', 'V',
              'W', 'L', 'I', 'G', 'S', 'E', 'G', '-', 'A', 'G', 'I', 'V', 'P',
              'L', 'N', 'I', 'E', 'T', 'L', 'L', 'F', 'M', 'V', 'L', 'D', 'V',
              'S', 'A', 'K', 'V', 'G', 'F', 'G', 'L', 'I', 'L', 'L', 'R', 'S',
              'R', 'A', 'I', 'F', 'G']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 100.00)
        self.assertAlmostEqual(alignment.annotations["E-value"], 1e-42)
        self.assertAlmostEqual(alignment.annotations["Score"], 291.62)
        self.assertAlmostEqual(alignment.annotations["Identities"], 57)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.005)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 7.800)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "CCcHHHHHHHHHHHHHHHHHHHHHHhcCCCChhhHHHHHHHHHHHHHHHHHHHHHHhCCCceeeccCCCCCcchHHHHHHHHHHHHHHHHHHHHHcCCCHHHHHHHHHHHHHHHHHHHHHHHcccchHHHHHHHHHHHHHHHHHHHHHHhHHHHHHhCCHHHHHHHHHHHHHHHHHHHHHHHHHHHccCCCCCCChhHHHHHHHHHHHHHHHHHHHHHHHhhhhhC",
        )
        self.assertEqual(alignment.query.id, "4Y9H:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[0:226],
            "GRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFG",
        )
        self.assertEqual(
            alignment[1],
            "GRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVP-FGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFG",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "~~~~~~~~~~~~~~m~~~~~~f~~~~~~~~~~~~r~~~~l~~~i~~ia~~~Y~~ma~~~g~~~~~~~~~~~~v~~~RYv~W~ittPlll~~l~~l~~~~~~~~~~~v~~~~~mi~~g~~g~~~~~~~~~~~~~~~s~~~~~~i~~~l~~~~~~~~~~~~~~~~~~~~~l~~~~~i~W~~YPi~w~l~~~g~~~i~~~~~~~~~~ilDi~~K~~f~~~l~~~~~~~~",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "               ~~~~~~~~~~~~~~m~~~~~~f~~~~~~~~~~~~r~~~~l~~~i~~vaa~~Y~~ma~g~g~~~~~~~~~~~r~v~~~RYv~W~iTtPlll~~l~~la~~~~~~~~~~v~~~~~Mi~~G~~g~~~~~~~~~~~~~~~g~~~f~~v~~~l~~~~~~~a~~~~~~~~~~~~~l~~~~~v~W~~YPivw~l~~~g~~~i~~~~e~i~~~ilDv~~K~~f~~~l~~~~~~~~                 ",
        )
        self.assertEqual(alignment.target.id, "3WQJ_A")
        self.assertEqual(
            alignment.target.seq[15:242],
            "GRPETLWLGIGTLLMLIGTFYFIARGWGVTDKEAREYYAITILVPGIASAAYLAMFFGIGVTEVELASGTVLDIYYARYADWLFTTPLLLLDLALLAKVDRVTIGTLIGVDALMIVTGLIGALSKTPLARYTWWLFSTIAFLFVLYYLLTSLRSAAAKRSEEVRSTFNTLTALVAVLWTAYPILWIVGTEGAGVVGLGIETLAFMVLDVTAKVGFGFVLLRSRAILG",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "3WQJ_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Archaerhodopsin-2; 7 trans-membrane helices, light-driven proton; HET: RET, L3P, L2P, 22B, SQL, L4P; 1.8A {Halobacterium}; Related PDB entries: 2Z55_B 2EI4_A 2Z55_D 1VGO_A 2Z55_E 1VGO_B 2Z55_A 1UAZ_A 1UAZ_B",
        )
        self.assertEqual(
            alignment[0],
            "GRPETLWLGIGTLLMLIGTFYFIARGWGVTDKEAREYYAITILVPGIASAAYLAMFFGIGVTEVELASGTVLDIYYARYADWLFTTPLLLLDLALLAKVDRVTIGTLIGVDALMIVTGLIGALSKTPLARYTWWLFSTIAFLFVLYYLLTSLRSAAAKRSEEVRSTFNTLTALVAVLWTAYPILWIVGTEGAGVVGLGIETLAFMVLDVTAKVGFGFVLLRSRAILG",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "               SCTTHHHHHHHHHHHHHHHHHHHHHHTTCCCHHHHHHHHHHHHHHHHHHHHHHHHHHSTTCEEEECTTSCEEEECTHHHHHHHHHHHHHHHHHHHHHTCCHHHHHHHHHHHHHHHHHHHHHHHCSSHHHHHHHHHHHHHHHHHHHHHCCCCCHHHHTTSCHHHHHHHHHHHHHHHHHHHHHHHHHHHSTTTTCCSCHHHHHHHHHHHHHCCCCCHHHHHHTCTTC                   ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "               CchHHHHHHHHHHHHHHHHHHHHHHhcCCCCHHHHHHHHHHHHHHHHHHHHHHHHHcCCCceeEecCCCCeeeecHHHHHHHHHHHHHHHHHHHHHhCCCHHHHHHHHHHHHHHHHHHHHHHhCCCchhHHHHHHHHHHHHHHHHHHHHHHHHHHHHhcCHHHHHHHHHHHHHHHHHHHHHHHHHHHCcCCCCCCCHHHHHHHHHHHHHHHHHHHHHHHHHHHHHhc                 ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 15,  80,  81, 242],
                             [  0,  65,  65, 226]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['G', 'R', 'P', 'E', 'T', 'L', 'W', 'L', 'G', 'I', 'G', 'T', 'L',
              'L', 'M', 'L', 'I', 'G', 'T', 'F', 'Y', 'F', 'I', 'A', 'R', 'G',
              'W', 'G', 'V', 'T', 'D', 'K', 'E', 'A', 'R', 'E', 'Y', 'Y', 'A',
              'I', 'T', 'I', 'L', 'V', 'P', 'G', 'I', 'A', 'S', 'A', 'A', 'Y',
              'L', 'A', 'M', 'F', 'F', 'G', 'I', 'G', 'V', 'T', 'E', 'V', 'E',
              'L', 'A', 'S', 'G', 'T', 'V', 'L', 'D', 'I', 'Y', 'Y', 'A', 'R',
              'Y', 'A', 'D', 'W', 'L', 'F', 'T', 'T', 'P', 'L', 'L', 'L', 'L',
              'D', 'L', 'A', 'L', 'L', 'A', 'K', 'V', 'D', 'R', 'V', 'T', 'I',
              'G', 'T', 'L', 'I', 'G', 'V', 'D', 'A', 'L', 'M', 'I', 'V', 'T',
              'G', 'L', 'I', 'G', 'A', 'L', 'S', 'K', 'T', 'P', 'L', 'A', 'R',
              'Y', 'T', 'W', 'W', 'L', 'F', 'S', 'T', 'I', 'A', 'F', 'L', 'F',
              'V', 'L', 'Y', 'Y', 'L', 'L', 'T', 'S', 'L', 'R', 'S', 'A', 'A',
              'A', 'K', 'R', 'S', 'E', 'E', 'V', 'R', 'S', 'T', 'F', 'N', 'T',
              'L', 'T', 'A', 'L', 'V', 'A', 'V', 'L', 'W', 'T', 'A', 'Y', 'P',
              'I', 'L', 'W', 'I', 'V', 'G', 'T', 'E', 'G', 'A', 'G', 'V', 'V',
              'G', 'L', 'G', 'I', 'E', 'T', 'L', 'A', 'F', 'M', 'V', 'L', 'D',
              'V', 'T', 'A', 'K', 'V', 'G', 'F', 'G', 'F', 'V', 'L', 'L', 'R',
              'S', 'R', 'A', 'I', 'L', 'G'],
             ['G', 'R', 'P', 'E', 'W', 'I', 'W', 'L', 'A', 'L', 'G', 'T', 'A',
              'L', 'M', 'G', 'L', 'G', 'T', 'L', 'Y', 'F', 'L', 'V', 'K', 'G',
              'M', 'G', 'V', 'S', 'D', 'P', 'D', 'A', 'K', 'K', 'F', 'Y', 'A',
              'I', 'T', 'T', 'L', 'V', 'P', 'A', 'I', 'A', 'F', 'T', 'M', 'Y',
              'L', 'S', 'M', 'L', 'L', 'G', 'Y', 'G', 'L', 'T', 'M', 'V', 'P',
              '-', 'F', 'G', 'G', 'E', 'Q', 'N', 'P', 'I', 'Y', 'W', 'A', 'R',
              'Y', 'A', 'D', 'W', 'L', 'F', 'T', 'T', 'P', 'L', 'L', 'L', 'L',
              'D', 'L', 'A', 'L', 'L', 'V', 'D', 'A', 'D', 'Q', 'G', 'T', 'I',
              'L', 'A', 'L', 'V', 'G', 'A', 'D', 'G', 'I', 'M', 'I', 'G', 'T',
              'G', 'L', 'V', 'G', 'A', 'L', 'T', 'K', 'V', 'Y', 'S', 'Y', 'R',
              'F', 'V', 'W', 'W', 'A', 'I', 'S', 'T', 'A', 'A', 'M', 'L', 'Y',
              'I', 'L', 'Y', 'V', 'L', 'F', 'F', 'G', 'F', 'T', 'S', 'K', 'A',
              'E', 'S', 'M', 'R', 'P', 'E', 'V', 'A', 'S', 'T', 'F', 'K', 'V',
              'L', 'R', 'N', 'V', 'T', 'V', 'V', 'L', 'W', 'S', 'A', 'Y', 'P',
              'V', 'V', 'W', 'L', 'I', 'G', 'S', 'E', 'G', 'A', 'G', 'I', 'V',
              'P', 'L', 'N', 'I', 'E', 'T', 'L', 'L', 'F', 'M', 'V', 'L', 'D',
              'V', 'S', 'A', 'K', 'V', 'G', 'F', 'G', 'L', 'I', 'L', 'L', 'R',
              'S', 'R', 'A', 'I', 'F', 'G']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 100.00)
        self.assertAlmostEqual(alignment.annotations["E-value"], 1.3e-42)
        self.assertAlmostEqual(alignment.annotations["Score"], 289.73)
        self.assertAlmostEqual(alignment.annotations["Identities"], 55)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.964)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 8.100)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            " CcHHHHHHHHHHHHHHHHHHHHHHhcCCCChhhHHHHHHHHHHHHHHHHHHHHHHhCCCceeeccCCCCCcchHHHHHHHHHHHHHHHHHHHHHcCCCHHHHHHHHHHHHHHHHHHHHHHHcccchHHHHHHHHHHHHHHHHHHHHHHhHHHHHHhCCHHHHHHHHHHHHHHHHHHHHHHHHHHHccCCCCCCChhHHHHHHHHHHHHHHHHHHHHHHHhhhhh ",
        )
        self.assertEqual(alignment.query.id, "4Y9H:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[1:225],
            "RPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIF",
        )
        self.assertEqual(
            alignment[1],
            "RPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKV------YSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIF",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            " ~~~~~~~~~~~~~m~~~~~~f~~~~~~~~~~~~r~~~~l~~~i~~ia~~~Y~~ma~~~g~~~~~~~~~~~~v~~~RYv~W~ittPlll~~l~~l~~~~~~~~~~~v~~~~~mi~~g~~g~~~~~~~~~~~~~~~s~~~~~~i~~~l~~~~~~~~~~~~~~~~~~~~~l~~~~~i~W~~YPi~w~l~~~g~~~i~~~~~~~~~~ilDi~~K~~f~~~l~~~~~~~ ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "    ~~~~~~lw~~~~~~~~~~~~f~~~~~~~~~~~~r~~~~l~~~i~~va~~~Y~~ma~~~g~~~~~~~~~~r~v~~~RYv~W~ittPlll~~l~~lag~~~~~~~~~~~~~~~mi~~g~~~~~~~~~~~~~~~~~~w~~~~~s~~~~~~i~~~l~~~~~~~~~~~~~~~~~~~~~l~~~~~v~W~~YPi~w~l~~~g~~~i~~~~e~i~~~ilDv~~K~~f~~~ll~~~~~~                ",
        )
        self.assertEqual(alignment.target.id, "4L35_A")
        self.assertEqual(
            alignment.target.seq[4:234],
            "EGEAIWLWLGTAGMFLGMLYFIARGWGETDSRRQKFYIATILITAIAFVNYLAMALGFGLTIVEIAGEQRPIYWARYSDWLFTTPLLLYDLGLLAGADRNTISSLVSLDVLMIGTGLVATLSAGSGVLSAGAERLVWWGISTAFLLVLLYFLFSSLSGRVADLPSDTRSTFKTLRNLVTVVWLVYPVWWLVGTEGIGLVGIGIETAGFMVIDLVAKVGFGIILLRSHGVL",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "4L35_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Cruxrhodopsin-3; protein-bacterioruberin complex, seven transmembrane alpha; HET: RET, 22B; 2.1A {Haloarcula vallismortis} SCOP: f.13.1.0; Related PDB entries: 4PXK_A 4JR8_A",
        )
        self.assertEqual(
            alignment[0],
            "EGEAIWLWLGTAGMFLGMLYFIARGWGETDSRRQKFYIATILITAIAFVNYLAMALGFGLTIVEIAGEQRPIYWARYSDWLFTTPLLLYDLGLLAGADRNTISSLVSLDVLMIGTGLVATLSAGSGVLSAGAERLVWWGISTAFLLVLLYFLFSSLSGRVADLPSDTRSTFKTLRNLVTVVWLVYPVWWLVGTEGIGLVGIGIETAGFMVIDLVAKVGFGIILLRSHGVL",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "    CSSHHHHHHHHHHHHHHHHHHHHHHTTCCCHHHHHHHHHHHHHHHHHHHHHHHHHTTTTEEEEEETTEEEEEETHHHHHHHHHHHHHHHHHHHHHTCCHHHHHHHHHHHHHHHHHHHHHHHCCCCSSSCHHHHHHHHHHHHHHHHHHHHHHHHHTTHHHHHTSCHHHHHHHHHHHHHHHHHHHHHHHHHHHSTTTTCSSCHHHHHHHHHHHHHCCCCCHHHHHTSCHHHH                ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "    hhHHHHHHHHHHHHHHHHHHHHHHHcCCCChhHHHHHHHHHHHHHHHHHHHHHHHcCCCeeEeecCCcceeeeHHHHHHHHHHHHHHHHHHHHHhCCCHHHHHHHHHHHHHHHHHHHHHHHcccCCccCcchHHHHHHHHHHHHHHHHHHHHHhhhhHHHhcCChhHHHHHHHHHHHHHHHHHHHHHHHHHccccccccCHHHHHHHHHHHHHHHHHHHHHHHHHHcchh                ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  4, 128, 134, 234],
                             [  1, 125, 125, 225]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['E', 'G', 'E', 'A', 'I', 'W', 'L', 'W', 'L', 'G', 'T', 'A', 'G',
              'M', 'F', 'L', 'G', 'M', 'L', 'Y', 'F', 'I', 'A', 'R', 'G', 'W',
              'G', 'E', 'T', 'D', 'S', 'R', 'R', 'Q', 'K', 'F', 'Y', 'I', 'A',
              'T', 'I', 'L', 'I', 'T', 'A', 'I', 'A', 'F', 'V', 'N', 'Y', 'L',
              'A', 'M', 'A', 'L', 'G', 'F', 'G', 'L', 'T', 'I', 'V', 'E', 'I',
              'A', 'G', 'E', 'Q', 'R', 'P', 'I', 'Y', 'W', 'A', 'R', 'Y', 'S',
              'D', 'W', 'L', 'F', 'T', 'T', 'P', 'L', 'L', 'L', 'Y', 'D', 'L',
              'G', 'L', 'L', 'A', 'G', 'A', 'D', 'R', 'N', 'T', 'I', 'S', 'S',
              'L', 'V', 'S', 'L', 'D', 'V', 'L', 'M', 'I', 'G', 'T', 'G', 'L',
              'V', 'A', 'T', 'L', 'S', 'A', 'G', 'S', 'G', 'V', 'L', 'S', 'A',
              'G', 'A', 'E', 'R', 'L', 'V', 'W', 'W', 'G', 'I', 'S', 'T', 'A',
              'F', 'L', 'L', 'V', 'L', 'L', 'Y', 'F', 'L', 'F', 'S', 'S', 'L',
              'S', 'G', 'R', 'V', 'A', 'D', 'L', 'P', 'S', 'D', 'T', 'R', 'S',
              'T', 'F', 'K', 'T', 'L', 'R', 'N', 'L', 'V', 'T', 'V', 'V', 'W',
              'L', 'V', 'Y', 'P', 'V', 'W', 'W', 'L', 'V', 'G', 'T', 'E', 'G',
              'I', 'G', 'L', 'V', 'G', 'I', 'G', 'I', 'E', 'T', 'A', 'G', 'F',
              'M', 'V', 'I', 'D', 'L', 'V', 'A', 'K', 'V', 'G', 'F', 'G', 'I',
              'I', 'L', 'L', 'R', 'S', 'H', 'G', 'V', 'L'],
             ['R', 'P', 'E', 'W', 'I', 'W', 'L', 'A', 'L', 'G', 'T', 'A', 'L',
              'M', 'G', 'L', 'G', 'T', 'L', 'Y', 'F', 'L', 'V', 'K', 'G', 'M',
              'G', 'V', 'S', 'D', 'P', 'D', 'A', 'K', 'K', 'F', 'Y', 'A', 'I',
              'T', 'T', 'L', 'V', 'P', 'A', 'I', 'A', 'F', 'T', 'M', 'Y', 'L',
              'S', 'M', 'L', 'L', 'G', 'Y', 'G', 'L', 'T', 'M', 'V', 'P', 'F',
              'G', 'G', 'E', 'Q', 'N', 'P', 'I', 'Y', 'W', 'A', 'R', 'Y', 'A',
              'D', 'W', 'L', 'F', 'T', 'T', 'P', 'L', 'L', 'L', 'L', 'D', 'L',
              'A', 'L', 'L', 'V', 'D', 'A', 'D', 'Q', 'G', 'T', 'I', 'L', 'A',
              'L', 'V', 'G', 'A', 'D', 'G', 'I', 'M', 'I', 'G', 'T', 'G', 'L',
              'V', 'G', 'A', 'L', 'T', 'K', 'V', '-', '-', '-', '-', '-', '-',
              'Y', 'S', 'Y', 'R', 'F', 'V', 'W', 'W', 'A', 'I', 'S', 'T', 'A',
              'A', 'M', 'L', 'Y', 'I', 'L', 'Y', 'V', 'L', 'F', 'F', 'G', 'F',
              'T', 'S', 'K', 'A', 'E', 'S', 'M', 'R', 'P', 'E', 'V', 'A', 'S',
              'T', 'F', 'K', 'V', 'L', 'R', 'N', 'V', 'T', 'V', 'V', 'L', 'W',
              'S', 'A', 'Y', 'P', 'V', 'V', 'W', 'L', 'I', 'G', 'S', 'E', 'G',
              'A', 'G', 'I', 'V', 'P', 'L', 'N', 'I', 'E', 'T', 'L', 'L', 'F',
              'M', 'V', 'L', 'D', 'V', 'S', 'A', 'K', 'V', 'G', 'F', 'G', 'L',
              'I', 'L', 'L', 'R', 'S', 'R', 'A', 'I', 'F']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 100.00)
        self.assertAlmostEqual(alignment.annotations["E-value"], 2.2e-42)
        self.assertAlmostEqual(alignment.annotations["Score"], 289.87)
        self.assertAlmostEqual(alignment.annotations["Identities"], 56)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.948)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 7.600)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            " CcHHHHHHHHHHHHHHHHHHHHHHhcCCCChhhHHHHHHHHHHHHHHHHHHHHHHhCCCceeeccCCCCCcchHHHHHHHHHHHHHHHHHHHHHcCCCHHHHHHHHHHHHHHHHHHHHHHHcccchHHHHHHHHHHHHHHHHHHHHHHhHHHHHHhCCHHHHHHHHHHHHHHHHHHHHHHHHHHHccCCCCCCChhHHHHHHHHHHHHHHHHHHHHHHHhhhhhC",
        )
        self.assertEqual(alignment.query.id, "4Y9H:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[1:226],
            "RPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFG",
        )
        self.assertEqual(
            alignment[1],
            "RPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPF-GGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFG",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            " ~~~~~~~~~~~~~m~~~~~~f~~~~~~~~~~~~r~~~~l~~~i~~ia~~~Y~~ma~~~g~~~~~~~~~~~~v~~~RYv~W~ittPlll~~l~~l~~~~~~~~~~~v~~~~~mi~~g~~g~~~~~~~~~~~~~~~s~~~~~~i~~~l~~~~~~~~~~~~~~~~~~~~~l~~~~~i~W~~YPi~w~l~~~g~~~i~~~~~~~~~~ilDi~~K~~f~~~l~~~~~~~~",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "             ~~~~~~~~~~~~~m~~~~l~f~~~~~~~~~~~~r~~~~l~~~i~~vaa~~Y~~ma~g~g~~~~~~~~g~~r~v~~~RYv~W~iTtPlll~~l~llag~~~~~~~~~i~~~~~Mi~~G~~g~~~~~~~~~w~~~~is~~~~~~v~~~l~~~~~~~~~~~~~~~~~~~~~L~~~~~v~W~~YPivw~l~~~G~~~is~~~e~i~~~vlDil~K~~fg~~ll~~~~~~~                       ",
        )
        self.assertEqual(alignment.target.id, "4QI1_A")
        self.assertEqual(
            alignment.target.seq[13:239],
            "EGEGIWLALGTIGMLLGMLYFIADGLDVQDPRQKEFYVITILIPAIAAASYLSMFFGFGLTEVSLANGRVVDVYWARYADWLFTTPLLLLDIGLLAGASQRDIGALVGIDAFMIVTGLVATLTKVVVARYAFWTISTISMVFLLYYLVAVFGEAVSDADEDTRSTFNALRNIILVTWAIYPVAWLVGTEGLALTGLYGETLLFMVLDLVAKVGFGFILLRSRAIMG",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "4QI1_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Bacteriorhodopsin-I; Bacteriorhodopsin, proton pump, Membrane, MEMBRANE; HET: MPG, RET; 1.85A {Haloquadratum walsbyi} SCOP: f.13.1.0; Related PDB entries: 5ITC_A 5KKH_A 5ITE_B 5KKH_C 5ITE_C 5ITC_C 5ITC_B 5KKH_B 5ITE_A 4QI1_B 4WAV_B 4WAV_A 4QID_B 4QID_A 4QI1_C",
        )
        self.assertEqual(
            alignment[0],
            "EGEGIWLALGTIGMLLGMLYFIADGLDVQDPRQKEFYVITILIPAIAAASYLSMFFGFGLTEVSLANGRVVDVYWARYADWLFTTPLLLLDIGLLAGASQRDIGALVGIDAFMIVTGLVATLTKVVVARYAFWTISTISMVFLLYYLVAVFGEAVSDADEDTRSTFNALRNIILVTWAIYPVAWLVGTEGLALTGLYGETLLFMVLDLVAKVGFGFILLRSRAIMG",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "             SSTHHHHHHHHHHHHHHHHHHHHHHTTCCCHHHHHHHHHHHHHHHHHHHHHHHHHTTTTEEEEECTTSCEEEEETHHHHHHHHHHHHHHHHHHHHHTCCHHHHHHHHHHHHHHHHHHHHHHHCCCHHHHHHHHHHHHHHHHHHHHHCCCCCHHHHTTSCHHHHHHHHHHHHHHHHHHHHHHHHHHHSTTTTCSSCHHHHHHHHHHHHCCCCCCHHHHHHTCGGGC                        ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "             hhHHHHHHHHHHHHHHHHHHHHHHhcCCCChhHHHHHHHHHHHHHHHHHHHHHHHcCCCceeeeccCCCeeeeeHHHHHHHHHHHHHHHHHHHHHcCCCHHHHHHHHHHHHHHHHHHHHHHHccchhhHHHHHHHHHHHHHHHHHHHHHHHHHHHhcCCHHHHHHHHHHHHHHHHHHHHHHHHHHHCcccccccCHHHHHHHHHHHHHHHHHHHHHHHHHHHHHhC                       ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 13,  78,  79, 239],
                             [  1,  66,  66, 226]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['E', 'G', 'E', 'G', 'I', 'W', 'L', 'A', 'L', 'G', 'T', 'I', 'G',
              'M', 'L', 'L', 'G', 'M', 'L', 'Y', 'F', 'I', 'A', 'D', 'G', 'L',
              'D', 'V', 'Q', 'D', 'P', 'R', 'Q', 'K', 'E', 'F', 'Y', 'V', 'I',
              'T', 'I', 'L', 'I', 'P', 'A', 'I', 'A', 'A', 'A', 'S', 'Y', 'L',
              'S', 'M', 'F', 'F', 'G', 'F', 'G', 'L', 'T', 'E', 'V', 'S', 'L',
              'A', 'N', 'G', 'R', 'V', 'V', 'D', 'V', 'Y', 'W', 'A', 'R', 'Y',
              'A', 'D', 'W', 'L', 'F', 'T', 'T', 'P', 'L', 'L', 'L', 'L', 'D',
              'I', 'G', 'L', 'L', 'A', 'G', 'A', 'S', 'Q', 'R', 'D', 'I', 'G',
              'A', 'L', 'V', 'G', 'I', 'D', 'A', 'F', 'M', 'I', 'V', 'T', 'G',
              'L', 'V', 'A', 'T', 'L', 'T', 'K', 'V', 'V', 'V', 'A', 'R', 'Y',
              'A', 'F', 'W', 'T', 'I', 'S', 'T', 'I', 'S', 'M', 'V', 'F', 'L',
              'L', 'Y', 'Y', 'L', 'V', 'A', 'V', 'F', 'G', 'E', 'A', 'V', 'S',
              'D', 'A', 'D', 'E', 'D', 'T', 'R', 'S', 'T', 'F', 'N', 'A', 'L',
              'R', 'N', 'I', 'I', 'L', 'V', 'T', 'W', 'A', 'I', 'Y', 'P', 'V',
              'A', 'W', 'L', 'V', 'G', 'T', 'E', 'G', 'L', 'A', 'L', 'T', 'G',
              'L', 'Y', 'G', 'E', 'T', 'L', 'L', 'F', 'M', 'V', 'L', 'D', 'L',
              'V', 'A', 'K', 'V', 'G', 'F', 'G', 'F', 'I', 'L', 'L', 'R', 'S',
              'R', 'A', 'I', 'M', 'G'],
             ['R', 'P', 'E', 'W', 'I', 'W', 'L', 'A', 'L', 'G', 'T', 'A', 'L',
              'M', 'G', 'L', 'G', 'T', 'L', 'Y', 'F', 'L', 'V', 'K', 'G', 'M',
              'G', 'V', 'S', 'D', 'P', 'D', 'A', 'K', 'K', 'F', 'Y', 'A', 'I',
              'T', 'T', 'L', 'V', 'P', 'A', 'I', 'A', 'F', 'T', 'M', 'Y', 'L',
              'S', 'M', 'L', 'L', 'G', 'Y', 'G', 'L', 'T', 'M', 'V', 'P', 'F',
              '-', 'G', 'G', 'E', 'Q', 'N', 'P', 'I', 'Y', 'W', 'A', 'R', 'Y',
              'A', 'D', 'W', 'L', 'F', 'T', 'T', 'P', 'L', 'L', 'L', 'L', 'D',
              'L', 'A', 'L', 'L', 'V', 'D', 'A', 'D', 'Q', 'G', 'T', 'I', 'L',
              'A', 'L', 'V', 'G', 'A', 'D', 'G', 'I', 'M', 'I', 'G', 'T', 'G',
              'L', 'V', 'G', 'A', 'L', 'T', 'K', 'V', 'Y', 'S', 'Y', 'R', 'F',
              'V', 'W', 'W', 'A', 'I', 'S', 'T', 'A', 'A', 'M', 'L', 'Y', 'I',
              'L', 'Y', 'V', 'L', 'F', 'F', 'G', 'F', 'T', 'S', 'K', 'A', 'E',
              'S', 'M', 'R', 'P', 'E', 'V', 'A', 'S', 'T', 'F', 'K', 'V', 'L',
              'R', 'N', 'V', 'T', 'V', 'V', 'L', 'W', 'S', 'A', 'Y', 'P', 'V',
              'V', 'W', 'L', 'I', 'G', 'S', 'E', 'G', 'A', 'G', 'I', 'V', 'P',
              'L', 'N', 'I', 'E', 'T', 'L', 'L', 'F', 'M', 'V', 'L', 'D', 'V',
              'S', 'A', 'K', 'V', 'G', 'F', 'G', 'L', 'I', 'L', 'L', 'R', 'S',
              'R', 'A', 'I', 'F', 'G']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 100.00)
        self.assertAlmostEqual(alignment.annotations["E-value"], 4.4e-42)
        self.assertAlmostEqual(alignment.annotations["Score"], 282.32)
        self.assertAlmostEqual(alignment.annotations["Identities"], 29)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.533)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 8.300)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            " CcHHHHHHHHHHHHHHHHHHHHHHhcCCCChhhHHHHHHHHHHHHHHHHHHHHHHhCCCceeeccCCCCCcchHHHHHHHHHHHHHHHHHHHHHcCCCHHHHHHHHHHHHHHHHHHHHHHHcccchHHHHHHHHHHHHHHHHHHHHHHhHHHHHHhCCHHHHHHHHHHHHHHHHHHHHHHHHHHHccCCCCCCChhHHHHHHHHHHHHHHHHHHHHHHHhhhh  ",
        )
        self.assertEqual(alignment.query.id, "4Y9H:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[1:224],
            "RPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAI",
        )
        self.assertEqual(
            alignment[1],
            "RPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAI",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            " ~~~~~~~~~~~~~m~~~~~~f~~~~~~~~~~~~r~~~~l~~~i~~ia~~~Y~~ma~~~g~~~~~~~~~~~~v~~~RYv~W~ittPlll~~l~~l~~~~~~~~~~~v~~~~~mi~~g~~g~~~~~~~~~~~~~~~s~~~~~~i~~~l~~~~~~~~~~~~~~~~~~~~~l~~~~~i~W~~YPi~w~l~~~g~~~i~~~~~~~~~~ilDi~~K~~f~~~l~~~~~~  ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            " ~~~~~~~w~~~~~~~~~~~~f~~~~~~~~~~~r~~~~l~~~i~~ia~~~Y~~ma~~~g~~~~~~~~~~~~RY~~W~vTtPlll~~l~~lag~~~~~~~~~i~~~~~mi~~G~~~~~~~~~~~~~~~~i~~~~~~~i~~~l~~~~~~~a~~~~~~~~~~~~~l~~~~~v~W~~yPi~w~l~~~g~~~i~~~~e~~~~~~~Di~~K~~f~~~l~~~~~~       ",
        )
        self.assertEqual(alignment.target.id, "1H2S_A")
        self.assertEqual(
            alignment.target.seq[1:218],
            "VGLTTLFWLGAIGMLVGTLAFAWAGRDAGSGERRYYVTLVGISGIAAVAYVVMALGVGWVPVAERTVFAPRYIDWILTTPLIVYFLGLLAGLDSREFGIVITLNTVVMLAGFAGAMVPGIERYALFGMGAVAFLGLVYYLVGPMTESASQRSSGIKSLYVRLRNLTVILWAIYPFIWLLGPPGVALLTPTVDVALIVYLDLVTKVGFGFIALDAAAT",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "1H2S_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "SENSORY RHODOPSIN II, SENSORY RHODOPSIN; MEMBRANE PROTEIN, MENBRANE PROTEIN COMPLEX; HET: BOG, RET; 1.93A {NATRONOMONAS PHARAONIS} SCOP: f.13.1.1; Related PDB entries: 1JGJ_A",
        )
        self.assertEqual(
            alignment[0],
            "VGLTTLFWLGAIGMLVGTLAFAWAGRDAG-SGERRYYVTLVGISGIAAVAYVVMALGVGWVPV----AERTVFAPRYIDWILTTPLIVYFLGLLAGLDSREFGIVITLNTVVMLAGFAGAMVPG-IERYALFGMGAVAFLGLVYYLVGPMTESASQRSSGIKSLYVRLRNLTVILWAIYPFIWLLGPPGVALLTPTVDVALIVYLDLVTKVGFGFIALDAAAT",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            " CSHHHHHHHHHHHHHHHHHHHHHHTTTCCTTTHHHHHHHHHHHHHHHHHHHHHHTTCSEECCTTSCEEHHHHHHHHHHHHHHHHHHHHHHTCCHHHHHHHHHHHHHHHHHHHHHHHCCSTHHHHHHHHHHHHHHHHHHHHHTHHHHHHTTSCHHHHHHHHHHHHHHHHHHTTHHHHHHHSTTTTCCSCHHHHHHHHHHHHCCCCCCHHHHHHHHHHH       ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            " chHHHHHHHHHHHHHHHHHHHHHHHcCCCCCcHHHHHHHHHHHHHHHHHHHHHHcCCCccccCcccccHHHHHHHHHHHHHHHHHHHHHcCCChHHHHHHHHHHHHHHHHHHHHHhCCcchHHHHHHHHHHHHHHHHHHHHchHHHHHhhcchhHHHHHHHHHHHHHHHHHHHHHHHHHCccccCCCChhHHHHHHHHHHHHHHHHHHHHHHHHHHH       ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  1,  30,  30,  63,  63, 120, 120, 218],
                             [  1,  30,  31,  64,  68, 125, 126, 224]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['V', 'G', 'L', 'T', 'T', 'L', 'F', 'W', 'L', 'G', 'A', 'I', 'G',
              'M', 'L', 'V', 'G', 'T', 'L', 'A', 'F', 'A', 'W', 'A', 'G', 'R',
              'D', 'A', 'G', '-', 'S', 'G', 'E', 'R', 'R', 'Y', 'Y', 'V', 'T',
              'L', 'V', 'G', 'I', 'S', 'G', 'I', 'A', 'A', 'V', 'A', 'Y', 'V',
              'V', 'M', 'A', 'L', 'G', 'V', 'G', 'W', 'V', 'P', 'V', '-', '-',
              '-', '-', 'A', 'E', 'R', 'T', 'V', 'F', 'A', 'P', 'R', 'Y', 'I',
              'D', 'W', 'I', 'L', 'T', 'T', 'P', 'L', 'I', 'V', 'Y', 'F', 'L',
              'G', 'L', 'L', 'A', 'G', 'L', 'D', 'S', 'R', 'E', 'F', 'G', 'I',
              'V', 'I', 'T', 'L', 'N', 'T', 'V', 'V', 'M', 'L', 'A', 'G', 'F',
              'A', 'G', 'A', 'M', 'V', 'P', 'G', '-', 'I', 'E', 'R', 'Y', 'A',
              'L', 'F', 'G', 'M', 'G', 'A', 'V', 'A', 'F', 'L', 'G', 'L', 'V',
              'Y', 'Y', 'L', 'V', 'G', 'P', 'M', 'T', 'E', 'S', 'A', 'S', 'Q',
              'R', 'S', 'S', 'G', 'I', 'K', 'S', 'L', 'Y', 'V', 'R', 'L', 'R',
              'N', 'L', 'T', 'V', 'I', 'L', 'W', 'A', 'I', 'Y', 'P', 'F', 'I',
              'W', 'L', 'L', 'G', 'P', 'P', 'G', 'V', 'A', 'L', 'L', 'T', 'P',
              'T', 'V', 'D', 'V', 'A', 'L', 'I', 'V', 'Y', 'L', 'D', 'L', 'V',
              'T', 'K', 'V', 'G', 'F', 'G', 'F', 'I', 'A', 'L', 'D', 'A', 'A',
              'A', 'T'],
             ['R', 'P', 'E', 'W', 'I', 'W', 'L', 'A', 'L', 'G', 'T', 'A', 'L',
              'M', 'G', 'L', 'G', 'T', 'L', 'Y', 'F', 'L', 'V', 'K', 'G', 'M',
              'G', 'V', 'S', 'D', 'P', 'D', 'A', 'K', 'K', 'F', 'Y', 'A', 'I',
              'T', 'T', 'L', 'V', 'P', 'A', 'I', 'A', 'F', 'T', 'M', 'Y', 'L',
              'S', 'M', 'L', 'L', 'G', 'Y', 'G', 'L', 'T', 'M', 'V', 'P', 'F',
              'G', 'G', 'E', 'Q', 'N', 'P', 'I', 'Y', 'W', 'A', 'R', 'Y', 'A',
              'D', 'W', 'L', 'F', 'T', 'T', 'P', 'L', 'L', 'L', 'L', 'D', 'L',
              'A', 'L', 'L', 'V', 'D', 'A', 'D', 'Q', 'G', 'T', 'I', 'L', 'A',
              'L', 'V', 'G', 'A', 'D', 'G', 'I', 'M', 'I', 'G', 'T', 'G', 'L',
              'V', 'G', 'A', 'L', 'T', 'K', 'V', 'Y', 'S', 'Y', 'R', 'F', 'V',
              'W', 'W', 'A', 'I', 'S', 'T', 'A', 'A', 'M', 'L', 'Y', 'I', 'L',
              'Y', 'V', 'L', 'F', 'F', 'G', 'F', 'T', 'S', 'K', 'A', 'E', 'S',
              'M', 'R', 'P', 'E', 'V', 'A', 'S', 'T', 'F', 'K', 'V', 'L', 'R',
              'N', 'V', 'T', 'V', 'V', 'L', 'W', 'S', 'A', 'Y', 'P', 'V', 'V',
              'W', 'L', 'I', 'G', 'S', 'E', 'G', 'A', 'G', 'I', 'V', 'P', 'L',
              'N', 'I', 'E', 'T', 'L', 'L', 'F', 'M', 'V', 'L', 'D', 'V', 'S',
              'A', 'K', 'V', 'G', 'F', 'G', 'L', 'I', 'L', 'L', 'R', 'S', 'R',
              'A', 'I']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 100.00)
        self.assertAlmostEqual(alignment.annotations["E-value"], 4.6e-42)
        self.assertAlmostEqual(alignment.annotations["Score"], 282.87)
        self.assertAlmostEqual(alignment.annotations["Identities"], 29)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.443)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 8.400)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "CCcHHHHHHHHHHHHHHHHHHHHHHhcCCCChhhHHHHHHHHHHHHHHHHHHHHHHhCCCceeeccCCCCCcchHHHHHHHHHHHHHHHHHHHHHcCCCHHHHHHHHHHHHHHHHHHHHHHHcccchHHHHHHHHHHHHHHHHHHHHHHhHHHHHHhCCHHHHHHHHHHHHHHHHHHHHHHHHHHHccCCCCCCChhHHHHHHHHHHHHHHHHHHHHHHHhhhhh ",
        )
        self.assertEqual(alignment.query.id, "4Y9H:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[0:225],
            "GRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIF",
        )
        self.assertEqual(
            alignment[1],
            "GRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIF",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "~~~~~~~~~~~~~~m~~~~~~f~~~~~~~~~~~~r~~~~l~~~i~~ia~~~Y~~ma~~~g~~~~~~~~~~~~v~~~RYv~W~ittPlll~~l~~l~~~~~~~~~~~v~~~~~mi~~g~~g~~~~~~~~~~~~~~~s~~~~~~i~~~l~~~~~~~~~~~~~~~~~~~~~l~~~~~i~W~~YPi~w~l~~~g~~~i~~~~~~~~~~ilDi~~K~~f~~~l~~~~~~~ ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "     ~~~g~~~lw~~~~~~~~~~~~f~~~~~r~~~~~R~~~~l~~~i~~va~~~Y~~ma~~~~~~~~~~~r~v~~~RYv~W~vTtPlll~~l~~l~g~~~~~~~~~v~~~~~mi~~g~~~~~~~~~~~~w~~~~~s~~~~~~i~~~l~~~~~~~~~~~~~~~~~~~~~l~~~~~v~W~~YPivw~l~~g~~~i~~~~e~i~~~vlDv~~K~~f~~~ll~~~~~~    ",
        )
        self.assertEqual(alignment.target.id, "3AM6_A")
        self.assertEqual(
            alignment.target.seq[5:225],
            "TETGMIAQWIVFAIMAAAAIAFGVAVHFRPSELKSAYYINIAICTIAATAYYAMAVNYQDLTMNGERQVVYARYIDWVLTTPLLLLDLIVMTKMGGVMISWVIGADIFMIVFGILGAFEDEHKFKWVYFIAGCVMQAVLTYGMYNATWKDDLKKSPEYHSSYVSLLVFLSILWVFYPVVWAFGSGSGVLSVDNEAILMGILDVLAKPLFGMGCLIAHETI",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "3AM6_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "rhodopsin-2; seven trans-membrane, Transport Protein; HET: RET, CLR; 3.2A {Acetabularia acetabulum}; Related PDB entries: 3AM6_B 3AM6_C 3AM6_D",
        )
        self.assertEqual(
            alignment[0],
            "TETGMIAQWIVFAIMAAAAIAFGVAVHFRP-SELKSAYYINIAICTIAATAYYAMAVNYQDLTMN---GERQVVYARYIDWVLTTPLLLLDLIVMTKMGGVMISWVIGADIFMIVFGILGAFEDEHKFKWVYFIAGCVMQAVLTYGMYNATWKDDLKKSPEYHSSYVSLLVFLSILWVFYPVVWAFG-SGSGVLSVDNEAILMGILDVLAKPLFGMGCLIAHETI",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "     CHHHHHHHHHHHHHHHHHHHHHHHHGGGCCCHHHHHHHHHHHHHHHHHHHHHHHTTTSCCBSSSBBCCTHHHHHHHHHHHHHHHHHHTTSCCCHHHHHHHHHHHHHHHHHHHHHHTCCCHHHHHHHHHHHHHHHHHHHHHHHHHHHSCCCCHHHHHHHHHHHHHHHHHHTHHHHHHCCCCCCSSSCHHHHHHHHHHHHHCCCCCHHHHHHHHHHHH        ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "     cHHHHHHHHHHHHHHHHHHHHHHHHHcCCChhhHHHHHHHHHHHHHHHHHHHHHHhcccccccCCcccchHHHHHHHHHHHHHHHHHHHHHhcCCHHHHHHHHHHHHHHHHHHHHHHhCCccchHHHHHHHHHHHHHHHHHHHHHhhHHHHHhcCHHHHHHHHHHHHHHHHHHHHHHHHHHHhcCcCCCChHHHHHHHHHHHHHhhHHHHHHHHHHHHHH    ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  5,  35,  35,  69,  69, 188, 188, 225],
                             [  0,  30,  31,  65,  68, 187, 188, 225]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['T', 'E', 'T', 'G', 'M', 'I', 'A', 'Q', 'W', 'I', 'V', 'F', 'A',
              'I', 'M', 'A', 'A', 'A', 'A', 'I', 'A', 'F', 'G', 'V', 'A', 'V',
              'H', 'F', 'R', 'P', '-', 'S', 'E', 'L', 'K', 'S', 'A', 'Y', 'Y',
              'I', 'N', 'I', 'A', 'I', 'C', 'T', 'I', 'A', 'A', 'T', 'A', 'Y',
              'Y', 'A', 'M', 'A', 'V', 'N', 'Y', 'Q', 'D', 'L', 'T', 'M', 'N',
              '-', '-', '-', 'G', 'E', 'R', 'Q', 'V', 'V', 'Y', 'A', 'R', 'Y',
              'I', 'D', 'W', 'V', 'L', 'T', 'T', 'P', 'L', 'L', 'L', 'L', 'D',
              'L', 'I', 'V', 'M', 'T', 'K', 'M', 'G', 'G', 'V', 'M', 'I', 'S',
              'W', 'V', 'I', 'G', 'A', 'D', 'I', 'F', 'M', 'I', 'V', 'F', 'G',
              'I', 'L', 'G', 'A', 'F', 'E', 'D', 'E', 'H', 'K', 'F', 'K', 'W',
              'V', 'Y', 'F', 'I', 'A', 'G', 'C', 'V', 'M', 'Q', 'A', 'V', 'L',
              'T', 'Y', 'G', 'M', 'Y', 'N', 'A', 'T', 'W', 'K', 'D', 'D', 'L',
              'K', 'K', 'S', 'P', 'E', 'Y', 'H', 'S', 'S', 'Y', 'V', 'S', 'L',
              'L', 'V', 'F', 'L', 'S', 'I', 'L', 'W', 'V', 'F', 'Y', 'P', 'V',
              'V', 'W', 'A', 'F', 'G', '-', 'S', 'G', 'S', 'G', 'V', 'L', 'S',
              'V', 'D', 'N', 'E', 'A', 'I', 'L', 'M', 'G', 'I', 'L', 'D', 'V',
              'L', 'A', 'K', 'P', 'L', 'F', 'G', 'M', 'G', 'C', 'L', 'I', 'A',
              'H', 'E', 'T', 'I'],
             ['G', 'R', 'P', 'E', 'W', 'I', 'W', 'L', 'A', 'L', 'G', 'T', 'A',
              'L', 'M', 'G', 'L', 'G', 'T', 'L', 'Y', 'F', 'L', 'V', 'K', 'G',
              'M', 'G', 'V', 'S', 'D', 'P', 'D', 'A', 'K', 'K', 'F', 'Y', 'A',
              'I', 'T', 'T', 'L', 'V', 'P', 'A', 'I', 'A', 'F', 'T', 'M', 'Y',
              'L', 'S', 'M', 'L', 'L', 'G', 'Y', 'G', 'L', 'T', 'M', 'V', 'P',
              'F', 'G', 'G', 'E', 'Q', 'N', 'P', 'I', 'Y', 'W', 'A', 'R', 'Y',
              'A', 'D', 'W', 'L', 'F', 'T', 'T', 'P', 'L', 'L', 'L', 'L', 'D',
              'L', 'A', 'L', 'L', 'V', 'D', 'A', 'D', 'Q', 'G', 'T', 'I', 'L',
              'A', 'L', 'V', 'G', 'A', 'D', 'G', 'I', 'M', 'I', 'G', 'T', 'G',
              'L', 'V', 'G', 'A', 'L', 'T', 'K', 'V', 'Y', 'S', 'Y', 'R', 'F',
              'V', 'W', 'W', 'A', 'I', 'S', 'T', 'A', 'A', 'M', 'L', 'Y', 'I',
              'L', 'Y', 'V', 'L', 'F', 'F', 'G', 'F', 'T', 'S', 'K', 'A', 'E',
              'S', 'M', 'R', 'P', 'E', 'V', 'A', 'S', 'T', 'F', 'K', 'V', 'L',
              'R', 'N', 'V', 'T', 'V', 'V', 'L', 'W', 'S', 'A', 'Y', 'P', 'V',
              'V', 'W', 'L', 'I', 'G', 'S', 'E', 'G', 'A', 'G', 'I', 'V', 'P',
              'L', 'N', 'I', 'E', 'T', 'L', 'L', 'F', 'M', 'V', 'L', 'D', 'V',
              'S', 'A', 'K', 'V', 'G', 'F', 'G', 'L', 'I', 'L', 'L', 'R', 'S',
              'R', 'A', 'I', 'F']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 100.00)
        self.assertAlmostEqual(alignment.annotations["E-value"], 1.8e-41)
        self.assertAlmostEqual(alignment.annotations["Score"], 282.46)
        self.assertAlmostEqual(alignment.annotations["Identities"], 30)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.544)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 7.900)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "   HHHHHHHHHHHHHHHHHHHHHHhcCCCChhhHHHHHHHHHHHHHHHHHHHHHHhCCCceeeccCCCCCcchHHHHHHHHHHHHHHHHHHHHHcCCCHHHHHHHHHHHHHHHHHHHHHHHcccchHHHHHHHHHHHHHHHHHHHHHHhHHHHHHhCCHHHHHHHHHHHHHHHHHHHHHHHHHHHccCCCCCCChhHHHHHHHHHHHHHHHHHHHHHHHhhhh  ",
        )
        self.assertEqual(alignment.query.id, "4Y9H:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[3:224],
            "EWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAI",
        )
        self.assertEqual(
            alignment[1],
            "EWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAI",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "   ~~~~~~~~~~~m~~~~~~f~~~~~~~~~~~~r~~~~l~~~i~~ia~~~Y~~ma~~~g~~~~~~~~~~~~v~~~RYv~W~ittPlll~~l~~l~~~~~~~~~~~v~~~~~mi~~g~~g~~~~~~~~~~~~~~~s~~~~~~i~~~l~~~~~~~~~~~~~~~~~~~~~l~~~~~i~W~~YPi~w~l~~~g~~~i~~~~~~~~~~ilDi~~K~~f~~~l~~~~~~  ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "  ~~~~lw~~~~~m~~~~~~f~~~~~~~~~~~r~~~~l~~~I~~ia~~~Y~~ma~g~g~~~~~~r~v~~~RYv~W~~TtPlll~~l~~lag~~~~~~~~~i~~~~~mi~~G~~~~~~~~~~kw~~f~~~~~~f~~v~~~l~~~~~~~a~~~~~~~~~~~~~l~~~~~v~W~~YPi~w~l~~~g~~~i~~~~e~~~~~vlDi~~K~~fg~~ll~~~~~                               ",
        )
        self.assertEqual(alignment.target.id, "5JJF_A")
        self.assertEqual(
            alignment.target.seq[2:217],
            "LTTLFWLGAIGMLVGTLAFAWAGRDAGSGERRYYVTLVGISGIAAVAYVVMALGVGWVPVAERTVFAPRYIDWILTTPLIVYFLGLLAGLDSREFGIVITLNTVVMLAGFAGAMVPGIERYALFGMGAVAFLGLVYYLVGPMTESASQRSSGIKSLYVRLRNLTVILWAIYPFIWLLGPPGVALLTPTVDVALIVYLDLVTKVGFGFIALDAAAT",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "5JJF_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "sensory rhodopsin; sensory rhodopsin II, transducer, membrane; HET: LFA, RET, BOG; 1.9A {Natronomonas pharaonis}; Related PDB entries: 5JJE_A 5JJN_C 2F93_A 5JJN_A 2F95_A 5JJJ_A 2KSY_A 4GYC_A 1GU8_A 1GUE_A 3QDC_A 3QAP_A 1H68_A",
        )
        self.assertEqual(
            alignment[0],
            "LTTLFWLGAIGMLVGTLAFAWAGRDAG-SGERRYYVTLVGISGIAAVAYVVMALGVGWVPV----AERTVFAPRYIDWILTTPLIVYFLGLLAGLDSREFGIVITLNTVVMLAGFAGAMVPG-IERYALFGMGAVAFLGLVYYLVGPMTESASQRSSGIKSLYVRLRNLTVILWAIYPFIWLLGPPGVALLTPTVDVALIVYLDLVTKVGFGFIALDAAAT",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "  HHHHHHHHHHHHHHHHHHHHHHHcCCCCCcHHHHHHHHHHHHHHHHHHHHHHcCCCceecCcccccHHHHHHHHHHHHHHHHHHHHHcCCCHHHHHHHHHHHHHHHHHHHHHHhCCcchHHHHHHHHHHHHHHHHHHHHcHHHHHHhcCCHHHHHHHHHHHHHHHHHHHHHHHHHHHCCcccCCCCHHHHHHHHHHHHHHHHHHHHHHHHHHHHH                               ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  2,  29,  29,  62,  62, 119, 119, 217],
                             [  3,  30,  31,  64,  68, 125, 126, 224]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['L', 'T', 'T', 'L', 'F', 'W', 'L', 'G', 'A', 'I', 'G', 'M', 'L',
              'V', 'G', 'T', 'L', 'A', 'F', 'A', 'W', 'A', 'G', 'R', 'D', 'A',
              'G', '-', 'S', 'G', 'E', 'R', 'R', 'Y', 'Y', 'V', 'T', 'L', 'V',
              'G', 'I', 'S', 'G', 'I', 'A', 'A', 'V', 'A', 'Y', 'V', 'V', 'M',
              'A', 'L', 'G', 'V', 'G', 'W', 'V', 'P', 'V', '-', '-', '-', '-',
              'A', 'E', 'R', 'T', 'V', 'F', 'A', 'P', 'R', 'Y', 'I', 'D', 'W',
              'I', 'L', 'T', 'T', 'P', 'L', 'I', 'V', 'Y', 'F', 'L', 'G', 'L',
              'L', 'A', 'G', 'L', 'D', 'S', 'R', 'E', 'F', 'G', 'I', 'V', 'I',
              'T', 'L', 'N', 'T', 'V', 'V', 'M', 'L', 'A', 'G', 'F', 'A', 'G',
              'A', 'M', 'V', 'P', 'G', '-', 'I', 'E', 'R', 'Y', 'A', 'L', 'F',
              'G', 'M', 'G', 'A', 'V', 'A', 'F', 'L', 'G', 'L', 'V', 'Y', 'Y',
              'L', 'V', 'G', 'P', 'M', 'T', 'E', 'S', 'A', 'S', 'Q', 'R', 'S',
              'S', 'G', 'I', 'K', 'S', 'L', 'Y', 'V', 'R', 'L', 'R', 'N', 'L',
              'T', 'V', 'I', 'L', 'W', 'A', 'I', 'Y', 'P', 'F', 'I', 'W', 'L',
              'L', 'G', 'P', 'P', 'G', 'V', 'A', 'L', 'L', 'T', 'P', 'T', 'V',
              'D', 'V', 'A', 'L', 'I', 'V', 'Y', 'L', 'D', 'L', 'V', 'T', 'K',
              'V', 'G', 'F', 'G', 'F', 'I', 'A', 'L', 'D', 'A', 'A', 'A', 'T'],
             ['E', 'W', 'I', 'W', 'L', 'A', 'L', 'G', 'T', 'A', 'L', 'M', 'G',
              'L', 'G', 'T', 'L', 'Y', 'F', 'L', 'V', 'K', 'G', 'M', 'G', 'V',
              'S', 'D', 'P', 'D', 'A', 'K', 'K', 'F', 'Y', 'A', 'I', 'T', 'T',
              'L', 'V', 'P', 'A', 'I', 'A', 'F', 'T', 'M', 'Y', 'L', 'S', 'M',
              'L', 'L', 'G', 'Y', 'G', 'L', 'T', 'M', 'V', 'P', 'F', 'G', 'G',
              'E', 'Q', 'N', 'P', 'I', 'Y', 'W', 'A', 'R', 'Y', 'A', 'D', 'W',
              'L', 'F', 'T', 'T', 'P', 'L', 'L', 'L', 'L', 'D', 'L', 'A', 'L',
              'L', 'V', 'D', 'A', 'D', 'Q', 'G', 'T', 'I', 'L', 'A', 'L', 'V',
              'G', 'A', 'D', 'G', 'I', 'M', 'I', 'G', 'T', 'G', 'L', 'V', 'G',
              'A', 'L', 'T', 'K', 'V', 'Y', 'S', 'Y', 'R', 'F', 'V', 'W', 'W',
              'A', 'I', 'S', 'T', 'A', 'A', 'M', 'L', 'Y', 'I', 'L', 'Y', 'V',
              'L', 'F', 'F', 'G', 'F', 'T', 'S', 'K', 'A', 'E', 'S', 'M', 'R',
              'P', 'E', 'V', 'A', 'S', 'T', 'F', 'K', 'V', 'L', 'R', 'N', 'V',
              'T', 'V', 'V', 'L', 'W', 'S', 'A', 'Y', 'P', 'V', 'V', 'W', 'L',
              'I', 'G', 'S', 'E', 'G', 'A', 'G', 'I', 'V', 'P', 'L', 'N', 'I',
              'E', 'T', 'L', 'L', 'F', 'M', 'V', 'L', 'D', 'V', 'S', 'A', 'K',
              'V', 'G', 'F', 'G', 'L', 'I', 'L', 'L', 'R', 'S', 'R', 'A', 'I']],
            dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 100.00)
        self.assertAlmostEqual(alignment.annotations["E-value"], 5e-41)
        self.assertAlmostEqual(alignment.annotations["Score"], 279.16)
        self.assertAlmostEqual(alignment.annotations["Identities"], 27)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.447)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 8.000)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "      HHHHHHHHHHHHHHHHHHHhcCCCChhhHHHHHHHHHHHHHHHHHHHHHHhCCCceeeccCCCCCcchHHHHHHHHHHHHHHHHHHHHHcCCCHHHHHHHHHHHHHHHHHHHHHHHcccchHHHHHHHHHHHHHHHHHHHHHHhHHHHHHhCCHHHHHHHHHHHHHHHHHHHHHHHHHHHccCCCCCCChhHHHHHHHHHHHHHHHHHHHHHHHhhhhh ",
        )
        self.assertEqual(alignment.query.id, "4Y9H:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[6:225],
            "WLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIF",
        )
        self.assertEqual(
            alignment[1],
            "WLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTM--VPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIF",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "      ~~~~~~~~m~~~~~~f~~~~~~~~~~~~r~~~~l~~~i~~ia~~~Y~~ma~~~g~~~~~~~~~~~~v~~~RYv~W~ittPlll~~l~~l~~~~~~~~~~~v~~~~~mi~~g~~g~~~~~~~~~~~~~~~s~~~~~~i~~~l~~~~~~~~~~~~~~~~~~~~~l~~~~~i~W~~YPi~w~l~~~g~~~i~~~~~~~~~~ilDi~~K~~f~~~l~~~~~~~ ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                     ~lw~~~~i~~~~~l~f~~~~~~~~~~~R~~~~l~~~i~~ia~~~Y~~ma~~~g~~~~~~~~~~~~r~v~~~RYv~W~iTtPlll~~l~~la~~~~~~~~~~~~~~~~mi~~G~~g~~~~~~~~w~~~~is~~~~~~~~~~l~~~~~~~~~~~~~~~~~~~~~l~~~~~v~W~~YPi~w~l~~~g~~i~~~~e~i~~~vlDi~~K~~f~~~ll~~~~~~     ",
        )
        self.assertEqual(alignment.target.id, "5AX0_A")
        self.assertEqual(
            alignment.target.seq[21:239],
            "AQWVVFAVMALAAIVFSIAVQFRPLPLRLTYYVNIAICTIAATAYYAMAVNGGDNKPTAGTGADERQVIYARYIDWVFTTPLLLLDLVLLTNMPATMIAWIMGADIAMIAFGIIGAFTVGSYKWFYFVVGCIMLAVLAWGMINPIFKEELQKHKEYTGAYTTLLIYLIVLWVIYPIVWGLGAGGHIIGVDVEIIAMGILDLLAKPLYAIGVLITVEVV",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "5AX0_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Rhodopsin I; PROTON TRANSPORT, MEMBRANE PROTEIN, RETINAL; HET: RET, D10, OLB, R16, D12, C14; 1.521A {Acetabularia acetabulum}; Related PDB entries: 5AX1_A 5AWZ_A",
        )
        self.assertEqual(
            alignment[0],
            "AQWVVFAVMALAAIVFSIAVQFRP-LPLRLTYYVNIAICTIAATAYYAMAVNGGDNKPTAGTGADERQVIYARYIDWVFTTPLLLLDLVLLTNMPATMIAWIMGADIAMIAFGIIGAFTVG-SYKWFYFVVGCIMLAVLAWGMINPIFKEELQKHKEYTGAYTTLLIYLIVLWVIYPIVWGL-GAGGHIIGVDVEIIAMGILDLLAKPLYAIGVLITVEVV",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "                     HHHHHHHHHHHHHHHHHHHHTTSCHHHHHHHHHHHHHHHHHHHHHHHHHHHGGGCCCEESCGGGCEECCCHHHHHHHHHHHHHHHHHHTTBCCCHHHHHHHHHHHHHHHHHHHHHHHCCSTTHHHHHHHHHHHHHHHHHHHHHHHHCGGGBSCGGGHHHHHHHHHHHHHHHHHHHHHHCCCCCTCCSCHHHHHHHHHHHHHCCCCCHHHHHHHHHHHH     ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                     HHHHHHHHHHHHHHHHHHHhcCCCcchhhHHHHHHHHHHHHHHHHHHHHhcCCCCCCCCCCCCCceeeeHHHHHHHHHHHHHHHHHHHHHcCCCHHHHHHHHHHHHHHHHHHHHHHHcccchHHHHHHHHHHHHHHHHHHHHHHHHHHHHHhChhhHHHHHHHHHHHHHHHHHHHHHHHHcCccccCChhHHHHHHHHHHHHHHHHHHHHHHHHHHHH     ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 21,  45,  45,  77,  79, 141, 141, 201, 201, 239],
                             [  6,  30,  31,  63,  63, 125, 126, 186, 187, 225]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['A', 'Q', 'W', 'V', 'V', 'F', 'A', 'V', 'M', 'A', 'L', 'A', 'A',
              'I', 'V', 'F', 'S', 'I', 'A', 'V', 'Q', 'F', 'R', 'P', '-', 'L',
              'P', 'L', 'R', 'L', 'T', 'Y', 'Y', 'V', 'N', 'I', 'A', 'I', 'C',
              'T', 'I', 'A', 'A', 'T', 'A', 'Y', 'Y', 'A', 'M', 'A', 'V', 'N',
              'G', 'G', 'D', 'N', 'K', 'P', 'T', 'A', 'G', 'T', 'G', 'A', 'D',
              'E', 'R', 'Q', 'V', 'I', 'Y', 'A', 'R', 'Y', 'I', 'D', 'W', 'V',
              'F', 'T', 'T', 'P', 'L', 'L', 'L', 'L', 'D', 'L', 'V', 'L', 'L',
              'T', 'N', 'M', 'P', 'A', 'T', 'M', 'I', 'A', 'W', 'I', 'M', 'G',
              'A', 'D', 'I', 'A', 'M', 'I', 'A', 'F', 'G', 'I', 'I', 'G', 'A',
              'F', 'T', 'V', 'G', '-', 'S', 'Y', 'K', 'W', 'F', 'Y', 'F', 'V',
              'V', 'G', 'C', 'I', 'M', 'L', 'A', 'V', 'L', 'A', 'W', 'G', 'M',
              'I', 'N', 'P', 'I', 'F', 'K', 'E', 'E', 'L', 'Q', 'K', 'H', 'K',
              'E', 'Y', 'T', 'G', 'A', 'Y', 'T', 'T', 'L', 'L', 'I', 'Y', 'L',
              'I', 'V', 'L', 'W', 'V', 'I', 'Y', 'P', 'I', 'V', 'W', 'G', 'L',
              '-', 'G', 'A', 'G', 'G', 'H', 'I', 'I', 'G', 'V', 'D', 'V', 'E',
              'I', 'I', 'A', 'M', 'G', 'I', 'L', 'D', 'L', 'L', 'A', 'K', 'P',
              'L', 'Y', 'A', 'I', 'G', 'V', 'L', 'I', 'T', 'V', 'E', 'V', 'V'],
             ['W', 'L', 'A', 'L', 'G', 'T', 'A', 'L', 'M', 'G', 'L', 'G', 'T',
              'L', 'Y', 'F', 'L', 'V', 'K', 'G', 'M', 'G', 'V', 'S', 'D', 'P',
              'D', 'A', 'K', 'K', 'F', 'Y', 'A', 'I', 'T', 'T', 'L', 'V', 'P',
              'A', 'I', 'A', 'F', 'T', 'M', 'Y', 'L', 'S', 'M', 'L', 'L', 'G',
              'Y', 'G', 'L', 'T', 'M', '-', '-', 'V', 'P', 'F', 'G', 'G', 'E',
              'Q', 'N', 'P', 'I', 'Y', 'W', 'A', 'R', 'Y', 'A', 'D', 'W', 'L',
              'F', 'T', 'T', 'P', 'L', 'L', 'L', 'L', 'D', 'L', 'A', 'L', 'L',
              'V', 'D', 'A', 'D', 'Q', 'G', 'T', 'I', 'L', 'A', 'L', 'V', 'G',
              'A', 'D', 'G', 'I', 'M', 'I', 'G', 'T', 'G', 'L', 'V', 'G', 'A',
              'L', 'T', 'K', 'V', 'Y', 'S', 'Y', 'R', 'F', 'V', 'W', 'W', 'A',
              'I', 'S', 'T', 'A', 'A', 'M', 'L', 'Y', 'I', 'L', 'Y', 'V', 'L',
              'F', 'F', 'G', 'F', 'T', 'S', 'K', 'A', 'E', 'S', 'M', 'R', 'P',
              'E', 'V', 'A', 'S', 'T', 'F', 'K', 'V', 'L', 'R', 'N', 'V', 'T',
              'V', 'V', 'L', 'W', 'S', 'A', 'Y', 'P', 'V', 'V', 'W', 'L', 'I',
              'G', 'S', 'E', 'G', 'A', 'G', 'I', 'V', 'P', 'L', 'N', 'I', 'E',
              'T', 'L', 'L', 'F', 'M', 'V', 'L', 'D', 'V', 'S', 'A', 'K', 'V',
              'G', 'F', 'G', 'L', 'I', 'L', 'L', 'R', 'S', 'R', 'A', 'I', 'F']],
            dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 100.00)
        self.assertAlmostEqual(alignment.annotations["E-value"], 1.2e-40)
        self.assertAlmostEqual(alignment.annotations["Score"], 279.52)
        self.assertAlmostEqual(alignment.annotations["Identities"], 28)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.470)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 7.900)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            " CcHHHHHHHHHHHHHHHHHHHHHHhcCCCChhhHHHHHHHHHHHHHHHHHHHHHHhCCCceeeccCCCCCcchHHHHHHHHHHHHHHHHHHHHHcCCCHHHHHHHHHHHHHHHHHHHHHHHcccchHHHHHHHHHHHHHHHHHHHHHHhHHHHHHhCCHHHHHHHHHHHHHHHHHHHHHHHHHHHccCCCCCCChhHHHHHHHHHHHHHHHHHHHHHHHhhhhhC",
        )
        self.assertEqual(alignment.query.id, "4Y9H:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[1:226],
            "RPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFG",
        )
        self.assertEqual(
            alignment[1],
            "RPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDA----DQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFG",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            " ~~~~~~~~~~~~~m~~~~~~f~~~~~~~~~~~~r~~~~l~~~i~~ia~~~Y~~ma~~~g~~~~~~~~~~~~v~~~RYv~W~ittPlll~~l~~l~~~~~~~~~~~v~~~~~mi~~g~~g~~~~~~~~~~~~~~~s~~~~~~i~~~l~~~~~~~~~~~~~~~~~~~~~l~~~~~i~W~~YPi~w~l~~~g~~~i~~~~~~~~~~ilDi~~K~~f~~~l~~~~~~~~",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            " ~~~~~~~w~~~~~~~~~~~~f~~~~~~~~~~r~~~~~~~~~i~~ia~~~Y~~ma~g~g~~~~~g~~v~~~RYv~W~vTtPlll~~l~~la~~~~~~~~~~~~~~i~~~~~mi~~G~~a~~~~~~~~kw~~~~~s~~~f~~il~~l~~~~~~~a~~~~~~~~~~~~~L~~~~~v~W~~YPivw~l~~~g~~~is~~~~~i~~~vlDi~~K~~f~~~ll~~~~~~~                                    ",
        )
        self.assertEqual(alignment.target.id, "1XIO_A")
        self.assertEqual(
            alignment.target.seq[1:225],
            "NLESLLHWIYVAGMTIGALHFWSLSRNPRGVPQYEYLVAMFIPIWSGLAYMAMAIDQGKVEAAGQIAHYARYIDWMVTTPLLLLSLSWTAMQFIKKDWTLIGFLMSTQIVVITSGLIADLSERDWVRYLWYICGVCAFLIILWGIWNPLRAKTRTQSSELANLYDKLVTYFTVLWIGYPIVWIIGPSGFGWINQTIDTFLFCLLPFFSKVGFSFLDLHGLRNLN",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "1XIO_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "ANABAENA SENSORY RHODOPSIN; signaling protein, photoreceptor; HET: PEE, RET; 2.0A {Nostoc sp. PCC 7120} SCOP: f.13.1.1; Related PDB entries: 2M3G_C 5UK6_B 5UK6_A 5UK6_C 2M3G_B 2M3G_A 4TL3_B 4TL3_A",
        )
        self.assertEqual(
            alignment[0],
            "NLESLLHWIYVAGMTIGALHFWSLSRN-PRGVPQYEYLVAMFIPIWSGLAYMAMAIDQGKVEA----AGQIAHYARYIDWMVTTPLLLLSLSWTAMQFIKKDWTLIGFLMSTQIVVITSGLIADLSERDWVRYLWYICGVCAFLIILWGIWNPLRAKTRTQSSELANLYDKLVTYFTVLWIGYPIVWIIGPSGFGWINQTIDTFLFCLLPFFSKVGFSFLDLHGLRNLN",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            " CHHHHHHHHHHHHHHHHHHHHHHHHTSCTTCCHHHHHHHHHHHHHHHHHHHHHHHCCTTHHHHHHHHHHHHHHHHHHHHHHTTTSCCCHHHHHHHHHHHHHHHHHHHHHHHCSSHHHHHHHHHHHHHHHHHHHHHHHTHHHHHHTTSCHHHHHHHHHHHHHHHHHHHHHHHHHHHSTTTTCSSCHHHHHHHHHHHHHCCCCCHHHHHHHHHHHTT                                             ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            " chHHHHHHHHHHHHHHHHHHHHHHhCCCCCCchHHHHHHHHHHHHHHHHHHHHHhccCceecccccccHHHHHHHHHHHHHHHHHHHHHHHhcccCCHHHHHHHHHHHHHHHHHHHHHHhcCCcchHHHHHHHHHHHHHHHHHHHHHHHHHHHHhcCHHHHHHHHHHHHHHHHHHHHHHHHHHHCCCCCCCCCHHHHHHHHHHHHHHHHHHHHHHHHHHHhhcc                                    ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  1,  28,  28,  63,  63,  93,  97, 225],
                             [  1,  28,  29,  64,  68,  98,  98, 226]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['N', 'L', 'E', 'S', 'L', 'L', 'H', 'W', 'I', 'Y', 'V', 'A', 'G',
              'M', 'T', 'I', 'G', 'A', 'L', 'H', 'F', 'W', 'S', 'L', 'S', 'R',
              'N', '-', 'P', 'R', 'G', 'V', 'P', 'Q', 'Y', 'E', 'Y', 'L', 'V',
              'A', 'M', 'F', 'I', 'P', 'I', 'W', 'S', 'G', 'L', 'A', 'Y', 'M',
              'A', 'M', 'A', 'I', 'D', 'Q', 'G', 'K', 'V', 'E', 'A', '-', '-',
              '-', '-', 'A', 'G', 'Q', 'I', 'A', 'H', 'Y', 'A', 'R', 'Y', 'I',
              'D', 'W', 'M', 'V', 'T', 'T', 'P', 'L', 'L', 'L', 'L', 'S', 'L',
              'S', 'W', 'T', 'A', 'M', 'Q', 'F', 'I', 'K', 'K', 'D', 'W', 'T',
              'L', 'I', 'G', 'F', 'L', 'M', 'S', 'T', 'Q', 'I', 'V', 'V', 'I',
              'T', 'S', 'G', 'L', 'I', 'A', 'D', 'L', 'S', 'E', 'R', 'D', 'W',
              'V', 'R', 'Y', 'L', 'W', 'Y', 'I', 'C', 'G', 'V', 'C', 'A', 'F',
              'L', 'I', 'I', 'L', 'W', 'G', 'I', 'W', 'N', 'P', 'L', 'R', 'A',
              'K', 'T', 'R', 'T', 'Q', 'S', 'S', 'E', 'L', 'A', 'N', 'L', 'Y',
              'D', 'K', 'L', 'V', 'T', 'Y', 'F', 'T', 'V', 'L', 'W', 'I', 'G',
              'Y', 'P', 'I', 'V', 'W', 'I', 'I', 'G', 'P', 'S', 'G', 'F', 'G',
              'W', 'I', 'N', 'Q', 'T', 'I', 'D', 'T', 'F', 'L', 'F', 'C', 'L',
              'L', 'P', 'F', 'F', 'S', 'K', 'V', 'G', 'F', 'S', 'F', 'L', 'D',
              'L', 'H', 'G', 'L', 'R', 'N', 'L', 'N'],
             ['R', 'P', 'E', 'W', 'I', 'W', 'L', 'A', 'L', 'G', 'T', 'A', 'L',
              'M', 'G', 'L', 'G', 'T', 'L', 'Y', 'F', 'L', 'V', 'K', 'G', 'M',
              'G', 'V', 'S', 'D', 'P', 'D', 'A', 'K', 'K', 'F', 'Y', 'A', 'I',
              'T', 'T', 'L', 'V', 'P', 'A', 'I', 'A', 'F', 'T', 'M', 'Y', 'L',
              'S', 'M', 'L', 'L', 'G', 'Y', 'G', 'L', 'T', 'M', 'V', 'P', 'F',
              'G', 'G', 'E', 'Q', 'N', 'P', 'I', 'Y', 'W', 'A', 'R', 'Y', 'A',
              'D', 'W', 'L', 'F', 'T', 'T', 'P', 'L', 'L', 'L', 'L', 'D', 'L',
              'A', 'L', 'L', 'V', 'D', 'A', '-', '-', '-', '-', 'D', 'Q', 'G',
              'T', 'I', 'L', 'A', 'L', 'V', 'G', 'A', 'D', 'G', 'I', 'M', 'I',
              'G', 'T', 'G', 'L', 'V', 'G', 'A', 'L', 'T', 'K', 'V', 'Y', 'S',
              'Y', 'R', 'F', 'V', 'W', 'W', 'A', 'I', 'S', 'T', 'A', 'A', 'M',
              'L', 'Y', 'I', 'L', 'Y', 'V', 'L', 'F', 'F', 'G', 'F', 'T', 'S',
              'K', 'A', 'E', 'S', 'M', 'R', 'P', 'E', 'V', 'A', 'S', 'T', 'F',
              'K', 'V', 'L', 'R', 'N', 'V', 'T', 'V', 'V', 'L', 'W', 'S', 'A',
              'Y', 'P', 'V', 'V', 'W', 'L', 'I', 'G', 'S', 'E', 'G', 'A', 'G',
              'I', 'V', 'P', 'L', 'N', 'I', 'E', 'T', 'L', 'L', 'F', 'M', 'V',
              'L', 'D', 'V', 'S', 'A', 'K', 'V', 'G', 'F', 'G', 'L', 'I', 'L',
              'L', 'R', 'S', 'R', 'A', 'I', 'F', 'G']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 100.00)
        self.assertAlmostEqual(alignment.annotations["E-value"], 1.4e-40)
        self.assertAlmostEqual(alignment.annotations["Score"], 282.56)
        self.assertAlmostEqual(alignment.annotations["Identities"], 32)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.598)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 7.300)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "      HHHHHHHHHHHHHHHHHHHhcCCCChhhHHHHHHHHHHHHHHHHHHHHHHhCCCceeeccCCCCCcchHHHHHHHHHHHHHHHHHHHHHcCCCHHHHHHHHHHHHHHHHHHHHHHHcccchHHHHHHHHHHHHHHHHHHHHHHhHHHHHHhCCHHHHHHHHHHHHHHHHHHHHHHHHHHHccCCCCCCChhHHHHHHHHHHHHHHHHHHHHHHHhhhh  ",
        )
        self.assertEqual(alignment.query.id, "4Y9H:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[6:224],
            "WLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAI",
        )
        self.assertEqual(
            alignment[1],
            "WLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVP----------------FGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVY-SYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAI",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "      ~~~~~~~~m~~~~~~f~~~~~~~~~~~~r~~~~l~~~i~~ia~~~Y~~ma~~~g~~~~~~~~~~~~v~~~RYv~W~ittPlll~~l~~l~~~~~~~~~~~v~~~~~mi~~g~~g~~~~~~~~~~~~~~~s~~~~~~i~~~l~~~~~~~~~~~~~~~~~~~~~l~~~~~i~W~~YPi~w~l~~~g~~~i~~~~~~~~~~ilDi~~K~~f~~~l~~~~~~  ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                                    ~lw~~~~i~~l~~~vf~~~s~~~~~~~~R~~~~l~~~i~~va~~aY~~Ma~g~G~~~v~~~~~~~~~~~~~~~g~~~~~~r~v~~~RYv~W~iTtPlll~~l~llag~~~~~i~~~v~~~~~mi~~Gl~gal~~~~~~~kw~~f~~g~~~~~~i~~~l~~~~~~~a~~~~~~~~~~~l~~~~~v~W~~YPi~w~l~~~G~~vi~~~~e~i~y~ilDi~~K~~f~~~ll~~~~~                      ",
        )
        self.assertEqual(alignment.target.id, "5B0W_G")
        self.assertEqual(
            alignment.target.seq[36:269],
            "SLYINIALAGLSILLFVFMTRGLDDPRAKLIAVSTILVPVVSIASYTGLASGLTISVLEMPAGHFAEGSSVMLGGEEVDGVVTMWGRYLTWALSTPMILLALGLLAGSNATKLFTAITFDIAMCVTGLAAALTTSSHLMRWFWYAISCACFIVVLYILLVEWAQDAKAAGTADIFSTLKLLTVVMWLGYPIVWALGVEGVAVLPVGYTSWAYSALDIVAKYIFAFLLLNYLTS",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "5B0W_G")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Halorhodopsin; seven trans-membrane helices, retinylidene protein; HET: RET, L3P, BNG, 22B; 1.7A {Natronomonas pharaonis DSM 2160}; Related PDB entries: 5B0W_D 3QBI_D 5ETZ_A 3ABW_D 3A7K_B 5ETZ_D 3QBK_D 3ABW_A 3VVK_C 4QRY_G 5B0W_F 3QBL_A 3ABW_B 3QBG_B 4QRY_E 3QBK_A 3QBI_B 3QBG_D 3VVK_B 3QBG_A 4QRY_C 5B0W_B 4QRY_F 3QBL_D 5B0W_A 3VVK_D 4QRY_A 5B0W_E 3VVK_A 3QBL_B 5ETZ_B 3A7K_A 3QBI_A 4QRY_B 3QBK_B 3VVK_F 3A7K_D 3VVK_E",
        )
        self.assertEqual(
            alignment[0],
            "SLYINIALAGLSILLFVFMTRGLDDPRAKLIAVSTILVPVVSIASYTGLASGLTISVLEMPAGHFAEGSSVMLGGEEVDGVVTMWGRYLTWALSTPMILLALGLLAGSNATKLFTAITFDIAMCVTGLAAALTTSSHLMRWFWYAISCACFIVVLYILLVEWAQDAKAAG--TADIFSTLKLLTVVMWLGYPIVWALGVEGVAVLPVGYTSWAYSALDIVAKYIFAFLLLNYLTS",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "                                    HHHHHHHHHHHHHHHHHHHHTTCCCHHHHHHHHHHHHHHHHHHHHHHHHHHTTTCEEECCCTTCTTTTCBCCBTTBCCBSEEECTHHHHHHHHHHHHHHHHHHHHHTCCHHHHHHHHHHHHHHHHHHHHHHHCCSCHHHHHHHHHHHHHHHHHHHHCCCCCHHHHHHHHTCHHHHHHHHHHHHHHHHHHHHHHHHBTTTTCSBCHHHHHHHHHHHHCCCCCCHHHHHHHHHTT                      ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                                    HHHHHHHHHHHHHHHHHHHHcCCCChHHHHHHHHHHHHHHHHHHHHHHHHcCCCceEEecccccccCCceeeecCCccCceeccHHHHHHHHHHHHHHHHHHHHHhCCCHHHHHHHHHHHHHHHHHHHHHHHcCCCcccHHHHHHHHHHHHHHHHHHHHHHHHHHHHHhChHHHHHHHHHHHHHHHHHHHHHHHHccCCCCCCChhHHHHHHHHHHHHHHHHHHHHHHHHHHh                      ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 36,  95, 111, 172, 173, 206, 206, 269],
                             [  6,  65,  65, 126, 126, 159, 161, 224]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['S', 'L', 'Y', 'I', 'N', 'I', 'A', 'L', 'A', 'G', 'L', 'S', 'I',
              'L', 'L', 'F', 'V', 'F', 'M', 'T', 'R', 'G', 'L', 'D', 'D', 'P',
              'R', 'A', 'K', 'L', 'I', 'A', 'V', 'S', 'T', 'I', 'L', 'V', 'P',
              'V', 'V', 'S', 'I', 'A', 'S', 'Y', 'T', 'G', 'L', 'A', 'S', 'G',
              'L', 'T', 'I', 'S', 'V', 'L', 'E', 'M', 'P', 'A', 'G', 'H', 'F',
              'A', 'E', 'G', 'S', 'S', 'V', 'M', 'L', 'G', 'G', 'E', 'E', 'V',
              'D', 'G', 'V', 'V', 'T', 'M', 'W', 'G', 'R', 'Y', 'L', 'T', 'W',
              'A', 'L', 'S', 'T', 'P', 'M', 'I', 'L', 'L', 'A', 'L', 'G', 'L',
              'L', 'A', 'G', 'S', 'N', 'A', 'T', 'K', 'L', 'F', 'T', 'A', 'I',
              'T', 'F', 'D', 'I', 'A', 'M', 'C', 'V', 'T', 'G', 'L', 'A', 'A',
              'A', 'L', 'T', 'T', 'S', 'S', 'H', 'L', 'M', 'R', 'W', 'F', 'W',
              'Y', 'A', 'I', 'S', 'C', 'A', 'C', 'F', 'I', 'V', 'V', 'L', 'Y',
              'I', 'L', 'L', 'V', 'E', 'W', 'A', 'Q', 'D', 'A', 'K', 'A', 'A',
              'G', '-', '-', 'T', 'A', 'D', 'I', 'F', 'S', 'T', 'L', 'K', 'L',
              'L', 'T', 'V', 'V', 'M', 'W', 'L', 'G', 'Y', 'P', 'I', 'V', 'W',
              'A', 'L', 'G', 'V', 'E', 'G', 'V', 'A', 'V', 'L', 'P', 'V', 'G',
              'Y', 'T', 'S', 'W', 'A', 'Y', 'S', 'A', 'L', 'D', 'I', 'V', 'A',
              'K', 'Y', 'I', 'F', 'A', 'F', 'L', 'L', 'L', 'N', 'Y', 'L', 'T',
              'S'],
             ['W', 'L', 'A', 'L', 'G', 'T', 'A', 'L', 'M', 'G', 'L', 'G', 'T',
              'L', 'Y', 'F', 'L', 'V', 'K', 'G', 'M', 'G', 'V', 'S', 'D', 'P',
              'D', 'A', 'K', 'K', 'F', 'Y', 'A', 'I', 'T', 'T', 'L', 'V', 'P',
              'A', 'I', 'A', 'F', 'T', 'M', 'Y', 'L', 'S', 'M', 'L', 'L', 'G',
              'Y', 'G', 'L', 'T', 'M', 'V', 'P', '-', '-', '-', '-', '-', '-',
              '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', 'F', 'G', 'G',
              'E', 'Q', 'N', 'P', 'I', 'Y', 'W', 'A', 'R', 'Y', 'A', 'D', 'W',
              'L', 'F', 'T', 'T', 'P', 'L', 'L', 'L', 'L', 'D', 'L', 'A', 'L',
              'L', 'V', 'D', 'A', 'D', 'Q', 'G', 'T', 'I', 'L', 'A', 'L', 'V',
              'G', 'A', 'D', 'G', 'I', 'M', 'I', 'G', 'T', 'G', 'L', 'V', 'G',
              'A', 'L', 'T', 'K', 'V', 'Y', '-', 'S', 'Y', 'R', 'F', 'V', 'W',
              'W', 'A', 'I', 'S', 'T', 'A', 'A', 'M', 'L', 'Y', 'I', 'L', 'Y',
              'V', 'L', 'F', 'F', 'G', 'F', 'T', 'S', 'K', 'A', 'E', 'S', 'M',
              'R', 'P', 'E', 'V', 'A', 'S', 'T', 'F', 'K', 'V', 'L', 'R', 'N',
              'V', 'T', 'V', 'V', 'L', 'W', 'S', 'A', 'Y', 'P', 'V', 'V', 'W',
              'L', 'I', 'G', 'S', 'E', 'G', 'A', 'G', 'I', 'V', 'P', 'L', 'N',
              'I', 'E', 'T', 'L', 'L', 'F', 'M', 'V', 'L', 'D', 'V', 'S', 'A',
              'K', 'V', 'G', 'F', 'G', 'L', 'I', 'L', 'L', 'R', 'S', 'R', 'A',
              'I']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 100.00)
        self.assertAlmostEqual(alignment.annotations["E-value"], 2.2e-40)
        self.assertAlmostEqual(alignment.annotations["Score"], 279.46)
        self.assertAlmostEqual(alignment.annotations["Identities"], 34)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.606)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 7.600)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "      HHHHHHHHHHHHHHHHHHHhcCCCChhhHHHHHHHHHHHHHHHHHHHHHHhCCCceeeccCCCCCcchHHHHHHHHHHHHHHHHHHHHHcCCCHHHHHHHHHHHHHHHHHHHHHHHcccchHHHHHHHHHHHHHHHHHHHHHHhHHHHHHhCCHHHHHHHHHHHHHHHHHHHHHHHHHHHccCCCCCCChhHHHHHHHHHHHHHHHHHHHHHHHhhhhhC",
        )
        self.assertEqual(alignment.query.id, "4Y9H:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[6:226],
            "WLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFG",
        )
        self.assertEqual(
            alignment[1],
            "WLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPF------GGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYS--YRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIV-PLNIETLLFMVLDVSAKVGFGLILLRSRAIFG",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "      ~~~~~~~~m~~~~~~f~~~~~~~~~~~~r~~~~l~~~i~~ia~~~Y~~ma~~~g~~~~~~~~~~~~v~~~RYv~W~ittPlll~~l~~l~~~~~~~~~~~v~~~~~mi~~g~~g~~~~~~~~~~~~~~~s~~~~~~i~~~l~~~~~~~~~~~~~~~~~~~~~l~~~~~i~W~~YPi~w~l~~~g~~~i~~~~~~~~~~ilDi~~K~~f~~~l~~~~~~~~",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                               ~lw~~~~~~~~~~l~f~~~~~~~~~~~~R~~~~l~~~i~~iaa~~Y~~ma~g~g~~~v~~~~~~~~~~~~~~~~~~RYvdW~vTtPlll~~l~~lag~~~~~~~~~v~~~~~mi~~g~~g~~~~~~~~~kw~~~~is~~~~~~i~~~l~~~~~~~a~~~~~~~~~~~l~~~~~v~W~~YPivw~l~~~G~~~i~s~~~e~i~y~ilDv~~K~~fg~~ll~~~~~~~                 ",
        )
        self.assertEqual(alignment.target.id, "2JAF_A")
        self.assertEqual(
            alignment.target.seq[31:257],
            "SLWVNVALAGIAILVFVYMGRTIRPGRPRLIWGATLMIPLVSISSYLGLLSGLTVGMIEMPAGHALAGEMVRSQWGRYLTWALSTPMILLALGLLADVDLGSLFTVIAADIGMCVTGLAAAMTTSALLFRWAFYAISCAFFVVVLSALVTDWAASASSAGTAEIFDTLRVLVVVLWLGYPIVWAVGVEGLALVQSVGATSWAYSVLDVFAKYVFAFILLRWVANNE",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "2JAF_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "HALORHODOPSIN; CHROMOPHORE, CHLORIDE PUMP, ION TRANSPORT; HET: BOG, RET, PLM; 1.7A {HALOBACTERIUM SALINARIUM} SCOP: f.13.1.1; Related PDB entries: 2JAG_A 5AHZ_A 5G36_A 5AHY_A 1E12_A",
        )
        self.assertEqual(
            alignment[0],
            "SLWVNVALAGIAILVFVYMGRTIRPGRPRLIWGATLMIPLVSISSYLGLLSGLTVGMIEMPAGHALAGEMVRSQWGRYLTWALSTPMILLALGLLADVDLGSLFTVIAADIGMCVTGLAAAMTTS-ALLFRWAFYAISCAFFVVVLSALVTDWAASAS--SAGTAEIFDTLRVLVVVLWLGYPIVWAVGVEGLALVQSVGATSWAYSVLDVFAKYVFAFILLRWVANNE",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "                               HHHHHHHHHHHHHHHHHHHTTTCCTTHHHHHHHHHHHHHHHHHHHHHHHHHTTTCEEEECCTTSTTTTSEEEECHHHHHHHHHHHHHHHHHHHHHHTCCHHHHHHHHHHHHHHHHHHHHHHHCCSCHHHHHHHHHHHHHHHHHHHHHHHTHHHHHHHHHTCHHHHHHHHHHHHHHHHHHHHHHHHSTTTTCSSCCHHHHHHHHHHHHCCCCCCHHHHHHHHHHHTH                 ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                               HHHHHHHHHHHHHHHHHHHHcCCCCCCccHHHHHHHHHHHHHHHHHHHHHcCCCceeeecccCCcccCCeEeeeHHHHHHHHhHHHHHHHHHHHHcCCCHHHHHHHHHHHHHHHHHHHHHHHcCCchhhHHHHHHHHHHHHHHHHHHHHhHHHHhcCChhHHHHHHHHHHHHHHHHHHHHHHHHHccccccccCChHHHHHHHHHHHHHHHHHHHHHHHHHHHhch                 ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 31,  91,  97, 156, 156, 157, 159, 188, 188, 224, 225, 257],
                             [  6,  66,  66, 125, 126, 127, 127, 156, 158, 194, 194, 226]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['S', 'L', 'W', 'V', 'N', 'V', 'A', 'L', 'A', 'G', 'I', 'A', 'I',
              'L', 'V', 'F', 'V', 'Y', 'M', 'G', 'R', 'T', 'I', 'R', 'P', 'G',
              'R', 'P', 'R', 'L', 'I', 'W', 'G', 'A', 'T', 'L', 'M', 'I', 'P',
              'L', 'V', 'S', 'I', 'S', 'S', 'Y', 'L', 'G', 'L', 'L', 'S', 'G',
              'L', 'T', 'V', 'G', 'M', 'I', 'E', 'M', 'P', 'A', 'G', 'H', 'A',
              'L', 'A', 'G', 'E', 'M', 'V', 'R', 'S', 'Q', 'W', 'G', 'R', 'Y',
              'L', 'T', 'W', 'A', 'L', 'S', 'T', 'P', 'M', 'I', 'L', 'L', 'A',
              'L', 'G', 'L', 'L', 'A', 'D', 'V', 'D', 'L', 'G', 'S', 'L', 'F',
              'T', 'V', 'I', 'A', 'A', 'D', 'I', 'G', 'M', 'C', 'V', 'T', 'G',
              'L', 'A', 'A', 'A', 'M', 'T', 'T', 'S', '-', 'A', 'L', 'L', 'F',
              'R', 'W', 'A', 'F', 'Y', 'A', 'I', 'S', 'C', 'A', 'F', 'F', 'V',
              'V', 'V', 'L', 'S', 'A', 'L', 'V', 'T', 'D', 'W', 'A', 'A', 'S',
              'A', 'S', '-', '-', 'S', 'A', 'G', 'T', 'A', 'E', 'I', 'F', 'D',
              'T', 'L', 'R', 'V', 'L', 'V', 'V', 'V', 'L', 'W', 'L', 'G', 'Y',
              'P', 'I', 'V', 'W', 'A', 'V', 'G', 'V', 'E', 'G', 'L', 'A', 'L',
              'V', 'Q', 'S', 'V', 'G', 'A', 'T', 'S', 'W', 'A', 'Y', 'S', 'V',
              'L', 'D', 'V', 'F', 'A', 'K', 'Y', 'V', 'F', 'A', 'F', 'I', 'L',
              'L', 'R', 'W', 'V', 'A', 'N', 'N', 'E'],
             ['W', 'L', 'A', 'L', 'G', 'T', 'A', 'L', 'M', 'G', 'L', 'G', 'T',
              'L', 'Y', 'F', 'L', 'V', 'K', 'G', 'M', 'G', 'V', 'S', 'D', 'P',
              'D', 'A', 'K', 'K', 'F', 'Y', 'A', 'I', 'T', 'T', 'L', 'V', 'P',
              'A', 'I', 'A', 'F', 'T', 'M', 'Y', 'L', 'S', 'M', 'L', 'L', 'G',
              'Y', 'G', 'L', 'T', 'M', 'V', 'P', 'F', '-', '-', '-', '-', '-',
              '-', 'G', 'G', 'E', 'Q', 'N', 'P', 'I', 'Y', 'W', 'A', 'R', 'Y',
              'A', 'D', 'W', 'L', 'F', 'T', 'T', 'P', 'L', 'L', 'L', 'L', 'D',
              'L', 'A', 'L', 'L', 'V', 'D', 'A', 'D', 'Q', 'G', 'T', 'I', 'L',
              'A', 'L', 'V', 'G', 'A', 'D', 'G', 'I', 'M', 'I', 'G', 'T', 'G',
              'L', 'V', 'G', 'A', 'L', 'T', 'K', 'V', 'Y', 'S', '-', '-', 'Y',
              'R', 'F', 'V', 'W', 'W', 'A', 'I', 'S', 'T', 'A', 'A', 'M', 'L',
              'Y', 'I', 'L', 'Y', 'V', 'L', 'F', 'F', 'G', 'F', 'T', 'S', 'K',
              'A', 'E', 'S', 'M', 'R', 'P', 'E', 'V', 'A', 'S', 'T', 'F', 'K',
              'V', 'L', 'R', 'N', 'V', 'T', 'V', 'V', 'L', 'W', 'S', 'A', 'Y',
              'P', 'V', 'V', 'W', 'L', 'I', 'G', 'S', 'E', 'G', 'A', 'G', 'I',
              'V', '-', 'P', 'L', 'N', 'I', 'E', 'T', 'L', 'L', 'F', 'M', 'V',
              'L', 'D', 'V', 'S', 'A', 'K', 'V', 'G', 'F', 'G', 'L', 'I', 'L',
              'L', 'R', 'S', 'R', 'A', 'I', 'F', 'G']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 100.00)
        self.assertAlmostEqual(alignment.annotations["E-value"], 2.6e-40)
        self.assertAlmostEqual(alignment.annotations["Score"], 266.44)
        self.assertAlmostEqual(alignment.annotations["Identities"], 99)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.599)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 8.100)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "CCcHHHHHHHHHHHHHHHHHHHHHHhcCCCChhhHHHHHHHHHHHHHHHHHHHHHHhCCCceeeccCCCCCcchHHHHHHHHHHHHHHHHHHHHHcCCCHHHHHHHHHHHHHHHHHHHHHHHcccchHHHHHHHHHHHHHHHHHHHHHHhHHHHHHhCCHHHHHHHHHHHHHHHHHHHHHHHHHHHccCCCCCCChhHHHHHHHHHHHHHHHHHHHH         ",
        )
        self.assertEqual(alignment.query.id, "4Y9H:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[0:217],
            "GRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLI",
        )
        self.assertEqual(
            alignment[1],
            "GRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLI",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "~~~~~~~~~~~~~~m~~~~~~f~~~~~~~~~~~~r~~~~l~~~i~~ia~~~Y~~ma~~~g~~~~~~~~~~~~v~~~RYv~W~ittPlll~~l~~l~~~~~~~~~~~v~~~~~mi~~g~~g~~~~~~~~~~~~~~~s~~~~~~i~~~l~~~~~~~~~~~~~~~~~~~~~l~~~~~i~W~~YPi~w~l~~~g~~~i~~~~~~~~~~ilDi~~K~~f~~~         ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            " ~~~~~~~~~~~~~~~~~~~~~f~~~~~~~~~~~~r~~~~l~~~i~~va~~~Y~~ma~~~g~~~~~~~~~~~~~~~~RYv~W~iTtPlll~~l~~l~g~~~~~~~~~v~~~~~mi~~g~~g~~~~~~~~~~~~~~~~~~~~~~~~~~l~~~~~v~W~~YPi~w~l~~~g~~~i~~~~~~~~~~~~Di~~K~~fg~~",
        )
        self.assertEqual(alignment.target.id, "1C8S_A")
        self.assertEqual(
            alignment.target.seq[1:196],
            "GRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLNLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLI",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "1C8S_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "BACTERIORHODOPSIN MUTANT/LIPID COMPLEX; ION PUMP, MEMBRANE PROTEIN, RETINAL; HET: LI1, SQU, RET; 2.0A {Halobacterium salinarum} SCOP: f.13.1.1",
        )
        self.assertEqual(
            alignment[0],
            "GRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLNLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLF----------------------NVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLI",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            " CCTTHHHHHHHHHHHHHHHHHHHHHHTTCCSHHHHHHHHHHHHHHHHHHHHHHHHHTTSSCBCCEETTEECCBCHHHHHHHHHHHHHHHHHHHHHTTCCHHHHHHHHHHHHHHHHHHHHHHHCSSHHHHHHHHHHHHHHHHHHHHHCCCHHHHHHHHHHHHHHHSTTTTCCSCHHHHHHHHHHHHHCCCCCHHHC",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            " CCcHHHHHHHHHHHHHHHHHHHHHHHhCCCChHHHHHHHHHHHHHHHHHHHHHHHHhCCCceeccCCCCCCcchhHHHHHHHHHHHHHHHHHHHHhCCCHHHHHHHHHHHHHHHHHHHHHHhccchhHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHccCCCCCCChhHHHHHHHHHHHHhHHhHhcC",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  1, 149, 149, 196],
                             [  0, 148, 170, 217]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['G', 'R', 'P', 'E', 'W', 'I', 'W', 'L', 'A', 'L', 'G', 'T', 'A',
              'L', 'M', 'G', 'L', 'G', 'T', 'L', 'Y', 'F', 'L', 'V', 'K', 'G',
              'M', 'G', 'V', 'S', 'D', 'P', 'D', 'A', 'K', 'K', 'F', 'Y', 'A',
              'I', 'T', 'T', 'L', 'V', 'P', 'A', 'I', 'A', 'F', 'T', 'M', 'Y',
              'L', 'S', 'M', 'L', 'L', 'G', 'Y', 'G', 'L', 'T', 'M', 'V', 'P',
              'F', 'G', 'G', 'E', 'Q', 'N', 'P', 'I', 'Y', 'W', 'A', 'R', 'Y',
              'A', 'D', 'W', 'L', 'F', 'T', 'T', 'P', 'L', 'L', 'L', 'L', 'N',
              'L', 'A', 'L', 'L', 'V', 'D', 'A', 'D', 'Q', 'G', 'T', 'I', 'L',
              'A', 'L', 'V', 'G', 'A', 'D', 'G', 'I', 'M', 'I', 'G', 'T', 'G',
              'L', 'V', 'G', 'A', 'L', 'T', 'K', 'V', 'Y', 'S', 'Y', 'R', 'F',
              'V', 'W', 'W', 'A', 'I', 'S', 'T', 'A', 'A', 'M', 'L', 'Y', 'I',
              'L', 'Y', 'V', 'L', 'F', '-', '-', '-', '-', '-', '-', '-', '-',
              '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',
              '-', 'N', 'V', 'T', 'V', 'V', 'L', 'W', 'S', 'A', 'Y', 'P', 'V',
              'V', 'W', 'L', 'I', 'G', 'S', 'E', 'G', 'A', 'G', 'I', 'V', 'P',
              'L', 'N', 'I', 'E', 'T', 'L', 'L', 'F', 'M', 'V', 'L', 'D', 'V',
              'S', 'A', 'K', 'V', 'G', 'F', 'G', 'L', 'I'],
             ['G', 'R', 'P', 'E', 'W', 'I', 'W', 'L', 'A', 'L', 'G', 'T', 'A',
              'L', 'M', 'G', 'L', 'G', 'T', 'L', 'Y', 'F', 'L', 'V', 'K', 'G',
              'M', 'G', 'V', 'S', 'D', 'P', 'D', 'A', 'K', 'K', 'F', 'Y', 'A',
              'I', 'T', 'T', 'L', 'V', 'P', 'A', 'I', 'A', 'F', 'T', 'M', 'Y',
              'L', 'S', 'M', 'L', 'L', 'G', 'Y', 'G', 'L', 'T', 'M', 'V', 'P',
              'F', 'G', 'G', 'E', 'Q', 'N', 'P', 'I', 'Y', 'W', 'A', 'R', 'Y',
              'A', 'D', 'W', 'L', 'F', 'T', 'T', 'P', 'L', 'L', 'L', 'L', 'D',
              'L', 'A', 'L', 'L', 'V', 'D', 'A', 'D', 'Q', 'G', 'T', 'I', 'L',
              'A', 'L', 'V', 'G', 'A', 'D', 'G', 'I', 'M', 'I', 'G', 'T', 'G',
              'L', 'V', 'G', 'A', 'L', 'T', 'K', 'V', 'Y', 'S', 'Y', 'R', 'F',
              'V', 'W', 'W', 'A', 'I', 'S', 'T', 'A', 'A', 'M', 'L', 'Y', 'I',
              'L', 'Y', 'V', 'L', 'F', 'F', 'G', 'F', 'T', 'S', 'K', 'A', 'E',
              'S', 'M', 'R', 'P', 'E', 'V', 'A', 'S', 'T', 'F', 'K', 'V', 'L',
              'R', 'N', 'V', 'T', 'V', 'V', 'L', 'W', 'S', 'A', 'Y', 'P', 'V',
              'V', 'W', 'L', 'I', 'G', 'S', 'E', 'G', 'A', 'G', 'I', 'V', 'P',
              'L', 'N', 'I', 'E', 'T', 'L', 'L', 'F', 'M', 'V', 'L', 'D', 'V',
              'S', 'A', 'K', 'V', 'G', 'F', 'G', 'L', 'I']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 100.00)
        self.assertAlmostEqual(alignment.annotations["E-value"], 4.7e-39)
        self.assertAlmostEqual(alignment.annotations["Score"], 273.19)
        self.assertAlmostEqual(alignment.annotations["Identities"], 22)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.327)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 7.600)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "CCcHHHHHHHHHHHHHHHHHHHHHHhcCCCChhhHHHHHHHHHHHHHHHHHHHHHHhCCCceeeccCCCCCcchHHHHHHHHHHHHHHHHHHHHHcCCCHHHHHHHHHHHHHHHHHHHHHHHcccchHHHHHHHHHHHHHHHHHHHHHHhHHHHHHhCCHHHHHHHHHHHHHHHHHHHHHHHHHHHccCCCCCCChhHHHHHHHHHHHHHHHHHHHHHHHhhhh  ",
        )
        self.assertEqual(alignment.query.id, "4Y9H:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[0:224],
            "GRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAI",
        )
        self.assertEqual(
            alignment[1],
            "GRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFG-------GEQNPI--YWARYADWLFTTPLLLLDLALLVDADQGTIL----ALVGADGIMIGTGLVGALTKV--YSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSE----GAGIV--PLNIETLLFMVLDVSAKVGFGLILLRSRAI",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "~~~~~~~~~~~~~~m~~~~~~f~~~~~~~~~~~~r~~~~l~~~i~~ia~~~Y~~ma~~~g~~~~~~~~~~~~v~~~RYv~W~ittPlll~~l~~l~~~~~~~~~~~v~~~~~mi~~g~~g~~~~~~~~~~~~~~~s~~~~~~i~~~l~~~~~~~~~~~~~~~~~~~~~l~~~~~i~W~~YPi~w~l~~~g~~~i~~~~~~~~~~ilDi~~K~~f~~~l~~~~~~  ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                        ~~~~s~~lw~~~~~~~~~~l~f~~~~~~~~~~~R~~~~l~~~i~~vaa~~Y~~ma~~~g~~~~~~~~~~~~~~~~~~~~~~~~~RYvdW~vTtPlll~~l~llag~~~~~~~~~~~~~v~~~~~Mi~~G~~g~~~~~~~~~~kw~~f~ig~~~~~~vl~~l~~~~~~~a~~~~~~~~~~~~~L~~~~~v~W~~YPivw~l~~~~~~~G~~li~~~~~~~~i~y~ilDv~aK~~fg~~ll~~~~~                    ",
        )
        self.assertEqual(alignment.target.id, "4XTL_A")
        self.assertEqual(
            alignment.target.seq[24:268],
            "YQFTSHILTLGYAVMLAGLLYFILTIKNVDKKFQMSNILSAVVMVSAFLLLYAQAQNWTSSFTFNEEVGRYFLDPSGDLFNNGYRYLNWLIDVPMLLFQILFVVSLTTSKFSSVRNQFWFSGAMMIITGYIGQFYEVSNLTAFLVWGAISSAFFFHILWVMKKVINEGKEGISPAGQKILSNIWILFLISWTLYPGAYLMPYLTGVDGFLYSEDGVMARQLVYTIADVSSXVIYGVLLGNLAIT",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "4XTL_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "KR2; ion pump, membrane protein; HET: GOL, LFA, LYR; 1.45A {Dokdonia eikasta}; Related PDB entries: 3X3B_A 3X3C_A 5JRF_A 4XTO_D 4XTN_C 4XTN_H 4XTN_A 4XTO_E 4XTN_B 4XTN_D 4XTO_A 4XTN_G 4XTO_B 4XTN_F 4XTO_C 4XTN_I 4XTN_J 4XTN_E",
        )
        self.assertEqual(
            alignment[0],
            "YQFTSHILTLGYAVMLAGLLYFILTIKNVD-KKFQMSNILSAVVMVSAFLLLYAQAQNWTSSFTFNEEVGRYFLDPSGDLFNNGYRYLNWLIDVPMLLFQILFVVSLTTSKFSSVRNQFWFSGAMMIITGYIGQFYEVSNLTAFLVWGAISSAFFFHILWVMKKVINEGKEGISPAGQKILSNIWILFLISWTLYPGAYLMPYLTGVDGFLYSEDGVMARQLVYTIADVSSXVIYGVLLGNLAIT",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "                        HHHHHHHHHHHHHHHHHHHHHHHHTGGGSCGGGHHHHHHHHHHHHHHHHHHHHHHHHHHHHEEEETTTTEEEECTTSCCCCTHHHHHHHHHHHHHHHHHGGGTSCCSSSCHHHHHHHHHHHHHHHHHHHHHHHHTTTTCHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHTTSCHHHHHHHHHHHHHHHHHHHHHHHHHHGGGTTCTTSTTSSHHHHHHHHHHHHHHHCCCCCCHHHHHHHHHHH                    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                        hHHHHHHHHHHHHHHHHHHHHHHHHHcCCCcchhHHHHHHHHHHHHHHHHHHHHHhchhhcccccchhcccccCCCCCccchHHHHHHHHHHHHHHHHHHHHHHcccccchHHHHHHHHHHHHHHHHHHHHHHHhccCcchhHHHHHHHHHHHHHHHHHHHHHHHHHhhcCCCHHHHHHHHHHHHHHHHHHHHHHHHHHHHHhcCCCCcccCCCCHHHHHHHHHHHHHHhHHHHHHHHHHHHHH                    ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 24,  54,  54,  90,  97, 103, 105, 136, 140, 161, 163, 227, 231,
                              236, 238, 268],
                             [  0,  30,  31,  67,  67,  73,  73, 104, 104, 125, 125, 189, 189,
                              194, 194, 224]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['Y', 'Q', 'F', 'T', 'S', 'H', 'I', 'L', 'T', 'L', 'G', 'Y', 'A',
              'V', 'M', 'L', 'A', 'G', 'L', 'L', 'Y', 'F', 'I', 'L', 'T', 'I',
              'K', 'N', 'V', 'D', '-', 'K', 'K', 'F', 'Q', 'M', 'S', 'N', 'I',
              'L', 'S', 'A', 'V', 'V', 'M', 'V', 'S', 'A', 'F', 'L', 'L', 'L',
              'Y', 'A', 'Q', 'A', 'Q', 'N', 'W', 'T', 'S', 'S', 'F', 'T', 'F',
              'N', 'E', 'E', 'V', 'G', 'R', 'Y', 'F', 'L', 'D', 'P', 'S', 'G',
              'D', 'L', 'F', 'N', 'N', 'G', 'Y', 'R', 'Y', 'L', 'N', 'W', 'L',
              'I', 'D', 'V', 'P', 'M', 'L', 'L', 'F', 'Q', 'I', 'L', 'F', 'V',
              'V', 'S', 'L', 'T', 'T', 'S', 'K', 'F', 'S', 'S', 'V', 'R', 'N',
              'Q', 'F', 'W', 'F', 'S', 'G', 'A', 'M', 'M', 'I', 'I', 'T', 'G',
              'Y', 'I', 'G', 'Q', 'F', 'Y', 'E', 'V', 'S', 'N', 'L', 'T', 'A',
              'F', 'L', 'V', 'W', 'G', 'A', 'I', 'S', 'S', 'A', 'F', 'F', 'F',
              'H', 'I', 'L', 'W', 'V', 'M', 'K', 'K', 'V', 'I', 'N', 'E', 'G',
              'K', 'E', 'G', 'I', 'S', 'P', 'A', 'G', 'Q', 'K', 'I', 'L', 'S',
              'N', 'I', 'W', 'I', 'L', 'F', 'L', 'I', 'S', 'W', 'T', 'L', 'Y',
              'P', 'G', 'A', 'Y', 'L', 'M', 'P', 'Y', 'L', 'T', 'G', 'V', 'D',
              'G', 'F', 'L', 'Y', 'S', 'E', 'D', 'G', 'V', 'M', 'A', 'R', 'Q',
              'L', 'V', 'Y', 'T', 'I', 'A', 'D', 'V', 'S', 'S', 'X', 'V', 'I',
              'Y', 'G', 'V', 'L', 'L', 'G', 'N', 'L', 'A', 'I', 'T'],
             ['G', 'R', 'P', 'E', 'W', 'I', 'W', 'L', 'A', 'L', 'G', 'T', 'A',
              'L', 'M', 'G', 'L', 'G', 'T', 'L', 'Y', 'F', 'L', 'V', 'K', 'G',
              'M', 'G', 'V', 'S', 'D', 'P', 'D', 'A', 'K', 'K', 'F', 'Y', 'A',
              'I', 'T', 'T', 'L', 'V', 'P', 'A', 'I', 'A', 'F', 'T', 'M', 'Y',
              'L', 'S', 'M', 'L', 'L', 'G', 'Y', 'G', 'L', 'T', 'M', 'V', 'P',
              'F', 'G', '-', '-', '-', '-', '-', '-', '-', 'G', 'E', 'Q', 'N',
              'P', 'I', '-', '-', 'Y', 'W', 'A', 'R', 'Y', 'A', 'D', 'W', 'L',
              'F', 'T', 'T', 'P', 'L', 'L', 'L', 'L', 'D', 'L', 'A', 'L', 'L',
              'V', 'D', 'A', 'D', 'Q', 'G', 'T', 'I', 'L', '-', '-', '-', '-',
              'A', 'L', 'V', 'G', 'A', 'D', 'G', 'I', 'M', 'I', 'G', 'T', 'G',
              'L', 'V', 'G', 'A', 'L', 'T', 'K', 'V', '-', '-', 'Y', 'S', 'Y',
              'R', 'F', 'V', 'W', 'W', 'A', 'I', 'S', 'T', 'A', 'A', 'M', 'L',
              'Y', 'I', 'L', 'Y', 'V', 'L', 'F', 'F', 'G', 'F', 'T', 'S', 'K',
              'A', 'E', 'S', 'M', 'R', 'P', 'E', 'V', 'A', 'S', 'T', 'F', 'K',
              'V', 'L', 'R', 'N', 'V', 'T', 'V', 'V', 'L', 'W', 'S', 'A', 'Y',
              'P', 'V', 'V', 'W', 'L', 'I', 'G', 'S', 'E', '-', '-', '-', '-',
              'G', 'A', 'G', 'I', 'V', '-', '-', 'P', 'L', 'N', 'I', 'E', 'T',
              'L', 'L', 'F', 'M', 'V', 'L', 'D', 'V', 'S', 'A', 'K', 'V', 'G',
              'F', 'G', 'L', 'I', 'L', 'L', 'R', 'S', 'R', 'A', 'I']],
            dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 100.00)
        self.assertAlmostEqual(alignment.annotations["E-value"], 1.3e-38)
        self.assertAlmostEqual(alignment.annotations["Score"], 275.60)
        self.assertAlmostEqual(alignment.annotations["Identities"], 16)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.231)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 6.600)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "CCcHHHHHHHHHHHHHHHHHHHHHHhcCCCChhhHHHHHHHHHHHHHHHHHHHHHHhCCCceeeccCCCCCcchHHHHHHHHHHHHHHHHHHHHHcCCCHHHHHHHHHHHHHHHHHHHHHHHcccchHHHHHHHHHHHHHHHHHHHHHHhHHHHHHhCCHHHHHHHHHHHHHHHHHHHHHHHHHHHccCCCCCCChhHHHHHHHHHHHHHHHHHHHHHHHhhhhhC",
        )
        self.assertEqual(alignment.query.id, "4Y9H:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[0:226],
            "GRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFG",
        )
        self.assertEqual(
            alignment[1],
            "GRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAIT---TLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQ---GTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESM-RPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFG",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "~~~~~~~~~~~~~~m~~~~~~f~~~~~~~~~~~~r~~~~l~~~i~~ia~~~Y~~ma~~~g~~~~~~~~~~~~v~~~RYv~W~ittPlll~~l~~l~~~~~~~~~~~v~~~~~mi~~g~~g~~~~~~~~~~~~~~~s~~~~~~i~~~l~~~~~~~~~~~~~~~~~~~~~l~~~~~i~W~~YPi~w~l~~~g~~~i~~~~~~~~~~ilDi~~K~~f~~~l~~~~~~~~",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                                                                                       ~~~~~~~lwv~f~~m~~~~l~f~~~~~r~~~~~r~~~~l~~~~~~I~~ia~~sY~~ma~~~~~~gr~v~~~RYvdWliTtPLLL~~L~lLag~~~~~~~~~~~lv~~d~lMI~~G~~gal~~~~~kw~~f~is~~~~l~l~~~l~~~~~~a~~~~~~~~~~~~~~~L~~l~~v~W~lYPIvw~L~~eG~~~is~~~e~i~y~ilDilsKv~fg~ll~~~~~~i~                                    ",
        )
        self.assertEqual(alignment.target.id, "5ZIH_B")
        self.assertEqual(
            alignment.target.seq[87:311],
            "KIGAQVCQWIAFSIAIALLTFYGFSAWKATCGWEEVYVCCVEVLFVTLEIFKEFSSPATVYLSTGNHAYCLRYFEWLLSCPVILIKLSNLSGLKNDYSKRTMGLIVSCVGMIVFGMAAGLATDWLKWLLYIVSCIYGGYMYFQAAKCYVEANHSVPKGHCRMVVKLMAYAYFASWGSYPILWAVGPEGLLKLSPYANSIGHSICDIIAXEFWTFLAHHLRIKIH",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "5ZIH_B")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Sensory opsin A,Chrimson; membrane protein, rhodopsin, ion channel; HET: OLC, LYR; 2.6A {Chlamydomonas reinhardtii}; Related PDB entries: 5ZIH_A",
        )
        self.assertEqual(
            alignment[0],
            "KIGAQVCQWIAFSIAIALLTFYGFSAWKAT-CGWEEVYVCCVEVLFVTLEIFKEFSSPATVYLS-------TGNHAYCLRYFEWLLSCPVILIKLSNLSGLKNDYSKRTMGLIVSCVGMIVFGMAAGLATD-WLKWLLYIVSCIYGGYMYFQAAKCYVEANHSVPKGHCRMVVKLMAYAYFASWGSYPILWAVGPEGLLKLSPYANSIGHSICDIIAXEFWTFLAHHLRIKIH",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "                                                                                       HHHHHHHHHHHHHHHHHHHHHHCCCCCTTHHHHHHHHHHHHHHHHHHHHTTSTTTSEEBTTSCEECHHHHHHHHHHHHHHHHHHHTCSSSSSCCHHHHHHHHHHHHHHHHHHHHHHHCCHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHTTSCTTHHHHHHHHHHHHHHHHHTHHHHHHHHSTTTTSCSCHHHHHHHHHHHHCCCCCCHHHHHHHHHHHHH                                       ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                                                                                       chhHHHHHHHHHHHHHHHHHHHHHHHccCCCCchhHHHHHHHHHHHHHHHHHHHHHHHHcCcCCCccccHHHHHHHHHHHHHHHHHHHHHhCCCCcchHHHHHHHHHHHHHHHHHHHHHHccchHHHHHHHHHHHHHHHHHHHHHHHHHHHhccCCchHHHHHHHHHHHHHHHHHHHHHHHHHHCcccCCCCCHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH                                    ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 87, 117, 117, 127, 130, 150, 150, 182, 185, 210, 210, 242, 243,
                              311],
                             [  0,  30,  31,  41,  41,  61,  68, 100, 100, 125, 126, 158, 158,
                              226]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['K', 'I', 'G', 'A', 'Q', 'V', 'C', 'Q', 'W', 'I', 'A', 'F', 'S',
              'I', 'A', 'I', 'A', 'L', 'L', 'T', 'F', 'Y', 'G', 'F', 'S', 'A',
              'W', 'K', 'A', 'T', '-', 'C', 'G', 'W', 'E', 'E', 'V', 'Y', 'V',
              'C', 'C', 'V', 'E', 'V', 'L', 'F', 'V', 'T', 'L', 'E', 'I', 'F',
              'K', 'E', 'F', 'S', 'S', 'P', 'A', 'T', 'V', 'Y', 'L', 'S', '-',
              '-', '-', '-', '-', '-', '-', 'T', 'G', 'N', 'H', 'A', 'Y', 'C',
              'L', 'R', 'Y', 'F', 'E', 'W', 'L', 'L', 'S', 'C', 'P', 'V', 'I',
              'L', 'I', 'K', 'L', 'S', 'N', 'L', 'S', 'G', 'L', 'K', 'N', 'D',
              'Y', 'S', 'K', 'R', 'T', 'M', 'G', 'L', 'I', 'V', 'S', 'C', 'V',
              'G', 'M', 'I', 'V', 'F', 'G', 'M', 'A', 'A', 'G', 'L', 'A', 'T',
              'D', '-', 'W', 'L', 'K', 'W', 'L', 'L', 'Y', 'I', 'V', 'S', 'C',
              'I', 'Y', 'G', 'G', 'Y', 'M', 'Y', 'F', 'Q', 'A', 'A', 'K', 'C',
              'Y', 'V', 'E', 'A', 'N', 'H', 'S', 'V', 'P', 'K', 'G', 'H', 'C',
              'R', 'M', 'V', 'V', 'K', 'L', 'M', 'A', 'Y', 'A', 'Y', 'F', 'A',
              'S', 'W', 'G', 'S', 'Y', 'P', 'I', 'L', 'W', 'A', 'V', 'G', 'P',
              'E', 'G', 'L', 'L', 'K', 'L', 'S', 'P', 'Y', 'A', 'N', 'S', 'I',
              'G', 'H', 'S', 'I', 'C', 'D', 'I', 'I', 'A', 'X', 'E', 'F', 'W',
              'T', 'F', 'L', 'A', 'H', 'H', 'L', 'R', 'I', 'K', 'I', 'H'],
             ['G', 'R', 'P', 'E', 'W', 'I', 'W', 'L', 'A', 'L', 'G', 'T', 'A',
              'L', 'M', 'G', 'L', 'G', 'T', 'L', 'Y', 'F', 'L', 'V', 'K', 'G',
              'M', 'G', 'V', 'S', 'D', 'P', 'D', 'A', 'K', 'K', 'F', 'Y', 'A',
              'I', 'T', '-', '-', '-', 'T', 'L', 'V', 'P', 'A', 'I', 'A', 'F',
              'T', 'M', 'Y', 'L', 'S', 'M', 'L', 'L', 'G', 'Y', 'G', 'L', 'T',
              'M', 'V', 'P', 'F', 'G', 'G', 'E', 'Q', 'N', 'P', 'I', 'Y', 'W',
              'A', 'R', 'Y', 'A', 'D', 'W', 'L', 'F', 'T', 'T', 'P', 'L', 'L',
              'L', 'L', 'D', 'L', 'A', 'L', 'L', 'V', 'D', 'A', 'D', 'Q', '-',
              '-', '-', 'G', 'T', 'I', 'L', 'A', 'L', 'V', 'G', 'A', 'D', 'G',
              'I', 'M', 'I', 'G', 'T', 'G', 'L', 'V', 'G', 'A', 'L', 'T', 'K',
              'V', 'Y', 'S', 'Y', 'R', 'F', 'V', 'W', 'W', 'A', 'I', 'S', 'T',
              'A', 'A', 'M', 'L', 'Y', 'I', 'L', 'Y', 'V', 'L', 'F', 'F', 'G',
              'F', 'T', 'S', 'K', 'A', 'E', 'S', 'M', '-', 'R', 'P', 'E', 'V',
              'A', 'S', 'T', 'F', 'K', 'V', 'L', 'R', 'N', 'V', 'T', 'V', 'V',
              'L', 'W', 'S', 'A', 'Y', 'P', 'V', 'V', 'W', 'L', 'I', 'G', 'S',
              'E', 'G', 'A', 'G', 'I', 'V', 'P', 'L', 'N', 'I', 'E', 'T', 'L',
              'L', 'F', 'M', 'V', 'L', 'D', 'V', 'S', 'A', 'K', 'V', 'G', 'F',
              'G', 'L', 'I', 'L', 'L', 'R', 'S', 'R', 'A', 'I', 'F', 'G']],
            dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 100.00)
        self.assertAlmostEqual(alignment.annotations["E-value"], 3.1e-38)
        self.assertAlmostEqual(alignment.annotations["Score"], 261.06)
        self.assertAlmostEqual(alignment.annotations["Identities"], 26)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.424)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 8.300)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "     HHHHHHHHHHHHHHHHHHHHhcCCCChhhHHHHHHHHHHHHHHHHHHHHHHhCCCceeeccCCCCCcchHHHHHHHHHHHHHHHHHHHHHcCCCHHHHHHHHHHHHHHHHHHHHHHHcccchHHHHHHHHHHHHHHHHHHHHHHhHHHHHHhCCHHHHHHHHHHHHHHHHHHHHHHHHHHHccCCCCCCChhHHHHHHHHHHHHHHHHHHHHHHHhhhh  ",
        )
        self.assertEqual(alignment.query.id, "4Y9H:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[5:224],
            "IWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAI",
        )
        self.assertEqual(
            alignment[1],
            "IWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDA----DQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAES--MRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAI",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "     ~~~~~~~~~m~~~~~~f~~~~~~~~~~~~r~~~~l~~~i~~ia~~~Y~~ma~~~g~~~~~~~~~~~~v~~~RYv~W~ittPlll~~l~~l~~~~~~~~~~~v~~~~~mi~~g~~g~~~~~~~~~~~~~~~s~~~~~~i~~~l~~~~~~~~~~~~~~~~~~~~~l~~~~~i~W~~YPi~w~l~~~g~~~i~~~~~~~~~~ilDi~~K~~f~~~l~~~~~~  ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "            ~~lw~~~~~~~~~~~~~~~~~~~~~~~~R~~~~l~~~i~~ia~~~Y~~m~~~~~~~g~~~~~~RYv~W~~TtPlll~~l~~l~~~~~~~~~~~~~~~~~~~~~mi~~G~~~~~~~~~~~~~~~~~~~~~~~v~~~l~~~~~~~~~~~~~~~~~~~~~~~l~~~~~v~W~~YPi~w~l~~~g~~i~~~~e~i~~~ilDv~~K~~f~~~ll~~~~~         ",
        )
        self.assertEqual(alignment.target.id, "4JQ6_C")
        self.assertEqual(
            alignment.target.seq[12:226],
            "ISFWLAAAIMLASTVFFFVERSDVPVKWKTSLTVAGLVTGVAFWHYLYMRGVWIYAGETPTVFRYIDWLITVPLQIIEFYLIIAAVTAISSAVFWKLLIASLVMLIGGFIGEAGLGDVVVWWIVGMIAWLYIIYEIFLGETAKANAGSGNAASQQAFNTIKWIVTVGWAIYPIGYAWGYFGDGLNEDALNIVYNLADLINKAAFGLAIWAAAMK",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "4JQ6_C")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Proteorhodopsin; retinylidene protein, ion transport, PROTON; HET: LI1, RET; 2.31A {uncultured bacterium} SCOP: f.13.1.0; Related PDB entries: 4JQ6_A 4JQ6_B",
        )
        self.assertEqual(
            alignment[0],
            "ISFWLAAAIMLASTVFFFVERSDVP-VKWKTSLTVAGLVTGVAFWHYLYMRGVWIY-------AGETPTVFRYIDWLITVPLQIIEFYLIIAAVTAISSAVFWKLLIASLVMLIGGFIGEAGLGD--VVVWWIVGMIAWLYIIYEIFLGETAKANAGSGNAASQQAFNTIKWIVTVGWAIYPIGYAW-GYFGDGLNEDALNIVYNLADLINKAAFGLAIWAAAMK",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "            HHHHHHHHHHHHHHHHHHHHGGGSCHHHHHHHHHHHHHHHHHHHHHHHHHHCCCCCCCCCHHHHHHHHHHHHHHHHHHHHHHHCCHHHHHHHHHHHHHHHHHHHHTSSCHHHHHHHHHHHHHHHHHTTCCTHHHHHHHHCCCCCHHHHHHHHHHHHSSSSCCHHHHHHHHHHHHCCCCCCHHHHHHHHHHH                                ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "            HHHHHHHHHHHHHHHHHHHHccCCChHHHHHHHHHHHHHHHHHHHHHHHHcchhhcCCCCchHHHHHHHhHHHHHHHHHHHHHHHHccCCHHHHHHHHHHHHHHHHHHHHHHhccccHHHHHHHHHHHHHHHHHHHHhhHHHHHhhCCCCHHHHHHHHHHHHHHHHHHHHHHHHHHHhhhccCCCHHHHHHHHHHHHHHHHHHHHHHHHHHHHH         ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 12,  37,  37,  67,  67,  97, 101, 129, 129, 158, 160, 189, 189,
                              226],
                             [  5,  30,  31,  61,  68,  98,  98, 126, 128, 157, 157, 186, 187,
                              224]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['I', 'S', 'F', 'W', 'L', 'A', 'A', 'A', 'I', 'M', 'L', 'A', 'S',
              'T', 'V', 'F', 'F', 'F', 'V', 'E', 'R', 'S', 'D', 'V', 'P', '-',
              'V', 'K', 'W', 'K', 'T', 'S', 'L', 'T', 'V', 'A', 'G', 'L', 'V',
              'T', 'G', 'V', 'A', 'F', 'W', 'H', 'Y', 'L', 'Y', 'M', 'R', 'G',
              'V', 'W', 'I', 'Y', '-', '-', '-', '-', '-', '-', '-', 'A', 'G',
              'E', 'T', 'P', 'T', 'V', 'F', 'R', 'Y', 'I', 'D', 'W', 'L', 'I',
              'T', 'V', 'P', 'L', 'Q', 'I', 'I', 'E', 'F', 'Y', 'L', 'I', 'I',
              'A', 'A', 'V', 'T', 'A', 'I', 'S', 'S', 'A', 'V', 'F', 'W', 'K',
              'L', 'L', 'I', 'A', 'S', 'L', 'V', 'M', 'L', 'I', 'G', 'G', 'F',
              'I', 'G', 'E', 'A', 'G', 'L', 'G', 'D', '-', '-', 'V', 'V', 'V',
              'W', 'W', 'I', 'V', 'G', 'M', 'I', 'A', 'W', 'L', 'Y', 'I', 'I',
              'Y', 'E', 'I', 'F', 'L', 'G', 'E', 'T', 'A', 'K', 'A', 'N', 'A',
              'G', 'S', 'G', 'N', 'A', 'A', 'S', 'Q', 'Q', 'A', 'F', 'N', 'T',
              'I', 'K', 'W', 'I', 'V', 'T', 'V', 'G', 'W', 'A', 'I', 'Y', 'P',
              'I', 'G', 'Y', 'A', 'W', '-', 'G', 'Y', 'F', 'G', 'D', 'G', 'L',
              'N', 'E', 'D', 'A', 'L', 'N', 'I', 'V', 'Y', 'N', 'L', 'A', 'D',
              'L', 'I', 'N', 'K', 'A', 'A', 'F', 'G', 'L', 'A', 'I', 'W', 'A',
              'A', 'A', 'M', 'K'],
             ['I', 'W', 'L', 'A', 'L', 'G', 'T', 'A', 'L', 'M', 'G', 'L', 'G',
              'T', 'L', 'Y', 'F', 'L', 'V', 'K', 'G', 'M', 'G', 'V', 'S', 'D',
              'P', 'D', 'A', 'K', 'K', 'F', 'Y', 'A', 'I', 'T', 'T', 'L', 'V',
              'P', 'A', 'I', 'A', 'F', 'T', 'M', 'Y', 'L', 'S', 'M', 'L', 'L',
              'G', 'Y', 'G', 'L', 'T', 'M', 'V', 'P', 'F', 'G', 'G', 'E', 'Q',
              'N', 'P', 'I', 'Y', 'W', 'A', 'R', 'Y', 'A', 'D', 'W', 'L', 'F',
              'T', 'T', 'P', 'L', 'L', 'L', 'L', 'D', 'L', 'A', 'L', 'L', 'V',
              'D', 'A', '-', '-', '-', '-', 'D', 'Q', 'G', 'T', 'I', 'L', 'A',
              'L', 'V', 'G', 'A', 'D', 'G', 'I', 'M', 'I', 'G', 'T', 'G', 'L',
              'V', 'G', 'A', 'L', 'T', 'K', 'V', 'Y', 'S', 'Y', 'R', 'F', 'V',
              'W', 'W', 'A', 'I', 'S', 'T', 'A', 'A', 'M', 'L', 'Y', 'I', 'L',
              'Y', 'V', 'L', 'F', 'F', 'G', 'F', 'T', 'S', 'K', 'A', 'E', 'S',
              '-', '-', 'M', 'R', 'P', 'E', 'V', 'A', 'S', 'T', 'F', 'K', 'V',
              'L', 'R', 'N', 'V', 'T', 'V', 'V', 'L', 'W', 'S', 'A', 'Y', 'P',
              'V', 'V', 'W', 'L', 'I', 'G', 'S', 'E', 'G', 'A', 'G', 'I', 'V',
              'P', 'L', 'N', 'I', 'E', 'T', 'L', 'L', 'F', 'M', 'V', 'L', 'D',
              'V', 'S', 'A', 'K', 'V', 'G', 'F', 'G', 'L', 'I', 'L', 'L', 'R',
              'S', 'R', 'A', 'I']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 100.00)
        self.assertAlmostEqual(alignment.annotations["E-value"], 3.9e-38)
        self.assertAlmostEqual(alignment.annotations["Score"], 265.93)
        self.assertAlmostEqual(alignment.annotations["Identities"], 24)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.298)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 7.600)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "      HHHHHHHHHHHHHHHHHHHhcCCCChhhHHHHHHHHHHHHHHHHHHHHHHhCCCceeeccCCCCCcchHHHHHHHHHHHHHHHHHHHHHcCCCHHHHHHHHHHHHHHHHHHHHHHHcccchHHHHHHHHHHHHHHHHHHHHHHhHHHHHHhCCHHHHHHHHHHHHHHHHHHHHHHHHHHHccCCCCCCChhHHHHHHHHHHHHHHHHHHHHHHHhhhh  ",
        )
        self.assertEqual(alignment.query.id, "4Y9H:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[6:224],
            "WLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAI",
        )
        self.assertEqual(
            alignment[1],
            "WLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPF------GGEQNPIYWARYADWLFTTPLLLLDLALLVDADQG----TILALVGADGIMIGTGLVGALTKVY--SYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAI",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "      ~~~~~~~~m~~~~~~f~~~~~~~~~~~~r~~~~l~~~i~~ia~~~Y~~ma~~~g~~~~~~~~~~~~v~~~RYv~W~ittPlll~~l~~l~~~~~~~~~~~v~~~~~mi~~g~~g~~~~~~~~~~~~~~~s~~~~~~i~~~l~~~~~~~~~~~~~~~~~~~~~l~~~~~i~W~~YPi~w~l~~~g~~~i~~~~~~~~~~ilDi~~K~~f~~~l~~~~~~  ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                       ~l~~~~~~~~~~~i~f~~~~~~~~~~~R~~~~l~~~v~~vaal~Y~~ma~~~~~~~~~~~~~~~~~~g~~~~~~RYvdW~lTtPlll~~L~llag~~~~~~~~~~~~~i~~~~~mi~~G~~g~~~~~~~~~~~~~w~~~s~~~f~~v~~~l~~~~~~~~~~~~~~~~~~~~~l~~~~~v~W~~YPi~w~l~~~g~~~is~~~~~i~~~ilDvl~K~~f~~~ll~~~~~                        ",
        )
        self.assertEqual(alignment.target.id, "5G28_A")
        self.assertEqual(
            alignment.target.seq[23:251],
            "LLTMGVGVHFAALIFFLVVSQFVAPKYRIATALSCIVMVSAGLILNSQAVMWTDAYAYVDGSYQLQDLTFSNGYRYVNWMATIPCLLLQLLIVLNLKGKELFSTATWLILAAWGMIITGYVGQLYEVDDIAQLMIWGAVSTAFFVVMNWIVGTKIFKNRATMLGGTDSTITKVFWLMMFAWTLYPIAYLVPAFMNNADGVVLRQLLFTIADISSKVIYGLMITYIAIQ",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "5G28_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "CHLORIDE PUMPING RHODOPSIN; SIGNALING PROTEIN; HET: OLA, RET; 1.57A {NONLABENS MARINUS S1-08}; Related PDB entries: 5B2N_A 5G2D_A 5G54_A 5G2C_A 5G2A_A 5FJG_A",
        )
        self.assertEqual(
            alignment[0],
            "LLTMGVGVHFAALIFFLVVSQFVA-PKYRIATALSCIVMVSAGLILNSQAVMWTDAYAYVDGSYQLQDLTFSN-GYRYVNWMATIPCLLLQLLIVLNLKGKELFSTATWLILAAWGMIITGYVGQLYEVDDIAQLMIWGAVSTAFFVVMNWIVGTKIFKNRATMLGGTDSTITKVFWLMMFAWTLYPIAYLVPAFMNNADGVVLRQLLFTIADISSKVIYGLMITYIAIQ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "                       HHHHHHHHHHHHHHHHHHHGGGSCTTTHHHHHHHHHHHHHHHHHHHHHHHHHHHHEEEETTEEEECSSCCCSHHHHHHHHHHHHHHHHHHHHHTTCCHHHHHHHHHHHHHHHHHHHHHHHHHHTTTTTCHHHHHHHHHHHHHHHHHHHHHHHHHHHHHGGGCCTTHHHHHHHHHHHHHHHHHHHHHHHHHHHHCCSHHHHHHHHHHHHHHHCCCCCCHHHHHHHHHHH                        ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                       HHHHHHHHHHHHHHHHHHHhccCCchhHHHHHHHHHHHHHHHHHHHHHHcccccceeeeCCccccCCcccccHHHHHHHHHHHHHHHHHHHHHHccCCcchHHHHHHHHHHHHHHHHHHHHHHhcCCChhHHHHHHHHHHHHHHHHHHHHHhhHHHHHHHhccCCchHHHHHHHHHHHHHHHHHHHHHHhhhhcCCCCcHHHHHHHHHHHHHHHHHHHHHHHHHHHHH                        ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 23,  47,  47,  82,  88,  95,  95, 122, 126, 151, 153, 251],
                             [  6,  30,  31,  66,  66,  73,  74, 101, 101, 126, 126, 224]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['L', 'L', 'T', 'M', 'G', 'V', 'G', 'V', 'H', 'F', 'A', 'A', 'L',
              'I', 'F', 'F', 'L', 'V', 'V', 'S', 'Q', 'F', 'V', 'A', '-', 'P',
              'K', 'Y', 'R', 'I', 'A', 'T', 'A', 'L', 'S', 'C', 'I', 'V', 'M',
              'V', 'S', 'A', 'G', 'L', 'I', 'L', 'N', 'S', 'Q', 'A', 'V', 'M',
              'W', 'T', 'D', 'A', 'Y', 'A', 'Y', 'V', 'D', 'G', 'S', 'Y', 'Q',
              'L', 'Q', 'D', 'L', 'T', 'F', 'S', 'N', '-', 'G', 'Y', 'R', 'Y',
              'V', 'N', 'W', 'M', 'A', 'T', 'I', 'P', 'C', 'L', 'L', 'L', 'Q',
              'L', 'L', 'I', 'V', 'L', 'N', 'L', 'K', 'G', 'K', 'E', 'L', 'F',
              'S', 'T', 'A', 'T', 'W', 'L', 'I', 'L', 'A', 'A', 'W', 'G', 'M',
              'I', 'I', 'T', 'G', 'Y', 'V', 'G', 'Q', 'L', 'Y', 'E', 'V', 'D',
              'D', 'I', 'A', 'Q', 'L', 'M', 'I', 'W', 'G', 'A', 'V', 'S', 'T',
              'A', 'F', 'F', 'V', 'V', 'M', 'N', 'W', 'I', 'V', 'G', 'T', 'K',
              'I', 'F', 'K', 'N', 'R', 'A', 'T', 'M', 'L', 'G', 'G', 'T', 'D',
              'S', 'T', 'I', 'T', 'K', 'V', 'F', 'W', 'L', 'M', 'M', 'F', 'A',
              'W', 'T', 'L', 'Y', 'P', 'I', 'A', 'Y', 'L', 'V', 'P', 'A', 'F',
              'M', 'N', 'N', 'A', 'D', 'G', 'V', 'V', 'L', 'R', 'Q', 'L', 'L',
              'F', 'T', 'I', 'A', 'D', 'I', 'S', 'S', 'K', 'V', 'I', 'Y', 'G',
              'L', 'M', 'I', 'T', 'Y', 'I', 'A', 'I', 'Q'],
             ['W', 'L', 'A', 'L', 'G', 'T', 'A', 'L', 'M', 'G', 'L', 'G', 'T',
              'L', 'Y', 'F', 'L', 'V', 'K', 'G', 'M', 'G', 'V', 'S', 'D', 'P',
              'D', 'A', 'K', 'K', 'F', 'Y', 'A', 'I', 'T', 'T', 'L', 'V', 'P',
              'A', 'I', 'A', 'F', 'T', 'M', 'Y', 'L', 'S', 'M', 'L', 'L', 'G',
              'Y', 'G', 'L', 'T', 'M', 'V', 'P', 'F', '-', '-', '-', '-', '-',
              '-', 'G', 'G', 'E', 'Q', 'N', 'P', 'I', 'Y', 'W', 'A', 'R', 'Y',
              'A', 'D', 'W', 'L', 'F', 'T', 'T', 'P', 'L', 'L', 'L', 'L', 'D',
              'L', 'A', 'L', 'L', 'V', 'D', 'A', 'D', 'Q', 'G', '-', '-', '-',
              '-', 'T', 'I', 'L', 'A', 'L', 'V', 'G', 'A', 'D', 'G', 'I', 'M',
              'I', 'G', 'T', 'G', 'L', 'V', 'G', 'A', 'L', 'T', 'K', 'V', 'Y',
              '-', '-', 'S', 'Y', 'R', 'F', 'V', 'W', 'W', 'A', 'I', 'S', 'T',
              'A', 'A', 'M', 'L', 'Y', 'I', 'L', 'Y', 'V', 'L', 'F', 'F', 'G',
              'F', 'T', 'S', 'K', 'A', 'E', 'S', 'M', 'R', 'P', 'E', 'V', 'A',
              'S', 'T', 'F', 'K', 'V', 'L', 'R', 'N', 'V', 'T', 'V', 'V', 'L',
              'W', 'S', 'A', 'Y', 'P', 'V', 'V', 'W', 'L', 'I', 'G', 'S', 'E',
              'G', 'A', 'G', 'I', 'V', 'P', 'L', 'N', 'I', 'E', 'T', 'L', 'L',
              'F', 'M', 'V', 'L', 'D', 'V', 'S', 'A', 'K', 'V', 'G', 'F', 'G',
              'L', 'I', 'L', 'L', 'R', 'S', 'R', 'A', 'I']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 100.00)
        self.assertAlmostEqual(alignment.annotations["E-value"], 4e-38)
        self.assertAlmostEqual(alignment.annotations["Score"], 263.73)
        self.assertAlmostEqual(alignment.annotations["Identities"], 28)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.412)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 7.900)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "CCcHHHHHHHHHHHHHHHHHHHHHHhcCCCChhhHHHHHHHHHHHHHHHHHHHHHHhCCCceeeccCCCCCcchHHHHHHHHHHHHHHHHHHHHHcCCCHHHHHHHHHHHHHHHHHHHHHHHcccchHHHHHHHHHHHHHHHHHHHHHHhHHHHHHhCCHHHHHHHHHHHHHHHHHHHHHHHHHHHccCCCCCCChhHHHHHHHHHHHHHHHHHHHHHHHhhhh  ",
        )
        self.assertEqual(alignment.query.id, "4Y9H:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[0:224],
            "GRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAI",
        )
        self.assertEqual(
            alignment[1],
            "GRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDA----DQGTILALVGADGIMIGTGLVGALTKV-----YSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAI",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "~~~~~~~~~~~~~~m~~~~~~f~~~~~~~~~~~~r~~~~l~~~i~~ia~~~Y~~ma~~~g~~~~~~~~~~~~v~~~RYv~W~ittPlll~~l~~l~~~~~~~~~~~v~~~~~mi~~g~~g~~~~~~~~~~~~~~~s~~~~~~i~~~l~~~~~~~~~~~~~~~~~~~~~l~~~~~i~W~~YPi~w~l~~~g~~~i~~~~~~~~~~ilDi~~K~~f~~~l~~~~~~  ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "       ~~~~~~~lw~~~~~~~~~~l~f~~~~~~~~~~~R~~~~l~~~i~~vaa~~Y~~ma~~~g~~~~~~~~~~~~~~~RYvdW~~TtPlll~~l~llag~~~~~~~~~~~~lv~~~~~mi~~G~~g~~~~~~~~~~~~~~~~~f~~s~~~~~~vl~~l~~~~~~~a~~~~~~~~~~~~~l~~~~~v~W~~YPivw~l~~~~~g~~i~~~~~i~~~vlDv~aK~~fg~~ll~~~~~                    ",
        )
        self.assertEqual(alignment.target.id, "4HYJ_A")
        self.assertEqual(
            alignment.target.seq[7:238],
            "VLATQYMFWVGFVGMAAGTLYFLVERNSLAPEYRSTATVAALVTFVAAIHYYFMKDAVGTSGLLSEIDGFPTEIRYIDWLVTTPLLLVKFPLLLGLKGRLGRPLLTKLVIADVIMIVGGYIGESSINIAGGFTQLGLWSYLIGCFAWIYIIYLLFTNVTKAAENKPAPIRDALLKMRLFILIGWAIYPIGYAVTLFAPGVEIQLVRELIYNFADLTNKVGFGLIAFFAVKT",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "4HYJ_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Rhodopsin; Seven-helical transmembrane protein, Proton pump; HET: RET, LFA; 2.3A {Exiguobacterium sibiricum} SCOP: f.13.1.0; Related PDB entries: 4HYJ_B",
        )
        self.assertEqual(
            alignment[0],
            "VLATQYMFWVGFVGMAAGTLYFLVERNSLA-PEYRSTATVAALVTFVAAIHYYFMKDAVGTSGLLSEIDGFPT-EIRYIDWLVTTPLLLVKFPLLLGLKGRLGRPLLTKLVIADVIMIVGGYIGESSINIAGGFTQLGLWSYLIGCFAWIYIIYLLFTNVTKAAENKPAPIRDALLKMRLFILIGWAIYPIGYAVTLFAPGVEIQLVRELIYNFADLTNKVGFGLIAFFAVKT",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "       HHHHHHHHHHHHHHHHHHHHHHHHSGGGSCTTTHHHHHHHHHHHHHHHHHHHHHHHHHCSSCCGGGCCSCHHHHHHHHHHHHHHHHHTHHHHHHCCCSHHHHHHHHHHHHHHHHHHHHHHHHHHHHTSCCHHHHHHHHHHHHHHHHHHHHCCCCCHHHHTTSCHHHHHHHHHHHHHHHTGGGHHHHHHHHHHHCCCHHHHHHHHHHHHHHHHCCCCCHHHHHHHHHHH                       ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "       hHHHHHHHHHHHHHHHHHHHHHHHHHcCCChhhHHHHHHHHHHHHHHHHHHHHHHhcccccccccccCCCCcHHHHhHHHhHHHHHHHHHHHHHccccccCHHHHHHHHHHHHHHHHHHHHHHHhhhccCCCcchHHHHHHHHHHHHHHHHHHHHHHHHHHHhCCCHHHHHHHHHHHHHHHHHHHHHHHHHHHhhcCCCcchhHHHHHHHHHHHHHHHHHHHHHHHHHHHh                    ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  7,  37,  37,  79,  79, 103, 107, 134, 139, 238],
                             [  0,  30,  31,  73,  74,  98,  98, 125, 125, 224]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['V', 'L', 'A', 'T', 'Q', 'Y', 'M', 'F', 'W', 'V', 'G', 'F', 'V',
              'G', 'M', 'A', 'A', 'G', 'T', 'L', 'Y', 'F', 'L', 'V', 'E', 'R',
              'N', 'S', 'L', 'A', '-', 'P', 'E', 'Y', 'R', 'S', 'T', 'A', 'T',
              'V', 'A', 'A', 'L', 'V', 'T', 'F', 'V', 'A', 'A', 'I', 'H', 'Y',
              'Y', 'F', 'M', 'K', 'D', 'A', 'V', 'G', 'T', 'S', 'G', 'L', 'L',
              'S', 'E', 'I', 'D', 'G', 'F', 'P', 'T', '-', 'E', 'I', 'R', 'Y',
              'I', 'D', 'W', 'L', 'V', 'T', 'T', 'P', 'L', 'L', 'L', 'V', 'K',
              'F', 'P', 'L', 'L', 'L', 'G', 'L', 'K', 'G', 'R', 'L', 'G', 'R',
              'P', 'L', 'L', 'T', 'K', 'L', 'V', 'I', 'A', 'D', 'V', 'I', 'M',
              'I', 'V', 'G', 'G', 'Y', 'I', 'G', 'E', 'S', 'S', 'I', 'N', 'I',
              'A', 'G', 'G', 'F', 'T', 'Q', 'L', 'G', 'L', 'W', 'S', 'Y', 'L',
              'I', 'G', 'C', 'F', 'A', 'W', 'I', 'Y', 'I', 'I', 'Y', 'L', 'L',
              'F', 'T', 'N', 'V', 'T', 'K', 'A', 'A', 'E', 'N', 'K', 'P', 'A',
              'P', 'I', 'R', 'D', 'A', 'L', 'L', 'K', 'M', 'R', 'L', 'F', 'I',
              'L', 'I', 'G', 'W', 'A', 'I', 'Y', 'P', 'I', 'G', 'Y', 'A', 'V',
              'T', 'L', 'F', 'A', 'P', 'G', 'V', 'E', 'I', 'Q', 'L', 'V', 'R',
              'E', 'L', 'I', 'Y', 'N', 'F', 'A', 'D', 'L', 'T', 'N', 'K', 'V',
              'G', 'F', 'G', 'L', 'I', 'A', 'F', 'F', 'A', 'V', 'K', 'T'],
             ['G', 'R', 'P', 'E', 'W', 'I', 'W', 'L', 'A', 'L', 'G', 'T', 'A',
              'L', 'M', 'G', 'L', 'G', 'T', 'L', 'Y', 'F', 'L', 'V', 'K', 'G',
              'M', 'G', 'V', 'S', 'D', 'P', 'D', 'A', 'K', 'K', 'F', 'Y', 'A',
              'I', 'T', 'T', 'L', 'V', 'P', 'A', 'I', 'A', 'F', 'T', 'M', 'Y',
              'L', 'S', 'M', 'L', 'L', 'G', 'Y', 'G', 'L', 'T', 'M', 'V', 'P',
              'F', 'G', 'G', 'E', 'Q', 'N', 'P', 'I', 'Y', 'W', 'A', 'R', 'Y',
              'A', 'D', 'W', 'L', 'F', 'T', 'T', 'P', 'L', 'L', 'L', 'L', 'D',
              'L', 'A', 'L', 'L', 'V', 'D', 'A', '-', '-', '-', '-', 'D', 'Q',
              'G', 'T', 'I', 'L', 'A', 'L', 'V', 'G', 'A', 'D', 'G', 'I', 'M',
              'I', 'G', 'T', 'G', 'L', 'V', 'G', 'A', 'L', 'T', 'K', 'V', '-',
              '-', '-', '-', '-', 'Y', 'S', 'Y', 'R', 'F', 'V', 'W', 'W', 'A',
              'I', 'S', 'T', 'A', 'A', 'M', 'L', 'Y', 'I', 'L', 'Y', 'V', 'L',
              'F', 'F', 'G', 'F', 'T', 'S', 'K', 'A', 'E', 'S', 'M', 'R', 'P',
              'E', 'V', 'A', 'S', 'T', 'F', 'K', 'V', 'L', 'R', 'N', 'V', 'T',
              'V', 'V', 'L', 'W', 'S', 'A', 'Y', 'P', 'V', 'V', 'W', 'L', 'I',
              'G', 'S', 'E', 'G', 'A', 'G', 'I', 'V', 'P', 'L', 'N', 'I', 'E',
              'T', 'L', 'L', 'F', 'M', 'V', 'L', 'D', 'V', 'S', 'A', 'K', 'V',
              'G', 'F', 'G', 'L', 'I', 'L', 'L', 'R', 'S', 'R', 'A', 'I']],
            dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 100.00)
        self.assertAlmostEqual(alignment.annotations["E-value"], 4.4e-38)
        self.assertAlmostEqual(alignment.annotations["Score"], 263.83)
        self.assertAlmostEqual(alignment.annotations["Identities"], 27)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.383)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 7.700)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "CCcHHHHHHHHHHHHHHHHHHHHHHhcCCCChhhHHHHHHHHHHHHHHHHHHHHHHhCCCceeeccCCCCCcchHHHHHHHHHHHHHHHHHHHHHcCCCHHHHHHHHHHHHHHHHHHHHHHHcccchHHHHHHHHHHHHHHHHHHHHHHhHHHHHHhCCHHHHHHHHHHHHHHHHHHHHHHHHHHHccCCCCCCChhHHHHHHHHHHHHHHHHHHHHHHHhhhh  ",
        )
        self.assertEqual(alignment.query.id, "4Y9H:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[0:224],
            "GRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAI",
        )
        self.assertEqual(
            alignment[1],
            "GRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDAD----QGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFG-FTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIG-SEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAI",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "~~~~~~~~~~~~~~m~~~~~~f~~~~~~~~~~~~r~~~~l~~~i~~ia~~~Y~~ma~~~g~~~~~~~~~~~~v~~~RYv~W~ittPlll~~l~~l~~~~~~~~~~~v~~~~~mi~~g~~g~~~~~~~~~~~~~~~s~~~~~~i~~~l~~~~~~~~~~~~~~~~~~~~~l~~~~~i~W~~YPi~w~l~~~g~~~i~~~~~~~~~~ilDi~~K~~f~~~l~~~~~~  ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                          ~~~~~~~lw~~~~~~~~~~~~f~~~~~~~~~~~R~~~~l~~~i~~ia~~aY~~ma~~~~~~g~~~~~~RYvdW~lTtPlll~~l~llag~~~~~~~~~~~~lv~~~~~Mi~~G~~g~~~~~~~w~~f~~g~~~f~~il~~l~~~~~~~~~~~~~~~~~~~~~~l~~~~~v~W~~YPi~w~l~~~eg~~~i~~~~e~i~y~ilDv~~K~~fg~~l~~~~~~               ",
        )
        self.assertEqual(alignment.target.id, "4KNF_D")
        self.assertEqual(
            alignment.target.seq[26:246],
            "SDTVGVSFWLVTAGMLAATVFFFVERDQVSAKWKTSLTVSGLITGIAFWHYLYMRGVWIDTGDTPTVFRYINWLLTVPLLVVEFYLILAACTSVAASLFKKLLAGSLVMLGAGFAGEAGLAPVLPAFIIGMAGWLYMIYELYMGEGKAAVSTASPAVNSAYNAMMMIIVVGWAIYPAGYAAGYLMGGEGVYASNLNLIYNLADFVNKILFGLIIWNVAVK",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "4KNF_D")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Blue-light absorbing proteorhodopsin; Membrane, Proton transport, ion transport; HET: RET; 2.6A {gamma proteobacterium 'Hot 75m4'} SCOP: l.1.1.1, f.13.1.0; Related PDB entries: 2L6X_A 4KNF_A 4KNF_B 4KNF_C 4KNF_E 4KLY_A 4KLY_B 4KLY_E 4KLY_C 4KLY_D",
        )
        self.assertEqual(
            alignment[0],
            "SDTVGVSFWLVTAGMLAATVFFFVERDQVS-AKWKTSLTVSGLITGIAFWHYLYMRGVWID-------TGDTPTVFRYINWLLTVPLLVVEFYLILAACTSVAASLFKKLLAGSLVMLGAGFAGEAGLA--PVLPAFIIGMAGWLYMIYELYMGEGKAAVSTASPAVNSAYNAMMMIIVVGWAIYPAGYAAGYLMGGEGVYASNLNLIYNLADFVNKILFGLIIWNVAVK",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "                          SCHHHHHHHHHHHHHHHHHHHHHHHGGGSCHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHCCCCHHHHHHHHHHHHHHHHHHHHHHHHHHCSCHHHHHHHHHHHHHHHHHHHHHHHTTSSCHHHHHHHHHHHHHHHHHHHHSTTSCCCCCSSHHHHHHHHHHCCCCCTTHHHHHHHHHHCCCCCCHHHHHHHHHHHHCCCCCCHHHHHHHHHHH                    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                          chHHHHHHHHHHHHHHHHHHHHHHHhcCCCcchHHHHHHHHHHHHHHHHHHHHHHhccccCCCCCcHHHHHHHHhHHHHHHHHHHHHHhhcccchHHHHHHHHHHHHHHHHHHHHHHhcccchHHHHHHHHHHHHHHHHHHHhhhHHHHHhcCCHHHHHHHHHHHHHHHHHHHHHHHHHHHHHhcCCCCCCHHHHHHHHHHHHHHHHHHHHHHHHHHHHH               ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 26,  56,  56,  86,  86, 117, 121, 147, 147, 170, 171, 208, 209,
                              246],
                             [  0,  30,  31,  61,  68,  99,  99, 125, 127, 150, 150, 187, 187,
                              224]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['S', 'D', 'T', 'V', 'G', 'V', 'S', 'F', 'W', 'L', 'V', 'T', 'A',
              'G', 'M', 'L', 'A', 'A', 'T', 'V', 'F', 'F', 'F', 'V', 'E', 'R',
              'D', 'Q', 'V', 'S', '-', 'A', 'K', 'W', 'K', 'T', 'S', 'L', 'T',
              'V', 'S', 'G', 'L', 'I', 'T', 'G', 'I', 'A', 'F', 'W', 'H', 'Y',
              'L', 'Y', 'M', 'R', 'G', 'V', 'W', 'I', 'D', '-', '-', '-', '-',
              '-', '-', '-', 'T', 'G', 'D', 'T', 'P', 'T', 'V', 'F', 'R', 'Y',
              'I', 'N', 'W', 'L', 'L', 'T', 'V', 'P', 'L', 'L', 'V', 'V', 'E',
              'F', 'Y', 'L', 'I', 'L', 'A', 'A', 'C', 'T', 'S', 'V', 'A', 'A',
              'S', 'L', 'F', 'K', 'K', 'L', 'L', 'A', 'G', 'S', 'L', 'V', 'M',
              'L', 'G', 'A', 'G', 'F', 'A', 'G', 'E', 'A', 'G', 'L', 'A', '-',
              '-', 'P', 'V', 'L', 'P', 'A', 'F', 'I', 'I', 'G', 'M', 'A', 'G',
              'W', 'L', 'Y', 'M', 'I', 'Y', 'E', 'L', 'Y', 'M', 'G', 'E', 'G',
              'K', 'A', 'A', 'V', 'S', 'T', 'A', 'S', 'P', 'A', 'V', 'N', 'S',
              'A', 'Y', 'N', 'A', 'M', 'M', 'M', 'I', 'I', 'V', 'V', 'G', 'W',
              'A', 'I', 'Y', 'P', 'A', 'G', 'Y', 'A', 'A', 'G', 'Y', 'L', 'M',
              'G', 'G', 'E', 'G', 'V', 'Y', 'A', 'S', 'N', 'L', 'N', 'L', 'I',
              'Y', 'N', 'L', 'A', 'D', 'F', 'V', 'N', 'K', 'I', 'L', 'F', 'G',
              'L', 'I', 'I', 'W', 'N', 'V', 'A', 'V', 'K'],
             ['G', 'R', 'P', 'E', 'W', 'I', 'W', 'L', 'A', 'L', 'G', 'T', 'A',
              'L', 'M', 'G', 'L', 'G', 'T', 'L', 'Y', 'F', 'L', 'V', 'K', 'G',
              'M', 'G', 'V', 'S', 'D', 'P', 'D', 'A', 'K', 'K', 'F', 'Y', 'A',
              'I', 'T', 'T', 'L', 'V', 'P', 'A', 'I', 'A', 'F', 'T', 'M', 'Y',
              'L', 'S', 'M', 'L', 'L', 'G', 'Y', 'G', 'L', 'T', 'M', 'V', 'P',
              'F', 'G', 'G', 'E', 'Q', 'N', 'P', 'I', 'Y', 'W', 'A', 'R', 'Y',
              'A', 'D', 'W', 'L', 'F', 'T', 'T', 'P', 'L', 'L', 'L', 'L', 'D',
              'L', 'A', 'L', 'L', 'V', 'D', 'A', 'D', '-', '-', '-', '-', 'Q',
              'G', 'T', 'I', 'L', 'A', 'L', 'V', 'G', 'A', 'D', 'G', 'I', 'M',
              'I', 'G', 'T', 'G', 'L', 'V', 'G', 'A', 'L', 'T', 'K', 'V', 'Y',
              'S', 'Y', 'R', 'F', 'V', 'W', 'W', 'A', 'I', 'S', 'T', 'A', 'A',
              'M', 'L', 'Y', 'I', 'L', 'Y', 'V', 'L', 'F', 'F', 'G', '-', 'F',
              'T', 'S', 'K', 'A', 'E', 'S', 'M', 'R', 'P', 'E', 'V', 'A', 'S',
              'T', 'F', 'K', 'V', 'L', 'R', 'N', 'V', 'T', 'V', 'V', 'L', 'W',
              'S', 'A', 'Y', 'P', 'V', 'V', 'W', 'L', 'I', 'G', '-', 'S', 'E',
              'G', 'A', 'G', 'I', 'V', 'P', 'L', 'N', 'I', 'E', 'T', 'L', 'L',
              'F', 'M', 'V', 'L', 'D', 'V', 'S', 'A', 'K', 'V', 'G', 'F', 'G',
              'L', 'I', 'L', 'L', 'R', 'S', 'R', 'A', 'I']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 100.00)
        self.assertAlmostEqual(alignment.annotations["E-value"], 6.6e-38)
        self.assertAlmostEqual(alignment.annotations["Score"], 263.77)
        self.assertAlmostEqual(alignment.annotations["Identities"], 27)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.402)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 7.900)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "       HHHHHHHHHHHHHHHHHHhcCCCChhhHHHHHHHHHHHHHHHHHHHHHHhCCCceeeccCCCCCcchHHHHHHHHHHHHHHHHHHHHHcCCCHHHHHHHHHHHHHHHHHHHHHHHcccchHHHHHHHHHHHHHHHHHHHHHHhHHHHHHhCCHHHHHHHHHHHHHHHHHHHHHHHHHHHccCCCCCCChhHHHHHHHHHHHHHHHHHHHHHHHhhhh  ",
        )
        self.assertEqual(alignment.query.id, "4Y9H:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[7:224],
            "LALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAI",
        )
        self.assertEqual(
            alignment[1],
            "LALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSML-LGYGLTMVPFG----GEQNPIYWARYADWLFTTPLLLLDLALLVDAD----QGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAES--MRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEG-AGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAI",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "       ~~~~~~~m~~~~~~f~~~~~~~~~~~~r~~~~l~~~i~~ia~~~Y~~ma~~~g~~~~~~~~~~~~v~~~RYv~W~ittPlll~~l~~l~~~~~~~~~~~v~~~~~mi~~g~~g~~~~~~~~~~~~~~~s~~~~~~i~~~l~~~~~~~~~~~~~~~~~~~~~l~~~~~i~W~~YPi~w~l~~~g~~~i~~~~~~~~~~ilDi~~K~~f~~~l~~~~~~  ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                  ~~~~~~~~~~~~~~f~~~~~~~~~~~r~~~~l~~~i~~va~~~Y~~ma~~~~g~~~~~~g~~~~~g~~~~~~~RYv~W~~TtPlll~~l~ll~g~~~~~~~~~~~~~i~~~~~mi~~G~~g~~~~~~~~~~~w~~is~~~f~~vl~~l~~~~~~~~~~~~~~~~~~~~~~~l~~~~~v~W~~YPi~w~l~~~~~~g~i~~~~~~i~y~vlDv~~K~~f~~~l~~~~~~                      ",
        )
        self.assertEqual(alignment.target.id, "5AZD_C")
        self.assertEqual(
            alignment.target.seq[18:246],
            "LSLTIAGMLAAFVFFLLARSYVAPRYHIALYLSALIVFIAGYHYLRIFESWVGAYQLQDGVYVPTGKPFNDFYRYADWLLTVPLLLLELILVLGLTAARTWNLSIKLVVASVLMLALGYVGEVNTEPGPRTLWGALSSIPFFYILYVLWVELGQAIREAKFGPRVLELLGATRLVLLMSWGFYPIAYALGTWLPGGAAQEVAIQIGYSLADLIAXPIYGLLVFAIARA",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "5AZD_C")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Thermophilic rhodopsin; Membrane protein, Retinal, Ion pump; HET: LYR; 2.8A {Thermus thermophilus JL-18}; Related PDB entries: 5AZD_D 5AZD_B 5AZD_A",
        )
        self.assertEqual(
            alignment[0],
            "LSLTIAGMLAAFVFFLLARSYVA-PRYHIALYLSALIVFIAGYHYLRIFESWVGAYQLQDGVYVPTGKPFNDFYRYADWLLTVPLLLLELILVLGLTAARTWNLSIKLVVASVLMLALGYVGEVNTEPGPRTLWGALSSIPFFYILYVLWVELGQAIREAKFGPRVLELLGATRLVLLMSWGFYPIAYALGTWLPGGAAQEVAIQIGYSLADLIAXPIYGLLVFAIARA",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "                  HHHHHHHHHHHHHHHHHGGGGSCGGGHHHHHHHHHHHHHHHHHHHHHHHHHHHHEEEETTEEEECSCCCCCCCHHHHHHHHHHHHHHHHHHHTTCCHHHHHHHHHHHHHHHHHHHHHHHHHHSCCSHHHHHHHHHHHHHHHHHHHHCCCCCHHHHHHHHTCCHHHHHHHHHHHHHHHHHHTHHHHHHHHHHSCCCCHHHHHHHHHHHHHHHHCCCCCHHHHHHHHHHH                      ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                  HHHHHHHHHHHHHHHHHHhccCCcchHHHHHHHHHHHHHHHHHHHHHHHhcccceeccCCcccCCCceeccHHHHHHHHHHHHHHHHHHHHHHccchhccHHHHHHHHHHHHHHHHhHHHHHhCCCCcchhHHHHHHHHHHHHHHHHHHHHHHHHHHHhhcCHHHHHHHHHHHHHHHHHHHHHHHHHHHhcCCCCCCchHHHHHHHHHHHHHhcHHHHHHHHHHHHHh                      ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 18,  41,  41,  66,  67,  78,  82, 114, 118, 176, 178, 211, 212,
                              246],
                             [  7,  30,  31,  56,  56,  67,  67,  99,  99, 157, 157, 190, 190,
                              224]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['L', 'S', 'L', 'T', 'I', 'A', 'G', 'M', 'L', 'A', 'A', 'F', 'V',
              'F', 'F', 'L', 'L', 'A', 'R', 'S', 'Y', 'V', 'A', '-', 'P', 'R',
              'Y', 'H', 'I', 'A', 'L', 'Y', 'L', 'S', 'A', 'L', 'I', 'V', 'F',
              'I', 'A', 'G', 'Y', 'H', 'Y', 'L', 'R', 'I', 'F', 'E', 'S', 'W',
              'V', 'G', 'A', 'Y', 'Q', 'L', 'Q', 'D', 'G', 'V', 'Y', 'V', 'P',
              'T', 'G', 'K', 'P', 'F', 'N', 'D', 'F', 'Y', 'R', 'Y', 'A', 'D',
              'W', 'L', 'L', 'T', 'V', 'P', 'L', 'L', 'L', 'L', 'E', 'L', 'I',
              'L', 'V', 'L', 'G', 'L', 'T', 'A', 'A', 'R', 'T', 'W', 'N', 'L',
              'S', 'I', 'K', 'L', 'V', 'V', 'A', 'S', 'V', 'L', 'M', 'L', 'A',
              'L', 'G', 'Y', 'V', 'G', 'E', 'V', 'N', 'T', 'E', 'P', 'G', 'P',
              'R', 'T', 'L', 'W', 'G', 'A', 'L', 'S', 'S', 'I', 'P', 'F', 'F',
              'Y', 'I', 'L', 'Y', 'V', 'L', 'W', 'V', 'E', 'L', 'G', 'Q', 'A',
              'I', 'R', 'E', 'A', 'K', 'F', 'G', 'P', 'R', 'V', 'L', 'E', 'L',
              'L', 'G', 'A', 'T', 'R', 'L', 'V', 'L', 'L', 'M', 'S', 'W', 'G',
              'F', 'Y', 'P', 'I', 'A', 'Y', 'A', 'L', 'G', 'T', 'W', 'L', 'P',
              'G', 'G', 'A', 'A', 'Q', 'E', 'V', 'A', 'I', 'Q', 'I', 'G', 'Y',
              'S', 'L', 'A', 'D', 'L', 'I', 'A', 'X', 'P', 'I', 'Y', 'G', 'L',
              'L', 'V', 'F', 'A', 'I', 'A', 'R', 'A'],
             ['L', 'A', 'L', 'G', 'T', 'A', 'L', 'M', 'G', 'L', 'G', 'T', 'L',
              'Y', 'F', 'L', 'V', 'K', 'G', 'M', 'G', 'V', 'S', 'D', 'P', 'D',
              'A', 'K', 'K', 'F', 'Y', 'A', 'I', 'T', 'T', 'L', 'V', 'P', 'A',
              'I', 'A', 'F', 'T', 'M', 'Y', 'L', 'S', 'M', 'L', '-', 'L', 'G',
              'Y', 'G', 'L', 'T', 'M', 'V', 'P', 'F', 'G', '-', '-', '-', '-',
              'G', 'E', 'Q', 'N', 'P', 'I', 'Y', 'W', 'A', 'R', 'Y', 'A', 'D',
              'W', 'L', 'F', 'T', 'T', 'P', 'L', 'L', 'L', 'L', 'D', 'L', 'A',
              'L', 'L', 'V', 'D', 'A', 'D', '-', '-', '-', '-', 'Q', 'G', 'T',
              'I', 'L', 'A', 'L', 'V', 'G', 'A', 'D', 'G', 'I', 'M', 'I', 'G',
              'T', 'G', 'L', 'V', 'G', 'A', 'L', 'T', 'K', 'V', 'Y', 'S', 'Y',
              'R', 'F', 'V', 'W', 'W', 'A', 'I', 'S', 'T', 'A', 'A', 'M', 'L',
              'Y', 'I', 'L', 'Y', 'V', 'L', 'F', 'F', 'G', 'F', 'T', 'S', 'K',
              'A', 'E', 'S', '-', '-', 'M', 'R', 'P', 'E', 'V', 'A', 'S', 'T',
              'F', 'K', 'V', 'L', 'R', 'N', 'V', 'T', 'V', 'V', 'L', 'W', 'S',
              'A', 'Y', 'P', 'V', 'V', 'W', 'L', 'I', 'G', 'S', 'E', 'G', '-',
              'A', 'G', 'I', 'V', 'P', 'L', 'N', 'I', 'E', 'T', 'L', 'L', 'F',
              'M', 'V', 'L', 'D', 'V', 'S', 'A', 'K', 'V', 'G', 'F', 'G', 'L',
              'I', 'L', 'L', 'R', 'S', 'R', 'A', 'I']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 100.00)
        self.assertAlmostEqual(alignment.annotations["E-value"], 7e-38)
        self.assertAlmostEqual(alignment.annotations["Score"], 264.24)
        self.assertAlmostEqual(alignment.annotations["Identities"], 22)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.378)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 7.800)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "      HHHHHHHHHHHHHHHHHHHhcCCCChhhHHHHHHHHHHHHHHHHHHHHHHhCCCceeeccCCCCCcchHHHHHHHHHHHHHHHHHHHHHcCCCHHHHHHHHHHHHHHHHHHHHHHHcccchHHHHHHHHHHHHHHHHHHHHHHhHHHHHHhCCHHHHHHHHHHHHHHHHHHHHHHHHHHHccCCCCCCChhHHHHHHHHHHHHHHHHHHHHHHHhhhh  ",
        )
        self.assertEqual(alignment.query.id, "4Y9H:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[6:224],
            "WLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAI",
        )
        self.assertEqual(
            alignment[1],
            "WLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGG---EQNPIYW--ARYADWLFTTPLLLLDLALLVDADQGT----ILALVGADGIMIGTGLVGALTKVY---SYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSE---GAGIVPL---NIETLLFMVLDVSAKVGFGLILLRSRAI",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "      ~~~~~~~~m~~~~~~f~~~~~~~~~~~~r~~~~l~~~i~~ia~~~Y~~ma~~~g~~~~~~~~~~~~v~~~RYv~W~ittPlll~~l~~l~~~~~~~~~~~v~~~~~mi~~g~~g~~~~~~~~~~~~~~~s~~~~~~i~~~l~~~~~~~~~~~~~~~~~~~~~l~~~~~i~W~~YPi~w~l~~~g~~~i~~~~~~~~~~ilDi~~K~~f~~~l~~~~~~  ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                  ~~w~~~~~m~~~~l~f~~~~~~~~~~~R~~~~l~~~i~~va~~~Y~~m~~~~~~~~~~~~g~~~~~~~~~~~~~RYvdW~iTtPlll~~l~llag~~~~~~~~~~~~~v~~~~~mi~~g~~~~~~~~~~~~~~~~~w~~is~~~~~~il~~l~~~~~~~a~~~~~~~~~~~~~l~~~~~v~W~~YPi~w~l~~~~~~~~~~~~~~~~~~e~i~y~ilDv~~K~~fg~~ll~~~~~                    ",
        )
        self.assertEqual(alignment.target.id, "3DDL_A")
        self.assertEqual(
            alignment.target.seq[18:253],
            "MFSFTVATMTASFVFFVLARNNVAPKYRISMMVSALVVFIAGYHYFRITSSWEAAYALQNGMYQPTGELFNDAYRYVDWLLTVPLLTVELVLVMGLPKNERGPLAAKLGFLAALMIVLGYPGEVSENAALFGTRGLWGFLSTIPFVWILYILFTQLGDTIQRQSSRVSTLLGNARLLLLATWGFYPIAYMIPMAFPEAFPSNTPGTIVALQVGYTIADVLAKAGYGVLIYNIAKA",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "3DDL_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Xanthorhodopsin; rhodopsin, carotenoid, ion pump, light-harvesting; HET: PX4, SXN, PCW, RET, UNL; 1.9A {Salinibacter ruber} SCOP: f.13.1.0; Related PDB entries: 3DDL_B",
        )
        self.assertEqual(
            alignment[0],
            "MFSFTVATMTASFVFFVLARNNVA-PKYRISMMVSALVVFIAGYHYFRITSSWEAAYALQNGMYQPTGELFNDAYRYVDWLLTVPLLTVELVLVMGLPKNERGPLAAKLGFLAALMIVLGYPGEVSENAALFGTRGLWGFLSTIPFVWILYILFTQLGDTIQRQSSRVSTLLGNARLLLLATWGFYPIAYMIPMAFPEAFPSNTPGTIVALQVGYTIADVLAKAGYGVLIYNIAKA",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "                  HHHHHHHHHHHHHHHHHHGGGGSCGGGHHHHHHHHHHHHHHHHHHHHHHHHHHHTEEECSSSEEECSCCCCCHHHHHHHHHHHHHHHHHHHHHSCCCHHHHHHHHHHHHHHHHHHHHHHHHHHTCSCCCTTSHHHHHHHHHHHHHHHHHHHHHHSCHHHHTTSCHHHHHHHHHHHHHHHHHHHHHHHHHHTTCCTTCHHHHHHHHHHHHHHHHCCCCCHHHHHHHHHHH                          ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                  HHHHHHHHHHHHHHHHHHHhcCCCcchhHHHHHHHHHHHHHHHHHHHHHhchHhHHHhhCCCCCCCCcccchHHHHHHHHHHHHHHHHHHHHHhCCCcccchhHHHHHHHHHHHHHHhchHHHhcccccccchhHHHHHHHHHHHHHHHHHHHHHHHHHHHhCCHHHHHHHHHHHHHHHHHHHHHHHHHHhhhcCcccCCCCCCCcHHHHHHHHHHHHHHHHHHHHHHHHHHHHH                    ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 18,  42,  42,  79,  82,  89,  91, 118, 122, 146, 149, 212, 215,
                              222, 225, 253],
                             [  6,  30,  31,  68,  68,  75,  75, 102, 102, 126, 126, 189, 189,
                              196, 196, 224]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['M', 'F', 'S', 'F', 'T', 'V', 'A', 'T', 'M', 'T', 'A', 'S', 'F',
              'V', 'F', 'F', 'V', 'L', 'A', 'R', 'N', 'N', 'V', 'A', '-', 'P',
              'K', 'Y', 'R', 'I', 'S', 'M', 'M', 'V', 'S', 'A', 'L', 'V', 'V',
              'F', 'I', 'A', 'G', 'Y', 'H', 'Y', 'F', 'R', 'I', 'T', 'S', 'S',
              'W', 'E', 'A', 'A', 'Y', 'A', 'L', 'Q', 'N', 'G', 'M', 'Y', 'Q',
              'P', 'T', 'G', 'E', 'L', 'F', 'N', 'D', 'A', 'Y', 'R', 'Y', 'V',
              'D', 'W', 'L', 'L', 'T', 'V', 'P', 'L', 'L', 'T', 'V', 'E', 'L',
              'V', 'L', 'V', 'M', 'G', 'L', 'P', 'K', 'N', 'E', 'R', 'G', 'P',
              'L', 'A', 'A', 'K', 'L', 'G', 'F', 'L', 'A', 'A', 'L', 'M', 'I',
              'V', 'L', 'G', 'Y', 'P', 'G', 'E', 'V', 'S', 'E', 'N', 'A', 'A',
              'L', 'F', 'G', 'T', 'R', 'G', 'L', 'W', 'G', 'F', 'L', 'S', 'T',
              'I', 'P', 'F', 'V', 'W', 'I', 'L', 'Y', 'I', 'L', 'F', 'T', 'Q',
              'L', 'G', 'D', 'T', 'I', 'Q', 'R', 'Q', 'S', 'S', 'R', 'V', 'S',
              'T', 'L', 'L', 'G', 'N', 'A', 'R', 'L', 'L', 'L', 'L', 'A', 'T',
              'W', 'G', 'F', 'Y', 'P', 'I', 'A', 'Y', 'M', 'I', 'P', 'M', 'A',
              'F', 'P', 'E', 'A', 'F', 'P', 'S', 'N', 'T', 'P', 'G', 'T', 'I',
              'V', 'A', 'L', 'Q', 'V', 'G', 'Y', 'T', 'I', 'A', 'D', 'V', 'L',
              'A', 'K', 'A', 'G', 'Y', 'G', 'V', 'L', 'I', 'Y', 'N', 'I', 'A',
              'K', 'A'],
             ['W', 'L', 'A', 'L', 'G', 'T', 'A', 'L', 'M', 'G', 'L', 'G', 'T',
              'L', 'Y', 'F', 'L', 'V', 'K', 'G', 'M', 'G', 'V', 'S', 'D', 'P',
              'D', 'A', 'K', 'K', 'F', 'Y', 'A', 'I', 'T', 'T', 'L', 'V', 'P',
              'A', 'I', 'A', 'F', 'T', 'M', 'Y', 'L', 'S', 'M', 'L', 'L', 'G',
              'Y', 'G', 'L', 'T', 'M', 'V', 'P', 'F', 'G', 'G', '-', '-', '-',
              'E', 'Q', 'N', 'P', 'I', 'Y', 'W', '-', '-', 'A', 'R', 'Y', 'A',
              'D', 'W', 'L', 'F', 'T', 'T', 'P', 'L', 'L', 'L', 'L', 'D', 'L',
              'A', 'L', 'L', 'V', 'D', 'A', 'D', 'Q', 'G', 'T', '-', '-', '-',
              '-', 'I', 'L', 'A', 'L', 'V', 'G', 'A', 'D', 'G', 'I', 'M', 'I',
              'G', 'T', 'G', 'L', 'V', 'G', 'A', 'L', 'T', 'K', 'V', 'Y', '-',
              '-', '-', 'S', 'Y', 'R', 'F', 'V', 'W', 'W', 'A', 'I', 'S', 'T',
              'A', 'A', 'M', 'L', 'Y', 'I', 'L', 'Y', 'V', 'L', 'F', 'F', 'G',
              'F', 'T', 'S', 'K', 'A', 'E', 'S', 'M', 'R', 'P', 'E', 'V', 'A',
              'S', 'T', 'F', 'K', 'V', 'L', 'R', 'N', 'V', 'T', 'V', 'V', 'L',
              'W', 'S', 'A', 'Y', 'P', 'V', 'V', 'W', 'L', 'I', 'G', 'S', 'E',
              '-', '-', '-', 'G', 'A', 'G', 'I', 'V', 'P', 'L', '-', '-', '-',
              'N', 'I', 'E', 'T', 'L', 'L', 'F', 'M', 'V', 'L', 'D', 'V', 'S',
              'A', 'K', 'V', 'G', 'F', 'G', 'L', 'I', 'L', 'L', 'R', 'S', 'R',
              'A', 'I']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 100.00)
        self.assertAlmostEqual(alignment.annotations["E-value"], 6.1e-37)
        self.assertAlmostEqual(alignment.annotations["Score"], 263.90)
        self.assertAlmostEqual(alignment.annotations["Identities"], 20)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.322)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 6.600)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "    HHHHHHHHHHHHHHHHHHHHHhcCCCChhhHHHHHHHHHHHHHHHHHHHHHHhCCCceeeccCCCCCcchHHHHHHHHHHHHHHHHHHHHHcCCCHHHHHHHHHHHHHHHHHHHHHHHcccchHHHHHHHHHHHHHHHHHHHHHHhHHHHHHhCCHHHHHHHHHHHHHHHHHHHHHHHHHHHccCCCCCCChhHHHHHHHHHHHHHHHHHHHHHHHhhhh  ",
        )
        self.assertEqual(alignment.query.id, "4Y9H:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[4:224],
            "WIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAI",
        )
        self.assertEqual(
            alignment[1],
            "WIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTL-VPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQG---TILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMR-PEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAI",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "    ~~~~~~~~~~m~~~~~~f~~~~~~~~~~~~r~~~~l~~~i~~ia~~~Y~~ma~~~g~~~~~~~~~~~~v~~~RYv~W~ittPlll~~l~~l~~~~~~~~~~~v~~~~~mi~~g~~g~~~~~~~~~~~~~~~s~~~~~~i~~~l~~~~~~~~~~~~~~~~~~~~~l~~~~~i~W~~YPi~w~l~~~g~~~i~~~~~~~~~~ilDi~~K~~f~~~l~~~~~~  ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                                                                    ~~~~W~~f~v~~l~~l~f~~~~~r~~~~~r~~y~~~~~~i~~iiaY~~ma~~~g~~~~~~gr~v~~~RYvdWllTtPLlL~~L~llag~~~~~~~~i~~li~~d~~MIvtGl~gal~~~~~kw~~f~ig~~~~l~~~~~~~~~~~~~a~~~~~~~~~~~~~~L~~~~~v~W~~YPIvwlL~~eG~~~is~~~e~i~y~ilDilaK~~fg~ll~~~r~~                                               ",
        )
        self.assertEqual(alignment.target.id, "3UG9_A")
        self.assertEqual(
            alignment.target.seq[68:286],
            "NILQWITFALSALCLMFYGYQTWKSTCGWEEIYVATIEMIKFIIEYFHEFDEPAVIYSSNGNKTVWLRYAEWLLTCPVILIHLSNLTGLANDYNKRTMGLLVSDIGTIVWGTTAALSKGYVRVIFFLMGLCYGIYTFFNAAKVYIEAYHTVPKGRCRQVVTGMAWLFFVSWGMFPILFILGPEGFGVLSVYGSTVGHTIIDLMSKNCWGLLGHYLRVL",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "3UG9_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Archaeal-type opsin 1, Archaeal-type opsin; microbialrhodopsin, seven-transmembrane, light-gated cation channel; HET: OLA, RET; 2.3A {Chlamydomonas reinhardtii}; Related PDB entries: 6CSO_A 6EIG_A 6EIG_B 6EID_B 6EID_A 4YZI_A 6CSN_A",
        )
        self.assertEqual(
            alignment[0],
            "NILQWITFALSALCLMFYGYQTWKST-CGWEEIYVATIEMIKFI--IEYFHEFDEPAVIYSS---NGNKTVWLRYAEWLLTCPVILIHLSNLTGLANDYNKRTMGLLVSDIGTIVWGTTAALSKG-YVRVIFFLMGLCYGIYTFFNAAKVYIEAYHTVPKGRCRQVVTGMAWLFFVSWGMFPILFILGPEGFGVLSVYGSTVGHTIIDLMSKNCWGLLGHYLRVL",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "                                                                    HHHHHHHHHHHHHHHHHCCCHHHHHHHHHHHHHHHHHHHHTTSTTTSEECTTSCEECHHHHHHHHHHHHHHHHHHTTTTSCCCCCCHHHHHHHHHHHHHHHHHHHHHHSCHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHSCSSHHHHHHHHHHHHHHHHHHHHHHHHHHSTTTTCSSCHHHHHHHHHHHHCCCCCCHHHHHHHHHHH                                                       ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                                                                    HHHHHHHHHHHHHHHHHHHHHHhcCCCCchHHHHHHHHHHHHHHHHHHHhcCCceeecCCCceeeHHHHHHHHhHHHHHHHHHHHHhCCCCcchHHHHHHHHHHHHHHHHHHHHHHcccchhHHHHHHHHHHHHHHHHHHHHHHHHHHhcCCChHHHHHHHHHHHHHHHHHHHHHHHHHHCccCCCCCCHHHHHHHHHHHHHHHHHHHHHHHHHHHHH                                               ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 68,  94,  94, 106, 107, 111, 111, 127, 127, 160, 163, 187, 187,
                              220, 221, 286],
                             [  4,  30,  31,  43,  43,  47,  49,  65,  68, 101, 101, 125, 126,
                              159, 159, 224]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['N', 'I', 'L', 'Q', 'W', 'I', 'T', 'F', 'A', 'L', 'S', 'A', 'L',
              'C', 'L', 'M', 'F', 'Y', 'G', 'Y', 'Q', 'T', 'W', 'K', 'S', 'T',
              '-', 'C', 'G', 'W', 'E', 'E', 'I', 'Y', 'V', 'A', 'T', 'I', 'E',
              'M', 'I', 'K', 'F', 'I', '-', '-', 'I', 'E', 'Y', 'F', 'H', 'E',
              'F', 'D', 'E', 'P', 'A', 'V', 'I', 'Y', 'S', 'S', '-', '-', '-',
              'N', 'G', 'N', 'K', 'T', 'V', 'W', 'L', 'R', 'Y', 'A', 'E', 'W',
              'L', 'L', 'T', 'C', 'P', 'V', 'I', 'L', 'I', 'H', 'L', 'S', 'N',
              'L', 'T', 'G', 'L', 'A', 'N', 'D', 'Y', 'N', 'K', 'R', 'T', 'M',
              'G', 'L', 'L', 'V', 'S', 'D', 'I', 'G', 'T', 'I', 'V', 'W', 'G',
              'T', 'T', 'A', 'A', 'L', 'S', 'K', 'G', '-', 'Y', 'V', 'R', 'V',
              'I', 'F', 'F', 'L', 'M', 'G', 'L', 'C', 'Y', 'G', 'I', 'Y', 'T',
              'F', 'F', 'N', 'A', 'A', 'K', 'V', 'Y', 'I', 'E', 'A', 'Y', 'H',
              'T', 'V', 'P', 'K', 'G', 'R', 'C', 'R', 'Q', 'V', 'V', 'T', 'G',
              'M', 'A', 'W', 'L', 'F', 'F', 'V', 'S', 'W', 'G', 'M', 'F', 'P',
              'I', 'L', 'F', 'I', 'L', 'G', 'P', 'E', 'G', 'F', 'G', 'V', 'L',
              'S', 'V', 'Y', 'G', 'S', 'T', 'V', 'G', 'H', 'T', 'I', 'I', 'D',
              'L', 'M', 'S', 'K', 'N', 'C', 'W', 'G', 'L', 'L', 'G', 'H', 'Y',
              'L', 'R', 'V', 'L'],
             ['W', 'I', 'W', 'L', 'A', 'L', 'G', 'T', 'A', 'L', 'M', 'G', 'L',
              'G', 'T', 'L', 'Y', 'F', 'L', 'V', 'K', 'G', 'M', 'G', 'V', 'S',
              'D', 'P', 'D', 'A', 'K', 'K', 'F', 'Y', 'A', 'I', 'T', 'T', 'L',
              '-', 'V', 'P', 'A', 'I', 'A', 'F', 'T', 'M', 'Y', 'L', 'S', 'M',
              'L', 'L', 'G', 'Y', 'G', 'L', 'T', 'M', 'V', 'P', 'F', 'G', 'G',
              'E', 'Q', 'N', 'P', 'I', 'Y', 'W', 'A', 'R', 'Y', 'A', 'D', 'W',
              'L', 'F', 'T', 'T', 'P', 'L', 'L', 'L', 'L', 'D', 'L', 'A', 'L',
              'L', 'V', 'D', 'A', 'D', 'Q', 'G', '-', '-', '-', 'T', 'I', 'L',
              'A', 'L', 'V', 'G', 'A', 'D', 'G', 'I', 'M', 'I', 'G', 'T', 'G',
              'L', 'V', 'G', 'A', 'L', 'T', 'K', 'V', 'Y', 'S', 'Y', 'R', 'F',
              'V', 'W', 'W', 'A', 'I', 'S', 'T', 'A', 'A', 'M', 'L', 'Y', 'I',
              'L', 'Y', 'V', 'L', 'F', 'F', 'G', 'F', 'T', 'S', 'K', 'A', 'E',
              'S', 'M', 'R', '-', 'P', 'E', 'V', 'A', 'S', 'T', 'F', 'K', 'V',
              'L', 'R', 'N', 'V', 'T', 'V', 'V', 'L', 'W', 'S', 'A', 'Y', 'P',
              'V', 'V', 'W', 'L', 'I', 'G', 'S', 'E', 'G', 'A', 'G', 'I', 'V',
              'P', 'L', 'N', 'I', 'E', 'T', 'L', 'L', 'F', 'M', 'V', 'L', 'D',
              'V', 'S', 'A', 'K', 'V', 'G', 'F', 'G', 'L', 'I', 'L', 'L', 'R',
              'S', 'R', 'A', 'I']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 100.00)
        self.assertAlmostEqual(alignment.annotations["E-value"], 1e-36)
        self.assertAlmostEqual(alignment.annotations["Score"], 249.63)
        self.assertAlmostEqual(alignment.annotations["Identities"], 23)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.407)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 8.300)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "     HHHHHHHHHHHHHHHHHHHHhcCCCChhhHHHHHHHHHHHHHHHHHHHHHHhCCCceeeccCCCCCcchHHHHHHHHHHHHHHHHHHHHHcCCCHHHHHHHHHHHHHHHHHHHHHHHcccchHHHHHHHHHHHHHHHHHHHHHHhHHHHHHhCCHHHHHHHHHHHHHHHHHHHHHHHHHHHccCCCCCCChhHHHHHHHHHHHHHHHHHHHHHHHhhhh  ",
        )
        self.assertEqual(alignment.query.id, "4Y9H:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[5:224],
            "IWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAI",
        )
        self.assertEqual(
            alignment[1],
            "IWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQG---TILALVGADGIMIGTGLVGALTKV-YSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAI",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "     ~~~~~~~~~m~~~~~~f~~~~~~~~~~~~r~~~~l~~~i~~ia~~~Y~~ma~~~g~~~~~~~~~~~~v~~~RYv~W~ittPlll~~l~~l~~~~~~~~~~~v~~~~~mi~~g~~g~~~~~~~~~~~~~~~s~~~~~~i~~~l~~~~~~~~~~~~~~~~~~~~~l~~~~~i~W~~YPi~w~l~~~g~~~i~~~~~~~~~~ilDi~~K~~f~~~l~~~~~~  ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "      ~~~w~~~~~~~~~~~~f~~~~~~~~~~~R~~~~l~~~i~~ia~~~Y~~m~~~~~~g~~~~~~RYv~W~~TtPlll~~l~~l~g~~~~~~~~~~~~~~~~~~mi~~g~~g~~~~~~~~~~w~~~~~s~~~~~~ll~~l~~~~~~~~~~~~~~~~~l~~~~~~~W~~YPivw~l~~~g~~~~~~~~~~~~~Di~~K~~f~~~ll~~~~~       ",
        )
        self.assertEqual(alignment.target.id, "5JSI_A")
        self.assertEqual(
            alignment.target.seq[6:213],
            "RLFMVATVGMLAGTVFLLASSREVKPEHRRGVYISALVCGIAWYHYQKMGASWESGSYDTGLRYVDWVLTVPLMFVEVLAVTRKGAAYNEAVRNWGIAATVMIGAGYYGETSAAGSNEYWTGFVIAMATYVWLMRNLQAEGEGLKGDQAVAFENIKNLILVGWIIYPLGYIAPVVGDFDAIREVLYTIADIINXVGLGVLVLQMARV",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "5JSI_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Nitrate/nitrite sensor protein NarQ (E.C.2.7.13.3); membrane protein, iodide, UNKNOWN FUNCTION; HET: LYR, OLC, LFA, IOD; 2.0A {Candidatus Actinomarina minuta}; Related PDB entries: 5JSI_B",
        )
        self.assertEqual(
            alignment[0],
            "RLFMVATVGMLAGTVFLLASSREVK-PEHRRGVYISALVCGIAWYHYQKMGASWE--------SGSYDTGLRYVDWVLTVPLMFVEVLAVTRKGAAYNEAVRNWGIAATVMIGAGYYGETSAAGSNEYWTGFVIAMATYVWLMRNLQA----EGEGLKGDQAVAFENIKNLILVGWIIYPLGYIAPVVG--DFD-AIREVLYTIADIINXVGLGVLVLQMARV",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "      HHHHHHHHHHHHHHHHHHHHGGGSCGGGCHHHHHHHHHHHHHHHHHHHHHHHHHHSSCCTHHHHHHHHHHHHHHHHHHHHHHCCHHHHHHHHHHHHHHHHHHHHHHHHHHTSCTTSHHHHHHHHHHHHHHHHHHHHHHHTTTTCCHHHHHHHHHHHHHHHHHTTHHHHHHHHHHHSCCHHHHHHHHHHHHHHHHHHHHHHHHHHHHH       ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "      HHHHHHHHHHHHHHHHHHHHcCCCChHHHHHHHHHHHHHHHHHHHHHHHHHHhhcCCCCcHHHHhHHHhHHHHHHHHHHHHHhcccchHHHHHHHHHHHHHHHHHHHHHHHhcCCccchHHHHHHHHHHHHHHHHHHHHhhccCCHHHHHHHHHHHHHHHHHHHHHHHHHHhhhcCChHHHHHHHHHHHHHHhHHHHHHHHHHHHHH       ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  6,  31,  31,  60,  60,  93,  96, 120, 121, 145, 145, 182, 182,
                              185, 185, 213],
                             [  5,  30,  31,  60,  68, 101, 101, 125, 125, 149, 153, 190, 192,
                              195, 196, 224]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['R', 'L', 'F', 'M', 'V', 'A', 'T', 'V', 'G', 'M', 'L', 'A', 'G',
              'T', 'V', 'F', 'L', 'L', 'A', 'S', 'S', 'R', 'E', 'V', 'K', '-',
              'P', 'E', 'H', 'R', 'R', 'G', 'V', 'Y', 'I', 'S', 'A', 'L', 'V',
              'C', 'G', 'I', 'A', 'W', 'Y', 'H', 'Y', 'Q', 'K', 'M', 'G', 'A',
              'S', 'W', 'E', '-', '-', '-', '-', '-', '-', '-', '-', 'S', 'G',
              'S', 'Y', 'D', 'T', 'G', 'L', 'R', 'Y', 'V', 'D', 'W', 'V', 'L',
              'T', 'V', 'P', 'L', 'M', 'F', 'V', 'E', 'V', 'L', 'A', 'V', 'T',
              'R', 'K', 'G', 'A', 'A', 'Y', 'N', 'E', 'A', 'V', 'R', 'N', 'W',
              'G', 'I', 'A', 'A', 'T', 'V', 'M', 'I', 'G', 'A', 'G', 'Y', 'Y',
              'G', 'E', 'T', 'S', 'A', 'A', 'G', 'S', 'N', 'E', 'Y', 'W', 'T',
              'G', 'F', 'V', 'I', 'A', 'M', 'A', 'T', 'Y', 'V', 'W', 'L', 'M',
              'R', 'N', 'L', 'Q', 'A', '-', '-', '-', '-', 'E', 'G', 'E', 'G',
              'L', 'K', 'G', 'D', 'Q', 'A', 'V', 'A', 'F', 'E', 'N', 'I', 'K',
              'N', 'L', 'I', 'L', 'V', 'G', 'W', 'I', 'I', 'Y', 'P', 'L', 'G',
              'Y', 'I', 'A', 'P', 'V', 'V', 'G', '-', '-', 'D', 'F', 'D', '-',
              'A', 'I', 'R', 'E', 'V', 'L', 'Y', 'T', 'I', 'A', 'D', 'I', 'I',
              'N', 'X', 'V', 'G', 'L', 'G', 'V', 'L', 'V', 'L', 'Q', 'M', 'A',
              'R', 'V'],
             ['I', 'W', 'L', 'A', 'L', 'G', 'T', 'A', 'L', 'M', 'G', 'L', 'G',
              'T', 'L', 'Y', 'F', 'L', 'V', 'K', 'G', 'M', 'G', 'V', 'S', 'D',
              'P', 'D', 'A', 'K', 'K', 'F', 'Y', 'A', 'I', 'T', 'T', 'L', 'V',
              'P', 'A', 'I', 'A', 'F', 'T', 'M', 'Y', 'L', 'S', 'M', 'L', 'L',
              'G', 'Y', 'G', 'L', 'T', 'M', 'V', 'P', 'F', 'G', 'G', 'E', 'Q',
              'N', 'P', 'I', 'Y', 'W', 'A', 'R', 'Y', 'A', 'D', 'W', 'L', 'F',
              'T', 'T', 'P', 'L', 'L', 'L', 'L', 'D', 'L', 'A', 'L', 'L', 'V',
              'D', 'A', 'D', 'Q', 'G', '-', '-', '-', 'T', 'I', 'L', 'A', 'L',
              'V', 'G', 'A', 'D', 'G', 'I', 'M', 'I', 'G', 'T', 'G', 'L', 'V',
              'G', 'A', 'L', 'T', 'K', 'V', '-', 'Y', 'S', 'Y', 'R', 'F', 'V',
              'W', 'W', 'A', 'I', 'S', 'T', 'A', 'A', 'M', 'L', 'Y', 'I', 'L',
              'Y', 'V', 'L', 'F', 'F', 'G', 'F', 'T', 'S', 'K', 'A', 'E', 'S',
              'M', 'R', 'P', 'E', 'V', 'A', 'S', 'T', 'F', 'K', 'V', 'L', 'R',
              'N', 'V', 'T', 'V', 'V', 'L', 'W', 'S', 'A', 'Y', 'P', 'V', 'V',
              'W', 'L', 'I', 'G', 'S', 'E', 'G', 'A', 'G', 'I', 'V', 'P', 'L',
              'N', 'I', 'E', 'T', 'L', 'L', 'F', 'M', 'V', 'L', 'D', 'V', 'S',
              'A', 'K', 'V', 'G', 'F', 'G', 'L', 'I', 'L', 'L', 'R', 'S', 'R',
              'A', 'I']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 100.00)
        self.assertAlmostEqual(alignment.annotations["E-value"], 2.7e-36)
        self.assertAlmostEqual(alignment.annotations["Score"], 255.00)
        self.assertAlmostEqual(alignment.annotations["Identities"], 21)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.325)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 7.700)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "CCcHHHHHHHHHHHHHHHHHHHHHHhcCCCChhhHHHHHHHHHHHHHHHHHHHHHHhCCCceeeccCCCCCcchHHHHHHHHHHHHHHHHHHHHHcCCCHHHHHHHHHHHHHHHHHHHHHHHcccchHHHHHHHHHHHHHHHHHHHHHHhHHHHHHhCCHHHHHHHHHHHHHHHHHHHHHHHHHHHccCCCCCCChhHHHHHHHHHHHHHHHHHHHHHHHhhhh  ",
        )
        self.assertEqual(alignment.query.id, "4Y9H:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[0:224],
            "GRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAI",
        )
        self.assertEqual(
            alignment[1],
            "GRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQ---GTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFG-------FTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAI",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "~~~~~~~~~~~~~~m~~~~~~f~~~~~~~~~~~~r~~~~l~~~i~~ia~~~Y~~ma~~~g~~~~~~~~~~~~v~~~RYv~W~ittPlll~~l~~l~~~~~~~~~~~v~~~~~mi~~g~~g~~~~~~~~~~~~~~~s~~~~~~i~~~l~~~~~~~~~~~~~~~~~~~~~l~~~~~i~W~~YPi~w~l~~~g~~~i~~~~~~~~~~ilDi~~K~~f~~~l~~~~~~  ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                          ~g~~~~w~v~~~~~~~~~~~~~sr~~~~~~l~~~i~~va~~~Y~~ma~g~g~~~~~~~r~v~~~RYv~W~vTtPlll~~l~llag~~~~~~~~~~~~i~~~~~mi~~G~~g~~~~~~~~kw~~~~is~~~fl~v~~~l~~~~~~~~~~~~~~a~~~~~~~~~~~~~l~~~~~v~W~~YPi~w~l~~~G~~~is~~~e~i~y~ilDil~K~~f~~~ll~~~~~                             ",
        )
        self.assertEqual(alignment.target.id, "6CSM_D")
        self.assertEqual(
            alignment.target.seq[26:248],
            "DGIKYVQLVMAVVSACQVFFMVTRAPKVPWEAIYLPTTEMITYSLAFTGNGYIRVANGKYLPWARMASWLCTCPIMLGLVSNMALVKYKSIPLNPMMIAASSICTVFGITASVVLDPLHVWLYCFISSIFFIFEMVVAFAIFAITIHDFQTIGSPMSLKVVERLKLMRIVFYVSWMAYPILWSFSSTGACIMSENTSSVLYLLGDALCKNTYGILLWATTWG",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "6CSM_D")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "GtACR1; rhodopsin, channelrhodopsin, anion channel, optogenetics; HET: RET, OLA; 2.9A {Guillardia theta CCMP2712}; Related PDB entries: 6CSM_B 6CSM_A 6CSM_C",
        )
        self.assertEqual(
            alignment[0],
            "DGIKYVQL---VMAVVSACQVFFMVTRAPK------VPWEAIYLPTTEMITYSLAFTGNGYIRVA---NGKYLPWARMASWLCTCPIMLGLVSNMALVKYKSIPLNPMMIAASSICTVFGITASVVLDPLHVWLYCFISSIFFIFEMVVAFAIFAITIHDFQTIGSPMSLKVVERLKLMRIVFYVSWMAYPILWSFSSTGACIMSENTSSVLYLLGDALCKNTYGILLWATTWG",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "                          HHHHHHHHHHHHHHHHHHHHHHHTCSSCCTHHHHHHHHHHHHHHHHHTTCCEEEBTTSCEEEHHHHHHHHHHHHHHHHHHHTTCCCEETTEECHHHHHHHHHHHHHHHHHHTTCSCHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHSCCHHHHHHHHHHHHHHHHHHHHHTHHHHHHHHSTTTTCCSCHHHHHHHHHHHHCCCCCCHHHHHHCCCCC                             ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                          HHHHHHHHHHHHHHHHHHHHHHhcCCCCcHHHHHHHHHHHHHHHHHHhCCCceecCCCccchHHHHHHHHHHHHHHHHHHHHHhcCCcccCChHHHHHHHHHHHHHHHHHHHhcCCchHHHHHHHHHHHHHHHHHHHHHHHHHHHHhhHHhhcccCCHHHHHHHHHHHHHHHHHHHHHHHHHHHCccCCCCCCHHHHHHHHHHHHHHHHHHHHHHHHHHHHh                             ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 26,  34,  34,  53,  53,  82,  82, 114, 117, 167, 174, 248],
                             [  0,   8,  11,  30,  36,  65,  68, 100, 100, 150, 150, 224]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['D', 'G', 'I', 'K', 'Y', 'V', 'Q', 'L', '-', '-', '-', 'V', 'M',
              'A', 'V', 'V', 'S', 'A', 'C', 'Q', 'V', 'F', 'F', 'M', 'V', 'T',
              'R', 'A', 'P', 'K', '-', '-', '-', '-', '-', '-', 'V', 'P', 'W',
              'E', 'A', 'I', 'Y', 'L', 'P', 'T', 'T', 'E', 'M', 'I', 'T', 'Y',
              'S', 'L', 'A', 'F', 'T', 'G', 'N', 'G', 'Y', 'I', 'R', 'V', 'A',
              '-', '-', '-', 'N', 'G', 'K', 'Y', 'L', 'P', 'W', 'A', 'R', 'M',
              'A', 'S', 'W', 'L', 'C', 'T', 'C', 'P', 'I', 'M', 'L', 'G', 'L',
              'V', 'S', 'N', 'M', 'A', 'L', 'V', 'K', 'Y', 'K', 'S', 'I', 'P',
              'L', 'N', 'P', 'M', 'M', 'I', 'A', 'A', 'S', 'S', 'I', 'C', 'T',
              'V', 'F', 'G', 'I', 'T', 'A', 'S', 'V', 'V', 'L', 'D', 'P', 'L',
              'H', 'V', 'W', 'L', 'Y', 'C', 'F', 'I', 'S', 'S', 'I', 'F', 'F',
              'I', 'F', 'E', 'M', 'V', 'V', 'A', 'F', 'A', 'I', 'F', 'A', 'I',
              'T', 'I', 'H', 'D', 'F', 'Q', 'T', 'I', 'G', 'S', 'P', 'M', 'S',
              'L', 'K', 'V', 'V', 'E', 'R', 'L', 'K', 'L', 'M', 'R', 'I', 'V',
              'F', 'Y', 'V', 'S', 'W', 'M', 'A', 'Y', 'P', 'I', 'L', 'W', 'S',
              'F', 'S', 'S', 'T', 'G', 'A', 'C', 'I', 'M', 'S', 'E', 'N', 'T',
              'S', 'S', 'V', 'L', 'Y', 'L', 'L', 'G', 'D', 'A', 'L', 'C', 'K',
              'N', 'T', 'Y', 'G', 'I', 'L', 'L', 'W', 'A', 'T', 'T', 'W', 'G'],
             ['G', 'R', 'P', 'E', 'W', 'I', 'W', 'L', 'A', 'L', 'G', 'T', 'A',
              'L', 'M', 'G', 'L', 'G', 'T', 'L', 'Y', 'F', 'L', 'V', 'K', 'G',
              'M', 'G', 'V', 'S', 'D', 'P', 'D', 'A', 'K', 'K', 'F', 'Y', 'A',
              'I', 'T', 'T', 'L', 'V', 'P', 'A', 'I', 'A', 'F', 'T', 'M', 'Y',
              'L', 'S', 'M', 'L', 'L', 'G', 'Y', 'G', 'L', 'T', 'M', 'V', 'P',
              'F', 'G', 'G', 'E', 'Q', 'N', 'P', 'I', 'Y', 'W', 'A', 'R', 'Y',
              'A', 'D', 'W', 'L', 'F', 'T', 'T', 'P', 'L', 'L', 'L', 'L', 'D',
              'L', 'A', 'L', 'L', 'V', 'D', 'A', 'D', 'Q', '-', '-', '-', 'G',
              'T', 'I', 'L', 'A', 'L', 'V', 'G', 'A', 'D', 'G', 'I', 'M', 'I',
              'G', 'T', 'G', 'L', 'V', 'G', 'A', 'L', 'T', 'K', 'V', 'Y', 'S',
              'Y', 'R', 'F', 'V', 'W', 'W', 'A', 'I', 'S', 'T', 'A', 'A', 'M',
              'L', 'Y', 'I', 'L', 'Y', 'V', 'L', 'F', 'F', 'G', '-', '-', '-',
              '-', '-', '-', '-', 'F', 'T', 'S', 'K', 'A', 'E', 'S', 'M', 'R',
              'P', 'E', 'V', 'A', 'S', 'T', 'F', 'K', 'V', 'L', 'R', 'N', 'V',
              'T', 'V', 'V', 'L', 'W', 'S', 'A', 'Y', 'P', 'V', 'V', 'W', 'L',
              'I', 'G', 'S', 'E', 'G', 'A', 'G', 'I', 'V', 'P', 'L', 'N', 'I',
              'E', 'T', 'L', 'L', 'F', 'M', 'V', 'L', 'D', 'V', 'S', 'A', 'K',
              'V', 'G', 'F', 'G', 'L', 'I', 'L', 'L', 'R', 'S', 'R', 'A', 'I']],
            dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 100.00)
        self.assertAlmostEqual(alignment.annotations["E-value"], 5.8e-36)
        self.assertAlmostEqual(alignment.annotations["Score"], 247.02)
        self.assertAlmostEqual(alignment.annotations["Identities"], 24)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.317)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 7.900)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "CCcHHHHHHHHHHHHHHHHHHHHHHhcCCCChhhHHHHHHHHHHHHHHHHHHHHHHhCCCceeeccCCCCCcchHHHHHHHHHHHHHHHHHHHHHcCCCHHHHHHHHHHHHHHHHHHHHHHHcccchHHHHHHHHHHHHHHHHHHHHHHhHHHHHHhCCHHHHHHHHHHHHHHHHHHHHHHHHHHHccCCCCCCChhHHHHHHHHHHHHHHHHHHHHHHHhhhh  ",
        )
        self.assertEqual(alignment.query.id, "4Y9H:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[0:224],
            "GRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAI",
        )
        self.assertEqual(
            alignment[1],
            "GRPEWIWLALGTAL-MGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGT-----ILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIG-----SEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAI",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "~~~~~~~~~~~~~~m~~~~~~f~~~~~~~~~~~~r~~~~l~~~i~~ia~~~Y~~ma~~~g~~~~~~~~~~~~v~~~RYv~W~ittPlll~~l~~l~~~~~~~~~~~v~~~~~mi~~g~~g~~~~~~~~~~~~~~~s~~~~~~i~~~l~~~~~~~~~~~~~~~~~~~~~l~~~~~i~W~~YPi~w~l~~~g~~~i~~~~~~~~~~ilDi~~K~~f~~~l~~~~~~  ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "        ~~~g~~~lw~~~~i~~~~~~~~f~~~~~~~~~~r~~~~l~~~i~~ia~~~Y~~m~~~g~~~~~~RYi~W~~TtPlll~~l~~l~~~~~~~~~~~~~~~l~~~~~~mi~~G~~~~~~~~~~~~~~~~~~f~~~~~~l~~~~~~~a~~~~~~~~~~~~~l~~~~~v~W~~YPi~w~l~~~~~~~~g~~~i~~~~~~~~~~ilDi~~K~~f~~~ll~~~~~      ",
        )
        self.assertEqual(alignment.target.id, "6EYU_B")
        self.assertEqual(
            alignment.target.seq[8:226],
            "GGFGSQPFILAYIITAMISGLLFLYLPRKLDVPQKFGIIHFFIVVWSGLMYTNFLNQSFLSDYAWYMDWMVSTPLILLALGLTAFHGADTKRYDLLGALLGAEFTLVITGLLAQAQGSITPYYVGVLLLLGVVYLLAKPFREIAEESSDGLARAYKILAGYIGIFFLSYPTVWYISGIDALPGSLNILDPTQTSIALVVLPFFCKQVYGFLDMYLIHK",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "6EYU_B")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Bacteriorhodopsin; Retinal protein, proton transport, MEMBRANE; HET: RET, LFA, MUN; 2.5A {Nanosalina sp. (strain J07AB43)}; Related PDB entries: 6EYU_C 6EYU_A",
        )
        self.assertEqual(
            alignment[0],
            "GGFGSQPFILAYIITAMISGLLFLYLPRKL--DVPQKFGIIHFFIVVWSGLMYTNFLNQ-----------SFLSDYAWYMDWMVSTPLILLALGLTAFHGADTKRYDLLGALLGAEFTLVITGLLAQAQGS----ITPYYVGVLLLLGVVYLLAKPFREIAEESSDGLARAYKILAGYIGIFFLSYPTVWYISGIDALPGSLNILDPTQTSIALVVLPFFCKQVYGFLDMYLIHK",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "        TTCTTHHHHHHHHHHHHHHHHHHHHHHHHHTCCHHHHHHHHHHHHHHHHHHTTTTSCCTTGGGHHHHHHHHHHHHHHHHHHHHHHHTCSCCCHHHHHHHHHHHHHHHHHHHHHHHTTCSHHHHHHHHHHHHHHHHHHTHHHHHHTTSCHHHHHHHHHHHHHHHHHHHHHHHHHHHBCCSSSCSSCCCBCHHHHHHHHHHHHHCCCCCHHHHHHHHHHH      ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "        cccchHHHHHHHHHHHHHHHHHHHHHHhCCCCChHHHHHHHHHHHHHHHHHHHHHccccchhHHHHHHHHHHHHHHHHHHHHHHhcCCCchhHHHHHHHHHHHHHHHHHHHHHHHcCCCHHHHHHHHHHHHHHHHHHchHHHHHHhcCHHHHHHHHHHHHHHHHHHHHHHHHHHHhccccCCCCcccCCHHHHHHHHHHHHHhhHHHHHHHHHHHHHH      ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  8,  22,  23,  38,  38,  65,  65,  98, 103, 126, 126, 184, 189,
                              226],
                             [  0,  14,  14,  29,  31,  58,  69, 102, 102, 125, 129, 187, 187,
                              224]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['G', 'G', 'F', 'G', 'S', 'Q', 'P', 'F', 'I', 'L', 'A', 'Y', 'I',
              'I', 'T', 'A', 'M', 'I', 'S', 'G', 'L', 'L', 'F', 'L', 'Y', 'L',
              'P', 'R', 'K', 'L', '-', '-', 'D', 'V', 'P', 'Q', 'K', 'F', 'G',
              'I', 'I', 'H', 'F', 'F', 'I', 'V', 'V', 'W', 'S', 'G', 'L', 'M',
              'Y', 'T', 'N', 'F', 'L', 'N', 'Q', '-', '-', '-', '-', '-', '-',
              '-', '-', '-', '-', '-', 'S', 'F', 'L', 'S', 'D', 'Y', 'A', 'W',
              'Y', 'M', 'D', 'W', 'M', 'V', 'S', 'T', 'P', 'L', 'I', 'L', 'L',
              'A', 'L', 'G', 'L', 'T', 'A', 'F', 'H', 'G', 'A', 'D', 'T', 'K',
              'R', 'Y', 'D', 'L', 'L', 'G', 'A', 'L', 'L', 'G', 'A', 'E', 'F',
              'T', 'L', 'V', 'I', 'T', 'G', 'L', 'L', 'A', 'Q', 'A', 'Q', 'G',
              'S', '-', '-', '-', '-', 'I', 'T', 'P', 'Y', 'Y', 'V', 'G', 'V',
              'L', 'L', 'L', 'L', 'G', 'V', 'V', 'Y', 'L', 'L', 'A', 'K', 'P',
              'F', 'R', 'E', 'I', 'A', 'E', 'E', 'S', 'S', 'D', 'G', 'L', 'A',
              'R', 'A', 'Y', 'K', 'I', 'L', 'A', 'G', 'Y', 'I', 'G', 'I', 'F',
              'F', 'L', 'S', 'Y', 'P', 'T', 'V', 'W', 'Y', 'I', 'S', 'G', 'I',
              'D', 'A', 'L', 'P', 'G', 'S', 'L', 'N', 'I', 'L', 'D', 'P', 'T',
              'Q', 'T', 'S', 'I', 'A', 'L', 'V', 'V', 'L', 'P', 'F', 'F', 'C',
              'K', 'Q', 'V', 'Y', 'G', 'F', 'L', 'D', 'M', 'Y', 'L', 'I', 'H',
              'K'],
             ['G', 'R', 'P', 'E', 'W', 'I', 'W', 'L', 'A', 'L', 'G', 'T', 'A',
              'L', '-', 'M', 'G', 'L', 'G', 'T', 'L', 'Y', 'F', 'L', 'V', 'K',
              'G', 'M', 'G', 'V', 'S', 'D', 'P', 'D', 'A', 'K', 'K', 'F', 'Y',
              'A', 'I', 'T', 'T', 'L', 'V', 'P', 'A', 'I', 'A', 'F', 'T', 'M',
              'Y', 'L', 'S', 'M', 'L', 'L', 'G', 'Y', 'G', 'L', 'T', 'M', 'V',
              'P', 'F', 'G', 'G', 'E', 'Q', 'N', 'P', 'I', 'Y', 'W', 'A', 'R',
              'Y', 'A', 'D', 'W', 'L', 'F', 'T', 'T', 'P', 'L', 'L', 'L', 'L',
              'D', 'L', 'A', 'L', 'L', 'V', 'D', 'A', 'D', 'Q', 'G', 'T', '-',
              '-', '-', '-', '-', 'I', 'L', 'A', 'L', 'V', 'G', 'A', 'D', 'G',
              'I', 'M', 'I', 'G', 'T', 'G', 'L', 'V', 'G', 'A', 'L', 'T', 'K',
              'V', 'Y', 'S', 'Y', 'R', 'F', 'V', 'W', 'W', 'A', 'I', 'S', 'T',
              'A', 'A', 'M', 'L', 'Y', 'I', 'L', 'Y', 'V', 'L', 'F', 'F', 'G',
              'F', 'T', 'S', 'K', 'A', 'E', 'S', 'M', 'R', 'P', 'E', 'V', 'A',
              'S', 'T', 'F', 'K', 'V', 'L', 'R', 'N', 'V', 'T', 'V', 'V', 'L',
              'W', 'S', 'A', 'Y', 'P', 'V', 'V', 'W', 'L', 'I', 'G', '-', '-',
              '-', '-', '-', 'S', 'E', 'G', 'A', 'G', 'I', 'V', 'P', 'L', 'N',
              'I', 'E', 'T', 'L', 'L', 'F', 'M', 'V', 'L', 'D', 'V', 'S', 'A',
              'K', 'V', 'G', 'F', 'G', 'L', 'I', 'L', 'L', 'R', 'S', 'R', 'A',
              'I']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 99.28)
        self.assertAlmostEqual(alignment.annotations["E-value"], 1.4e-15)
        self.assertAlmostEqual(alignment.annotations["Score"], 102.00)
        self.assertAlmostEqual(alignment.annotations["Identities"], 100)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.503)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 7.200)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                                                                                                                                                              CHHHHHHHHHHHHHHHHHHHHHHHHHHHccCCCCCCChhHHHHHHHHHHHHHHHHHHHHHHHhhhhhC",
        )
        self.assertEqual(alignment.query.id, "4Y9H:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[158:226],
            "RPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFG",
        )
        self.assertEqual(
            alignment[1],
            "RPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFG",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                                                                                              ~~~~~~~~~~l~~~~~i~W~~YPi~w~l~~~g~~~i~~~~~~~~~~ilDi~~K~~f~~~l~~~~~~~~",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            " ~~~~~~~~~~l~~~~~v~W~~YPi~w~lg~~g~~~i~~~~~~i~y~ilDv~~K~~fg~~~l~~~~~~~",
        )
        self.assertEqual(alignment.target.id, "1BCT_A")
        self.assertEqual(
            alignment.target.seq[1:69],
            "RPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFG",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "1BCT_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "BACTERIORHODOPSIN (FRAGMENT 163-231) (NMR, 14; PHOTORECEPTOR; NMR {Halobacterium salinarum} SCOP: j.35.1.1",
        )
        self.assertEqual(
            alignment[0],
            "RPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFG",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            " CCCHHHHHHHHHHHHHHHHHHHHHHHHHTTTTTTHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHC",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            " CHHHHHHHHHHHHHHHHHHHHHHHHHHHccCCCCCCCHHHHHHHHHHHHHHHHHHHHHHHHHHHHHhC",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  1,  69],
                             [158, 226]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['R', 'P', 'E', 'V', 'A', 'S', 'T', 'F', 'K', 'V', 'L', 'R', 'N',
              'V', 'T', 'V', 'V', 'L', 'W', 'S', 'A', 'Y', 'P', 'V', 'V', 'W',
              'L', 'I', 'G', 'S', 'E', 'G', 'A', 'G', 'I', 'V', 'P', 'L', 'N',
              'I', 'E', 'T', 'L', 'L', 'F', 'M', 'V', 'L', 'D', 'V', 'S', 'A',
              'K', 'V', 'G', 'F', 'G', 'L', 'I', 'L', 'L', 'R', 'S', 'R', 'A',
              'I', 'F', 'G'],
             ['R', 'P', 'E', 'V', 'A', 'S', 'T', 'F', 'K', 'V', 'L', 'R', 'N',
              'V', 'T', 'V', 'V', 'L', 'W', 'S', 'A', 'Y', 'P', 'V', 'V', 'W',
              'L', 'I', 'G', 'S', 'E', 'G', 'A', 'G', 'I', 'V', 'P', 'L', 'N',
              'I', 'E', 'T', 'L', 'L', 'F', 'M', 'V', 'L', 'D', 'V', 'S', 'A',
              'K', 'V', 'G', 'F', 'G', 'L', 'I', 'L', 'L', 'R', 'S', 'R', 'A',
              'I', 'F', 'G']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 98.44)
        self.assertAlmostEqual(alignment.annotations["E-value"], 2e-10)
        self.assertAlmostEqual(alignment.annotations["Score"], 77.33)
        self.assertAlmostEqual(alignment.annotations["Identities"], 97)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.598)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 7.100)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "CCcHHHHHHHHHHHHHHHHHHHHHHhcCCCChhhHHHHHHHHHHHHHHHHHHHHHHhCCCceeecc                                                                                                                                                                ",
        )
        self.assertEqual(alignment.query.id, "4Y9H:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[0:66],
            "GRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPF",
        )
        self.assertEqual(
            alignment[1],
            "GRPEWIWLALGTALMGLGTLYFLVKGMG-VSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPF",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "~~~~~~~~~~~~~~m~~~~~~f~~~~~~~~~~~~r~~~~l~~~i~~ia~~~Y~~ma~~~g~~~~~~                                                                                                                                                                ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "     ~~~~~~~lw~~~~~m~~~~~~f~~~s~~~~~~~~R~~~~~~~~i~~iaaiaY~~MA~g~G~~~v~~",
        )
        self.assertEqual(alignment.target.id, "1BHA_A")
        self.assertEqual(
            alignment.target.seq[5:71],
            "GRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPF",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "1BHA_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "BACTERIORHODOPSIN (PROTEOLYTIC FRAGMENT 1 -; PHOTORECEPTOR; NMR {Halobacterium salinarum} SCOP: j.35.1.1; Related PDB entries: 1BHB_A",
        )
        self.assertEqual(
            alignment[0],
            "GRPEWIWLALGTALMGLGTLYFLVKGMGVSD-PDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPF",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "     SCHHHHHHHHHHHHHHHHHHHHHHHHHHCCHHHHHHHHHHSHHHHHHHHHHHHHTCSSSSCC    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "     ChhHHHHHHHHHHHHHHHHHHHHHHHhCCCChhhHHHHHHHHHHHHHHHHHHHHHHhcCCcccCCC",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 5, 33, 34, 36, 36, 71],
                             [ 0, 28, 28, 30, 31, 66]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['G', 'R', 'P', 'E', 'W', 'I', 'W', 'L', 'A', 'L', 'G', 'T', 'A',
              'L', 'M', 'G', 'L', 'G', 'T', 'L', 'Y', 'F', 'L', 'V', 'K', 'G',
              'M', 'G', 'V', 'S', 'D', '-', 'P', 'D', 'A', 'K', 'K', 'F', 'Y',
              'A', 'I', 'T', 'T', 'L', 'V', 'P', 'A', 'I', 'A', 'F', 'T', 'M',
              'Y', 'L', 'S', 'M', 'L', 'L', 'G', 'Y', 'G', 'L', 'T', 'M', 'V',
              'P', 'F'],
             ['G', 'R', 'P', 'E', 'W', 'I', 'W', 'L', 'A', 'L', 'G', 'T', 'A',
              'L', 'M', 'G', 'L', 'G', 'T', 'L', 'Y', 'F', 'L', 'V', 'K', 'G',
              'M', 'G', '-', 'V', 'S', 'D', 'P', 'D', 'A', 'K', 'K', 'F', 'Y',
              'A', 'I', 'T', 'T', 'L', 'V', 'P', 'A', 'I', 'A', 'F', 'T', 'M',
              'Y', 'L', 'S', 'M', 'L', 'L', 'G', 'Y', 'G', 'L', 'T', 'M', 'V',
              'P', 'F']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 96.55)
        self.assertAlmostEqual(alignment.annotations["E-value"], 3.3e-05)
        self.assertAlmostEqual(alignment.annotations["Score"], 51.24)
        self.assertAlmostEqual(alignment.annotations["Identities"], 31)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.326)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 4.800)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "       HHHHHHHHHHHHHHHHHHhcCCCChhhHHHHHHHHHHHHHHHHHHHHHHhCC                                                                                                                                                                       ",
        )
        self.assertEqual(alignment.query.id, "4Y9H:A|PDBID|CHAIN|SEQUENCE")
        self.assertEqual(
            alignment.query.seq[7:59],
            "LALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGY",
        )
        self.assertEqual(
            alignment[1], "LALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGY"
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "       ~~~~~~~m~~~~~~f~~~~~~~~~~~~r~~~~l~~~i~~ia~~~Y~~ma~~~                                                                                                                                                                       ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "              fs~g~aaM~AatvfF~l~~~~V~~~yr~s~~lsalVt~iA~~hY~~m~~~W    ",
        )
        self.assertEqual(alignment.target.id, "5ABB_Z")
        self.assertEqual(
            alignment.target.seq[14:65],
            "FWLVTAALLASTVFFFVERDRVSAKWKTSLTVSGLVTGIAFWHYMYMRGVW",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "5ABB_Z")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "PROTEIN TRANSLOCASE SUBUNIT SECY, PROTEIN; TRANSLATION, RIBOSOME, MEMBRANE PROTEIN, TRANSLOCON; 8.0A {ESCHERICHIA COLI}",
        )
        self.assertEqual(
            alignment[0], "FWLVTAALLASTVFFFVERDRVS-AKWKTSLTVSGLVTGIAFWHYMYMRGVW"
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "              HHHHHHHHHHHHHHHHTTTSSCCSCCCCSSSCHHHHHHHHHHHHHHHHHHH    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "              HHHHHHHHHHHHHHHHHHhhccChhhhHHHHHHHHHHHHHHHHHHHHHHHH    ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[14, 37, 37, 65],
                             [ 7, 30, 31, 59]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['F', 'W', 'L', 'V', 'T', 'A', 'A', 'L', 'L', 'A', 'S', 'T', 'V',
              'F', 'F', 'F', 'V', 'E', 'R', 'D', 'R', 'V', 'S', '-', 'A', 'K',
              'W', 'K', 'T', 'S', 'L', 'T', 'V', 'S', 'G', 'L', 'V', 'T', 'G',
              'I', 'A', 'F', 'W', 'H', 'Y', 'M', 'Y', 'M', 'R', 'G', 'V', 'W'],
             ['L', 'A', 'L', 'G', 'T', 'A', 'L', 'M', 'G', 'L', 'G', 'T', 'L',
              'Y', 'F', 'L', 'V', 'K', 'G', 'M', 'G', 'V', 'S', 'D', 'P', 'D',
              'A', 'K', 'K', 'F', 'Y', 'A', 'I', 'T', 'T', 'L', 'V', 'P', 'A',
              'I', 'A', 'F', 'T', 'M', 'Y', 'L', 'S', 'M', 'L', 'L', 'G', 'Y']],
            dtype='U')
                # fmt: on
            )
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_length(self):
        """Test getting the number of alignments without parsing the file."""
        stream = open(self.path)
        alignments = Align.parse(stream, "hhr")
        stream.close()
        self.assertEqual(len(alignments), 29)


class Align_hhr_hhpred_9590198(unittest.TestCase):
    path = os.path.join("HHsuite", "hhpred_9590198.hhr")

    def test_reading(self):
        alignments = Align.parse(self.path, "hhr")
        self.assertEqual(alignments.metadata["No_of_seqs"], (157, 584))
        self.assertAlmostEqual(alignments.metadata["Neff"], 6.82639)
        self.assertEqual(alignments.metadata["Searched_HMMs"], 64707)
        self.assertEqual(alignments.metadata["Rundate"], "Fri Feb  1 15:48:30 2019")
        self.assertEqual(
            alignments.metadata["Command line"],
            "hhsearch -cpu 8 -i ../results/full.a3m -d /cluster/toolkit/production/databases/hh-suite/mmcif70/pdb70 -d /cluster/toolkit/production/databases/hh-suite/pfama/pfama -o ../results/9590198.hhr -oa3m ../results/9590198.a3m -p 20 -Z 250 -loc -z 1 -b 1 -B 250 -ssm 2 -sc 1 -seq 1 -dbstrlen 10000 -norealign -maxres 32000 -contxt /cluster/toolkit/production/bioprogs/tools/hh-suite-build/data/context_data.crf",
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 100.00)
        self.assertAlmostEqual(alignment.annotations["E-value"], 9.9e-102)
        self.assertAlmostEqual(alignment.annotations["Score"], 792.76)
        self.assertAlmostEqual(alignment.annotations["Identities"], 53)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.957)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 7.100)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                     CCcHHHHHHHHHHccCCcCcEEEEECCCCCcCCCEEEEeCCCCeEEEecCccceeeEEEEecCCCCEEEECCEEecCCCCCCcHHHHHHhhcCCCCCCccCcccEEEEEeCcEEEEeeCcccCCCCccCCCcccccccCCCCCCCeEEEEEEEeCCChHhcCCCCCCcccccCcceeeEEEEeeCCCCCceEEEEEEecCCCCccccccceeEEEEEEEecCCHHHHHHHhcCCcEEEEeccccccccCCCCCcCCCCcCCccchhhHhhccceeecCCCCEEEEEEecCCCCCCccccccccCcEEeeccccccCCCCCcceeccCcCHHHHHHHhCCCCCCCeEeeCCCCCCCCCCCCceEEeccCCEEEEEecCCCEEEEEEe               ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(
            alignment.query.seq[21:407],
            "GMPLAQAVAILQKHCRIIKNVQVLYSEQSPLSHDLILNLTQDGIKLMFDAFNQRLKVIEVCDLTKVKLKYCGVHFNSQAIAPTIEQIDQSFGATHPGVYNSAEQLFHLNFRGLSFSFQLDSWTEAPKYEPNFAHGLASLQIPHGATVKRMYIYSGNSLQDTKAPMMPLSCFLGNVYAESVDVLRDGTGPAGLRLRLLAAGCGPGLLADAKMRVFERSVYFGDSCQDVLSMLGSPHKVFYKSEDKMKIHSPSPHKQVPSKCNDYFFNYFTLGVDILFDANTHKVKKFVLHTNYPGHYNFNIYHRCEFKIPLAIKKENADGQTETCTTYSKWDNIQELLGHPVEKPVVLHRSSSPNNTNPFGSTFCFGLQRMIFEVMQNNHIASVTLY",
        )
        self.assertEqual(
            alignment[1],
            "GMPLAQAVAILQKHCRIIKNVQVLYSEQSPLSHDLILNLTQDGIKLMFDAFNQRLKVIEVCDLTKVKLKYCGVHFNSQAIAPTIEQIDQSFGATHPGVYNSAEQLFHLNFRGLSFSFQLDSWTEAPKYEPNFAHGLASLQIPHGA--TVKRMYIYSGNSLQ---------DTKA-PMMPLSCFLGNVYAESVDVLRDGTGPAGLRLRLLAAGCGPGLLADAKMRVFERSVYFGDSCQDVLSMLGSPHKVFYKSEDKMKIHSPSPHKQVPSKCNDYFFNYFTLGVDILFDANTHKVKKFVLHTNYPGHYNFNIYHRCEFKIPLAIKKENADG------QTETCTTYSKWDNIQELLGHPVEKPVVLHRSSSPNNTNPFGSTFCFGLQRMIFEVMQNNHIASVTLY",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                     G~sL~~vi~~Lr~~~~~~~~vel~Ys~~~pl~~~Ivi~L~~~GirL~Fd~~~QrL~lIEv~d~~~~~L~Y~~~~~~~~~~~ptf~~I~~~FGPTyPG~yd~~~~~Y~LsYpGisF~Fpi~~~~~~~~~~~~~~~~l~~l~~~~~~~~s~i~If~G~s~~~~~~~~~p~~~~~~~~~~~~v~v~~~~~~~~gl~~~~~~~~~~~~~~~~~~~~~~~~~I~~G~T~QDVl~~LG~P~~~f~K~ddrm~IH~~~~~~~~~~~~~dyF~NYF~lG~DiLfd~~t~~v~KiILHtN~PG~~~F~~Y~RC~w~i~~~~~~~~~~~~~~~it~~~~~~~i~~~l~~~~~~pvvlnR~~s~~~~~~fg~T~~yG~~g~IfEVm~ng~IasvTl~               ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "G~~L~~vi~~L~~~~~~~~~Vel~Ys~~~p~~~~Ivl~L~~~glrL~Fd~~~QrL~lIEv~d~~~~~L~Y~~~~~~~~~~~~t~~~Iy~~FGPTyPG~~d~~~~~y~LsYpGisF~F~~~~~~~~~~~~~~~~l~~~~~~~~~~~~m~If~G~~~~~~~~~~~~~~~~~~p~~p~~~~~~~~~~~~v~v~~~~~~~~gl~l~f~~~~~~~~~~~~~~~~~~~I~~G~t~QDVl~~LG~P~~~f~K~ddrm~Ih~~~~~~~~~~~~~yF~NYF~lG~DiLfd~~th~v~KiILHtN~Pg~~~F~~Y~RC~~~l~~~~~~~~~~~~~~~~~~~~i~~~~~~~~i~~~l~~~~~p~vlnR~~~~~~~~~~g~T~lyg~~g~ifEV~~ng~Iasvtlf",
        )
        self.assertEqual(alignment.target.id, "H9J4A9_BOMMO/5")
        self.assertEqual(
            alignment.target.seq[0:394],
            "GMHFSQSVAIIQSQVGTIRGVQVLYSDQNPLSVDLVINMPQDGMRLIFDPVAQRLKIIEIYNMKLVKLRYSGMCFNSPEITPSIEQVEHCFGATHPGLYDSQRHLFALNFRGLSFYFPVDSKFEPGYAHGLGSLQFPNGGSPVVSRTTIYYGSQHQLSSNTSSRVSGVPLPDLPLSCYRQQLHLRRCDVLRNTTSTMGLRLHMFTEGTSRALEPSQVALVRVVRFGDSCQGVARALGAPARLYYKADDKMRIHRPTARRRPPPASDYLFNYFTLGLDVLFDARTNQVKKFVLHTNYPGHYNFNMYHRCEFELTVQPDKSEAHSLVESGGGVAVTAYSKWEVVSRALRVCERPVVLNRASSTNTTNPFGSTFCYGYQDIIFEVMSNNYIASITLY",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "PF03676.14")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "; UPF0183 ; Uncharacterised protein family (UPF0183)",
        )
        self.assertEqual(
            alignment[0],
            "GMHFSQSVAIIQSQVGTIRGVQVLYSDQNPLSVDLVINMPQDGMRLIFDPVAQRLKIIEIYNMKLVKLRYSGMCFNSPEITPSIEQVEHCFGATHPGLYDSQRHLFALNFRGLSFYFPVDS-----KFEPGYAHGLGSLQFPNGGSPVVSRTTIYYGSQHQLSSNTSSRVSGVPLPDLPLSCYRQQLHLRRCDVLRNTTSTMGLRLHMFTEGT--SRALEPSQVALVRVVRFGDSCQGVARALGAPARLYYKADDKMRIHRPTARRR-PPPASDYLFNYFTLGLDVLFDARTNQVKKFVLHTNYPGHYNFNMYHRCEFELTVQPD-KSEAHSLVESGGGVAVTAYSKWEVVSRAL-RVCERPVVLNRASSTNTTNPFGSTFCYGYQDIIFEVMSNNYIASITLY",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "CCcHHHHHHHHHhcCCccceEEEEECCCCCCCCcEEEEcCCCceEEEEcCccCeEEEEEEEcCCCcEEEEcCeeecCCCCCCcHHHHHHHhcCCCCCcccCCCCeEEEEeCcEEEEeeCCcccCCCCccccccccCCCCCCCceEEEEEEeCCccccccccccccCCCCCCCCCchhcCCCCceEEEEEEEcCCCCceEEEEEEcCCCCCcCCccceeceEEEEecCCHHHHHHHhCCCceEEEcCcccccccCCCCCCCCCCcchhchhHHHhCcceeecCCCcEEEEEEeeCCCCCCccccccccCeEEEEeCCCcchhcceecCCCCceecCcccHHHHHHHhccCCCCEEEecCCCCCCCCCCCcEEEEecCCEEEEEecCCcEEEEEeC",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  0, 121, 121, 140, 142, 156, 165, 169, 170, 208, 208, 260, 260,
                              317, 317, 322, 328, 346, 346, 394],
                             [ 21, 142, 147, 166, 166, 180, 180, 184, 184, 222, 224, 276, 277,
                              334, 335, 340, 340, 358, 359, 407]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['G', 'M', 'H', 'F', 'S', 'Q', 'S', 'V', 'A', 'I', 'I', 'Q', 'S',
              'Q', 'V', 'G', 'T', 'I', 'R', 'G', 'V', 'Q', 'V', 'L', 'Y', 'S',
              'D', 'Q', 'N', 'P', 'L', 'S', 'V', 'D', 'L', 'V', 'I', 'N', 'M',
              'P', 'Q', 'D', 'G', 'M', 'R', 'L', 'I', 'F', 'D', 'P', 'V', 'A',
              'Q', 'R', 'L', 'K', 'I', 'I', 'E', 'I', 'Y', 'N', 'M', 'K', 'L',
              'V', 'K', 'L', 'R', 'Y', 'S', 'G', 'M', 'C', 'F', 'N', 'S', 'P',
              'E', 'I', 'T', 'P', 'S', 'I', 'E', 'Q', 'V', 'E', 'H', 'C', 'F',
              'G', 'A', 'T', 'H', 'P', 'G', 'L', 'Y', 'D', 'S', 'Q', 'R', 'H',
              'L', 'F', 'A', 'L', 'N', 'F', 'R', 'G', 'L', 'S', 'F', 'Y', 'F',
              'P', 'V', 'D', 'S', '-', '-', '-', '-', '-', 'K', 'F', 'E', 'P',
              'G', 'Y', 'A', 'H', 'G', 'L', 'G', 'S', 'L', 'Q', 'F', 'P', 'N',
              'G', 'G', 'S', 'P', 'V', 'V', 'S', 'R', 'T', 'T', 'I', 'Y', 'Y',
              'G', 'S', 'Q', 'H', 'Q', 'L', 'S', 'S', 'N', 'T', 'S', 'S', 'R',
              'V', 'S', 'G', 'V', 'P', 'L', 'P', 'D', 'L', 'P', 'L', 'S', 'C',
              'Y', 'R', 'Q', 'Q', 'L', 'H', 'L', 'R', 'R', 'C', 'D', 'V', 'L',
              'R', 'N', 'T', 'T', 'S', 'T', 'M', 'G', 'L', 'R', 'L', 'H', 'M',
              'F', 'T', 'E', 'G', 'T', '-', '-', 'S', 'R', 'A', 'L', 'E', 'P',
              'S', 'Q', 'V', 'A', 'L', 'V', 'R', 'V', 'V', 'R', 'F', 'G', 'D',
              'S', 'C', 'Q', 'G', 'V', 'A', 'R', 'A', 'L', 'G', 'A', 'P', 'A',
              'R', 'L', 'Y', 'Y', 'K', 'A', 'D', 'D', 'K', 'M', 'R', 'I', 'H',
              'R', 'P', 'T', 'A', 'R', 'R', 'R', '-', 'P', 'P', 'P', 'A', 'S',
              'D', 'Y', 'L', 'F', 'N', 'Y', 'F', 'T', 'L', 'G', 'L', 'D', 'V',
              'L', 'F', 'D', 'A', 'R', 'T', 'N', 'Q', 'V', 'K', 'K', 'F', 'V',
              'L', 'H', 'T', 'N', 'Y', 'P', 'G', 'H', 'Y', 'N', 'F', 'N', 'M',
              'Y', 'H', 'R', 'C', 'E', 'F', 'E', 'L', 'T', 'V', 'Q', 'P', 'D',
              '-', 'K', 'S', 'E', 'A', 'H', 'S', 'L', 'V', 'E', 'S', 'G', 'G',
              'G', 'V', 'A', 'V', 'T', 'A', 'Y', 'S', 'K', 'W', 'E', 'V', 'V',
              'S', 'R', 'A', 'L', '-', 'R', 'V', 'C', 'E', 'R', 'P', 'V', 'V',
              'L', 'N', 'R', 'A', 'S', 'S', 'T', 'N', 'T', 'T', 'N', 'P', 'F',
              'G', 'S', 'T', 'F', 'C', 'Y', 'G', 'Y', 'Q', 'D', 'I', 'I', 'F',
              'E', 'V', 'M', 'S', 'N', 'N', 'Y', 'I', 'A', 'S', 'I', 'T', 'L',
              'Y'],
             ['G', 'M', 'P', 'L', 'A', 'Q', 'A', 'V', 'A', 'I', 'L', 'Q', 'K',
              'H', 'C', 'R', 'I', 'I', 'K', 'N', 'V', 'Q', 'V', 'L', 'Y', 'S',
              'E', 'Q', 'S', 'P', 'L', 'S', 'H', 'D', 'L', 'I', 'L', 'N', 'L',
              'T', 'Q', 'D', 'G', 'I', 'K', 'L', 'M', 'F', 'D', 'A', 'F', 'N',
              'Q', 'R', 'L', 'K', 'V', 'I', 'E', 'V', 'C', 'D', 'L', 'T', 'K',
              'V', 'K', 'L', 'K', 'Y', 'C', 'G', 'V', 'H', 'F', 'N', 'S', 'Q',
              'A', 'I', 'A', 'P', 'T', 'I', 'E', 'Q', 'I', 'D', 'Q', 'S', 'F',
              'G', 'A', 'T', 'H', 'P', 'G', 'V', 'Y', 'N', 'S', 'A', 'E', 'Q',
              'L', 'F', 'H', 'L', 'N', 'F', 'R', 'G', 'L', 'S', 'F', 'S', 'F',
              'Q', 'L', 'D', 'S', 'W', 'T', 'E', 'A', 'P', 'K', 'Y', 'E', 'P',
              'N', 'F', 'A', 'H', 'G', 'L', 'A', 'S', 'L', 'Q', 'I', 'P', 'H',
              'G', 'A', '-', '-', 'T', 'V', 'K', 'R', 'M', 'Y', 'I', 'Y', 'S',
              'G', 'N', 'S', 'L', 'Q', '-', '-', '-', '-', '-', '-', '-', '-',
              '-', 'D', 'T', 'K', 'A', '-', 'P', 'M', 'M', 'P', 'L', 'S', 'C',
              'F', 'L', 'G', 'N', 'V', 'Y', 'A', 'E', 'S', 'V', 'D', 'V', 'L',
              'R', 'D', 'G', 'T', 'G', 'P', 'A', 'G', 'L', 'R', 'L', 'R', 'L',
              'L', 'A', 'A', 'G', 'C', 'G', 'P', 'G', 'L', 'L', 'A', 'D', 'A',
              'K', 'M', 'R', 'V', 'F', 'E', 'R', 'S', 'V', 'Y', 'F', 'G', 'D',
              'S', 'C', 'Q', 'D', 'V', 'L', 'S', 'M', 'L', 'G', 'S', 'P', 'H',
              'K', 'V', 'F', 'Y', 'K', 'S', 'E', 'D', 'K', 'M', 'K', 'I', 'H',
              'S', 'P', 'S', 'P', 'H', 'K', 'Q', 'V', 'P', 'S', 'K', 'C', 'N',
              'D', 'Y', 'F', 'F', 'N', 'Y', 'F', 'T', 'L', 'G', 'V', 'D', 'I',
              'L', 'F', 'D', 'A', 'N', 'T', 'H', 'K', 'V', 'K', 'K', 'F', 'V',
              'L', 'H', 'T', 'N', 'Y', 'P', 'G', 'H', 'Y', 'N', 'F', 'N', 'I',
              'Y', 'H', 'R', 'C', 'E', 'F', 'K', 'I', 'P', 'L', 'A', 'I', 'K',
              'K', 'E', 'N', 'A', 'D', 'G', '-', '-', '-', '-', '-', '-', 'Q',
              'T', 'E', 'T', 'C', 'T', 'T', 'Y', 'S', 'K', 'W', 'D', 'N', 'I',
              'Q', 'E', 'L', 'L', 'G', 'H', 'P', 'V', 'E', 'K', 'P', 'V', 'V',
              'L', 'H', 'R', 'S', 'S', 'S', 'P', 'N', 'N', 'T', 'N', 'P', 'F',
              'G', 'S', 'T', 'F', 'C', 'F', 'G', 'L', 'Q', 'R', 'M', 'I', 'F',
              'E', 'V', 'M', 'Q', 'N', 'N', 'H', 'I', 'A', 'S', 'V', 'T', 'L',
              'Y']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 71.77)
        self.assertAlmostEqual(alignment.annotations["E-value"], 67)
        self.assertAlmostEqual(alignment.annotations["Score"], 28.65)
        self.assertAlmostEqual(alignment.annotations["Identities"], 14)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.172)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 9.600)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                                                                                                                                                                                                                                              EEecCCHHHHHHHhcCCcEEEEeccccccccCCCCCcCCCCcCCccchhhHhhccceeecCCCCEEEEEEecCCCCCCccccccccCcEEeeccccccCCCCCcceeccCcCHHHHHHHhCCCCCCCeEeeCCCCCCCCCCCCceEEeccCCEEEEEecCCCEEEEEEeCC             ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(
            alignment.query.seq[238:409],
            "VYFGDSCQDVLSMLGSPHKVFYKSEDKMKIHSPSPHKQVPSKCNDYFFNYFTLGVDILFDANTHKVKKFVLHTNYPGHYNFNIYHRCEFKIPLAIKKENADGQTETCTTYSKWDNIQELLGHPVEKPVVLHRSSSPNNTNPFGSTFCFGLQRMIFEVMQNNHIASVTLYGP",
        )
        self.assertEqual(
            alignment[1],
            "VYF---------GDSCQDVLSMLGSPHKVFYKSEDKMKIHSPSPHKQVPSKCNDYFFNYFTLGVDILFDANTHKVKKFVL-----HTNYPGHYNFNIYHRCEFKIPLAIKKENADGQTETCTTYSKWDNIQELLGHPVEKPVVLHRSSSPNNTNPFGSTFCFGLQ------RMIFEVM-QNNHIASVTLYGP",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                                                                                                                                                                              I~~G~T~QDVl~~LG~P~~~f~K~ddrm~IH~~~~~~~~~~~~~dyF~NYF~lG~DiLfd~~t~~v~KiILHtN~PG~~~F~~Y~RC~w~i~~~~~~~~~~~~~~~it~~~~~~~i~~~l~~~~~~pvvlnR~~s~~~~~~fg~T~~yG~~g~IfEVm~ng~IasvTl~~~             ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "           i~~~~~~~~~~~G~t~~eV~~~lG~p~~~~~~~~~~~~~~~~~w~~~~~~~~v~~~~~~~~~~~~~~~~~~~~~~~t~~~~~~i~~G~t~~eV~~~lG~p~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~F~~g~l~~~~~~g~",
        )
        self.assertEqual(alignment.target.id, "Q8DTX1_STRMU/4")
        self.assertEqual(
            alignment.target.seq[11:159],
            "IKVTTDQNHFSGGTSIEQLKQWFGDPNKSEQRNAGNITLDSYTWVKDGAVINAQLYKNSTVARSISNFSFSREAKIGKEDYDELKIGESYKKVVEKLGEPDVLSQSMSSDKEEMQTVWSSGIKTKSSSATIELYFENGLLKNKTQKDL",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "PF12978.7")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "; DUF3862 ; Domain of Unknown Function with PDB structure (DUF3862)",
        )
        self.assertEqual(
            alignment[0],
            "IKVTTDQNHFSGGTSIEQLKQWFGDPNKSEQRNAG-------------NITLDSYTWVKDGAVINAQLY-KNSTVARSISNFSFSREAKIGKEDYD-----------------------ELKIGESYKKVVEKL----GEPDVLSQSMS---SDKEEMQTVWSSGIKTKSSSATIELYFENGLLKNKTQKDL",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "           eeecCCcccCCCCCCHHHHHHHHCCCceEEeeeeCCeEEEEEEEEeCCeEEEEEEECCEEEEEEEEeceecCCCCCCHHHHHhcCCCCCHHHHHHHHCCCCeeEEEeeCCceEEEEEEEecccCCCCCcEEEEEEECCeEEEeEecCC",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 11,  14,  23,  46,  46,  67,  67,  77,  82,  93,  93, 108, 108,
                              119, 119, 132, 138, 145, 146, 159],
                             [238, 241, 241, 264, 277, 298, 299, 309, 309, 320, 343, 358, 362,
                              373, 376, 389, 389, 396, 396, 409]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['I', 'K', 'V', 'T', 'T', 'D', 'Q', 'N', 'H', 'F', 'S', 'G', 'G',
              'T', 'S', 'I', 'E', 'Q', 'L', 'K', 'Q', 'W', 'F', 'G', 'D', 'P',
              'N', 'K', 'S', 'E', 'Q', 'R', 'N', 'A', 'G', '-', '-', '-', '-',
              '-', '-', '-', '-', '-', '-', '-', '-', '-', 'N', 'I', 'T', 'L',
              'D', 'S', 'Y', 'T', 'W', 'V', 'K', 'D', 'G', 'A', 'V', 'I', 'N',
              'A', 'Q', 'L', 'Y', '-', 'K', 'N', 'S', 'T', 'V', 'A', 'R', 'S',
              'I', 'S', 'N', 'F', 'S', 'F', 'S', 'R', 'E', 'A', 'K', 'I', 'G',
              'K', 'E', 'D', 'Y', 'D', '-', '-', '-', '-', '-', '-', '-', '-',
              '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',
              '-', '-', 'E', 'L', 'K', 'I', 'G', 'E', 'S', 'Y', 'K', 'K', 'V',
              'V', 'E', 'K', 'L', '-', '-', '-', '-', 'G', 'E', 'P', 'D', 'V',
              'L', 'S', 'Q', 'S', 'M', 'S', '-', '-', '-', 'S', 'D', 'K', 'E',
              'E', 'M', 'Q', 'T', 'V', 'W', 'S', 'S', 'G', 'I', 'K', 'T', 'K',
              'S', 'S', 'S', 'A', 'T', 'I', 'E', 'L', 'Y', 'F', 'E', 'N', 'G',
              'L', 'L', 'K', 'N', 'K', 'T', 'Q', 'K', 'D', 'L'],
             ['V', 'Y', 'F', '-', '-', '-', '-', '-', '-', '-', '-', '-', 'G',
              'D', 'S', 'C', 'Q', 'D', 'V', 'L', 'S', 'M', 'L', 'G', 'S', 'P',
              'H', 'K', 'V', 'F', 'Y', 'K', 'S', 'E', 'D', 'K', 'M', 'K', 'I',
              'H', 'S', 'P', 'S', 'P', 'H', 'K', 'Q', 'V', 'P', 'S', 'K', 'C',
              'N', 'D', 'Y', 'F', 'F', 'N', 'Y', 'F', 'T', 'L', 'G', 'V', 'D',
              'I', 'L', 'F', 'D', 'A', 'N', 'T', 'H', 'K', 'V', 'K', 'K', 'F',
              'V', 'L', '-', '-', '-', '-', '-', 'H', 'T', 'N', 'Y', 'P', 'G',
              'H', 'Y', 'N', 'F', 'N', 'I', 'Y', 'H', 'R', 'C', 'E', 'F', 'K',
              'I', 'P', 'L', 'A', 'I', 'K', 'K', 'E', 'N', 'A', 'D', 'G', 'Q',
              'T', 'E', 'T', 'C', 'T', 'T', 'Y', 'S', 'K', 'W', 'D', 'N', 'I',
              'Q', 'E', 'L', 'L', 'G', 'H', 'P', 'V', 'E', 'K', 'P', 'V', 'V',
              'L', 'H', 'R', 'S', 'S', 'S', 'P', 'N', 'N', 'T', 'N', 'P', 'F',
              'G', 'S', 'T', 'F', 'C', 'F', 'G', 'L', 'Q', '-', '-', '-', '-',
              '-', '-', 'R', 'M', 'I', 'F', 'E', 'V', 'M', '-', 'Q', 'N', 'N',
              'H', 'I', 'A', 'S', 'V', 'T', 'L', 'Y', 'G', 'P']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 70.25)
        self.assertAlmostEqual(alignment.annotations["E-value"], 4.9)
        self.assertAlmostEqual(alignment.annotations["Score"], 30.39)
        self.assertAlmostEqual(alignment.annotations["Identities"], 30)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.363)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 10.200)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                                                                                                                                                                                                                                              EEecCCHHHHHHHhcCCcEE                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(alignment.query.seq[238:258], "VYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(alignment[1], "VYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                                                                                                                                                                              I~~G~T~QDVl~~LG~P~~~                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "           i~~Gmt~~eV~~~lG~P~~~                                          ",
        )
        self.assertEqual(alignment.target.id, "A6W2D6_MARMS/3")
        self.assertEqual(alignment.target.seq[11:31], "LQIGMSESQVTYLLGNPMLR")
        self.assertEqual(alignment.target.annotations["hmm_name"], "PF04355.13")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "; SmpA_OmlA ; SmpA / OmlA family",
        )
        self.assertEqual(alignment[0], "LQIGMSESQVTYLLGNPMLR")
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "           CCCCCCHHHHHHHHCCCcee                                          ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 11,  31],
                             [238, 258]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['L', 'Q', 'I', 'G', 'M', 'S', 'E', 'S', 'Q', 'V', 'T', 'Y', 'L',
              'L', 'G', 'N', 'P', 'M', 'L', 'R'],
             ['V', 'Y', 'F', 'G', 'D', 'S', 'C', 'Q', 'D', 'V', 'L', 'S', 'M',
              'L', 'G', 'S', 'P', 'H', 'K', 'V']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 68.06)
        self.assertAlmostEqual(alignment.annotations["E-value"], 5.9)
        self.assertAlmostEqual(alignment.annotations["Score"], 32.02)
        self.assertAlmostEqual(alignment.annotations["Identities"], 30)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.380)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 9.600)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                                                                                                                                                                                                                                              EEecCCHHHHHHHhcCCcEE                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(alignment.query.seq[238:258], "VYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(alignment[1], "VYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                                                                                                                                                                              I~~G~T~QDVl~~LG~P~~~                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                          l~~G~t~~eV~~~lG~P~~~                                                ",
        )
        self.assertEqual(alignment.target.id, "5EKQ_E")
        self.assertEqual(alignment.target.seq[26:46], "IRVGMTQQQVAYALGTPLMS")
        self.assertEqual(alignment.target.annotations["hmm_name"], "5EKQ_E")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Outer membrane protein assembly factor; membrane protein, insertase, beta-barrel, outer; 3.392A {Escherichia coli (strain K12)}; Related PDB entries: 2KM7_A 2KXX_A 5AYW_E 5LJO_E",
        )
        self.assertEqual(alignment[0], "IRVGMTQQQVAYALGTPLMS")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "                          CCTTCCSTTTTTSSCCCSEE                                                ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                          CCCCCCHHHHHHHhCCCcee                                                ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 26,  46],
                             [238, 258]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['I', 'R', 'V', 'G', 'M', 'T', 'Q', 'Q', 'Q', 'V', 'A', 'Y', 'A',
              'L', 'G', 'T', 'P', 'L', 'M', 'S'],
             ['V', 'Y', 'F', 'G', 'D', 'S', 'C', 'Q', 'D', 'V', 'L', 'S', 'M',
              'L', 'G', 'S', 'P', 'H', 'K', 'V']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 58.72)
        self.assertAlmostEqual(alignment.annotations["E-value"], 13)
        self.assertAlmostEqual(alignment.annotations["Score"], 30.74)
        self.assertAlmostEqual(alignment.annotations["Identities"], 15)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.388)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 9.500)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                                                                                                                                                                                                                                              EEecCCHHHHHHHhcCCcEE                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(alignment.query.seq[238:258], "VYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(alignment[1], "VYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                                                                                                                                                                              I~~G~T~QDVl~~LG~P~~~                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                                 i~~Gmtk~eV~~~lG~P~~~                                                   ",
        )
        self.assertEqual(alignment.target.id, "A3QGB8_SHELP/1")
        self.assertEqual(alignment.target.seq[33:53], "LSLGMTRDQVMTLMGTADFN")
        self.assertEqual(alignment.target.annotations["hmm_name"], "PF11399.8")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "; DUF3192 ; Protein of unknown function (DUF3192)",
        )
        self.assertEqual(alignment[0], "LSLGMTRDQVMTLMGTADFN")
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                                 cCCCCCHHHHHHHhCCCcee                                                   ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 33,  53],
                             [238, 258]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['L', 'S', 'L', 'G', 'M', 'T', 'R', 'D', 'Q', 'V', 'M', 'T', 'L',
              'M', 'G', 'T', 'A', 'D', 'F', 'N'],
             ['V', 'Y', 'F', 'G', 'D', 'S', 'C', 'Q', 'D', 'V', 'L', 'S', 'M',
              'L', 'G', 'S', 'P', 'H', 'K', 'V']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 58.46)
        self.assertAlmostEqual(alignment.annotations["E-value"], 11)
        self.assertAlmostEqual(alignment.annotations["Score"], 31.36)
        self.assertAlmostEqual(alignment.annotations["Identities"], 35)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.376)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 9.400)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                                                                                                                                                                                                                                              EEecCCHHHHHHHhcCCcEE                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(alignment.query.seq[238:258], "VYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(alignment[1], "VYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                                                                                                                                                                              I~~G~T~QDVl~~LG~P~~~                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                             l~~G~t~~eV~~~lG~P~~~                                                           ",
        )
        self.assertEqual(alignment.target.id, "5WAM_B")
        self.assertEqual(alignment.target.seq[29:49], "LRPGMTKDQVLLLLGSPILR")
        self.assertEqual(alignment.target.annotations["hmm_name"], "5WAM_B")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Outer membrane protein assembly factor; beta-barrel assembly machinery, BAM complex; 2.45A {Neisseria gonorrhoeae (strain ATCC 700825 / FA 1090)}; Related PDB entries: 5WAM_A",
        )
        self.assertEqual(alignment[0], "LRPGMTKDQVLLLLGSPILR")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "                             CCTTCBHHHHHHHHCCCSCC                                                           ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                             CCCCCCHHHHHHHHCCCcee                                                           ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 29,  49],
                             [238, 258]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['L', 'R', 'P', 'G', 'M', 'T', 'K', 'D', 'Q', 'V', 'L', 'L', 'L',
              'L', 'G', 'S', 'P', 'I', 'L', 'R'],
             ['V', 'Y', 'F', 'G', 'D', 'S', 'C', 'Q', 'D', 'V', 'L', 'S', 'M',
              'L', 'G', 'S', 'P', 'H', 'K', 'V']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 57.72)
        self.assertAlmostEqual(alignment.annotations["E-value"], 1.6e02)
        self.assertAlmostEqual(alignment.annotations["Score"], 26.08)
        self.assertAlmostEqual(alignment.annotations["Identities"], 11)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.028)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 9.400)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                  EEcCCcHHHHHHHHHHccCCcCcEEEEECCCCCcCCCEEEEeCCCCeEEEecCccceeeEEEEecCCCCEEEECCEEecCCCCCCcHHHHHHhhcCCCCCCccCcccEEEEEeCcEEEEee                                                                                                                                                                                                                                                                                           ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(
            alignment.query.seq[18:139],
            "FTLGMPLAQAVAILQKHCRIIKNVQVLYSEQSPLSHDLILNLTQDGIKLMFDAFNQRLKVIEVCDLTKVKLKYCGVHFNSQAIAPTIEQIDQSFGATHPGVYNSAEQLFHLNFRGLSFSFQ",
        )
        self.assertEqual(
            alignment[1],
            "FTLGMPLAQAVAILQKHCRIIKN--VQVLYSEQSPLSHDLILNLTQDGIKLMFDAFNQRLKVIEVCDLT-KVKLKYCGVHFNSQAIAPTIEQIDQSFGATHPGVYNS----------AEQLFHLNFR--------GLSFSFQ",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                  f~LG~sL~~vi~~Lr~~~~~~~~vel~Ys~~~pl~~~Ivi~L~~~GirL~Fd~~~QrL~lIEv~d~~~~~L~Y~~~~~~~~~~~ptf~~I~~~FGPTyPG~yd~~~~~Y~LsYpGisF~Fp                                                                                                                                                                                                                                                                                           ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "           I~~G~s~~eV~~~lG~p~~~~~~~~~~~~W~~~~~~~~~~~v~f~~~~kv~~k~~~~~~~~~~~~~t~~~~~~i~~Gmt~~~V~~~lG~p~~~~~~~~~~~~~~~~~~~~~y~~~~~~~~~~~~~~~~F~             ",
        )
        self.assertEqual(alignment.target.id, "3GMX_B")
        self.assertEqual(
            alignment.target.seq[11:141],
            "IQFGMDRTLVWQLAGADQSCSDQVERIICYNNPDHYGPQGHFFFNAADKLIHKRQMELFPAPKPTMRLATYNKTQTGMTEAQFWAAVPSDTCSALAEQYPNWPATNGNLREYVCPSKAERFAPSAYFTFT",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "3GMX_B")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "BLP; 2-layer alpha/beta sandwich, PROTEIN BINDING; HET: ACT; 1.05A {Streptomyces clavuligerus} SCOP: d.98.1.0; Related PDB entries: 3GMX_A 3GMY_A 3GMY_B",
        )
        self.assertEqual(
            alignment[0],
            "IQFGMDRTLVWQLAGADQSCSDQVERIICYNNPDHYGPQGHFFFNAAD----------KLIHKRQMELFPAPKPTMRLATYNKTQTGMTEAQFWAAVPS--DTCSALAEQYPNWPATNGNLREYVCPSKAERFAPSAYFTFT",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "           CCTTCBHHHHHHHHTHHHHEEECSSCEEEESSSCTTSSEEEEEECTTSBEEEEEEESSSCCSSCCCCHHHHTTCCTTCBHHHHHHHSCGGGCEEEEEECTTTTSCTTCEEEEEEESSSSTTCCEEEEEEE             ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "           CCCCCCHHHHHHHhCCCceecCCceeEEEEcCCCCCcceEEEEEccCCEEEEEEecccCCCCCCCCCHHHHhhcCCCCCHHHHHHHcCCccceeceeccCCCCCCCCeEEEEEeCCCCCCCCCeEEEEEE             ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 11,  34,  36,  59,  59,  70,  71, 100, 100, 106, 116, 126, 134,
                              141],
                             [ 18,  41,  41,  64,  74,  85,  85, 114, 116, 122, 122, 132, 132,
                              139]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['I', 'Q', 'F', 'G', 'M', 'D', 'R', 'T', 'L', 'V', 'W', 'Q', 'L',
              'A', 'G', 'A', 'D', 'Q', 'S', 'C', 'S', 'D', 'Q', 'V', 'E', 'R',
              'I', 'I', 'C', 'Y', 'N', 'N', 'P', 'D', 'H', 'Y', 'G', 'P', 'Q',
              'G', 'H', 'F', 'F', 'F', 'N', 'A', 'A', 'D', '-', '-', '-', '-',
              '-', '-', '-', '-', '-', '-', 'K', 'L', 'I', 'H', 'K', 'R', 'Q',
              'M', 'E', 'L', 'F', 'P', 'A', 'P', 'K', 'P', 'T', 'M', 'R', 'L',
              'A', 'T', 'Y', 'N', 'K', 'T', 'Q', 'T', 'G', 'M', 'T', 'E', 'A',
              'Q', 'F', 'W', 'A', 'A', 'V', 'P', 'S', '-', '-', 'D', 'T', 'C',
              'S', 'A', 'L', 'A', 'E', 'Q', 'Y', 'P', 'N', 'W', 'P', 'A', 'T',
              'N', 'G', 'N', 'L', 'R', 'E', 'Y', 'V', 'C', 'P', 'S', 'K', 'A',
              'E', 'R', 'F', 'A', 'P', 'S', 'A', 'Y', 'F', 'T', 'F', 'T'],
             ['F', 'T', 'L', 'G', 'M', 'P', 'L', 'A', 'Q', 'A', 'V', 'A', 'I',
              'L', 'Q', 'K', 'H', 'C', 'R', 'I', 'I', 'K', 'N', '-', '-', 'V',
              'Q', 'V', 'L', 'Y', 'S', 'E', 'Q', 'S', 'P', 'L', 'S', 'H', 'D',
              'L', 'I', 'L', 'N', 'L', 'T', 'Q', 'D', 'G', 'I', 'K', 'L', 'M',
              'F', 'D', 'A', 'F', 'N', 'Q', 'R', 'L', 'K', 'V', 'I', 'E', 'V',
              'C', 'D', 'L', 'T', '-', 'K', 'V', 'K', 'L', 'K', 'Y', 'C', 'G',
              'V', 'H', 'F', 'N', 'S', 'Q', 'A', 'I', 'A', 'P', 'T', 'I', 'E',
              'Q', 'I', 'D', 'Q', 'S', 'F', 'G', 'A', 'T', 'H', 'P', 'G', 'V',
              'Y', 'N', 'S', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',
              'A', 'E', 'Q', 'L', 'F', 'H', 'L', 'N', 'F', 'R', '-', '-', '-',
              '-', '-', '-', '-', '-', 'G', 'L', 'S', 'F', 'S', 'F', 'Q']],
            dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 57.07)
        self.assertAlmostEqual(alignment.annotations["E-value"], 24)
        self.assertAlmostEqual(alignment.annotations["Score"], 30.56)
        self.assertAlmostEqual(alignment.annotations["Identities"], 7)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.027)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 10.000)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                                                                                                                                                                                                                                              EEecCCHHHHHHHhcCCcEEEEeccccccccCCCCCcCCCCcCCccchhhHhhccceeecCCCCEEEEEEecCC                                                                                                              ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(
            alignment.query.seq[238:312],
            "VYFGDSCQDVLSMLGSPHKVFYKSEDKMKIHSPSPHKQVPSKCNDYFFNYFTLGVDILFDANTHKVKKFVLHTN",
        )
        self.assertEqual(
            alignment[1],
            "VYFGDSCQDVLSMLGSPHKVFYKSEDKMKIHSPSPHKQVPSKCNDYFFNYFTLGVDILFDANTHKVKKFVLHTN",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                                                                                                                                                                              I~~G~T~QDVl~~LG~P~~~f~K~ddrm~IH~~~~~~~~~~~~~dyF~NYF~lG~DiLfd~~t~~v~KiILHtN                                                                                                              ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                                                             i~vG~t~~~v~~~~G~p~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~fd~~~~~v~~I~~~~~         ",
        )
        self.assertEqual(alignment.target.id, "U5LAW6_9BACI/1")
        self.assertEqual(
            alignment.target.seq[61:131],
            "FHIGQPVSEIYSSVFIDTNINFQYKGSSYRFELSEDDLNTRPLIKAGNIYAQLYIDRFTGELSSIRYMDA",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "PF14504.6")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "; CAP_assoc_N ; CAP-associated N-terminal",
        )
        self.assertEqual(
            alignment[0],
            "FHIGQPVSEIYSSVFIDTNINFQYKGSSYRFELSEDD----LNTRPLIKAGNIYAQLYIDRFTGELSSIRYMDA",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                                                             ccCCCcHHHHHHhhcCCceEEEEEcCcEEEEEeCcchhccceeEEECCEEEEEEEECCCCEEEEEEEecH         ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 61,  98,  98, 131],
                             [238, 275, 279, 312]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['F', 'H', 'I', 'G', 'Q', 'P', 'V', 'S', 'E', 'I', 'Y', 'S', 'S',
              'V', 'F', 'I', 'D', 'T', 'N', 'I', 'N', 'F', 'Q', 'Y', 'K', 'G',
              'S', 'S', 'Y', 'R', 'F', 'E', 'L', 'S', 'E', 'D', 'D', '-', '-',
              '-', '-', 'L', 'N', 'T', 'R', 'P', 'L', 'I', 'K', 'A', 'G', 'N',
              'I', 'Y', 'A', 'Q', 'L', 'Y', 'I', 'D', 'R', 'F', 'T', 'G', 'E',
              'L', 'S', 'S', 'I', 'R', 'Y', 'M', 'D', 'A'],
             ['V', 'Y', 'F', 'G', 'D', 'S', 'C', 'Q', 'D', 'V', 'L', 'S', 'M',
              'L', 'G', 'S', 'P', 'H', 'K', 'V', 'F', 'Y', 'K', 'S', 'E', 'D',
              'K', 'M', 'K', 'I', 'H', 'S', 'P', 'S', 'P', 'H', 'K', 'Q', 'V',
              'P', 'S', 'K', 'C', 'N', 'D', 'Y', 'F', 'F', 'N', 'Y', 'F', 'T',
              'L', 'G', 'V', 'D', 'I', 'L', 'F', 'D', 'A', 'N', 'T', 'H', 'K',
              'V', 'K', 'K', 'F', 'V', 'L', 'H', 'T', 'N']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 56.37)
        self.assertAlmostEqual(alignment.annotations["E-value"], 14)
        self.assertAlmostEqual(alignment.annotations["Score"], 29.98)
        self.assertAlmostEqual(alignment.annotations["Identities"], 30)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.380)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 8.700)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                                                                                                                                                                                                                                              EEecCCHHHHHHHhcCCcEE                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(alignment.query.seq[238:258], "VYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(alignment[1], "VYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                                                                                                                                                                              I~~G~T~QDVl~~LG~P~~~                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "            l~~Gmtk~eV~~lLG~P~~~                                                        ",
        )
        self.assertEqual(alignment.target.id, "2YH9_C")
        self.assertEqual(alignment.target.seq[12:32], "IRVGMTQQQVAYALGTPLMS")
        self.assertEqual(alignment.target.annotations["hmm_name"], "2YH9_C")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "SMALL PROTEIN A; LIPOPROTEIN, 3D DOMAIN SWAP, MEMBRANE; HET: MSE; 1.8A {ESCHERICHIA COLI}; Related PDB entries: 2YH9_B 2YH9_A",
        )
        self.assertEqual(alignment[0], "IRVGMTQQQVAYALGTPLMS")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "            CCTTCBHHHHHHHHCSCSEE                                                        ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "            CCCCCCHHHHHHHhCCCcee                                                        ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 12,  32],
                             [238, 258]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['I', 'R', 'V', 'G', 'M', 'T', 'Q', 'Q', 'Q', 'V', 'A', 'Y', 'A',
              'L', 'G', 'T', 'P', 'L', 'M', 'S'],
             ['V', 'Y', 'F', 'G', 'D', 'S', 'C', 'Q', 'D', 'V', 'L', 'S', 'M',
              'L', 'G', 'S', 'P', 'H', 'K', 'V']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 56.24)
        self.assertAlmostEqual(alignment.annotations["E-value"], 30)
        self.assertAlmostEqual(alignment.annotations["Score"], 29.96)
        self.assertAlmostEqual(alignment.annotations["Identities"], 11)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.142)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 10.000)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                                                                                                                                                                                                                                                ecCCHHHHHHHhcCCcEEEEeccccccccCCCCCcCCCCcCCccchhhHhhccceeecCCCCEEEEEEecCC                                                                                                              ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(
            alignment.query.seq[240:312],
            "FGDSCQDVLSMLGSPHKVFYKSEDKMKIHSPSPHKQVPSKCNDYFFNYFTLGVDILFDANTHKVKKFVLHTN",
        )
        self.assertEqual(
            alignment[1],
            "FGDSCQDVLSMLGSPHKVFYKSEDKMKIHSPSPHKQVPSKCNDYFFNYFT-LGVDILFDANTHKVKKFVLHTN",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                                                                                                                                                                                ~G~T~QDVl~~LG~P~~~f~K~ddrm~IH~~~~~~~~~~~~~dyF~NYF~lG~DiLfd~~t~~v~KiILHtN                                                                                                              ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "lG~~~~~v~~~lG~P~~~~~~~~g~~~~~Y~~~~~~~~~v~~~~~~v~~i~~~~~                                                                                     ",
        )
        self.assertEqual(alignment.target.id, "U5LAW6_9BACI/1")
        self.assertEqual(
            alignment.target.seq[0:55],
            "IGKNASDLQVLLGDPERKDPSEYGYEWWIYKKGTSQYVQAGVLDGRIVTLFATGP",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "PF14504.6")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "; CAP_assoc_N ; CAP-associated N-terminal",
        )
        self.assertEqual(
            alignment[0],
            "IGKNASDLQVLLGDPERKDPSEYG------------------YEWWIYKKGTSQYVQAGVLDGRIVTLFATGP",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "CCCCHHHHHHHHCCCcEEeccccCCeEEEEecCCCcEEEEEEECCEEEEEEEcCC                                                                                     ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  0,  24,  24,  32,  33,  55],
                             [240, 264, 282, 290, 290, 312]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['I', 'G', 'K', 'N', 'A', 'S', 'D', 'L', 'Q', 'V', 'L', 'L', 'G',
              'D', 'P', 'E', 'R', 'K', 'D', 'P', 'S', 'E', 'Y', 'G', '-', '-',
              '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',
              '-', '-', '-', 'Y', 'E', 'W', 'W', 'I', 'Y', 'K', 'K', 'G', 'T',
              'S', 'Q', 'Y', 'V', 'Q', 'A', 'G', 'V', 'L', 'D', 'G', 'R', 'I',
              'V', 'T', 'L', 'F', 'A', 'T', 'G', 'P'],
             ['F', 'G', 'D', 'S', 'C', 'Q', 'D', 'V', 'L', 'S', 'M', 'L', 'G',
              'S', 'P', 'H', 'K', 'V', 'F', 'Y', 'K', 'S', 'E', 'D', 'K', 'M',
              'K', 'I', 'H', 'S', 'P', 'S', 'P', 'H', 'K', 'Q', 'V', 'P', 'S',
              'K', 'C', 'N', 'D', 'Y', 'F', 'F', 'N', 'Y', 'F', 'T', '-', 'L',
              'G', 'V', 'D', 'I', 'L', 'F', 'D', 'A', 'N', 'T', 'H', 'K', 'V',
              'K', 'K', 'F', 'V', 'L', 'H', 'T', 'N']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 56.21)
        self.assertAlmostEqual(alignment.annotations["E-value"], 14)
        self.assertAlmostEqual(alignment.annotations["Score"], 31.47)
        self.assertAlmostEqual(alignment.annotations["Identities"], 40)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.433)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 7.400)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                                                                                                                                                                                                                                              EEecCCHHHHHHHhcCCcEE                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(alignment.query.seq[238:258], "VYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(alignment[1], "VYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                                                                                                                                                                              I~~G~T~QDVl~~LG~P~~~                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                       I~~GmTk~qV~~~lG~P~~~                                                            ",
        )
        self.assertEqual(alignment.target.id, "4DM5_A")
        self.assertEqual(alignment.target.seq[23:43], "VEKGMSQQEVLRIGGTPSGT")
        self.assertEqual(alignment.target.annotations["hmm_name"], "4DM5_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Osmotically inducible lipoprotein OsmE; lipoprotein, UNKNOWN FUNCTION; HET: MSE; 1.5A {Pseudomonas aeruginosa PAO1}; Related PDB entries: 4DM5_C 4DM5_D 4DM5_B",
        )
        self.assertEqual(alignment[0], "VEKGMSQQEVLRIGGTPSGT")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "                       CCTTCBHHHHHHHHCSCSEE                                                            ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                       CCCCCCHHHHHHHhCCCCee                                                            ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 23,  43],
                             [238, 258]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['V', 'E', 'K', 'G', 'M', 'S', 'Q', 'Q', 'E', 'V', 'L', 'R', 'I',
              'G', 'G', 'T', 'P', 'S', 'G', 'T'],
             ['V', 'Y', 'F', 'G', 'D', 'S', 'C', 'Q', 'D', 'V', 'L', 'S', 'M',
              'L', 'G', 'S', 'P', 'H', 'K', 'V']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 56.09)
        self.assertAlmostEqual(alignment.annotations["E-value"], 77)
        self.assertAlmostEqual(alignment.annotations["Score"], 27.35)
        self.assertAlmostEqual(alignment.annotations["Identities"], 22)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.238)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 10.000)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "            CCCCccEEcCCcHHHHHHHHHHccCCcCcEEEEECCCCCcCCCEEEEeCCCCeEEEecCccceeeEEEEec                                                                                                                                                                                                                                                                                                                                                   ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(
            alignment.query.seq[12:83],
            "GNEQWEFTLGMPLAQAVAILQKHCRIIKNVQVLYSEQSPLSHDLILNLTQDGIKLMFDAFNQRLKVIEVCD",
        )
        self.assertEqual(
            alignment[1],
            "GNEQWEFTLGMPLAQAVAILQKHCRI-----IKNVQVLYSEQSPLSHDLILNLTQDGIKLMFDAFNQRLKVIEVCD",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "            G~~~~~f~LG~sL~~vi~~Lr~~~~~~~~vel~Ys~~~pl~~~Ivi~L~~~GirL~Fd~~~QrL~lIEv~d                                                                                                                                                                                                                                                                                                                                                   ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                                                            gi~vG~t~~~v~~~~G~p~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~fd~~~~~v~~I~~~~          ",
        )
        self.assertEqual(alignment.target.id, "U5LAW6_9BACI/1")
        self.assertEqual(
            alignment.target.seq[60:130],
            "PFHIGQPVSEIYSSVFIDTNINFQYKGSSYRFELSEDDLNTRPLIKAGNIYAQLYIDRFTGELSSIRYMD",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "PF14504.6")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "; CAP_assoc_N ; CAP-associated N-terminal",
        )
        self.assertEqual(
            alignment[0],
            "P-----FHIGQPVSEIYSSVFIDTNINFQYKGSSYRFELSEDDLNTRPLIKA-GNIYAQLYIDRFTGELSSIRYMD",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                                                            CccCCCcHHHHHHhhcCCceEEEEEcCcEEEEEeCcchhccceeEEECCEEEEEEEECCCCEEEEEEEec          ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 60,  61,  61,  81,  86, 107, 107, 130],
                             [ 12,  13,  18,  38,  38,  59,  60,  83]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['P', '-', '-', '-', '-', '-', 'F', 'H', 'I', 'G', 'Q', 'P', 'V',
              'S', 'E', 'I', 'Y', 'S', 'S', 'V', 'F', 'I', 'D', 'T', 'N', 'I',
              'N', 'F', 'Q', 'Y', 'K', 'G', 'S', 'S', 'Y', 'R', 'F', 'E', 'L',
              'S', 'E', 'D', 'D', 'L', 'N', 'T', 'R', 'P', 'L', 'I', 'K', 'A',
              '-', 'G', 'N', 'I', 'Y', 'A', 'Q', 'L', 'Y', 'I', 'D', 'R', 'F',
              'T', 'G', 'E', 'L', 'S', 'S', 'I', 'R', 'Y', 'M', 'D'],
             ['G', 'N', 'E', 'Q', 'W', 'E', 'F', 'T', 'L', 'G', 'M', 'P', 'L',
              'A', 'Q', 'A', 'V', 'A', 'I', 'L', 'Q', 'K', 'H', 'C', 'R', 'I',
              '-', '-', '-', '-', '-', 'I', 'K', 'N', 'V', 'Q', 'V', 'L', 'Y',
              'S', 'E', 'Q', 'S', 'P', 'L', 'S', 'H', 'D', 'L', 'I', 'L', 'N',
              'L', 'T', 'Q', 'D', 'G', 'I', 'K', 'L', 'M', 'F', 'D', 'A', 'F',
              'N', 'Q', 'R', 'L', 'K', 'V', 'I', 'E', 'V', 'C', 'D']],
            dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 53.82)
        self.assertAlmostEqual(alignment.annotations["E-value"], 2.4e02)
        self.assertAlmostEqual(alignment.annotations["Score"], 24.83)
        self.assertAlmostEqual(alignment.annotations["Identities"], 8)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.090)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 8.400)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "        ccccCCCCccEEcCCcHHHHHHHHHHccCCcCcEEEEECCCCCcCCCEEEEeCCCCeEEEecCccceeeEEEEecCCCCEEEECCEEecCCCCCCcHHHHHHhhcCCCCCCccCcccEEEEEeCcEEEEeeCcccCCCCccCCCcccccccCCCCCCCeEEEEEEE                                                                                                                                                                                                                                                        ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(
            alignment.query.seq[8:174],
            "ERSLGNEQWEFTLGMPLAQAVAILQKHCRIIKNVQVLYSEQSPLSHDLILNLTQDGIKLMFDAFNQRLKVIEVCDLTKVKLKYCGVHFNSQAIAPTIEQIDQSFGATHPGVYNSAEQLFHLNFRGLSFSFQLDSWTEAPKYEPNFAHGLASLQIPHGATVKRMYIY",
        )
        self.assertEqual(
            alignment[1],
            "ERSLGNEQWEFTLGMP-LAQAVAILQKHCRIIKNVQVLYSEQSPLSHDLILNLTQDGIKLMFDAFNQRLKVIEVCDLTKVKLKYCGVHFNSQAIAPTIEQIDQSFGATHPGVYNSAEQLFHLNF-RG----LSFSFQLDSWTEAPKYEPNFAHGLASLQIPHGATVKRMYIY",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "        g~~LG~~~~~f~LG~sL~~vi~~Lr~~~~~~~~vel~Ys~~~pl~~~Ivi~L~~~GirL~Fd~~~QrL~lIEv~d~~~~~L~Y~~~~~~~~~~~ptf~~I~~~FGPTyPG~yd~~~~~Y~LsYpGisF~Fpi~~~~~~~~~~~~~~~~l~~l~~~~~~~~s~i~If                                                                                                                                                                                                                                                        ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "          Gk~~g~~~~lG~s~~~~V~~~~G~P~~~~~~~~~~g~~~~Y~~~~~~f~~~~~~~v~~i~~~~~~~~~lt~~~v~~~LG~P~~~~~~~~~~~~~Y~~g~~y~L~f~~~~~~~~~~~~~v~~i~v~",
        )
        self.assertEqual(alignment.target.id, "R4K5N0_CLOPA/8")
        self.assertEqual(
            alignment.target.seq[10:135],
            "GKVFNSDFPAKDTNIDSVESKWGKADNSEWVASAKGLYSTYSKHNIVFGSNKGGQIFEVRSLDKQLGNIYLSMVKDKLGTPQHDVKVNGEEIIGYKMGNDFKILFVFPEPTNQHANPIMSHYSVL",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "PF14172.6")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "; DUF4309 ; Domain of unknown function (DUF4309)",
        )
        self.assertEqual(
            alignment[0],
            "GKVFNS---DFPAKDTNIDSVESKWGK------ADNSEWVASA---KGLYSTYSKHNIVFGSNKGG-QIFEVRSLDKQLGN--------------IYLSMVKDKLG--TPQHDVKVNGEEIIGYKMGNDFKILFVFPEPTNQ------------------HANPIMSHYSVL",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "          CCCCCCCCCCCCCcHHHHHHHhCCCCceehhhhcCceeEEecCceEEEEECCCCeEEEEEEechhhcCCCHHHHHHHhCCCcEEEEECCEEEEEEEeCCCEEEEEEeCCCCCCCCCCeeeEEEeC",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 10,  16,  16,  23,  24,  34,  34,  44,  44,  64,  64,  78,  78,
                               89,  89, 105, 106, 108, 112, 123, 123, 135],
                             [  8,  14,  17,  24,  24,  34,  40,  50,  53,  73,  74,  88, 102,
                              113, 115, 131, 131, 133, 133, 144, 162, 174]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['G', 'K', 'V', 'F', 'N', 'S', '-', '-', '-', 'D', 'F', 'P', 'A',
              'K', 'D', 'T', 'N', 'I', 'D', 'S', 'V', 'E', 'S', 'K', 'W', 'G',
              'K', '-', '-', '-', '-', '-', '-', 'A', 'D', 'N', 'S', 'E', 'W',
              'V', 'A', 'S', 'A', '-', '-', '-', 'K', 'G', 'L', 'Y', 'S', 'T',
              'Y', 'S', 'K', 'H', 'N', 'I', 'V', 'F', 'G', 'S', 'N', 'K', 'G',
              'G', '-', 'Q', 'I', 'F', 'E', 'V', 'R', 'S', 'L', 'D', 'K', 'Q',
              'L', 'G', 'N', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',
              '-', '-', '-', '-', 'I', 'Y', 'L', 'S', 'M', 'V', 'K', 'D', 'K',
              'L', 'G', '-', '-', 'T', 'P', 'Q', 'H', 'D', 'V', 'K', 'V', 'N',
              'G', 'E', 'E', 'I', 'I', 'G', 'Y', 'K', 'M', 'G', 'N', 'D', 'F',
              'K', 'I', 'L', 'F', 'V', 'F', 'P', 'E', 'P', 'T', 'N', 'Q', '-',
              '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',
              '-', '-', '-', '-', 'H', 'A', 'N', 'P', 'I', 'M', 'S', 'H', 'Y',
              'S', 'V', 'L'],
             ['E', 'R', 'S', 'L', 'G', 'N', 'E', 'Q', 'W', 'E', 'F', 'T', 'L',
              'G', 'M', 'P', '-', 'L', 'A', 'Q', 'A', 'V', 'A', 'I', 'L', 'Q',
              'K', 'H', 'C', 'R', 'I', 'I', 'K', 'N', 'V', 'Q', 'V', 'L', 'Y',
              'S', 'E', 'Q', 'S', 'P', 'L', 'S', 'H', 'D', 'L', 'I', 'L', 'N',
              'L', 'T', 'Q', 'D', 'G', 'I', 'K', 'L', 'M', 'F', 'D', 'A', 'F',
              'N', 'Q', 'R', 'L', 'K', 'V', 'I', 'E', 'V', 'C', 'D', 'L', 'T',
              'K', 'V', 'K', 'L', 'K', 'Y', 'C', 'G', 'V', 'H', 'F', 'N', 'S',
              'Q', 'A', 'I', 'A', 'P', 'T', 'I', 'E', 'Q', 'I', 'D', 'Q', 'S',
              'F', 'G', 'A', 'T', 'H', 'P', 'G', 'V', 'Y', 'N', 'S', 'A', 'E',
              'Q', 'L', 'F', 'H', 'L', 'N', 'F', '-', 'R', 'G', '-', '-', '-',
              '-', 'L', 'S', 'F', 'S', 'F', 'Q', 'L', 'D', 'S', 'W', 'T', 'E',
              'A', 'P', 'K', 'Y', 'E', 'P', 'N', 'F', 'A', 'H', 'G', 'L', 'A',
              'S', 'L', 'Q', 'I', 'P', 'H', 'G', 'A', 'T', 'V', 'K', 'R', 'M',
              'Y', 'I', 'Y']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 52.88)
        self.assertAlmostEqual(alignment.annotations["E-value"], 16)
        self.assertAlmostEqual(alignment.annotations["Score"], 31.41)
        self.assertAlmostEqual(alignment.annotations["Identities"], 20)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.320)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 8.800)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                                                                                                                                                                                                                                              EEecCCHHHHHHHhcCCcEE                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(alignment.query.seq[238:258], "VYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(alignment[1], "VYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                                                                                                                                                                              I~~G~T~QDVl~~LG~P~~~                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "       I~~G~t~~~V~~~lG~P~~~                                                                                              ",
        )
        self.assertEqual(alignment.target.id, "B5GLC0_STRC2/3")
        self.assertEqual(alignment.target.seq[7:27], "IQFGMDRTLVWQLAGADQSC")
        self.assertEqual(alignment.target.annotations["hmm_name"], "PF07467.11")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "; BLIP ; Beta-lactamase inhibitor (BLIP)",
        )
        self.assertEqual(alignment[0], "IQFGMDRTLVWQLAGADQSC")
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "       CCCCCCHHHHHHHHCCCcee                                                                                              ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  7,  27],
                             [238, 258]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['I', 'Q', 'F', 'G', 'M', 'D', 'R', 'T', 'L', 'V', 'W', 'Q', 'L',
              'A', 'G', 'A', 'D', 'Q', 'S', 'C'],
             ['V', 'Y', 'F', 'G', 'D', 'S', 'C', 'Q', 'D', 'V', 'L', 'S', 'M',
              'L', 'G', 'S', 'P', 'H', 'K', 'V']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 50.88)
        self.assertAlmostEqual(alignment.annotations["E-value"], 27)
        self.assertAlmostEqual(alignment.annotations["Score"], 28.63)
        self.assertAlmostEqual(alignment.annotations["Identities"], 19)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.145)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 8.700)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                                                                                                                                                                                                                                              EEecCCHHHHHHHhcCCcEEEEecccc                                                                                                                                                             ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(alignment.query.seq[238:265], "VYFGDSCQDVLSMLGSPHKVFYKSEDK")
        self.assertEqual(alignment[1], "VYFGDSCQDVLSMLGSPHKVFYKSEDK")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                                                                                                                                                                              I~~G~T~QDVl~~LG~P~~~f~K~ddr                                                                                                                                                             ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "          V~~Gds~~~Vl~~cG~P~~~~~~~~~~                                                         ",
        )
        self.assertEqual(alignment.target.id, "B3PEX3_CELJU/2")
        self.assertEqual(alignment.target.seq[10:37], "TQTGDTKAEVIAKCGDPVFTDHYCAPM")
        self.assertEqual(alignment.target.annotations["hmm_name"], "PF11006.8")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "; DUF2845 ; Protein of unknown function (DUF2845)",
        )
        self.assertEqual(alignment[0], "TQTGDTKAEVIAKCGDPVFTDHYCAPM")
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "          cCCCCCHHHHHHHhCCCcEEEeEeeec                                                         ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 10,  37],
                             [238, 265]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['T', 'Q', 'T', 'G', 'D', 'T', 'K', 'A', 'E', 'V', 'I', 'A', 'K',
              'C', 'G', 'D', 'P', 'V', 'F', 'T', 'D', 'H', 'Y', 'C', 'A', 'P',
              'M'],
             ['V', 'Y', 'F', 'G', 'D', 'S', 'C', 'Q', 'D', 'V', 'L', 'S', 'M',
              'L', 'G', 'S', 'P', 'H', 'K', 'V', 'F', 'Y', 'K', 'S', 'E', 'D',
              'K']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 49.44)
        self.assertAlmostEqual(alignment.annotations["E-value"], 1.1e02)
        self.assertAlmostEqual(alignment.annotations["Score"], 27.28)
        self.assertAlmostEqual(alignment.annotations["Identities"], 9)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.155)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 9.600)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "    EeccccccCCCCccEEcCCcHHHHHHHHHHccCCcCcEEEEECCCCCcCCCEEEEeCCCCeEEEecCccceeeEEEEecCCCCEEEECCEEecCCCCCCcHHHHHHhhcCCCCCCccCcccEEEEEeCcEEEEe                                                                                                                                                                                                                                                                                            ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(
            alignment.query.seq[4:138],
            "EVVPERSLGNEQWEFTLGMPLAQAVAILQKHCRIIKNVQVLYSEQSPLSHDLILNLTQDGIKLMFDAFNQRLKVIEVCDLTKVKLKYCGVHFNSQAIAPTIEQIDQSFGATHPGVYNSAEQLFHLNFRGLSFSF",
        )
        self.assertEqual(
            alignment[1],
            "EVVPERSLGNEQWEFTLGMPLAQAVAILQKHCRIIKNV-------QVLYSEQSPLSHDLILNLTQDGIKLMFDAFNQRLKVIEVCDLTKV-KLKYCGVHFNSQAIAPTIEQIDQSFGATHPGVYNSAE----QLFHLNFR----------GLSFSF",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "    ~i~Pg~~LG~~~~~f~LG~sL~~vi~~Lr~~~~~~~~vel~Ys~~~pl~~~Ivi~L~~~GirL~Fd~~~QrL~lIEv~d~~~~~L~Y~~~~~~~~~~~ptf~~I~~~FGPTyPG~yd~~~~~Y~LsYpGisF~F                                                                                                                                                                                                                                                                                            ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "          ~i~~~~~~~~~~~G~t~~eV~~~lG~p~~~~~~~~~~~~~~~~~w~~~~~~~~v~~~~~~~~~~~~~~~~~~~~~~~t~~~~~~i~~G~t~~eV~~~lG~p~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~F             ",
        )
        self.assertEqual(alignment.target.id, "Q8DTX1_STRMU/4")
        self.assertEqual(
            alignment.target.seq[10:146],
            "KIKVTTDQNHFSGGTSIEQLKQWFGDPNKSEQRNAGNITLDSYTWVKDGAVINAQLYKNSTVARSISNFSFSREAKIGKEDYDELKIGESYKKVVEKLGEPDVLSQSMSSDKEEMQTVWSSGIKTKSSSATIELYF",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "PF12978.7")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "; DUF3862 ; Domain of Unknown Function with PDB structure (DUF3862)",
        )
        self.assertEqual(
            alignment[0],
            "KIKVTTDQNH----FSGGTSIEQLKQWFGDPNKSEQRNAGNITLDSYTWVK------------DGAVINAQLY--KNSTVARSISNFSFSREAKIGKEDYDELKIGESYKKVVEKLG--EPDVLSQSMSSDKEEMQTVWSSGIKTKSSSATIELYF",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "          ceeecCCcccCCCCCCHHHHHHHHCCCceEEeeeeCCeEEEEEEEEeCCeEEEEEEECCEEEEEEEEeceecCCCCCCHHHHHhcCCCCCHHHHHHHHCCCCeeEEEeeCCceEEEEEEEecccCCCCCcEEEEEE             ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 10,  20,  20,  44,  51,  57,  57,  67,  67,  82,  83, 109, 109,
                              118, 122, 130, 140, 146],
                             [  4,  14,  18,  42,  42,  48,  60,  70,  72,  87,  87, 113, 115,
                              124, 124, 132, 132, 138]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['K', 'I', 'K', 'V', 'T', 'T', 'D', 'Q', 'N', 'H', '-', '-', '-',
              '-', 'F', 'S', 'G', 'G', 'T', 'S', 'I', 'E', 'Q', 'L', 'K', 'Q',
              'W', 'F', 'G', 'D', 'P', 'N', 'K', 'S', 'E', 'Q', 'R', 'N', 'A',
              'G', 'N', 'I', 'T', 'L', 'D', 'S', 'Y', 'T', 'W', 'V', 'K', '-',
              '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', 'D', 'G',
              'A', 'V', 'I', 'N', 'A', 'Q', 'L', 'Y', '-', '-', 'K', 'N', 'S',
              'T', 'V', 'A', 'R', 'S', 'I', 'S', 'N', 'F', 'S', 'F', 'S', 'R',
              'E', 'A', 'K', 'I', 'G', 'K', 'E', 'D', 'Y', 'D', 'E', 'L', 'K',
              'I', 'G', 'E', 'S', 'Y', 'K', 'K', 'V', 'V', 'E', 'K', 'L', 'G',
              '-', '-', 'E', 'P', 'D', 'V', 'L', 'S', 'Q', 'S', 'M', 'S', 'S',
              'D', 'K', 'E', 'E', 'M', 'Q', 'T', 'V', 'W', 'S', 'S', 'G', 'I',
              'K', 'T', 'K', 'S', 'S', 'S', 'A', 'T', 'I', 'E', 'L', 'Y', 'F'],
             ['E', 'V', 'V', 'P', 'E', 'R', 'S', 'L', 'G', 'N', 'E', 'Q', 'W',
              'E', 'F', 'T', 'L', 'G', 'M', 'P', 'L', 'A', 'Q', 'A', 'V', 'A',
              'I', 'L', 'Q', 'K', 'H', 'C', 'R', 'I', 'I', 'K', 'N', 'V', '-',
              '-', '-', '-', '-', '-', '-', 'Q', 'V', 'L', 'Y', 'S', 'E', 'Q',
              'S', 'P', 'L', 'S', 'H', 'D', 'L', 'I', 'L', 'N', 'L', 'T', 'Q',
              'D', 'G', 'I', 'K', 'L', 'M', 'F', 'D', 'A', 'F', 'N', 'Q', 'R',
              'L', 'K', 'V', 'I', 'E', 'V', 'C', 'D', 'L', 'T', 'K', 'V', '-',
              'K', 'L', 'K', 'Y', 'C', 'G', 'V', 'H', 'F', 'N', 'S', 'Q', 'A',
              'I', 'A', 'P', 'T', 'I', 'E', 'Q', 'I', 'D', 'Q', 'S', 'F', 'G',
              'A', 'T', 'H', 'P', 'G', 'V', 'Y', 'N', 'S', 'A', 'E', '-', '-',
              '-', '-', 'Q', 'L', 'F', 'H', 'L', 'N', 'F', 'R', '-', '-', '-',
              '-', '-', '-', '-', '-', '-', '-', 'G', 'L', 'S', 'F', 'S', 'F']],
            dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 48.69)
        self.assertAlmostEqual(alignment.annotations["E-value"], 1.6e02)
        self.assertAlmostEqual(alignment.annotations["Score"], 25.16)
        self.assertAlmostEqual(alignment.annotations["Identities"], 18)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.331)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 4.200)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                      CcHHHHHHHHHHccCCcCcEEEEECCCCCcCCCEEEEeCCCCeEEEecCccceeeEEEEe                                                                                                                                                                                                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(
            alignment.query.seq[22:82],
            "MPLAQAVAILQKHCRIIKNVQVLYSEQSPLSHDLILNLTQDGIKLMFDAFNQRLKVIEVC",
        )
        self.assertEqual(
            alignment[1],
            "MPLAQAVAILQKHCRIIKNVQVLYSEQSPLSHDLILNL----TQDGIKLMFDAFNQRLK-VIEVC",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                      ~sL~~vi~~Lr~~~~~~~~vel~Ys~~~pl~~~Ivi~L~~~GirL~Fd~~~QrL~lIEv~                                                                                                                                                                                                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "    m~~dE~ieYik~nVk~~D~lEisynRvf~pGeVl~i~~~~~~~~~~~~v~l~l~ge~~~~~veiD                       ",
        )
        self.assertEqual(alignment.target.id, "A5UN43_METS3/2")
        self.assertEqual(
            alignment.target.seq[4:69],
            "LTPDKAVEYLKDNVKIHDNLEISYNRIFGSGEVLNMDFSEYFGKPGFKMLLSLDGDSINPTIEID",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "PF09870.9")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "; DUF2097 ; Uncharacterized protein conserved in archaea (DUF2097)",
        )
        self.assertEqual(
            alignment[0],
            "LTPDKAVEYLKDNVKIHDNLEISYNRIFGSGEVLNMDFSEYFGKPGFKMLLSLDGDSINPTIEID",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "    cCHHHHHHHHHHhCCccCeEEEEeceeecceeEEEEeeceecCCCeEEEEEEeCCCeecceEEEe                       ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 4, 42, 46, 63, 64, 69],
                             [22, 60, 60, 77, 77, 82]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['L', 'T', 'P', 'D', 'K', 'A', 'V', 'E', 'Y', 'L', 'K', 'D', 'N',
              'V', 'K', 'I', 'H', 'D', 'N', 'L', 'E', 'I', 'S', 'Y', 'N', 'R',
              'I', 'F', 'G', 'S', 'G', 'E', 'V', 'L', 'N', 'M', 'D', 'F', 'S',
              'E', 'Y', 'F', 'G', 'K', 'P', 'G', 'F', 'K', 'M', 'L', 'L', 'S',
              'L', 'D', 'G', 'D', 'S', 'I', 'N', 'P', 'T', 'I', 'E', 'I', 'D'],
             ['M', 'P', 'L', 'A', 'Q', 'A', 'V', 'A', 'I', 'L', 'Q', 'K', 'H',
              'C', 'R', 'I', 'I', 'K', 'N', 'V', 'Q', 'V', 'L', 'Y', 'S', 'E',
              'Q', 'S', 'P', 'L', 'S', 'H', 'D', 'L', 'I', 'L', 'N', 'L', '-',
              '-', '-', '-', 'T', 'Q', 'D', 'G', 'I', 'K', 'L', 'M', 'F', 'D',
              'A', 'F', 'N', 'Q', 'R', 'L', 'K', '-', 'V', 'I', 'E', 'V', 'C']],
            dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 43.68)
        self.assertAlmostEqual(alignment.annotations["E-value"], 24)
        self.assertAlmostEqual(alignment.annotations["Score"], 30.46)
        self.assertAlmostEqual(alignment.annotations["Identities"], 35)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.526)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 8.300)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                                                                                                                                                                                                                                              EEecCCHHHHHHHhcCCcEE                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(alignment.query.seq[238:258], "VYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(alignment[1], "VYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                                                                                                                                                                              I~~G~T~QDVl~~LG~P~~~                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                         I~~G~Tk~eV~~~LG~P~~~                                                                         ",
        )
        self.assertEqual(alignment.target.id, "2PXG_A")
        self.assertEqual(alignment.target.seq[25:45], "LQVGQSKQQVSALLGTPSIP")
        self.assertEqual(alignment.target.annotations["hmm_name"], "2PXG_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Outer membrane protein; two layer alpha/beta plait, two; NMR {Xanthomonas axonopodis pv. citri}",
        )
        self.assertEqual(alignment[0], "LQVGQSKQQVSALLGTPSIP")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "                         TCCCSCHHHHHHHHTSCCCC                                                                         ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                         CCCCCCHHHHHHHhCCCccC                                                                         ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 25,  45],
                             [238, 258]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['L', 'Q', 'V', 'G', 'Q', 'S', 'K', 'Q', 'Q', 'V', 'S', 'A', 'L',
              'L', 'G', 'T', 'P', 'S', 'I', 'P'],
             ['V', 'Y', 'F', 'G', 'D', 'S', 'C', 'Q', 'D', 'V', 'L', 'S', 'M',
              'L', 'G', 'S', 'P', 'H', 'K', 'V']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 42.97)
        self.assertAlmostEqual(alignment.annotations["E-value"], 29)
        self.assertAlmostEqual(alignment.annotations["Score"], 30.94)
        self.assertAlmostEqual(alignment.annotations["Identities"], 20)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.320)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 9.400)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                                                                                                                                                                                                                                              EEecCCHHHHHHHhcCCcEE                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(alignment.query.seq[238:258], "VYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(alignment[1], "VYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                                                                                                                                                                              I~~G~T~QDVl~~LG~P~~~                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "           I~~G~s~~eV~~~lG~p~~~                                                                                                                           ",
        )
        self.assertEqual(alignment.target.id, "3GMX_B")
        self.assertEqual(alignment.target.seq[11:31], "IQFGMDRTLVWQLAGADQSC")
        self.assertEqual(alignment.target.annotations["hmm_name"], "3GMX_B")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "BLP; 2-layer alpha/beta sandwich, PROTEIN BINDING; HET: ACT; 1.05A {Streptomyces clavuligerus} SCOP: d.98.1.0; Related PDB entries: 3GMX_A 3GMY_A 3GMY_B",
        )
        self.assertEqual(alignment[0], "IQFGMDRTLVWQLAGADQSC")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "           CCTTCBHHHHHHHHTHHHHE                                                                                                                           ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "           CCCCCCHHHHHHHhCCCcee                                                                                                                           ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 11,  31],
                             [238, 258]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['I', 'Q', 'F', 'G', 'M', 'D', 'R', 'T', 'L', 'V', 'W', 'Q', 'L',
              'A', 'G', 'A', 'D', 'Q', 'S', 'C'],
             ['V', 'Y', 'F', 'G', 'D', 'S', 'C', 'Q', 'D', 'V', 'L', 'S', 'M',
              'L', 'G', 'S', 'P', 'H', 'K', 'V']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 41.38)
        self.assertAlmostEqual(alignment.annotations["E-value"], 33)
        self.assertAlmostEqual(alignment.annotations["Score"], 29.00)
        self.assertAlmostEqual(alignment.annotations["Identities"], 30)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.380)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 9.900)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                                                                                                                                                                                                                                              EEecCCHHHHHHHhcCCcEE                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(alignment.query.seq[238:258], "VYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(alignment[1], "VYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                                                                                                                                                                              I~~G~T~QDVl~~LG~P~~~                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                                             l~~G~t~~eV~~~lG~P~~~                                                          ",
        )
        self.assertEqual(alignment.target.id, "5D0O_E")
        self.assertEqual(alignment.target.seq[45:65], "IRVGMTQQQVAYALGTPLMS")
        self.assertEqual(alignment.target.annotations["hmm_name"], "5D0O_E")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Outer membrane protein assembly factor; E.coli, Bacterial outer membrane beta; 2.9A {Escherichia coli}; Related PDB entries: 5D0Q_E 5D0Q_I",
        )
        self.assertEqual(alignment[0], "IRVGMTQQQVAYALGTPLMS")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "                                             CCSSCBHHHHHHHHCCCSBC                                                          ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                                             CCCCCCHHHHHHHhCCCcee                                                          ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 45,  65],
                             [238, 258]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['I', 'R', 'V', 'G', 'M', 'T', 'Q', 'Q', 'Q', 'V', 'A', 'Y', 'A',
              'L', 'G', 'T', 'P', 'L', 'M', 'S'],
             ['V', 'Y', 'F', 'G', 'D', 'S', 'C', 'Q', 'D', 'V', 'L', 'S', 'M',
              'L', 'G', 'S', 'P', 'H', 'K', 'V']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 39.07)
        self.assertAlmostEqual(alignment.annotations["E-value"], 83)
        self.assertAlmostEqual(alignment.annotations["Score"], 29.23)
        self.assertAlmostEqual(alignment.annotations["Identities"], 12)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.176)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 8.500)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                                                                                                                                                                                                                                              EEecCCHHHHHHHhcCCcEEEEeccccccccCCCCCcCCCCcCCccchhhHhhccceeecCCCCEEEEEEec                                                                                                                ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(
            alignment.query.seq[238:310],
            "VYFGDSCQDVLSMLGSPHKVFYKSEDKMKIHSPSPHKQVPSKCNDYFFNYFTLGVDILFDANTHKVKKFVLH",
        )
        self.assertEqual(
            alignment[1],
            "VYF---------GDSCQDVLSMLGSPHKVFYKSEDKMKIHSPSPHKQVPSKCNDYFFNYFTLGVDILFDANTHKVKKFVLH",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                                                                                                                                                                              I~~G~T~QDVl~~LG~P~~~f~K~ddrm~IH~~~~~~~~~~~~~dyF~NYF~lG~DiLfd~~t~~v~KiILH                                                                                                                ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                              I~~g~~~~~~~~Gmt~~eV~~ilG~P~~~~~~~~~~~~~~~~~W~~~~~~i~v~F~~~kv~~k~~~~                                                                                  ",
        )
        self.assertEqual(alignment.target.id, "3D4E_A")
        self.assertEqual(
            alignment.target.seq[30:97],
            "IKVTTDQNHFSGGTSIEQLKQWFGDPNKSEQRNAGNITLDSYTWVKDGAVINAQLYKNSTVARSISN",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "3D4E_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "putative beta-lactamase inhibitor protein; NP_721579.1, putative beta-lactamase inhibitor protein; HET: EDO, MSE; 1.4A {Streptococcus mutans}",
        )
        self.assertEqual(
            alignment[0],
            "IKVTTDQNHFSGGTSIEQLKQWFGDPNKSEQRNAG-------------NITLDSYTWVKDGAVINAQLY-KNSTVARSISN",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "                              CCCCCCTTTCCSSCCHHHHHHHHCSCSEEEEEEETTEEEEEEEEEETTEEEEEEEETTEEEEEEEES                                                                                  ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                              eeecCCccCCCCCCCHHHHHHHHCCCceeEEEeeCCEEEEEEEEEeCCeEEEEEEECCeEEEeeeec                                                                                  ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 30,  33,  42,  65,  65,  86,  86,  97],
                             [238, 241, 241, 264, 277, 298, 299, 310]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['I', 'K', 'V', 'T', 'T', 'D', 'Q', 'N', 'H', 'F', 'S', 'G', 'G',
              'T', 'S', 'I', 'E', 'Q', 'L', 'K', 'Q', 'W', 'F', 'G', 'D', 'P',
              'N', 'K', 'S', 'E', 'Q', 'R', 'N', 'A', 'G', '-', '-', '-', '-',
              '-', '-', '-', '-', '-', '-', '-', '-', '-', 'N', 'I', 'T', 'L',
              'D', 'S', 'Y', 'T', 'W', 'V', 'K', 'D', 'G', 'A', 'V', 'I', 'N',
              'A', 'Q', 'L', 'Y', '-', 'K', 'N', 'S', 'T', 'V', 'A', 'R', 'S',
              'I', 'S', 'N'],
             ['V', 'Y', 'F', '-', '-', '-', '-', '-', '-', '-', '-', '-', 'G',
              'D', 'S', 'C', 'Q', 'D', 'V', 'L', 'S', 'M', 'L', 'G', 'S', 'P',
              'H', 'K', 'V', 'F', 'Y', 'K', 'S', 'E', 'D', 'K', 'M', 'K', 'I',
              'H', 'S', 'P', 'S', 'P', 'H', 'K', 'Q', 'V', 'P', 'S', 'K', 'C',
              'N', 'D', 'Y', 'F', 'F', 'N', 'Y', 'F', 'T', 'L', 'G', 'V', 'D',
              'I', 'L', 'F', 'D', 'A', 'N', 'T', 'H', 'K', 'V', 'K', 'K', 'F',
              'V', 'L', 'H']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 36.91)
        self.assertAlmostEqual(alignment.annotations["E-value"], 3e02)
        self.assertAlmostEqual(alignment.annotations["Score"], 25.55)
        self.assertAlmostEqual(alignment.annotations["Identities"], 10)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.142)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 8.500)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                  EEcCCcHHHHHHHHHHccCCcCcEEEEECCCCCcCCCEEEEeCCCCeEEEecCccceeeEEEEecCCCCEEEECCEEecCCCCCCcHHHHHHhhcCCCCCCccCcccEEEEEeCcEEEEee                                                                                                                                                                                                                                                                                           ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(
            alignment.query.seq[18:139],
            "FTLGMPLAQAVAILQKHCRIIKNVQVLYSEQSPLSHDLILNLTQDGIKLMFDAFNQRLKVIEVCDLTKVKLKYCGVHFNSQAIAPTIEQIDQSFGATHPGVYNSAEQLFHLNFRGLSFSFQ",
        )
        self.assertEqual(
            alignment[1],
            "FTL---------GMPLAQAVAILQKHCRII-------KNVQVLYSEQSPLSHDLILNLTQDGIKLMFDAFNQRLKVIEVCDLTK-VKLKYCGVHFNSQAIAPTIEQIDQSFGATHPGVYNSAEQ----LFHLNF----------RGLSFSFQ",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                  f~LG~sL~~vi~~Lr~~~~~~~~vel~Ys~~~pl~~~Ivi~L~~~GirL~Fd~~~QrL~lIEv~d~~~~~L~Y~~~~~~~~~~~ptf~~I~~~FGPTyPG~yd~~~~~Y~LsYpGisF~Fp                                                                                                                                                                                                                                                                                           ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                              I~~g~~~~~~~~Gmt~~eV~~ilG~P~~~~~~~~~~~~~~~~~W~~~~~~i~v~F~~~kv~~k~~~~~~~~~~~~it~~~~~~I~~Gmt~~eV~~~lG~p~~~~~~~~~~~~~~~~~w~~~~~~~~~~~~i~~~F~             ",
        )
        self.assertEqual(alignment.target.id, "3D4E_A")
        self.assertEqual(
            alignment.target.seq[30:166],
            "IKVTTDQNHFSGGTSIEQLKQWFGDPNKSEQRNAGNITLDSYTWVKDGAVINAQLYKNSTVARSISNFSFSREAKIGKEDYDELKIGESYKKIVEKLGEPDVLSQSMSSDKEEMQTVWSSGIKTKSSSATIELYFE",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "3D4E_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "putative beta-lactamase inhibitor protein; NP_721579.1, putative beta-lactamase inhibitor protein; HET: EDO, MSE; 1.4A {Streptococcus mutans}",
        )
        self.assertEqual(
            alignment[0],
            "IKVTTDQNHFSGGTSIEQLKQWFGDPNKSEQRNAGNITLDSYTWVK------------DGAVINAQLY--KNSTVARSISNFSFSREAKIGKEDYDELKIGESYKKIVEKLGE--PDVLSQSMSSDKEEMQTVWSSGIKTKSSSATIELYFE",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "                              CCCCCCTTTCCSSCCHHHHHHHHCSCSEEEEEEETTEEEEEEEEEETTEEEEEEEETTEEEEEEEESCCCCCCCCBCHHHHHHCCTTCBHHHHHHHHCCCSEEEEEECSSCEEEEEEECSSBCSSCTTCEEEEEEE             ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                              eeecCCccCCCCCCCHHHHHHHHCCCceeEEEeeCCEEEEEEEEEeCCeEEEEEEECCeEEEeeeecCccccCCCCCHHHHHHCCCCCCHHHHHHHHCCCceeeeecCCCceEEEEEEEcCCcCCCCccEEEEEEE             ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 30,  33,  42,  60,  67,  76,  76,  86,  86, 100, 101, 129, 129,
                              138, 142, 148, 158, 166],
                             [ 18,  21,  21,  39,  39,  48,  60,  70,  72,  86,  86, 114, 116,
                              125, 125, 131, 131, 139]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['I', 'K', 'V', 'T', 'T', 'D', 'Q', 'N', 'H', 'F', 'S', 'G', 'G',
              'T', 'S', 'I', 'E', 'Q', 'L', 'K', 'Q', 'W', 'F', 'G', 'D', 'P',
              'N', 'K', 'S', 'E', 'Q', 'R', 'N', 'A', 'G', 'N', 'I', 'T', 'L',
              'D', 'S', 'Y', 'T', 'W', 'V', 'K', '-', '-', '-', '-', '-', '-',
              '-', '-', '-', '-', '-', '-', 'D', 'G', 'A', 'V', 'I', 'N', 'A',
              'Q', 'L', 'Y', '-', '-', 'K', 'N', 'S', 'T', 'V', 'A', 'R', 'S',
              'I', 'S', 'N', 'F', 'S', 'F', 'S', 'R', 'E', 'A', 'K', 'I', 'G',
              'K', 'E', 'D', 'Y', 'D', 'E', 'L', 'K', 'I', 'G', 'E', 'S', 'Y',
              'K', 'K', 'I', 'V', 'E', 'K', 'L', 'G', 'E', '-', '-', 'P', 'D',
              'V', 'L', 'S', 'Q', 'S', 'M', 'S', 'S', 'D', 'K', 'E', 'E', 'M',
              'Q', 'T', 'V', 'W', 'S', 'S', 'G', 'I', 'K', 'T', 'K', 'S', 'S',
              'S', 'A', 'T', 'I', 'E', 'L', 'Y', 'F', 'E'],
             ['F', 'T', 'L', '-', '-', '-', '-', '-', '-', '-', '-', '-', 'G',
              'M', 'P', 'L', 'A', 'Q', 'A', 'V', 'A', 'I', 'L', 'Q', 'K', 'H',
              'C', 'R', 'I', 'I', '-', '-', '-', '-', '-', '-', '-', 'K', 'N',
              'V', 'Q', 'V', 'L', 'Y', 'S', 'E', 'Q', 'S', 'P', 'L', 'S', 'H',
              'D', 'L', 'I', 'L', 'N', 'L', 'T', 'Q', 'D', 'G', 'I', 'K', 'L',
              'M', 'F', 'D', 'A', 'F', 'N', 'Q', 'R', 'L', 'K', 'V', 'I', 'E',
              'V', 'C', 'D', 'L', 'T', 'K', '-', 'V', 'K', 'L', 'K', 'Y', 'C',
              'G', 'V', 'H', 'F', 'N', 'S', 'Q', 'A', 'I', 'A', 'P', 'T', 'I',
              'E', 'Q', 'I', 'D', 'Q', 'S', 'F', 'G', 'A', 'T', 'H', 'P', 'G',
              'V', 'Y', 'N', 'S', 'A', 'E', 'Q', '-', '-', '-', '-', 'L', 'F',
              'H', 'L', 'N', 'F', '-', '-', '-', '-', '-', '-', '-', '-', '-',
              '-', 'R', 'G', 'L', 'S', 'F', 'S', 'F', 'Q']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 36.62)
        self.assertAlmostEqual(alignment.annotations["E-value"], 88)
        self.assertAlmostEqual(alignment.annotations["Score"], 24.38)
        self.assertAlmostEqual(alignment.annotations["Identities"], 23)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.342)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 5.900)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                                                                                       EEEECCEEecCCCCCCcHHHHHHhhcCCCC                                                                                                                                                                                                                                                                                                                 ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(alignment.query.seq[87:117], "KLKYCGVHFNSQAIAPTIEQIDQSFGATHP")
        self.assertEqual(alignment[1], "KLKYCGVHFNSQAIAPTIEQIDQSFGATHP")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                       ~L~Y~~~~~~~~~~~ptf~~I~~~FGPTyP                                                                                                                                                                                                                                                                                                                 ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "    ~F~~~~~~l~DP~P~~s~e~v~~~~a~~yP                            ",
        )
        self.assertEqual(alignment.target.id, "H8H183_DEIGI/2")
        self.assertEqual(alignment.target.seq[4:34], "VFKFDGKVLDDPNPKSTPEQVKTFYAPTYP")
        self.assertEqual(alignment.target.annotations["hmm_name"], "PF14454.6")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "; Prok_Ub ; Prokaryotic Ubiquitin",
        )
        self.assertEqual(alignment[0], "VFKFDGKVLDDPNPKSTPEQVKTFYAPTYP")
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "    EEEECCEEecCCCCCCCHHHHHHHHHHhCH                            ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  4,  34],
                             [ 87, 117]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['V', 'F', 'K', 'F', 'D', 'G', 'K', 'V', 'L', 'D', 'D', 'P', 'N',
              'P', 'K', 'S', 'T', 'P', 'E', 'Q', 'V', 'K', 'T', 'F', 'Y', 'A',
              'P', 'T', 'Y', 'P'],
             ['K', 'L', 'K', 'Y', 'C', 'G', 'V', 'H', 'F', 'N', 'S', 'Q', 'A',
              'I', 'A', 'P', 'T', 'I', 'E', 'Q', 'I', 'D', 'Q', 'S', 'F', 'G',
              'A', 'T', 'H', 'P']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 36.36)
        self.assertAlmostEqual(alignment.annotations["E-value"], 43)
        self.assertAlmostEqual(alignment.annotations["Score"], 30.67)
        self.assertAlmostEqual(alignment.annotations["Identities"], 30)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.464)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 8.500)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                                                                                                                                                                                                                                              EEecCCHHHHHHHhcCCcEE                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(alignment.query.seq[238:258], "VYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(alignment[1], "VYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                                                                                                                                                                              I~~G~T~QDVl~~LG~P~~~                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "           I~~Gmt~~eV~~~lG~p~~~                                                                                                                                      ",
        )
        self.assertEqual(alignment.target.id, "3N4I_B")
        self.assertEqual(alignment.target.seq[11:31], "IQFGMTRQQVLDIAGAENCE")
        self.assertEqual(alignment.target.annotations["hmm_name"], "3N4I_B")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Crystal structure of the SHV-1; beta-lactamase, protein-protein complex, BLIP, SHV-1; 1.56A {Klebsiella pneumoniae} SCOP: d.98.1.1; Related PDB entries: 3E2K_D 1XXM_C 2G2W_B 3E2L_C 3GMU_B 1XXM_D 3C4P_B 1S0W_C 3C7U_D 3C4O_B 1JTG_B 3C7V_B 3C7V_D 3E2L_D 3E2K_C 2G2U_B 3C7U_B 2B5R_C 1JTG_D 2B5R_D 1S0W_D",
        )
        self.assertEqual(alignment[0], "IQFGMTRQQVLDIAGAENCE")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "           SCTTCCHHHHHHHHCGGGEE                                                                                                                                      ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "           CCCCCCHHHHHHHHCCCCee                                                                                                                                      ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 11,  31],
                             [238, 258]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['I', 'Q', 'F', 'G', 'M', 'T', 'R', 'Q', 'Q', 'V', 'L', 'D', 'I',
              'A', 'G', 'A', 'E', 'N', 'C', 'E'],
             ['V', 'Y', 'F', 'G', 'D', 'S', 'C', 'Q', 'D', 'V', 'L', 'S', 'M',
              'L', 'G', 'S', 'P', 'H', 'K', 'V']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 34.69)
        self.assertAlmostEqual(alignment.annotations["E-value"], 44)
        self.assertAlmostEqual(alignment.annotations["Score"], 29.77)
        self.assertAlmostEqual(alignment.annotations["Identities"], 10)
        self.assertAlmostEqual(alignment.annotations["Similarity"], -0.059)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 9.400)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                                                                                                                                                                                                                                              EEecCCHHHHHHHhcCCcEE                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(alignment.query.seq[238:258], "VYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(alignment[1], "VYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                                                                                                                                                                              I~~G~T~QDVl~~LG~P~~~                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                                                                                    i~~Gmt~~~V~~~lG~p~~~                                                  ",
        )
        self.assertEqual(alignment.target.id, "3GMX_B")
        self.assertEqual(alignment.target.seq[84:104], "TQTGMTEAQFWAAVPSDTCS")
        self.assertEqual(alignment.target.annotations["hmm_name"], "3GMX_B")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "BLP; 2-layer alpha/beta sandwich, PROTEIN BINDING; HET: ACT; 1.05A {Streptomyces clavuligerus} SCOP: d.98.1.0; Related PDB entries: 3GMX_A 3GMY_A 3GMY_B",
        )
        self.assertEqual(alignment[0], "TQTGMTEAQFWAAVPSDTCS")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "                                                                                    CCTTCBHHHHHHHSCGGGCE                                                  ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                                                                                    cCCCCCHHHHHHHcCCccce                                                  ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 84, 104],
                             [238, 258]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['T', 'Q', 'T', 'G', 'M', 'T', 'E', 'A', 'Q', 'F', 'W', 'A', 'A',
              'V', 'P', 'S', 'D', 'T', 'C', 'S'],
             ['V', 'Y', 'F', 'G', 'D', 'S', 'C', 'Q', 'D', 'V', 'L', 'S', 'M',
              'L', 'G', 'S', 'P', 'H', 'K', 'V']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 34.59)
        self.assertAlmostEqual(alignment.annotations["E-value"], 71)
        self.assertAlmostEqual(alignment.annotations["Score"], 28.22)
        self.assertAlmostEqual(alignment.annotations["Identities"], 12)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.088)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 8.400)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                                                                                                                                                                                                                                              EEecCCHHHHHHHhcCCcEEEEeccccccccCCCCCcCCCCcCCccchhhHhhccceeecCCCCEEEEE                                                                                                                   ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(
            alignment.query.seq[238:307],
            "VYFGDSCQDVLSMLGSPHKVFYKSEDKMKIHSPSPHKQVPSKCNDYFFNYFTLGVDILFDANTHKVKKF",
        )
        self.assertEqual(
            alignment[1],
            "VYFGDS-CQDVLSMLGSPHKVFYKSEDKMKIHSPSPHKQVPSKCNDYFFNYFTLGVDILFDANTHKVKKF",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                                                                                                                                                                              I~~G~T~QDVl~~LG~P~~~f~K~ddrm~IH~~~~~~~~~~~~~dyF~NYF~lG~DiLfd~~t~~v~Ki                                                                                                                   ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                 ~~lG~s~~~~V~~~~G~P~~~~~~~~~~g~~~~Y~~~~~~f~~~~~~~v~~i                                                                  ",
        )
        self.assertEqual(alignment.target.id, "R4K5N0_CLOPA/8")
        self.assertEqual(
            alignment.target.seq[17:69],
            "FPAKDTNIDSVESKWGKADNSEWVASAKGLYSTYSKHNIVFGSNKGGQIFEV",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "PF14172.6")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "; DUF4309 ; Domain of unknown function (DUF4309)",
        )
        self.assertEqual(
            alignment[0],
            "FPAKDTNIDSVESKWGKADNSEWVASAK-----------------GLYSTYSKHNIVFGSNK-GGQIFEV",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                 CCCCCCcHHHHHHHhCCCCceehhhhcCceeEEecCceEEEEECCCCeEEEE                                                                  ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 17,  23,  24,  45,  45,  62,  62,  69],
                             [238, 244, 244, 265, 282, 299, 300, 307]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['F', 'P', 'A', 'K', 'D', 'T', 'N', 'I', 'D', 'S', 'V', 'E', 'S',
              'K', 'W', 'G', 'K', 'A', 'D', 'N', 'S', 'E', 'W', 'V', 'A', 'S',
              'A', 'K', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',
              '-', '-', '-', '-', '-', '-', 'G', 'L', 'Y', 'S', 'T', 'Y', 'S',
              'K', 'H', 'N', 'I', 'V', 'F', 'G', 'S', 'N', 'K', '-', 'G', 'G',
              'Q', 'I', 'F', 'E', 'V'],
             ['V', 'Y', 'F', 'G', 'D', 'S', '-', 'C', 'Q', 'D', 'V', 'L', 'S',
              'M', 'L', 'G', 'S', 'P', 'H', 'K', 'V', 'F', 'Y', 'K', 'S', 'E',
              'D', 'K', 'M', 'K', 'I', 'H', 'S', 'P', 'S', 'P', 'H', 'K', 'Q',
              'V', 'P', 'S', 'K', 'C', 'N', 'D', 'Y', 'F', 'F', 'N', 'Y', 'F',
              'T', 'L', 'G', 'V', 'D', 'I', 'L', 'F', 'D', 'A', 'N', 'T', 'H',
              'K', 'V', 'K', 'K', 'F']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 34.21)
        self.assertAlmostEqual(alignment.annotations["E-value"], 54)
        self.assertAlmostEqual(alignment.annotations["Score"], 29.27)
        self.assertAlmostEqual(alignment.annotations["Identities"], 30)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.610)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 9.600)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                                                                                                                                                                                                                                              EEecCCHHHHHHHhcCCcEE                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(alignment.query.seq[238:258], "VYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(alignment[1], "VYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                                                                                                                                                                              I~~G~T~QDVl~~LG~P~~~                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                                                                                              i~~G~t~~eV~~~lG~p~~~                                             ",
        )
        self.assertEqual(alignment.target.id, "Q8DTX1_STRMU/4")
        self.assertEqual(alignment.target.seq[94:114], "LKIGESYKKVVEKLGEPDVL")
        self.assertEqual(alignment.target.annotations["hmm_name"], "PF12978.7")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "; DUF3862 ; Domain of Unknown Function with PDB structure (DUF3862)",
        )
        self.assertEqual(alignment[0], "LKIGESYKKVVEKLGEPDVL")
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                                                                                              cCCCCCHHHHHHHHCCCCee                                             ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 94, 114],
                             [238, 258]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['L', 'K', 'I', 'G', 'E', 'S', 'Y', 'K', 'K', 'V', 'V', 'E', 'K',
              'L', 'G', 'E', 'P', 'D', 'V', 'L'],
             ['V', 'Y', 'F', 'G', 'D', 'S', 'C', 'Q', 'D', 'V', 'L', 'S', 'M',
              'L', 'G', 'S', 'P', 'H', 'K', 'V']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 30.93)
        self.assertAlmostEqual(alignment.annotations["E-value"], 61)
        self.assertAlmostEqual(alignment.annotations["Score"], 33.05)
        self.assertAlmostEqual(alignment.annotations["Identities"], 19)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.498)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 9.000)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                                                                                                                                                                                                                                              EEecCCHHHHHHHhcCCcEEE                                                                                                                                                                   ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(alignment.query.seq[238:259], "VYFGDSCQDVLSMLGSPHKVF")
        self.assertEqual(alignment[1], "VYFGDSCQDVLSMLGSPHKVF")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                                                                                                                                                                              I~~G~T~QDVl~~LG~P~~~f                                                                                                                                                                   ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                                          i~iG~s~~~v~~~~G~P~r~~                                                                                                                                                                                                                                                                    ",
        )
        self.assertEqual(alignment.target.id, "4H0A_A")
        self.assertEqual(alignment.target.seq[42:63], "TWVGKDIKVLTSKFGQADRVY")
        self.assertEqual(alignment.target.annotations["hmm_name"], "4H0A_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Uncharacterized protein; CAP protein family, cysteine-rich secretory; HET: EDO, MSE; 1.9A {Staphylococcus aureus subsp. aureus}; Related PDB entries: 4H0A_B",
        )
        self.assertEqual(alignment[0], "TWVGKDIKVLTSKFGQADRVY")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "                                          GGTTSBTHHHHHHHCSCSEEE                                                                                                                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                                          HhcCCcHHHHHHHhCCCcEEE                                                                                                                                                                                                                                                                    ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 42,  63],
                             [238, 259]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['T', 'W', 'V', 'G', 'K', 'D', 'I', 'K', 'V', 'L', 'T', 'S', 'K',
              'F', 'G', 'Q', 'A', 'D', 'R', 'V', 'Y'],
             ['V', 'Y', 'F', 'G', 'D', 'S', 'C', 'Q', 'D', 'V', 'L', 'S', 'M',
              'L', 'G', 'S', 'P', 'H', 'K', 'V', 'F']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 30.07)
        self.assertAlmostEqual(alignment.annotations["E-value"], 3.4e02)
        self.assertAlmostEqual(alignment.annotations["Score"], 25.76)
        self.assertAlmostEqual(alignment.annotations["Identities"], 20)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.191)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 6.200)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "            CCCCccEEcCCcHHHHHHHHHHccCCcCcEEEEECCCCCcCCCEEEEeCCCCeEEEecCcc                                                                                                                                                                                                                                                                                                                                                             ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(
            alignment.query.seq[12:73],
            "GNEQWEFTLGMPLAQAVAILQKHCRIIKNVQVLYSEQSPLSHDLILNLTQDGIKLMFDAFN",
        )
        self.assertEqual(
            alignment[1],
            "GNEQWEFTLGMPLAQAVAILQKHCRII--KNVQVLYSEQ----SPLSHDLILNLTQDGIKLMFDAFN",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "            G~~~~~f~LG~sL~~vi~~Lr~~~~~~~~vel~Ys~~~pl~~~Ivi~L~~~GirL~Fd~~~                                                                                                                                                                                                                                                                                                                                                             ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                                                                                                     G~EHiEfVi~~~~~~~~~~~~~~~~~~~~~~i~~k~s~~~~~~~r~~NP~~~~~~~~~~VKfH~~s        ",
        )
        self.assertEqual(alignment.target.id, "Q7N560_PHOLL/9")
        self.assertEqual(
            alignment.target.seq[101:167],
            "GWEHVELVLPVAPEILSTAAKALLPQPLPAGFSVKESQPKGEQERLPNPTLAITDGEITVKFHPFT",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "PF06185.12")
        self.assertEqual(
            alignment.target.annotations["hmm_description"], "; YecM ; YecM protein"
        )
        self.assertEqual(
            alignment[0],
            "GWEHVELVLPVAPEILSTAAKALLPQPLPAGFSVKESQPKGEQERLPNPT-LAITDGEITVKFHPFT",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                                                                                                     ceeEEEEEcCCCHHHHHHHHHHhCcCcCCCCCeeEecCCCCccccCCCceEEeecCCEEEEEccCC        ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[101, 128, 130, 140, 144, 151, 151, 167],
                             [ 12,  39,  39,  49,  49,  56,  57,  73]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['G', 'W', 'E', 'H', 'V', 'E', 'L', 'V', 'L', 'P', 'V', 'A', 'P',
              'E', 'I', 'L', 'S', 'T', 'A', 'A', 'K', 'A', 'L', 'L', 'P', 'Q',
              'P', 'L', 'P', 'A', 'G', 'F', 'S', 'V', 'K', 'E', 'S', 'Q', 'P',
              'K', 'G', 'E', 'Q', 'E', 'R', 'L', 'P', 'N', 'P', 'T', '-', 'L',
              'A', 'I', 'T', 'D', 'G', 'E', 'I', 'T', 'V', 'K', 'F', 'H', 'P',
              'F', 'T'],
             ['G', 'N', 'E', 'Q', 'W', 'E', 'F', 'T', 'L', 'G', 'M', 'P', 'L',
              'A', 'Q', 'A', 'V', 'A', 'I', 'L', 'Q', 'K', 'H', 'C', 'R', 'I',
              'I', '-', '-', 'K', 'N', 'V', 'Q', 'V', 'L', 'Y', 'S', 'E', 'Q',
              '-', '-', '-', '-', 'S', 'P', 'L', 'S', 'H', 'D', 'L', 'I', 'L',
              'N', 'L', 'T', 'Q', 'D', 'G', 'I', 'K', 'L', 'M', 'F', 'D', 'A',
              'F', 'N']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 29.78)
        self.assertAlmostEqual(alignment.annotations["E-value"], 6.8e02)
        self.assertAlmostEqual(alignment.annotations["Score"], 22.72)
        self.assertAlmostEqual(alignment.annotations["Identities"], 21)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.270)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 1.300)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "      ccccccCCCCccEEcCCcHHHHHHHHHHccCCcCcEEEEECCCCCcCCCEEEEeCCCCeEEEecCccceeeEEEEecCCCCEEEECCEEecCCCCCCcHHHHHHhhcCCCCCCccCcccEEEEEeCcEEEEeeCcccCCCCccCCCcccccccCCCCCCCeEEEEEEEe                                                                                                                                                                                                                                                       ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(
            alignment.query.seq[6:175],
            "VPERSLGNEQWEFTLGMPLAQAVAILQKHCRIIKNVQVLYSEQSPLSHDLILNLTQDGIKLMFDAFNQRLKVIEVCDLTKVKLKYCGVHFNSQAIAPTIEQIDQSFGATHPGVYNSAEQLFHLNFRGLSFSFQLDSWTEAPKYEPNFAHGLASLQIPHGATVKRMYIYS",
        )
        self.assertEqual(
            alignment[1],
            "VPERSLGNEQWEFTLGMPLAQAVAILQKHCRIIKNVQVLYSEQSPLSHDLILNLTQDGIKLMFDAFNQRLKVIEVCDLTKVKLKYCGVHFNSQAIAPTIEQIDQSFGATHPGVYNSAE----------------QLFHLNFRGLSFSFQLDSWTEAPKYEPNFAHGLASLQIPHGATVKRMYIYS",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "      ~Pg~~LG~~~~~f~LG~sL~~vi~~Lr~~~~~~~~vel~Ys~~~pl~~~Ivi~L~~~GirL~Fd~~~QrL~lIEv~d~~~~~L~Y~~~~~~~~~~~ptf~~I~~~FGPTyPG~yd~~~~~Y~LsYpGisF~Fpi~~~~~~~~~~~~~~~~l~~l~~~~~~~~s~i~If~                                                                                                                                                                                                                                                       ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "        ipgegtgI~Ls~GQiLkf~N~PlaEI~~~Ydpsn~S~vSS~VklKGtIHPLFEvPsqisien~Q~nEN~LI~sgfgtslpq~y~~pangyl~~~itnt~tgnigqitl~~g~~tmtfnlqtgenkipviagtqit~~~lts~sailiye         ",
        )
        self.assertEqual(alignment.target.id, "4IL7_A")
        self.assertEqual(
            alignment.target.seq[8:157],
            "IPGEGTGIQLSAGQILKFYNVPIAEIIVEYDPSNVSGVSSNVKLKGTIHPLFEVPSQISIENFQPTENYLIYSGFGTSLPQTYTIPANGYLIISITNTSTGNIGQITLTIGSTTMTFNLQTGENKIPVIAGTQITNMTLTSSSAILIYE",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "4IL7_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Putative uncharacterized protein; partial jelly roll fold, hypothetical; 1.4A {Sulfolobus turreted icosahedral virus}",
        )
        self.assertEqual(
            alignment[0],
            "IPGEGTG-----IQLSAG-----QILKFYNVPIAEIIVEYDPSNVSGVSSNVKLK-GTIHPLFEVPSQ--ISIENFQPTENYLIYSG------------------FGTSLPQTYTIPANGYLIISITNTSTGNIGQITLTIGSTTMTFNLQTGENKIPVIAGTQITNMTL-----TSSSAILIYE",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "        CEEEEEEEEEEECSSCEECCSSEEEEEEEECSSCSEEEEEEEEETTEEEEEEEESEEEEEEECTTCEEEEEEEEEEEEEEEEE                                                                           ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "        CCCCCcceEeccccEEEEecccHheEEEEeCCCCcCCccceeEeeeeeeeceecceeeeeeccCCCCcEEEEecCCCCCCceEEEcCCcEEEEEEEeCCCCceeEEEEEEeeEEEEEEeccCCccccEEeeceeEeeeecCcceEeeee         ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  8,  15,  15,  21,  21,  53,  53,  65,  65,  82,  82,  95, 111,
                              147, 147, 157],
                             [  6,  13,  18,  24,  29,  61,  62,  74,  76,  93, 111, 124, 124,
                              160, 165, 175]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['I', 'P', 'G', 'E', 'G', 'T', 'G', '-', '-', '-', '-', '-', 'I',
              'Q', 'L', 'S', 'A', 'G', '-', '-', '-', '-', '-', 'Q', 'I', 'L',
              'K', 'F', 'Y', 'N', 'V', 'P', 'I', 'A', 'E', 'I', 'I', 'V', 'E',
              'Y', 'D', 'P', 'S', 'N', 'V', 'S', 'G', 'V', 'S', 'S', 'N', 'V',
              'K', 'L', 'K', '-', 'G', 'T', 'I', 'H', 'P', 'L', 'F', 'E', 'V',
              'P', 'S', 'Q', '-', '-', 'I', 'S', 'I', 'E', 'N', 'F', 'Q', 'P',
              'T', 'E', 'N', 'Y', 'L', 'I', 'Y', 'S', 'G', '-', '-', '-', '-',
              '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',
              '-', 'F', 'G', 'T', 'S', 'L', 'P', 'Q', 'T', 'Y', 'T', 'I', 'P',
              'A', 'N', 'G', 'Y', 'L', 'I', 'I', 'S', 'I', 'T', 'N', 'T', 'S',
              'T', 'G', 'N', 'I', 'G', 'Q', 'I', 'T', 'L', 'T', 'I', 'G', 'S',
              'T', 'T', 'M', 'T', 'F', 'N', 'L', 'Q', 'T', 'G', 'E', 'N', 'K',
              'I', 'P', 'V', 'I', 'A', 'G', 'T', 'Q', 'I', 'T', 'N', 'M', 'T',
              'L', '-', '-', '-', '-', '-', 'T', 'S', 'S', 'S', 'A', 'I', 'L',
              'I', 'Y', 'E'],
             ['V', 'P', 'E', 'R', 'S', 'L', 'G', 'N', 'E', 'Q', 'W', 'E', 'F',
              'T', 'L', 'G', 'M', 'P', 'L', 'A', 'Q', 'A', 'V', 'A', 'I', 'L',
              'Q', 'K', 'H', 'C', 'R', 'I', 'I', 'K', 'N', 'V', 'Q', 'V', 'L',
              'Y', 'S', 'E', 'Q', 'S', 'P', 'L', 'S', 'H', 'D', 'L', 'I', 'L',
              'N', 'L', 'T', 'Q', 'D', 'G', 'I', 'K', 'L', 'M', 'F', 'D', 'A',
              'F', 'N', 'Q', 'R', 'L', 'K', 'V', 'I', 'E', 'V', 'C', 'D', 'L',
              'T', 'K', 'V', 'K', 'L', 'K', 'Y', 'C', 'G', 'V', 'H', 'F', 'N',
              'S', 'Q', 'A', 'I', 'A', 'P', 'T', 'I', 'E', 'Q', 'I', 'D', 'Q',
              'S', 'F', 'G', 'A', 'T', 'H', 'P', 'G', 'V', 'Y', 'N', 'S', 'A',
              'E', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',
              '-', '-', '-', '-', 'Q', 'L', 'F', 'H', 'L', 'N', 'F', 'R', 'G',
              'L', 'S', 'F', 'S', 'F', 'Q', 'L', 'D', 'S', 'W', 'T', 'E', 'A',
              'P', 'K', 'Y', 'E', 'P', 'N', 'F', 'A', 'H', 'G', 'L', 'A', 'S',
              'L', 'Q', 'I', 'P', 'H', 'G', 'A', 'T', 'V', 'K', 'R', 'M', 'Y',
              'I', 'Y', 'S']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 26.30)
        self.assertAlmostEqual(alignment.annotations["E-value"], 89)
        self.assertAlmostEqual(alignment.annotations["Score"], 29.01)
        self.assertAlmostEqual(alignment.annotations["Identities"], 25)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.605)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 8.500)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                                                                                                                                                                                                                                              EEecCCHHHHHHHhcCCcEE                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(alignment.query.seq[238:258], "VYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(alignment[1], "VYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                                                                                                                                                                              I~~G~T~QDVl~~LG~P~~~                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                                                                                                                 I~~Gmt~~eV~~~lG~p~~~                                              ",
        )
        self.assertEqual(alignment.target.id, "3D4E_A")
        self.assertEqual(alignment.target.seq[113:133], "LKIGESYKKIVEKLGEPDVL")
        self.assertEqual(alignment.target.annotations["hmm_name"], "3D4E_A")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "putative beta-lactamase inhibitor protein; NP_721579.1, putative beta-lactamase inhibitor protein; HET: EDO, MSE; 1.4A {Streptococcus mutans}",
        )
        self.assertEqual(alignment[0], "LKIGESYKKIVEKLGEPDVL")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "                                                                                                                 CCTTCBHHHHHHHHCCCSEE                                              ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                                                                                                                 CCCCCCHHHHHHHHCCCcee                                              ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[113, 133],
                             [238, 258]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['L', 'K', 'I', 'G', 'E', 'S', 'Y', 'K', 'K', 'I', 'V', 'E', 'K',
              'L', 'G', 'E', 'P', 'D', 'V', 'L'],
             ['V', 'Y', 'F', 'G', 'D', 'S', 'C', 'Q', 'D', 'V', 'L', 'S', 'M',
              'L', 'G', 'S', 'P', 'H', 'K', 'V']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 24.37)
        self.assertAlmostEqual(alignment.annotations["E-value"], 1.1e02)
        self.assertAlmostEqual(alignment.annotations["Score"], 28.12)
        self.assertAlmostEqual(alignment.annotations["Identities"], 25)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.287)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 8.500)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                                                                                                                                                                                                                                              EEecCCHHHHHHHhcCCcEE                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(alignment.query.seq[238:258], "VYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(alignment[1], "VYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                                                                                                                                                                              I~~G~T~QDVl~~LG~P~~~                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                                                                                          i~~Gmt~~eV~~~lG~p~~~                                                       ",
        )
        self.assertEqual(alignment.target.id, "3N4I_B")
        self.assertEqual(alignment.target.seq[90:110], "VTVGMTRAQVLATVGQGSCT")
        self.assertEqual(alignment.target.annotations["hmm_name"], "3N4I_B")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "Crystal structure of the SHV-1; beta-lactamase, protein-protein complex, BLIP, SHV-1; 1.56A {Klebsiella pneumoniae} SCOP: d.98.1.1; Related PDB entries: 3E2K_D 1XXM_C 2G2W_B 3E2L_C 3GMU_B 1XXM_D 3C4P_B 1S0W_C 3C7U_D 3C4O_B 1JTG_B 3C7V_B 3C7V_D 3E2L_D 3E2K_C 2G2U_B 3C7U_B 2B5R_C 1JTG_D 2B5R_D 1S0W_D",
        )
        self.assertEqual(alignment[0], "VTVGMTRAQVLATVGQGSCT")
        self.assertEqual(
            alignment.target.letter_annotations["ss_dssp"],
            "                                                                                          CCTTCBHHHHHHHHCTBSEE                                                       ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                                                                                          ccCCCCHHHHHHHHCCCcee                                                       ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 90, 110],
                             [238, 258]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['V', 'T', 'V', 'G', 'M', 'T', 'R', 'A', 'Q', 'V', 'L', 'A', 'T',
              'V', 'G', 'Q', 'G', 'S', 'C', 'T'],
             ['V', 'Y', 'F', 'G', 'D', 'S', 'C', 'Q', 'D', 'V', 'L', 'S', 'M',
              'L', 'G', 'S', 'P', 'H', 'K', 'V']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 23.84)
        self.assertAlmostEqual(alignment.annotations["E-value"], 89)
        self.assertAlmostEqual(alignment.annotations["Score"], 26.82)
        self.assertAlmostEqual(alignment.annotations["Identities"], 10)
        self.assertAlmostEqual(alignment.annotations["Similarity"], -0.092)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 8.800)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                                                                                                                                                                                                                                              EEecCCHHHHHHHhcCCcEEE                                                                                                                                                                   ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(alignment.query.seq[238:259], "VYFGDSCQDVLSMLGSPHKVF")
        self.assertEqual(alignment[1], "VYFGDSCQDVLSMLGSPHKVF")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                                                                                                                                                                              I~~G~T~QDVl~~LG~P~~~f                                                                                                                                                                   ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                                                                                i~~G~t~~~V~~~lG~p~~~~                    ",
        )
        self.assertEqual(alignment.target.id, "B5GLC0_STRC2/3")
        self.assertEqual(alignment.target.seq[80:101], "TQTGMTEAQFWAAVPSDTCSA")
        self.assertEqual(alignment.target.annotations["hmm_name"], "PF07467.11")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "; BLIP ; Beta-lactamase inhibitor (BLIP)",
        )
        self.assertEqual(alignment[0], "TQTGMTEAQFWAAVPSDTCSA")
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                                                                                cCCCCCHHHHHHhcCCccccc                    ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 80, 101],
                             [238, 259]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['T', 'Q', 'T', 'G', 'M', 'T', 'E', 'A', 'Q', 'F', 'W', 'A', 'A',
              'V', 'P', 'S', 'D', 'T', 'C', 'S', 'A'],
             ['V', 'Y', 'F', 'G', 'D', 'S', 'C', 'Q', 'D', 'V', 'L', 'S', 'M',
              'L', 'G', 'S', 'P', 'H', 'K', 'V', 'F']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 21.56)
        self.assertAlmostEqual(alignment.annotations["E-value"], 3.9e02)
        self.assertAlmostEqual(alignment.annotations["Score"], 22.84)
        self.assertAlmostEqual(alignment.annotations["Identities"], 13)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.078)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 8.800)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                  EEcCCcHHHHHHHHHHccCCcCcEEEEECCCCCcCCCEEEEeCCCCeEEEecCccceeeEEEEecCCCCEEEECCEEecCCCCCCcHHHHHHhhcC                                                                                                                                                                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(
            alignment.query.seq[18:114],
            "FTLGMPLAQAVAILQKHCRIIKNVQVLYSEQSPLSHDLILNLTQDGIKLMFDAFNQRLKVIEVCDLTKVKLKYCGVHFNSQAIAPTIEQIDQSFGA",
        )
        self.assertEqual(
            alignment[1],
            "FTLGMPLAQAVAILQKHCRIIKNVQVLYSEQSPLSHDLILNLTQDGIKLMFDAFNQRLKVIEVCDLTKVKLKYCGVH-FNSQAIAPTIEQIDQSFGA",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                  f~LG~sL~~vi~~Lr~~~~~~~~vel~Ys~~~pl~~~Ivi~L~~~GirL~Fd~~~QrL~lIEv~d~~~~~L~Y~~~~~~~~~~~ptf~~I~~~FGP                                                                                                                                                                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "       I~~G~t~~~V~~~lG~P~~~~~~~~~~~~y~~~~~~~~~~~v~F~~~gkv~~k~~~~~~~~~~~~~t~~~~~~i~~G~t~~~V~~~lG~                         ",
        )
        self.assertEqual(alignment.target.id, "B5GLC0_STRC2/3")
        self.assertEqual(
            alignment.target.seq[7:96],
            "IQFGMDRTLVWQLAGADQSCSDQVERIICYNNPDHYGPQGHFFFNAADKLIHKRQMELFPAPKPTMRLATYNKTQTGMTEAQFWAAVPS",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "PF07467.11")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "; BLIP ; Beta-lactamase inhibitor (BLIP)",
        )
        self.assertEqual(
            alignment[0],
            "IQFGMDRTLVWQLAGADQSCSDQVERIICYNNPDH-------YGPQGHFFFNA-ADKLIHKRQMELFPAPKPTMRLATYNKTQTGMTEAQFWAAVPS",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "       CCCCCCHHHHHHHHCCCceecCCcceEEEEeCCCCCCCEEEEEEcCCCEEEEEEECCCCCCCCCCCCHHHHHhcCCCCCHHHHHHhcCC                         ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  7,  42,  42,  53,  53,  76,  77,  96],
                             [ 18,  53,  60,  71,  72,  95,  95, 114]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['I', 'Q', 'F', 'G', 'M', 'D', 'R', 'T', 'L', 'V', 'W', 'Q', 'L',
              'A', 'G', 'A', 'D', 'Q', 'S', 'C', 'S', 'D', 'Q', 'V', 'E', 'R',
              'I', 'I', 'C', 'Y', 'N', 'N', 'P', 'D', 'H', '-', '-', '-', '-',
              '-', '-', '-', 'Y', 'G', 'P', 'Q', 'G', 'H', 'F', 'F', 'F', 'N',
              'A', '-', 'A', 'D', 'K', 'L', 'I', 'H', 'K', 'R', 'Q', 'M', 'E',
              'L', 'F', 'P', 'A', 'P', 'K', 'P', 'T', 'M', 'R', 'L', 'A', 'T',
              'Y', 'N', 'K', 'T', 'Q', 'T', 'G', 'M', 'T', 'E', 'A', 'Q', 'F',
              'W', 'A', 'A', 'V', 'P', 'S'],
             ['F', 'T', 'L', 'G', 'M', 'P', 'L', 'A', 'Q', 'A', 'V', 'A', 'I',
              'L', 'Q', 'K', 'H', 'C', 'R', 'I', 'I', 'K', 'N', 'V', 'Q', 'V',
              'L', 'Y', 'S', 'E', 'Q', 'S', 'P', 'L', 'S', 'H', 'D', 'L', 'I',
              'L', 'N', 'L', 'T', 'Q', 'D', 'G', 'I', 'K', 'L', 'M', 'F', 'D',
              'A', 'F', 'N', 'Q', 'R', 'L', 'K', 'V', 'I', 'E', 'V', 'C', 'D',
              'L', 'T', 'K', 'V', 'K', 'L', 'K', 'Y', 'C', 'G', 'V', 'H', '-',
              'F', 'N', 'S', 'Q', 'A', 'I', 'A', 'P', 'T', 'I', 'E', 'Q', 'I',
              'D', 'Q', 'S', 'F', 'G', 'A']], dtype='U')
                # fmt: on
            )
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_length(self):
        """Test getting the number of alignments without parsing the file."""
        stream = open(self.path)
        alignments = Align.parse(stream, "hhr")
        stream.close()
        self.assertEqual(len(alignments), 34)


class Align_hhr_hhsearch_q9bsu1_uniclust_w_ss_pfamA_30(unittest.TestCase):
    path = os.path.join("HHsuite", "hhsearch_q9bsu1_uniclust_w_ss_pfamA_30.hhr")

    def test_reading(self):
        alignments = Align.parse(self.path, "hhr")
        self.assertEqual(alignments.metadata["No_of_seqs"], (149, 573))
        self.assertAlmostEqual(alignments.metadata["Neff"], 6.62119)
        self.assertEqual(alignments.metadata["Searched_HMMs"], 16712)
        self.assertEqual(alignments.metadata["Rundate"], "Wed Feb 13 09:26:07 2019")
        self.assertEqual(
            alignments.metadata["Command line"],
            "/home/shah/hh-suite/build/bin/hhsearch -i /home/shah/seq/q9bsu1/hhblits/q9bsu1_uniclust_w_ss.a3m -d /home/shah/db/pfamA_30/pfam -o /home/shah/seq/q9bsu1/hhsearch_q9bsu1_uniclust_w_ss_pfamA_30.hhr -p 20 -Z 250 -loc -z 1 -b 1 -B 250 -ssm 2 -sc 1 -seq 1 -dbstrlen 10000 -norealign -maxres 32000",
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 100.00)
        self.assertAlmostEqual(alignment.annotations["E-value"], 2e-106)
        self.assertAlmostEqual(alignment.annotations["Score"], 822.75)
        self.assertAlmostEqual(alignment.annotations["Identities"], 88)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 1.401)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 6.600)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "          ccCCCCceEEcCCcHHHHHHHHHHccCccccEEEEECCCCCcCCCeEEEecCCCeEEEeccccCceeEEEEecCCCCEEEECCEEecCCCCCCcHhHHHhhcCCCCCCCccCccceEEeeeCeEEEEeeCcccCCCCccCCCcccccccCCCCCCCeEEEEEEEeCCChHhcCCCCCCcchhcCcceeEEEEEeeCCCCCceEEEEEEecCCCCccccccceeEEEEEEEeCCchHHHHHhhcCCcEEEEcCcccCccCCCCCCcCCCCcccccchhHHHhccceeecCCCCeEEEEEecCCCCCCchhhcccCCcEEeeccccccCCCCCcccccCCCCHHHHHHHhCCCCCCCeEeeCCCCCCCCCCCCceEEeccCCEEEEEecCCCEEEEEEe               ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(
            alignment.query.seq[10:407],
            "SLGNEQWEFTLGMPLAQAVAILQKHCRIIKNVQVLYSEQSPLSHDLILNLTQDGIKLMFDAFNQRLKVIEVCDLTKVKLKYCGVHFNSQAIAPTIEQIDQSFGATHPGVYNSAEQLFHLNFRGLSFSFQLDSWTEAPKYEPNFAHGLASLQIPHGATVKRMYIYSGNSLQDTKAPMMPLSCFLGNVYAESVDVLRDGTGPAGLRLRLLAAGCGPGLLADAKMRVFERSVYFGDSCQDVLSMLGSPHKVFYKSEDKMKIHSPSPHKQVPSKCNDYFFNYFTLGVDILFDANTHKVKKFVLHTNYPGHYNFNIYHRCEFKIPLAIKKENADGQTETCTTYSKWDNIQELLGHPVEKPVVLHRSSSPNNTNPFGSTFCFGLQRMIFEVMQNNHIASVTLY",
        )
        self.assertEqual(
            alignment[1],
            "SLGNEQWEFTLGMPLAQAVAILQKHCRIIKNVQVLYSEQSPLSHDLILNLTQDGIKLMFDAFNQRLKVIEVCDLTKVKLKYCGVHFNSQAIAPTIEQIDQSFGATHPGVYNSAEQLFHLNFRGLSFSFQLDSWTEAPKYEPNFAHGLASLQIPHGATVKRMYIYSGNSLQDTKAPMMPLSCFLGNVYAESVDVLRDGTGPAGLRLRLLAAGCGPGLLADAKMRVFERSVYFGDSCQDVLSMLGSPHKVFYKSEDKMKIHSPSPHKQVPSKCNDYFFNYFTLGVDILFDANTHKVKKFVLHTNYPGHYNFNIYHRCEFKIPLAIKKENADGQTE--TCTTYSKWDNIQELLGHPVEKPVVLHRSSSPNNTNPFGSTFCFGLQRMIFEVMQNNHIASVTLY",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "          ~LG~~~~~f~LG~sL~~vi~~Lr~~~~~~~~Vel~Ys~~~Pl~~dIvi~L~~~GirL~Fd~~~QrL~lIEV~d~s~v~L~Y~g~~~~~~~~~pT~~~I~~~FGPTyPG~yd~~~~~Y~LsYpGisF~Fpi~~~~~~~~~~~~~~~~l~~~~~~~~~~~~~m~If~G~s~~~~~~~~~p~~~~~~~~~~~~v~v~~~~~~~~gl~~~~~~~~~~~~~~~~~~~~~~~~~I~fG~T~QDVl~~LG~P~~vf~K~ddrm~IH~~~~~~~~~~~~~dyF~NYF~lGlDiLfd~~th~v~KiILHtN~PGh~~F~~Y~RC~f~i~~~~~~~~~~~~~~~it~~~k~~~i~~~l~~~~~~pvvlnR~~s~~~~npfg~T~lyG~~~~ifEVm~n~~IasVTlf               ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "~lG~f~LG~sL~~vl~~Lr~~~~~~~~v~i~Ys~~~P~~~~Ivi~Lp~~GirL~Fd~~~QrL~lIEv~d~~~~~L~Y~~~~~~~~~~~~tf~~Iy~~FGPTyPG~~d~~~~~y~LsYpGlsF~Fp~~~~~~~~~~~~~~~~~l~~~~~~~~~~~~~~~If~G~~~~~~~~p~~p~~~~~~~~~~~~v~v~~~~~~~~g~~l~f~~~~~~~~~~~~~~~~~~~~~I~~G~TtQDVl~~LG~P~~vf~K~ddrm~IH~~~~~~~~~~~~~~yF~NYF~lGlDiLfd~~th~v~KiILHtN~Pg~~~F~~Y~RC~w~i~~~~~~~~~~~~~~~~~i~~~~~~~~i~~~l~~~~~~p~vLnR~~~~~~~~~~g~t~lyG~~g~IfEV~~ng~IasVTlf",
        )
        self.assertEqual(alignment.target.id, "Q1LV04_DANRE/1")
        self.assertEqual(
            alignment.target.seq[0:395],
            "EQWEFALGMPLAQAISILQKHCRIIKNVQVLYSEQMPLSHDLILNLTQDGIKLLFDACNQRLKVIEVYDLTKVKLKYCGVHFNSQAIAPTIEQIDQSFGATHPGVYNAAEQLFHLNFRGLSFSFQLDSWSEAPKYEPNFAHGLASLQIPHGATVKRMYIYSGNNLQETKAPAMPLACFLGNVYAECVEVLRDGAGPLGLKLRLLTAGCGPGVLADTKVRAVERSIYFGDSCQDVLSALGSPHKVFYKSEDKMKIHSPSPHKQVPSKCNDYFFNYYILGVDILFDSTTHLVKKFVLHTNFPGHYNFNIYHRCDFKIPLIIKKDGADAHSEDCILTTYSKWDQIQELLGHPMEKPVVLHRSSSANNTNPFGSTFCFGLQRMIFEVMQNNHIASVTLY",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "PF03676.13")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "; UPF0183 ; Uncharacterised protein family (UPF0183)",
        )
        self.assertEqual(
            alignment[0],
            "EQWE----FALGMPLAQAISILQKHCRIIKNVQVLYSEQMPLSHDLILNLTQDGIKLLFDACNQRLKVIEVYDLTKVKLKYCGVHFNSQAIAPTIEQIDQSFGATHPGVYNAAEQLFHLNFRGLSFSFQLDSWSEAPKYEPNFAHGLASLQIPHGATVKRMYIYSGNNLQETKAPAMPLACFLGNVYAECVEVLRDGAGPLGLKLRLLTAGCGPGVLADTKVRAVERSIYFGDSCQDVLSALGSPHKVFYKSEDKMKIHSPSPHKQVPSKCNDYFFNYYILGVDILFDSTTHLVKKFVLHTNFPGHYNFNIYHRCDFKIPLIIKKDGADAHSEDCILTTYSKWDQIQELLGHPMEKPVVLHRSSSANNTNPFGSTFCFGLQRMIFEVMQNNHIASVTLY",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "CCccEEeCCcHHHHHHHHHhcccccccEEEEECCCCccCCcEEEEeCCCcEEEEEcCCCcEEEEEEEEcCCccEEEECCEEecCCCCCCCHHHHHHHhcCCCCCccccccCEEEEeeCCEEEEeeCCccccCCccCCCcccccccccCCCCCeeEEEEEeeCCCccccCCCCCCchhcCCceeeEEEEEEeCCCCCcceEEEEEeCCCCCCcccccccccceEEEEecCCHHHHHHHhCCCceEEEcCcccccccCCCCCCCCCccccccchhhHHhceeeeeeCCCceeEEEEEecCCCCccccCccccCcEEEEecCCcCCCccCcccceecCCCcHHHHHHHhcCCCCCCEEEecCCCCCCCCCCccEEEEeeCCEEEEEecCCcEEEEEeC",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  0,   4,   4, 329, 331, 395],
                             [ 10,  14,  18, 343, 343, 407]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['E', 'Q', 'W', 'E', '-', '-', '-', '-', 'F', 'A', 'L', 'G', 'M',
              'P', 'L', 'A', 'Q', 'A', 'I', 'S', 'I', 'L', 'Q', 'K', 'H', 'C',
              'R', 'I', 'I', 'K', 'N', 'V', 'Q', 'V', 'L', 'Y', 'S', 'E', 'Q',
              'M', 'P', 'L', 'S', 'H', 'D', 'L', 'I', 'L', 'N', 'L', 'T', 'Q',
              'D', 'G', 'I', 'K', 'L', 'L', 'F', 'D', 'A', 'C', 'N', 'Q', 'R',
              'L', 'K', 'V', 'I', 'E', 'V', 'Y', 'D', 'L', 'T', 'K', 'V', 'K',
              'L', 'K', 'Y', 'C', 'G', 'V', 'H', 'F', 'N', 'S', 'Q', 'A', 'I',
              'A', 'P', 'T', 'I', 'E', 'Q', 'I', 'D', 'Q', 'S', 'F', 'G', 'A',
              'T', 'H', 'P', 'G', 'V', 'Y', 'N', 'A', 'A', 'E', 'Q', 'L', 'F',
              'H', 'L', 'N', 'F', 'R', 'G', 'L', 'S', 'F', 'S', 'F', 'Q', 'L',
              'D', 'S', 'W', 'S', 'E', 'A', 'P', 'K', 'Y', 'E', 'P', 'N', 'F',
              'A', 'H', 'G', 'L', 'A', 'S', 'L', 'Q', 'I', 'P', 'H', 'G', 'A',
              'T', 'V', 'K', 'R', 'M', 'Y', 'I', 'Y', 'S', 'G', 'N', 'N', 'L',
              'Q', 'E', 'T', 'K', 'A', 'P', 'A', 'M', 'P', 'L', 'A', 'C', 'F',
              'L', 'G', 'N', 'V', 'Y', 'A', 'E', 'C', 'V', 'E', 'V', 'L', 'R',
              'D', 'G', 'A', 'G', 'P', 'L', 'G', 'L', 'K', 'L', 'R', 'L', 'L',
              'T', 'A', 'G', 'C', 'G', 'P', 'G', 'V', 'L', 'A', 'D', 'T', 'K',
              'V', 'R', 'A', 'V', 'E', 'R', 'S', 'I', 'Y', 'F', 'G', 'D', 'S',
              'C', 'Q', 'D', 'V', 'L', 'S', 'A', 'L', 'G', 'S', 'P', 'H', 'K',
              'V', 'F', 'Y', 'K', 'S', 'E', 'D', 'K', 'M', 'K', 'I', 'H', 'S',
              'P', 'S', 'P', 'H', 'K', 'Q', 'V', 'P', 'S', 'K', 'C', 'N', 'D',
              'Y', 'F', 'F', 'N', 'Y', 'Y', 'I', 'L', 'G', 'V', 'D', 'I', 'L',
              'F', 'D', 'S', 'T', 'T', 'H', 'L', 'V', 'K', 'K', 'F', 'V', 'L',
              'H', 'T', 'N', 'F', 'P', 'G', 'H', 'Y', 'N', 'F', 'N', 'I', 'Y',
              'H', 'R', 'C', 'D', 'F', 'K', 'I', 'P', 'L', 'I', 'I', 'K', 'K',
              'D', 'G', 'A', 'D', 'A', 'H', 'S', 'E', 'D', 'C', 'I', 'L', 'T',
              'T', 'Y', 'S', 'K', 'W', 'D', 'Q', 'I', 'Q', 'E', 'L', 'L', 'G',
              'H', 'P', 'M', 'E', 'K', 'P', 'V', 'V', 'L', 'H', 'R', 'S', 'S',
              'S', 'A', 'N', 'N', 'T', 'N', 'P', 'F', 'G', 'S', 'T', 'F', 'C',
              'F', 'G', 'L', 'Q', 'R', 'M', 'I', 'F', 'E', 'V', 'M', 'Q', 'N',
              'N', 'H', 'I', 'A', 'S', 'V', 'T', 'L', 'Y'],
             ['S', 'L', 'G', 'N', 'E', 'Q', 'W', 'E', 'F', 'T', 'L', 'G', 'M',
              'P', 'L', 'A', 'Q', 'A', 'V', 'A', 'I', 'L', 'Q', 'K', 'H', 'C',
              'R', 'I', 'I', 'K', 'N', 'V', 'Q', 'V', 'L', 'Y', 'S', 'E', 'Q',
              'S', 'P', 'L', 'S', 'H', 'D', 'L', 'I', 'L', 'N', 'L', 'T', 'Q',
              'D', 'G', 'I', 'K', 'L', 'M', 'F', 'D', 'A', 'F', 'N', 'Q', 'R',
              'L', 'K', 'V', 'I', 'E', 'V', 'C', 'D', 'L', 'T', 'K', 'V', 'K',
              'L', 'K', 'Y', 'C', 'G', 'V', 'H', 'F', 'N', 'S', 'Q', 'A', 'I',
              'A', 'P', 'T', 'I', 'E', 'Q', 'I', 'D', 'Q', 'S', 'F', 'G', 'A',
              'T', 'H', 'P', 'G', 'V', 'Y', 'N', 'S', 'A', 'E', 'Q', 'L', 'F',
              'H', 'L', 'N', 'F', 'R', 'G', 'L', 'S', 'F', 'S', 'F', 'Q', 'L',
              'D', 'S', 'W', 'T', 'E', 'A', 'P', 'K', 'Y', 'E', 'P', 'N', 'F',
              'A', 'H', 'G', 'L', 'A', 'S', 'L', 'Q', 'I', 'P', 'H', 'G', 'A',
              'T', 'V', 'K', 'R', 'M', 'Y', 'I', 'Y', 'S', 'G', 'N', 'S', 'L',
              'Q', 'D', 'T', 'K', 'A', 'P', 'M', 'M', 'P', 'L', 'S', 'C', 'F',
              'L', 'G', 'N', 'V', 'Y', 'A', 'E', 'S', 'V', 'D', 'V', 'L', 'R',
              'D', 'G', 'T', 'G', 'P', 'A', 'G', 'L', 'R', 'L', 'R', 'L', 'L',
              'A', 'A', 'G', 'C', 'G', 'P', 'G', 'L', 'L', 'A', 'D', 'A', 'K',
              'M', 'R', 'V', 'F', 'E', 'R', 'S', 'V', 'Y', 'F', 'G', 'D', 'S',
              'C', 'Q', 'D', 'V', 'L', 'S', 'M', 'L', 'G', 'S', 'P', 'H', 'K',
              'V', 'F', 'Y', 'K', 'S', 'E', 'D', 'K', 'M', 'K', 'I', 'H', 'S',
              'P', 'S', 'P', 'H', 'K', 'Q', 'V', 'P', 'S', 'K', 'C', 'N', 'D',
              'Y', 'F', 'F', 'N', 'Y', 'F', 'T', 'L', 'G', 'V', 'D', 'I', 'L',
              'F', 'D', 'A', 'N', 'T', 'H', 'K', 'V', 'K', 'K', 'F', 'V', 'L',
              'H', 'T', 'N', 'Y', 'P', 'G', 'H', 'Y', 'N', 'F', 'N', 'I', 'Y',
              'H', 'R', 'C', 'E', 'F', 'K', 'I', 'P', 'L', 'A', 'I', 'K', 'K',
              'E', 'N', 'A', 'D', 'G', 'Q', 'T', 'E', '-', '-', 'T', 'C', 'T',
              'T', 'Y', 'S', 'K', 'W', 'D', 'N', 'I', 'Q', 'E', 'L', 'L', 'G',
              'H', 'P', 'V', 'E', 'K', 'P', 'V', 'V', 'L', 'H', 'R', 'S', 'S',
              'S', 'P', 'N', 'N', 'T', 'N', 'P', 'F', 'G', 'S', 'T', 'F', 'C',
              'F', 'G', 'L', 'Q', 'R', 'M', 'I', 'F', 'E', 'V', 'M', 'Q', 'N',
              'N', 'H', 'I', 'A', 'S', 'V', 'T', 'L', 'Y']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 71.64)
        self.assertAlmostEqual(alignment.annotations["E-value"], 0.97)
        self.assertAlmostEqual(alignment.annotations["Score"], 31.45)
        self.assertAlmostEqual(alignment.annotations["Identities"], 30)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.363)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 9.700)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                                                                                                                                                                                                                                              EEeCCchHHHHHhhcCCcEE                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(alignment.query.seq[238:258], "VYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(alignment[1], "VYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                                                                                                                                                                              I~fG~T~QDVl~~LG~P~~v                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "           i~~G~Tk~eV~~~lG~P~~~                                          ",
        )
        self.assertEqual(alignment.target.id, "A6W2D6_MARMS/3")
        self.assertEqual(alignment.target.seq[11:31], "LQIGMSESQVTYLLGNPMLR")
        self.assertEqual(alignment.target.annotations["hmm_name"], "PF04355.12")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "; SmpA_OmlA ; SmpA / OmlA family",
        )
        self.assertEqual(alignment[0], "LQIGMSESQVTYLLGNPMLR")
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "           CCCCCCHHHHHHHHCCCCcc                                          ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 11,  31],
                             [238, 258]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['L', 'Q', 'I', 'G', 'M', 'S', 'E', 'S', 'Q', 'V', 'T', 'Y', 'L',
              'L', 'G', 'N', 'P', 'M', 'L', 'R'],
             ['V', 'Y', 'F', 'G', 'D', 'S', 'C', 'Q', 'D', 'V', 'L', 'S', 'M',
              'L', 'G', 'S', 'P', 'H', 'K', 'V']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 53.47)
        self.assertAlmostEqual(alignment.annotations["E-value"], 4.1)
        self.assertAlmostEqual(alignment.annotations["Score"], 30.99)
        self.assertAlmostEqual(alignment.annotations["Identities"], 15)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.388)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 8.100)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                                                                                                                                                                                                                                              EEeCCchHHHHHhhcCCcEE                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(alignment.query.seq[238:258], "VYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(alignment[1], "VYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                                                                                                                                                                              I~fG~T~QDVl~~LG~P~~v                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                                 l~~GMTk~qV~~iLG~P~~~                                                   ",
        )
        self.assertEqual(alignment.target.id, "A3QGB8_SHELP/1")
        self.assertEqual(alignment.target.seq[33:53], "LSLGMTRDQVMTLMGTADFN")
        self.assertEqual(alignment.target.annotations["hmm_name"], "PF11399.7")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "; DUF3192 ; Protein of unknown function (DUF3192)",
        )
        self.assertEqual(alignment[0], "LSLGMTRDQVMTLMGTADFN")
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                                 CCCCCCHHHHHHHhCCCcee                                                   ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 33,  53],
                             [238, 258]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['L', 'S', 'L', 'G', 'M', 'T', 'R', 'D', 'Q', 'V', 'M', 'T', 'L',
              'M', 'G', 'T', 'A', 'D', 'F', 'N'],
             ['V', 'Y', 'F', 'G', 'D', 'S', 'C', 'Q', 'D', 'V', 'L', 'S', 'M',
              'L', 'G', 'S', 'P', 'H', 'K', 'V']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 52.40)
        self.assertAlmostEqual(alignment.annotations["E-value"], 7.5)
        self.assertAlmostEqual(alignment.annotations["Score"], 31.10)
        self.assertAlmostEqual(alignment.annotations["Identities"], 12)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.062)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 7.600)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                                                                                                                                                                                                                                              EEeCCchHHHHHhhcCCcEEEEcCcccCccCCCCCCcCCCCcccccchhHHHhccceeecCCCCeEEEEEecCCCCC                                                                                                           ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(
            alignment.query.seq[238:315],
            "VYFGDSCQDVLSMLGSPHKVFYKSEDKMKIHSPSPHKQVPSKCNDYFFNYFTLGVDILFDANTHKVKKFVLHTNYPG",
        )
        self.assertEqual(
            alignment[1],
            "VYFGDS-CQDVLSMLGSPHKVFYKSEDKMKIHSPSPHKQVPSKCNDYFFNYFTLGVDILFDANTHKVKKFVLHTNYPG",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                                                                                                                                                                              I~fG~T~QDVl~~LG~P~~vf~K~ddrm~IH~~~~~~~~~~~~~dyF~NYF~lGlDiLfd~~th~v~KiILHtN~PG                                                                                                           ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                 ~~ig~s~~~eV~~~~G~P~~~~~~~~~~g~y~~y~~~~~~fg~~~~~~i~~ir~~~~~l~~                                                         ",
        )
        self.assertEqual(alignment.target.id, "R4K5N0_CLOPA/8")
        self.assertEqual(
            alignment.target.seq[17:78],
            "FPAKDTNIDSVESKWGKADNSEWVASAKGLYSTYSKHNIVFGSNKGGQIFEVRSLDKQLGN",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "PF14172.5")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "; DUF4309 ; Domain of unknown function (DUF4309)",
        )
        self.assertEqual(
            alignment[0],
            "FPAKDTNIDSVESKWGKADNSEWVASAK-----------------GLYSTYSKHNIVFGSNKGGQIFEVRSLDKQLGN",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                 CCCCCCcHHHHHHHhCCCCceehhhhcCCeeEEEeCCeEEEEECCCCcEEEEEEeCccccc                                                         ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 17,  23,  24,  45,  45,  78],
                             [238, 244, 244, 265, 282, 315]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['F', 'P', 'A', 'K', 'D', 'T', 'N', 'I', 'D', 'S', 'V', 'E', 'S',
              'K', 'W', 'G', 'K', 'A', 'D', 'N', 'S', 'E', 'W', 'V', 'A', 'S',
              'A', 'K', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',
              '-', '-', '-', '-', '-', '-', 'G', 'L', 'Y', 'S', 'T', 'Y', 'S',
              'K', 'H', 'N', 'I', 'V', 'F', 'G', 'S', 'N', 'K', 'G', 'G', 'Q',
              'I', 'F', 'E', 'V', 'R', 'S', 'L', 'D', 'K', 'Q', 'L', 'G', 'N'],
             ['V', 'Y', 'F', 'G', 'D', 'S', '-', 'C', 'Q', 'D', 'V', 'L', 'S',
              'M', 'L', 'G', 'S', 'P', 'H', 'K', 'V', 'F', 'Y', 'K', 'S', 'E',
              'D', 'K', 'M', 'K', 'I', 'H', 'S', 'P', 'S', 'P', 'H', 'K', 'Q',
              'V', 'P', 'S', 'K', 'C', 'N', 'D', 'Y', 'F', 'F', 'N', 'Y', 'F',
              'T', 'L', 'G', 'V', 'D', 'I', 'L', 'F', 'D', 'A', 'N', 'T', 'H',
              'K', 'V', 'K', 'K', 'F', 'V', 'L', 'H', 'T', 'N', 'Y', 'P', 'G']],
            dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 48.79)
        self.assertAlmostEqual(alignment.annotations["E-value"], 5.3)
        self.assertAlmostEqual(alignment.annotations["Score"], 29.67)
        self.assertAlmostEqual(alignment.annotations["Identities"], 27)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.282)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 8.000)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                                                                                                                                                                                                                                    ccceeEEEEEEEeCCchHHHHHhhcCCcEE                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(alignment.query.seq[228:258], "DAKMRVFERSVYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(alignment[1], "DAKMRVFERSVYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                                                                                                                                                                    ~~~~~~~~~~I~fG~T~QDVl~~LG~P~~v                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "a~s~RCg~~lV~~Gds~~~Vl~~cGeP~~~                                                               ",
        )
        self.assertEqual(alignment.target.id, "H8MKT6_CORCM/2")
        self.assertEqual(alignment.target.seq[0:30], "ASALRCDNKIVSEGASQVDALAKCGQPVTK")
        self.assertEqual(alignment.target.annotations["hmm_name"], "PF11006.7")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "; DUF2845 ; Protein of unknown function (DUF2845)",
        )
        self.assertEqual(alignment[0], "ASALRCDNKIVSEGASQVDALAKCGQPVTK")
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "CCCcccCCcccCCCCCHHHHHHHhCCCeee                                                               ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  0,  30],
                             [228, 258]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['A', 'S', 'A', 'L', 'R', 'C', 'D', 'N', 'K', 'I', 'V', 'S', 'E',
              'G', 'A', 'S', 'Q', 'V', 'D', 'A', 'L', 'A', 'K', 'C', 'G', 'Q',
              'P', 'V', 'T', 'K'],
             ['D', 'A', 'K', 'M', 'R', 'V', 'F', 'E', 'R', 'S', 'V', 'Y', 'F',
              'G', 'D', 'S', 'C', 'Q', 'D', 'V', 'L', 'S', 'M', 'L', 'G', 'S',
              'P', 'H', 'K', 'V']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 41.13)
        self.assertAlmostEqual(alignment.annotations["E-value"], 9.7)
        self.assertAlmostEqual(alignment.annotations["Score"], 29.55)
        self.assertAlmostEqual(alignment.annotations["Identities"], 26)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.562)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 9.700)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                                                                                                                                                                                                                                                eCCchHHHHHhhcCCcEEE                                                                                                                                                                   ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(alignment.query.seq[240:259], "FGDSCQDVLSMLGSPHKVF")
        self.assertEqual(alignment[1], "FGDSCQDVLSMLGSPHKVF")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                                                                                                                                                                                fG~T~QDVl~~LG~P~~vf                                                                                                                                                                   ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "iG~~~~~v~~~lG~P~~~~                                                                                                                         ",
        )
        self.assertEqual(alignment.target.id, "U5LAW6_9BACI/1")
        self.assertEqual(alignment.target.seq[0:19], "IGKNASDLQVLLGDPERKD")
        self.assertEqual(alignment.target.annotations["hmm_name"], "PF14504.5")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "; CAP_assoc_N ; CAP-associated N-terminal",
        )
        self.assertEqual(alignment[0], "IGKNASDLQVLLGDPERKD")
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "CCCCHHHHHHHHCCCceeC                                                                                                                         ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  0,  19],
                             [240, 259]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['I', 'G', 'K', 'N', 'A', 'S', 'D', 'L', 'Q', 'V', 'L', 'L', 'G',
              'D', 'P', 'E', 'R', 'K', 'D'],
             ['F', 'G', 'D', 'S', 'C', 'Q', 'D', 'V', 'L', 'S', 'M', 'L', 'G',
              'S', 'P', 'H', 'K', 'V', 'F']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 41.12)
        self.assertAlmostEqual(alignment.annotations["E-value"], 17)
        self.assertAlmostEqual(alignment.annotations["Score"], 28.06)
        self.assertAlmostEqual(alignment.annotations["Identities"], 7)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.077)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 9.700)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                                                                                                                                                                                                                                              EEeCCchHHHHHhhcCCcEEEEcCcccCccCCCCCCcCCCCcccccchhHHHhccceeecCCCCeEEEEEe                                                                                                                 ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(
            alignment.query.seq[238:309],
            "VYFGDSCQDVLSMLGSPHKVFYKSEDKMKIHSPSPHKQVPSKCNDYFFNYFTLGVDILFDANTHKVKKFVL",
        )
        self.assertEqual(
            alignment[1],
            "VYFGDSCQDVLSMLGSPHKV---------FYKSEDKMKIHSPSPHKQVPSKCNDYFFNYFTLGVDILFDANTHKVKKFVL",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                                                                                                                                                                              I~fG~T~QDVl~~LG~P~~vf~K~ddrm~IH~~~~~~~~~~~~~dyF~NYF~lGlDiLfd~~th~v~KiIL                                                                                                                 ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                                                             i~iG~s~~~v~~~~G~p~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~d~~~g~v~~i~~            ",
        )
        self.assertEqual(alignment.target.id, "U5LAW6_9BACI/1")
        self.assertEqual(
            alignment.target.seq[61:128],
            "FHIGQPVSEIYSSVFIDTNINFQYKGSSYRFELSEDDLNTRPLIKAGNIYAQLYIDRFTGELSSIRY",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "PF14504.5")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "; CAP_assoc_N ; CAP-associated N-terminal",
        )
        self.assertEqual(
            alignment[0],
            "FHIGQPVSEIYSSVFIDTNINFQYKGSSYRFELSE-------------DDLNTRPLIKAGNIYAQLYIDRFTGELSSIRY",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                                                             ccCCCcHHHHHHhhCCCceEEEEEcCeEEEEEECcchhceeeeEEECCEEEEEEEeCCCCEEEEEEE            ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 61,  81,  90,  96,  96, 128],
                             [238, 258, 258, 264, 277, 309]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['F', 'H', 'I', 'G', 'Q', 'P', 'V', 'S', 'E', 'I', 'Y', 'S', 'S',
              'V', 'F', 'I', 'D', 'T', 'N', 'I', 'N', 'F', 'Q', 'Y', 'K', 'G',
              'S', 'S', 'Y', 'R', 'F', 'E', 'L', 'S', 'E', '-', '-', '-', '-',
              '-', '-', '-', '-', '-', '-', '-', '-', '-', 'D', 'D', 'L', 'N',
              'T', 'R', 'P', 'L', 'I', 'K', 'A', 'G', 'N', 'I', 'Y', 'A', 'Q',
              'L', 'Y', 'I', 'D', 'R', 'F', 'T', 'G', 'E', 'L', 'S', 'S', 'I',
              'R', 'Y'],
             ['V', 'Y', 'F', 'G', 'D', 'S', 'C', 'Q', 'D', 'V', 'L', 'S', 'M',
              'L', 'G', 'S', 'P', 'H', 'K', 'V', '-', '-', '-', '-', '-', '-',
              '-', '-', '-', 'F', 'Y', 'K', 'S', 'E', 'D', 'K', 'M', 'K', 'I',
              'H', 'S', 'P', 'S', 'P', 'H', 'K', 'Q', 'V', 'P', 'S', 'K', 'C',
              'N', 'D', 'Y', 'F', 'F', 'N', 'Y', 'F', 'T', 'L', 'G', 'V', 'D',
              'I', 'L', 'F', 'D', 'A', 'N', 'T', 'H', 'K', 'V', 'K', 'K', 'F',
              'V', 'L']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 36.58)
        self.assertAlmostEqual(alignment.annotations["E-value"], 21)
        self.assertAlmostEqual(alignment.annotations["Score"], 24.91)
        self.assertAlmostEqual(alignment.annotations["Identities"], 23)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.342)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 4.700)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                                                                                       EEEECCEEecCCCCCCcHhHHHhhcCCCCC                                                                                                                                                                                                                                                                                                                 ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(alignment.query.seq[87:117], "KLKYCGVHFNSQAIAPTIEQIDQSFGATHP")
        self.assertEqual(alignment[1], "KLKYCGVHFNSQAIAPTIEQIDQSFGATHP")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                       ~L~Y~g~~~~~~~~~pT~~~I~~~FGPTyP                                                                                                                                                                                                                                                                                                                 ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "    ~F~~~g~~L~DP~P~~spe~V~~~ya~~yP                            ",
        )
        self.assertEqual(alignment.target.id, "H8H183_DEIGI/2")
        self.assertEqual(alignment.target.seq[4:34], "VFKFDGKVLDDPNPKSTPEQVKTFYAPTYP")
        self.assertEqual(alignment.target.annotations["hmm_name"], "PF14454.5")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "; Prok_Ub ; Prokaryotic Ubiquitin",
        )
        self.assertEqual(alignment[0], "VFKFDGKVLDDPNPKSTPEQVKTFYAPTYP")
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "    EEEECCEEccCCCCCCCHHHHHHHHHhhCH                            ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[  4,  34],
                             [ 87, 117]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['V', 'F', 'K', 'F', 'D', 'G', 'K', 'V', 'L', 'D', 'D', 'P', 'N',
              'P', 'K', 'S', 'T', 'P', 'E', 'Q', 'V', 'K', 'T', 'F', 'Y', 'A',
              'P', 'T', 'Y', 'P'],
             ['K', 'L', 'K', 'Y', 'C', 'G', 'V', 'H', 'F', 'N', 'S', 'Q', 'A',
              'I', 'A', 'P', 'T', 'I', 'E', 'Q', 'I', 'D', 'Q', 'S', 'F', 'G',
              'A', 'T', 'H', 'P']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 35.95)
        self.assertAlmostEqual(alignment.annotations["E-value"], 1.3e02)
        self.assertAlmostEqual(alignment.annotations["Score"], 23.17)
        self.assertAlmostEqual(alignment.annotations["Identities"], 12)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.155)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 7.600)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "        ccccCCCCceEEcCCcHHHHHHHHHHccCccccEEEEECCCCCcCCCeEEEecCCCeEEEeccccCceeEEEEecCCCCEEEECCEEecCCCCCCcHhHHHhhcCCCCCCCccCccceEEeeeCeEEEEeeCcccCCCCccCCCcccccccCCCCCCCeEEEEEEE                                                                                                                                                                                                                                                        ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(
            alignment.query.seq[8:174],
            "ERSLGNEQWEFTLGMPLAQAVAILQKHCRIIKNVQVLYSEQSPLSHDLILNLTQDGIKLMFDAFNQRLKVIEVCDLTKVKLKYCGVHFNSQAIAPTIEQIDQSFGATHPGVYNSAEQLFHLNFRGLSFSFQLDSWTEAPKYEPNFAHGLASLQIPHGATVKRMYIY",
        )
        self.assertEqual(
            alignment[1],
            "ERSLGNEQWEFTLGMP-LAQAVAILQKHCRIIKNVQVLYSEQSPLSHD-------LILNLTQDGIKLMFDAFNQRLKVIEVCDLTKVKLKYCGVHFNSQAIAPTIEQIDQSFGATHPGVYNSAEQLFHLNFR-G----LSFSFQLDSWTEAPKYEPNFAHGLASLQIPHGATVKRMYIY",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "        g~~LG~~~~~f~LG~sL~~vi~~Lr~~~~~~~~Vel~Ys~~~Pl~~dIvi~L~~~GirL~Fd~~~QrL~lIEV~d~s~v~L~Y~g~~~~~~~~~pT~~~I~~~FGPTyPG~yd~~~~~Y~LsYpGisF~Fpi~~~~~~~~~~~~~~~~l~~~~~~~~~~~~~m~If                                                                                                                                                                                                                                                        ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "          G~~~~~~~~ig~s~~~eV~~~~G~P~~~~~~~~~~g~y~~y~~~~~~fg~~~~~~i~~ir~~~~~l~~lt~~~ikk~LG~P~~~~~~~g~~~~~Y~~g~~y~L~f~~~~~~~~~~~~~v~hi~l~",
        )
        self.assertEqual(alignment.target.id, "R4K5N0_CLOPA/8")
        self.assertEqual(
            alignment.target.seq[10:135],
            "GKVFNSDFPAKDTNIDSVESKWGKADNSEWVASAKGLYSTYSKHNIVFGSNKGGQIFEVRSLDKQLGNIYLSMVKDKLGTPQHDVKVNGEEIIGYKMGNDFKILFVFPEPTNQHANPIMSHYSVL",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "PF14172.5")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "; DUF4309 ; Domain of unknown function (DUF4309)",
        )
        self.assertEqual(
            alignment[0],
            "GKVFNS---DFPAKDTNIDSVESKW----------------GKADNSEWVASAKGLYSTYSKHNIVFGSNKGGQ---IFEVRSL------------DKQLGNIYLSMVKDKLGT--PQHDVKVNGEEIIGYKMGNDFKILFVFPEPTNQ------------------HANPIMSHYSVL",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "          CCCCCCCCCCCCCcHHHHHHHhCCCCceehhhhcCCeeEEEeCCeEEEEECCCCcEEEEEEeCccccccCHHHHHHhhCCCceEEEECCEEEEEEEcCCCEEEEEEeCCCCCCCCCCceeEEEeC",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 10,  16,  16,  23,  24,  32,  32,  39,  46,  65,  65,  72,  72,
                               90,  90, 106, 107, 108, 112, 123, 123, 135],
                             [  8,  14,  17,  24,  24,  32,  48,  55,  55,  74,  77,  84,  96,
                              114, 116, 132, 132, 133, 133, 144, 162, 174]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['G', 'K', 'V', 'F', 'N', 'S', '-', '-', '-', 'D', 'F', 'P', 'A',
              'K', 'D', 'T', 'N', 'I', 'D', 'S', 'V', 'E', 'S', 'K', 'W', '-',
              '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',
              '-', '-', 'G', 'K', 'A', 'D', 'N', 'S', 'E', 'W', 'V', 'A', 'S',
              'A', 'K', 'G', 'L', 'Y', 'S', 'T', 'Y', 'S', 'K', 'H', 'N', 'I',
              'V', 'F', 'G', 'S', 'N', 'K', 'G', 'G', 'Q', '-', '-', '-', 'I',
              'F', 'E', 'V', 'R', 'S', 'L', '-', '-', '-', '-', '-', '-', '-',
              '-', '-', '-', '-', '-', 'D', 'K', 'Q', 'L', 'G', 'N', 'I', 'Y',
              'L', 'S', 'M', 'V', 'K', 'D', 'K', 'L', 'G', 'T', '-', '-', 'P',
              'Q', 'H', 'D', 'V', 'K', 'V', 'N', 'G', 'E', 'E', 'I', 'I', 'G',
              'Y', 'K', 'M', 'G', 'N', 'D', 'F', 'K', 'I', 'L', 'F', 'V', 'F',
              'P', 'E', 'P', 'T', 'N', 'Q', '-', '-', '-', '-', '-', '-', '-',
              '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', 'H', 'A',
              'N', 'P', 'I', 'M', 'S', 'H', 'Y', 'S', 'V', 'L'],
             ['E', 'R', 'S', 'L', 'G', 'N', 'E', 'Q', 'W', 'E', 'F', 'T', 'L',
              'G', 'M', 'P', '-', 'L', 'A', 'Q', 'A', 'V', 'A', 'I', 'L', 'Q',
              'K', 'H', 'C', 'R', 'I', 'I', 'K', 'N', 'V', 'Q', 'V', 'L', 'Y',
              'S', 'E', 'Q', 'S', 'P', 'L', 'S', 'H', 'D', '-', '-', '-', '-',
              '-', '-', '-', 'L', 'I', 'L', 'N', 'L', 'T', 'Q', 'D', 'G', 'I',
              'K', 'L', 'M', 'F', 'D', 'A', 'F', 'N', 'Q', 'R', 'L', 'K', 'V',
              'I', 'E', 'V', 'C', 'D', 'L', 'T', 'K', 'V', 'K', 'L', 'K', 'Y',
              'C', 'G', 'V', 'H', 'F', 'N', 'S', 'Q', 'A', 'I', 'A', 'P', 'T',
              'I', 'E', 'Q', 'I', 'D', 'Q', 'S', 'F', 'G', 'A', 'T', 'H', 'P',
              'G', 'V', 'Y', 'N', 'S', 'A', 'E', 'Q', 'L', 'F', 'H', 'L', 'N',
              'F', 'R', '-', 'G', '-', '-', '-', '-', 'L', 'S', 'F', 'S', 'F',
              'Q', 'L', 'D', 'S', 'W', 'T', 'E', 'A', 'P', 'K', 'Y', 'E', 'P',
              'N', 'F', 'A', 'H', 'G', 'L', 'A', 'S', 'L', 'Q', 'I', 'P', 'H',
              'G', 'A', 'T', 'V', 'K', 'R', 'M', 'Y', 'I', 'Y']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 33.21)
        self.assertAlmostEqual(alignment.annotations["E-value"], 69)
        self.assertAlmostEqual(alignment.annotations["Score"], 26.69)
        self.assertAlmostEqual(alignment.annotations["Identities"], 17)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.109)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 5.600)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "            CCCCceEEcCCcHHHHHHHHHHccCccccEEEEECCCCCcCCCeEEEecCCCeEEEecccc                                                                                                                                                                                                                                                                                                                                                             ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(
            alignment.query.seq[12:73],
            "GNEQWEFTLGMPLAQAVAILQKHCRIIKNVQVLYSEQSPLSHDLILNLTQDGIKLMFDAFN",
        )
        self.assertEqual(
            alignment[1],
            "GNEQWEFTLGMPLAQAVAILQKHCRIIK--NVQVLYSE-----QSPLSHDLILNLTQDGIKLMFDAFN",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "            G~~~~~f~LG~sL~~vi~~Lr~~~~~~~~Vel~Ys~~~Pl~~dIvi~L~~~GirL~Fd~~~                                                                                                                                                                                                                                                                                                                                                             ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                                                                                                     GwEHiE~Vvp~~~~~~~~~~~~~~~~~~~~~i~~k~s~pk~~~erl~NPtva~~~~~~~iKfHp~s        ",
        )
        self.assertEqual(alignment.target.id, "Q7N560_PHOLL/9")
        self.assertEqual(
            alignment.target.seq[101:167],
            "GWEHVELVLPVAPEILSTAAKALLPQPLPAGFSVKESQPKGEQERLPNPTLAITDGEITVKFHPFT",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "PF06185.11")
        self.assertEqual(
            alignment.target.annotations["hmm_description"], "; YecM ; YecM protein"
        )
        self.assertEqual(
            alignment[0],
            "GWEHVELVLPVAPEILSTAAKALLPQPLPAGFSVKESQPKGEQERLPNPTLAITD--GEITVKFHPFT",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                                                                                                     CeeEEEEEecCChHHHHHHHHHHCCCCCCCCcEEEEecCCCcccCCCCceEEEEcCCEEEEEccCC        ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[101, 129, 131, 139, 144, 156, 156, 167],
                             [ 12,  40,  40,  48,  48,  60,  62,  73]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['G', 'W', 'E', 'H', 'V', 'E', 'L', 'V', 'L', 'P', 'V', 'A', 'P',
              'E', 'I', 'L', 'S', 'T', 'A', 'A', 'K', 'A', 'L', 'L', 'P', 'Q',
              'P', 'L', 'P', 'A', 'G', 'F', 'S', 'V', 'K', 'E', 'S', 'Q', 'P',
              'K', 'G', 'E', 'Q', 'E', 'R', 'L', 'P', 'N', 'P', 'T', 'L', 'A',
              'I', 'T', 'D', '-', '-', 'G', 'E', 'I', 'T', 'V', 'K', 'F', 'H',
              'P', 'F', 'T'],
             ['G', 'N', 'E', 'Q', 'W', 'E', 'F', 'T', 'L', 'G', 'M', 'P', 'L',
              'A', 'Q', 'A', 'V', 'A', 'I', 'L', 'Q', 'K', 'H', 'C', 'R', 'I',
              'I', 'K', '-', '-', 'N', 'V', 'Q', 'V', 'L', 'Y', 'S', 'E', '-',
              '-', '-', '-', '-', 'Q', 'S', 'P', 'L', 'S', 'H', 'D', 'L', 'I',
              'L', 'N', 'L', 'T', 'Q', 'D', 'G', 'I', 'K', 'L', 'M', 'F', 'D',
              'A', 'F', 'N']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 30.70)
        self.assertAlmostEqual(alignment.annotations["E-value"], 18)
        self.assertAlmostEqual(alignment.annotations["Score"], 28.76)
        self.assertAlmostEqual(alignment.annotations["Identities"], 30)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.610)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 9.300)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                                                                                                                                                                                                                                              EEeCCchHHHHHhhcCCcEE                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(alignment.query.seq[238:258], "VYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(alignment[1], "VYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                                                                                                                                                                              I~fG~T~QDVl~~LG~P~~v                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                                                                                              i~~Gmt~~eV~~~lG~P~~~                                             ",
        )
        self.assertEqual(alignment.target.id, "Q8DTX1_STRMU/4")
        self.assertEqual(alignment.target.seq[94:114], "LKIGESYKKVVEKLGEPDVL")
        self.assertEqual(alignment.target.annotations["hmm_name"], "PF12978.6")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "; DUF3862 ; Domain of Unknown Function with PDB structure (DUF3862)",
        )
        self.assertEqual(alignment[0], "LKIGESYKKVVEKLGEPDVL")
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                                                                                              cCCCCCHHHHHHHHCCCCee                                             ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 94, 114],
                             [238, 258]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['L', 'K', 'I', 'G', 'E', 'S', 'Y', 'K', 'K', 'V', 'V', 'E', 'K',
              'L', 'G', 'E', 'P', 'D', 'V', 'L'],
             ['V', 'Y', 'F', 'G', 'D', 'S', 'C', 'Q', 'D', 'V', 'L', 'S', 'M',
              'L', 'G', 'S', 'P', 'H', 'K', 'V']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 29.70)
        self.assertAlmostEqual(alignment.annotations["E-value"], 83)
        self.assertAlmostEqual(alignment.annotations["Score"], 23.66)
        self.assertAlmostEqual(alignment.annotations["Identities"], 22)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.238)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 9.700)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "            CCCCceEEcCCcHHHHHHHHHHccCccccEEEEECCCCCcCCCeEEEecCCCeEEEeccccCceeEEEEec                                                                                                                                                                                                                                                                                                                                                   ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(
            alignment.query.seq[12:83],
            "GNEQWEFTLGMPLAQAVAILQKHCRIIKNVQVLYSEQSPLSHDLILNLTQDGIKLMFDAFNQRLKVIEVCD",
        )
        self.assertEqual(
            alignment[1],
            "GNEQWEFTLGMPLAQAVAILQKHCRI-----IKNVQVLYSEQSPLSHDLILNLTQDGIKLMFDAFNQRLKVIEVCD",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "            G~~~~~f~LG~sL~~vi~~Lr~~~~~~~~Vel~Ys~~~Pl~~dIvi~L~~~GirL~Fd~~~QrL~lIEV~d                                                                                                                                                                                                                                                                                                                                                   ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                                                            gi~iG~s~~~v~~~~G~p~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~d~~~g~v~~i~~~~          ",
        )
        self.assertEqual(alignment.target.id, "U5LAW6_9BACI/1")
        self.assertEqual(
            alignment.target.seq[60:130],
            "PFHIGQPVSEIYSSVFIDTNINFQYKGSSYRFELSEDDLNTRPLIKAGNIYAQLYIDRFTGELSSIRYMD",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "PF14504.5")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "; CAP_assoc_N ; CAP-associated N-terminal",
        )
        self.assertEqual(
            alignment[0],
            "P-----FHIGQPVSEIYSSVFIDTNINFQYKGSSYRFELSEDDLNTRPLIKA-GNIYAQLYIDRFTGELSSIRYMD",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                                                            CccCCCcHHHHHHhhCCCceEEEEEcCeEEEEEECcchhceeeeEEECCEEEEEEEeCCCCEEEEEEEec          ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 60,  61,  61,  81,  86, 107, 107, 130],
                             [ 12,  13,  18,  38,  38,  59,  60,  83]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['P', '-', '-', '-', '-', '-', 'F', 'H', 'I', 'G', 'Q', 'P', 'V',
              'S', 'E', 'I', 'Y', 'S', 'S', 'V', 'F', 'I', 'D', 'T', 'N', 'I',
              'N', 'F', 'Q', 'Y', 'K', 'G', 'S', 'S', 'Y', 'R', 'F', 'E', 'L',
              'S', 'E', 'D', 'D', 'L', 'N', 'T', 'R', 'P', 'L', 'I', 'K', 'A',
              '-', 'G', 'N', 'I', 'Y', 'A', 'Q', 'L', 'Y', 'I', 'D', 'R', 'F',
              'T', 'G', 'E', 'L', 'S', 'S', 'I', 'R', 'Y', 'M', 'D'],
             ['G', 'N', 'E', 'Q', 'W', 'E', 'F', 'T', 'L', 'G', 'M', 'P', 'L',
              'A', 'Q', 'A', 'V', 'A', 'I', 'L', 'Q', 'K', 'H', 'C', 'R', 'I',
              '-', '-', '-', '-', '-', 'I', 'K', 'N', 'V', 'Q', 'V', 'L', 'Y',
              'S', 'E', 'Q', 'S', 'P', 'L', 'S', 'H', 'D', 'L', 'I', 'L', 'N',
              'L', 'T', 'Q', 'D', 'G', 'I', 'K', 'L', 'M', 'F', 'D', 'A', 'F',
              'N', 'Q', 'R', 'L', 'K', 'V', 'I', 'E', 'V', 'C', 'D']],
            dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 28.54)
        self.assertAlmostEqual(alignment.annotations["E-value"], 1.2e02)
        self.assertAlmostEqual(alignment.annotations["Score"], 22.65)
        self.assertAlmostEqual(alignment.annotations["Identities"], 18)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.331)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 4.100)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                      CcHHHHHHHHHHccCccccEEEEECCCCCcCCCeEEEecCCCeEEEeccccCceeEEEEe                                                                                                                                                                                                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(
            alignment.query.seq[22:82],
            "MPLAQAVAILQKHCRIIKNVQVLYSEQSPLSHDLILNLTQDGIKLMFDAFNQRLKVIEVC",
        )
        self.assertEqual(
            alignment[1],
            "MPLAQAVAILQKHCRIIKNVQVLYSEQSPLSHDLILNL----TQDGIKLMFDAFNQRLK-VIEVC",
        )
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                      ~sL~~vi~~Lr~~~~~~~~Vel~Ys~~~Pl~~dIvi~L~~~GirL~Fd~~~QrL~lIEV~                                                                                                                                                                                                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "    m~~eE~ieYi~~nV~~gD~lEisy~Rv~vpGeVi~i~~~~~~~~~~~~~~l~l~ge~~~~~vevD                       ",
        )
        self.assertEqual(alignment.target.id, "A5UN43_METS3/2")
        self.assertEqual(
            alignment.target.seq[4:69],
            "LTPDKAVEYLKDNVKIHDNLEISYNRIFGSGEVLNMDFSEYFGKPGFKMLLSLDGDSINPTIEID",
        )
        self.assertEqual(alignment.target.annotations["hmm_name"], "PF09870.8")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "; DUF2097 ; Uncharacterized protein conserved in archaea (DUF2097)",
        )
        self.assertEqual(
            alignment[0],
            "LTPDKAVEYLKDNVKIHDNLEISYNRIFGSGEVLNMDFSEYFGKPGFKMLLSLDGDSINPTIEID",
        )
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "    CCHHHHHHHHHHhCCcCCeEEEEeceeeeeeeEEEEecccccCCCcEEEEEEECCCcccceEEEe                       ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 4, 42, 46, 63, 64, 69],
                             [22, 60, 60, 77, 77, 82]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['L', 'T', 'P', 'D', 'K', 'A', 'V', 'E', 'Y', 'L', 'K', 'D', 'N',
              'V', 'K', 'I', 'H', 'D', 'N', 'L', 'E', 'I', 'S', 'Y', 'N', 'R',
              'I', 'F', 'G', 'S', 'G', 'E', 'V', 'L', 'N', 'M', 'D', 'F', 'S',
              'E', 'Y', 'F', 'G', 'K', 'P', 'G', 'F', 'K', 'M', 'L', 'L', 'S',
              'L', 'D', 'G', 'D', 'S', 'I', 'N', 'P', 'T', 'I', 'E', 'I', 'D'],
             ['M', 'P', 'L', 'A', 'Q', 'A', 'V', 'A', 'I', 'L', 'Q', 'K', 'H',
              'C', 'R', 'I', 'I', 'K', 'N', 'V', 'Q', 'V', 'L', 'Y', 'S', 'E',
              'Q', 'S', 'P', 'L', 'S', 'H', 'D', 'L', 'I', 'L', 'N', 'L', '-',
              '-', '-', '-', 'T', 'Q', 'D', 'G', 'I', 'K', 'L', 'M', 'F', 'D',
              'A', 'F', 'N', 'Q', 'R', 'L', 'K', '-', 'V', 'I', 'E', 'V', 'C']],
            dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 26.31)
        self.assertAlmostEqual(alignment.annotations["E-value"], 22)
        self.assertAlmostEqual(alignment.annotations["Score"], 28.97)
        self.assertAlmostEqual(alignment.annotations["Identities"], 20)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.378)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 9.400)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                                                                                                                                                                                                                                              EEeCCchHHHHHhhcCCcEE                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(alignment.query.seq[238:258], "VYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(alignment[1], "VYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                                                                                                                                                                              I~fG~T~QDVl~~LG~P~~v                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                                       i~~G~t~~eV~~~lG~p~~~                                                                                                                             ",
        )
        self.assertEqual(alignment.target.id, "Q9KJ90_STREX/2")
        self.assertEqual(alignment.target.seq[39:59], "IQFGMTFDEVWEIGGGEAAC")
        self.assertEqual(alignment.target.annotations["hmm_name"], "PF07467.10")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "; BLIP ; Beta-lactamase inhibitor (BLIP)",
        )
        self.assertEqual(alignment[0], "IQFGMTFDEVWEIGGGEAAC")
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                                       cCCCCCHHHHHHHhCCCCee                                                                                                                             ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 39,  59],
                             [238, 258]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['I', 'Q', 'F', 'G', 'M', 'T', 'F', 'D', 'E', 'V', 'W', 'E', 'I',
              'G', 'G', 'G', 'E', 'A', 'A', 'C'],
             ['V', 'Y', 'F', 'G', 'D', 'S', 'C', 'Q', 'D', 'V', 'L', 'S', 'M',
              'L', 'G', 'S', 'P', 'H', 'K', 'V']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 25.89)
        self.assertAlmostEqual(alignment.annotations["E-value"], 32)
        self.assertAlmostEqual(alignment.annotations["Score"], 27.15)
        self.assertAlmostEqual(alignment.annotations["Identities"], 20)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.380)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 9.300)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                                                                                                                                                                                                                                              EEeCCchHHHHHhhcCCcEE                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(alignment.query.seq[238:258], "VYFGDSCQDVLSMLGSPHKV")
        self.assertEqual(alignment[1], "VYFG---------DSCQDVLSMLGSPHKV")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                                                                                                                                                                                                              I~fG~T~QDVl~~LG~P~~v                                                                                                                                                                    ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "           i~~g~~~~~~~~g~t~~~V~~~~G~p~~~                                                                                                                       ",
        )
        self.assertEqual(alignment.target.id, "Q8DTX1_STRMU/4")
        self.assertEqual(alignment.target.seq[11:40], "IKVTTDQNHFSGGTSIEQLKQWFGDPNKS")
        self.assertEqual(alignment.target.annotations["hmm_name"], "PF12978.6")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "; DUF3862 ; Domain of Unknown Function with PDB structure (DUF3862)",
        )
        self.assertEqual(alignment[0], "IKVTTDQNHFSGGTSIEQLKQWFGDPNKS")
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "           eeecCcccCCCCCCCHHHHHHHHCCCceE                                                                                                                       ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 11,  15,  24,  40],
                             [238, 242, 242, 258]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['I', 'K', 'V', 'T', 'T', 'D', 'Q', 'N', 'H', 'F', 'S', 'G', 'G',
              'T', 'S', 'I', 'E', 'Q', 'L', 'K', 'Q', 'W', 'F', 'G', 'D', 'P',
              'N', 'K', 'S'],
             ['V', 'Y', 'F', 'G', '-', '-', '-', '-', '-', '-', '-', '-', '-',
              'D', 'S', 'C', 'Q', 'D', 'V', 'L', 'S', 'M', 'L', 'G', 'S', 'P',
              'H', 'K', 'V']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertAlmostEqual(alignment.annotations["Probab"], 20.88)
        self.assertAlmostEqual(alignment.annotations["E-value"], 78)
        self.assertAlmostEqual(alignment.annotations["Score"], 19.81)
        self.assertAlmostEqual(alignment.annotations["Identities"], 26)
        self.assertAlmostEqual(alignment.annotations["Similarity"], 0.368)
        self.assertAlmostEqual(alignment.annotations["Sum_probs"], 0.0)
        self.assertAlmostEqual(alignment.annotations["Template_Neff"], 7.400)
        self.assertEqual(
            alignment.query.letter_annotations["ss_pred"],
            "                                                             CCCeEEEeccccCceeEEEEecCC                                                                                                                                                                                                                                                                                                                                                 ",
        )
        self.assertEqual(
            alignment.query.id,
            "sp|Q9BSU1|CP070_HUMAN UPF0183 protein C16orf70 OS=Homo sapiens OX=9606 GN=C16orf70 PE=1 SV=1",
        )
        self.assertEqual(alignment.query.seq[61:85], "QDGIKLMFDAFNQRLKVIEVCDLT")
        self.assertEqual(alignment[1], "QDGIKLMFDAFNQRLKVIEVCDLT")
        self.assertEqual(
            alignment.query.letter_annotations["Consensus"],
            "                                                             ~~GirL~Fd~~~QrL~lIEV~d~s                                                                                                                                                                                                                                                                                                                                                 ",
        )
        self.assertEqual(
            alignment.target.letter_annotations["Consensus"],
            "                         ~~~i~ld~d~~g~ivgiEil~~s",
        )
        self.assertEqual(alignment.target.id, "Q21BS2_RHOPB/5")
        self.assertEqual(alignment.target.seq[25:48], "APNVIFDYDAEGRIVGIELLDAR")
        self.assertEqual(alignment.target.annotations["hmm_name"], "PF10049.8")
        self.assertEqual(
            alignment.target.annotations["hmm_description"],
            "; DUF2283 ; Protein of unknown function (DUF2283)",
        )
        self.assertEqual(alignment[0], "APNVIFDYDA-EGRIVGIELLDAR")
        self.assertEqual(
            alignment.target.letter_annotations["ss_pred"],
            "                         cCCEEEEECCCCCEEEEEEEecC",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[25, 35, 35, 48],
                             [61, 71, 72, 85]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['A', 'P', 'N', 'V', 'I', 'F', 'D', 'Y', 'D', 'A', '-', 'E', 'G',
              'R', 'I', 'V', 'G', 'I', 'E', 'L', 'L', 'D', 'A', 'R'],
             ['Q', 'D', 'G', 'I', 'K', 'L', 'M', 'F', 'D', 'A', 'F', 'N', 'Q',
              'R', 'L', 'K', 'V', 'I', 'E', 'V', 'C', 'D', 'L', 'T']],
            dtype='U')
                # fmt: on
            )
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_length(self):
        """Test getting the number of alignments without parsing the file."""
        stream = open(self.path)
        alignments = Align.parse(stream, "hhr")
        stream.close()
        self.assertEqual(len(alignments), 16)


class Align_hhr_2uvo_hhblits_emptytable(unittest.TestCase):
    path = os.path.join("HHsuite", "2uvo_hhblits_emptytable.hhr")

    def test_reading(self):
        alignments = Align.parse(self.path, "hhr")
        self.assertEqual(alignments.metadata["Match_columns"], 171)
        self.assertEqual(alignments.metadata["No_of_seqs"], (1560, 4005))
        self.assertAlmostEqual(alignments.metadata["Neff"], 8.3)
        self.assertEqual(alignments.metadata["Searched_HMMs"], 34)
        self.assertEqual(alignments.metadata["Rundate"], "Fri Feb 15 16:34:13 2019")
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_length(self):
        """Test getting the number of alignments without parsing the file."""
        stream = open(self.path)
        alignments = Align.parse(stream, "hhr")
        stream.close()
        self.assertEqual(len(alignments), 0)


class Align_hhr_2uvo_hhblits_onlyheader(unittest.TestCase):
    path = os.path.join("HHsuite", "2uvo_hhblits_onlyheader.hhr")

    def test_reading(self):
        with self.assertRaises(ValueError) as cm:
            alignments = Align.parse(self.path, "hhr")
        self.assertEqual(str(cm.exception), "Truncated file.")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
