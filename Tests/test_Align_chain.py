# Copyright 2023 by Michiel de Hoon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Align.chain module."""
import unittest
from io import StringIO


from Bio.Align import Alignment
from Bio.Seq import Seq, reverse_complement
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SimpleLocation, ExactPosition, CompoundLocation
from Bio import SeqIO
from Bio import Align


try:
    import numpy as np
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install numpy if you want to use Bio.Align.chain."
    ) from None


class TestAlign_dna_rna(unittest.TestCase):
    # The chain file dna_rna.chain was generated from the PSL file using:
    # pslToChain dna_rna.psl dna_rna.chain

    def setUp(self):
        data = {}
        records = SeqIO.parse("Blat/dna.fa", "fasta")
        for record in records:
            name, start_end = record.id.split(":")
            assert name == "chr3"
            start, end = start_end.split("-")
            start = int(start)
            end = int(end)
            sequence = str(record.seq)
            assert len(sequence) == end - start
            data[start] = sequence
        self.dna = data
        records = SeqIO.parse("Blat/rna.fa", "fasta")
        self.rna = {record.id: record.seq for record in records}

    def test_reading(self):
        """Test parsing dna_rna.chain."""
        path = "Blat/dna_rna.chain"
        alignments = Align.parse(path, "chain")
        self.check_alignments(alignments)
        alignments.rewind()
        self.check_alignments(alignments)
        with Align.parse(path, "chain") as alignments:
            self.check_alignments(alignments)
        with self.assertRaises(AttributeError):
            alignments._stream
        with Align.parse(path, "chain") as alignments:
            pass
        with self.assertRaises(AttributeError):
            alignments._stream

    def check_alignments(self, alignments):
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 1711))
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.score, 176)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr3")
        self.assertEqual(alignment.query.id, "NR_046654.1")
        self.assertEqual(len(alignment.target.seq), 198295559)
        self.assertEqual(len(alignment.query.seq), 181)
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[42530895, 42530958, 42532020, 42532095, 42532563, 42532606],
                          [     181,      118,      118,       43,       43,        0]])
                # fmt: on
            )
        )
        dna = Seq(self.dna, length=len(alignment.target))
        alignment.target.seq = dna
        alignment.query.seq = self.rna[alignment.query.id]
        self.assertTrue(
            np.array_equal(
                alignment.substitutions,
                # fmt: off
            np.array([[36.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0., 40.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0., 57.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0., 42.,  0.,  0.,  0.,  0.],
                      [ 2.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  1.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  3.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                     ])
                # fmt: on
            )
        )
        self.assertEqual(alignment.substitutions.alphabet, "ACGTacgt")
        self.assertEqual(
            str(alignment),
            """\
chr3       42530895 CGGAAGTACTTCTGGGGGTACATACTCATCGGCTGGGGTATGGTACCAGGGAGGGCTTCC
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
NR_046654       181 CGGAAGTACTTCTGGGGGTACATACTCATCGGCTGGGGTATGGTACCAGGGAGGGCTTCC

chr3       42530955 AGGCTGGGGACAGAGGGGGCAAGGCCTGGAGAACTCCCTAGGGGGAGGGTGCCAACCCAG
                 60 |||---------------------------------------------------------
NR_046654       121 AGG---------------------------------------------------------

chr3       42531015 CTTGCAGTCCTACGTCTTGCTTAGCTGCAGGTCCTGCCTGCAAGGATATCAGCCAAGGGT
                120 ------------------------------------------------------------
NR_046654       118 ------------------------------------------------------------

chr3       42531075 CAAGAAAGTCCTCAAAATGTCTGATCCCAGGACAAGTCCCTCAGGTTGCAGCTGCACCTA
                180 ------------------------------------------------------------
NR_046654       118 ------------------------------------------------------------

chr3       42531135 GGGCTGACCTGTGGGACAGATTTTGTGAACATCTTTCCATTTCCCTTTAGTTCCCGAAAT
                240 ------------------------------------------------------------
NR_046654       118 ------------------------------------------------------------

chr3       42531195 ACAcagggccactgctaatctataaagggcctctgtcacaattagaaagagaatgtccgt
                300 ------------------------------------------------------------
NR_046654       118 ------------------------------------------------------------

chr3       42531255 ctaggtagacacagcccttcaggcatacagcttCACCCCCTCAGTGGAGCATCCCTCCGT
                360 ------------------------------------------------------------
NR_046654       118 ------------------------------------------------------------

chr3       42531315 GGTGAACAACCTATGCAACCAAAGACAGCAGACTGACAACCCACCCTTTtctctctccct
                420 ------------------------------------------------------------
NR_046654       118 ------------------------------------------------------------

chr3       42531375 ccctctccctgcttttctccaaaatctctccctcatgccctctacccctgcttcctgtgc
                480 ------------------------------------------------------------
NR_046654       118 ------------------------------------------------------------

chr3       42531435 cctctctgctctttcactctccctGGGCCTGACAGGGGTACCCAGCACATTCACCATGGT
                540 ------------------------------------------------------------
NR_046654       118 ------------------------------------------------------------

chr3       42531495 GTGGACCATCGCCAGGATCCATTTTGAGGATTATGGGTGAGCTGCTGCCCCACACACTCC
                600 ------------------------------------------------------------
NR_046654       118 ------------------------------------------------------------

chr3       42531555 CCCGGCCGCCATCACTTGGGCAGGCCCCCTGGGTGGGATGATAATGCCATCTGGCCTTGG
                660 ------------------------------------------------------------
NR_046654       118 ------------------------------------------------------------

chr3       42531615 TGAGTGGACAAAAACCACAGCTCTCGGGCCAGAGGGGAGGCTGGAGGAGGACCTGGGGAG
                720 ------------------------------------------------------------
NR_046654       118 ------------------------------------------------------------

chr3       42531675 CAACAGACTCTGGGCCCGGGGTTGCTAAAGTGCTCAGGAGCAGAGCTGGGGACAACTGGG
                780 ------------------------------------------------------------
NR_046654       118 ------------------------------------------------------------

chr3       42531735 GGAGGTGCTGCTGAGTCTCTCTCTGGCTGAGGACAATCCCTCTCATTCCTCCCCACGGTC
                840 ------------------------------------------------------------
NR_046654       118 ------------------------------------------------------------

chr3       42531795 TGCTCAGGTGCTGGGACACCATCAACTCCTCACTGTGGTGGATCATAAAGGGCCCCATCC
                900 ------------------------------------------------------------
NR_046654       118 ------------------------------------------------------------

chr3       42531855 TCACCTCCATCTTGGTAAGATaccctcccaccacctagagatggggaaacaggcccaaag
                960 ------------------------------------------------------------
NR_046654       118 ------------------------------------------------------------

chr3       42531915 ggcaggcaacttagcccaaggtcacatgggaaattagtatctaggtcagaactgaaacgt
               1020 ------------------------------------------------------------
NR_046654       118 ------------------------------------------------------------

chr3       42531975 agcttcctaatgcccaatgcaggatcatccccacccctgtcctaccagTTCTTCCTTGAG
               1080 ---------------------------------------------...||||||||||||
NR_046654       118 ---------------------------------------------CAGTTCTTCCTTGAG

chr3       42532035 CGTAAGCGGATTGGGAGCACAGTCCTTAGGGATTTGAAGGAGGTAGAGTTCCCGGATGAC
               1140 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
NR_046654       103 CGTAAGCGGATTGGGAGCACAGTCCTTAGGGATTTGAAGGAGGTAGAGTTCCCGGATGAC

chr3       42532095 CTGCCCAAAGGGGAAATGCCAGAGGAGAGGTAAGATAGAGAGAGGGGCAGCAGGACCCTG
               1200 ------------------------------------------------------------
NR_046654        43 ------------------------------------------------------------

chr3       42532155 GGAAAGAAGACAGGCCAGCAGTCAAGGGGCCTGAACACCTCAGCCTTCCCGCTCTGACTG
               1260 ------------------------------------------------------------
NR_046654        43 ------------------------------------------------------------

chr3       42532215 CCCGAACTCGGGTCCCCACCCACTAGGTAAACTTCATCCTGTTTATTTGCATCATCCGAA
               1320 ------------------------------------------------------------
NR_046654        43 ------------------------------------------------------------

chr3       42532275 TCCTGCTTCAGAAACTGCGGCCCCCAGATATCAGGAAGAGTGACAGCAGTCCATACTCGT
               1380 ------------------------------------------------------------
NR_046654        43 ------------------------------------------------------------

chr3       42532335 GAGTGTGGGCCTAGTGCCTCAGCCCCCAGTACCTCCATCCCCAGTCCTCAAATCATCCCA
               1440 ------------------------------------------------------------
NR_046654        43 ------------------------------------------------------------

chr3       42532395 CATCTCCTTGAAGTCCTCCCACCCCAAACATCCAGAGTCACCAAAGAGCCACATTGTTCT
               1500 ------------------------------------------------------------
NR_046654        43 ------------------------------------------------------------

chr3       42532455 TTCCCACCTCCACCATGGCCTGGCTcagcccaccaccatcccctgctccagccccaccct
               1560 ------------------------------------------------------------
NR_046654        43 ------------------------------------------------------------

chr3       42532515 caCCAGGCTGCACTCAGAGCCCTGCATGCTTCTCCTGCCCACACTCACCTAGCATCCTTC
               1620 ------------------------------------------------||||||||||||
NR_046654        43 ------------------------------------------------CTAGCATCCTTC

chr3       42532575 CCAGGTATGCATCTGCTGCCAAGCCAGGgag 42532606
               1680 ||||||||||||||||||||||||||||...     1711
NR_046654        31 CCAGGTATGCATCTGCTGCCAAGCCAGGGAG        0
""",
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 176 chr3 198295559 + 42530895 42532606 NR_046654.1 181 - 0 181 1
63	1062	0
75	468	0
43

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 1714))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr3")
        self.assertEqual(alignment.query.id, "NR_046654.1_modified")
        self.assertEqual(len(alignment.target.seq), 198295559)
        self.assertEqual(len(alignment.query.seq), 190)
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[42530895, 42530922, 42530922, 42530958, 42532020,
                           42532037, 42532039, 42532095, 42532563, 42532606],
                          [     185,      158,      155,      119,      119,
                                102,      102,       46,       46,        3],
                         ])
                # fmt: on
            )
        )
        dna = Seq(self.dna, length=len(alignment.target))
        alignment.target.seq = dna
        alignment.query.seq = self.rna[alignment.query.id]
        self.assertTrue(
            np.array_equal(
                alignment.substitutions,
                # fmt: off
            np.array([[34.,  0.,  0.,  1.,  0.,  0.,  0.,  0.],
                      [ 0., 40.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0., 57.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0., 41.,  0.,  0.,  0.,  0.],
                      [ 2.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  1.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  3.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                     ]),
                # fmt: on
            )
        )
        self.assertEqual(alignment.substitutions.alphabet, "ACGTacgt")
        self.assertEqual(
            str(alignment),
            """\
chr3       42530895 CGGAAGTACTTCTGGGGGTACATACTC---ATCGGCTGGGGTATGGTACCAGGGAGGGCT
                  0 |||||||||||||||||||||||||||---||||||||||||||||||||||||||||||
NR_046654       185 CGGAAGTACTTCTGGGGGTACATACTCCCCATCGGCTGGGGTATGGTACCAGGGAGGGCT

chr3       42530952 TCCAGGCTGGGGACAGAGGGGGCAAGGCCTGGAGAACTCCCTAGGGGGAGGGTGCCAACC
                 60 ||||||------------------------------------------------------
NR_046654       125 TCCAGG------------------------------------------------------

chr3       42531012 CAGCTTGCAGTCCTACGTCTTGCTTAGCTGCAGGTCCTGCCTGCAAGGATATCAGCCAAG
                120 ------------------------------------------------------------
NR_046654       119 ------------------------------------------------------------

chr3       42531072 GGTCAAGAAAGTCCTCAAAATGTCTGATCCCAGGACAAGTCCCTCAGGTTGCAGCTGCAC
                180 ------------------------------------------------------------
NR_046654       119 ------------------------------------------------------------

chr3       42531132 CTAGGGCTGACCTGTGGGACAGATTTTGTGAACATCTTTCCATTTCCCTTTAGTTCCCGA
                240 ------------------------------------------------------------
NR_046654       119 ------------------------------------------------------------

chr3       42531192 AATACAcagggccactgctaatctataaagggcctctgtcacaattagaaagagaatgtc
                300 ------------------------------------------------------------
NR_046654       119 ------------------------------------------------------------

chr3       42531252 cgtctaggtagacacagcccttcaggcatacagcttCACCCCCTCAGTGGAGCATCCCTC
                360 ------------------------------------------------------------
NR_046654       119 ------------------------------------------------------------

chr3       42531312 CGTGGTGAACAACCTATGCAACCAAAGACAGCAGACTGACAACCCACCCTTTtctctctc
                420 ------------------------------------------------------------
NR_046654       119 ------------------------------------------------------------

chr3       42531372 cctccctctccctgcttttctccaaaatctctccctcatgccctctacccctgcttcctg
                480 ------------------------------------------------------------
NR_046654       119 ------------------------------------------------------------

chr3       42531432 tgccctctctgctctttcactctccctGGGCCTGACAGGGGTACCCAGCACATTCACCAT
                540 ------------------------------------------------------------
NR_046654       119 ------------------------------------------------------------

chr3       42531492 GGTGTGGACCATCGCCAGGATCCATTTTGAGGATTATGGGTGAGCTGCTGCCCCACACAC
                600 ------------------------------------------------------------
NR_046654       119 ------------------------------------------------------------

chr3       42531552 TCCCCCGGCCGCCATCACTTGGGCAGGCCCCCTGGGTGGGATGATAATGCCATCTGGCCT
                660 ------------------------------------------------------------
NR_046654       119 ------------------------------------------------------------

chr3       42531612 TGGTGAGTGGACAAAAACCACAGCTCTCGGGCCAGAGGGGAGGCTGGAGGAGGACCTGGG
                720 ------------------------------------------------------------
NR_046654       119 ------------------------------------------------------------

chr3       42531672 GAGCAACAGACTCTGGGCCCGGGGTTGCTAAAGTGCTCAGGAGCAGAGCTGGGGACAACT
                780 ------------------------------------------------------------
NR_046654       119 ------------------------------------------------------------

chr3       42531732 GGGGGAGGTGCTGCTGAGTCTCTCTCTGGCTGAGGACAATCCCTCTCATTCCTCCCCACG
                840 ------------------------------------------------------------
NR_046654       119 ------------------------------------------------------------

chr3       42531792 GTCTGCTCAGGTGCTGGGACACCATCAACTCCTCACTGTGGTGGATCATAAAGGGCCCCA
                900 ------------------------------------------------------------
NR_046654       119 ------------------------------------------------------------

chr3       42531852 TCCTCACCTCCATCTTGGTAAGATaccctcccaccacctagagatggggaaacaggccca
                960 ------------------------------------------------------------
NR_046654       119 ------------------------------------------------------------

chr3       42531912 aagggcaggcaacttagcccaaggtcacatgggaaattagtatctaggtcagaactgaaa
               1020 ------------------------------------------------------------
NR_046654       119 ------------------------------------------------------------

chr3       42531972 cgtagcttcctaatgcccaatgcaggatcatccccacccctgtcctaccagTTCTTCCTT
               1080 ------------------------------------------------...|||||||||
NR_046654       119 ------------------------------------------------CAGTTCTTCCTT

chr3       42532032 GAGCGTAAGCGGATTGGGAGCACAGTCCTTAGGGATTTGAAGGAGGTAGAGTTCCCGGAT
               1140 |||||--|||||||||||.|||||||||||||||||||||||||||||||||||||||||
NR_046654       107 GAGCG--AGCGGATTGGGTGCACAGTCCTTAGGGATTTGAAGGAGGTAGAGTTCCCGGAT

chr3       42532092 GACCTGCCCAAAGGGGAAATGCCAGAGGAGAGGTAAGATAGAGAGAGGGGCAGCAGGACC
               1200 |||---------------------------------------------------------
NR_046654        49 GAC---------------------------------------------------------

chr3       42532152 CTGGGAAAGAAGACAGGCCAGCAGTCAAGGGGCCTGAACACCTCAGCCTTCCCGCTCTGA
               1260 ------------------------------------------------------------
NR_046654        46 ------------------------------------------------------------

chr3       42532212 CTGCCCGAACTCGGGTCCCCACCCACTAGGTAAACTTCATCCTGTTTATTTGCATCATCC
               1320 ------------------------------------------------------------
NR_046654        46 ------------------------------------------------------------

chr3       42532272 GAATCCTGCTTCAGAAACTGCGGCCCCCAGATATCAGGAAGAGTGACAGCAGTCCATACT
               1380 ------------------------------------------------------------
NR_046654        46 ------------------------------------------------------------

chr3       42532332 CGTGAGTGTGGGCCTAGTGCCTCAGCCCCCAGTACCTCCATCCCCAGTCCTCAAATCATC
               1440 ------------------------------------------------------------
NR_046654        46 ------------------------------------------------------------

chr3       42532392 CCACATCTCCTTGAAGTCCTCCCACCCCAAACATCCAGAGTCACCAAAGAGCCACATTGT
               1500 ------------------------------------------------------------
NR_046654        46 ------------------------------------------------------------

chr3       42532452 TCTTTCCCACCTCCACCATGGCCTGGCTcagcccaccaccatcccctgctccagccccac
               1560 ------------------------------------------------------------
NR_046654        46 ------------------------------------------------------------

chr3       42532512 cctcaCCAGGCTGCACTCAGAGCCCTGCATGCTTCTCCTGCCCACACTCACCTAGCATCC
               1620 ---------------------------------------------------|||||||||
NR_046654        46 ---------------------------------------------------CTAGCATCC

chr3       42532572 TTCCCAGGTATGCATCTGCTGCCAAGCCAGGgag 42532606
               1680 |||||||||||||||||||||||||||||||...     1714
NR_046654        37 TTCCCAGGTATGCATCTGCTGCCAAGCCAGGGAG        3
""",
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 170 chr3 198295559 + 42530895 42532606 NR_046654.1_modified 190 - 5 187 2
27	0	3
36	1062	0
17	2	0
56	468	0
43

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 5407))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr3")
        self.assertEqual(alignment.query.id, "NR_111921.1")
        self.assertEqual(len(alignment.target.seq), 198295559)
        self.assertEqual(len(alignment.query.seq), 216)
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array( [[48663767, 48663813, 48665640, 48665722, 48669098, 48669174],
                              [       0,        46,      46,      128,      128,      204]]),
                # fmt: on
            )
        )
        dna = Seq(self.dna, length=len(alignment.target.seq))
        alignment.target.seq = dna
        alignment.query.seq = self.rna[alignment.query.id]
        self.assertTrue(
            np.array_equal(
                alignment.substitutions,
                # fmt: off
            np.array([[53.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0., 35.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0., 50.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0., 27.,  0.,  0.,  0.,  0.],
                      [ 9.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  7.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0., 16.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0.,  7.,  0.,  0.,  0.,  0.],
                     ])
                # fmt: on
            )
        )
        self.assertEqual(alignment.substitutions.alphabet, "ACGTacgt")
        self.assertEqual(
            str(alignment),
            """\
chr3       48663767 CACGAGAGGAGCGGAGGCGAGGGGTGAACGCGGAGCACTCCAATCGGTCAGTCATTGTTT
                  0 ||||||||||||||||||||||||||||||||||||||||||||||--------------
NR_111921         0 CACGAGAGGAGCGGAGGCGAGGGGTGAACGCGGAGCACTCCAATCG--------------

chr3       48663827 CTATTGGCACAATGGGAGGCCCCGCCCCTCACGGCGGACTCATCGCATGGGGGAGGGGGC
                 60 ------------------------------------------------------------
NR_111921        46 ------------------------------------------------------------

chr3       48663887 TCCGCGGGTTGCCGGCTAACCGTGAGAGAGTCCGGGAGGTACACTATACGGACCGGCCTC
                120 ------------------------------------------------------------
NR_111921        46 ------------------------------------------------------------

chr3       48663947 CAAAGGCGGAATCGATAACGAGCTGCAGCGCCGGGTGCAGAGGACGCGGGCATCCCGAAG
                180 ------------------------------------------------------------
NR_111921        46 ------------------------------------------------------------

chr3       48664007 CCCAGGAAGAGGTCAGGGCCGGGACCCCAGAACGCTCCACAGGGTGCGGCTCCCGCGATG
                240 ------------------------------------------------------------
NR_111921        46 ------------------------------------------------------------

chr3       48664067 GGGTGGATCCTGGTTCTAACAGGCGAGGAACTCCTGGCCAAGGCCTCTGGCCCGCCCCGA
                300 ------------------------------------------------------------
NR_111921        46 ------------------------------------------------------------

chr3       48664127 ACGGTCCCTATGACATCACCATCAACCAATCAGTCGGCGCATCCTTTCGCCCCTTGACTG
                360 ------------------------------------------------------------
NR_111921        46 ------------------------------------------------------------

chr3       48664187 CTCCGCTTCCGGGAGGCGGGGCTTCTGCGGGTTCCACCTCCCGAGCGCCCCTTGTGGCTA
                420 ------------------------------------------------------------
NR_111921        46 ------------------------------------------------------------

chr3       48664247 CCAAGGTCAGGCAACAGGTGTCCAGTTGTTCCCTCTCCTGTCTACGAATCTGAGGACCTC
                480 ------------------------------------------------------------
NR_111921        46 ------------------------------------------------------------

chr3       48664307 CCCAGGATCAGAGCTCTGGGCCTGATACACGGCCGGGGTTCCTACGGGTTTGTGAGTGGG
                540 ------------------------------------------------------------
NR_111921        46 ------------------------------------------------------------

chr3       48664367 GGTGGAAGATCTGCAGAGGCACTTAGGGCTGAACTCCTTTGAATGGGAGCCAATCGGTGC
                600 ------------------------------------------------------------
NR_111921        46 ------------------------------------------------------------

chr3       48664427 AGGGGCTGGAGGAGCGAGTCCCCCAAAGTAGttttatttatctatttagagacaaggtct
                660 ------------------------------------------------------------
NR_111921        46 ------------------------------------------------------------

chr3       48664487 cactctttcggagtgcagtggtgatcacagctcaccgtagcctcgaactccccaggcgat
                720 ------------------------------------------------------------
NR_111921        46 ------------------------------------------------------------

chr3       48664547 tctctcacctcagcctcccgagtagctgggactacgggtacatgtcatcacacttggcta
                780 ------------------------------------------------------------
NR_111921        46 ------------------------------------------------------------

chr3       48664607 atttttgcattttttatagagacagggtctcaccatgtaggccagattagtcttgaactc
                840 ------------------------------------------------------------
NR_111921        46 ------------------------------------------------------------

chr3       48664667 ctgggctcaagcaatccgcccatcttggcctcccaaagtgctgggattataggtgtgagc
                900 ------------------------------------------------------------
NR_111921        46 ------------------------------------------------------------

chr3       48664727 caccgcgcccggcAACCCAGAAGTGGTTTTGACAGCAccagcgctttctgtgtccacaat
                960 ------------------------------------------------------------
NR_111921        46 ------------------------------------------------------------

chr3       48664787 ctagtgagtagagggcacaaaacctgacaccacggaggcagacaggcaggggctctgccg
               1020 ------------------------------------------------------------
NR_111921        46 ------------------------------------------------------------

chr3       48664847 gggaagggtgttggagtcccaaaggaggcgtctgagtcaccttcgcaacctgggacgcct
               1080 ------------------------------------------------------------
NR_111921        46 ------------------------------------------------------------

chr3       48664907 tcttgcataagatgcctgagcagtgccttgaatgaccaaggggagatccgcatctgcaaa
               1140 ------------------------------------------------------------
NR_111921        46 ------------------------------------------------------------

chr3       48664967 ggaagggcagggagggatagggattgggggtgggcatcctaggtcttggagactgtgtgg
               1200 ------------------------------------------------------------
NR_111921        46 ------------------------------------------------------------

chr3       48665027 gcaaatgtgcagagacataaagggactatggctgagggaaatcaagCCCTGCCCTCTCAC
               1260 ------------------------------------------------------------
NR_111921        46 ------------------------------------------------------------

chr3       48665087 CAATAGGGCTGGCGCTGGTCCCAGCTAACACTCCTTTTGGAGAGCAAAGCTCCTCACTTC
               1320 ------------------------------------------------------------
NR_111921        46 ------------------------------------------------------------

chr3       48665147 TGAGTAGTGAGATTGATTGCGGATCACTCTCCATGTTGCTGCCTGCTGTGTGTCATCCCA
               1380 ------------------------------------------------------------
NR_111921        46 ------------------------------------------------------------

chr3       48665207 CTGTCATCCTCCCTTTGTGGCTGTTCTGTGGAGCCCCTCTCCCTCAATCTGCACTCACCT
               1440 ------------------------------------------------------------
NR_111921        46 ------------------------------------------------------------

chr3       48665267 CTATGCCCCAGCCCCATTGGCAGCTCCTAATGCACTCCCGGTaaaaaaaaaaaaacaaaa
               1500 ------------------------------------------------------------
NR_111921        46 ------------------------------------------------------------

chr3       48665327 aCCAGATGTTAGTGATAGTGGTGGTAGTTCTTCTCTCCACCTCCAAATCTTGCCCTTGCC
               1560 ------------------------------------------------------------
NR_111921        46 ------------------------------------------------------------

chr3       48665387 TCCTAATAAGACCCCTATGTGGTTTAACCTCAttttttttttttttttttttttttgaga
               1620 ------------------------------------------------------------
NR_111921        46 ------------------------------------------------------------

chr3       48665447 tggagtttcactctgtcacccaggctggagtgaagtggtgtgatGGGGCTTCACCATGTg
               1680 ------------------------------------------------------------
NR_111921        46 ------------------------------------------------------------

chr3       48665507 atggggcttcaccatgttggccaggctggtatcaaactcctgacctctagtgatctgccc
               1740 ------------------------------------------------------------
NR_111921        46 ------------------------------------------------------------

chr3       48665567 gcctcagcctcccaaagtgctgggattaccggcatgaggcaccgtgcccagccTATCCTC
               1800 ------------------------------------------------------------
NR_111921        46 ------------------------------------------------------------

chr3       48665627 CTTCTCTTATCAGCTCCCAACTAGAGGTCCACCCAGGACCCAGAGACCTGGATTTGAGGC
               1860 -------------|||||||||||||||||||||||||||||||||||||||||||||||
NR_111921        46 -------------CTCCCAACTAGAGGTCCACCCAGGACCCAGAGACCTGGATTTGAGGC

chr3       48665687 TGCTGGGCGGCAGATGGAGCGATCAGAAGACCAGGGTAAGGGTGTGGCAGATACTGCCAC
               1920 |||||||||||||||||||||||||||||||||||-------------------------
NR_111921        93 TGCTGGGCGGCAGATGGAGCGATCAGAAGACCAGG-------------------------

chr3       48665747 TAACACTTCTCAGCCTTTCCTTCTCCTGCCTTTTCCACCCCACCCTGTGTTTGTCTACTC
               1980 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48665807 CCAGCCAGGTGTACCTTTCCAGGGGAAGACCTGGCCAACCTGTCCAGCTCAATTAtccag
               2040 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48665867 cagttctttgacctcactgagatctcgagtccattgttcatcacctcagctattgacctg
               2100 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48665927 tgtcattagccttatagagttcagtgccacggaaactccctgccctgttctttttctttt
               2160 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48665987 tctttttttttttttttttgagacagagccttgctctgtcgcccaggctggagtgcagtg
               2220 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48666047 gcgcgatctcggctcactgcaagctctgcctcccaggttcacaccattctcctgactcag
               2280 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48666107 cctcccgagtagctgggactacaggcgtccaccaccatgcccagctaatttttttttttg
               2340 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48666167 tatttttagtagagacggcgtttcaccgtgttagccaggctggtctcgatctcctgacct
               2400 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48666227 tgtgatgctcccgcctcggcctcccaaagtgctgggattacaggcatgagccattgtgcc
               2460 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48666287 cggcctgccctgttcttcttagacaaacttgctgggctaaaatctaaccccgttaaaata
               2520 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48666347 gactatttacgtattgtttgcctctagcgcagcagaacattgctggagaaaaacaaacaa
               2580 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48666407 ccgtgctaattggtctcattttatattcatgaccacaagcctcagtattatatcggaggg
               2640 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48666467 cctatccagtgcagtagggcaagaaaaataataagttatgaagattggaagggaaaaaaa
               2700 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48666527 actaattcacaagcagtaggattgtatatgtaaaaatttcaaaggaacctataggtaagt
               2760 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48666587 tgttagaatgagttcagcaaagttgttggacacaagatcaatatataaaaatcagttgca
               2820 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48666647 atttctatatgtcaccaacagttagaaaataaatttcttgcctgggcatgttggctcaag
               2880 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48666707 cctgtaatcccagcactttgggtggccaaggcgggcagatcacctgaggtcaggagtttg
               2940 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48666767 agaccagtttggccaacatggtgaaatcccgtctctactaaaaatacagaaattagccgg
               3000 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48666827 gcgtggtggtgggcacctgtagtcccagctactgaggaggctgaggcaggagaatcactt
               3060 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48666887 gaacctgggaggcagaggttgcagtgaacgagaaaaaaaaattttttttcttaaaaacaa
               3120 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48666947 tgatgtttacaatagcatcaagtaatatcaaatgctgaggaataaacctaatgaaagatg
               3180 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48667007 tgcaaagactacatacacacatacaaaaaaactataaaacattattgagggaaataaaga
               3240 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48667067 cataggcctggcattggtggctcatgcctgaaatctcagcactttggagggccaaggtgg
               3300 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48667127 gtggatcatttgaggtcaggagttagagatcagtccggccaacatggtgaaacctcatct
               3360 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48667187 ctactaaaaatacaaaaaaattagcttggccaggtgcagtggctcacacctgtaatccca
               3420 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48667247 gcactttgggaggctgaggcgggcggatcatgaggtcaggagatcgagaccatcctggct
               3480 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48667307 aacacggtgaaaccctgtctctactaaaaatacaaaaaaaaattagccgggcctgatggc
               3540 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48667367 gggcgcccgtagtcccagctactcgggaggctgaggtagcagaatggcgtgaacctggga
               3600 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48667427 ggtgcagcttgcagtgagcctaaattgcgccactgcactccagcctgggtaacagagcga
               3660 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48667487 gactccgtttcaaaaaaaaaaaaaaaaattagctgggcatgctgttgtgcacctgcaatc
               3720 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48667547 ccagctactctggaggatgaggcagaagtgcctgaacctgggacacagaggttgcagtga
               3780 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48667607 gccaagatcatgccattgcactccagcctggacaacacagccagacgctatctgaaaaaa
               3840 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48667667 aaaaaaaaaaaaaaagtaaaaaaaatgagaaataaagacataaataaagtgaaaaattgt
               3900 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48667727 tccaatattggaaaagtcaatattataaaggtgccaattttcccaaattgatatatggat
               3960 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48667787 tcgatgcaacttcagttaaaaatcccactaaattttggctgggtgcggtggctcacacct
               4020 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48667847 gtaatcccagcactttgggaggctgaggcgggcggatcacaaggtcaggagatcgagacc
               4080 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48667907 atcttggctaacatggtgaaaccgtctctactaaaaatacaaaagttagccgggtgtggt
               4140 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48667967 ggcgggcacctgtagtcccagctacttgggaggctgagacagaatggcgtgaacctgggg
               4200 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48668027 aggcggagcttgcagtgagccaagttgacgccactgcactccagcctgggcgacagagca
               4260 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48668087 agactctgtctcaaaaaaaaaaaaaaaaaaaTCCCACTAGATTTTGTGTGTGTGTAAACT
               4320 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48668147 GACAAACTAGATTTAGcagcctgagcaacacagcaaaaccccatctctacaaaaaataca
               4380 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48668207 aaaattttgcacatgcctgtatagtcccagctacttgggaggctgaagtgggaggatcat
               4440 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48668267 gtgagctctggggaggtcgaggctgtagtgagctatgatcacatgctgcactctagcctg
               4500 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48668327 ggcaacagagcaagagaccctgtatctaaaaaaagaatgaaaattaaaaaataaaaaGAa
               4560 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48668387 accaagattgtgtggtactggtacgaggataggaagactaaaggaacgaaatccagagac
               4620 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48668447 aggcctgaagatgtgtggaaacttgaattttgacaagggtgGTTCTTCAGAGCTAACATG
               4680 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48668507 AAGAAAGGGTTGTTTTCTTTTTTTTGTTTCCCcaggagcaactctattaactgaaagaat
               4740 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48668567 aggcttttcaataaatgatgctgggtcagttggatatccatatagaaaaaattaaatgag
               4800 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48668627 atctctatttcacactgcttgcataatcaattccatataaatttgacatctgaaaatata
               4860 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48668687 cagtttctagaaaacagtatTAAGACCttgttttgttttttgttgttgttgttttttgtt
               4920 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48668747 ttgttttttgttttttgagacagagtctcgctctgtcgccaggctggaatacagtggtgc
               4980 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48668807 aaccttggctcactgcaacctctgactccctagttcaagcaattctcctgcctcagcctc
               5040 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48668867 ccgagtagctgcgattacaggcacatgccaccacgcccagctaatttttgtatttttagt
               5100 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48668927 agagatgggggtttcaccatgttggccaggatggtctcgatctcctgaccctgtaatccg
               5160 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48668987 cccacctcggcctcccaaagtgctgggattacaggcgtgagccactgcacctggccAAGA
               5220 ------------------------------------------------------------
NR_111921       128 ------------------------------------------------------------

chr3       48669047 GAAGATCTTAAAGGTGACTTTAAGCAAACttttttttttttttttttacagagacgggag
               5280 ---------------------------------------------------.........
NR_111921       128 ---------------------------------------------------AGACGGGAG

chr3       48669107 ctggagtgcagtggctgttcacaagcgtgaAAGCAAAGATTAAAAAATTTGTTTTTATAT
               5340 ..............................||||||||||||||||||||||||||||||
NR_111921       137 CTGGAGTGCAGTGGCTGTTCACAAGCGTGAAAGCAAAGATTAAAAAATTTGTTTTTATAT

chr3       48669167 TAAAAAA 48669174
               5400 |||||||     5407
NR_111921       197 TAAAAAA      204
""",
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 182 chr3 198295559 + 48663767 48669174 NR_111921.1 216 + 0 204 3
46	1827	0
82	3376	0
76

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 5409))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr3")
        self.assertEqual(alignment.query.id, "NR_111921.1_modified")
        self.assertEqual(len(alignment.target.seq), 198295559)
        self.assertEqual(len(alignment.query.seq), 220)
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[48663767, 48663795, 48663796, 48663813, 48665640,
                           48665716, 48665716, 48665722, 48669098, 48669174],
                          [       3,       31,       31,       48,       48,
                                124,      126,      132,      132,      208]
                         ])
                # fmt: on
            )
        )
        dna = Seq(self.dna, length=len(alignment.target))
        alignment.target.seq = dna
        alignment.query.seq = self.rna[alignment.query.id]
        self.assertTrue(
            np.array_equal(
                alignment.substitutions,
                # fmt: off
            np.array([[53.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0., 34.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  2., 48.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0., 27.,  0.,  0.,  0.,  0.],
                      [ 9.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  7.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0., 16.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0.,  7.,  0.,  0.,  0.,  0.],
                     ]),
                # fmt: on
            )
        )
        self.assertEqual(alignment.substitutions.alphabet, "ACGTacgt")
        self.assertEqual(
            str(alignment),
            """\
chr3       48663767 CACGAGAGGAGCGGAGGCGAGGGGTGAACGCGGAGCACTCCAATCGGTCAGTCATTGTTT
                  0 ||||||||||||||||||||||||||||-|||||||||||||||||--------------
NR_111921         3 CACGAGAGGAGCGGAGGCGAGGGGTGAA-GCGGAGCACTCCAATCG--------------

chr3       48663827 CTATTGGCACAATGGGAGGCCCCGCCCCTCACGGCGGACTCATCGCATGGGGGAGGGGGC
                 60 ------------------------------------------------------------
NR_111921        48 ------------------------------------------------------------

chr3       48663887 TCCGCGGGTTGCCGGCTAACCGTGAGAGAGTCCGGGAGGTACACTATACGGACCGGCCTC
                120 ------------------------------------------------------------
NR_111921        48 ------------------------------------------------------------

chr3       48663947 CAAAGGCGGAATCGATAACGAGCTGCAGCGCCGGGTGCAGAGGACGCGGGCATCCCGAAG
                180 ------------------------------------------------------------
NR_111921        48 ------------------------------------------------------------

chr3       48664007 CCCAGGAAGAGGTCAGGGCCGGGACCCCAGAACGCTCCACAGGGTGCGGCTCCCGCGATG
                240 ------------------------------------------------------------
NR_111921        48 ------------------------------------------------------------

chr3       48664067 GGGTGGATCCTGGTTCTAACAGGCGAGGAACTCCTGGCCAAGGCCTCTGGCCCGCCCCGA
                300 ------------------------------------------------------------
NR_111921        48 ------------------------------------------------------------

chr3       48664127 ACGGTCCCTATGACATCACCATCAACCAATCAGTCGGCGCATCCTTTCGCCCCTTGACTG
                360 ------------------------------------------------------------
NR_111921        48 ------------------------------------------------------------

chr3       48664187 CTCCGCTTCCGGGAGGCGGGGCTTCTGCGGGTTCCACCTCCCGAGCGCCCCTTGTGGCTA
                420 ------------------------------------------------------------
NR_111921        48 ------------------------------------------------------------

chr3       48664247 CCAAGGTCAGGCAACAGGTGTCCAGTTGTTCCCTCTCCTGTCTACGAATCTGAGGACCTC
                480 ------------------------------------------------------------
NR_111921        48 ------------------------------------------------------------

chr3       48664307 CCCAGGATCAGAGCTCTGGGCCTGATACACGGCCGGGGTTCCTACGGGTTTGTGAGTGGG
                540 ------------------------------------------------------------
NR_111921        48 ------------------------------------------------------------

chr3       48664367 GGTGGAAGATCTGCAGAGGCACTTAGGGCTGAACTCCTTTGAATGGGAGCCAATCGGTGC
                600 ------------------------------------------------------------
NR_111921        48 ------------------------------------------------------------

chr3       48664427 AGGGGCTGGAGGAGCGAGTCCCCCAAAGTAGttttatttatctatttagagacaaggtct
                660 ------------------------------------------------------------
NR_111921        48 ------------------------------------------------------------

chr3       48664487 cactctttcggagtgcagtggtgatcacagctcaccgtagcctcgaactccccaggcgat
                720 ------------------------------------------------------------
NR_111921        48 ------------------------------------------------------------

chr3       48664547 tctctcacctcagcctcccgagtagctgggactacgggtacatgtcatcacacttggcta
                780 ------------------------------------------------------------
NR_111921        48 ------------------------------------------------------------

chr3       48664607 atttttgcattttttatagagacagggtctcaccatgtaggccagattagtcttgaactc
                840 ------------------------------------------------------------
NR_111921        48 ------------------------------------------------------------

chr3       48664667 ctgggctcaagcaatccgcccatcttggcctcccaaagtgctgggattataggtgtgagc
                900 ------------------------------------------------------------
NR_111921        48 ------------------------------------------------------------

chr3       48664727 caccgcgcccggcAACCCAGAAGTGGTTTTGACAGCAccagcgctttctgtgtccacaat
                960 ------------------------------------------------------------
NR_111921        48 ------------------------------------------------------------

chr3       48664787 ctagtgagtagagggcacaaaacctgacaccacggaggcagacaggcaggggctctgccg
               1020 ------------------------------------------------------------
NR_111921        48 ------------------------------------------------------------

chr3       48664847 gggaagggtgttggagtcccaaaggaggcgtctgagtcaccttcgcaacctgggacgcct
               1080 ------------------------------------------------------------
NR_111921        48 ------------------------------------------------------------

chr3       48664907 tcttgcataagatgcctgagcagtgccttgaatgaccaaggggagatccgcatctgcaaa
               1140 ------------------------------------------------------------
NR_111921        48 ------------------------------------------------------------

chr3       48664967 ggaagggcagggagggatagggattgggggtgggcatcctaggtcttggagactgtgtgg
               1200 ------------------------------------------------------------
NR_111921        48 ------------------------------------------------------------

chr3       48665027 gcaaatgtgcagagacataaagggactatggctgagggaaatcaagCCCTGCCCTCTCAC
               1260 ------------------------------------------------------------
NR_111921        48 ------------------------------------------------------------

chr3       48665087 CAATAGGGCTGGCGCTGGTCCCAGCTAACACTCCTTTTGGAGAGCAAAGCTCCTCACTTC
               1320 ------------------------------------------------------------
NR_111921        48 ------------------------------------------------------------

chr3       48665147 TGAGTAGTGAGATTGATTGCGGATCACTCTCCATGTTGCTGCCTGCTGTGTGTCATCCCA
               1380 ------------------------------------------------------------
NR_111921        48 ------------------------------------------------------------

chr3       48665207 CTGTCATCCTCCCTTTGTGGCTGTTCTGTGGAGCCCCTCTCCCTCAATCTGCACTCACCT
               1440 ------------------------------------------------------------
NR_111921        48 ------------------------------------------------------------

chr3       48665267 CTATGCCCCAGCCCCATTGGCAGCTCCTAATGCACTCCCGGTaaaaaaaaaaaaacaaaa
               1500 ------------------------------------------------------------
NR_111921        48 ------------------------------------------------------------

chr3       48665327 aCCAGATGTTAGTGATAGTGGTGGTAGTTCTTCTCTCCACCTCCAAATCTTGCCCTTGCC
               1560 ------------------------------------------------------------
NR_111921        48 ------------------------------------------------------------

chr3       48665387 TCCTAATAAGACCCCTATGTGGTTTAACCTCAttttttttttttttttttttttttgaga
               1620 ------------------------------------------------------------
NR_111921        48 ------------------------------------------------------------

chr3       48665447 tggagtttcactctgtcacccaggctggagtgaagtggtgtgatGGGGCTTCACCATGTg
               1680 ------------------------------------------------------------
NR_111921        48 ------------------------------------------------------------

chr3       48665507 atggggcttcaccatgttggccaggctggtatcaaactcctgacctctagtgatctgccc
               1740 ------------------------------------------------------------
NR_111921        48 ------------------------------------------------------------

chr3       48665567 gcctcagcctcccaaagtgctgggattaccggcatgaggcaccgtgcccagccTATCCTC
               1800 ------------------------------------------------------------
NR_111921        48 ------------------------------------------------------------

chr3       48665627 CTTCTCTTATCAGCTCCCAACTAGAGGTCCACCCAGGACCCAGAGACCTGGATTTGAGGC
               1860 -------------|||||||||||||||||||||||||||||||||||||||||||||||
NR_111921        48 -------------CTCCCAACTAGAGGTCCACCCAGGACCCAGAGACCTGGATTTGAGGC

chr3       48665687 TGCTGGGCGGCAGATGGAGCGATCAGAAG--ACCAGGGTAAGGGTGTGGCAGATACTGCC
               1920 |||||..||||||||||||||||||||||--||||||-----------------------
NR_111921        95 TGCTGCCCGGCAGATGGAGCGATCAGAAGCCACCAGG-----------------------

chr3       48665745 ACTAACACTTCTCAGCCTTTCCTTCTCCTGCCTTTTCCACCCCACCCTGTGTTTGTCTAC
               1980 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48665805 TCCCAGCCAGGTGTACCTTTCCAGGGGAAGACCTGGCCAACCTGTCCAGCTCAATTAtcc
               2040 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48665865 agcagttctttgacctcactgagatctcgagtccattgttcatcacctcagctattgacc
               2100 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48665925 tgtgtcattagccttatagagttcagtgccacggaaactccctgccctgttctttttctt
               2160 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48665985 tttctttttttttttttttttgagacagagccttgctctgtcgcccaggctggagtgcag
               2220 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48666045 tggcgcgatctcggctcactgcaagctctgcctcccaggttcacaccattctcctgactc
               2280 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48666105 agcctcccgagtagctgggactacaggcgtccaccaccatgcccagctaatttttttttt
               2340 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48666165 tgtatttttagtagagacggcgtttcaccgtgttagccaggctggtctcgatctcctgac
               2400 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48666225 cttgtgatgctcccgcctcggcctcccaaagtgctgggattacaggcatgagccattgtg
               2460 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48666285 cccggcctgccctgttcttcttagacaaacttgctgggctaaaatctaaccccgttaaaa
               2520 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48666345 tagactatttacgtattgtttgcctctagcgcagcagaacattgctggagaaaaacaaac
               2580 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48666405 aaccgtgctaattggtctcattttatattcatgaccacaagcctcagtattatatcggag
               2640 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48666465 ggcctatccagtgcagtagggcaagaaaaataataagttatgaagattggaagggaaaaa
               2700 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48666525 aaactaattcacaagcagtaggattgtatatgtaaaaatttcaaaggaacctataggtaa
               2760 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48666585 gttgttagaatgagttcagcaaagttgttggacacaagatcaatatataaaaatcagttg
               2820 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48666645 caatttctatatgtcaccaacagttagaaaataaatttcttgcctgggcatgttggctca
               2880 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48666705 agcctgtaatcccagcactttgggtggccaaggcgggcagatcacctgaggtcaggagtt
               2940 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48666765 tgagaccagtttggccaacatggtgaaatcccgtctctactaaaaatacagaaattagcc
               3000 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48666825 gggcgtggtggtgggcacctgtagtcccagctactgaggaggctgaggcaggagaatcac
               3060 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48666885 ttgaacctgggaggcagaggttgcagtgaacgagaaaaaaaaattttttttcttaaaaac
               3120 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48666945 aatgatgtttacaatagcatcaagtaatatcaaatgctgaggaataaacctaatgaaaga
               3180 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48667005 tgtgcaaagactacatacacacatacaaaaaaactataaaacattattgagggaaataaa
               3240 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48667065 gacataggcctggcattggtggctcatgcctgaaatctcagcactttggagggccaaggt
               3300 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48667125 gggtggatcatttgaggtcaggagttagagatcagtccggccaacatggtgaaacctcat
               3360 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48667185 ctctactaaaaatacaaaaaaattagcttggccaggtgcagtggctcacacctgtaatcc
               3420 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48667245 cagcactttgggaggctgaggcgggcggatcatgaggtcaggagatcgagaccatcctgg
               3480 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48667305 ctaacacggtgaaaccctgtctctactaaaaatacaaaaaaaaattagccgggcctgatg
               3540 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48667365 gcgggcgcccgtagtcccagctactcgggaggctgaggtagcagaatggcgtgaacctgg
               3600 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48667425 gaggtgcagcttgcagtgagcctaaattgcgccactgcactccagcctgggtaacagagc
               3660 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48667485 gagactccgtttcaaaaaaaaaaaaaaaaattagctgggcatgctgttgtgcacctgcaa
               3720 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48667545 tcccagctactctggaggatgaggcagaagtgcctgaacctgggacacagaggttgcagt
               3780 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48667605 gagccaagatcatgccattgcactccagcctggacaacacagccagacgctatctgaaaa
               3840 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48667665 aaaaaaaaaaaaaaaaagtaaaaaaaatgagaaataaagacataaataaagtgaaaaatt
               3900 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48667725 gttccaatattggaaaagtcaatattataaaggtgccaattttcccaaattgatatatgg
               3960 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48667785 attcgatgcaacttcagttaaaaatcccactaaattttggctgggtgcggtggctcacac
               4020 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48667845 ctgtaatcccagcactttgggaggctgaggcgggcggatcacaaggtcaggagatcgaga
               4080 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48667905 ccatcttggctaacatggtgaaaccgtctctactaaaaatacaaaagttagccgggtgtg
               4140 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48667965 gtggcgggcacctgtagtcccagctacttgggaggctgagacagaatggcgtgaacctgg
               4200 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48668025 ggaggcggagcttgcagtgagccaagttgacgccactgcactccagcctgggcgacagag
               4260 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48668085 caagactctgtctcaaaaaaaaaaaaaaaaaaaTCCCACTAGATTTTGTGTGTGTGTAAA
               4320 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48668145 CTGACAAACTAGATTTAGcagcctgagcaacacagcaaaaccccatctctacaaaaaata
               4380 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48668205 caaaaattttgcacatgcctgtatagtcccagctacttgggaggctgaagtgggaggatc
               4440 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48668265 atgtgagctctggggaggtcgaggctgtagtgagctatgatcacatgctgcactctagcc
               4500 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48668325 tgggcaacagagcaagagaccctgtatctaaaaaaagaatgaaaattaaaaaataaaaaG
               4560 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48668385 Aaaccaagattgtgtggtactggtacgaggataggaagactaaaggaacgaaatccagag
               4620 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48668445 acaggcctgaagatgtgtggaaacttgaattttgacaagggtgGTTCTTCAGAGCTAACA
               4680 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48668505 TGAAGAAAGGGTTGTTTTCTTTTTTTTGTTTCCCcaggagcaactctattaactgaaaga
               4740 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48668565 ataggcttttcaataaatgatgctgggtcagttggatatccatatagaaaaaattaaatg
               4800 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48668625 agatctctatttcacactgcttgcataatcaattccatataaatttgacatctgaaaata
               4860 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48668685 tacagtttctagaaaacagtatTAAGACCttgttttgttttttgttgttgttgttttttg
               4920 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48668745 ttttgttttttgttttttgagacagagtctcgctctgtcgccaggctggaatacagtggt
               4980 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48668805 gcaaccttggctcactgcaacctctgactccctagttcaagcaattctcctgcctcagcc
               5040 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48668865 tcccgagtagctgcgattacaggcacatgccaccacgcccagctaatttttgtattttta
               5100 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48668925 gtagagatgggggtttcaccatgttggccaggatggtctcgatctcctgaccctgtaatc
               5160 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48668985 cgcccacctcggcctcccaaagtgctgggattacaggcgtgagccactgcacctggccAA
               5220 ------------------------------------------------------------
NR_111921       132 ------------------------------------------------------------

chr3       48669045 GAGAAGATCTTAAAGGTGACTTTAAGCAAACttttttttttttttttttacagagacggg
               5280 -----------------------------------------------------.......
NR_111921       132 -----------------------------------------------------AGACGGG

chr3       48669105 agctggagtgcagtggctgttcacaagcgtgaAAGCAAAGATTAAAAAATTTGTTTTTAT
               5340 ................................||||||||||||||||||||||||||||
NR_111921       139 AGCTGGAGTGCAGTGGCTGTTCACAAGCGTGAAAGCAAAGATTAAAAAATTTGTTTTTAT

chr3       48669165 ATTAAAAAA 48669174
               5400 |||||||||     5409
NR_111921       199 ATTAAAAAA      208
""",
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 175 chr3 198295559 + 48663767 48669174 NR_111921.1_modified 220 + 3 208 4
28	1	0
17	1827	0
76	0	2
6	3376	0
76

""",
        )
        self.assertRaises(StopIteration, next, alignments)

    def test_writing(self):
        """Test writing the alignments in dna_rna.chain."""
        path = "Blat/dna_rna.chain"
        alignments = Align.parse(path, "chain")
        stream = StringIO()
        n = Align.write(alignments, stream, "chain")
        self.assertEqual(n, 4)
        stream.seek(0)
        alignments = Align.parse(stream, "chain")
        self.check_alignments(alignments)


class TestAlign_dna(unittest.TestCase):
    queries = {
        record.id: str(record.seq)
        for record in SeqIO.parse("Blat/fasta_34.fa", "fasta")
    }

    def test_reading_psl_34_001(self):
        """Test parsing psl_34_001.chain."""
        # The chain file psl_34_001.chain was generated from the PSL file using:
        # pslToChain psl_34_001.psl psl_34_001.chain
        path = "Blat/psl_34_001.chain"
        alignments = Align.parse(path, "chain")
        self.check_reading_psl_34_001(alignments)

    def check_reading_psl_34_001(self, alignments):
        """Check parsing psl_34_001.chain."""
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 16))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertEqual(len(alignment.target.seq), 191154276)
        self.assertEqual(len(alignment.query.seq), 33)
        self.assertEqual(
            str(alignment),
            """\
chr4       61646095 ???????????????? 61646111
                  0 ||||||||||||||||       16
hg18_dna         11 ????????????????       27
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[61646095, 61646111],
                          [      11,       27]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 16 chr4 191154276 + 61646095 61646111 hg18_dna 33 + 11 27 1
16

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 33))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 33)
        self.assertEqual(
            str(alignment),
            """\
chr1       10271783 ????????????????????????????????? 10271816
                  0 |||||||||||||||||||||||||||||||||       33
hg18_dna          0 ?????????????????????????????????       33
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[10271783, 10271816],
                          [       0,       33]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 33 chr1 249250621 + 10271783 10271816 hg18_dna 33 + 0 33 2
33

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 17))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 33)
        self.assertEqual(
            str(alignment),
            """\
chr2       53575980 ????????????????? 53575997
                  0 |||||||||||||||||       17
hg18_dna         25 ?????????????????        8
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[53575980, 53575997],
                          [      25,        8]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 17 chr2 243199373 + 53575980 53575997 hg18_dna 33 - 8 25 3
17

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 41))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr9")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 141213431)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr9       85737865 ????????????????????????????????????????? 85737906
                  0 |||||||||||||||||||||||||||||||||||||||||       41
hg19_dna          9 ?????????????????????????????????????????       50
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[85737865, 85737906],
                          [       9,       50]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 35 chr9 141213431 + 85737865 85737906 hg19_dna 50 + 9 50 4
41

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 41))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr8")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 146364022)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr8       95160479 ????????????????????????????????????????? 95160520
                  0 |||||||||||||||||||||||||||||||||||||||||       41
hg19_dna          8 ?????????????????????????????????????????       49
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[95160479, 95160520],
                          [       8,       49]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 41 chr8 146364022 + 95160479 95160520 hg19_dna 50 + 8 49 5
41

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr22")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 51304566)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr22      42144400 ???????????????????????????????????? 42144436
                  0 ||||||||||||||||||||||||||||||||||||       36
hg19_dna         11 ????????????????????????????????????       47
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[42144400, 42144436],
                          [      11,       47]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 30 chr22 51304566 + 42144400 42144436 hg19_dna 50 + 11 47 6
36

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 48))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr2      183925984 ??????----?????????????????????????????????????? 183926028
                  0 ||||||----||||||||||||||||||||||||||||||||||||||        48
hg19_dna          1 ????????????????????????????????????????????????        49
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[183925984, 183925990, 183925990, 183926028],
                          [        1,         7,        11,        49]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 41 chr2 243199373 + 183925984 183926028 hg19_dna 50 + 1 49 7
6	0	4
38

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 170))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr19      35483340 ????????????????????????????????????????????????????????????
                  0 |||||||||||||||||||||||||-----------------------------------
hg19_dna         10 ?????????????????????????-----------------------------------

chr19      35483400 ????????????????????????????????????????????????????????????
                 60 ------------------------------------------------------------
hg19_dna         35 ------------------------------------------------------------

chr19      35483460 ?????????????????????????????????????????????????? 35483510
                120 ---------------------------------------|||||||||||      170
hg19_dna         35 ---------------------------------------???????????       46
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[35483340, 35483365, 35483499, 35483510],
                          [      10,       35,       35,       46]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 31 chr19 59128983 + 35483340 35483510 hg19_dna 50 + 10 46 8
25	134	0
11

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr18")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 78077248)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr18      23891310 ??????????????????????????????????????? 23891349
                  0 |||||||||||||||||||||||||||||||||||||||       39
hg19_dna         10 ???????????????????????????????????????       49
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[23891310, 23891349],
                          [      10,       49]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 39 chr18 78077248 + 23891310 23891349 hg19_dna 50 + 10 49 9
39

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 28))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr18")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 78077248)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr18      43252217 ???????????????????????????? 43252245
                  0 ||||||||||||||||||||||||||||       28
hg19_dna         21 ????????????????????????????       49
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[43252217, 43252245],
                          [      21,       49]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 26 chr18 78077248 + 43252217 43252245 hg19_dna 50 + 21 49 10
28

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 54))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 115169878)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr13      52759147 ?????????????---??????????????????????????????????????
                  0 |||||||---------||||||||||||||||||||||||||||||||||||||
hg19_dna          1 ???????------?????????????????????????????????????????

chr13      52759198
                 54
hg19_dna         49
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[52759147, 52759154, 52759160, 52759160, 52759198],
                          [       1,        8,        8,       11,       49]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 41 chr13 115169878 + 52759147 52759198 hg19_dna 50 + 1 49 11
7	6	3
38

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 50))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr1        1207056 ?????????????????????????????????????????????????? 1207106
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||      50
hg19_dna          0 ??????????????????????????????????????????????????      50
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[1207056, 1207106],
                          [      0,      50]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 50 chr1 249250621 + 1207056 1207106 hg19_dna 50 + 0 50 12
50

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr1       61700837 ?????????????????????????????????? 61700871
                  0 ||||||||||||||||||||||||||||||||||       34
hg19_dna          1 ??????????????????????????????????       35
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[61700837, 61700871],
                          [       1,       35]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 28 chr1 249250621 + 61700837 61700871 hg19_dna 50 + 1 35 13
34

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 44))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 191154276)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr4       37558157 ????????????????----------?????????????????? 37558191
                  0 ||||||||||----------------||||||||||||||||||       44
hg19_dna         49 ??????????------????????????????????????????       11
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[37558157, 37558167, 37558173, 37558173, 37558191],
                          [      49,       39,       39,       29,       11]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
28	0	0	0	1	10	1	6	-	hg19_dna	50	11	49	chr4	191154276	37558157	37558191	2	10,18,	1,21,	37558157,37558173,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 37))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr22")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 51304566)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr22      48997405 ????????????????????????????????????? 48997442
                  0 |||||||||||||||||||||||||||||||||||||       37
hg19_dna         49 ?????????????????????????????????????       12
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[48997405, 48997442],
                          [      49,       12]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 33 chr22 51304566 + 48997405 48997442 hg19_dna 50 - 1 38 15
37

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr2      120641740 ???????????????????????????????????? 120641776
                  0 ||||||||||||||||||||||||||||||||||||        36
hg19_dna         49 ????????????????????????????????????        13
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[120641740, 120641776],
                          [       49,        13]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 34 chr2 243199373 + 120641740 120641776 hg19_dna 50 - 1 37 16
36

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr19      54017130 ??????????????????????????????????????? 54017169
                  0 |||||||||||||||||||||||||||||||||||||||       39
hg19_dna         49 ???????????????????????????????????????       10
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[54017130, 54017169],
                          [      49,       10]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 39 chr19 59128983 + 54017130 54017169 hg19_dna 50 - 1 40 17
39

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr19        553742 ??????????????????????????????????????? 553781
                  0 |||||||||||||||||||||||||||||||||||||||     39
hg19_dna         49 ???????????????????????????????????????     10
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[553742, 553781],
                          [    49,     10]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 33 chr19 59128983 + 553742 553781 hg19_dna 50 - 1 40 18
39

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr10")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 135534747)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr10      99388555 ???????????????????????????????????? 99388591
                  0 ||||||||||||||||||||||||||||||||||||       36
hg19_dna         49 ????????????????????????????????????       13
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[99388555, 99388591],
                          [      49,       13]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 30 chr10 135534747 + 99388555 99388591 hg19_dna 50 - 1 37 19
36

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 25))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr10")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 135534747)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr10     112178171 ????????????????????????? 112178196
                  0 |||||||||||||||||||||||||        25
hg19_dna         35 ?????????????????????????        10
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[112178171, 112178196],
                          [       35,        10]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 23 chr10 135534747 + 112178171 112178196 hg19_dna 50 - 15 40 20
25

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr1       39368490 ???????????????????????????????????? 39368526
                  0 ||||||||||||||||||||||||||||||||||||       36
hg19_dna         49 ????????????????????????????????????       13
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[39368490, 39368526],
                          [      49,       13]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 34 chr1 249250621 + 39368490 39368526 hg19_dna 50 - 1 37 21
36

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr1      220325687 ?????????????????????????????????? 220325721
                  0 ||||||||||||||||||||||||||||||||||        34
hg19_dna         47 ??????????????????????????????????        13
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[220325687, 220325721],
                          [       47,        13]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 32 chr1 249250621 + 220325687 220325721 hg19_dna 50 - 3 37 22
34

""",
        )
        self.assertRaises(StopIteration, next, alignments)

    def test_writing_chain_34_001(self):
        """Test writing the alignments in psl_34_001.chain."""
        path = "Blat/psl_34_001.chain"
        with open(path) as stream:
            original_data = stream.read()
        alignments = Align.parse(path, "chain")
        stream = StringIO()
        n = Align.write(alignments, stream, "chain")
        self.assertEqual(n, 22)
        stream.seek(0)
        alignments = Align.parse(stream, "chain")
        self.check_reading_psl_34_001(alignments)

    def test_reading_chain_34_002(self):
        """Test parsing psl_34_002.chain."""
        # The chain file psl_34_002.chain was generated from the PSL file using:
        # pslToChain psl_34_002.psl psl_34_002.chain
        path = "Blat/psl_34_002.chain"
        alignments = Align.parse(path, "chain")
        self.assertRaises(StopIteration, next, alignments)

    def test_writing_psl_34_002(self):
        """Test writing the alignments in psl_34_002.chain."""
        path = "Blat/psl_34_002.chain"
        alignments = Align.parse(path, "chain")
        stream = StringIO()
        n = Align.write(alignments, stream, "chain")
        self.assertEqual(n, 0)
        stream.seek(0)
        alignments = Align.parse(stream, "chain")
        self.assertRaises(StopIteration, next, alignments)

    def test_reading_psl_34_003(self):
        """Test parsing psl_34_003.chain."""
        # The chain file psl_34_003.chain was generated from the PSL file using:
        # pslToChain psl_34_003.psl psl_34_003.chain
        path = "Blat/psl_34_003.chain"
        alignments = Align.parse(path, "chain")
        self.check_reading_psl_34_003(alignments)

    def check_reading_psl_34_003(self, alignments):
        """Check parsing psl_34_003.chain."""
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 16))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertEqual(len(alignment.target.seq), 191154276)
        self.assertEqual(len(alignment.query.seq), 33)
        self.assertEqual(
            str(alignment),
            """\
chr4       61646095 ???????????????? 61646111
                  0 ||||||||||||||||       16
hg18_dna         11 ????????????????       27
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[61646095, 61646111],
                          [      11,       27]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 16 chr4 191154276 + 61646095 61646111 hg18_dna 33 + 11 27 1
16

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 33))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 33)
        self.assertEqual(
            str(alignment),
            """\
chr1       10271783 ????????????????????????????????? 10271816
                  0 |||||||||||||||||||||||||||||||||       33
hg18_dna          0 ?????????????????????????????????       33
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[10271783, 10271816],
                          [       0,       33]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 33 chr1 249250621 + 10271783 10271816 hg18_dna 33 + 0 33 2
33

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 17))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 33)
        self.assertEqual(
            str(alignment),
            """\
chr2       53575980 ????????????????? 53575997
                  0 |||||||||||||||||       17
hg18_dna         25 ?????????????????        8
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[53575980, 53575997],
                          [      25,        8]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 17 chr2 243199373 + 53575980 53575997 hg18_dna 33 - 8 25 3
17

""",
        )
        self.assertRaises(StopIteration, next, alignments)

    def test_writing_psl_34_003(self):
        """Test writing the alignments in psl_34_003.chain."""
        path = "Blat/psl_34_003.chain"
        alignments = Align.parse(path, "chain")
        stream = StringIO()
        n = Align.write(alignments, stream, "chain")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "chain")
        self.check_reading_psl_34_003(alignments)

    def test_reading_psl_34_004(self):
        """Test parsing psl_34_004.chain."""
        # The chain file psl_34_004.chain was generated from the PSL file using:
        # pslToChain psl_34_004.psl psl_34_004.chain
        path = "Blat/psl_34_004.chain"
        alignments = Align.parse(path, "chain")
        self.check_reading_psl_34_004(alignments)

    def check_reading_psl_34_004(self, alignments):
        """Check parsing psl_34_004.chain."""
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 41))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr9")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 141213431)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr9       85737865 ????????????????????????????????????????? 85737906
                  0 |||||||||||||||||||||||||||||||||||||||||       41
hg19_dna          9 ?????????????????????????????????????????       50
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[85737865, 85737906],
                          [       9,       50]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 35 chr9 141213431 + 85737865 85737906 hg19_dna 50 + 9 50 1
41

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 41))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr8")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 146364022)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr8       95160479 ????????????????????????????????????????? 95160520
                  0 |||||||||||||||||||||||||||||||||||||||||       41
hg19_dna          8 ?????????????????????????????????????????       49
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[95160479, 95160520],
                          [       8,       49]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 41 chr8 146364022 + 95160479 95160520 hg19_dna 50 + 8 49 2
41

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr22")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 51304566)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr22      42144400 ???????????????????????????????????? 42144436
                  0 ||||||||||||||||||||||||||||||||||||       36
hg19_dna         11 ????????????????????????????????????       47
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[42144400, 42144436],
                          [      11,       47]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 30 chr22 51304566 + 42144400 42144436 hg19_dna 50 + 11 47 3
36

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 48))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr2      183925984 ??????----?????????????????????????????????????? 183926028
                  0 ||||||----||||||||||||||||||||||||||||||||||||||        48
hg19_dna          1 ????????????????????????????????????????????????        49
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[183925984, 183925990, 183925990, 183926028],
                          [        1,         7,        11,        49]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 41 chr2 243199373 + 183925984 183926028 hg19_dna 50 + 1 49 4
6	0	4
38

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 170))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr19      35483340 ????????????????????????????????????????????????????????????
                  0 |||||||||||||||||||||||||-----------------------------------
hg19_dna         10 ?????????????????????????-----------------------------------

chr19      35483400 ????????????????????????????????????????????????????????????
                 60 ------------------------------------------------------------
hg19_dna         35 ------------------------------------------------------------

chr19      35483460 ?????????????????????????????????????????????????? 35483510
                120 ---------------------------------------|||||||||||      170
hg19_dna         35 ---------------------------------------???????????       46
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[35483340, 35483365, 35483499, 35483510],
                          [      10,       35,       35,       46]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 31 chr19 59128983 + 35483340 35483510 hg19_dna 50 + 10 46 5
25	134	0
11

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr18")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 78077248)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr18      23891310 ??????????????????????????????????????? 23891349
                  0 |||||||||||||||||||||||||||||||||||||||       39
hg19_dna         10 ???????????????????????????????????????       49
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[23891310, 23891349],
                          [      10,       49]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 39 chr18 78077248 + 23891310 23891349 hg19_dna 50 + 10 49 6
39

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 28))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr18")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 78077248)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr18      43252217 ???????????????????????????? 43252245
                  0 ||||||||||||||||||||||||||||       28
hg19_dna         21 ????????????????????????????       49
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[43252217, 43252245],
                          [      21,       49]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 26 chr18 78077248 + 43252217 43252245 hg19_dna 50 + 21 49 7
28

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 54))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 115169878)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr13      52759147 ?????????????---??????????????????????????????????????
                  0 |||||||---------||||||||||||||||||||||||||||||||||||||
hg19_dna          1 ???????------?????????????????????????????????????????

chr13      52759198
                 54
hg19_dna         49
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[52759147, 52759154, 52759160, 52759160, 52759198],
                          [       1,        8,        8,       11,       49]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 41 chr13 115169878 + 52759147 52759198 hg19_dna 50 + 1 49 8
7	6	3
38

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 50))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr1        1207056 ?????????????????????????????????????????????????? 1207106
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||      50
hg19_dna          0 ??????????????????????????????????????????????????      50
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[1207056, 1207106],
                          [      0,      50]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 50 chr1 249250621 + 1207056 1207106 hg19_dna 50 + 0 50 9
50

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr1       61700837 ?????????????????????????????????? 61700871
                  0 ||||||||||||||||||||||||||||||||||       34
hg19_dna          1 ??????????????????????????????????       35
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[61700837, 61700871],
                          [       1,       35]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 28 chr1 249250621 + 61700837 61700871 hg19_dna 50 + 1 35 10
34

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 44))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 191154276)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr4       37558157 ????????????????----------?????????????????? 37558191
                  0 ||||||||||----------------||||||||||||||||||       44
hg19_dna         49 ??????????------????????????????????????????       11
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[37558157, 37558167, 37558173, 37558173, 37558191],
                          [      49,       39,       39,       29,       11]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 26 chr4 191154276 + 37558157 37558191 hg19_dna 50 - 1 39 11
10	6	10
18

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 37))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr22")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 51304566)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr22      48997405 ????????????????????????????????????? 48997442
                  0 |||||||||||||||||||||||||||||||||||||       37
hg19_dna         49 ?????????????????????????????????????       12
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[48997405, 48997442],
                          [      49,       12]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 33 chr22 51304566 + 48997405 48997442 hg19_dna 50 - 1 38 12
37

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr2      120641740 ???????????????????????????????????? 120641776
                  0 ||||||||||||||||||||||||||||||||||||        36
hg19_dna         49 ????????????????????????????????????        13
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[120641740, 120641776],
                          [       49,        13]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 34 chr2 243199373 + 120641740 120641776 hg19_dna 50 - 1 37 13
36

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr19      54017130 ??????????????????????????????????????? 54017169
                  0 |||||||||||||||||||||||||||||||||||||||       39
hg19_dna         49 ???????????????????????????????????????       10
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[54017130, 54017169],
                          [      49,       10]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 39 chr19 59128983 + 54017130 54017169 hg19_dna 50 - 1 40 14
39

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr19        553742 ??????????????????????????????????????? 553781
                  0 |||||||||||||||||||||||||||||||||||||||     39
hg19_dna         49 ???????????????????????????????????????     10
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[553742, 553781],
                          [    49,     10]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 33 chr19 59128983 + 553742 553781 hg19_dna 50 - 1 40 15
39

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr10")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 135534747)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr10      99388555 ???????????????????????????????????? 99388591
                  0 ||||||||||||||||||||||||||||||||||||       36
hg19_dna         49 ????????????????????????????????????       13
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[99388555, 99388591],
                          [      49,       13]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 30 chr10 135534747 + 99388555 99388591 hg19_dna 50 - 1 37 16
36

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 25))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr10")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 135534747)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr10     112178171 ????????????????????????? 112178196
                  0 |||||||||||||||||||||||||        25
hg19_dna         35 ?????????????????????????        10
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[112178171, 112178196],
                          [       35,        10]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 23 chr10 135534747 + 112178171 112178196 hg19_dna 50 - 15 40 17
25

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr1       39368490 ???????????????????????????????????? 39368526
                  0 ||||||||||||||||||||||||||||||||||||       36
hg19_dna         49 ????????????????????????????????????       13
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[39368490, 39368526],
                          [      49,       13]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 34 chr1 249250621 + 39368490 39368526 hg19_dna 50 - 1 37 18
36

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr1      220325687 ?????????????????????????????????? 220325721
                  0 ||||||||||||||||||||||||||||||||||        34
hg19_dna         47 ??????????????????????????????????        13
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[220325687, 220325721],
                          [       47,        13]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 32 chr1 249250621 + 220325687 220325721 hg19_dna 50 - 3 37 19
34

""",
        )
        self.assertRaises(StopIteration, next, alignments)

    def test_writing_psl_34_004(self):
        """Test writing the alignments in psl_34_004.chain."""
        path = "Blat/psl_34_004.chain"
        with open(path) as stream:
            original_data = stream.read()
        alignments = Align.parse(path, "chain")
        stream = StringIO()
        n = Align.write(alignments, stream, "chain")
        self.assertEqual(n, 19)
        stream.seek(0)
        alignments = Align.parse(stream, "chain")
        self.check_reading_psl_34_004(alignments)

    def test_reading_psl_34_005(self):
        """Test parsing psl_34_005.chain."""
        # The chain file psl_34_005.chain was generated from the PSL file using:
        # pslToChain psl_34_005.psl psl_34_005.chain
        path = "Blat/psl_34_005.chain"
        alignments = Align.parse(path, "chain")
        self.check_reading_psl_34_005(alignments)

    def check_reading_psl_34_005(self, alignments):
        """Check parsing psl_34_005.chain."""
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 16))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertEqual(len(alignment.target.seq), 191154276)
        self.assertEqual(len(alignment.query.seq), 33)
        self.assertEqual(
            str(alignment),
            """\
chr4       61646095 ???????????????? 61646111
                  0 ||||||||||||||||       16
hg18_dna         11 ????????????????       27
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[61646095, 61646111],
                          [      11,       27]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 16 chr4 191154276 + 61646095 61646111 hg18_dna 33 + 11 27 1
16

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 33))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 33)
        self.assertEqual(
            str(alignment),
            """\
chr1       10271783 ????????????????????????????????? 10271816
                  0 |||||||||||||||||||||||||||||||||       33
hg18_dna          0 ?????????????????????????????????       33
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[10271783, 10271816],
                          [       0,       33]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 33 chr1 249250621 + 10271783 10271816 hg18_dna 33 + 0 33 2
33

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 17))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 33)
        self.assertEqual(
            str(alignment),
            """\
chr2       53575980 ????????????????? 53575997
                  0 |||||||||||||||||       17
hg18_dna         25 ?????????????????        8
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[53575980, 53575997],
                          [      25,        8]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 17 chr2 243199373 + 53575980 53575997 hg18_dna 33 - 8 25 3
17

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 41))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr9")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 141213431)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr9       85737865 ????????????????????????????????????????? 85737906
                  0 |||||||||||||||||||||||||||||||||||||||||       41
hg19_dna          9 ?????????????????????????????????????????       50
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[85737865, 85737906],
                          [       9,       50]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 35 chr9 141213431 + 85737865 85737906 hg19_dna 50 + 9 50 4
41

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 41))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr8")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 146364022)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr8       95160479 ????????????????????????????????????????? 95160520
                  0 |||||||||||||||||||||||||||||||||||||||||       41
hg19_dna          8 ?????????????????????????????????????????       49
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[95160479, 95160520],
                          [       8,       49]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 41 chr8 146364022 + 95160479 95160520 hg19_dna 50 + 8 49 5
41

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr22")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 51304566)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr22      42144400 ???????????????????????????????????? 42144436
                  0 ||||||||||||||||||||||||||||||||||||       36
hg19_dna         11 ????????????????????????????????????       47
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[42144400, 42144436],
                          [      11,       47]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 30 chr22 51304566 + 42144400 42144436 hg19_dna 50 + 11 47 6
36

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 48))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr2      183925984 ??????----?????????????????????????????????????? 183926028
                  0 ||||||----||||||||||||||||||||||||||||||||||||||        48
hg19_dna          1 ????????????????????????????????????????????????        49
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[183925984, 183925990, 183925990, 183926028],
                          [        1,         7,        11,        49]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 41 chr2 243199373 + 183925984 183926028 hg19_dna 50 + 1 49 7
6	0	4
38

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 170))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr19      35483340 ????????????????????????????????????????????????????????????
                  0 |||||||||||||||||||||||||-----------------------------------
hg19_dna         10 ?????????????????????????-----------------------------------

chr19      35483400 ????????????????????????????????????????????????????????????
                 60 ------------------------------------------------------------
hg19_dna         35 ------------------------------------------------------------

chr19      35483460 ?????????????????????????????????????????????????? 35483510
                120 ---------------------------------------|||||||||||      170
hg19_dna         35 ---------------------------------------???????????       46
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[35483340, 35483365, 35483499, 35483510],
                          [      10,       35,       35,       46]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 31 chr19 59128983 + 35483340 35483510 hg19_dna 50 + 10 46 8
25	134	0
11

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr18")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 78077248)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr18      23891310 ??????????????????????????????????????? 23891349
                  0 |||||||||||||||||||||||||||||||||||||||       39
hg19_dna         10 ???????????????????????????????????????       49
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[23891310, 23891349],
                          [      10,       49]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 39 chr18 78077248 + 23891310 23891349 hg19_dna 50 + 10 49 9
39

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 28))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr18")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 78077248)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr18      43252217 ???????????????????????????? 43252245
                  0 ||||||||||||||||||||||||||||       28
hg19_dna         21 ????????????????????????????       49
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[43252217, 43252245],
                          [      21,       49]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 26 chr18 78077248 + 43252217 43252245 hg19_dna 50 + 21 49 10
28

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 54))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 115169878)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr13      52759147 ?????????????---??????????????????????????????????????
                  0 |||||||---------||||||||||||||||||||||||||||||||||||||
hg19_dna          1 ???????------?????????????????????????????????????????

chr13      52759198
                 54
hg19_dna         49
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[52759147, 52759154, 52759160, 52759160, 52759198],
                          [       1,        8,        8,       11,       49]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 41 chr13 115169878 + 52759147 52759198 hg19_dna 50 + 1 49 11
7	6	3
38

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 50))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr1        1207056 ?????????????????????????????????????????????????? 1207106
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||      50
hg19_dna          0 ??????????????????????????????????????????????????      50
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[1207056, 1207106],
                          [      0,      50]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 50 chr1 249250621 + 1207056 1207106 hg19_dna 50 + 0 50 12
50

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr1       61700837 ?????????????????????????????????? 61700871
                  0 ||||||||||||||||||||||||||||||||||       34
hg19_dna          1 ??????????????????????????????????       35
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[61700837, 61700871],
                          [       1,       35]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 28 chr1 249250621 + 61700837 61700871 hg19_dna 50 + 1 35 13
34

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 44))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 191154276)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr4       37558157 ????????????????----------?????????????????? 37558191
                  0 ||||||||||----------------||||||||||||||||||       44
hg19_dna         49 ??????????------????????????????????????????       11
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[37558157, 37558167, 37558173, 37558173, 37558191],
                          [      49,       39,       39,       29,       11]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 26 chr4 191154276 + 37558157 37558191 hg19_dna 50 - 1 39 14
10	6	10
18

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 37))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr22")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 51304566)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr22      48997405 ????????????????????????????????????? 48997442
                  0 |||||||||||||||||||||||||||||||||||||       37
hg19_dna         49 ?????????????????????????????????????       12
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[48997405, 48997442],
                          [      49,       12]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 33 chr22 51304566 + 48997405 48997442 hg19_dna 50 - 1 38 15
37

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr2      120641740 ???????????????????????????????????? 120641776
                  0 ||||||||||||||||||||||||||||||||||||        36
hg19_dna         49 ????????????????????????????????????        13
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[120641740, 120641776],
                          [       49,        13]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 34 chr2 243199373 + 120641740 120641776 hg19_dna 50 - 1 37 16
36

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr19      54017130 ??????????????????????????????????????? 54017169
                  0 |||||||||||||||||||||||||||||||||||||||       39
hg19_dna         49 ???????????????????????????????????????       10
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[54017130, 54017169],
                          [      49,       10]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 39 chr19 59128983 + 54017130 54017169 hg19_dna 50 - 1 40 17
39

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr19        553742 ??????????????????????????????????????? 553781
                  0 |||||||||||||||||||||||||||||||||||||||     39
hg19_dna         49 ???????????????????????????????????????     10
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[553742, 553781],
                          [    49,     10]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 33 chr19 59128983 + 553742 553781 hg19_dna 50 - 1 40 18
39

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr10")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 135534747)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr10      99388555 ???????????????????????????????????? 99388591
                  0 ||||||||||||||||||||||||||||||||||||       36
hg19_dna         49 ????????????????????????????????????       13
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[99388555, 99388591],
                          [      49,       13]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 30 chr10 135534747 + 99388555 99388591 hg19_dna 50 - 1 37 19
36

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 25))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr10")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 135534747)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr10     112178171 ????????????????????????? 112178196
                  0 |||||||||||||||||||||||||        25
hg19_dna         35 ?????????????????????????        10
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[112178171, 112178196],
                          [       35,        10]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 23 chr10 135534747 + 112178171 112178196 hg19_dna 50 - 15 40 20
25

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr1       39368490 ???????????????????????????????????? 39368526
                  0 ||||||||||||||||||||||||||||||||||||       36
hg19_dna         49 ????????????????????????????????????       13
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[39368490, 39368526],
                             [      49,       13]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 34 chr1 249250621 + 39368490 39368526 hg19_dna 50 - 1 37 21
36

""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertEqual(
            str(alignment),
            """\
chr1      220325687 ?????????????????????????????????? 220325721
                  0 ||||||||||||||||||||||||||||||||||        34
hg19_dna         47 ??????????????????????????????????        13
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[220325687, 220325721],
                          [       47,        13]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "chain"),
            """\
chain 32 chr1 249250621 + 220325687 220325721 hg19_dna 50 - 3 37 22
34

""",
        )
        self.assertRaises(StopIteration, next, alignments)

    def test_writing_psl_34_005(self):
        """Test writing the alignments in psl_34_005.chain."""
        path = "Blat/psl_34_005.chain"
        alignments = Align.parse(path, "chain")
        stream = StringIO()
        n = Align.write(alignments, stream, "chain")
        self.assertEqual(n, 22)
        stream.seek(0)
        alignments = Align.parse(stream, "chain")
        self.check_reading_psl_34_005(alignments)


class TestAlign_strand(unittest.TestCase):
    def test_format(self):
        """Test alignment with the target on the opposite strand."""
        sequences = ["AACAGCAGCGTGTCG", "CAGCTAGCGAA"]
        coordinates = np.array(
            [[0, 2, 2, 3, 4, 6, 6, 9, 10, 12, 15], [11, 11, 9, 8, 8, 6, 5, 2, 2, 0, 0]]
        )
        score = 8
        alignment = Alignment(sequences, coordinates)
        alignment.score = score
        chain1 = """\
chain 8 target 15 + 0 15 query 11 - 0 11
0	2	2
1	1	0
2	0	1
3	1	0
0	3	0
2

"""
        chain2 = """\
chain 8 target 15 - 0 15 query 11 + 0 11
2	3	0
3	1	0
2	0	1
1	1	0
0	2	2
0

"""
        chain3 = """\
chain 8 target 15 + 0 15 query 11 - 0 11
2	3	0
3	1	0
2	0	1
1	1	0
0	2	2
0

"""
        chain4 = """\
chain 8 target 15 - 0 15 query 11 + 0 11
0	2	2
1	1	0
2	0	1
3	1	0
0	3	0
2

"""
        self.assertEqual(
            str(alignment),
            """\
target            0 AA--CAGC-AGCGTGTCG 15
                  0 ----|-||-|||-||--- 18
query            11 --TTC-GCTAGC-TG---  0
""",
        )
        self.assertEqual(format(alignment, "chain"), chain1)
        alignment.coordinates = alignment.coordinates[:, ::-1]
        self.assertEqual(
            str(alignment),
            """\
target           15 CGACACGCT-GCTG--TT  0
                  0 ---||-|||-||-|---- 18
query             0 ---CA-GCTAGC-GAA-- 11
""",
        )
        self.assertEqual(format(alignment, "chain"), chain2)
        alignment.coordinates = alignment.coordinates[:, ::-1]
        alignment = alignment.reverse_complement()
        alignment.score = score
        self.assertEqual(
            str(alignment),
            """\
target            0 CGACACGCT-GCTG--TT 15
                  0 ---||-|||-||-|---- 18
query            11 ---CA-GCTAGC-GAA--  0
""",
        )
        self.assertEqual(format(alignment, "chain"), chain3)
        alignment.coordinates = alignment.coordinates[:, ::-1]
        self.assertEqual(
            str(alignment),
            """\
target           15 AA--CAGC-AGCGTGTCG  0
                  0 ----|-||-|||-||--- 18
query             0 --TTC-GCTAGC-TG--- 11
""",
        )
        self.assertEqual(format(alignment, "chain"), chain4)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
