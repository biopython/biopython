# Copyright 2022 by Michiel de Hoon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Align.psl module."""
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
        "Install numpy if you want to use Bio.Align.psl."
    ) from None


class TestAlign_dna_rna(unittest.TestCase):
    # The PSL file dna_rna.psl was generated using this command:
    # blat -mask=lower hg38.2bit rna.fa dna_rna.unsorted.psl
    # pslSort dirs dna_rna.psl . dna_rna.unsorted.psl

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
        """Test parsing dna_rna.psl."""
        path = "Blat/dna_rna.psl"
        alignments = Align.parse(path, "psl")
        self.check_alignments(alignments)
        alignments.rewind()
        self.check_alignments(alignments)
        with Align.parse(path, "psl") as alignments:
            self.check_alignments(alignments)
        with self.assertRaises(AttributeError):
            alignments._stream
        with Align.parse(path, "psl") as alignments:
            pass
        with self.assertRaises(AttributeError):
            alignments._stream

    def check_alignments(self, alignments):
        self.assertEqual(alignments.metadata["psLayout version"], "3")
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 175)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 6)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 1711))
        self.assertEqual(len(alignment), 2)
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
        matches = sum(
            alignment.substitutions[c, c] for c in alignment.substitutions.alphabet
        )
        repMatches = sum(
            alignment.substitutions[c, c.swapcase()]
            for c in alignment.substitutions.alphabet
        )
        self.assertEqual(matches, alignment.matches)
        self.assertEqual(repMatches, alignment.repMatches)
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
            format(alignment, "psl"),
            """\
175	0	6	0	0	0	2	1530	-	NR_046654.1	181	0	181	chr3	198295559	42530895	42532606	3	63,75,43,	0,63,138,	42530895,42532020,42532563,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 172)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 6)
        self.assertEqual(alignment.nCount, 0)
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
        matches = sum(
            alignment.substitutions[c, c] for c in alignment.substitutions.alphabet
        )
        repMatches = sum(
            alignment.substitutions[c, c.swapcase()]
            for c in alignment.substitutions.alphabet
            if c != "X"
        )
        self.assertEqual(matches, alignment.matches)
        self.assertEqual(repMatches, alignment.repMatches)
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
            format(alignment, "psl"),
            """\
172	1	6	0	1	3	3	1532	-	NR_046654.1_modified	190	3	185	chr3	198295559	42530895	42532606	5	27,36,17,56,43,	5,35,71,88,144,	42530895,42530922,42532020,42532039,42532563,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 165)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 39)
        self.assertEqual(alignment.nCount, 0)
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
        matches = sum(
            alignment.substitutions[c, c] for c in alignment.substitutions.alphabet
        )
        repMatches = sum(
            alignment.substitutions[c, c.swapcase()]
            for c in alignment.substitutions.alphabet
        )
        self.assertEqual(matches, alignment.matches)
        self.assertEqual(repMatches, alignment.repMatches)
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
            format(alignment, "psl"),
            """\
165	0	39	0	0	0	2	5203	+	NR_111921.1	216	0	204	chr3	198295559	48663767	48669174	3	46,82,76,	0,46,128,	48663767,48665640,48669098,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 162)
        self.assertEqual(alignment.misMatches, 2)
        self.assertEqual(alignment.repMatches, 39)
        self.assertEqual(alignment.nCount, 0)
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
        matches = sum(
            alignment.substitutions[c, c] for c in alignment.substitutions.alphabet
        )
        repMatches = sum(
            alignment.substitutions[c, c.swapcase()]
            for c in alignment.substitutions.alphabet
            if c != "X"
        )
        self.assertEqual(matches, alignment.matches)
        self.assertEqual(repMatches, alignment.repMatches)
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
            format(alignment, "psl"),
            """\
162	2	39	0	1	2	3	5204	+	NR_111921.1_modified	220	3	208	chr3	198295559	48663767	48669174	5	28,17,76,6,76,	3,31,48,126,132,	48663767,48663796,48665640,48665716,48669098,
""",
        )
        self.assertRaises(StopIteration, next, alignments)

    def test_writing(self):
        """Test writing the alignments in dna_rna.psl."""
        path = "Blat/dna_rna.psl"
        with open(path) as stream:
            original_data = stream.read()
        alignments = Align.parse(path, "psl")
        stream = StringIO()
        n = Align.write(alignments, stream, "psl")
        self.assertEqual(n, 4)
        stream.seek(0)
        written_data = stream.read()
        stream.close()
        self.assertEqual(original_data, written_data)
        # Try this again. This time, we first strip the matches, misMatches,
        # repMatches, and nCount attributes from each alignment, and insert the
        # appropriate sequence data in each alignment. The writer will then
        # recalculate the matches, misMatches, repMatches, and nCount values
        # from the sequence data and the alignment, and store those values in
        # the PSL file.
        alignments = []
        for alignment in Align.parse(path, "psl"):
            del alignment.matches
            del alignment.misMatches
            del alignment.repMatches
            del alignment.nCount
            dna = Seq(self.dna, length=len(alignment.target))
            alignment.target.seq = dna
            alignment.query.seq = self.rna[alignment.sequences[1].id]
            alignments.append(alignment)
        stream = StringIO()
        n = Align.write(alignments, stream, "psl", mask="lower")
        self.assertEqual(n, 4)
        stream.seek(0)
        written_data = stream.read()
        stream.close()
        self.assertEqual(original_data, written_data)


class TestAlign_dna(unittest.TestCase):
    queries = {
        record.id: str(record.seq)
        for record in SeqIO.parse("Blat/fasta_34.fa", "fasta")
    }

    def test_reading_psl_34_001(self):
        """Test parsing psl_34_001.psl and pslx_34_001.pslx."""
        self.check_reading_psl_34_001("psl")
        self.check_reading_psl_34_001("pslx")

    def check_reading_psl_34_001(self, fmt):
        """Check parsing psl_34_001.psl or pslx_34_001.pslx."""
        path = "Blat/%s_34_001.%s" % (fmt, fmt)
        alignments = Align.parse(path, "psl")
        self.assertEqual(alignments.metadata["psLayout version"], "3")
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 16)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 16))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertEqual(len(alignment.target.seq), 191154276)
        self.assertEqual(len(alignment.query.seq), 33)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr4       61646095 ???????????????? 61646111
                  0 ||||||||||||||||       16
hg18_dna         11 ????????????????       27
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[61646095:61646111], "aggtaaactgccttca"
            )
            self.assertEqual(
                alignment.query.seq[11:27], self.queries[alignment.query.id][11:27]
            )
            self.assertEqual(
                str(alignment),
                """\
chr4       61646095 aggtaaactgccttca 61646111
                  0 ||||||||||||||||       16
hg18_dna         11 aggtaaactgccttca       27
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['a', 'g', 'g', 't', 'a', 'a', 'a', 'c', 't', 'g', 'c', 'c', 't',
           't', 'c', 'a'],
          ['a', 'g', 'g', 't', 'a', 'a', 'a', 'c', 't', 'g', 'c', 'c', 't',
           't', 'c', 'a']], dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
16	0	0	0	0	0	0	0	+	hg18_dna	33	11	27	chr4	191154276	61646095	61646111	1	16,	11,	61646095,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 33)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 33))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 33)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr1       10271783 ????????????????????????????????? 10271816
                  0 |||||||||||||||||||||||||||||||||       33
hg18_dna          0 ?????????????????????????????????       33
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[10271783:10271816],
                "atgagcttccaaggtaaactgccttcaagattc",
            )
            self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
            self.assertEqual(
                str(alignment),
                """\
chr1       10271783 atgagcttccaaggtaaactgccttcaagattc 10271816
                  0 |||||||||||||||||||||||||||||||||       33
hg18_dna          0 atgagcttccaaggtaaactgccttcaagattc       33
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['a', 't', 'g', 'a', 'g', 'c', 't', 't', 'c', 'c', 'a', 'a', 'g',
           'g', 't', 'a', 'a', 'a', 'c', 't', 'g', 'c', 'c', 't', 't', 'c',
           'a', 'a', 'g', 'a', 't', 't', 'c'],
          ['a', 't', 'g', 'a', 'g', 'c', 't', 't', 'c', 'c', 'a', 'a', 'g',
           'g', 't', 'a', 'a', 'a', 'c', 't', 'g', 'c', 'c', 't', 't', 'c',
           'a', 'a', 'g', 'a', 't', 't', 'c']], dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
33	0	0	0	0	0	0	0	+	hg18_dna	33	0	33	chr1	249250621	10271783	10271816	1	33,	0,	10271783,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 17)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 17))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 33)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr2       53575980 ????????????????? 53575997
                  0 |||||||||||||||||       17
hg18_dna         25 ?????????????????        8
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[53575980:53575997], "aaggcagtttaccttgg"
            )
            self.assertEqual(
                alignment.query.seq[8:25], self.queries[alignment.query.id][8:25]
            )
            self.assertEqual(
                str(alignment),
                """\
chr2       53575980 aaggcagtttaccttgg 53575997
                  0 |||||||||||||||||       17
hg18_dna         25 aaggcagtttaccttgg        8
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['a', 'a', 'g', 'g', 'c', 'a', 'g', 't', 't', 't', 'a', 'c', 'c',
           't', 't', 'g', 'g'],
          ['a', 'a', 'g', 'g', 'c', 'a', 'g', 't', 't', 't', 'a', 'c', 'c',
           't', 't', 'g', 'g']], dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
17	0	0	0	0	0	0	0	-	hg18_dna	33	8	25	chr2	243199373	53575980	53575997	1	17,	8,	53575980,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 38)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 41))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr9")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 141213431)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr9       85737865 ????????????????????????????????????????? 85737906
                  0 |||||||||||||||||||||||||||||||||||||||||       41
hg19_dna          9 ?????????????????????????????????????????       50
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[85737865:85737906],
                "acaaaggggctgggcgcagtggctcacgcctgtaatcccaa",
            )
            self.assertEqual(
                alignment.query.seq[9:50], self.queries[alignment.query.id][9:50]
            )
            self.assertEqual(
                str(alignment),
                """\
chr9       85737865 acaaaggggctgggcgcagtggctcacgcctgtaatcccaa 85737906
                  0 ||||||||||||||||..|||||||||.|||||||||||||       41
hg19_dna          9 acaaaggggctgggcgtggtggctcacacctgtaatcccaa       50
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['a', 'c', 'a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g',
           'g', 'c', 'g', 'c', 'a', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a',
           'c', 'g', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c', 'c',
           'a', 'a'],
          ['a', 'c', 'a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g',
           'g', 'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a',
           'c', 'a', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c', 'c',
           'a', 'a']], dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
38	3	0	0	0	0	0	0	+	hg19_dna	50	9	50	chr9	141213431	85737865	85737906	1	41,	9,	85737865,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 41)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 41))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr8")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 146364022)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr8       95160479 ????????????????????????????????????????? 95160520
                  0 |||||||||||||||||||||||||||||||||||||||||       41
hg19_dna          8 ?????????????????????????????????????????       49
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[95160479:95160520],
                "cacaaaggggctgggcgtggtggctcacacctgtaatccca",
            )
            self.assertEqual(
                alignment.query.seq[8:49], self.queries[alignment.query.id][8:49]
            )
            self.assertEqual(
                str(alignment),
                """\
chr8       95160479 cacaaaggggctgggcgtggtggctcacacctgtaatccca 95160520
                  0 |||||||||||||||||||||||||||||||||||||||||       41
hg19_dna          8 cacaaaggggctgggcgtggtggctcacacctgtaatccca       49
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['c', 'a', 'c', 'a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g',
           'g', 'g', 'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c',
           'a', 'c', 'a', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c',
           'c', 'a'],
          ['c', 'a', 'c', 'a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g',
           'g', 'g', 'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c',
           'a', 'c', 'a', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c',
           'c', 'a']], dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
41	0	0	0	0	0	0	0	+	hg19_dna	50	8	49	chr8	146364022	95160479	95160520	1	41,	8,	95160479,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 33)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr22")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 51304566)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr22      42144400 ???????????????????????????????????? 42144436
                  0 ||||||||||||||||||||||||||||||||||||       36
hg19_dna         11 ????????????????????????????????????       47
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[42144400:42144436],
                "aaaggggctgggcgtggtagctcatgcctgtaatcc",
            )
            self.assertEqual(
                alignment.query.seq[11:47], self.queries[alignment.query.id][11:47]
            )
            self.assertEqual(
                str(alignment),
                """\
chr22      42144400 aaaggggctgggcgtggtagctcatgcctgtaatcc 42144436
                  0 ||||||||||||||||||.|||||..||||||||||       36
hg19_dna         11 aaaggggctgggcgtggtggctcacacctgtaatcc       47
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g', 'c',
           'g', 't', 'g', 'g', 't', 'a', 'g', 'c', 't', 'c', 'a', 't', 'g',
           'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c'],
          ['a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g', 'c',
           'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a', 'c', 'a',
           'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c']], dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
33	3	0	0	0	0	0	0	+	hg19_dna	50	11	47	chr22	51304566	42144400	42144436	1	36,	11,	42144400,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 43)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 48))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr2      183925984 ??????----?????????????????????????????????????? 183926028
                  0 ||||||----||||||||||||||||||||||||||||||||||||||        48
hg19_dna          1 ????????????????????????????????????????????????        49
""",
            )
        elif fmt == "pslx":
            self.assertEqual(alignment.target.seq[183925984:183925990], "aaaaat")
            self.assertEqual(
                alignment.target.seq[183925990:183926028],
                "aaaggggctgggcgtggtggctcacgcctgtaatccca",
            )
            self.assertEqual(
                alignment.query.seq[1:7], self.queries[alignment.query.id][1:7]
            )
            self.assertEqual(
                alignment.query.seq[38:49], self.queries[alignment.query.id][38:49]
            )
            self.assertEqual(
                str(alignment),
                """\
chr2      183925984 aaaaat----aaaggggctgggcgtggtggctcacgcctgtaatccca 183926028
                  0 ||||||----|||||||||||||||||||||||||.||||||||||||        48
hg19_dna          1 aaaaat????aaaggggctgggcgtggtggctcacacctgtaatccca        49
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
            format(alignment, "psl"),
            """\
43	1	0	0	1	4	0	0	+	hg19_dna	50	1	49	chr2	243199373	183925984	183926028	2	6,38,	1,11,	183925984,183925990,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 34)
        self.assertEqual(alignment.misMatches, 2)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 170))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
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
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[35483340:35483365], "caaaggggctgggcgtagtggctga"
            )
            self.assertEqual(alignment.target.seq[35483499:35483510], "cacctgtaatc")
            self.assertEqual(
                alignment.query.seq[10:46], self.queries[alignment.query.id][10:46]
            )
            self.assertEqual(
                str(alignment),
                """\
chr19      35483340 caaaggggctgggcgtagtggctga???????????????????????????????????
                  0 ||||||||||||||||.||||||.|-----------------------------------
hg19_dna         10 caaaggggctgggcgtggtggctca-----------------------------------

chr19      35483400 ????????????????????????????????????????????????????????????
                 60 ------------------------------------------------------------
hg19_dna         35 ------------------------------------------------------------

chr19      35483460 ???????????????????????????????????????cacctgtaatc 35483510
                120 ---------------------------------------|||||||||||      170
hg19_dna         35 ---------------------------------------cacctgtaatc       46
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
            format(alignment, "psl"),
            """\
34	2	0	0	0	0	1	134	+	hg19_dna	50	10	46	chr19	59128983	35483340	35483510	2	25,11,	10,35,	35483340,35483499,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 39)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr18")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 78077248)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr18      23891310 ??????????????????????????????????????? 23891349
                  0 |||||||||||||||||||||||||||||||||||||||       39
hg19_dna         10 ???????????????????????????????????????       49
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[23891310:23891349],
                "caaaggggctgggcgtggtggctcacacctgtaatccca",
            )
            self.assertEqual(
                alignment.query.seq[10:49], self.queries[alignment.query.id][10:49]
            )
            self.assertEqual(
                str(alignment),
                """\
chr18      23891310 caaaggggctgggcgtggtggctcacacctgtaatccca 23891349
                  0 |||||||||||||||||||||||||||||||||||||||       39
hg19_dna         10 caaaggggctgggcgtggtggctcacacctgtaatccca       49
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['c', 'a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g',
           'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a', 'c',
           'a', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c', 'c', 'a'],
          ['c', 'a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g',
           'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a', 'c',
           'a', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c', 'c', 'a']],
         dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
39	0	0	0	0	0	0	0	+	hg19_dna	50	10	49	chr18	78077248	23891310	23891349	1	39,	10,	23891310,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 27)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 28))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr18")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 78077248)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr18      43252217 ???????????????????????????? 43252245
                  0 ||||||||||||||||||||||||||||       28
hg19_dna         21 ????????????????????????????       49
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[43252217:43252245], "ggcgtggtggctcacgcctgtaatccca"
            )
            self.assertEqual(
                alignment.query.seq[21:49], self.queries[alignment.query.id][21:49]
            )
            self.assertEqual(
                str(alignment),
                """\
chr18      43252217 ggcgtggtggctcacgcctgtaatccca 43252245
                  0 |||||||||||||||.||||||||||||       28
hg19_dna         21 ggcgtggtggctcacacctgtaatccca       49
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['g', 'g', 'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c',
           'a', 'c', 'g', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c',
           'c', 'a'],
          ['g', 'g', 'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c',
           'a', 'c', 'a', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c',
           'c', 'a']], dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
27	1	0	0	0	0	0	0	+	hg19_dna	50	21	49	chr18	78077248	43252217	43252245	1	28,	21,	43252217,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 44)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 54))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 115169878)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
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
        elif fmt == "pslx":
            self.assertEqual(alignment.target.seq[52759147:52759154], "aaaaatt")
            self.assertEqual(
                alignment.target.seq[52759160:52759198],
                "aaaggggctgggcgtggtggctcacgcctgtaatccca",
            )
            self.assertEqual(
                alignment.query.seq[1:8], self.queries[alignment.query.id][1:8]
            )
            self.assertEqual(
                alignment.query.seq[38:49], self.queries[alignment.query.id][38:49]
            )
            self.assertEqual(
                str(alignment),
                """\
chr13      52759147 aaaaatt??????---aaaggggctgggcgtggtggctcacgcctgtaatccca
                  0 |||||||---------|||||||||||||||||||||||||.||||||||||||
hg19_dna          1 aaaaatt------???aaaggggctgggcgtggtggctcacacctgtaatccca

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
            format(alignment, "psl"),
            """\
44	1	0	0	1	3	1	6	+	hg19_dna	50	1	49	chr13	115169878	52759147	52759198	2	7,38,	1,11,	52759147,52759160,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 50)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 50))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr1        1207056 ?????????????????????????????????????????????????? 1207106
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||      50
hg19_dna          0 ??????????????????????????????????????????????????      50
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[1207056:1207106],
                "caaaaattcacaaaggggctgggcgtggtggctcacacctgtaatcccaa",
            )
            self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
            self.assertEqual(
                str(alignment),
                """\
chr1        1207056 caaaaattcacaaaggggctgggcgtggtggctcacacctgtaatcccaa 1207106
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||      50
hg19_dna          0 caaaaattcacaaaggggctgggcgtggtggctcacacctgtaatcccaa      50
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['c', 'a', 'a', 'a', 'a', 'a', 't', 't', 'c', 'a', 'c', 'a', 'a',
           'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g', 'c', 'g', 't',
           'g', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a', 'c', 'a', 'c', 'c',
           't', 'g', 't', 'a', 'a', 't', 'c', 'c', 'c', 'a', 'a'],
          ['c', 'a', 'a', 'a', 'a', 'a', 't', 't', 'c', 'a', 'c', 'a', 'a',
           'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g', 'c', 'g', 't',
           'g', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a', 'c', 'a', 'c', 'c',
           't', 'g', 't', 'a', 'a', 't', 'c', 'c', 'c', 'a', 'a']],
         dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
50	0	0	0	0	0	0	0	+	hg19_dna	50	0	50	chr1	249250621	1207056	1207106	1	50,	0,	1207056,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 31)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr1       61700837 ?????????????????????????????????? 61700871
                  0 ||||||||||||||||||||||||||||||||||       34
hg19_dna          1 ??????????????????????????????????       35
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[61700837:61700871],
                "aaaaatgaacaaaggggctgggcgcggtggctca",
            )
            self.assertEqual(
                alignment.query.seq[1:35], self.queries[alignment.query.id][1:35]
            )
            self.assertEqual(
                str(alignment),
                """\
chr1       61700837 aaaaatgaacaaaggggctgggcgcggtggctca 61700871
                  0 ||||||..||||||||||||||||.|||||||||       34
hg19_dna          1 aaaaattcacaaaggggctgggcgtggtggctca       35
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['a', 'a', 'a', 'a', 'a', 't', 'g', 'a', 'a', 'c', 'a', 'a', 'a',
           'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g', 'c', 'g', 'c', 'g',
           'g', 't', 'g', 'g', 'c', 't', 'c', 'a'],
          ['a', 'a', 'a', 'a', 'a', 't', 't', 'c', 'a', 'c', 'a', 'a', 'a',
           'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g', 'c', 'g', 't', 'g',
           'g', 't', 'g', 'g', 'c', 't', 'c', 'a']], dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
31	3	0	0	0	0	0	0	+	hg19_dna	50	1	35	chr1	249250621	61700837	61700871	1	34,	1,	61700837,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 28)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 44))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 191154276)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr4       37558157 ????????????????----------?????????????????? 37558191
                  0 ||||||||||----------------||||||||||||||||||       44
hg19_dna         49 ??????????------????????????????????????????       11
""",
            )
        elif fmt == "pslx":
            self.assertEqual(alignment.target.seq[37558157:37558167], "tgggattaca")
            self.assertEqual(
                alignment.target.seq[37558173:37558191], "accacgcccagccccttt"
            )
            self.assertEqual(
                alignment.query.seq[11:29], self.queries[alignment.query.id][11:29]
            )
            self.assertEqual(
                alignment.query.seq[39:49], self.queries[alignment.query.id][39:49]
            )
            self.assertEqual(
                str(alignment),
                """\
chr4       37558157 tgggattaca??????----------accacgcccagccccttt 37558191
                  0 ||||||||||----------------||||||||||||||||||       44
hg19_dna         49 tgggattaca------??????????accacgcccagccccttt       11
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
        self.assertEqual(alignment.matches, 35)
        self.assertEqual(alignment.misMatches, 2)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 37))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr22")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 51304566)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr22      48997405 ????????????????????????????????????? 48997442
                  0 |||||||||||||||||||||||||||||||||||||       37
hg19_dna         49 ?????????????????????????????????????       12
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[48997405:48997442],
                "tgggattacaggcgggagccaccacgcccagcccctt",
            )
            self.assertEqual(
                alignment.query.seq[12:49], self.queries[alignment.query.id][12:49]
            )
            self.assertEqual(
                str(alignment),
                """\
chr22      48997405 tgggattacaggcgggagccaccacgcccagcccctt 48997442
                  0 ||||||||||||.|.||||||||||||||||||||||       37
hg19_dna         49 tgggattacaggtgtgagccaccacgcccagcccctt       12
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 'c',
           'g', 'g', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
           'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't'],
          ['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
           'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
           'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't']],
         dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
35	2	0	0	0	0	0	0	-	hg19_dna	50	12	49	chr22	51304566	48997405	48997442	1	37,	1,	48997405,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 35)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr2      120641740 ???????????????????????????????????? 120641776
                  0 ||||||||||||||||||||||||||||||||||||        36
hg19_dna         49 ????????????????????????????????????        13
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[120641740:120641776],
                "tgggattacaggcgtgagccaccacgcccagcccct",
            )
            self.assertEqual(
                alignment.query.seq[13:49], self.queries[alignment.query.id][13:49]
            )
            self.assertEqual(
                str(alignment),
                """\
chr2      120641740 tgggattacaggcgtgagccaccacgcccagcccct 120641776
                  0 ||||||||||||.|||||||||||||||||||||||        36
hg19_dna         49 tgggattacaggtgtgagccaccacgcccagcccct        13
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 'c',
           'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
           'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't'],
          ['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
           'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
           'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't']], dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
35	1	0	0	0	0	0	0	-	hg19_dna	50	13	49	chr2	243199373	120641740	120641776	1	36,	1,	120641740,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 39)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr19      54017130 ??????????????????????????????????????? 54017169
                  0 |||||||||||||||||||||||||||||||||||||||       39
hg19_dna         49 ???????????????????????????????????????       10
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[54017130:54017169],
                "tgggattacaggtgtgagccaccacgcccagcccctttg",
            )
            self.assertEqual(
                alignment.query.seq[10:49], self.queries[alignment.query.id][10:49]
            )
            self.assertEqual(
                str(alignment),
                """\
chr19      54017130 tgggattacaggtgtgagccaccacgcccagcccctttg 54017169
                  0 |||||||||||||||||||||||||||||||||||||||       39
hg19_dna         49 tgggattacaggtgtgagccaccacgcccagcccctttg       10
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
           'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
           'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't', 't', 'g'],
          ['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
           'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
           'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't', 't', 'g']],
         dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
39	0	0	0	0	0	0	0	-	hg19_dna	50	10	49	chr19	59128983	54017130	54017169	1	39,	1,	54017130,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 36)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr19        553742 ??????????????????????????????????????? 553781
                  0 |||||||||||||||||||||||||||||||||||||||     39
hg19_dna         49 ???????????????????????????????????????     10
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[553742:553781],
                "tgggatgacaggggtgaggcaccacgcccagcccctttg",
            )
            self.assertEqual(
                alignment.query.seq[10:49], self.queries[alignment.query.id][10:49]
            )
            self.assertEqual(
                str(alignment),
                """\
chr19        553742 tgggatgacaggggtgaggcaccacgcccagcccctttg 553781
                  0 ||||||.|||||.|||||.||||||||||||||||||||     39
hg19_dna         49 tgggattacaggtgtgagccaccacgcccagcccctttg     10
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['t', 'g', 'g', 'g', 'a', 't', 'g', 'a', 'c', 'a', 'g', 'g', 'g',
           'g', 't', 'g', 'a', 'g', 'g', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
           'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't', 't', 'g'],
          ['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
           'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
           'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't', 't', 'g']],
         dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
36	3	0	0	0	0	0	0	-	hg19_dna	50	10	49	chr19	59128983	553742	553781	1	39,	1,	553742,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 33)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr10")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 135534747)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr10      99388555 ???????????????????????????????????? 99388591
                  0 ||||||||||||||||||||||||||||||||||||       36
hg19_dna         49 ????????????????????????????????????       13
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[99388555:99388591],
                "tgggattataggcatgagccaccacgcccagcccct",
            )
            self.assertEqual(
                alignment.query.seq[13:49], self.queries[alignment.query.id][13:49]
            )
            self.assertEqual(
                str(alignment),
                """\
chr10      99388555 tgggattataggcatgagccaccacgcccagcccct 99388591
                  0 ||||||||.|||..||||||||||||||||||||||       36
hg19_dna         49 tgggattacaggtgtgagccaccacgcccagcccct       13
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 't', 'a', 'g', 'g', 'c',
           'a', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
           'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't'],
          ['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
           'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
           'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't']], dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
33	3	0	0	0	0	0	0	-	hg19_dna	50	13	49	chr10	135534747	99388555	99388591	1	36,	1,	99388555,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 24)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 25))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr10")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 135534747)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr10     112178171 ????????????????????????? 112178196
                  0 |||||||||||||||||||||||||        25
hg19_dna         35 ?????????????????????????        10
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[112178171:112178196], "tgagtcaccacgcccagcccctttg"
            )
            self.assertEqual(
                alignment.query.seq[10:35], self.queries[alignment.query.id][10:35]
            )
            self.assertEqual(
                str(alignment),
                """\
chr10     112178171 tgagtcaccacgcccagcccctttg 112178196
                  0 ||||.||||||||||||||||||||        25
hg19_dna         35 tgagccaccacgcccagcccctttg        10
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['t', 'g', 'a', 'g', 't', 'c', 'a', 'c', 'c', 'a', 'c', 'g', 'c',
           'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't', 't', 'g'],
          ['t', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g', 'c',
           'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't', 't', 'g']],
         dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
24	1	0	0	0	0	0	0	-	hg19_dna	50	10	35	chr10	135534747	112178171	112178196	1	25,	15,	112178171,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 35)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr1       39368490 ???????????????????????????????????? 39368526
                  0 ||||||||||||||||||||||||||||||||||||       36
hg19_dna         49 ????????????????????????????????????       13
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[39368490:39368526],
                "tgggattacaggcgtgagccaccacgcccagcccct",
            )
            self.assertEqual(
                alignment.query.seq[13:49], self.queries[alignment.query.id][13:49]
            )
            self.assertEqual(
                str(alignment),
                """\
chr1       39368490 tgggattacaggcgtgagccaccacgcccagcccct 39368526
                  0 ||||||||||||.|||||||||||||||||||||||       36
hg19_dna         49 tgggattacaggtgtgagccaccacgcccagcccct       13
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 'c',
           'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
           'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't'],
          ['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
           'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
           'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't']], dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
35	1	0	0	0	0	0	0	-	hg19_dna	50	13	49	chr1	249250621	39368490	39368526	1	36,	1,	39368490,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 33)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr1      220325687 ?????????????????????????????????? 220325721
                  0 ||||||||||||||||||||||||||||||||||        34
hg19_dna         47 ??????????????????????????????????        13
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[220325687:220325721],
                "ggattacaggcgtgagccaccacgcccagcccct",
            )
            self.assertEqual(
                alignment.query.seq[13:47], self.queries[alignment.query.id][13:47]
            )
            self.assertEqual(
                str(alignment),
                """\
chr1      220325687 ggattacaggcgtgagccaccacgcccagcccct 220325721
                  0 ||||||||||.|||||||||||||||||||||||        34
hg19_dna         47 ggattacaggtgtgagccaccacgcccagcccct        13
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 'c', 'g', 't',
           'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g', 'c', 'c',
           'c', 'a', 'g', 'c', 'c', 'c', 'c', 't'],
          ['g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't', 'g', 't',
           'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g', 'c', 'c',
           'c', 'a', 'g', 'c', 'c', 'c', 'c', 't']], dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
33	1	0	0	0	0	0	0	-	hg19_dna	50	13	47	chr1	249250621	220325687	220325721	1	34,	3,	220325687,
""",
        )
        self.assertRaises(StopIteration, next, alignments)

    def test_writing_psl_34_001(self):
        """Test writing the alignments in psl_34_001.psl."""
        path = "Blat/psl_34_001.psl"
        with open(path) as stream:
            original_data = stream.read()
        alignments = Align.parse(path, "psl")
        stream = StringIO()
        n = Align.write(alignments, stream, "psl")
        self.assertEqual(n, 22)
        stream.seek(0)
        written_data = stream.read()
        stream.close()
        self.assertEqual(original_data, written_data)

    def test_reading_psl_34_002(self):
        """Test parsing psl_34_002.psl and pslx_34_002.pslx."""
        path = "Blat/psl_34_002.psl"
        self.check_reading_psl_34_002(path)
        path = "Blat/pslx_34_002.pslx"
        self.check_reading_psl_34_002(path)

    def check_reading_psl_34_002(self, path):
        """Check parsing psl_34_002.psl or pslx_34_002.pslx."""
        alignments = Align.parse(path, "psl")
        self.assertEqual(alignments.metadata["psLayout version"], "3")
        self.assertRaises(StopIteration, next, alignments)

    def test_writing_psl_34_002(self):
        """Test writing the alignments in psl_34_002.psl."""
        path = "Blat/psl_34_002.psl"
        with open(path) as stream:
            original_data = stream.read()
        alignments = Align.parse(path, "psl")
        stream = StringIO()
        n = Align.write(alignments, stream, "psl")
        self.assertEqual(n, 0)
        stream.seek(0)
        written_data = stream.read()
        stream.close()
        self.assertEqual(original_data, written_data)

    def test_reading_psl_34_003(self):
        """Test parsing psl_34_003.psl and pslx_34_003.pslx."""
        self.check_reading_psl_34_003("psl")
        self.check_reading_psl_34_003("pslx")

    def check_reading_psl_34_003(self, fmt):
        """Check parsing psl_34_003.psl or pslx_34_003.pslx."""
        path = "Blat/%s_34_003.%s" % (fmt, fmt)
        alignments = Align.parse(path, "psl")
        self.assertEqual(alignments.metadata["psLayout version"], "3")
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 16)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 16))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertEqual(len(alignment.target.seq), 191154276)
        self.assertEqual(len(alignment.query.seq), 33)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr4       61646095 ???????????????? 61646111
                  0 ||||||||||||||||       16
hg18_dna         11 ????????????????       27
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[61646095:61646111], "aggtaaactgccttca"
            )
            self.assertEqual(
                alignment.query.seq[11:27], self.queries[alignment.query.id][11:27]
            )
            self.assertEqual(
                str(alignment),
                """\
chr4       61646095 aggtaaactgccttca 61646111
                  0 ||||||||||||||||       16
hg18_dna         11 aggtaaactgccttca       27
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['a', 'g', 'g', 't', 'a', 'a', 'a', 'c', 't', 'g', 'c', 'c', 't',
           't', 'c', 'a'],
          ['a', 'g', 'g', 't', 'a', 'a', 'a', 'c', 't', 'g', 'c', 'c', 't',
           't', 'c', 'a']], dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
16	0	0	0	0	0	0	0	+	hg18_dna	33	11	27	chr4	191154276	61646095	61646111	1	16,	11,	61646095,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 33)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 33))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 33)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr1       10271783 ????????????????????????????????? 10271816
                  0 |||||||||||||||||||||||||||||||||       33
hg18_dna          0 ?????????????????????????????????       33
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[10271783:10271816],
                "atgagcttccaaggtaaactgccttcaagattc",
            )
            self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
            self.assertEqual(
                str(alignment),
                """\
chr1       10271783 atgagcttccaaggtaaactgccttcaagattc 10271816
                  0 |||||||||||||||||||||||||||||||||       33
hg18_dna          0 atgagcttccaaggtaaactgccttcaagattc       33
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['a', 't', 'g', 'a', 'g', 'c', 't', 't', 'c', 'c', 'a', 'a', 'g',
           'g', 't', 'a', 'a', 'a', 'c', 't', 'g', 'c', 'c', 't', 't', 'c',
           'a', 'a', 'g', 'a', 't', 't', 'c'],
          ['a', 't', 'g', 'a', 'g', 'c', 't', 't', 'c', 'c', 'a', 'a', 'g',
           'g', 't', 'a', 'a', 'a', 'c', 't', 'g', 'c', 'c', 't', 't', 'c',
           'a', 'a', 'g', 'a', 't', 't', 'c']], dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
33	0	0	0	0	0	0	0	+	hg18_dna	33	0	33	chr1	249250621	10271783	10271816	1	33,	0,	10271783,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 17)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 17))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 33)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr2       53575980 ????????????????? 53575997
                  0 |||||||||||||||||       17
hg18_dna         25 ?????????????????        8
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[53575980:53575997], "aaggcagtttaccttgg"
            )
            self.assertEqual(
                alignment.query.seq[8:25], self.queries[alignment.query.id][8:25]
            )
            self.assertEqual(
                str(alignment),
                """\
chr2       53575980 aaggcagtttaccttgg 53575997
                  0 |||||||||||||||||       17
hg18_dna         25 aaggcagtttaccttgg        8
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['a', 'a', 'g', 'g', 'c', 'a', 'g', 't', 't', 't', 'a', 'c', 'c',
           't', 't', 'g', 'g'],
          ['a', 'a', 'g', 'g', 'c', 'a', 'g', 't', 't', 't', 'a', 'c', 'c',
           't', 't', 'g', 'g']], dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
17	0	0	0	0	0	0	0	-	hg18_dna	33	8	25	chr2	243199373	53575980	53575997	1	17,	8,	53575980,
""",
        )
        self.assertRaises(StopIteration, next, alignments)

    def test_writing_psl_34_003(self):
        """Test writing the alignments in psl_34_003.psl."""
        path = "Blat/psl_34_003.psl"
        with open(path) as stream:
            original_data = stream.read()
        alignments = Align.parse(path, "psl")
        stream = StringIO()
        n = Align.write(alignments, stream, "psl")
        self.assertEqual(n, 3)
        stream.seek(0)
        written_data = stream.read()
        stream.close()
        self.assertEqual(original_data, written_data)

    def test_reading_psl_34_004(self):
        """Test parsing psl_34_004.psl and pslx_34_004.pslx."""
        self.check_reading_psl_34_004("psl")
        self.check_reading_psl_34_004("pslx")

    def check_reading_psl_34_004(self, fmt):
        """Check parsing psl_34_004.psl or pslx_34_004.pslx."""
        path = "Blat/%s_34_004.%s" % (fmt, fmt)
        alignments = Align.parse(path, "psl")
        self.assertEqual(alignments.metadata["psLayout version"], "3")
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 38)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 41))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr9")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 141213431)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr9       85737865 ????????????????????????????????????????? 85737906
                  0 |||||||||||||||||||||||||||||||||||||||||       41
hg19_dna          9 ?????????????????????????????????????????       50
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[85737865:85737906],
                "acaaaggggctgggcgcagtggctcacgcctgtaatcccaa",
            )
            self.assertEqual(
                alignment.query.seq[9:50], self.queries[alignment.query.id][9:50]
            )
            self.assertEqual(
                str(alignment),
                """\
chr9       85737865 acaaaggggctgggcgcagtggctcacgcctgtaatcccaa 85737906
                  0 ||||||||||||||||..|||||||||.|||||||||||||       41
hg19_dna          9 acaaaggggctgggcgtggtggctcacacctgtaatcccaa       50
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['a', 'c', 'a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g',
           'g', 'c', 'g', 'c', 'a', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a',
           'c', 'g', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c', 'c',
           'a', 'a'],
          ['a', 'c', 'a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g',
           'g', 'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a',
           'c', 'a', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c', 'c',
           'a', 'a']], dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
38	3	0	0	0	0	0	0	+	hg19_dna	50	9	50	chr9	141213431	85737865	85737906	1	41,	9,	85737865,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 41)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 41))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr8")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 146364022)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr8       95160479 ????????????????????????????????????????? 95160520
                  0 |||||||||||||||||||||||||||||||||||||||||       41
hg19_dna          8 ?????????????????????????????????????????       49
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[95160479:95160520],
                "cacaaaggggctgggcgtggtggctcacacctgtaatccca",
            )
            self.assertEqual(
                alignment.query.seq[8:49], self.queries[alignment.query.id][8:49]
            )
            self.assertEqual(
                str(alignment),
                """\
chr8       95160479 cacaaaggggctgggcgtggtggctcacacctgtaatccca 95160520
                  0 |||||||||||||||||||||||||||||||||||||||||       41
hg19_dna          8 cacaaaggggctgggcgtggtggctcacacctgtaatccca       49
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['c', 'a', 'c', 'a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g',
           'g', 'g', 'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c',
           'a', 'c', 'a', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c',
           'c', 'a'],
          ['c', 'a', 'c', 'a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g',
           'g', 'g', 'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c',
           'a', 'c', 'a', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c',
           'c', 'a']], dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
41	0	0	0	0	0	0	0	+	hg19_dna	50	8	49	chr8	146364022	95160479	95160520	1	41,	8,	95160479,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 33)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr22")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 51304566)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr22      42144400 ???????????????????????????????????? 42144436
                  0 ||||||||||||||||||||||||||||||||||||       36
hg19_dna         11 ????????????????????????????????????       47
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[42144400:42144436],
                "aaaggggctgggcgtggtagctcatgcctgtaatcc",
            )
            self.assertEqual(
                alignment.query.seq[11:47], self.queries[alignment.query.id][11:47]
            )
            self.assertEqual(
                str(alignment),
                """\
chr22      42144400 aaaggggctgggcgtggtagctcatgcctgtaatcc 42144436
                  0 ||||||||||||||||||.|||||..||||||||||       36
hg19_dna         11 aaaggggctgggcgtggtggctcacacctgtaatcc       47
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g', 'c',
           'g', 't', 'g', 'g', 't', 'a', 'g', 'c', 't', 'c', 'a', 't', 'g',
           'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c'],
          ['a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g', 'c',
           'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a', 'c', 'a',
           'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c']], dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
33	3	0	0	0	0	0	0	+	hg19_dna	50	11	47	chr22	51304566	42144400	42144436	1	36,	11,	42144400,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 43)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 48))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr2      183925984 ??????----?????????????????????????????????????? 183926028
                  0 ||||||----||||||||||||||||||||||||||||||||||||||        48
hg19_dna          1 ????????????????????????????????????????????????        49
""",
            )
        elif fmt == "pslx":
            self.assertEqual(alignment.target.seq[183925984:183925990], "aaaaat")
            self.assertEqual(
                alignment.target.seq[183925990:183926028],
                "aaaggggctgggcgtggtggctcacgcctgtaatccca",
            )
            self.assertEqual(
                str(alignment),
                """\
chr2      183925984 aaaaat----aaaggggctgggcgtggtggctcacgcctgtaatccca 183926028
                  0 ||||||----|||||||||||||||||||||||||.||||||||||||        48
hg19_dna          1 aaaaat????aaaggggctgggcgtggtggctcacacctgtaatccca        49
""",
            )
            self.assertEqual(
                alignment.query.seq[1:7], self.queries[alignment.query.id][1:7]
            )
            self.assertEqual(
                alignment.query.seq[38:49], self.queries[alignment.query.id][38:49]
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
            format(alignment, "psl"),
            """\
43	1	0	0	1	4	0	0	+	hg19_dna	50	1	49	chr2	243199373	183925984	183926028	2	6,38,	1,11,	183925984,183925990,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 34)
        self.assertEqual(alignment.misMatches, 2)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 170))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
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
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[35483340:35483365], "caaaggggctgggcgtagtggctga"
            )
            self.assertEqual(alignment.target.seq[35483499:35483510], "cacctgtaatc")
            self.assertEqual(
                alignment.query.seq[10:46], self.queries[alignment.query.id][10:46]
            )
            self.assertEqual(
                str(alignment),
                """\
chr19      35483340 caaaggggctgggcgtagtggctga???????????????????????????????????
                  0 ||||||||||||||||.||||||.|-----------------------------------
hg19_dna         10 caaaggggctgggcgtggtggctca-----------------------------------

chr19      35483400 ????????????????????????????????????????????????????????????
                 60 ------------------------------------------------------------
hg19_dna         35 ------------------------------------------------------------

chr19      35483460 ???????????????????????????????????????cacctgtaatc 35483510
                120 ---------------------------------------|||||||||||      170
hg19_dna         35 ---------------------------------------cacctgtaatc       46
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
            format(alignment, "psl"),
            """\
34	2	0	0	0	0	1	134	+	hg19_dna	50	10	46	chr19	59128983	35483340	35483510	2	25,11,	10,35,	35483340,35483499,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 39)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr18")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 78077248)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr18      23891310 ??????????????????????????????????????? 23891349
                  0 |||||||||||||||||||||||||||||||||||||||       39
hg19_dna         10 ???????????????????????????????????????       49
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[23891310:23891349],
                "caaaggggctgggcgtggtggctcacacctgtaatccca",
            )
            self.assertEqual(
                alignment.query.seq[10:49], self.queries[alignment.query.id][10:49]
            )
            self.assertEqual(
                str(alignment),
                """\
chr18      23891310 caaaggggctgggcgtggtggctcacacctgtaatccca 23891349
                  0 |||||||||||||||||||||||||||||||||||||||       39
hg19_dna         10 caaaggggctgggcgtggtggctcacacctgtaatccca       49
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['c', 'a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g',
           'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a', 'c',
           'a', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c', 'c', 'a'],
          ['c', 'a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g',
           'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a', 'c',
           'a', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c', 'c', 'a']],
         dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
39	0	0	0	0	0	0	0	+	hg19_dna	50	10	49	chr18	78077248	23891310	23891349	1	39,	10,	23891310,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 27)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 28))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr18")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 78077248)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr18      43252217 ???????????????????????????? 43252245
                  0 ||||||||||||||||||||||||||||       28
hg19_dna         21 ????????????????????????????       49
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[43252217:43252245], "ggcgtggtggctcacgcctgtaatccca"
            )
            self.assertEqual(
                alignment.query.seq[21:49], self.queries[alignment.query.id][21:49]
            )
            self.assertEqual(
                str(alignment),
                """\
chr18      43252217 ggcgtggtggctcacgcctgtaatccca 43252245
                  0 |||||||||||||||.||||||||||||       28
hg19_dna         21 ggcgtggtggctcacacctgtaatccca       49
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['g', 'g', 'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c',
           'a', 'c', 'g', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c',
           'c', 'a'],
          ['g', 'g', 'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c',
           'a', 'c', 'a', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c',
           'c', 'a']], dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
27	1	0	0	0	0	0	0	+	hg19_dna	50	21	49	chr18	78077248	43252217	43252245	1	28,	21,	43252217,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 44)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 54))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 115169878)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
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
        elif fmt == "pslx":
            self.assertEqual(alignment.target.seq[52759147:52759154], "aaaaatt")
            self.assertEqual(
                alignment.target.seq[52759160:52759198],
                "aaaggggctgggcgtggtggctcacgcctgtaatccca",
            )
            self.assertEqual(
                alignment.query.seq[1:8], self.queries[alignment.query.id][1:8]
            )
            self.assertEqual(
                alignment.query.seq[38:49], self.queries[alignment.query.id][38:49]
            )
            self.assertEqual(
                str(alignment),
                """\
chr13      52759147 aaaaatt??????---aaaggggctgggcgtggtggctcacgcctgtaatccca
                  0 |||||||---------|||||||||||||||||||||||||.||||||||||||
hg19_dna          1 aaaaatt------???aaaggggctgggcgtggtggctcacacctgtaatccca

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
            format(alignment, "psl"),
            """\
44	1	0	0	1	3	1	6	+	hg19_dna	50	1	49	chr13	115169878	52759147	52759198	2	7,38,	1,11,	52759147,52759160,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 50)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 50))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr1        1207056 ?????????????????????????????????????????????????? 1207106
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||      50
hg19_dna          0 ??????????????????????????????????????????????????      50
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[1207056:1207106],
                "caaaaattcacaaaggggctgggcgtggtggctcacacctgtaatcccaa",
            )
            self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
            self.assertEqual(
                str(alignment),
                """\
chr1        1207056 caaaaattcacaaaggggctgggcgtggtggctcacacctgtaatcccaa 1207106
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||      50
hg19_dna          0 caaaaattcacaaaggggctgggcgtggtggctcacacctgtaatcccaa      50
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['c', 'a', 'a', 'a', 'a', 'a', 't', 't', 'c', 'a', 'c', 'a', 'a',
           'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g', 'c', 'g', 't',
           'g', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a', 'c', 'a', 'c', 'c',
           't', 'g', 't', 'a', 'a', 't', 'c', 'c', 'c', 'a', 'a'],
          ['c', 'a', 'a', 'a', 'a', 'a', 't', 't', 'c', 'a', 'c', 'a', 'a',
           'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g', 'c', 'g', 't',
           'g', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a', 'c', 'a', 'c', 'c',
           't', 'g', 't', 'a', 'a', 't', 'c', 'c', 'c', 'a', 'a']],
         dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
50	0	0	0	0	0	0	0	+	hg19_dna	50	0	50	chr1	249250621	1207056	1207106	1	50,	0,	1207056,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 31)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr1       61700837 ?????????????????????????????????? 61700871
                  0 ||||||||||||||||||||||||||||||||||       34
hg19_dna          1 ??????????????????????????????????       35
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[61700837:61700871],
                "aaaaatgaacaaaggggctgggcgcggtggctca",
            )
            self.assertEqual(
                alignment.query.seq[1:35], self.queries[alignment.query.id][1:35]
            )
            self.assertEqual(
                str(alignment),
                """\
chr1       61700837 aaaaatgaacaaaggggctgggcgcggtggctca 61700871
                  0 ||||||..||||||||||||||||.|||||||||       34
hg19_dna          1 aaaaattcacaaaggggctgggcgtggtggctca       35
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['a', 'a', 'a', 'a', 'a', 't', 'g', 'a', 'a', 'c', 'a', 'a', 'a',
           'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g', 'c', 'g', 'c', 'g',
           'g', 't', 'g', 'g', 'c', 't', 'c', 'a'],
          ['a', 'a', 'a', 'a', 'a', 't', 't', 'c', 'a', 'c', 'a', 'a', 'a',
           'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g', 'c', 'g', 't', 'g',
           'g', 't', 'g', 'g', 'c', 't', 'c', 'a']], dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
31	3	0	0	0	0	0	0	+	hg19_dna	50	1	35	chr1	249250621	61700837	61700871	1	34,	1,	61700837,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 28)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 44))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 191154276)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr4       37558157 ????????????????----------?????????????????? 37558191
                  0 ||||||||||----------------||||||||||||||||||       44
hg19_dna         49 ??????????------????????????????????????????       11
""",
            )
        elif fmt == "pslx":
            self.assertEqual(alignment.target.seq[37558157:37558167], "tgggattaca")
            self.assertEqual(
                alignment.target.seq[37558173:37558191], "accacgcccagccccttt"
            )
            self.assertEqual(
                alignment.query.seq[11:29], self.queries[alignment.query.id][11:29]
            )
            self.assertEqual(
                alignment.query.seq[39:49], self.queries[alignment.query.id][39:49]
            )
            self.assertEqual(
                str(alignment),
                """\
chr4       37558157 tgggattaca??????----------accacgcccagccccttt 37558191
                  0 ||||||||||----------------||||||||||||||||||       44
hg19_dna         49 tgggattaca------??????????accacgcccagccccttt       11
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
        self.assertEqual(alignment.matches, 35)
        self.assertEqual(alignment.misMatches, 2)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 37))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr22")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 51304566)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr22      48997405 ????????????????????????????????????? 48997442
                  0 |||||||||||||||||||||||||||||||||||||       37
hg19_dna         49 ?????????????????????????????????????       12
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[48997405:48997442],
                "tgggattacaggcgggagccaccacgcccagcccctt",
            )
            self.assertEqual(
                alignment.query.seq[12:49], self.queries[alignment.query.id][12:49]
            )
            self.assertEqual(
                str(alignment),
                """\
chr22      48997405 tgggattacaggcgggagccaccacgcccagcccctt 48997442
                  0 ||||||||||||.|.||||||||||||||||||||||       37
hg19_dna         49 tgggattacaggtgtgagccaccacgcccagcccctt       12
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 'c',
           'g', 'g', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
           'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't'],
          ['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
           'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
           'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't']],
         dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
35	2	0	0	0	0	0	0	-	hg19_dna	50	12	49	chr22	51304566	48997405	48997442	1	37,	1,	48997405,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 35)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr2      120641740 ???????????????????????????????????? 120641776
                  0 ||||||||||||||||||||||||||||||||||||        36
hg19_dna         49 ????????????????????????????????????        13
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[120641740:120641776],
                "tgggattacaggcgtgagccaccacgcccagcccct",
            )
            self.assertEqual(
                alignment.query.seq[13:49], self.queries[alignment.query.id][13:49]
            )
            self.assertEqual(
                str(alignment),
                """\
chr2      120641740 tgggattacaggcgtgagccaccacgcccagcccct 120641776
                  0 ||||||||||||.|||||||||||||||||||||||        36
hg19_dna         49 tgggattacaggtgtgagccaccacgcccagcccct        13
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 'c',
           'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
           'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't'],
          ['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
           'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
           'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't']], dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
35	1	0	0	0	0	0	0	-	hg19_dna	50	13	49	chr2	243199373	120641740	120641776	1	36,	1,	120641740,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 39)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr19      54017130 ??????????????????????????????????????? 54017169
                  0 |||||||||||||||||||||||||||||||||||||||       39
hg19_dna         49 ???????????????????????????????????????       10
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[54017130:54017169],
                "tgggattacaggtgtgagccaccacgcccagcccctttg",
            )
            self.assertEqual(
                alignment.query.seq[10:49], self.queries[alignment.query.id][10:49]
            )
            self.assertEqual(
                str(alignment),
                """\
chr19      54017130 tgggattacaggtgtgagccaccacgcccagcccctttg 54017169
                  0 |||||||||||||||||||||||||||||||||||||||       39
hg19_dna         49 tgggattacaggtgtgagccaccacgcccagcccctttg       10
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
           'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
           'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't', 't', 'g'],
          ['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
           'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
           'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't', 't', 'g']],
         dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
39	0	0	0	0	0	0	0	-	hg19_dna	50	10	49	chr19	59128983	54017130	54017169	1	39,	1,	54017130,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 36)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr19        553742 ??????????????????????????????????????? 553781
                  0 |||||||||||||||||||||||||||||||||||||||     39
hg19_dna         49 ???????????????????????????????????????     10
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[553742:553781],
                "tgggatgacaggggtgaggcaccacgcccagcccctttg",
            )
            self.assertEqual(
                alignment.query.seq[10:49], self.queries[alignment.query.id][10:49]
            )
            self.assertEqual(
                str(alignment),
                """\
chr19        553742 tgggatgacaggggtgaggcaccacgcccagcccctttg 553781
                  0 ||||||.|||||.|||||.||||||||||||||||||||     39
hg19_dna         49 tgggattacaggtgtgagccaccacgcccagcccctttg     10
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['t', 'g', 'g', 'g', 'a', 't', 'g', 'a', 'c', 'a', 'g', 'g', 'g',
           'g', 't', 'g', 'a', 'g', 'g', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
           'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't', 't', 'g'],
          ['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
           'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
           'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't', 't', 'g']],
         dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
36	3	0	0	0	0	0	0	-	hg19_dna	50	10	49	chr19	59128983	553742	553781	1	39,	1,	553742,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 33)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr10")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 135534747)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr10      99388555 ???????????????????????????????????? 99388591
                  0 ||||||||||||||||||||||||||||||||||||       36
hg19_dna         49 ????????????????????????????????????       13
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[99388555:99388591],
                "tgggattataggcatgagccaccacgcccagcccct",
            )
            self.assertEqual(
                alignment.query.seq[13:49], self.queries[alignment.query.id][13:49]
            )
            self.assertEqual(
                str(alignment),
                """\
chr10      99388555 tgggattataggcatgagccaccacgcccagcccct 99388591
                  0 ||||||||.|||..||||||||||||||||||||||       36
hg19_dna         49 tgggattacaggtgtgagccaccacgcccagcccct       13
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 't', 'a', 'g', 'g', 'c',
           'a', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
           'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't'],
          ['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
           'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
           'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't']], dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
33	3	0	0	0	0	0	0	-	hg19_dna	50	13	49	chr10	135534747	99388555	99388591	1	36,	1,	99388555,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 24)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 25))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr10")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 135534747)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr10     112178171 ????????????????????????? 112178196
                  0 |||||||||||||||||||||||||        25
hg19_dna         35 ?????????????????????????        10
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[112178171:112178196], "tgagtcaccacgcccagcccctttg"
            )
            self.assertEqual(
                alignment.query.seq[10:35], self.queries[alignment.query.id][10:35]
            )
            self.assertEqual(
                str(alignment),
                """\
chr10     112178171 tgagtcaccacgcccagcccctttg 112178196
                  0 ||||.||||||||||||||||||||        25
hg19_dna         35 tgagccaccacgcccagcccctttg        10
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['t', 'g', 'a', 'g', 't', 'c', 'a', 'c', 'c', 'a', 'c', 'g', 'c',
           'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't', 't', 'g'],
          ['t', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g', 'c',
           'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't', 't', 'g']],
         dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
24	1	0	0	0	0	0	0	-	hg19_dna	50	10	35	chr10	135534747	112178171	112178196	1	25,	15,	112178171,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 35)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr1       39368490 ???????????????????????????????????? 39368526
                  0 ||||||||||||||||||||||||||||||||||||       36
hg19_dna         49 ????????????????????????????????????       13
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[39368490:39368526],
                "tgggattacaggcgtgagccaccacgcccagcccct",
            )
            self.assertEqual(
                alignment.query.seq[13:49], self.queries[alignment.query.id][13:49]
            )
            self.assertEqual(
                str(alignment),
                """\
chr1       39368490 tgggattacaggcgtgagccaccacgcccagcccct 39368526
                  0 ||||||||||||.|||||||||||||||||||||||       36
hg19_dna         49 tgggattacaggtgtgagccaccacgcccagcccct       13
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 'c',
           'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
           'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't'],
          ['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
           'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
           'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't']], dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
35	1	0	0	0	0	0	0	-	hg19_dna	50	13	49	chr1	249250621	39368490	39368526	1	36,	1,	39368490,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 33)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr1      220325687 ?????????????????????????????????? 220325721
                  0 ||||||||||||||||||||||||||||||||||        34
hg19_dna         47 ??????????????????????????????????        13
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[220325687:220325721],
                "ggattacaggcgtgagccaccacgcccagcccct",
            )
            self.assertEqual(
                alignment.query.seq[13:47], self.queries[alignment.query.id][13:47]
            )
            self.assertEqual(
                str(alignment),
                """\
chr1      220325687 ggattacaggcgtgagccaccacgcccagcccct 220325721
                  0 ||||||||||.|||||||||||||||||||||||        34
hg19_dna         47 ggattacaggtgtgagccaccacgcccagcccct        13
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 'c', 'g', 't',
           'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g', 'c', 'c',
           'c', 'a', 'g', 'c', 'c', 'c', 'c', 't'],
          ['g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't', 'g', 't',
           'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g', 'c', 'c',
           'c', 'a', 'g', 'c', 'c', 'c', 'c', 't']], dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
33	1	0	0	0	0	0	0	-	hg19_dna	50	13	47	chr1	249250621	220325687	220325721	1	34,	3,	220325687,
""",
        )
        self.assertRaises(StopIteration, next, alignments)

    def test_writing_psl_34_004(self):
        """Test writing the alignments in psl_34_004.psl."""
        path = "Blat/psl_34_004.psl"
        with open(path) as stream:
            original_data = stream.read()
        alignments = Align.parse(path, "psl")
        stream = StringIO()
        n = Align.write(alignments, stream, "psl")
        self.assertEqual(n, 19)
        stream.seek(0)
        written_data = stream.read()
        stream.close()
        self.assertEqual(original_data, written_data)

    def test_reading_psl_34_005(self):
        """Test parsing psl_34_005.psl and pslx_34_005.pslx."""
        self.check_reading_psl_34_005("psl")
        self.check_reading_psl_34_005("pslx")

    def check_reading_psl_34_005(self, fmt):
        """Check parsing psl_34_005.psl or pslx_34_005.pslx."""
        path = "Blat/%s_34_005.%s" % (fmt, fmt)
        alignments = Align.parse(path, "psl")
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 16)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 16))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertEqual(len(alignment.target.seq), 191154276)
        self.assertEqual(len(alignment.query.seq), 33)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr4       61646095 ???????????????? 61646111
                  0 ||||||||||||||||       16
hg18_dna         11 ????????????????       27
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[61646095:61646111], "aggtaaactgccttca"
            )
            self.assertEqual(
                alignment.query.seq[11:27], self.queries[alignment.query.id][11:27]
            )
            self.assertEqual(
                str(alignment),
                """\
chr4       61646095 aggtaaactgccttca 61646111
                  0 ||||||||||||||||       16
hg18_dna         11 aggtaaactgccttca       27
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['a', 'g', 'g', 't', 'a', 'a', 'a', 'c', 't', 'g', 'c', 'c', 't',
           't', 'c', 'a'],
          ['a', 'g', 'g', 't', 'a', 'a', 'a', 'c', 't', 'g', 'c', 'c', 't',
           't', 'c', 'a']], dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
16	0	0	0	0	0	0	0	+	hg18_dna	33	11	27	chr4	191154276	61646095	61646111	1	16,	11,	61646095,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 33)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 33))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 33)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr1       10271783 ????????????????????????????????? 10271816
                  0 |||||||||||||||||||||||||||||||||       33
hg18_dna          0 ?????????????????????????????????       33
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[10271783:10271816],
                "atgagcttccaaggtaaactgccttcaagattc",
            )
            self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
            self.assertEqual(
                str(alignment),
                """\
chr1       10271783 atgagcttccaaggtaaactgccttcaagattc 10271816
                  0 |||||||||||||||||||||||||||||||||       33
hg18_dna          0 atgagcttccaaggtaaactgccttcaagattc       33
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['a', 't', 'g', 'a', 'g', 'c', 't', 't', 'c', 'c', 'a', 'a', 'g',
           'g', 't', 'a', 'a', 'a', 'c', 't', 'g', 'c', 'c', 't', 't', 'c',
           'a', 'a', 'g', 'a', 't', 't', 'c'],
          ['a', 't', 'g', 'a', 'g', 'c', 't', 't', 'c', 'c', 'a', 'a', 'g',
           'g', 't', 'a', 'a', 'a', 'c', 't', 'g', 'c', 'c', 't', 't', 'c',
           'a', 'a', 'g', 'a', 't', 't', 'c']], dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
33	0	0	0	0	0	0	0	+	hg18_dna	33	0	33	chr1	249250621	10271783	10271816	1	33,	0,	10271783,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 17)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 17))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 33)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr2       53575980 ????????????????? 53575997
                  0 |||||||||||||||||       17
hg18_dna         25 ?????????????????        8
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[53575980:53575997], "aaggcagtttaccttgg"
            )
            self.assertEqual(
                alignment.query.seq[8:25], self.queries[alignment.query.id][8:25]
            )
            self.assertEqual(
                str(alignment),
                """\
chr2       53575980 aaggcagtttaccttgg 53575997
                  0 |||||||||||||||||       17
hg18_dna         25 aaggcagtttaccttgg        8
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['a', 'a', 'g', 'g', 'c', 'a', 'g', 't', 't', 't', 'a', 'c', 'c',
           't', 't', 'g', 'g'],
          ['a', 'a', 'g', 'g', 'c', 'a', 'g', 't', 't', 't', 'a', 'c', 'c',
           't', 't', 'g', 'g']], dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
17	0	0	0	0	0	0	0	-	hg18_dna	33	8	25	chr2	243199373	53575980	53575997	1	17,	8,	53575980,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 38)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 41))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr9")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 141213431)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr9       85737865 ????????????????????????????????????????? 85737906
                  0 |||||||||||||||||||||||||||||||||||||||||       41
hg19_dna          9 ?????????????????????????????????????????       50
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[85737865:85737906],
                "acaaaggggctgggcgcagtggctcacgcctgtaatcccaa",
            )
            self.assertEqual(
                alignment.query.seq[9:], self.queries[alignment.query.id][9:]
            )
            self.assertEqual(
                str(alignment),
                """\
chr9       85737865 acaaaggggctgggcgcagtggctcacgcctgtaatcccaa 85737906
                  0 ||||||||||||||||..|||||||||.|||||||||||||       41
hg19_dna          9 acaaaggggctgggcgtggtggctcacacctgtaatcccaa       50
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['a', 'c', 'a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g',
           'g', 'c', 'g', 'c', 'a', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a',
           'c', 'g', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c', 'c',
           'a', 'a'],
          ['a', 'c', 'a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g',
           'g', 'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a',
           'c', 'a', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c', 'c',
           'a', 'a']], dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
38	3	0	0	0	0	0	0	+	hg19_dna	50	9	50	chr9	141213431	85737865	85737906	1	41,	9,	85737865,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 41)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 41))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr8")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 146364022)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr8       95160479 ????????????????????????????????????????? 95160520
                  0 |||||||||||||||||||||||||||||||||||||||||       41
hg19_dna          8 ?????????????????????????????????????????       49
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[95160479:95160520],
                "cacaaaggggctgggcgtggtggctcacacctgtaatccca",
            )
            self.assertEqual(
                alignment.query.seq[8:49], self.queries[alignment.query.id][8:49]
            )
            self.assertEqual(
                str(alignment),
                """\
chr8       95160479 cacaaaggggctgggcgtggtggctcacacctgtaatccca 95160520
                  0 |||||||||||||||||||||||||||||||||||||||||       41
hg19_dna          8 cacaaaggggctgggcgtggtggctcacacctgtaatccca       49
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['c', 'a', 'c', 'a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g',
           'g', 'g', 'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c',
           'a', 'c', 'a', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c',
           'c', 'a'],
          ['c', 'a', 'c', 'a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g',
           'g', 'g', 'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c',
           'a', 'c', 'a', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c',
           'c', 'a']], dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
41	0	0	0	0	0	0	0	+	hg19_dna	50	8	49	chr8	146364022	95160479	95160520	1	41,	8,	95160479,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 33)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr22")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 51304566)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr22      42144400 ???????????????????????????????????? 42144436
                  0 ||||||||||||||||||||||||||||||||||||       36
hg19_dna         11 ????????????????????????????????????       47
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[42144400:42144436],
                "aaaggggctgggcgtggtagctcatgcctgtaatcc",
            )
            self.assertEqual(
                alignment.query.seq[11:47], self.queries[alignment.query.id][11:47]
            )
            self.assertEqual(
                str(alignment),
                """\
chr22      42144400 aaaggggctgggcgtggtagctcatgcctgtaatcc 42144436
                  0 ||||||||||||||||||.|||||..||||||||||       36
hg19_dna         11 aaaggggctgggcgtggtggctcacacctgtaatcc       47
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g', 'c',
           'g', 't', 'g', 'g', 't', 'a', 'g', 'c', 't', 'c', 'a', 't', 'g',
           'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c'],
          ['a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g', 'c',
           'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a', 'c', 'a',
           'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c']], dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
33	3	0	0	0	0	0	0	+	hg19_dna	50	11	47	chr22	51304566	42144400	42144436	1	36,	11,	42144400,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 43)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 48))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr2      183925984 ??????----?????????????????????????????????????? 183926028
                  0 ||||||----||||||||||||||||||||||||||||||||||||||        48
hg19_dna          1 ????????????????????????????????????????????????        49
""",
            )
        elif fmt == "pslx":
            self.assertEqual(alignment.target.seq[183925984:183925990], "aaaaat")
            self.assertEqual(
                alignment.target.seq[183925990:183926028],
                "aaaggggctgggcgtggtggctcacgcctgtaatccca",
            )
            self.assertEqual(
                alignment.query.seq[1:7], self.queries[alignment.query.id][1:7]
            )
            self.assertEqual(
                alignment.query.seq[11:49], self.queries[alignment.query.id][11:49]
            )
            self.assertEqual(
                str(alignment),
                """\
chr2      183925984 aaaaat----aaaggggctgggcgtggtggctcacgcctgtaatccca 183926028
                  0 ||||||----|||||||||||||||||||||||||.||||||||||||        48
hg19_dna          1 aaaaat????aaaggggctgggcgtggtggctcacacctgtaatccca        49
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
            format(alignment, "psl"),
            """\
43	1	0	0	1	4	0	0	+	hg19_dna	50	1	49	chr2	243199373	183925984	183926028	2	6,38,	1,11,	183925984,183925990,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 34)
        self.assertEqual(alignment.misMatches, 2)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 170))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
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
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[35483340:35483365], "caaaggggctgggcgtagtggctga"
            )
            self.assertEqual(alignment.target.seq[35483499:35483510], "cacctgtaatc")
            self.assertEqual(
                alignment.query.seq[10:46], self.queries[alignment.query.id][10:46]
            )
            self.assertEqual(
                str(alignment),
                """\
chr19      35483340 caaaggggctgggcgtagtggctga???????????????????????????????????
                  0 ||||||||||||||||.||||||.|-----------------------------------
hg19_dna         10 caaaggggctgggcgtggtggctca-----------------------------------

chr19      35483400 ????????????????????????????????????????????????????????????
                 60 ------------------------------------------------------------
hg19_dna         35 ------------------------------------------------------------

chr19      35483460 ???????????????????????????????????????cacctgtaatc 35483510
                120 ---------------------------------------|||||||||||      170
hg19_dna         35 ---------------------------------------cacctgtaatc       46
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
            format(alignment, "psl"),
            """\
34	2	0	0	0	0	1	134	+	hg19_dna	50	10	46	chr19	59128983	35483340	35483510	2	25,11,	10,35,	35483340,35483499,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 39)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr18")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 78077248)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr18      23891310 ??????????????????????????????????????? 23891349
                  0 |||||||||||||||||||||||||||||||||||||||       39
hg19_dna         10 ???????????????????????????????????????       49
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[23891310:23891349],
                "caaaggggctgggcgtggtggctcacacctgtaatccca",
            )
            self.assertEqual(
                alignment.query.seq[10:49], self.queries[alignment.query.id][10:49]
            )
            self.assertEqual(
                str(alignment),
                """\
chr18      23891310 caaaggggctgggcgtggtggctcacacctgtaatccca 23891349
                  0 |||||||||||||||||||||||||||||||||||||||       39
hg19_dna         10 caaaggggctgggcgtggtggctcacacctgtaatccca       49
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['c', 'a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g',
           'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a', 'c',
           'a', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c', 'c', 'a'],
          ['c', 'a', 'a', 'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g',
           'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a', 'c',
           'a', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c', 'c', 'a']],
         dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
39	0	0	0	0	0	0	0	+	hg19_dna	50	10	49	chr18	78077248	23891310	23891349	1	39,	10,	23891310,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 27)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 28))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr18")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 78077248)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr18      43252217 ???????????????????????????? 43252245
                  0 ||||||||||||||||||||||||||||       28
hg19_dna         21 ????????????????????????????       49
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[43252217:43252245], "ggcgtggtggctcacgcctgtaatccca"
            )
            self.assertEqual(
                alignment.query.seq[21:49], self.queries[alignment.query.id][21:49]
            )
            self.assertEqual(
                str(alignment),
                """\
chr18      43252217 ggcgtggtggctcacgcctgtaatccca 43252245
                  0 |||||||||||||||.||||||||||||       28
hg19_dna         21 ggcgtggtggctcacacctgtaatccca       49
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['g', 'g', 'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c',
           'a', 'c', 'g', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c',
           'c', 'a'],
          ['g', 'g', 'c', 'g', 't', 'g', 'g', 't', 'g', 'g', 'c', 't', 'c',
           'a', 'c', 'a', 'c', 'c', 't', 'g', 't', 'a', 'a', 't', 'c', 'c',
           'c', 'a']], dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
27	1	0	0	0	0	0	0	+	hg19_dna	50	21	49	chr18	78077248	43252217	43252245	1	28,	21,	43252217,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 44)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 54))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 115169878)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
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
        elif fmt == "pslx":
            self.assertEqual(alignment.target.seq[52759147:52759154], "aaaaatt")
            self.assertEqual(
                alignment.target.seq[52759160:52759198],
                "aaaggggctgggcgtggtggctcacgcctgtaatccca",
            )
            self.assertEqual(
                alignment.query.seq[1:8], self.queries[alignment.query.id][1:8]
            )
            self.assertEqual(
                alignment.query.seq[38:49], self.queries[alignment.query.id][38:49]
            )
            self.assertEqual(
                str(alignment),
                """\
chr13      52759147 aaaaatt??????---aaaggggctgggcgtggtggctcacgcctgtaatccca
                  0 |||||||---------|||||||||||||||||||||||||.||||||||||||
hg19_dna          1 aaaaatt------???aaaggggctgggcgtggtggctcacacctgtaatccca

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
            format(alignment, "psl"),
            """\
44	1	0	0	1	3	1	6	+	hg19_dna	50	1	49	chr13	115169878	52759147	52759198	2	7,38,	1,11,	52759147,52759160,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 50)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 50))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr1        1207056 ?????????????????????????????????????????????????? 1207106
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||      50
hg19_dna          0 ??????????????????????????????????????????????????      50
""",
            )
        if fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[1207056:1207106],
                "caaaaattcacaaaggggctgggcgtggtggctcacacctgtaatcccaa",
            )
            self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
            self.assertEqual(
                str(alignment),
                """\
chr1        1207056 caaaaattcacaaaggggctgggcgtggtggctcacacctgtaatcccaa 1207106
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||      50
hg19_dna          0 caaaaattcacaaaggggctgggcgtggtggctcacacctgtaatcccaa      50
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['c', 'a', 'a', 'a', 'a', 'a', 't', 't', 'c', 'a', 'c', 'a', 'a',
           'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g', 'c', 'g', 't',
           'g', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a', 'c', 'a', 'c', 'c',
           't', 'g', 't', 'a', 'a', 't', 'c', 'c', 'c', 'a', 'a'],
          ['c', 'a', 'a', 'a', 'a', 'a', 't', 't', 'c', 'a', 'c', 'a', 'a',
           'a', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g', 'c', 'g', 't',
           'g', 'g', 't', 'g', 'g', 'c', 't', 'c', 'a', 'c', 'a', 'c', 'c',
           't', 'g', 't', 'a', 'a', 't', 'c', 'c', 'c', 'a', 'a']],
         dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
50	0	0	0	0	0	0	0	+	hg19_dna	50	0	50	chr1	249250621	1207056	1207106	1	50,	0,	1207056,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 31)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr1       61700837 ?????????????????????????????????? 61700871
                  0 ||||||||||||||||||||||||||||||||||       34
hg19_dna          1 ??????????????????????????????????       35
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[61700837:61700871],
                "aaaaatgaacaaaggggctgggcgcggtggctca",
            )
            self.assertEqual(
                alignment.query.seq[1:35], self.queries[alignment.query.id][1:35]
            )
            self.assertEqual(
                str(alignment),
                """\
chr1       61700837 aaaaatgaacaaaggggctgggcgcggtggctca 61700871
                  0 ||||||..||||||||||||||||.|||||||||       34
hg19_dna          1 aaaaattcacaaaggggctgggcgtggtggctca       35
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['a', 'a', 'a', 'a', 'a', 't', 'g', 'a', 'a', 'c', 'a', 'a', 'a',
           'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g', 'c', 'g', 'c', 'g',
           'g', 't', 'g', 'g', 'c', 't', 'c', 'a'],
          ['a', 'a', 'a', 'a', 'a', 't', 't', 'c', 'a', 'c', 'a', 'a', 'a',
           'g', 'g', 'g', 'g', 'c', 't', 'g', 'g', 'g', 'c', 'g', 't', 'g',
           'g', 't', 'g', 'g', 'c', 't', 'c', 'a']], dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
31	3	0	0	0	0	0	0	+	hg19_dna	50	1	35	chr1	249250621	61700837	61700871	1	34,	1,	61700837,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 28)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 44))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 191154276)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr4       37558157 ????????????????----------?????????????????? 37558191
                  0 ||||||||||----------------||||||||||||||||||       44
hg19_dna         49 ??????????------????????????????????????????       11
""",
            )
        elif fmt == "pslx":
            self.assertEqual(alignment.target.seq[37558157:37558167], "tgggattaca")
            self.assertEqual(
                alignment.target.seq[37558173:37558191], "accacgcccagccccttt"
            )
            self.assertEqual(
                alignment.query.seq[11:29], self.queries[alignment.query.id][11:29]
            )
            self.assertEqual(
                alignment.query.seq[39:49], self.queries[alignment.query.id][39:49]
            )
            self.assertEqual(
                str(alignment),
                """\
chr4       37558157 tgggattaca??????----------accacgcccagccccttt 37558191
                  0 ||||||||||----------------||||||||||||||||||       44
hg19_dna         49 tgggattaca------??????????accacgcccagccccttt       11
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
        self.assertEqual(alignment.matches, 35)
        self.assertEqual(alignment.misMatches, 2)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 37))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr22")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 51304566)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr22      48997405 ????????????????????????????????????? 48997442
                  0 |||||||||||||||||||||||||||||||||||||       37
hg19_dna         49 ?????????????????????????????????????       12
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[48997405:48997442],
                "tgggattacaggcgggagccaccacgcccagcccctt",
            )
            self.assertEqual(
                alignment.query.seq[12:49], self.queries[alignment.query.id][12:49]
            )
            self.assertEqual(
                str(alignment),
                """\
chr22      48997405 tgggattacaggcgggagccaccacgcccagcccctt 48997442
                  0 ||||||||||||.|.||||||||||||||||||||||       37
hg19_dna         49 tgggattacaggtgtgagccaccacgcccagcccctt       12
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 'c',
           'g', 'g', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
           'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't'],
          ['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
           'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
           'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't']],
         dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
35	2	0	0	0	0	0	0	-	hg19_dna	50	12	49	chr22	51304566	48997405	48997442	1	37,	1,	48997405,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 35)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr2      120641740 ???????????????????????????????????? 120641776
                  0 ||||||||||||||||||||||||||||||||||||        36
hg19_dna         49 ????????????????????????????????????        13
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[120641740:120641776],
                "tgggattacaggcgtgagccaccacgcccagcccct",
            )
            self.assertEqual(
                alignment.query.seq[13:49], self.queries[alignment.query.id][13:49]
            )
            self.assertEqual(
                str(alignment),
                """\
chr2      120641740 tgggattacaggcgtgagccaccacgcccagcccct 120641776
                  0 ||||||||||||.|||||||||||||||||||||||        36
hg19_dna         49 tgggattacaggtgtgagccaccacgcccagcccct        13
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 'c',
           'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
           'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't'],
          ['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
           'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
           'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't']], dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
35	1	0	0	0	0	0	0	-	hg19_dna	50	13	49	chr2	243199373	120641740	120641776	1	36,	1,	120641740,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 39)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr19      54017130 ??????????????????????????????????????? 54017169
                  0 |||||||||||||||||||||||||||||||||||||||       39
hg19_dna         49 ???????????????????????????????????????       10
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[54017130:54017169],
                "tgggattacaggtgtgagccaccacgcccagcccctttg",
            )
            self.assertEqual(
                alignment.query.seq[10:49], self.queries[alignment.query.id][10:49]
            )
            self.assertEqual(
                str(alignment),
                """\
chr19      54017130 tgggattacaggtgtgagccaccacgcccagcccctttg 54017169
                  0 |||||||||||||||||||||||||||||||||||||||       39
hg19_dna         49 tgggattacaggtgtgagccaccacgcccagcccctttg       10
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
           'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
           'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't', 't', 'g'],
          ['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
           'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
           'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't', 't', 'g']],
         dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
39	0	0	0	0	0	0	0	-	hg19_dna	50	10	49	chr19	59128983	54017130	54017169	1	39,	1,	54017130,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 36)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr19        553742 ??????????????????????????????????????? 553781
                  0 |||||||||||||||||||||||||||||||||||||||     39
hg19_dna         49 ???????????????????????????????????????     10
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[553742:553781],
                "tgggatgacaggggtgaggcaccacgcccagcccctttg",
            )
            self.assertEqual(
                alignment.query.seq[10:49], self.queries[alignment.query.id][10:49]
            )
            self.assertEqual(
                str(alignment),
                """\
chr19        553742 tgggatgacaggggtgaggcaccacgcccagcccctttg 553781
                  0 ||||||.|||||.|||||.||||||||||||||||||||     39
hg19_dna         49 tgggattacaggtgtgagccaccacgcccagcccctttg     10
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['t', 'g', 'g', 'g', 'a', 't', 'g', 'a', 'c', 'a', 'g', 'g', 'g',
           'g', 't', 'g', 'a', 'g', 'g', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
           'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't', 't', 'g'],
          ['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
           'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
           'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't', 't', 'g']],
         dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
36	3	0	0	0	0	0	0	-	hg19_dna	50	10	49	chr19	59128983	553742	553781	1	39,	1,	553742,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 33)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr10")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 135534747)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr10      99388555 ???????????????????????????????????? 99388591
                  0 ||||||||||||||||||||||||||||||||||||       36
hg19_dna         49 ????????????????????????????????????       13
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[99388555:99388591],
                "tgggattataggcatgagccaccacgcccagcccct",
            )
            self.assertEqual(
                alignment.query.seq[13:49], self.queries[alignment.query.id][13:49]
            )
            self.assertEqual(
                str(alignment),
                """\
chr10      99388555 tgggattataggcatgagccaccacgcccagcccct 99388591
                  0 ||||||||.|||..||||||||||||||||||||||       36
hg19_dna         49 tgggattacaggtgtgagccaccacgcccagcccct       13
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 't', 'a', 'g', 'g', 'c',
           'a', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
           'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't'],
          ['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
           'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
           'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't']], dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
33	3	0	0	0	0	0	0	-	hg19_dna	50	13	49	chr10	135534747	99388555	99388591	1	36,	1,	99388555,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 24)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 25))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr10")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 135534747)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr10     112178171 ????????????????????????? 112178196
                  0 |||||||||||||||||||||||||        25
hg19_dna         35 ?????????????????????????        10
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[112178171:112178196], "tgagtcaccacgcccagcccctttg"
            )
            self.assertEqual(
                alignment.query.seq[10:35], self.queries[alignment.query.id][10:35]
            )
            self.assertEqual(
                str(alignment),
                """\
chr10     112178171 tgagtcaccacgcccagcccctttg 112178196
                  0 ||||.||||||||||||||||||||        25
hg19_dna         35 tgagccaccacgcccagcccctttg        10
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['t', 'g', 'a', 'g', 't', 'c', 'a', 'c', 'c', 'a', 'c', 'g', 'c',
           'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't', 't', 'g'],
          ['t', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g', 'c',
           'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't', 't', 't', 'g']],
         dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
24	1	0	0	0	0	0	0	-	hg19_dna	50	10	35	chr10	135534747	112178171	112178196	1	25,	15,	112178171,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 35)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr1       39368490 ???????????????????????????????????? 39368526
                  0 ||||||||||||||||||||||||||||||||||||       36
hg19_dna         49 ????????????????????????????????????       13
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[39368490:39368526],
                "tgggattacaggcgtgagccaccacgcccagcccct",
            )
            self.assertEqual(
                alignment.query.seq[13:49], self.queries[alignment.query.id][13:49]
            )
            self.assertEqual(
                str(alignment),
                """\
chr1       39368490 tgggattacaggcgtgagccaccacgcccagcccct 39368526
                  0 ||||||||||||.|||||||||||||||||||||||       36
hg19_dna         49 tgggattacaggtgtgagccaccacgcccagcccct       13
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 'c',
           'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
           'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't'],
          ['t', 'g', 'g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't',
           'g', 't', 'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g',
           'c', 'c', 'c', 'a', 'g', 'c', 'c', 'c', 'c', 't']], dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
35	1	0	0	0	0	0	0	-	hg19_dna	50	13	49	chr1	249250621	39368490	39368526	1	36,	1,	39368490,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 33)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        if fmt == "psl":
            self.assertEqual(
                str(alignment),
                """\
chr1      220325687 ?????????????????????????????????? 220325721
                  0 ||||||||||||||||||||||||||||||||||        34
hg19_dna         47 ??????????????????????????????????        13
""",
            )
        elif fmt == "pslx":
            self.assertEqual(
                alignment.target.seq[220325687:220325721],
                "ggattacaggcgtgagccaccacgcccagcccct",
            )
            self.assertEqual(
                alignment.query.seq[13:47], self.queries[alignment.query.id][13:47]
            )
            self.assertEqual(
                str(alignment),
                """\
chr1      220325687 ggattacaggcgtgagccaccacgcccagcccct 220325721
                  0 ||||||||||.|||||||||||||||||||||||        34
hg19_dna         47 ggattacaggtgtgagccaccacgcccagcccct        13
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 'c', 'g', 't',
           'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g', 'c', 'c',
           'c', 'a', 'g', 'c', 'c', 'c', 'c', 't'],
          ['g', 'g', 'a', 't', 't', 'a', 'c', 'a', 'g', 'g', 't', 'g', 't',
           'g', 'a', 'g', 'c', 'c', 'a', 'c', 'c', 'a', 'c', 'g', 'c', 'c',
           'c', 'a', 'g', 'c', 'c', 'c', 'c', 't']], dtype='U')
                    # fmt: on
                )
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
            format(alignment, "psl"),
            """\
33	1	0	0	0	0	0	0	-	hg19_dna	50	13	47	chr1	249250621	220325687	220325721	1	34,	3,	220325687,
""",
        )
        self.assertRaises(StopIteration, next, alignments)

    def test_writing_psl_34_005(self):
        """Test writing the alignments in psl_34_005.psl."""
        path = "Blat/psl_34_005.psl"
        with open(path) as stream:
            original_data = stream.read()
        alignments = Align.parse(path, "psl")
        stream = StringIO()
        n = Align.write(alignments, stream, "psl", header=False)
        self.assertEqual(n, 22)
        stream.seek(0)
        written_data = stream.read()
        stream.close()
        self.assertEqual(original_data, written_data)


class TestAlign_dnax_prot(unittest.TestCase):
    @classmethod
    def read_dna(cls, assembly, sequence):
        path = "Blat/%s.fa" % assembly
        records = SeqIO.parse(path, "fasta")
        for record in records:
            name, start_end = record.id.split(":")
            if name == sequence.id:
                break
        else:
            raise Exception("Failed to find DNA sequence")
        start, end = start_end.split("-")
        start = int(start)
        end = int(end)
        length = len(sequence)
        sequence = str(record.seq)
        dna = Seq({start: sequence}, length=length)
        return dna

    def test_reading_psl_35_001(self):
        """Test parsing psl_35_001.psl and pslx_35_001.pslx."""
        self.check_reading_psl_35_001("psl")
        self.check_reading_psl_35_001("pslx")

    def check_reading_psl_35_001(self, fmt):
        """Check parsing psl_35_001.psl or pslx_35_001.pslx."""
        path = "Blat/%s_35_001.%s" % (fmt, fmt)
        alignments = Align.parse(path, "psl")
        self.assertEqual(alignments.metadata["psLayout version"], "3")
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 52)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertEqual(len(alignment.target.seq), 114364328)
        self.assertEqual(len(alignment.query.seq), 230)
        if fmt == "pslx":
            self.assertEqual(
                alignment.query.seq[61:113],
                "YEVFRTEEEEKIKSQGQDVTSSVYFMKQTISNACGTIGLIHAIANNKDKMHF",
            )
            self.assertEqual(len(alignment.target.seq.defined_ranges), 0)
            self.assertEqual(len(alignment.target.features), 1)
            feature = alignment.target.features[0]
            self.assertEqual(
                feature.location,
                SimpleLocation(
                    ExactPosition(75566694), ExactPosition(75566850), strand=+1
                ),
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(len(feature.qualifiers["translation"]), 1)
            self.assertEqual(
                feature.qualifiers["translation"][0],
                "YEVFRTEEEEKIKSQGQDVTSSVYFMKQTISNACGTIGLIHAIANNKDKMHF",
            )
            alignment.target.seq = TestAlign_dnax_prot.read_dna(
                "hg38", alignment.target
            )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[75566694, 75566850],
                          [      61,      113]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
52	0	0	0	0	0	0	0	++	CAG33136.1	230	61	113	chr13	114364328	75566694	75566850	1	52,	61,	75566694,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 44)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertEqual(len(alignment.target.seq), 114364328)
        self.assertEqual(len(alignment.query.seq), 230)
        if fmt == "pslx":
            self.assertEqual(
                alignment.query.seq[17:61],
                "QFLKQLGLHPNWQFVDVYGMDPELLSMVPRPVCAVLLLFPITEK",
            )
            self.assertEqual(len(alignment.target.seq.defined_ranges), 0)
            self.assertEqual(len(alignment.target.features), 1)
            feature = alignment.target.features[0]
            self.assertEqual(
                feature.location,
                SimpleLocation(
                    ExactPosition(75560749), ExactPosition(75560881), strand=+1
                ),
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(len(feature.qualifiers["translation"]), 1)
            self.assertEqual(
                feature.qualifiers["translation"][0],
                "QFLKQLGLHPNWQFVDVYGMDPELLSMVPRPVCAVLLLFPITEK",
            )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[75560749, 75560881],
                          [      17,       61]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
44	0	0	0	0	0	0	0	++	CAG33136.1	230	17	61	chr13	114364328	75560749	75560881	1	44,	17,	75560749,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 44)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertEqual(len(alignment.target.seq), 114364328)
        self.assertEqual(len(alignment.query.seq), 230)
        if fmt == "pslx":
            self.assertEqual(alignment.query.seq[:15], "MEGQRWLPLEANPEV")
            self.assertEqual(
                alignment.query.seq[113:142], "ESGSTLKKFLEESVSMSPEERARYLENYD"
            )
            self.assertEqual(len(alignment.target.seq.defined_ranges), 0)
            self.assertEqual(len(alignment.target.features), 1)
            feature = alignment.target.features[0]
            self.assertEqual(
                feature.location,
                CompoundLocation(
                    [
                        SimpleLocation(
                            ExactPosition(75549820), ExactPosition(75549865), strand=+1
                        ),
                        SimpleLocation(
                            ExactPosition(75567225), ExactPosition(75567312), strand=+1
                        ),
                    ],
                    operator="join",
                ),
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(len(feature.qualifiers["translation"]), 1)
            self.assertEqual(
                feature.qualifiers["translation"][0],
                "MEGQRWLPLEANPEVESGSTLKKFLEESVSMSPEERARYLENYD",
            )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[75549820, 75549865, 75567225, 75567225, 75567312],
                          [       0,       15,       15,      113,      142]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
44	0	0	0	1	98	1	17360	++	CAG33136.1	230	0	142	chr13	114364328	75549820	75567312	2	15,29,	0,113,	75549820,75567225,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 47)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertEqual(len(alignment.target.seq), 114364328)
        self.assertEqual(len(alignment.query.seq), 230)
        if fmt == "pslx":
            self.assertEqual(
                alignment.query.seq[183:],
                "DGRKPFPINHGETSDETLLEDAIEVCKKFMERDPDELRFNAIALSAA",
            )
            self.assertEqual(len(alignment.target.seq.defined_ranges), 0)
            self.assertEqual(len(alignment.target.features), 1)
            feature = alignment.target.features[0]
            self.assertEqual(
                feature.location,
                CompoundLocation(
                    [
                        SimpleLocation(
                            ExactPosition(75604767), ExactPosition(75604827), strand=+1
                        ),
                        SimpleLocation(
                            ExactPosition(75605728), ExactPosition(75605809), strand=+1
                        ),
                    ],
                    operator="join",
                ),
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(len(feature.qualifiers["translation"]), 1)
            self.assertEqual(
                feature.qualifiers["translation"][0],
                "DGRKPFPINHGETSDETLLEDAIEVCKKFMERDPDELRFNAIALSAA",
            )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[75604767, 75604827, 75605728, 75605809],
                          [     183,      203,      203,      230]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
47	0	0	0	0	0	1	901	++	CAG33136.1	230	183	230	chr13	114364328	75604767	75605809	2	20,27,	183,203,	75604767,75605728,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 25)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertEqual(len(alignment.target.seq), 114364328)
        self.assertEqual(len(alignment.query.seq), 230)
        if fmt == "pslx":
            self.assertEqual(alignment.query.seq[158:183], "APSIDEKVDLHFIALVHVDGHLYEL")
            self.assertEqual(len(alignment.target.seq.defined_ranges), 0)
            self.assertEqual(len(alignment.target.features), 1)
            feature = alignment.target.features[0]
            self.assertEqual(
                feature.location,
                SimpleLocation(
                    ExactPosition(75594914), ExactPosition(75594989), strand=+1
                ),
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(len(feature.qualifiers["translation"]), 1)
            self.assertEqual(
                feature.qualifiers["translation"][0], "APSIDEKVDLHFIALVHVDGHLYEL"
            )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[75594914, 75594989],
                          [     158,      183]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
25	0	0	0	0	0	0	0	++	CAG33136.1	230	158	183	chr13	114364328	75594914	75594989	1	25,	158,	75594914,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 16)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertEqual(len(alignment.target.seq), 114364328)
        self.assertEqual(len(alignment.query.seq), 230)
        if fmt == "pslx":
            self.assertEqual(alignment.query.seq[142:158], "AIRVTHETSAHEGQTE")
            self.assertEqual(len(alignment.target.seq.defined_ranges), 0)
            self.assertEqual(len(alignment.target.features), 1)
            feature = alignment.target.features[0]
            self.assertEqual(
                feature.location,
                SimpleLocation(
                    ExactPosition(75569459), ExactPosition(75569507), strand=+1
                ),
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(len(feature.qualifiers["translation"]), 1)
            self.assertEqual(feature.qualifiers["translation"][0], "AIRVTHETSAHEGQTE")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[75569459, 75569507],
                          [     142,      158]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
16	0	0	0	0	0	0	0	++	CAG33136.1	230	142	158	chr13	114364328	75569459	75569507	1	16,	142,	75569459,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 26)
        self.assertEqual(alignment.misMatches, 8)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertEqual(len(alignment.target.seq), 190214555)
        self.assertEqual(len(alignment.query.seq), 230)
        if fmt == "pslx":
            self.assertEqual(
                alignment.query.seq[76:110], "GQDVTSSVYFMKQTISNACGTIGLIHAIANNKDK"
            )
            self.assertEqual(len(alignment.target.seq.defined_ranges), 0)
            self.assertEqual(len(alignment.target.features), 1)
            feature = alignment.target.features[0]
            self.assertEqual(
                feature.location,
                SimpleLocation(
                    ExactPosition(41260685), ExactPosition(41260787), strand=+1
                ),
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(len(feature.qualifiers["translation"]), 1)
            self.assertEqual(
                feature.qualifiers["translation"][0],
                "GQEVSPKVYFMKQTIGNSCGTIGLIHAVANNQDK",
            )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[41260685, 41260787],
                          [      76,      110]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
26	8	0	0	0	0	0	0	++	CAG33136.1	230	76	110	chr4	190214555	41260685	41260787	1	34,	76,	41260685,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 37)
        self.assertEqual(alignment.misMatches, 26)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertEqual(len(alignment.target.seq), 190214555)
        self.assertEqual(len(alignment.query.seq), 230)
        if fmt == "pslx":
            self.assertEqual(
                alignment.query.seq[17:59], "QFLKQLGLHPNWQFVDVYGMDPELLSMVPRPVCAVLLLFPIT"
            )
            self.assertEqual(alignment.query.seq[162:183], "DEKVDLHFIALVHVDGHLYEL")
            self.assertEqual(len(alignment.target.seq.defined_ranges), 0)
            self.assertEqual(len(alignment.target.features), 1)
            feature = alignment.target.features[0]
            self.assertEqual(
                feature.location,
                CompoundLocation(
                    [
                        SimpleLocation(
                            ExactPosition(41257605), ExactPosition(41257731), strand=+1
                        ),
                        SimpleLocation(
                            ExactPosition(41263227), ExactPosition(41263290), strand=+1
                        ),
                    ],
                    operator="join",
                ),
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(len(feature.qualifiers["translation"]), 1)
            self.assertEqual(
                feature.qualifiers["translation"][0],
                "QVLSRLGVAGQWRFVDVLGLEEESLGSVPAPACALLLLFPLTDDKVNFHFILFNNVDGHLYEL",
            )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[41257605, 41257731, 41263227, 41263227, 41263290],
                          [      17,       59,       59,      162,      183]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
37	26	0	0	1	103	1	5496	++	CAG33136.1	230	17	183	chr4	190214555	41257605	41263290	2	42,21,	17,162,	41257605,41263227,
""",
        )
        self.assertRaises(StopIteration, next, alignments)

    def test_writing_psl_35_001(self):
        """Test writing the alignments in psl_35_001.psl."""
        path = "Blat/psl_35_001.psl"
        with open(path) as stream:
            original_data = stream.read()
        alignments = Align.parse(path, "psl")
        stream = StringIO()
        n = Align.write(alignments, stream, "psl")
        self.assertEqual(n, 8)
        stream.seek(0)
        written_data = stream.read()
        stream.close()
        self.assertEqual(original_data, written_data)
        # Convert the alignment to a protein alignment and insert the
        # appropriate sequence data. Write this alignment in a PSL file;
        # the writer will recalculate the values for matches, misMatches,
        # repMatches, and nCount from the sequence data and the alignment.
        #
        # The alignments were generated using
        # blat -t=dnax -q=prot hg38.2bit CAG33136.1.fasta psl_35_001.psl
        #
        # To save disk space, we extracted the necessary sequence data using
        #
        # twoBitToFa hg38.2bit:chr13:75549820-75605809 stdout
        # twoBitToFa hg38.2bit:chr4:41257605-41263290 stdout
        #
        # and concatenating the results into file hg38.fa. We will use this
        # file below, and create partially defined Seq objects.
        #
        # Load the protein sequence:
        protein = SeqIO.read("Blat/CAG33136.1.fasta", "fasta")
        protein_alignments = []
        alignments = Align.parse(path, "psl")
        for i, alignment in enumerate(alignments):
            alignment.sequences[0].seq = TestAlign_dnax_prot.read_dna(
                "hg38", alignment.sequences[0]
            )
            self.assertEqual(alignment.sequences[1].id, protein.id)
            alignment.sequences[1].seq = protein.seq
            # The alignment is on the forward strand of the DNA sequence:
            self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
            # The protein alignment is also in the forward orientation:
            self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
            # Now extract the aligned sequences:
            aligned_dna = ""
            aligned_protein = ""
            for start, end in alignment.aligned[0]:
                aligned_dna += alignment.sequences[0].seq[start:end]
            for start, end in alignment.aligned[1]:
                aligned_protein += alignment.sequences[1].seq[start:end]
            # Translate the aligned DNA sequence:
            aligned_dna = Seq(aligned_dna)
            aligned_dna_translated = Seq(aligned_dna.translate())
            aligned_protein = Seq(aligned_protein)
            # Create a new alignment including the aligned sequences only:
            records = [
                SeqRecord(aligned_dna_translated, id=alignment.sequences[0].id),
                SeqRecord(aligned_protein, id=alignment.sequences[1].id),
            ]
            coordinates = np.array(
                [[0, len(aligned_dna_translated)], [0, len(aligned_protein)]]
            )
            protein_alignment = Alignment(records, coordinates)
            protein_alignments.append(protein_alignment)
            if i == 0:
                self.assertEqual(
                    str(protein_alignment),
                    """\
chr13             0 YEVFRTEEEEKIKSQGQDVTSSVYFMKQTISNACGTIGLIHAIANNKDKMHF 52
                  0 |||||||||||||||||||||||||||||||||||||||||||||||||||| 52
CAG33136.         0 YEVFRTEEEEKIKSQGQDVTSSVYFMKQTISNACGTIGLIHAIANNKDKMHF 52
""",
                )
            elif i == 1:
                self.assertEqual(
                    str(protein_alignment),
                    """\
chr13             0 QFLKQLGLHPNWQFVDVYGMDPELLSMVPRPVCAVLLLFPITEK 44
                  0 |||||||||||||||||||||||||||||||||||||||||||| 44
CAG33136.         0 QFLKQLGLHPNWQFVDVYGMDPELLSMVPRPVCAVLLLFPITEK 44
""",
                )
            elif i == 2:
                self.assertEqual(
                    str(protein_alignment),
                    """\
chr13             0 MEGQRWLPLEANPEVESGSTLKKFLEESVSMSPEERARYLENYD 44
                  0 |||||||||||||||||||||||||||||||||||||||||||| 44
CAG33136.         0 MEGQRWLPLEANPEVESGSTLKKFLEESVSMSPEERARYLENYD 44
""",
                )
            elif i == 3:
                self.assertEqual(
                    str(protein_alignment),
                    """\
chr13             0 DGRKPFPINHGETSDETLLEDAIEVCKKFMERDPDELRFNAIALSAA 47
                  0 ||||||||||||||||||||||||||||||||||||||||||||||| 47
CAG33136.         0 DGRKPFPINHGETSDETLLEDAIEVCKKFMERDPDELRFNAIALSAA 47
""",
                )
            elif i == 4:
                self.assertEqual(
                    str(protein_alignment),
                    """\
chr13             0 APSIDEKVDLHFIALVHVDGHLYEL 25
                  0 ||||||||||||||||||||||||| 25
CAG33136.         0 APSIDEKVDLHFIALVHVDGHLYEL 25
""",
                )
            elif i == 5:
                self.assertEqual(
                    str(protein_alignment),
                    """\
chr13             0 AIRVTHETSAHEGQTE 16
                  0 |||||||||||||||| 16
CAG33136.         0 AIRVTHETSAHEGQTE 16
""",
                )
            elif i == 6:
                self.assertEqual(
                    str(protein_alignment),
                    """\
chr4              0 GQEVSPKVYFMKQTIGNSCGTIGLIHAVANNQDK 34
                  0 ||.|...||||||||.|.|||||||||.|||.|| 34
CAG33136.         0 GQDVTSSVYFMKQTISNACGTIGLIHAIANNKDK 34
""",
                )
            elif i == 7:
                self.assertEqual(
                    str(protein_alignment),
                    """\
chr4              0 QVLSRLGVAGQWRFVDVLGLEEESLGSVPAPACALLLLFPLTDDKVNFHFILFNNVDGHL
                  0 |.|..||....|.||||.|...|.|..||.|.||.|||||.||.||..|||....|||||
CAG33136.         0 QFLKQLGLHPNWQFVDVYGMDPELLSMVPRPVCAVLLLFPITDEKVDLHFIALVHVDGHL

chr4             60 YEL 63
                 60 ||| 63
CAG33136.        60 YEL 63
""",
                )
        # Write the protein alignments to a PSL file:
        stream = StringIO()
        n = Align.write(protein_alignments, stream, "psl", wildcard="X")
        self.assertEqual(n, 8)
        # Read the alignments back in:
        alignments = Align.parse(path, "psl")
        stream.seek(0)
        protein_alignments = Align.parse(stream, "psl")
        for alignment, protein_alignment in zip(alignments, protein_alignments):
            # Confirm that the recalculated values for matches, misMatches,
            # repMatches, and nCount are correct:
            self.assertEqual(alignment.matches, protein_alignment.matches)
            self.assertEqual(alignment.misMatches, protein_alignment.misMatches)
            self.assertEqual(alignment.repMatches, protein_alignment.repMatches)
            self.assertEqual(alignment.nCount, protein_alignment.nCount)

    def test_reading_psl_35_002(self):
        """Test parsing psl_35_002.psl."""
        # See below for a description of the file balAcu1.fa.
        # We use this file here so we can check the SeqFeatures.
        records = SeqIO.parse("Blat/balAcu1.fa", "fasta")
        self.dna = {}
        for record in records:
            name, start_end = record.id.split(":")
            start, end = start_end.split("-")
            start = int(start)
            end = int(end)
            sequence = str(record.seq)
            self.dna[name] = Seq({start: sequence}, length=end)
        self.check_reading_psl_35_002("psl")
        self.check_reading_psl_35_002("pslx")

    def check_reading_psl_35_002(self, fmt):
        """Check parsing psl_35_002.psl or pslx_35_002.pslx."""
        path = "Blat/%s_35_002.%s" % (fmt, fmt)
        alignments = Align.parse(path, "psl")
        self.assertEqual(alignments.metadata["psLayout version"], "3")
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 210)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "KI537979")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertEqual(len(alignment.target.seq), 14052872)
        self.assertEqual(len(alignment.query.seq), 230)
        if fmt == "pslx":
            self.assertEqual(
                alignment.query.seq[17:],
                "QFLKQLGLHPNWQFVDVYGMDPELLSMVPRPVCAVLLLFPITEKYEVFRTEEEEKIKSQGQDVTSSVYFMKQTISNACGTIGLIHAIANNKDKMHFESGSTLKKFLEESVSMSPEERARYLENYDAIRVTHETSAHEGQTEAPSIDEKVDLHFIALVHVDGHLYELDGRKPFPINHGETSDETLLEDAIEVCKKFMERDPDELRFNAIALSAA",
            )
            self.assertEqual(len(alignment.target.seq.defined_ranges), 0)
            self.assertEqual(len(alignment.target.features), 1)
            feature = alignment.target.features[0]
            self.assertEqual(
                feature.location,
                CompoundLocation(
                    [
                        SimpleLocation(
                            ExactPosition(9712654), ExactPosition(9712786), strand=+1
                        ),
                        SimpleLocation(
                            ExactPosition(9715941), ExactPosition(9716097), strand=+1
                        ),
                        SimpleLocation(
                            ExactPosition(9716445), ExactPosition(9716532), strand=+1
                        ),
                        SimpleLocation(
                            ExactPosition(9718374), ExactPosition(9718422), strand=+1
                        ),
                        SimpleLocation(
                            ExactPosition(9739264), ExactPosition(9739339), strand=+1
                        ),
                        SimpleLocation(
                            ExactPosition(9743706), ExactPosition(9743766), strand=+1
                        ),
                        SimpleLocation(
                            ExactPosition(9744511), ExactPosition(9744592), strand=+1
                        ),
                    ],
                    operator="join",
                ),
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(len(feature.qualifiers["translation"]), 1)
            self.assertEqual(
                feature.qualifiers["translation"][0],
                "QFLKQLGLHPNWQFVDVYGMDPELLSMVPRPVCAVLLLFPITEKYEIFRTEEEEKIKSQGQDVTSSVYFMKQTISNACGTIGLIHAIANNKDKMHFESGSTLKKFLEESASMSPEERARYLENYDAIRVTHETSAHEGQTEAPNIDEKVDLHFIALVHVDGHLYELDGRKPFPINHGETSDETLLEDAIEVCKKFMERDPDELRFNAIALSAA",
            )
            # confirm that the feature coordinates are correct by extracting
            # the feature sequence from the target sequence and translating it.
            cds = feature.extract(self.dna[alignment.target.id]).translate()
            self.assertEqual(feature.qualifiers["translation"][0], cds)
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[9712654, 9712786, 9715941, 9716097, 9716445, 9716532, 9718374,
                           9718422, 9739264, 9739339, 9743706, 9743766, 9744511, 9744592],
                          [     17,      61,      61,     113,     113,     142,     142,
                               158,     158,     183,     183,     203,     203,     230]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
210	3	0	0	0	0	6	31299	++	CAG33136.1	230	17	230	KI537979	14052872	9712654	9744592	7	44,52,29,16,25,20,27,	17,61,113,142,158,183,203,	9712654,9715941,9716445,9718374,9739264,9743706,9744511,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 207)
        self.assertEqual(alignment.misMatches, 22)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "KI538594")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertEqual(len(alignment.target.seq), 7819582)
        self.assertEqual(len(alignment.query.seq), 230)
        if fmt == "pslx":
            self.assertEqual(alignment.query.seq[:20], "MEGQRWLPLEANPEVTNQFL")
            self.assertEqual(
                alignment.query.seq[21:],
                "QLGLHPNWQFVDVYGMDPELLSMVPRPVCAVLLLFPITEKYEVFRTEEEEKIKSQGQDVTSSVYFMKQTISNACGTIGLIHAIANNKDKMHFESGSTLKKFLEESVSMSPEERARYLENYDAIRVTHETSAHEGQTEAPSIDEKVDLHFIALVHVDGHLYELDGRKPFPINHGETSDETLLEDAIEVCKKFMERDPDELRFNAIALSAA",
            )
            self.assertEqual(len(alignment.target.seq.defined_ranges), 0)
            self.assertEqual(len(alignment.target.features), 1)
            feature = alignment.target.features[0]
            self.assertEqual(
                feature.location,
                CompoundLocation(
                    [
                        SimpleLocation(
                            ExactPosition(2103463), ExactPosition(2103523), strand=+1
                        ),
                        SimpleLocation(
                            ExactPosition(2103522), ExactPosition(2104149), strand=+1
                        ),
                    ],
                    operator="join",
                ),
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(len(feature.qualifiers["translation"]), 1)
            self.assertEqual(
                feature.qualifiers["translation"][0],
                "MEGQCWLPLEANPEVTNQLLQLGLHPNWQFVDVYGMDPELLSMVPRPVCAVLLLFPITEKYEVFRTEEEEKIKSQGQNITSSGYFMRQTISSACGTIGLIHAIANNKDKMHFESGSTLKKFLEESASLSPEERAIYLENYDSIRVTHKTSDHEGQTEAQNIDEKVDLHFIALVHVDGHLYELDGWKPFPINHGETSDATLLRDAIEVFKKFRERDPDERRFNVIALSAA",
            )
            # confirm that the feature coordinates are correct by extracting
            # the feature sequence from the target sequence and translating it.
            cds = feature.extract(self.dna[alignment.target.id]).translate()
            self.assertEqual(feature.qualifiers["translation"][0], cds)
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[2103463, 2103523, 2103522, 2103522, 2104149],
                          [      0,      20,      20,      21,     230]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
207	22	0	0	1	1	1	-1	++	CAG33136.1	230	0	230	KI538594	7819582	2103463	2104149	2	20,209,	0,21,	2103463,2103522,
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 204)
        self.assertEqual(alignment.misMatches, 6)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "KI537194")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertEqual(len(alignment.target.seq), 37111980)
        self.assertEqual(len(alignment.query.seq), 230)
        if fmt == "pslx":
            self.assertEqual(
                alignment.query.seq[:183],
                "MEGQRWLPLEANPEVTNQFLKQLGLHPNWQFVDVYGMDPELLSMVPRPVCAVLLLFPITEKYEVFRTEEEEKIKSQGQDVTSSVYFMKQTISNACGTIGLIHAIANNKDKMHFESGSTLKKFLEESVSMSPEERARYLENYDAIRVTHETSAHEGQTEAPSIDEKVDLHFIALVHVDGHLYEL",
            )
            self.assertEqual(alignment.query.seq[203:], "DAIEVCKKFMERDPDELRFNAIALSAA")
            self.assertEqual(len(alignment.target.seq.defined_ranges), 0)
            self.assertEqual(len(alignment.target.features), 1)
            feature = alignment.target.features[0]
            self.assertEqual(
                feature.location,
                CompoundLocation(
                    [
                        SimpleLocation(
                            ExactPosition(20872472), ExactPosition(20873021), strand=-1
                        ),
                        SimpleLocation(
                            ExactPosition(20872390), ExactPosition(20872471), strand=-1
                        ),
                    ],
                    operator="join",
                ),
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(len(feature.qualifiers["translation"]), 1)
            self.assertEqual(
                feature.qualifiers["translation"][0],
                "MESQRWLPLEANPEVTNQFLKQLGLHPNWQCVDVYGMDPELLSMVPRPVCAVLLLFPITEKYEIFRTEEEEKTKSQGQDVTSSVYFMKQTISNACGTIGLIHAIANNKDKMHFESGSTLKKFLEESASMSPEERARYLENYDAIRVTHETSAHEGQTEAPNIDEKVDLHFIALVHVDGHLYELDAIEVCKKFMERDPDELRFNAIALSAA",
            )
            # confirm that the feature coordinates are correct by extracting
            # the feature sequence from the target sequence and translating it.
            cds = feature.extract(self.dna[alignment.target.id]).translate()
            self.assertEqual(feature.qualifiers["translation"][0], cds)
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[20873021, 20872472, 20872471, 20872471, 20872390],
                          [       0,      183,      183,      203,      230]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "psl"),
            """\
204	6	0	0	1	20	1	1	+-	CAG33136.1	230	0	230	KI537194	37111980	20872390	20873021	2	183,27,	0,203,	16238959,16239509,
""",
        )
        self.assertRaises(StopIteration, next, alignments)

    def test_writing_psl_35_002(self):
        """Test writing the alignments in psl_35_002.psl."""
        path = "Blat/psl_35_002.psl"
        with open(path) as stream:
            original_data = stream.read()
        alignments = Align.parse(path, "psl")
        stream = StringIO()
        n = Align.write(alignments, stream, "psl")
        self.assertEqual(n, 3)
        stream.seek(0)
        written_data = stream.read()
        stream.close()
        self.assertEqual(original_data, written_data)
        # Convert the alignment to a protein alignment and insert the
        # appropriate sequence data. Write this alignment in a PSL file;
        # the writer will recalculate the values for matches, misMatches,
        # repMatches, and nCount from the sequence data and the alignment.
        #
        # The alignments were generated using
        # blat -t=dnax -q=prot balAcu1.2bit CAG33136.1.fasta psl_35_001.psl
        #
        # To save disk space, we extracted the necessary sequence data using
        #
        # twoBitToFa balAcu1.2bit:KI537979:9712654-9744592 stdout
        # twoBitToFa balAcu1.2bit:KI538594:2103463-2104149 stdout
        # twoBitToFa balAcu1.2bit:KI537194:20872390-20873021 stdout
        #
        # and concatenating the results into file balAcu1.fa. We will use this
        # file below, and create partially defined Seq objects.
        #
        # Load the protein sequence:
        protein = SeqIO.read("Blat/CAG33136.1.fasta", "fasta")
        protein_alignments = []
        alignments = Align.parse(path, "psl")
        for i, alignment in enumerate(alignments):
            alignment.sequences[0].seq = TestAlign_dnax_prot.read_dna(
                "balAcu1", alignment.sequences[0]
            )
            self.assertEqual(alignment.sequences[1].id, protein.id)
            alignment.sequences[1].seq = protein.seq
            if i == 0 or i == 1:
                # The alignment is on the forward strand of the DNA sequence:
                self.assertLess(
                    alignment.coordinates[0, 0], alignment.coordinates[0, -1]
                )
            elif i == 2:
                # The alignment is on the reverse strand of the DNA sequence:
                self.assertGreater(
                    alignment.coordinates[0, 0], alignment.coordinates[0, -1]
                )
                # so we take the reverse complement:
                alignment.coordinates[0, :] = (
                    len(alignment.sequences[0].seq) - alignment.coordinates[0, :]
                )
                alignment.sequences[0].seq = alignment.sequences[
                    0
                ].seq.reverse_complement()
            # The protein alignment is always in the forward orientation:
            self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
            # Now extract the aligned sequences:
            aligned_dna = ""
            aligned_protein = ""
            for start, end in alignment.aligned[0]:
                aligned_dna += alignment.sequences[0].seq[start:end]
            for start, end in alignment.aligned[1]:
                aligned_protein += alignment.sequences[1].seq[start:end]
            # Translate the aligned DNA sequence:
            aligned_dna = Seq(aligned_dna)
            aligned_dna_translated = Seq(aligned_dna.translate())
            aligned_protein = Seq(aligned_protein)
            # Create a new alignment including the aligned sequences only:
            records = [
                SeqRecord(aligned_dna_translated, id=alignment.sequences[0].id),
                SeqRecord(aligned_protein, id=alignment.sequences[1].id),
            ]
            coordinates = np.array(
                [[0, len(aligned_dna_translated)], [0, len(aligned_protein)]]
            )
            protein_alignment = Alignment(records, coordinates)
            protein_alignments.append(protein_alignment)
            if i == 0:
                self.assertEqual(
                    str(protein_alignment),
                    """\
KI537979          0 QFLKQLGLHPNWQFVDVYGMDPELLSMVPRPVCAVLLLFPITEKYEIFRTEEEEKIKSQG
                  0 ||||||||||||||||||||||||||||||||||||||||||||||.|||||||||||||
CAG33136.         0 QFLKQLGLHPNWQFVDVYGMDPELLSMVPRPVCAVLLLFPITEKYEVFRTEEEEKIKSQG

KI537979         60 QDVTSSVYFMKQTISNACGTIGLIHAIANNKDKMHFESGSTLKKFLEESASMSPEERARY
                 60 |||||||||||||||||||||||||||||||||||||||||||||||||.||||||||||
CAG33136.        60 QDVTSSVYFMKQTISNACGTIGLIHAIANNKDKMHFESGSTLKKFLEESVSMSPEERARY

KI537979        120 LENYDAIRVTHETSAHEGQTEAPNIDEKVDLHFIALVHVDGHLYELDGRKPFPINHGETS
                120 |||||||||||||||||||||||.||||||||||||||||||||||||||||||||||||
CAG33136.       120 LENYDAIRVTHETSAHEGQTEAPSIDEKVDLHFIALVHVDGHLYELDGRKPFPINHGETS

KI537979        180 DETLLEDAIEVCKKFMERDPDELRFNAIALSAA 213
                180 ||||||||||||||||||||||||||||||||| 213
CAG33136.       180 DETLLEDAIEVCKKFMERDPDELRFNAIALSAA 213
""",
                )
            elif i == 1:
                self.assertEqual(
                    str(protein_alignment),
                    """\
KI538594          0 MEGQCWLPLEANPEVTNQLLQLGLHPNWQFVDVYGMDPELLSMVPRPVCAVLLLFPITEK
                  0 ||||.|||||||||||||.|||||||||||||||||||||||||||||||||||||||||
CAG33136.         0 MEGQRWLPLEANPEVTNQFLQLGLHPNWQFVDVYGMDPELLSMVPRPVCAVLLLFPITEK

KI538594         60 YEVFRTEEEEKIKSQGQNITSSGYFMRQTISSACGTIGLIHAIANNKDKMHFESGSTLKK
                 60 |||||||||||||||||..|||.|||.||||.||||||||||||||||||||||||||||
CAG33136.        60 YEVFRTEEEEKIKSQGQDVTSSVYFMKQTISNACGTIGLIHAIANNKDKMHFESGSTLKK

KI538594        120 FLEESASLSPEERAIYLENYDSIRVTHKTSDHEGQTEAQNIDEKVDLHFIALVHVDGHLY
                120 |||||.|.||||||.||||||.|||||.||.|||||||..||||||||||||||||||||
CAG33136.       120 FLEESVSMSPEERARYLENYDAIRVTHETSAHEGQTEAPSIDEKVDLHFIALVHVDGHLY

KI538594        180 ELDGWKPFPINHGETSDATLLRDAIEVFKKFRERDPDERRFNVIALSAA 229
                180 ||||.||||||||||||.|||.|||||.|||.||||||.|||.|||||| 229
CAG33136.       180 ELDGRKPFPINHGETSDETLLEDAIEVCKKFMERDPDELRFNAIALSAA 229
""",
                )
            elif i == 2:
                self.assertEqual(
                    str(protein_alignment),
                    """\
KI537194          0 MESQRWLPLEANPEVTNQFLKQLGLHPNWQCVDVYGMDPELLSMVPRPVCAVLLLFPITE
                  0 ||.|||||||||||||||||||||||||||.|||||||||||||||||||||||||||||
CAG33136.         0 MEGQRWLPLEANPEVTNQFLKQLGLHPNWQFVDVYGMDPELLSMVPRPVCAVLLLFPITE

KI537194         60 KYEIFRTEEEEKTKSQGQDVTSSVYFMKQTISNACGTIGLIHAIANNKDKMHFESGSTLK
                 60 |||.||||||||.|||||||||||||||||||||||||||||||||||||||||||||||
CAG33136.        60 KYEVFRTEEEEKIKSQGQDVTSSVYFMKQTISNACGTIGLIHAIANNKDKMHFESGSTLK

KI537194        120 KFLEESASMSPEERARYLENYDAIRVTHETSAHEGQTEAPNIDEKVDLHFIALVHVDGHL
                120 ||||||.|||||||||||||||||||||||||||||||||.|||||||||||||||||||
CAG33136.       120 KFLEESVSMSPEERARYLENYDAIRVTHETSAHEGQTEAPSIDEKVDLHFIALVHVDGHL

KI537194        180 YELDAIEVCKKFMERDPDELRFNAIALSAA 210
                180 |||||||||||||||||||||||||||||| 210
CAG33136.       180 YELDAIEVCKKFMERDPDELRFNAIALSAA 210
""",
                )
        # Write the protein alignments to a PSL file:
        stream = StringIO()
        n = Align.write(protein_alignments, stream, "psl", wildcard="X")
        self.assertEqual(n, 3)
        # Read the alignments back in:
        alignments = Align.parse(path, "psl")
        stream.seek(0)
        protein_alignments = Align.parse(stream, "psl")
        for alignment, protein_alignment in zip(alignments, protein_alignments):
            # Confirm that the recalculated values for matches, misMatches,
            # repMatches, and nCount are correct:
            self.assertEqual(alignment.matches, protein_alignment.matches)
            self.assertEqual(alignment.misMatches, protein_alignment.misMatches)
            self.assertEqual(alignment.repMatches, protein_alignment.repMatches)
            self.assertEqual(alignment.nCount, protein_alignment.nCount)


class TestAlign_strand(unittest.TestCase):
    def test_format(self):
        """Test alignment with the target on the opposite strand."""
        sequences = ["AACAGCAGCGTGTCG", "CAGCTAGCGAA"]
        coordinates = np.array(
            [[0, 2, 2, 3, 4, 6, 6, 9, 10, 12, 15], [11, 11, 9, 8, 8, 6, 5, 2, 2, 0, 0]]
        )
        alignment = Alignment(sequences, coordinates)
        alignment.score = 8
        line = """\
8	0	0	0	1	1	2	2	-	query	11	0	9	target	15	2	12	4	1,2,3,2,	2,3,6,9,	2,4,6,10,
"""
        self.assertEqual(
            str(alignment),
            """\
target            0 AA--CAGC-AGCGTGTCG 15
                  0 ----|-||-|||-||--- 18
query            11 --TTC-GCTAGC-TG---  0
""",
        )
        self.assertEqual(format(alignment, "psl"), line)
        alignment.coordinates = alignment.coordinates[:, ::-1]
        self.assertEqual(
            str(alignment),
            """\
target           15 CGACACGCT-GCTG--TT  0
                  0 ---||-|||-||-|---- 18
query             0 ---CA-GCTAGC-GAA-- 11
""",
        )
        self.assertEqual(format(alignment, "psl"), line)
        alignment.coordinates = alignment.coordinates[:, ::-1]
        line = """\
8	0	0	0	1	1	2	2	-	query	11	2	11	target	15	3	13	4	2,3,2,1,	0,2,6,8,	3,6,9,12,
"""
        alignment = alignment.reverse_complement()
        self.assertEqual(
            str(alignment),
            """\
target            0 CGACACGCT-GCTG--TT 15
                  0 ---||-|||-||-|---- 18
query            11 ---CA-GCTAGC-GAA--  0
""",
        )
        self.assertEqual(format(alignment, "psl"), line)
        alignment.coordinates = alignment.coordinates[:, ::-1]
        self.assertEqual(
            str(alignment),
            """\
target           15 AA--CAGC-AGCGTGTCG  0
                  0 ----|-||-|||-||--- 18
query             0 --TTC-GCTAGC-TG--- 11
""",
        )
        self.assertEqual(format(alignment, "psl"), line)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
