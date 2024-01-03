# Copyright 2008 by Bartek Wilczynski.  All rights reserved.
# Revisions copyright 2019 by Victor Lin.
# Adapted from test_Mymodule.py by Jeff Chang.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Tests for motifs module."""

import tempfile
import unittest
import math

try:
    import numpy as np
except ImportError:
    from Bio import MissingExternalDependencyError

    raise MissingExternalDependencyError(
        "Install numpy if you want to use Bio.motifs."
    ) from None

from Bio import BiopythonDeprecationWarning
from Bio import motifs
from Bio.Seq import Seq


class TestBasic(unittest.TestCase):
    """Basic motif tests."""

    def test_format(self):
        m = motifs.create([Seq("ATATA")])
        m.name = "Foo"
        s1 = format(m, "pfm")
        expected_pfm = """  1.00   0.00   1.00   0.00  1.00
  0.00   0.00   0.00   0.00  0.00
  0.00   0.00   0.00   0.00  0.00
  0.00   1.00   0.00   1.00  0.00
"""
        s2 = format(m, "jaspar")
        expected_jaspar = """>None Foo
A [  1.00   0.00   1.00   0.00   1.00]
C [  0.00   0.00   0.00   0.00   0.00]
G [  0.00   0.00   0.00   0.00   0.00]
T [  0.00   1.00   0.00   1.00   0.00]
"""
        self.assertEqual(s2, expected_jaspar)
        s3 = format(m, "transfac")
        expected_transfac = """P0      A      C      G      T
01      1      0      0      0      A
02      0      0      0      1      T
03      1      0      0      0      A
04      0      0      0      1      T
05      1      0      0      0      A
XX
//
"""
        self.assertEqual(s3, expected_transfac)
        self.assertRaises(ValueError, format, m, "foo_bar")

    def test_relative_entropy(self):
        m = motifs.create([Seq("ATATA"), Seq("ATCTA"), Seq("TTGTA")])
        self.assertEqual(len(m.alignment), 3)
        self.assertEqual(m.background, {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25})
        self.assertEqual(m.pseudocounts, {"A": 0.0, "C": 0.0, "G": 0.0, "T": 0.0})
        self.assertTrue(
            np.allclose(
                m.relative_entropy,
                np.array([1.0817041659455104, 2.0, 0.4150374992788437, 2.0, 2.0]),
            )
        )
        m.background = {"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3}
        self.assertTrue(
            np.allclose(
                m.relative_entropy,
                np.array(
                    [
                        0.8186697601117167,
                        1.7369655941662063,
                        0.5419780939258206,
                        1.7369655941662063,
                        1.7369655941662063,
                    ]
                ),
            )
        )
        m.background = None
        self.assertEqual(m.background, {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25})
        pseudocounts = math.sqrt(len(m.alignment))
        m.pseudocounts = {
            letter: m.background[letter] * pseudocounts for letter in "ACGT"
        }
        self.assertTrue(
            np.allclose(
                m.relative_entropy,
                np.array(
                    [
                        0.3532586861097656,
                        0.7170228827697498,
                        0.11859369972847714,
                        0.7170228827697498,
                        0.7170228827697499,
                    ]
                ),
            )
        )
        m.background = {"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3}
        self.assertTrue(
            np.allclose(
                m.relative_entropy,
                np.array(
                    [
                        0.19727984803857979,
                        0.561044044698564,
                        0.20984910512125132,
                        0.561044044698564,
                        0.5610440446985638,
                    ]
                ),
            )
        )

    def test_reverse_complement(self):
        """Test if motifs can be reverse-complemented."""
        background = {"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3}
        pseudocounts = 0.5
        m = motifs.create([Seq("ATATA")])
        m.background = background
        m.pseudocounts = pseudocounts
        received_forward = format(m, "transfac")
        expected_forward = """\
P0      A      C      G      T
01      1      0      0      0      A
02      0      0      0      1      T
03      1      0      0      0      A
04      0      0      0      1      T
05      1      0      0      0      A
XX
//
"""
        self.assertEqual(received_forward, expected_forward)
        expected_forward_pwm = """\
        0      1      2      3      4
A:   0.50   0.17   0.50   0.17   0.50
C:   0.17   0.17   0.17   0.17   0.17
G:   0.17   0.17   0.17   0.17   0.17
T:   0.17   0.50   0.17   0.50   0.17
"""
        self.assertEqual(str(m.pwm), expected_forward_pwm)
        m = m.reverse_complement()
        received_reverse = format(m, "transfac")
        expected_reverse = """\
P0      A      C      G      T
01      0      0      0      1      T
02      1      0      0      0      A
03      0      0      0      1      T
04      1      0      0      0      A
05      0      0      0      1      T
XX
//
"""
        self.assertEqual(received_reverse, expected_reverse)
        expected_reverse_pwm = """\
        0      1      2      3      4
A:   0.17   0.50   0.17   0.50   0.17
C:   0.17   0.17   0.17   0.17   0.17
G:   0.17   0.17   0.17   0.17   0.17
T:   0.50   0.17   0.50   0.17   0.50
"""
        self.assertEqual(str(m.pwm), expected_reverse_pwm)
        # Same thing, but now start with a motif calculated from a count matrix
        m = motifs.create([Seq("ATATA")])
        counts = m.counts
        m = motifs.Motif(counts=counts)
        m.background = background
        m.pseudocounts = pseudocounts
        received_forward = format(m, "transfac")
        self.assertEqual(received_forward, expected_forward)
        self.assertEqual(str(m.pwm), expected_forward_pwm)
        m = m.reverse_complement()
        received_reverse = format(m, "transfac")
        self.assertEqual(received_reverse, expected_reverse)
        self.assertEqual(str(m.pwm), expected_reverse_pwm)


class TestAlignAce(unittest.TestCase):
    """Testing parsing AlignAce output files."""

    def test_alignace_parsing(self):
        """Test if Bio.motifs can parse AlignAce output files."""
        with open("motifs/alignace.out") as stream:
            record = motifs.parse(stream, "AlignAce")
        self.assertEqual(record.version, "AlignACE 4.0 05/13/04")
        self.assertEqual(record.command, "./AlignACE -i test.fa")
        self.assertEqual(len(record.parameters), 7)
        self.assertEqual(record.parameters["expect"], "10")
        self.assertEqual(record.parameters["gcback"], "0.38")
        self.assertEqual(record.parameters["minpass"], "200")
        self.assertEqual(record.parameters["seed"], "1227623309")
        self.assertEqual(record.parameters["numcols"], "10")
        self.assertEqual(record.parameters["undersample"], "1")
        self.assertEqual(record.parameters["oversample"], "1")
        self.assertEqual(len(record.sequences), 10)
        self.assertEqual(record.sequences[0], "SEQ1; M: CTCAATCGTAGA at 52")
        self.assertEqual(record.sequences[1], "SEQ2; M: CTCAATCGTAGA at 172")
        self.assertEqual(record.sequences[2], "SEQ3; M: CTCAATCGTAGA at 112")
        self.assertEqual(record.sequences[3], "SEQ4; M: CTCAATCGTAGA at 173")
        self.assertEqual(record.sequences[4], "SEQ5; M: CTCAATCGTAGA at 185")
        self.assertEqual(record.sequences[5], "SEQ6; M: CTCAATCGTAGA at 105")
        self.assertEqual(record.sequences[6], "SEQ7; M: CTCAATCGTAGA at 177")
        self.assertEqual(record.sequences[7], "SEQ8; M: CTCAATCGTAGA at 172")
        self.assertEqual(record.sequences[8], "SEQ9; M: CTCAATCGTAGA at 93")
        self.assertEqual(record.sequences[9], "SEQ10; M: CTCAATCGTAGA at 3")
        self.assertEqual(len(record), 16)
        self.assertEqual(record[0].alphabet, "ACGT")
        # using the old instances property:
        with self.assertWarns(BiopythonDeprecationWarning):
            self.assertEqual(len(record[0].instances), 11)
            self.assertEqual(record[0].instances[0], "TCTACGATTGAG")
            self.assertEqual(record[0].instances[1], "TCTACGATTGAG")
            self.assertEqual(record[0].instances[2], "TCTACGATTGAG")
            self.assertEqual(record[0].instances[3], "TCTACGATTGAG")
            self.assertEqual(record[0].instances[4], "TCTACGATTGAG")
            self.assertEqual(record[0].instances[5], "TCTACGATTGAG")
            self.assertEqual(record[0].instances[6], "TCTACGATTGAG")
            self.assertEqual(record[0].instances[7], "TCTACGATTGAG")
            self.assertEqual(record[0].instances[8], "TCTACGATTGAG")
            self.assertEqual(record[0].instances[9], "TCAAAGATAGAG")
            self.assertEqual(record[0].instances[10], "TCTACGATTGAG")
        self.assertEqual(len(record[0].alignment.sequences), 11)
        self.assertEqual(record[0].alignment.sequences[0], "TCTACGATTGAG")
        self.assertEqual(record[0].alignment.sequences[1], "TCTACGATTGAG")
        self.assertEqual(record[0].alignment.sequences[2], "TCTACGATTGAG")
        self.assertEqual(record[0].alignment.sequences[3], "TCTACGATTGAG")
        self.assertEqual(record[0].alignment.sequences[4], "TCTACGATTGAG")
        self.assertEqual(record[0].alignment.sequences[5], "TCTACGATTGAG")
        self.assertEqual(record[0].alignment.sequences[6], "TCTACGATTGAG")
        self.assertEqual(record[0].alignment.sequences[7], "TCTACGATTGAG")
        self.assertEqual(record[0].alignment.sequences[8], "TCTACGATTGAG")
        self.assertEqual(record[0].alignment.sequences[9], "TCAAAGATAGAG")
        self.assertEqual(record[0].alignment.sequences[10], "TCTACGATTGAG")
        self.assertEqual(record[0].mask, (1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1))
        self.assertAlmostEqual(record[0].score, 57.9079)
        self.assertEqual(
            str(record[0]),
            """\
TCTACGATTGAG
TCTACGATTGAG
TCTACGATTGAG
TCTACGATTGAG
TCTACGATTGAG
TCTACGATTGAG
TCTACGATTGAG
TCTACGATTGAG
TCTACGATTGAG
TCAAAGATAGAG
TCTACGATTGAG""",
        )
        motif = record[0][2:-1]
        self.assertEqual(motif.alphabet, "ACGT")
        self.assertEqual(
            str(motif),
            """\
TACGATTGA
TACGATTGA
TACGATTGA
TACGATTGA
TACGATTGA
TACGATTGA
TACGATTGA
TACGATTGA
TACGATTGA
AAAGATAGA
TACGATTGA""",
        )
        self.assertEqual(motif.mask, (0, 1, 1, 1, 1, 1, 0, 1, 1))
        self.assertEqual(record[1].alphabet, "ACGT")
        # using the old instances property:
        with self.assertWarns(BiopythonDeprecationWarning):
            self.assertEqual(len(record[1].instances), 22)
            self.assertEqual(record[1].instances[0], "GCGAAGGAAGCAGCGCGTGTG")
            self.assertEqual(record[1].instances[1], "GGCACCGCCTCTACGATTGAG")
            self.assertEqual(record[1].instances[2], "CAGAGCTTAGCATTGAACGCG")
            self.assertEqual(record[1].instances[3], "CTAATGAAAGCAATGAGAGTG")
            self.assertEqual(record[1].instances[4], "CTTGTGCCCTCTAAGCGTCCG")
            self.assertEqual(record[1].instances[5], "GAGCACGACGCTTTGTACCTG")
            self.assertEqual(record[1].instances[6], "CGGCACTTAGCAGCGTATCGT")
            self.assertEqual(record[1].instances[7], "CTGGTTTCATCTACGATTGAG")
            self.assertEqual(record[1].instances[8], "GGGCCAATAGCGGCGCCGGAG")
            self.assertEqual(record[1].instances[9], "GTGGAGTTATCTTAGTGCGCG")
            self.assertEqual(record[1].instances[10], "GAGAGGTTATCTACGATTGAG")
            self.assertEqual(record[1].instances[11], "CTGCTCCCCGCATACAGCGCG")
            self.assertEqual(record[1].instances[12], "CAGAACCGAGGTCCGGTACGG")
            self.assertEqual(record[1].instances[13], "GTGCCCCAAGCTTACCCAGGG")
            self.assertEqual(record[1].instances[14], "CGCCTCTGATCTACGATTGAG")
            self.assertEqual(record[1].instances[15], "GTGCTCATAGGGACGTCGCGG")
            self.assertEqual(record[1].instances[16], "CTGCCCCCCGCATAGTAGGGG")
            self.assertEqual(record[1].instances[17], "GTAAAGAAATCGATGTGCCAG")
            self.assertEqual(record[1].instances[18], "CACCTGCAATTGCTGGCAGCG")
            self.assertEqual(record[1].instances[19], "GGCGGGCCATCCCTGTATGAA")
            self.assertEqual(record[1].instances[20], "CTCCAGGTCGCATGGAGAGAG")
            self.assertEqual(record[1].instances[21], "CCTCGGATCGCTTGGGAAGAG")
        self.assertEqual(len(record[1].alignment.sequences), 22)
        self.assertEqual(record[1].alignment.sequences[0], "GCGAAGGAAGCAGCGCGTGTG")
        self.assertEqual(record[1].alignment.sequences[1], "GGCACCGCCTCTACGATTGAG")
        self.assertEqual(record[1].alignment.sequences[2], "CAGAGCTTAGCATTGAACGCG")
        self.assertEqual(record[1].alignment.sequences[3], "CTAATGAAAGCAATGAGAGTG")
        self.assertEqual(record[1].alignment.sequences[4], "CTTGTGCCCTCTAAGCGTCCG")
        self.assertEqual(record[1].alignment.sequences[5], "GAGCACGACGCTTTGTACCTG")
        self.assertEqual(record[1].alignment.sequences[6], "CGGCACTTAGCAGCGTATCGT")
        self.assertEqual(record[1].alignment.sequences[7], "CTGGTTTCATCTACGATTGAG")
        self.assertEqual(record[1].alignment.sequences[8], "GGGCCAATAGCGGCGCCGGAG")
        self.assertEqual(record[1].alignment.sequences[9], "GTGGAGTTATCTTAGTGCGCG")
        self.assertEqual(record[1].alignment.sequences[10], "GAGAGGTTATCTACGATTGAG")
        self.assertEqual(record[1].alignment.sequences[11], "CTGCTCCCCGCATACAGCGCG")
        self.assertEqual(record[1].alignment.sequences[12], "CAGAACCGAGGTCCGGTACGG")
        self.assertEqual(record[1].alignment.sequences[13], "GTGCCCCAAGCTTACCCAGGG")
        self.assertEqual(record[1].alignment.sequences[14], "CGCCTCTGATCTACGATTGAG")
        self.assertEqual(record[1].alignment.sequences[15], "GTGCTCATAGGGACGTCGCGG")
        self.assertEqual(record[1].alignment.sequences[16], "CTGCCCCCCGCATAGTAGGGG")
        self.assertEqual(record[1].alignment.sequences[17], "GTAAAGAAATCGATGTGCCAG")
        self.assertEqual(record[1].alignment.sequences[18], "CACCTGCAATTGCTGGCAGCG")
        self.assertEqual(record[1].alignment.sequences[19], "GGCGGGCCATCCCTGTATGAA")
        self.assertEqual(record[1].alignment.sequences[20], "CTCCAGGTCGCATGGAGAGAG")
        self.assertEqual(record[1].alignment.sequences[21], "CCTCGGATCGCTTGGGAAGAG")
        self.assertEqual(
            record[1].mask,
            (1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1),
        )
        self.assertAlmostEqual(record[1].score, 19.6235)
        self.assertEqual(
            str(record[1]),
            """\
GCGAAGGAAGCAGCGCGTGTG
GGCACCGCCTCTACGATTGAG
CAGAGCTTAGCATTGAACGCG
CTAATGAAAGCAATGAGAGTG
CTTGTGCCCTCTAAGCGTCCG
GAGCACGACGCTTTGTACCTG
CGGCACTTAGCAGCGTATCGT
CTGGTTTCATCTACGATTGAG
GGGCCAATAGCGGCGCCGGAG
GTGGAGTTATCTTAGTGCGCG
GAGAGGTTATCTACGATTGAG
CTGCTCCCCGCATACAGCGCG
CAGAACCGAGGTCCGGTACGG
GTGCCCCAAGCTTACCCAGGG
CGCCTCTGATCTACGATTGAG
GTGCTCATAGGGACGTCGCGG
CTGCCCCCCGCATAGTAGGGG
GTAAAGAAATCGATGTGCCAG
CACCTGCAATTGCTGGCAGCG
GGCGGGCCATCCCTGTATGAA
CTCCAGGTCGCATGGAGAGAG
CCTCGGATCGCTTGGGAAGAG""",
        )
        motif = record[1][2:-1]
        self.assertEqual(motif.alphabet, "ACGT")
        self.assertEqual(
            str(motif),
            """\
GAAGGAAGCAGCGCGTGT
CACCGCCTCTACGATTGA
GAGCTTAGCATTGAACGC
AATGAAAGCAATGAGAGT
TGTGCCCTCTAAGCGTCC
GCACGACGCTTTGTACCT
GCACTTAGCAGCGTATCG
GGTTTCATCTACGATTGA
GCCAATAGCGGCGCCGGA
GGAGTTATCTTAGTGCGC
GAGGTTATCTACGATTGA
GCTCCCCGCATACAGCGC
GAACCGAGGTCCGGTACG
GCCCCAAGCTTACCCAGG
CCTCTGATCTACGATTGA
GCTCATAGGGACGTCGCG
GCCCCCCGCATAGTAGGG
AAAGAAATCGATGTGCCA
CCTGCAATTGCTGGCAGC
CGGGCCATCCCTGTATGA
CCAGGTCGCATGGAGAGA
TCGGATCGCTTGGGAAGA""",
        )

        self.assertEqual(record[2].alphabet, "ACGT")
        # using the old instances property:
        with self.assertWarns(BiopythonDeprecationWarning):
            self.assertEqual(len(record[2].instances), 18)
            self.assertEqual(record[2].instances[0], "GTGCGCGAAGGAAGCAGCGCG")
            self.assertEqual(record[2].instances[1], "CAGAGCTTAGCATTGAACGCG")
            self.assertEqual(record[2].instances[2], "GTGCCCGATGACCACCCGTCG")
            self.assertEqual(record[2].instances[3], "GCCCTCTAAGCGTCCGCGGAT")
            self.assertEqual(record[2].instances[4], "GAGCACGACGCTTTGTACCTG")
            self.assertEqual(record[2].instances[5], "CGGCACTTAGCAGCGTATCGT")
            self.assertEqual(record[2].instances[6], "GGGCCAATAGCGGCGCCGGAG")
            self.assertEqual(record[2].instances[7], "GCGCACTAAGATAACTCCACG")
            self.assertEqual(record[2].instances[8], "CGGCCCGTTGTCCAGCAGACG")
            self.assertEqual(record[2].instances[9], "CTGCTCCCCGCATACAGCGCG")
            self.assertEqual(record[2].instances[10], "GTGCCCCAAGCTTACCCAGGG")
            self.assertEqual(record[2].instances[11], "GTGCTCATAGGGACGTCGCGG")
            self.assertEqual(record[2].instances[12], "CTGCCCCCCGCATAGTAGGGG")
            self.assertEqual(record[2].instances[13], "CGCCGCCATGCGACGCAGAGG")
            self.assertEqual(record[2].instances[14], "AACCTCTAAGCATACTCTACG")
            self.assertEqual(record[2].instances[15], "GACCTGGAGGCTTAGACTTGG")
            self.assertEqual(record[2].instances[16], "GCGCTCTTCCCAAGCGATCCG")
            self.assertEqual(record[2].instances[17], "GGGCCGTCAGCTCTCAAGTCT")
        self.assertEqual(len(record[2].alignment.sequences), 18)
        self.assertEqual(record[2].alignment.sequences[0], "GTGCGCGAAGGAAGCAGCGCG")
        self.assertEqual(record[2].alignment.sequences[1], "CAGAGCTTAGCATTGAACGCG")
        self.assertEqual(record[2].alignment.sequences[2], "GTGCCCGATGACCACCCGTCG")
        self.assertEqual(record[2].alignment.sequences[3], "GCCCTCTAAGCGTCCGCGGAT")
        self.assertEqual(record[2].alignment.sequences[4], "GAGCACGACGCTTTGTACCTG")
        self.assertEqual(record[2].alignment.sequences[5], "CGGCACTTAGCAGCGTATCGT")
        self.assertEqual(record[2].alignment.sequences[6], "GGGCCAATAGCGGCGCCGGAG")
        self.assertEqual(record[2].alignment.sequences[7], "GCGCACTAAGATAACTCCACG")
        self.assertEqual(record[2].alignment.sequences[8], "CGGCCCGTTGTCCAGCAGACG")
        self.assertEqual(record[2].alignment.sequences[9], "CTGCTCCCCGCATACAGCGCG")
        self.assertEqual(record[2].alignment.sequences[10], "GTGCCCCAAGCTTACCCAGGG")
        self.assertEqual(record[2].alignment.sequences[11], "GTGCTCATAGGGACGTCGCGG")
        self.assertEqual(record[2].alignment.sequences[12], "CTGCCCCCCGCATAGTAGGGG")
        self.assertEqual(record[2].alignment.sequences[13], "CGCCGCCATGCGACGCAGAGG")
        self.assertEqual(record[2].alignment.sequences[14], "AACCTCTAAGCATACTCTACG")
        self.assertEqual(record[2].alignment.sequences[15], "GACCTGGAGGCTTAGACTTGG")
        self.assertEqual(record[2].alignment.sequences[16], "GCGCTCTTCCCAAGCGATCCG")
        self.assertEqual(record[2].alignment.sequences[17], "GGGCCGTCAGCTCTCAAGTCT")
        self.assertEqual(
            record[2].mask,
            (1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1),
        )
        self.assertAlmostEqual(record[2].score, 19.1804)
        self.assertEqual(
            str(record[2]),
            """\
GTGCGCGAAGGAAGCAGCGCG
CAGAGCTTAGCATTGAACGCG
GTGCCCGATGACCACCCGTCG
GCCCTCTAAGCGTCCGCGGAT
GAGCACGACGCTTTGTACCTG
CGGCACTTAGCAGCGTATCGT
GGGCCAATAGCGGCGCCGGAG
GCGCACTAAGATAACTCCACG
CGGCCCGTTGTCCAGCAGACG
CTGCTCCCCGCATACAGCGCG
GTGCCCCAAGCTTACCCAGGG
GTGCTCATAGGGACGTCGCGG
CTGCCCCCCGCATAGTAGGGG
CGCCGCCATGCGACGCAGAGG
AACCTCTAAGCATACTCTACG
GACCTGGAGGCTTAGACTTGG
GCGCTCTTCCCAAGCGATCCG
GGGCCGTCAGCTCTCAAGTCT""",
        )
        motif = record[2][2:-1]
        self.assertEqual(motif.alphabet, "ACGT")
        self.assertEqual(
            str(motif),
            """\
GCGCGAAGGAAGCAGCGC
GAGCTTAGCATTGAACGC
GCCCGATGACCACCCGTC
CCTCTAAGCGTCCGCGGA
GCACGACGCTTTGTACCT
GCACTTAGCAGCGTATCG
GCCAATAGCGGCGCCGGA
GCACTAAGATAACTCCAC
GCCCGTTGTCCAGCAGAC
GCTCCCCGCATACAGCGC
GCCCCAAGCTTACCCAGG
GCTCATAGGGACGTCGCG
GCCCCCCGCATAGTAGGG
CCGCCATGCGACGCAGAG
CCTCTAAGCATACTCTAC
CCTGGAGGCTTAGACTTG
GCTCTTCCCAAGCGATCC
GCCGTCAGCTCTCAAGTC""",
        )

        self.assertEqual(record[3].alphabet, "ACGT")
        # using the old instances property:
        with self.assertWarns(BiopythonDeprecationWarning):
            self.assertEqual(len(record[3].instances), 16)
            self.assertEqual(record[3].instances[0], "GCCCCAAGCTTACCCAGGGAC")
            self.assertEqual(record[3].instances[1], "GCCGTCTGCTGGACAACGGGC")
            self.assertEqual(record[3].instances[2], "GCCGACGGGTGGTCATCGGGC")
            self.assertEqual(record[3].instances[3], "GCCAATAGCGGCGCCGGAGTC")
            self.assertEqual(record[3].instances[4], "GCCCCCCGCATAGTAGGGGGA")
            self.assertEqual(record[3].instances[5], "GCCCGTACCGGACCTCGGTTC")
            self.assertEqual(record[3].instances[6], "GCCTCATGTACCGGAAGGGAC")
            self.assertEqual(record[3].instances[7], "GACACGCGCCTGGGAGGGTTC")
            self.assertEqual(record[3].instances[8], "GCCTTTGGCCTTGGATGAGAA")
            self.assertEqual(record[3].instances[9], "GGCCCTCGGATCGCTTGGGAA")
            self.assertEqual(record[3].instances[10], "GCATGTTGGGAATCCGCGGAC")
            self.assertEqual(record[3].instances[11], "GACACGCGCTGTATGCGGGGA")
            self.assertEqual(record[3].instances[12], "GCCAGGTACAAAGCGTCGTGC")
            self.assertEqual(record[3].instances[13], "GCGATCAGCTTGTGGGCGTGC")
            self.assertEqual(record[3].instances[14], "GACAAATCGGATACTGGGGCA")
            self.assertEqual(record[3].instances[15], "GCACTTAGCAGCGTATCGTTA")
        self.assertEqual(len(record[3].alignment.sequences), 16)
        self.assertEqual(record[3].alignment.sequences[0], "GCCCCAAGCTTACCCAGGGAC")
        self.assertEqual(record[3].alignment.sequences[1], "GCCGTCTGCTGGACAACGGGC")
        self.assertEqual(record[3].alignment.sequences[2], "GCCGACGGGTGGTCATCGGGC")
        self.assertEqual(record[3].alignment.sequences[3], "GCCAATAGCGGCGCCGGAGTC")
        self.assertEqual(record[3].alignment.sequences[4], "GCCCCCCGCATAGTAGGGGGA")
        self.assertEqual(record[3].alignment.sequences[5], "GCCCGTACCGGACCTCGGTTC")
        self.assertEqual(record[3].alignment.sequences[6], "GCCTCATGTACCGGAAGGGAC")
        self.assertEqual(record[3].alignment.sequences[7], "GACACGCGCCTGGGAGGGTTC")
        self.assertEqual(record[3].alignment.sequences[8], "GCCTTTGGCCTTGGATGAGAA")
        self.assertEqual(record[3].alignment.sequences[9], "GGCCCTCGGATCGCTTGGGAA")
        self.assertEqual(record[3].alignment.sequences[10], "GCATGTTGGGAATCCGCGGAC")
        self.assertEqual(record[3].alignment.sequences[11], "GACACGCGCTGTATGCGGGGA")
        self.assertEqual(record[3].alignment.sequences[12], "GCCAGGTACAAAGCGTCGTGC")
        self.assertEqual(record[3].alignment.sequences[13], "GCGATCAGCTTGTGGGCGTGC")
        self.assertEqual(record[3].alignment.sequences[14], "GACAAATCGGATACTGGGGCA")
        self.assertEqual(record[3].alignment.sequences[15], "GCACTTAGCAGCGTATCGTTA")
        self.assertEqual(
            record[3].mask,
            (1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1),
        )
        self.assertAlmostEqual(record[3].score, 18.0097)
        self.assertEqual(
            str(record[3]),
            """\
GCCCCAAGCTTACCCAGGGAC
GCCGTCTGCTGGACAACGGGC
GCCGACGGGTGGTCATCGGGC
GCCAATAGCGGCGCCGGAGTC
GCCCCCCGCATAGTAGGGGGA
GCCCGTACCGGACCTCGGTTC
GCCTCATGTACCGGAAGGGAC
GACACGCGCCTGGGAGGGTTC
GCCTTTGGCCTTGGATGAGAA
GGCCCTCGGATCGCTTGGGAA
GCATGTTGGGAATCCGCGGAC
GACACGCGCTGTATGCGGGGA
GCCAGGTACAAAGCGTCGTGC
GCGATCAGCTTGTGGGCGTGC
GACAAATCGGATACTGGGGCA
GCACTTAGCAGCGTATCGTTA""",
        )
        motif = record[3][2:-1]
        self.assertEqual(motif.alphabet, "ACGT")
        self.assertEqual(
            str(motif),
            """\
CCCAAGCTTACCCAGGGA
CGTCTGCTGGACAACGGG
CGACGGGTGGTCATCGGG
CAATAGCGGCGCCGGAGT
CCCCCGCATAGTAGGGGG
CCGTACCGGACCTCGGTT
CTCATGTACCGGAAGGGA
CACGCGCCTGGGAGGGTT
CTTTGGCCTTGGATGAGA
CCCTCGGATCGCTTGGGA
ATGTTGGGAATCCGCGGA
CACGCGCTGTATGCGGGG
CAGGTACAAAGCGTCGTG
GATCAGCTTGTGGGCGTG
CAAATCGGATACTGGGGC
ACTTAGCAGCGTATCGTT""",
        )
        self.assertEqual(record[4].alphabet, "ACGT")
        # using the old instances property:
        with self.assertWarns(BiopythonDeprecationWarning):
            self.assertEqual(len(record[4].instances), 15)
            self.assertEqual(record[4].instances[0], "CGGCACAGAGCTT")
            self.assertEqual(record[4].instances[1], "ATCCGCGGACGCT")
            self.assertEqual(record[4].instances[2], "CGCCTGGGAGGGT")
            self.assertEqual(record[4].instances[3], "CGGAAGGGACGTT")
            self.assertEqual(record[4].instances[4], "ACACACAGACGGT")
            self.assertEqual(record[4].instances[5], "TGCCAGAGAGGTT")
            self.assertEqual(record[4].instances[6], "AGACTGAGACGTT")
            self.assertEqual(record[4].instances[7], "AATCGTAGAGGAT")
            self.assertEqual(record[4].instances[8], "CGTCTCGTAGGGT")
            self.assertEqual(record[4].instances[9], "CGTCGCGGAGGAT")
            self.assertEqual(record[4].instances[10], "CTTCTTAGACGCT")
            self.assertEqual(record[4].instances[11], "CGACGCAGAGGAT")
            self.assertEqual(record[4].instances[12], "ATGCTTAGAGGTT")
            self.assertEqual(record[4].instances[13], "AGACTTGGGCGAT")
            self.assertEqual(record[4].instances[14], "CGACCTGGAGGCT")
        self.assertEqual(len(record[4].alignment.sequences), 15)
        self.assertEqual(record[4].alignment.sequences[0], "CGGCACAGAGCTT")
        self.assertEqual(record[4].alignment.sequences[1], "ATCCGCGGACGCT")
        self.assertEqual(record[4].alignment.sequences[2], "CGCCTGGGAGGGT")
        self.assertEqual(record[4].alignment.sequences[3], "CGGAAGGGACGTT")
        self.assertEqual(record[4].alignment.sequences[4], "ACACACAGACGGT")
        self.assertEqual(record[4].alignment.sequences[5], "TGCCAGAGAGGTT")
        self.assertEqual(record[4].alignment.sequences[6], "AGACTGAGACGTT")
        self.assertEqual(record[4].alignment.sequences[7], "AATCGTAGAGGAT")
        self.assertEqual(record[4].alignment.sequences[8], "CGTCTCGTAGGGT")
        self.assertEqual(record[4].alignment.sequences[9], "CGTCGCGGAGGAT")
        self.assertEqual(record[4].alignment.sequences[10], "CTTCTTAGACGCT")
        self.assertEqual(record[4].alignment.sequences[11], "CGACGCAGAGGAT")
        self.assertEqual(record[4].alignment.sequences[12], "ATGCTTAGAGGTT")
        self.assertEqual(record[4].alignment.sequences[13], "AGACTTGGGCGAT")
        self.assertEqual(record[4].alignment.sequences[14], "CGACCTGGAGGCT")
        self.assertEqual(record[4].mask, (1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1))
        self.assertAlmostEqual(record[4].score, 16.8287)
        self.assertEqual(
            str(record[4]),
            """\
CGGCACAGAGCTT
ATCCGCGGACGCT
CGCCTGGGAGGGT
CGGAAGGGACGTT
ACACACAGACGGT
TGCCAGAGAGGTT
AGACTGAGACGTT
AATCGTAGAGGAT
CGTCTCGTAGGGT
CGTCGCGGAGGAT
CTTCTTAGACGCT
CGACGCAGAGGAT
ATGCTTAGAGGTT
AGACTTGGGCGAT
CGACCTGGAGGCT""",
        )
        motif = record[4][2:-1]
        self.assertEqual(motif.alphabet, "ACGT")
        self.assertEqual(
            str(motif),
            """\
GCACAGAGCT
CCGCGGACGC
CCTGGGAGGG
GAAGGGACGT
ACACAGACGG
CCAGAGAGGT
ACTGAGACGT
TCGTAGAGGA
TCTCGTAGGG
TCGCGGAGGA
TCTTAGACGC
ACGCAGAGGA
GCTTAGAGGT
ACTTGGGCGA
ACCTGGAGGC""",
        )
        self.assertEqual(record[5].alphabet, "ACGT")
        # using the old instances property:
        with self.assertWarns(BiopythonDeprecationWarning):
            self.assertEqual(len(record[5].instances), 18)
            self.assertEqual(record[5].instances[0], "GTGCGCGAAGGAAGCAGCGCGTG")
            self.assertEqual(record[5].instances[1], "TTGAGCCGAGTAAAGGGCTGGTG")
            self.assertEqual(record[5].instances[2], "CAATGCTAAGCTCTGTGCCGACG")
            self.assertEqual(record[5].instances[3], "CAACTCTCTATGTAGTGCCCGAG")
            self.assertEqual(record[5].instances[4], "CGACGCTTTGTACCTGGCTTGCG")
            self.assertEqual(record[5].instances[5], "CGAGTCAATGACACGCGCCTGGG")
            self.assertEqual(record[5].instances[6], "CGATACGCTGCTAAGTGCCGTCC")
            self.assertEqual(record[5].instances[7], "CCGGGCCAATAGCGGCGCCGGAG")
            self.assertEqual(record[5].instances[8], "CCACGCTTCGACACGTGGTATAG")
            self.assertEqual(record[5].instances[9], "CCGAGCCTCATGTACCGGAAGGG")
            self.assertEqual(record[5].instances[10], "CTGCTCCCCGCATACAGCGCGTG")
            self.assertEqual(record[5].instances[11], "CCGAGGTCCGGTACGGGCAAGCC")
            self.assertEqual(record[5].instances[12], "GTGCTCATAGGGACGTCGCGGAG")
            self.assertEqual(record[5].instances[13], "CCCTACTATGCGGGGGGCAGGTC")
            self.assertEqual(record[5].instances[14], "GCCAGCAATTGCAGGTGGTCGTG")
            self.assertEqual(record[5].instances[15], "CTCTGCGTCGCATGGCGGCGTGG")
            self.assertEqual(record[5].instances[16], "GGAGGCTTAGACTTGGGCGATAC")
            self.assertEqual(record[5].instances[17], "GCATGGAGAGAGATCCGGAGGAG")
        self.assertEqual(len(record[5].alignment.sequences), 18)
        self.assertEqual(record[5].alignment.sequences[0], "GTGCGCGAAGGAAGCAGCGCGTG")
        self.assertEqual(record[5].alignment.sequences[1], "TTGAGCCGAGTAAAGGGCTGGTG")
        self.assertEqual(record[5].alignment.sequences[2], "CAATGCTAAGCTCTGTGCCGACG")
        self.assertEqual(record[5].alignment.sequences[3], "CAACTCTCTATGTAGTGCCCGAG")
        self.assertEqual(record[5].alignment.sequences[4], "CGACGCTTTGTACCTGGCTTGCG")
        self.assertEqual(record[5].alignment.sequences[5], "CGAGTCAATGACACGCGCCTGGG")
        self.assertEqual(record[5].alignment.sequences[6], "CGATACGCTGCTAAGTGCCGTCC")
        self.assertEqual(record[5].alignment.sequences[7], "CCGGGCCAATAGCGGCGCCGGAG")
        self.assertEqual(record[5].alignment.sequences[8], "CCACGCTTCGACACGTGGTATAG")
        self.assertEqual(record[5].alignment.sequences[9], "CCGAGCCTCATGTACCGGAAGGG")
        self.assertEqual(record[5].alignment.sequences[10], "CTGCTCCCCGCATACAGCGCGTG")
        self.assertEqual(record[5].alignment.sequences[11], "CCGAGGTCCGGTACGGGCAAGCC")
        self.assertEqual(record[5].alignment.sequences[12], "GTGCTCATAGGGACGTCGCGGAG")
        self.assertEqual(record[5].alignment.sequences[13], "CCCTACTATGCGGGGGGCAGGTC")
        self.assertEqual(record[5].alignment.sequences[14], "GCCAGCAATTGCAGGTGGTCGTG")
        self.assertEqual(record[5].alignment.sequences[15], "CTCTGCGTCGCATGGCGGCGTGG")
        self.assertEqual(record[5].alignment.sequences[16], "GGAGGCTTAGACTTGGGCGATAC")
        self.assertEqual(record[5].alignment.sequences[17], "GCATGGAGAGAGATCCGGAGGAG")
        self.assertEqual(
            record[5].mask,
            (1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1),
        )
        self.assertAlmostEqual(record[5].score, 15.0441)
        self.assertEqual(
            str(record[5]),
            """\
GTGCGCGAAGGAAGCAGCGCGTG
TTGAGCCGAGTAAAGGGCTGGTG
CAATGCTAAGCTCTGTGCCGACG
CAACTCTCTATGTAGTGCCCGAG
CGACGCTTTGTACCTGGCTTGCG
CGAGTCAATGACACGCGCCTGGG
CGATACGCTGCTAAGTGCCGTCC
CCGGGCCAATAGCGGCGCCGGAG
CCACGCTTCGACACGTGGTATAG
CCGAGCCTCATGTACCGGAAGGG
CTGCTCCCCGCATACAGCGCGTG
CCGAGGTCCGGTACGGGCAAGCC
GTGCTCATAGGGACGTCGCGGAG
CCCTACTATGCGGGGGGCAGGTC
GCCAGCAATTGCAGGTGGTCGTG
CTCTGCGTCGCATGGCGGCGTGG
GGAGGCTTAGACTTGGGCGATAC
GCATGGAGAGAGATCCGGAGGAG""",
        )
        motif = record[5][2:-1]
        self.assertEqual(motif.alphabet, "ACGT")
        self.assertEqual(
            str(motif),
            """\
GCGCGAAGGAAGCAGCGCGT
GAGCCGAGTAAAGGGCTGGT
ATGCTAAGCTCTGTGCCGAC
ACTCTCTATGTAGTGCCCGA
ACGCTTTGTACCTGGCTTGC
AGTCAATGACACGCGCCTGG
ATACGCTGCTAAGTGCCGTC
GGGCCAATAGCGGCGCCGGA
ACGCTTCGACACGTGGTATA
GAGCCTCATGTACCGGAAGG
GCTCCCCGCATACAGCGCGT
GAGGTCCGGTACGGGCAAGC
GCTCATAGGGACGTCGCGGA
CTACTATGCGGGGGGCAGGT
CAGCAATTGCAGGTGGTCGT
CTGCGTCGCATGGCGGCGTG
AGGCTTAGACTTGGGCGATA
ATGGAGAGAGATCCGGAGGA""",
        )
        self.assertEqual(record[6].alphabet, "ACGT")
        # using the old instances property:
        with self.assertWarns(BiopythonDeprecationWarning):
            self.assertEqual(len(record[6].instances), 20)
            self.assertEqual(record[6].instances[0], "GCGCGTGTGTGTAAC")
            self.assertEqual(record[6].instances[1], "GCACAGAGCTTAGCA")
            self.assertEqual(record[6].instances[2], "GGTGGTCATCGGGCA")
            self.assertEqual(record[6].instances[3], "GCGCGTGTCATTGAC")
            self.assertEqual(record[6].instances[4], "GGACGGCACTTAGCA")
            self.assertEqual(record[6].instances[5], "GCGCGTCCCGGGCCA")
            self.assertEqual(record[6].instances[6], "GCTCGGCCCGTTGTC")
            self.assertEqual(record[6].instances[7], "GCGCGTGTCCTTTAA")
            self.assertEqual(record[6].instances[8], "GCTGATCGCTGCTCC")
            self.assertEqual(record[6].instances[9], "GCCCGTACCGGACCT")
            self.assertEqual(record[6].instances[10], "GGACGTCGCGGAGGA")
            self.assertEqual(record[6].instances[11], "GCGGGGGGCAGGTCA")
            self.assertEqual(record[6].instances[12], "GGACGTACTGGCACA")
            self.assertEqual(record[6].instances[13], "GCAGGTGGTCGTGCA")
            self.assertEqual(record[6].instances[14], "GCGCATACCTTAACA")
            self.assertEqual(record[6].instances[15], "GCACGGGACTTCAAC")
            self.assertEqual(record[6].instances[16], "GCACGTAGCTGGTAA")
            self.assertEqual(record[6].instances[17], "GCTCGTCTATGGTCA")
            self.assertEqual(record[6].instances[18], "GCGCATGCTGGATCC")
            self.assertEqual(record[6].instances[19], "GGCCGTCAGCTCTCA")
        self.assertEqual(len(record[6].alignment.sequences), 20)
        self.assertEqual(record[6].alignment.sequences[0], "GCGCGTGTGTGTAAC")
        self.assertEqual(record[6].alignment.sequences[1], "GCACAGAGCTTAGCA")
        self.assertEqual(record[6].alignment.sequences[2], "GGTGGTCATCGGGCA")
        self.assertEqual(record[6].alignment.sequences[3], "GCGCGTGTCATTGAC")
        self.assertEqual(record[6].alignment.sequences[4], "GGACGGCACTTAGCA")
        self.assertEqual(record[6].alignment.sequences[5], "GCGCGTCCCGGGCCA")
        self.assertEqual(record[6].alignment.sequences[6], "GCTCGGCCCGTTGTC")
        self.assertEqual(record[6].alignment.sequences[7], "GCGCGTGTCCTTTAA")
        self.assertEqual(record[6].alignment.sequences[8], "GCTGATCGCTGCTCC")
        self.assertEqual(record[6].alignment.sequences[9], "GCCCGTACCGGACCT")
        self.assertEqual(record[6].alignment.sequences[10], "GGACGTCGCGGAGGA")
        self.assertEqual(record[6].alignment.sequences[11], "GCGGGGGGCAGGTCA")
        self.assertEqual(record[6].alignment.sequences[12], "GGACGTACTGGCACA")
        self.assertEqual(record[6].alignment.sequences[13], "GCAGGTGGTCGTGCA")
        self.assertEqual(record[6].alignment.sequences[14], "GCGCATACCTTAACA")
        self.assertEqual(record[6].alignment.sequences[15], "GCACGGGACTTCAAC")
        self.assertEqual(record[6].alignment.sequences[16], "GCACGTAGCTGGTAA")
        self.assertEqual(record[6].alignment.sequences[17], "GCTCGTCTATGGTCA")
        self.assertEqual(record[6].alignment.sequences[18], "GCGCATGCTGGATCC")
        self.assertEqual(record[6].alignment.sequences[19], "GGCCGTCAGCTCTCA")
        self.assertEqual(record[6].mask, (1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1))
        self.assertAlmostEqual(record[6].score, 13.3145)
        self.assertEqual(
            str(record[6]),
            """\
GCGCGTGTGTGTAAC
GCACAGAGCTTAGCA
GGTGGTCATCGGGCA
GCGCGTGTCATTGAC
GGACGGCACTTAGCA
GCGCGTCCCGGGCCA
GCTCGGCCCGTTGTC
GCGCGTGTCCTTTAA
GCTGATCGCTGCTCC
GCCCGTACCGGACCT
GGACGTCGCGGAGGA
GCGGGGGGCAGGTCA
GGACGTACTGGCACA
GCAGGTGGTCGTGCA
GCGCATACCTTAACA
GCACGGGACTTCAAC
GCACGTAGCTGGTAA
GCTCGTCTATGGTCA
GCGCATGCTGGATCC
GGCCGTCAGCTCTCA""",
        )
        motif = record[6][2:-1]
        self.assertEqual(motif.alphabet, "ACGT")
        self.assertEqual(
            str(motif),
            """\
GCGTGTGTGTAA
ACAGAGCTTAGC
TGGTCATCGGGC
GCGTGTCATTGA
ACGGCACTTAGC
GCGTCCCGGGCC
TCGGCCCGTTGT
GCGTGTCCTTTA
TGATCGCTGCTC
CCGTACCGGACC
ACGTCGCGGAGG
GGGGGGCAGGTC
ACGTACTGGCAC
AGGTGGTCGTGC
GCATACCTTAAC
ACGGGACTTCAA
ACGTAGCTGGTA
TCGTCTATGGTC
GCATGCTGGATC
CCGTCAGCTCTC""",
        )
        self.assertEqual(record[7].alphabet, "ACGT")
        # using the old instances property:
        with self.assertWarns(BiopythonDeprecationWarning):
            self.assertEqual(len(record[7].instances), 20)
            self.assertEqual(record[7].instances[0], "GAACCGAGGTCCGGTACGGGC")
            self.assertEqual(record[7].instances[1], "GCCCCCCGCATAGTAGGGGGA")
            self.assertEqual(record[7].instances[2], "GTCCCTGGGTAAGCTTGGGGC")
            self.assertEqual(record[7].instances[3], "ACTCCACGCTTCGACACGTGG")
            self.assertEqual(record[7].instances[4], "ATCCTCTGCGTCGCATGGCGG")
            self.assertEqual(record[7].instances[5], "GTTCAATGCTAAGCTCTGTGC")
            self.assertEqual(record[7].instances[6], "GCTCATAGGGACGTCGCGGAG")
            self.assertEqual(record[7].instances[7], "GTCCCGGGCCAATAGCGGCGC")
            self.assertEqual(record[7].instances[8], "GCACTTAGCAGCGTATCGTTA")
            self.assertEqual(record[7].instances[9], "GGCCCTCGGATCGCTTGGGAA")
            self.assertEqual(record[7].instances[10], "CTGCTGGACAACGGGCCGAGC")
            self.assertEqual(record[7].instances[11], "GGGCACTACATAGAGAGTTGC")
            self.assertEqual(record[7].instances[12], "AGCCTCCAGGTCGCATGGAGA")
            self.assertEqual(record[7].instances[13], "AATCGTAGATCAGAGGCGAGA")
            self.assertEqual(record[7].instances[14], "GAACTCCACTAAGACTTGAGA")
            self.assertEqual(record[7].instances[15], "GAGCAGCGATCAGCTTGTGGG")
            self.assertEqual(record[7].instances[16], "GCCAGGTACAAAGCGTCGTGC")
            self.assertEqual(record[7].instances[17], "AGTCAATGACACGCGCCTGGG")
            self.assertEqual(record[7].instances[18], "GGTCATGGAATCTTATGTAGC")
            self.assertEqual(record[7].instances[19], "GTAGATAACAGAGGTCGGGGG")
        self.assertEqual(len(record[7].alignment.sequences), 20)
        self.assertEqual(record[7].alignment.sequences[0], "GAACCGAGGTCCGGTACGGGC")
        self.assertEqual(record[7].alignment.sequences[1], "GCCCCCCGCATAGTAGGGGGA")
        self.assertEqual(record[7].alignment.sequences[2], "GTCCCTGGGTAAGCTTGGGGC")
        self.assertEqual(record[7].alignment.sequences[3], "ACTCCACGCTTCGACACGTGG")
        self.assertEqual(record[7].alignment.sequences[4], "ATCCTCTGCGTCGCATGGCGG")
        self.assertEqual(record[7].alignment.sequences[5], "GTTCAATGCTAAGCTCTGTGC")
        self.assertEqual(record[7].alignment.sequences[6], "GCTCATAGGGACGTCGCGGAG")
        self.assertEqual(record[7].alignment.sequences[7], "GTCCCGGGCCAATAGCGGCGC")
        self.assertEqual(record[7].alignment.sequences[8], "GCACTTAGCAGCGTATCGTTA")
        self.assertEqual(record[7].alignment.sequences[9], "GGCCCTCGGATCGCTTGGGAA")
        self.assertEqual(record[7].alignment.sequences[10], "CTGCTGGACAACGGGCCGAGC")
        self.assertEqual(record[7].alignment.sequences[11], "GGGCACTACATAGAGAGTTGC")
        self.assertEqual(record[7].alignment.sequences[12], "AGCCTCCAGGTCGCATGGAGA")
        self.assertEqual(record[7].alignment.sequences[13], "AATCGTAGATCAGAGGCGAGA")
        self.assertEqual(record[7].alignment.sequences[14], "GAACTCCACTAAGACTTGAGA")
        self.assertEqual(record[7].alignment.sequences[15], "GAGCAGCGATCAGCTTGTGGG")
        self.assertEqual(record[7].alignment.sequences[16], "GCCAGGTACAAAGCGTCGTGC")
        self.assertEqual(record[7].alignment.sequences[17], "AGTCAATGACACGCGCCTGGG")
        self.assertEqual(record[7].alignment.sequences[18], "GGTCATGGAATCTTATGTAGC")
        self.assertEqual(record[7].alignment.sequences[19], "GTAGATAACAGAGGTCGGGGG")
        self.assertEqual(
            record[7].mask,
            (1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1),
        )
        self.assertAlmostEqual(record[7].score, 11.6098)
        self.assertEqual(
            str(record[7]),
            """\
GAACCGAGGTCCGGTACGGGC
GCCCCCCGCATAGTAGGGGGA
GTCCCTGGGTAAGCTTGGGGC
ACTCCACGCTTCGACACGTGG
ATCCTCTGCGTCGCATGGCGG
GTTCAATGCTAAGCTCTGTGC
GCTCATAGGGACGTCGCGGAG
GTCCCGGGCCAATAGCGGCGC
GCACTTAGCAGCGTATCGTTA
GGCCCTCGGATCGCTTGGGAA
CTGCTGGACAACGGGCCGAGC
GGGCACTACATAGAGAGTTGC
AGCCTCCAGGTCGCATGGAGA
AATCGTAGATCAGAGGCGAGA
GAACTCCACTAAGACTTGAGA
GAGCAGCGATCAGCTTGTGGG
GCCAGGTACAAAGCGTCGTGC
AGTCAATGACACGCGCCTGGG
GGTCATGGAATCTTATGTAGC
GTAGATAACAGAGGTCGGGGG""",
        )
        motif = record[7][2:-1]
        self.assertEqual(motif.alphabet, "ACGT")
        self.assertEqual(
            str(motif),
            """\
ACCGAGGTCCGGTACGGG
CCCCCGCATAGTAGGGGG
CCCTGGGTAAGCTTGGGG
TCCACGCTTCGACACGTG
CCTCTGCGTCGCATGGCG
TCAATGCTAAGCTCTGTG
TCATAGGGACGTCGCGGA
CCCGGGCCAATAGCGGCG
ACTTAGCAGCGTATCGTT
CCCTCGGATCGCTTGGGA
GCTGGACAACGGGCCGAG
GCACTACATAGAGAGTTG
CCTCCAGGTCGCATGGAG
TCGTAGATCAGAGGCGAG
ACTCCACTAAGACTTGAG
GCAGCGATCAGCTTGTGG
CAGGTACAAAGCGTCGTG
TCAATGACACGCGCCTGG
TCATGGAATCTTATGTAG
AGATAACAGAGGTCGGGG""",
        )
        self.assertEqual(record[8].alphabet, "ACGT")
        # using the old instances property:
        with self.assertWarns(BiopythonDeprecationWarning):
            self.assertEqual(len(record[8].instances), 14)
            self.assertEqual(record[8].instances[0], "CCGAGTAAAGGGCTG")
            self.assertEqual(record[8].instances[1], "GTGGTCATCGGGCAC")
            self.assertEqual(record[8].instances[2], "GATAACAGAGGTCGG")
            self.assertEqual(record[8].instances[3], "CGGCGCCGGAGTCTG")
            self.assertEqual(record[8].instances[4], "GCGCGTCCCGGGCCA")
            self.assertEqual(record[8].instances[5], "CTGGACAACGGGCCG")
            self.assertEqual(record[8].instances[6], "CGGATACTGGGGCAG")
            self.assertEqual(record[8].instances[7], "GGGAGCAGCGATCAG")
            self.assertEqual(record[8].instances[8], "CAGAACCGAGGTCCG")
            self.assertEqual(record[8].instances[9], "GGGTCCCTGGGTAAG")
            self.assertEqual(record[8].instances[10], "GTGCTCATAGGGACG")
            self.assertEqual(record[8].instances[11], "GAGATCCGGAGGAGG")
            self.assertEqual(record[8].instances[12], "GCGATCCGAGGGCCG")
            self.assertEqual(record[8].instances[13], "GAGTTCACATGGCTG")
        self.assertEqual(len(record[8].alignment.sequences), 14)
        self.assertEqual(record[8].alignment.sequences[0], "CCGAGTAAAGGGCTG")
        self.assertEqual(record[8].alignment.sequences[1], "GTGGTCATCGGGCAC")
        self.assertEqual(record[8].alignment.sequences[2], "GATAACAGAGGTCGG")
        self.assertEqual(record[8].alignment.sequences[3], "CGGCGCCGGAGTCTG")
        self.assertEqual(record[8].alignment.sequences[4], "GCGCGTCCCGGGCCA")
        self.assertEqual(record[8].alignment.sequences[5], "CTGGACAACGGGCCG")
        self.assertEqual(record[8].alignment.sequences[6], "CGGATACTGGGGCAG")
        self.assertEqual(record[8].alignment.sequences[7], "GGGAGCAGCGATCAG")
        self.assertEqual(record[8].alignment.sequences[8], "CAGAACCGAGGTCCG")
        self.assertEqual(record[8].alignment.sequences[9], "GGGTCCCTGGGTAAG")
        self.assertEqual(record[8].alignment.sequences[10], "GTGCTCATAGGGACG")
        self.assertEqual(record[8].alignment.sequences[11], "GAGATCCGGAGGAGG")
        self.assertEqual(record[8].alignment.sequences[12], "GCGATCCGAGGGCCG")
        self.assertEqual(record[8].alignment.sequences[13], "GAGTTCACATGGCTG")
        self.assertEqual(record[8].mask, (1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1))
        self.assertAlmostEqual(record[8].score, 11.2943)
        self.assertEqual(
            str(record[8]),
            """\
CCGAGTAAAGGGCTG
GTGGTCATCGGGCAC
GATAACAGAGGTCGG
CGGCGCCGGAGTCTG
GCGCGTCCCGGGCCA
CTGGACAACGGGCCG
CGGATACTGGGGCAG
GGGAGCAGCGATCAG
CAGAACCGAGGTCCG
GGGTCCCTGGGTAAG
GTGCTCATAGGGACG
GAGATCCGGAGGAGG
GCGATCCGAGGGCCG
GAGTTCACATGGCTG""",
        )
        motif = record[8][2:-1]
        self.assertEqual(motif.alphabet, "ACGT")
        self.assertEqual(
            str(motif),
            """\
GAGTAAAGGGCT
GGTCATCGGGCA
TAACAGAGGTCG
GCGCCGGAGTCT
GCGTCCCGGGCC
GGACAACGGGCC
GATACTGGGGCA
GAGCAGCGATCA
GAACCGAGGTCC
GTCCCTGGGTAA
GCTCATAGGGAC
GATCCGGAGGAG
GATCCGAGGGCC
GTTCACATGGCT""",
        )
        self.assertEqual(record[9].alphabet, "ACGT")
        # using the old instances property:
        with self.assertWarns(BiopythonDeprecationWarning):
            self.assertEqual(len(record[9].instances), 18)
            self.assertEqual(record[9].instances[0], "TAGAGGCGGTG")
            self.assertEqual(record[9].instances[1], "GCTAAGCTCTG")
            self.assertEqual(record[9].instances[2], "TGGAAGCAGTG")
            self.assertEqual(record[9].instances[3], "GCGAGGCTGTG")
            self.assertEqual(record[9].instances[4], "ACGACGCTTTG")
            self.assertEqual(record[9].instances[5], "GGGACGCGCAC")
            self.assertEqual(record[9].instances[6], "TCGAAGCGTGG")
            self.assertEqual(record[9].instances[7], "TGTATGCGGGG")
            self.assertEqual(record[9].instances[8], "GGTAAGCTTGG")
            self.assertEqual(record[9].instances[9], "TGTACGCTGGG")
            self.assertEqual(record[9].instances[10], "ACTATGCGGGG")
            self.assertEqual(record[9].instances[11], "GGTATGCGCTG")
            self.assertEqual(record[9].instances[12], "GGTACCCGGAG")
            self.assertEqual(record[9].instances[13], "GCGACGCAGAG")
            self.assertEqual(record[9].instances[14], "TGGCGGCGTGG")
            self.assertEqual(record[9].instances[15], "TCTAGGCGGGC")
            self.assertEqual(record[9].instances[16], "AGTATGCTTAG")
            self.assertEqual(record[9].instances[17], "TGGAGGCTTAG")
        self.assertEqual(len(record[9].alignment.sequences), 18)
        self.assertEqual(record[9].alignment.sequences[0], "TAGAGGCGGTG")
        self.assertEqual(record[9].alignment.sequences[1], "GCTAAGCTCTG")
        self.assertEqual(record[9].alignment.sequences[2], "TGGAAGCAGTG")
        self.assertEqual(record[9].alignment.sequences[3], "GCGAGGCTGTG")
        self.assertEqual(record[9].alignment.sequences[4], "ACGACGCTTTG")
        self.assertEqual(record[9].alignment.sequences[5], "GGGACGCGCAC")
        self.assertEqual(record[9].alignment.sequences[6], "TCGAAGCGTGG")
        self.assertEqual(record[9].alignment.sequences[7], "TGTATGCGGGG")
        self.assertEqual(record[9].alignment.sequences[8], "GGTAAGCTTGG")
        self.assertEqual(record[9].alignment.sequences[9], "TGTACGCTGGG")
        self.assertEqual(record[9].alignment.sequences[10], "ACTATGCGGGG")
        self.assertEqual(record[9].alignment.sequences[11], "GGTATGCGCTG")
        self.assertEqual(record[9].alignment.sequences[12], "GGTACCCGGAG")
        self.assertEqual(record[9].alignment.sequences[13], "GCGACGCAGAG")
        self.assertEqual(record[9].alignment.sequences[14], "TGGCGGCGTGG")
        self.assertEqual(record[9].alignment.sequences[15], "TCTAGGCGGGC")
        self.assertEqual(record[9].alignment.sequences[16], "AGTATGCTTAG")
        self.assertEqual(record[9].alignment.sequences[17], "TGGAGGCTTAG")
        self.assertEqual(record[9].mask, (1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1))
        self.assertAlmostEqual(record[9].score, 9.7924)
        self.assertEqual(
            str(record[9]),
            """\
TAGAGGCGGTG
GCTAAGCTCTG
TGGAAGCAGTG
GCGAGGCTGTG
ACGACGCTTTG
GGGACGCGCAC
TCGAAGCGTGG
TGTATGCGGGG
GGTAAGCTTGG
TGTACGCTGGG
ACTATGCGGGG
GGTATGCGCTG
GGTACCCGGAG
GCGACGCAGAG
TGGCGGCGTGG
TCTAGGCGGGC
AGTATGCTTAG
TGGAGGCTTAG""",
        )
        motif = record[9][2:-1]
        self.assertEqual(motif.alphabet, "ACGT")
        self.assertEqual(
            str(motif),
            """\
GAGGCGGT
TAAGCTCT
GAAGCAGT
GAGGCTGT
GACGCTTT
GACGCGCA
GAAGCGTG
TATGCGGG
TAAGCTTG
TACGCTGG
TATGCGGG
TATGCGCT
TACCCGGA
GACGCAGA
GCGGCGTG
TAGGCGGG
TATGCTTA
GAGGCTTA""",
        )
        self.assertEqual(record[10].alphabet, "ACGT")
        # using the old instances property:
        with self.assertWarns(BiopythonDeprecationWarning):
            self.assertEqual(len(record[10].instances), 13)
            self.assertEqual(record[10].instances[0], "GCACAGAGCTTAGCATTGAAC")
            self.assertEqual(record[10].instances[1], "GTCCGCGGATTCCCAACATGC")
            self.assertEqual(record[10].instances[2], "ATACACAGCCTCGCAAGCCAG")
            self.assertEqual(record[10].instances[3], "GGCCCGGGACGCGCACTAAGA")
            self.assertEqual(record[10].instances[4], "GCCCGTTGTCCAGCAGACGGC")
            self.assertEqual(record[10].instances[5], "GAGCAGCGATCAGCTTGTGGG")
            self.assertEqual(record[10].instances[6], "GAACCGAGGTCCGGTACGGGC")
            self.assertEqual(record[10].instances[7], "GTCCCTGGGTAAGCTTGGGGC")
            self.assertEqual(record[10].instances[8], "GACCTGCCCCCCGCATAGTAG")
            self.assertEqual(record[10].instances[9], "AACCAGCGCATACCTTAACAG")
            self.assertEqual(record[10].instances[10], "ATCCTCTGCGTCGCATGGCGG")
            self.assertEqual(record[10].instances[11], "GACCATAGACGAGCATCAAAG")
            self.assertEqual(record[10].instances[12], "GGCCCTCGGATCGCTTGGGAA")
        self.assertEqual(len(record[10].alignment.sequences), 13)
        self.assertEqual(record[10].alignment.sequences[0], "GCACAGAGCTTAGCATTGAAC")
        self.assertEqual(record[10].alignment.sequences[1], "GTCCGCGGATTCCCAACATGC")
        self.assertEqual(record[10].alignment.sequences[2], "ATACACAGCCTCGCAAGCCAG")
        self.assertEqual(record[10].alignment.sequences[3], "GGCCCGGGACGCGCACTAAGA")
        self.assertEqual(record[10].alignment.sequences[4], "GCCCGTTGTCCAGCAGACGGC")
        self.assertEqual(record[10].alignment.sequences[5], "GAGCAGCGATCAGCTTGTGGG")
        self.assertEqual(record[10].alignment.sequences[6], "GAACCGAGGTCCGGTACGGGC")
        self.assertEqual(record[10].alignment.sequences[7], "GTCCCTGGGTAAGCTTGGGGC")
        self.assertEqual(record[10].alignment.sequences[8], "GACCTGCCCCCCGCATAGTAG")
        self.assertEqual(record[10].alignment.sequences[9], "AACCAGCGCATACCTTAACAG")
        self.assertEqual(record[10].alignment.sequences[10], "ATCCTCTGCGTCGCATGGCGG")
        self.assertEqual(record[10].alignment.sequences[11], "GACCATAGACGAGCATCAAAG")
        self.assertEqual(record[10].alignment.sequences[12], "GGCCCTCGGATCGCTTGGGAA")
        self.assertEqual(
            record[10].mask,
            (1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1),
        )
        self.assertAlmostEqual(record[10].score, 9.01393)
        self.assertEqual(
            str(record[10]),
            """\
GCACAGAGCTTAGCATTGAAC
GTCCGCGGATTCCCAACATGC
ATACACAGCCTCGCAAGCCAG
GGCCCGGGACGCGCACTAAGA
GCCCGTTGTCCAGCAGACGGC
GAGCAGCGATCAGCTTGTGGG
GAACCGAGGTCCGGTACGGGC
GTCCCTGGGTAAGCTTGGGGC
GACCTGCCCCCCGCATAGTAG
AACCAGCGCATACCTTAACAG
ATCCTCTGCGTCGCATGGCGG
GACCATAGACGAGCATCAAAG
GGCCCTCGGATCGCTTGGGAA""",
        )
        motif = record[10][2:-1]
        self.assertEqual(motif.alphabet, "ACGT")
        self.assertEqual(
            str(motif),
            """\
ACAGAGCTTAGCATTGAA
CCGCGGATTCCCAACATG
ACACAGCCTCGCAAGCCA
CCCGGGACGCGCACTAAG
CCGTTGTCCAGCAGACGG
GCAGCGATCAGCTTGTGG
ACCGAGGTCCGGTACGGG
CCCTGGGTAAGCTTGGGG
CCTGCCCCCCGCATAGTA
CCAGCGCATACCTTAACA
CCTCTGCGTCGCATGGCG
CCATAGACGAGCATCAAA
CCCTCGGATCGCTTGGGA""",
        )
        self.assertEqual(record[11].alphabet, "ACGT")
        # using the old instances property:
        with self.assertWarns(BiopythonDeprecationWarning):
            self.assertEqual(len(record[11].instances), 16)
            self.assertEqual(record[11].instances[0], "GCCGTCCGTC")
            self.assertEqual(record[11].instances[1], "GGCGTGCGCG")
            self.assertEqual(record[11].instances[2], "GGCGCGTGTC")
            self.assertEqual(record[11].instances[3], "AGCGCGTGTG")
            self.assertEqual(record[11].instances[4], "GCGGTGCGTG")
            self.assertEqual(record[11].instances[5], "AGCGCGTGTC")
            self.assertEqual(record[11].instances[6], "AGCGTCCGCG")
            self.assertEqual(record[11].instances[7], "ACCGTCTGTG")
            self.assertEqual(record[11].instances[8], "GCCATGCGAC")
            self.assertEqual(record[11].instances[9], "ACCACCCGTC")
            self.assertEqual(record[11].instances[10], "GGCGCCGGAG")
            self.assertEqual(record[11].instances[11], "ACCACGTGTC")
            self.assertEqual(record[11].instances[12], "GGCTTGCGAG")
            self.assertEqual(record[11].instances[13], "GCGATCCGAG")
            self.assertEqual(record[11].instances[14], "AGTGCGCGTC")
            self.assertEqual(record[11].instances[15], "AGTGCCCGAG")
        self.assertEqual(len(record[11].alignment.sequences), 16)
        self.assertEqual(record[11].alignment.sequences[0], "GCCGTCCGTC")
        self.assertEqual(record[11].alignment.sequences[1], "GGCGTGCGCG")
        self.assertEqual(record[11].alignment.sequences[2], "GGCGCGTGTC")
        self.assertEqual(record[11].alignment.sequences[3], "AGCGCGTGTG")
        self.assertEqual(record[11].alignment.sequences[4], "GCGGTGCGTG")
        self.assertEqual(record[11].alignment.sequences[5], "AGCGCGTGTC")
        self.assertEqual(record[11].alignment.sequences[6], "AGCGTCCGCG")
        self.assertEqual(record[11].alignment.sequences[7], "ACCGTCTGTG")
        self.assertEqual(record[11].alignment.sequences[8], "GCCATGCGAC")
        self.assertEqual(record[11].alignment.sequences[9], "ACCACCCGTC")
        self.assertEqual(record[11].alignment.sequences[10], "GGCGCCGGAG")
        self.assertEqual(record[11].alignment.sequences[11], "ACCACGTGTC")
        self.assertEqual(record[11].alignment.sequences[12], "GGCTTGCGAG")
        self.assertEqual(record[11].alignment.sequences[13], "GCGATCCGAG")
        self.assertEqual(record[11].alignment.sequences[14], "AGTGCGCGTC")
        self.assertEqual(record[11].alignment.sequences[15], "AGTGCCCGAG")
        self.assertEqual(record[11].mask, (1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
        self.assertAlmostEqual(record[11].score, 7.51121)
        self.assertEqual(
            str(record[11]),
            """\
GCCGTCCGTC
GGCGTGCGCG
GGCGCGTGTC
AGCGCGTGTG
GCGGTGCGTG
AGCGCGTGTC
AGCGTCCGCG
ACCGTCTGTG
GCCATGCGAC
ACCACCCGTC
GGCGCCGGAG
ACCACGTGTC
GGCTTGCGAG
GCGATCCGAG
AGTGCGCGTC
AGTGCCCGAG""",
        )
        motif = record[11][2:-1]
        self.assertEqual(motif.alphabet, "ACGT")
        self.assertEqual(
            str(motif),
            """\
CGTCCGT
CGTGCGC
CGCGTGT
CGCGTGT
GGTGCGT
CGCGTGT
CGTCCGC
CGTCTGT
CATGCGA
CACCCGT
CGCCGGA
CACGTGT
CTTGCGA
GATCCGA
TGCGCGT
TGCCCGA""",
        )
        self.assertEqual(record[12].alphabet, "ACGT")
        # using the old instances property:
        with self.assertWarns(BiopythonDeprecationWarning):
            self.assertEqual(len(record[12].instances), 16)
            self.assertEqual(record[12].instances[0], "GCCGACGGGTGGTCATCGGG")
            self.assertEqual(record[12].instances[1], "GCACGACGCTTTGTACCTGG")
            self.assertEqual(record[12].instances[2], "CCTGGGAGGGTTCAATAACG")
            self.assertEqual(record[12].instances[3], "GCGCGTCCCGGGCCAATAGC")
            self.assertEqual(record[12].instances[4], "GCCGTCTGCTGGACAACGGG")
            self.assertEqual(record[12].instances[5], "GTCCCTTCCGGTACATGAGG")
            self.assertEqual(record[12].instances[6], "GCTGCTCCCCGCATACAGCG")
            self.assertEqual(record[12].instances[7], "GCCCCAAGCTTACCCAGGGA")
            self.assertEqual(record[12].instances[8], "ACCGGCTGACGCTAATACGG")
            self.assertEqual(record[12].instances[9], "GCGGGGGGCAGGTCATTACA")
            self.assertEqual(record[12].instances[10], "GCTGGCAGCGTCTAAGAAGG")
            self.assertEqual(record[12].instances[11], "GCAGGTGGTCGTGCAATACG")
            self.assertEqual(record[12].instances[12], "GCTGGTTGAAGTCCCGTGCG")
            self.assertEqual(record[12].instances[13], "GCACGTAGCTGGTAAATAGG")
            self.assertEqual(record[12].instances[14], "GCGGCGTGGATTTCATACAG")
            self.assertEqual(record[12].instances[15], "CCTGGAGGCTTAGACTTGGG")
        self.assertEqual(len(record[12].alignment.sequences), 16)
        self.assertEqual(record[12].alignment.sequences[0], "GCCGACGGGTGGTCATCGGG")
        self.assertEqual(record[12].alignment.sequences[1], "GCACGACGCTTTGTACCTGG")
        self.assertEqual(record[12].alignment.sequences[2], "CCTGGGAGGGTTCAATAACG")
        self.assertEqual(record[12].alignment.sequences[3], "GCGCGTCCCGGGCCAATAGC")
        self.assertEqual(record[12].alignment.sequences[4], "GCCGTCTGCTGGACAACGGG")
        self.assertEqual(record[12].alignment.sequences[5], "GTCCCTTCCGGTACATGAGG")
        self.assertEqual(record[12].alignment.sequences[6], "GCTGCTCCCCGCATACAGCG")
        self.assertEqual(record[12].alignment.sequences[7], "GCCCCAAGCTTACCCAGGGA")
        self.assertEqual(record[12].alignment.sequences[8], "ACCGGCTGACGCTAATACGG")
        self.assertEqual(record[12].alignment.sequences[9], "GCGGGGGGCAGGTCATTACA")
        self.assertEqual(record[12].alignment.sequences[10], "GCTGGCAGCGTCTAAGAAGG")
        self.assertEqual(record[12].alignment.sequences[11], "GCAGGTGGTCGTGCAATACG")
        self.assertEqual(record[12].alignment.sequences[12], "GCTGGTTGAAGTCCCGTGCG")
        self.assertEqual(record[12].alignment.sequences[13], "GCACGTAGCTGGTAAATAGG")
        self.assertEqual(record[12].alignment.sequences[14], "GCGGCGTGGATTTCATACAG")
        self.assertEqual(record[12].alignment.sequences[15], "CCTGGAGGCTTAGACTTGGG")
        self.assertEqual(
            record[12].mask,
            (1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1),
        )
        self.assertAlmostEqual(record[12].score, 5.63667)
        self.assertEqual(
            str(record[12]),
            """\
GCCGACGGGTGGTCATCGGG
GCACGACGCTTTGTACCTGG
CCTGGGAGGGTTCAATAACG
GCGCGTCCCGGGCCAATAGC
GCCGTCTGCTGGACAACGGG
GTCCCTTCCGGTACATGAGG
GCTGCTCCCCGCATACAGCG
GCCCCAAGCTTACCCAGGGA
ACCGGCTGACGCTAATACGG
GCGGGGGGCAGGTCATTACA
GCTGGCAGCGTCTAAGAAGG
GCAGGTGGTCGTGCAATACG
GCTGGTTGAAGTCCCGTGCG
GCACGTAGCTGGTAAATAGG
GCGGCGTGGATTTCATACAG
CCTGGAGGCTTAGACTTGGG""",
        )
        motif = record[12][2:-1]
        self.assertEqual(motif.alphabet, "ACGT")
        self.assertEqual(
            str(motif),
            """\
CGACGGGTGGTCATCGG
ACGACGCTTTGTACCTG
TGGGAGGGTTCAATAAC
GCGTCCCGGGCCAATAG
CGTCTGCTGGACAACGG
CCCTTCCGGTACATGAG
TGCTCCCCGCATACAGC
CCCAAGCTTACCCAGGG
CGGCTGACGCTAATACG
GGGGGGCAGGTCATTAC
TGGCAGCGTCTAAGAAG
AGGTGGTCGTGCAATAC
TGGTTGAAGTCCCGTGC
ACGTAGCTGGTAAATAG
GGCGTGGATTTCATACA
TGGAGGCTTAGACTTGG""",
        )
        self.assertEqual(record[13].alphabet, "ACGT")
        # using the old instances property:
        with self.assertWarns(BiopythonDeprecationWarning):
            self.assertEqual(len(record[13].instances), 15)
            self.assertEqual(record[13].instances[0], "GCCGACGGGTGGTCATCGGG")
            self.assertEqual(record[13].instances[1], "ATCCGCGGACGCTTAGAGGG")
            self.assertEqual(record[13].instances[2], "ACGCTTTGTACCTGGCTTGC")
            self.assertEqual(record[13].instances[3], "ACGGACGGCACTTAGCAGCG")
            self.assertEqual(record[13].instances[4], "GCCGTCTGCTGGACAACGGG")
            self.assertEqual(record[13].instances[5], "ACACACAGACGGTTGAAAGG")
            self.assertEqual(record[13].instances[6], "GCCGATAGTGCTTAAGTTCG")
            self.assertEqual(record[13].instances[7], "CTTGCCCGTACCGGACCTCG")
            self.assertEqual(record[13].instances[8], "ACCGGCTGACGCTAATACGG")
            self.assertEqual(record[13].instances[9], "GCCCCCCGCATAGTAGGGGG")
            self.assertEqual(record[13].instances[10], "GCTGGCAGCGTCTAAGAAGG")
            self.assertEqual(record[13].instances[11], "GCAGGTGGTCGTGCAATACG")
            self.assertEqual(record[13].instances[12], "ACGCACGGGACTTCAACCAG")
            self.assertEqual(record[13].instances[13], "GCACGTAGCTGGTAAATAGG")
            self.assertEqual(record[13].instances[14], "ATCCTCTGCGTCGCATGGCG")
        self.assertEqual(len(record[13].alignment.sequences), 15)
        self.assertEqual(record[13].alignment.sequences[0], "GCCGACGGGTGGTCATCGGG")
        self.assertEqual(record[13].alignment.sequences[1], "ATCCGCGGACGCTTAGAGGG")
        self.assertEqual(record[13].alignment.sequences[2], "ACGCTTTGTACCTGGCTTGC")
        self.assertEqual(record[13].alignment.sequences[3], "ACGGACGGCACTTAGCAGCG")
        self.assertEqual(record[13].alignment.sequences[4], "GCCGTCTGCTGGACAACGGG")
        self.assertEqual(record[13].alignment.sequences[5], "ACACACAGACGGTTGAAAGG")
        self.assertEqual(record[13].alignment.sequences[6], "GCCGATAGTGCTTAAGTTCG")
        self.assertEqual(record[13].alignment.sequences[7], "CTTGCCCGTACCGGACCTCG")
        self.assertEqual(record[13].alignment.sequences[8], "ACCGGCTGACGCTAATACGG")
        self.assertEqual(record[13].alignment.sequences[9], "GCCCCCCGCATAGTAGGGGG")
        self.assertEqual(record[13].alignment.sequences[10], "GCTGGCAGCGTCTAAGAAGG")
        self.assertEqual(record[13].alignment.sequences[11], "GCAGGTGGTCGTGCAATACG")
        self.assertEqual(record[13].alignment.sequences[12], "ACGCACGGGACTTCAACCAG")
        self.assertEqual(record[13].alignment.sequences[13], "GCACGTAGCTGGTAAATAGG")
        self.assertEqual(record[13].alignment.sequences[14], "ATCCTCTGCGTCGCATGGCG")
        self.assertEqual(
            record[13].mask,
            (1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1),
        )
        self.assertAlmostEqual(record[13].score, 3.89842)
        self.assertEqual(
            str(record[13]),
            """\
GCCGACGGGTGGTCATCGGG
ATCCGCGGACGCTTAGAGGG
ACGCTTTGTACCTGGCTTGC
ACGGACGGCACTTAGCAGCG
GCCGTCTGCTGGACAACGGG
ACACACAGACGGTTGAAAGG
GCCGATAGTGCTTAAGTTCG
CTTGCCCGTACCGGACCTCG
ACCGGCTGACGCTAATACGG
GCCCCCCGCATAGTAGGGGG
GCTGGCAGCGTCTAAGAAGG
GCAGGTGGTCGTGCAATACG
ACGCACGGGACTTCAACCAG
GCACGTAGCTGGTAAATAGG
ATCCTCTGCGTCGCATGGCG""",
        )
        motif = record[13][2:-1]
        self.assertEqual(motif.alphabet, "ACGT")
        self.assertEqual(
            str(motif),
            """\
CGACGGGTGGTCATCGG
CCGCGGACGCTTAGAGG
GCTTTGTACCTGGCTTG
GGACGGCACTTAGCAGC
CGTCTGCTGGACAACGG
ACACAGACGGTTGAAAG
CGATAGTGCTTAAGTTC
TGCCCGTACCGGACCTC
CGGCTGACGCTAATACG
CCCCCGCATAGTAGGGG
TGGCAGCGTCTAAGAAG
AGGTGGTCGTGCAATAC
GCACGGGACTTCAACCA
ACGTAGCTGGTAAATAG
CCTCTGCGTCGCATGGC""",
        )
        self.assertEqual(record[14].alphabet, "ACGT")
        # using the old instances property:
        with self.assertWarns(BiopythonDeprecationWarning):
            self.assertEqual(len(record[14].instances), 14)
            self.assertEqual(record[14].instances[0], "GAGGCTGTGTAT")
            self.assertEqual(record[14].instances[1], "GAGGTCGGGGGT")
            self.assertEqual(record[14].instances[2], "GACGGACGGCAC")
            self.assertEqual(record[14].instances[3], "TTGGCCCGGGAC")
            self.assertEqual(record[14].instances[4], "GAGGCTCGGCCC")
            self.assertEqual(record[14].instances[5], "CACGCGCTGTAT")
            self.assertEqual(record[14].instances[6], "TAGGCCAGGTAT")
            self.assertEqual(record[14].instances[7], "GAGGTCCGGTAC")
            self.assertEqual(record[14].instances[8], "TACGCTGGGGAT")
            self.assertEqual(record[14].instances[9], "GTCGCGGAGGAT")
            self.assertEqual(record[14].instances[10], "TACGCACGGGAC")
            self.assertEqual(record[14].instances[11], "TACTCCGGGTAC")
            self.assertEqual(record[14].instances[12], "GACGCAGAGGAT")
            self.assertEqual(record[14].instances[13], "TAGGCGGGCCAT")
        self.assertEqual(len(record[14].alignment.sequences), 14)
        self.assertEqual(record[14].alignment.sequences[0], "GAGGCTGTGTAT")
        self.assertEqual(record[14].alignment.sequences[1], "GAGGTCGGGGGT")
        self.assertEqual(record[14].alignment.sequences[2], "GACGGACGGCAC")
        self.assertEqual(record[14].alignment.sequences[3], "TTGGCCCGGGAC")
        self.assertEqual(record[14].alignment.sequences[4], "GAGGCTCGGCCC")
        self.assertEqual(record[14].alignment.sequences[5], "CACGCGCTGTAT")
        self.assertEqual(record[14].alignment.sequences[6], "TAGGCCAGGTAT")
        self.assertEqual(record[14].alignment.sequences[7], "GAGGTCCGGTAC")
        self.assertEqual(record[14].alignment.sequences[8], "TACGCTGGGGAT")
        self.assertEqual(record[14].alignment.sequences[9], "GTCGCGGAGGAT")
        self.assertEqual(record[14].alignment.sequences[10], "TACGCACGGGAC")
        self.assertEqual(record[14].alignment.sequences[11], "TACTCCGGGTAC")
        self.assertEqual(record[14].alignment.sequences[12], "GACGCAGAGGAT")
        self.assertEqual(record[14].alignment.sequences[13], "TAGGCGGGCCAT")
        self.assertEqual(record[14].mask, (1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1))
        self.assertAlmostEqual(record[14].score, 3.33444)
        self.assertEqual(
            str(record[14]),
            """\
GAGGCTGTGTAT
GAGGTCGGGGGT
GACGGACGGCAC
TTGGCCCGGGAC
GAGGCTCGGCCC
CACGCGCTGTAT
TAGGCCAGGTAT
GAGGTCCGGTAC
TACGCTGGGGAT
GTCGCGGAGGAT
TACGCACGGGAC
TACTCCGGGTAC
GACGCAGAGGAT
TAGGCGGGCCAT""",
        )
        motif = record[14][2:-1]
        self.assertEqual(motif.alphabet, "ACGT")
        self.assertEqual(
            str(motif),
            """\
GGCTGTGTA
GGTCGGGGG
CGGACGGCA
GGCCCGGGA
GGCTCGGCC
CGCGCTGTA
GGCCAGGTA
GGTCCGGTA
CGCTGGGGA
CGCGGAGGA
CGCACGGGA
CTCCGGGTA
CGCAGAGGA
GGCGGGCCA""",
        )
        self.assertEqual(record[15].alphabet, "ACGT")
        # using the old instances property:
        with self.assertWarns(BiopythonDeprecationWarning):
            self.assertEqual(len(record[15].instances), 21)
            self.assertEqual(record[15].instances[0], "CGGCTCAATCGTAGAGGC")
            self.assertEqual(record[15].instances[1], "CGACGGGTGGTCATCGGG")
            self.assertEqual(record[15].instances[2], "CGCTTAGAGGGCACAAGC")
            self.assertEqual(record[15].instances[3], "TGACACGCGCCTGGGAGG")
            self.assertEqual(record[15].instances[4], "CGATACGCTGCTAAGTGC")
            self.assertEqual(record[15].instances[5], "CGTCCCGGGCCAATAGCG")
            self.assertEqual(record[15].instances[6], "CCACGCTTCGACACGTGG")
            self.assertEqual(record[15].instances[7], "CGTCTGCTGGACAACGGG")
            self.assertEqual(record[15].instances[8], "ACACAGACGGTTGAAAGG")
            self.assertEqual(record[15].instances[9], "TGCTCCCCGCATACAGCG")
            self.assertEqual(record[15].instances[10], "TGAGGCTTGCCCGTACCG")
            self.assertEqual(record[15].instances[11], "TGCCCCAAGCTTACCCAG")
            self.assertEqual(record[15].instances[12], "CGGCTGACGCTAATACGG")
            self.assertEqual(record[15].instances[13], "CGCGACGTCCCTATGAGC")
            self.assertEqual(record[15].instances[14], "TGCCCCCCGCATAGTAGG")
            self.assertEqual(record[15].instances[15], "CGTTGCCTTCTTAGACGC")
            self.assertEqual(record[15].instances[16], "TGACTCAATCGTAGACCC")
            self.assertEqual(record[15].instances[17], "AGTCCCGTGCGTATGTGG")
            self.assertEqual(record[15].instances[18], "AGGCTCGCACGTAGCTGG")
            self.assertEqual(record[15].instances[19], "CCACGCCGCCATGCGACG")
            self.assertEqual(record[15].instances[20], "AGCCTCCAGGTCGCATGG")
        self.assertEqual(len(record[15].alignment.sequences), 21)
        self.assertEqual(record[15].alignment.sequences[0], "CGGCTCAATCGTAGAGGC")
        self.assertEqual(record[15].alignment.sequences[1], "CGACGGGTGGTCATCGGG")
        self.assertEqual(record[15].alignment.sequences[2], "CGCTTAGAGGGCACAAGC")
        self.assertEqual(record[15].alignment.sequences[3], "TGACACGCGCCTGGGAGG")
        self.assertEqual(record[15].alignment.sequences[4], "CGATACGCTGCTAAGTGC")
        self.assertEqual(record[15].alignment.sequences[5], "CGTCCCGGGCCAATAGCG")
        self.assertEqual(record[15].alignment.sequences[6], "CCACGCTTCGACACGTGG")
        self.assertEqual(record[15].alignment.sequences[7], "CGTCTGCTGGACAACGGG")
        self.assertEqual(record[15].alignment.sequences[8], "ACACAGACGGTTGAAAGG")
        self.assertEqual(record[15].alignment.sequences[9], "TGCTCCCCGCATACAGCG")
        self.assertEqual(record[15].alignment.sequences[10], "TGAGGCTTGCCCGTACCG")
        self.assertEqual(record[15].alignment.sequences[11], "TGCCCCAAGCTTACCCAG")
        self.assertEqual(record[15].alignment.sequences[12], "CGGCTGACGCTAATACGG")
        self.assertEqual(record[15].alignment.sequences[13], "CGCGACGTCCCTATGAGC")
        self.assertEqual(record[15].alignment.sequences[14], "TGCCCCCCGCATAGTAGG")
        self.assertEqual(record[15].alignment.sequences[15], "CGTTGCCTTCTTAGACGC")
        self.assertEqual(record[15].alignment.sequences[16], "TGACTCAATCGTAGACCC")
        self.assertEqual(record[15].alignment.sequences[17], "AGTCCCGTGCGTATGTGG")
        self.assertEqual(record[15].alignment.sequences[18], "AGGCTCGCACGTAGCTGG")
        self.assertEqual(record[15].alignment.sequences[19], "CCACGCCGCCATGCGACG")
        self.assertEqual(record[15].alignment.sequences[20], "AGCCTCCAGGTCGCATGG")
        self.assertEqual(
            record[15].mask, (1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1)
        )
        self.assertAlmostEqual(record[15].score, 1.0395)
        self.assertEqual(
            str(record[15]),
            """\
CGGCTCAATCGTAGAGGC
CGACGGGTGGTCATCGGG
CGCTTAGAGGGCACAAGC
TGACACGCGCCTGGGAGG
CGATACGCTGCTAAGTGC
CGTCCCGGGCCAATAGCG
CCACGCTTCGACACGTGG
CGTCTGCTGGACAACGGG
ACACAGACGGTTGAAAGG
TGCTCCCCGCATACAGCG
TGAGGCTTGCCCGTACCG
TGCCCCAAGCTTACCCAG
CGGCTGACGCTAATACGG
CGCGACGTCCCTATGAGC
TGCCCCCCGCATAGTAGG
CGTTGCCTTCTTAGACGC
TGACTCAATCGTAGACCC
AGTCCCGTGCGTATGTGG
AGGCTCGCACGTAGCTGG
CCACGCCGCCATGCGACG
AGCCTCCAGGTCGCATGG""",
        )
        motif = record[15][2:-1]
        self.assertEqual(motif.alphabet, "ACGT")
        self.assertEqual(
            str(motif),
            """\
GCTCAATCGTAGAGG
ACGGGTGGTCATCGG
CTTAGAGGGCACAAG
ACACGCGCCTGGGAG
ATACGCTGCTAAGTG
TCCCGGGCCAATAGC
ACGCTTCGACACGTG
TCTGCTGGACAACGG
ACAGACGGTTGAAAG
CTCCCCGCATACAGC
AGGCTTGCCCGTACC
CCCCAAGCTTACCCA
GCTGACGCTAATACG
CGACGTCCCTATGAG
CCCCCCGCATAGTAG
TTGCCTTCTTAGACG
ACTCAATCGTAGACC
TCCCGTGCGTATGTG
GCTCGCACGTAGCTG
ACGCCGCCATGCGAC
CCTCCAGGTCGCATG""",
        )


class TestClusterBuster(unittest.TestCase):
    """Testing parsing Cluster-Buster output files."""

    def test_clusterbuster_parsing_and_output(self):
        """Test if Bio.motifs can parse and output Cluster-Buster PFM files."""
        with open("motifs/clusterbuster.pfm") as stream:
            record = motifs.parse(stream, "clusterbuster")
            self.assertEqual(len(record), 3)
            motif = record[0]
            self.assertEqual(motif.name, "MA0004.1")
            self.assertEqual(motif.alphabet, "GATC")
            self.assertEqual(motif.consensus, "CACGTG")
            self.assertEqual(motif.degenerate_consensus, "CACGTG")
            self.assertTrue(
                np.allclose(
                    motif.relative_entropy,
                    np.array(
                        [1.278071905112638, 1.7136030428840439, 2.0, 2.0, 2.0, 2.0]
                    ),
                )
            )
            self.assertEqual(motif[1:-2].consensus, "ACG")
            self.assertEqual(motif.length, 6)
            self.assertAlmostEqual(motif.counts["G", 0], 0.0)
            self.assertAlmostEqual(motif.counts["G", 1], 1.0)
            self.assertAlmostEqual(motif.counts["G", 2], 0.0)
            self.assertAlmostEqual(motif.counts["G", 3], 20.0)
            self.assertAlmostEqual(motif.counts["G", 4], 0.0)
            self.assertAlmostEqual(motif.counts["G", 5], 20.0)
            self.assertAlmostEqual(motif.counts["A", 0], 4.0)
            self.assertAlmostEqual(motif.counts["A", 1], 19.0)
            self.assertAlmostEqual(motif.counts["A", 2], 0.0)
            self.assertAlmostEqual(motif.counts["A", 3], 0.0)
            self.assertAlmostEqual(motif.counts["A", 4], 0.0)
            self.assertAlmostEqual(motif.counts["A", 5], 0.0)
            self.assertAlmostEqual(motif.counts["T", 0], 0.0)
            self.assertAlmostEqual(motif.counts["T", 1], 0.0)
            self.assertAlmostEqual(motif.counts["T", 2], 0.0)
            self.assertAlmostEqual(motif.counts["T", 3], 0.0)
            self.assertAlmostEqual(motif.counts["T", 4], 20.0)
            self.assertAlmostEqual(motif.counts["T", 5], 0.0)
            self.assertAlmostEqual(motif.counts["C", 0], 16.0)
            self.assertAlmostEqual(motif.counts["C", 1], 0.0)
            self.assertAlmostEqual(motif.counts["C", 2], 20.0)
            self.assertAlmostEqual(motif.counts["C", 3], 0.0)
            self.assertAlmostEqual(motif.counts["C", 4], 0.0)
            self.assertAlmostEqual(motif.counts["C", 5], 0.0)
            motif = record[1]
            self.assertEqual(motif.name, "MA0006.1")
            self.assertEqual(motif.alphabet, "GATC")
            self.assertEqual(motif.consensus, "TGCGTG")
            self.assertEqual(motif.degenerate_consensus, "YGCGTG")
            self.assertTrue(
                np.allclose(
                    motif.relative_entropy,
                    np.array(
                        [
                            0.28206397041108283,
                            1.7501177071668148,
                            1.7501177071668148,
                            1.7501177071668148,
                            2.0,
                            2.0,
                        ]
                    ),
                )
            )
            self.assertEqual(motif[1:-2].consensus, "GCG")
            self.assertEqual(motif.length, 6)
            self.assertAlmostEqual(motif.counts["G", 0], 2.0)
            self.assertAlmostEqual(motif.counts["G", 1], 23.0)
            self.assertAlmostEqual(motif.counts["G", 2], 0.0)
            self.assertAlmostEqual(motif.counts["G", 3], 23.0)
            self.assertAlmostEqual(motif.counts["G", 4], 0.0)
            self.assertAlmostEqual(motif.counts["G", 5], 24.0)
            self.assertAlmostEqual(motif.counts["A", 0], 3.0)
            self.assertAlmostEqual(motif.counts["A", 1], 0.0)
            self.assertAlmostEqual(motif.counts["A", 2], 0.0)
            self.assertAlmostEqual(motif.counts["A", 3], 0.0)
            self.assertAlmostEqual(motif.counts["A", 4], 0.0)
            self.assertAlmostEqual(motif.counts["A", 5], 0.0)
            self.assertAlmostEqual(motif.counts["T", 0], 11.0)
            self.assertAlmostEqual(motif.counts["T", 1], 1.0)
            self.assertAlmostEqual(motif.counts["T", 2], 1.0)
            self.assertAlmostEqual(motif.counts["T", 3], 1.0)
            self.assertAlmostEqual(motif.counts["T", 4], 24.0)
            self.assertAlmostEqual(motif.counts["T", 5], 0.0)
            self.assertAlmostEqual(motif.counts["C", 0], 8.0)
            self.assertAlmostEqual(motif.counts["C", 1], 0.0)
            self.assertAlmostEqual(motif.counts["C", 2], 23.0)
            self.assertAlmostEqual(motif.counts["C", 3], 0.0)
            self.assertAlmostEqual(motif.counts["C", 4], 0.0)
            self.assertAlmostEqual(motif.counts["C", 5], 0.0)
            motif = record[2]
            self.assertEqual(motif.name, "MA0008.1")
            self.assertEqual(motif.alphabet, "GATC")
            self.assertEqual(motif.consensus, "CAATTATT")
            self.assertEqual(motif.degenerate_consensus, "CAATTATT")
            self.assertTrue(
                np.allclose(
                    motif.relative_entropy,
                    np.array(
                        [
                            0.2549535827226545,
                            1.2358859454459725,
                            2.0,
                            2.0,
                            1.278071905112638,
                            1.7577078109175852,
                            1.7577078109175852,
                            1.5978208097977271,
                        ]
                    ),
                )
            )
            self.assertEqual(motif[1:-2].consensus, "AATTA")
            self.assertEqual(motif.length, 8)
            self.assertAlmostEqual(motif.counts["G", 0], 4.0)
            self.assertAlmostEqual(motif.counts["G", 1], 0.0)
            self.assertAlmostEqual(motif.counts["G", 2], 0.0)
            self.assertAlmostEqual(motif.counts["G", 3], 0.0)
            self.assertAlmostEqual(motif.counts["G", 4], 0.0)
            self.assertAlmostEqual(motif.counts["G", 5], 1.0)
            self.assertAlmostEqual(motif.counts["G", 6], 0.0)
            self.assertAlmostEqual(motif.counts["G", 7], 2.0)
            self.assertAlmostEqual(motif.counts["A", 0], 3.0)
            self.assertAlmostEqual(motif.counts["A", 1], 21.0)
            self.assertAlmostEqual(motif.counts["A", 2], 25.0)
            self.assertAlmostEqual(motif.counts["A", 3], 0.0)
            self.assertAlmostEqual(motif.counts["A", 4], 0.0)
            self.assertAlmostEqual(motif.counts["A", 5], 24.0)
            self.assertAlmostEqual(motif.counts["A", 6], 1.0)
            self.assertAlmostEqual(motif.counts["A", 7], 0.0)
            self.assertAlmostEqual(motif.counts["T", 0], 5.0)
            self.assertAlmostEqual(motif.counts["T", 1], 3.0)
            self.assertAlmostEqual(motif.counts["T", 2], 0.0)
            self.assertAlmostEqual(motif.counts["T", 3], 25.0)
            self.assertAlmostEqual(motif.counts["T", 4], 20.0)
            self.assertAlmostEqual(motif.counts["T", 5], 0.0)
            self.assertAlmostEqual(motif.counts["T", 6], 24.0)
            self.assertAlmostEqual(motif.counts["T", 7], 23.0)
            self.assertAlmostEqual(motif.counts["C", 0], 13.0)
            self.assertAlmostEqual(motif.counts["C", 1], 1.0)
            self.assertAlmostEqual(motif.counts["C", 2], 0.0)
            self.assertAlmostEqual(motif.counts["C", 3], 0.0)
            self.assertAlmostEqual(motif.counts["C", 4], 5.0)
            self.assertAlmostEqual(motif.counts["C", 5], 0.0)
            self.assertAlmostEqual(motif.counts["C", 6], 0.0)
            self.assertAlmostEqual(motif.counts["C", 7], 0.0)
            stream.seek(0)
            self.assertEqual(
                motifs.write(record, "clusterbuster").split(),
                stream.read().split(),
            )


class TestXMS(unittest.TestCase):
    """Testing parsing xms output files."""

    def test_xms_parsing(self):
        """Test if Bio.motifs can parse and output xms PFM files."""
        with open("motifs/abdb.xms") as stream:
            record = motifs.parse(stream, "xms")
        self.assertEqual(len(record), 1)
        motif = record[0]
        self.assertEqual(motif.name, "Abd-B")
        self.assertEqual(motif.length, 14)
        self.assertEqual(motif.alphabet, "GATC")
        self.assertAlmostEqual(motif.counts["G", 0], 0.333333333)
        self.assertAlmostEqual(motif.counts["G", 1], 0.379310345)
        self.assertAlmostEqual(motif.counts["G", 2], 0.264705882)
        self.assertAlmostEqual(motif.counts["G", 3], 0.194444444)
        self.assertAlmostEqual(motif.counts["G", 4], 0.102564103)
        self.assertAlmostEqual(motif.counts["G", 5], 0.177777778)
        self.assertAlmostEqual(motif.counts["G", 6], 0.000000000)
        self.assertAlmostEqual(motif.counts["G", 7], 0.022222222)
        self.assertAlmostEqual(motif.counts["G", 8], 0.697674419)
        self.assertAlmostEqual(motif.counts["G", 9], 0.571428571)
        self.assertAlmostEqual(motif.counts["G", 10], 0.150000000)
        self.assertAlmostEqual(motif.counts["G", 11], 0.305555556)
        self.assertAlmostEqual(motif.counts["G", 12], 0.258064516)
        self.assertAlmostEqual(motif.counts["G", 13], 0.259259259)
        self.assertAlmostEqual(motif.counts["A", 0], 0.333333333)
        self.assertAlmostEqual(motif.counts["A", 1], 0.103448276)
        self.assertAlmostEqual(motif.counts["A", 2], 0.264705882)
        self.assertAlmostEqual(motif.counts["A", 3], 0.000000000)
        self.assertAlmostEqual(motif.counts["A", 4], 0.102564103)
        self.assertAlmostEqual(motif.counts["A", 5], 0.244444444)
        self.assertAlmostEqual(motif.counts["A", 6], 0.800000000)
        self.assertAlmostEqual(motif.counts["A", 7], 0.133333333)
        self.assertAlmostEqual(motif.counts["A", 8], 0.046511628)
        self.assertAlmostEqual(motif.counts["A", 9], 0.238095238)
        self.assertAlmostEqual(motif.counts["A", 10], 0.025000000)
        self.assertAlmostEqual(motif.counts["A", 11], 0.222222222)
        self.assertAlmostEqual(motif.counts["A", 12], 0.354838710)
        self.assertAlmostEqual(motif.counts["A", 13], 0.185185185)
        self.assertAlmostEqual(motif.counts["T", 0], 0.125000000)
        self.assertAlmostEqual(motif.counts["T", 1], 0.103448276)
        self.assertAlmostEqual(motif.counts["T", 2], 0.205882353)
        self.assertAlmostEqual(motif.counts["T", 3], 0.777777778)
        self.assertAlmostEqual(motif.counts["T", 4], 0.743589744)
        self.assertAlmostEqual(motif.counts["T", 5], 0.533333333)
        self.assertAlmostEqual(motif.counts["T", 6], 0.155555556)
        self.assertAlmostEqual(motif.counts["T", 7], 0.688888889)
        self.assertAlmostEqual(motif.counts["T", 8], 0.209302326)
        self.assertAlmostEqual(motif.counts["T", 9], 0.095238095)
        self.assertAlmostEqual(motif.counts["T", 10], 0.025000000)
        self.assertAlmostEqual(motif.counts["T", 11], 0.194444444)
        self.assertAlmostEqual(motif.counts["T", 12], 0.129032258)
        self.assertAlmostEqual(motif.counts["T", 13], 0.222222222)
        self.assertAlmostEqual(motif.counts["C", 0], 0.208333333)
        self.assertAlmostEqual(motif.counts["C", 1], 0.413793103)
        self.assertAlmostEqual(motif.counts["C", 2], 0.264705882)
        self.assertAlmostEqual(motif.counts["C", 3], 0.027777778)
        self.assertAlmostEqual(motif.counts["C", 4], 0.051282051)
        self.assertAlmostEqual(motif.counts["C", 5], 0.044444444)
        self.assertAlmostEqual(motif.counts["C", 6], 0.044444444)
        self.assertAlmostEqual(motif.counts["C", 7], 0.155555556)
        self.assertAlmostEqual(motif.counts["C", 8], 0.046511628)
        self.assertAlmostEqual(motif.counts["C", 9], 0.095238095)
        self.assertAlmostEqual(motif.counts["C", 10], 0.800000000)
        self.assertAlmostEqual(motif.counts["C", 11], 0.277777778)
        self.assertAlmostEqual(motif.counts["C", 12], 0.258064516)
        self.assertAlmostEqual(motif.counts["C", 13], 0.333333333)
        self.assertEqual(motif.consensus, "GCGTTTATGGCGAC")
        self.assertEqual(motif.degenerate_consensus, "NSNTTTATGGCNNN")
        self.assertTrue(
            np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        0.09689283163718865,
                        0.26557323997556864,
                        0.007815379142180268,
                        1.1150033950025815,
                        0.78848108520697,
                        0.3768768552773923,
                        1.125231003810913,
                        0.7023990165752877,
                        0.7536432192801433,
                        0.3995487907017483,
                        1.0658162802208113,
                        0.022422587676774776,
                        0.07979555429087543,
                        0.03400971806422712,
                    ]
                ),
            )
        )
        self.assertEqual(motif[3::2].consensus, "TTTGGC")
        self.assertEqual(motif[3::2].degenerate_consensus, "TTTGNN")
        self.assertTrue(
            np.allclose(
                motif[3::2].relative_entropy,
                np.array(
                    [
                        1.1150033950025815,
                        0.3768768552773923,
                        0.7023990165752877,
                        0.3995487907017483,
                        0.022422587676774776,
                        0.03400971806422712,
                    ]
                ),
            )
        )


class TestJASPAR(unittest.TestCase):
    """Testing parsing JASPAR files."""

    def test_pfm_parsing(self):
        """Test if Bio.motifs can parse JASPAR-style pfm files."""
        with open("motifs/SRF.pfm") as stream:
            m = motifs.read(stream, "pfm")
        self.assertEqual(m.length, 12)

    def test_pfm_four_columns_parsing(self):
        """Test if Bio.motifs.pfm can parse motifs in position frequency matrix format (4 columns)."""
        with open("motifs/fourcolumns.pfm") as stream:
            record = motifs.parse(stream, "pfm-four-columns")
        self.assertEqual(len(record), 8)
        motif = record[0]
        self.assertIsNone(motif.name)
        self.assertEqual(motif.length, 8)
        self.assertEqual(motif.alphabet, "GATC")
        self.assertAlmostEqual(motif.counts["G", 0], 0.009615385)
        self.assertAlmostEqual(motif.counts["G", 1], 0.009615385)
        self.assertAlmostEqual(motif.counts["G", 2], 0.009615385)
        self.assertAlmostEqual(motif.counts["G", 3], 0.009615385)
        self.assertAlmostEqual(motif.counts["G", 4], 0.009615385)
        self.assertAlmostEqual(motif.counts["G", 5], 0.009615385)
        self.assertAlmostEqual(motif.counts["G", 6], 0.009615385)
        self.assertAlmostEqual(motif.counts["G", 7], 0.009615385)
        self.assertAlmostEqual(motif.counts["A", 0], 0.009615385)
        self.assertAlmostEqual(motif.counts["A", 1], 0.009615385)
        self.assertAlmostEqual(motif.counts["A", 2], 0.971153846)
        self.assertAlmostEqual(motif.counts["A", 3], 0.009615385)
        self.assertAlmostEqual(motif.counts["A", 4], 0.009615385)
        self.assertAlmostEqual(motif.counts["A", 5], 0.971153846)
        self.assertAlmostEqual(motif.counts["A", 6], 0.009615385)
        self.assertAlmostEqual(motif.counts["A", 7], 0.009615385)
        self.assertAlmostEqual(motif.counts["T", 0], 0.971153846)
        self.assertAlmostEqual(motif.counts["T", 1], 0.971153846)
        self.assertAlmostEqual(motif.counts["T", 2], 0.009615385)
        self.assertAlmostEqual(motif.counts["T", 3], 0.971153846)
        self.assertAlmostEqual(motif.counts["T", 4], 0.009615385)
        self.assertAlmostEqual(motif.counts["T", 5], 0.009615385)
        self.assertAlmostEqual(motif.counts["T", 6], 0.009615385)
        self.assertAlmostEqual(motif.counts["T", 7], 0.971153846)
        self.assertAlmostEqual(motif.counts["C", 0], 0.009615385)
        self.assertAlmostEqual(motif.counts["C", 1], 0.009615385)
        self.assertAlmostEqual(motif.counts["C", 2], 0.009615385)
        self.assertAlmostEqual(motif.counts["C", 3], 0.009615385)
        self.assertAlmostEqual(motif.counts["C", 4], 0.971153846)
        self.assertAlmostEqual(motif.counts["C", 5], 0.009615385)
        self.assertAlmostEqual(motif.counts["C", 6], 0.971153846)
        self.assertAlmostEqual(motif.counts["C", 7], 0.009615385)
        self.assertEqual(motif.consensus, "TTATCACT")
        self.assertEqual(motif.degenerate_consensus, "TTATCACT")
        self.assertTrue(
            np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        1.765707971839016,
                        1.765707971839016,
                        1.7657079718390165,
                        1.765707971839016,
                        1.7657079718390158,
                        1.7657079718390165,
                        1.7657079718390158,
                        1.765707971839016,
                    ]
                ),
            )
        )
        self.assertEqual(motif[1:-2].consensus, "TATCA")
        motif = record[1]
        self.assertEqual(motif.name, "ENSG00000197372")
        self.assertEqual(motif.length, 20)
        self.assertEqual(motif.alphabet, "GATC")
        self.assertAlmostEqual(motif.counts["G", 0], 0.117054000)
        self.assertAlmostEqual(motif.counts["G", 1], 0.364552000)
        self.assertAlmostEqual(motif.counts["G", 2], 0.310520000)
        self.assertAlmostEqual(motif.counts["G", 3], 0.131007000)
        self.assertAlmostEqual(motif.counts["G", 4], 0.176504000)
        self.assertAlmostEqual(motif.counts["G", 5], 0.197793000)
        self.assertAlmostEqual(motif.counts["G", 6], 0.926202000)
        self.assertAlmostEqual(motif.counts["G", 7], 0.983797000)
        self.assertAlmostEqual(motif.counts["G", 8], 0.002387000)
        self.assertAlmostEqual(motif.counts["G", 9], 0.002418000)
        self.assertAlmostEqual(motif.counts["G", 10], 0.001991000)
        self.assertAlmostEqual(motif.counts["G", 11], 0.002868000)
        self.assertAlmostEqual(motif.counts["G", 12], 0.350783000)
        self.assertAlmostEqual(motif.counts["G", 13], 1.000000000)
        self.assertAlmostEqual(motif.counts["G", 14], 0.000000000)
        self.assertAlmostEqual(motif.counts["G", 15], 1.000000000)
        self.assertAlmostEqual(motif.counts["G", 16], 1.000000000)
        self.assertAlmostEqual(motif.counts["G", 17], 0.000000000)
        self.assertAlmostEqual(motif.counts["G", 18], 0.000000000)
        self.assertAlmostEqual(motif.counts["G", 19], 0.000000000)
        self.assertAlmostEqual(motif.counts["A", 0], 0.341303000)
        self.assertAlmostEqual(motif.counts["A", 1], 0.283785000)
        self.assertAlmostEqual(motif.counts["A", 2], 0.491055000)
        self.assertAlmostEqual(motif.counts["A", 3], 0.492621000)
        self.assertAlmostEqual(motif.counts["A", 4], 0.250645000)
        self.assertAlmostEqual(motif.counts["A", 5], 0.276694000)
        self.assertAlmostEqual(motif.counts["A", 6], 0.056317000)
        self.assertAlmostEqual(motif.counts["A", 7], 0.004470000)
        self.assertAlmostEqual(motif.counts["A", 8], 0.936213000)
        self.assertAlmostEqual(motif.counts["A", 9], 0.004352000)
        self.assertAlmostEqual(motif.counts["A", 10], 0.013277000)
        self.assertAlmostEqual(motif.counts["A", 11], 0.968132000)
        self.assertAlmostEqual(motif.counts["A", 12], 0.397623000)
        self.assertAlmostEqual(motif.counts["A", 13], 0.000000000)
        self.assertAlmostEqual(motif.counts["A", 14], 1.000000000)
        self.assertAlmostEqual(motif.counts["A", 15], 0.000000000)
        self.assertAlmostEqual(motif.counts["A", 16], 0.000000000)
        self.assertAlmostEqual(motif.counts["A", 17], 1.000000000)
        self.assertAlmostEqual(motif.counts["A", 18], 0.000000000)
        self.assertAlmostEqual(motif.counts["A", 19], 1.000000000)
        self.assertAlmostEqual(motif.counts["T", 0], 0.409215000)
        self.assertAlmostEqual(motif.counts["T", 1], 0.274597000)
        self.assertAlmostEqual(motif.counts["T", 2], 0.120217000)
        self.assertAlmostEqual(motif.counts["T", 3], 0.300256000)
        self.assertAlmostEqual(motif.counts["T", 4], 0.211387000)
        self.assertAlmostEqual(motif.counts["T", 5], 0.027444000)
        self.assertAlmostEqual(motif.counts["T", 6], 0.002850000)
        self.assertAlmostEqual(motif.counts["T", 7], 0.003964000)
        self.assertAlmostEqual(motif.counts["T", 8], 0.002613000)
        self.assertAlmostEqual(motif.counts["T", 9], 0.989200000)
        self.assertAlmostEqual(motif.counts["T", 10], 0.976567000)
        self.assertAlmostEqual(motif.counts["T", 11], 0.026737000)
        self.assertAlmostEqual(motif.counts["T", 12], 0.199577000)
        self.assertAlmostEqual(motif.counts["T", 13], 0.000000000)
        self.assertAlmostEqual(motif.counts["T", 14], 0.000000000)
        self.assertAlmostEqual(motif.counts["T", 15], 0.000000000)
        self.assertAlmostEqual(motif.counts["T", 16], 0.000000000)
        self.assertAlmostEqual(motif.counts["T", 17], 0.000000000)
        self.assertAlmostEqual(motif.counts["T", 18], 0.000000000)
        self.assertAlmostEqual(motif.counts["T", 19], 0.000000000)
        self.assertAlmostEqual(motif.counts["C", 0], 0.132427000)
        self.assertAlmostEqual(motif.counts["C", 1], 0.077066000)
        self.assertAlmostEqual(motif.counts["C", 2], 0.078208000)
        self.assertAlmostEqual(motif.counts["C", 3], 0.076117000)
        self.assertAlmostEqual(motif.counts["C", 4], 0.361464000)
        self.assertAlmostEqual(motif.counts["C", 5], 0.498070000)
        self.assertAlmostEqual(motif.counts["C", 6], 0.014631000)
        self.assertAlmostEqual(motif.counts["C", 7], 0.007769000)
        self.assertAlmostEqual(motif.counts["C", 8], 0.058787000)
        self.assertAlmostEqual(motif.counts["C", 9], 0.004030000)
        self.assertAlmostEqual(motif.counts["C", 10], 0.008165000)
        self.assertAlmostEqual(motif.counts["C", 11], 0.002263000)
        self.assertAlmostEqual(motif.counts["C", 12], 0.052017000)
        self.assertAlmostEqual(motif.counts["C", 13], 0.000000000)
        self.assertAlmostEqual(motif.counts["C", 14], 0.000000000)
        self.assertAlmostEqual(motif.counts["C", 15], 0.000000000)
        self.assertAlmostEqual(motif.counts["C", 16], 0.000000000)
        self.assertAlmostEqual(motif.counts["C", 17], 0.000000000)
        self.assertAlmostEqual(motif.counts["C", 18], 1.000000000)
        self.assertAlmostEqual(motif.counts["C", 19], 0.000000000)
        self.assertEqual(motif.consensus, "TGAACCGGATTAAGAGGACA")
        self.assertEqual(motif.degenerate_consensus, "WNRWNMGGATTANGAGGACA")
        self.assertTrue(
            np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        0.1946677220077018,
                        0.1566211351816578,
                        0.31728135119311995,
                        0.3086747573287918,
                        0.053393542701508756,
                        0.381471417197324,
                        1.5505596169174871,
                        1.8558501430757017,
                        1.6274200195132635,
                        1.8972899364737197,
                        1.809312450637467,
                        1.7709547585539227,
                        0.2549373240046801,
                        2.0,
                        2.0,
                        2.0,
                        2.0,
                        2.0,
                        2.0,
                        2.0,
                    ]
                ),
            )
        )
        self.assertEqual(motif[1:-2].consensus, "GAACCGGATTAAGAGGA")
        motif = record[2]
        self.assertAlmostEqual(motif.counts["G", 0], 0.083333300)
        self.assertAlmostEqual(motif.counts["G", 1], 0.083333300)
        self.assertAlmostEqual(motif.counts["G", 2], 0.000000000)
        self.assertAlmostEqual(motif.counts["G", 3], 0.000000000)
        self.assertAlmostEqual(motif.counts["G", 4], 0.083333300)
        self.assertAlmostEqual(motif.counts["G", 5], 0.000000000)
        self.assertAlmostEqual(motif.counts["G", 6], 0.000000000)
        self.assertAlmostEqual(motif.counts["G", 7], 0.333333000)
        self.assertAlmostEqual(motif.counts["G", 8], 0.166667000)
        self.assertAlmostEqual(motif.counts["G", 9], 0.166667000)
        self.assertAlmostEqual(motif.counts["G", 10], 0.416667000)
        self.assertAlmostEqual(motif.counts["A", 0], 0.250000000)
        self.assertAlmostEqual(motif.counts["A", 1], 0.750000000)
        self.assertAlmostEqual(motif.counts["A", 2], 0.833333000)
        self.assertAlmostEqual(motif.counts["A", 3], 1.000000000)
        self.assertAlmostEqual(motif.counts["A", 4], 0.000000000)
        self.assertAlmostEqual(motif.counts["A", 5], 0.333333000)
        self.assertAlmostEqual(motif.counts["A", 6], 0.833333000)
        self.assertAlmostEqual(motif.counts["A", 7], 0.500000000)
        self.assertAlmostEqual(motif.counts["A", 8], 0.500000000)
        self.assertAlmostEqual(motif.counts["A", 9], 0.333333000)
        self.assertAlmostEqual(motif.counts["A", 10], 0.166667000)
        self.assertAlmostEqual(motif.counts["T", 0], 0.583333000)
        self.assertAlmostEqual(motif.counts["T", 1], 0.000000000)
        self.assertAlmostEqual(motif.counts["T", 2], 0.166667000)
        self.assertAlmostEqual(motif.counts["T", 3], 0.000000000)
        self.assertAlmostEqual(motif.counts["T", 4], 0.083333300)
        self.assertAlmostEqual(motif.counts["T", 5], 0.666667000)
        self.assertAlmostEqual(motif.counts["T", 6], 0.166667000)
        self.assertAlmostEqual(motif.counts["T", 7], 0.166667000)
        self.assertAlmostEqual(motif.counts["T", 8], 0.250000000)
        self.assertAlmostEqual(motif.counts["T", 9], 0.250000000)
        self.assertAlmostEqual(motif.counts["T", 10], 0.166667000)
        self.assertAlmostEqual(motif.counts["C", 0], 0.083333300)
        self.assertAlmostEqual(motif.counts["C", 1], 0.166667000)
        self.assertAlmostEqual(motif.counts["C", 2], 0.000000000)
        self.assertAlmostEqual(motif.counts["C", 3], 0.000000000)
        self.assertAlmostEqual(motif.counts["C", 4], 0.833333000)
        self.assertAlmostEqual(motif.counts["C", 5], 0.000000000)
        self.assertAlmostEqual(motif.counts["C", 6], 0.000000000)
        self.assertAlmostEqual(motif.counts["C", 7], 0.000000000)
        self.assertAlmostEqual(motif.counts["C", 8], 0.083333300)
        self.assertAlmostEqual(motif.counts["C", 9], 0.250000000)
        self.assertAlmostEqual(motif.counts["C", 10], 0.250000000)
        self.assertEqual(motif.name, "M1734_0.90")
        self.assertEqual(motif.length, 11)
        self.assertEqual(motif.alphabet, "GATC")
        self.assertEqual(motif.consensus, "TAAACTAAAAG")
        self.assertEqual(motif.degenerate_consensus, "TAAACTARNNN")
        self.assertTrue(
            np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        0.4489017067534855,
                        0.9591474871280075,
                        1.3499768043761913,
                        2.0,
                        1.1833109116849791,
                        1.0817044992792044,
                        1.3499768043761913,
                        0.5408517496401433,
                        0.2704258182036411,
                        0.04085174964014324,
                        0.1120812409282564,
                    ]
                ),
            )
        )
        self.assertEqual(motif[1:-2].consensus, "AAACTAAA")
        motif = record[3]
        self.assertEqual(motif.name, "AbdA_Cell_FBgn0000014")
        self.assertEqual(motif.length, 7)
        self.assertEqual(motif.alphabet, "GATC")
        self.assertAlmostEqual(motif.counts["G", 0], 0.000000000)
        self.assertAlmostEqual(motif.counts["G", 1], 0.000000000)
        self.assertAlmostEqual(motif.counts["G", 2], 0.000000000)
        self.assertAlmostEqual(motif.counts["G", 3], 0.000000000)
        self.assertAlmostEqual(motif.counts["G", 4], 0.000000000)
        self.assertAlmostEqual(motif.counts["G", 5], 6.000000000)
        self.assertAlmostEqual(motif.counts["G", 6], 2.000000000)
        self.assertAlmostEqual(motif.counts["A", 0], 1.000000000)
        self.assertAlmostEqual(motif.counts["A", 1], 0.000000000)
        self.assertAlmostEqual(motif.counts["A", 2], 16.000000000)
        self.assertAlmostEqual(motif.counts["A", 3], 18.000000000)
        self.assertAlmostEqual(motif.counts["A", 4], 1.000000000)
        self.assertAlmostEqual(motif.counts["A", 5], 0.000000000)
        self.assertAlmostEqual(motif.counts["A", 6], 15.000000000)
        self.assertAlmostEqual(motif.counts["T", 0], 14.000000000)
        self.assertAlmostEqual(motif.counts["T", 1], 18.000000000)
        self.assertAlmostEqual(motif.counts["T", 2], 2.000000000)
        self.assertAlmostEqual(motif.counts["T", 3], 0.000000000)
        self.assertAlmostEqual(motif.counts["T", 4], 17.000000000)
        self.assertAlmostEqual(motif.counts["T", 5], 12.000000000)
        self.assertAlmostEqual(motif.counts["T", 6], 0.000000000)
        self.assertAlmostEqual(motif.counts["C", 0], 3.000000000)
        self.assertAlmostEqual(motif.counts["C", 1], 0.000000000)
        self.assertAlmostEqual(motif.counts["C", 2], 0.000000000)
        self.assertAlmostEqual(motif.counts["C", 3], 0.000000000)
        self.assertAlmostEqual(motif.counts["C", 4], 0.000000000)
        self.assertAlmostEqual(motif.counts["C", 5], 0.000000000)
        self.assertAlmostEqual(motif.counts["C", 6], 1.000000000)
        self.assertEqual(motif.consensus, "TTAATTA")
        self.assertEqual(motif.degenerate_consensus, "TTAATKA")
        self.assertTrue(
            np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        1.0555114658337947,
                        2.0,
                        1.4967416652243541,
                        2.0,
                        1.6904565708496748,
                        1.0817041659455104,
                        1.1969282726758976,
                    ]
                ),
            )
        )
        self.assertEqual(motif[1:-2].consensus, "TAAT")
        motif = record[4]
        self.assertEqual(
            motif.name,
            "ATGACTCATC AP-1(bZIP)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer    6.049537    -1.782996e+03   0   9805.3,5781.0,3085.1,2715.0,0.00e+00",
        )
        self.assertEqual(motif.length, 10)
        self.assertEqual(motif.alphabet, "GATC")
        self.assertAlmostEqual(motif.counts["G", 0], 0.277000000)
        self.assertAlmostEqual(motif.counts["G", 1], 0.001000000)
        self.assertAlmostEqual(motif.counts["G", 2], 0.965000000)
        self.assertAlmostEqual(motif.counts["G", 3], 0.001000000)
        self.assertAlmostEqual(motif.counts["G", 4], 0.305000000)
        self.assertAlmostEqual(motif.counts["G", 5], 0.001000000)
        self.assertAlmostEqual(motif.counts["G", 6], 0.001000000)
        self.assertAlmostEqual(motif.counts["G", 7], 0.001000000)
        self.assertAlmostEqual(motif.counts["G", 8], 0.307000000)
        self.assertAlmostEqual(motif.counts["G", 9], 0.211000000)
        self.assertAlmostEqual(motif.counts["A", 0], 0.419000000)
        self.assertAlmostEqual(motif.counts["A", 1], 0.001000000)
        self.assertAlmostEqual(motif.counts["A", 2], 0.010000000)
        self.assertAlmostEqual(motif.counts["A", 3], 0.984000000)
        self.assertAlmostEqual(motif.counts["A", 4], 0.062000000)
        self.assertAlmostEqual(motif.counts["A", 5], 0.026000000)
        self.assertAlmostEqual(motif.counts["A", 6], 0.043000000)
        self.assertAlmostEqual(motif.counts["A", 7], 0.980000000)
        self.assertAlmostEqual(motif.counts["A", 8], 0.050000000)
        self.assertAlmostEqual(motif.counts["A", 9], 0.149000000)
        self.assertAlmostEqual(motif.counts["T", 0], 0.028000000)
        self.assertAlmostEqual(motif.counts["T", 1], 0.997000000)
        self.assertAlmostEqual(motif.counts["T", 2], 0.023000000)
        self.assertAlmostEqual(motif.counts["T", 3], 0.012000000)
        self.assertAlmostEqual(motif.counts["T", 4], 0.054000000)
        self.assertAlmostEqual(motif.counts["T", 5], 0.972000000)
        self.assertAlmostEqual(motif.counts["T", 6], 0.012000000)
        self.assertAlmostEqual(motif.counts["T", 7], 0.014000000)
        self.assertAlmostEqual(motif.counts["T", 8], 0.471000000)
        self.assertAlmostEqual(motif.counts["T", 9], 0.195000000)
        self.assertAlmostEqual(motif.counts["C", 0], 0.275000000)
        self.assertAlmostEqual(motif.counts["C", 1], 0.001000000)
        self.assertAlmostEqual(motif.counts["C", 2], 0.002000000)
        self.assertAlmostEqual(motif.counts["C", 3], 0.003000000)
        self.assertAlmostEqual(motif.counts["C", 4], 0.579000000)
        self.assertAlmostEqual(motif.counts["C", 5], 0.001000000)
        self.assertAlmostEqual(motif.counts["C", 6], 0.943000000)
        self.assertAlmostEqual(motif.counts["C", 7], 0.005000000)
        self.assertAlmostEqual(motif.counts["C", 8], 0.172000000)
        self.assertAlmostEqual(motif.counts["C", 9], 0.444000000)
        self.assertEqual(motif.consensus, "ATGACTCATC")
        self.assertEqual(motif.degenerate_consensus, "NTGASTCAKN")
        self.assertTrue(
            np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        0.30427230622817475,
                        1.9657810606529142,
                        1.7408585738061,
                        1.8654244261025423,
                        0.5449286810918202,
                        1.8033449015144003,
                        1.639502374827662,
                        1.8370335049436752,
                        0.3124728907316759,
                        0.13671828556764112,
                    ]
                ),
            )
        )
        self.assertEqual(motif[1:-2].consensus, "TGACTCA")
        motif = record[5]
        self.assertEqual(motif.name, "AHR_si")
        self.assertEqual(motif.length, 9)
        self.assertEqual(motif.alphabet, "GATC")
        self.assertAlmostEqual(motif.counts["G", 0], 56.412537571)
        self.assertAlmostEqual(motif.counts["G", 1], 34.663129823)
        self.assertAlmostEqual(motif.counts["G", 2], 20.706746562)
        self.assertAlmostEqual(motif.counts["G", 3], 145.863705132)
        self.assertAlmostEqual(motif.counts["G", 4], 1.492783630)
        self.assertAlmostEqual(motif.counts["G", 5], 149.376137203)
        self.assertAlmostEqual(motif.counts["G", 6], 0.702486414)
        self.assertAlmostEqual(motif.counts["G", 7], 153.958717377)
        self.assertAlmostEqual(motif.counts["G", 8], 16.159862547)
        self.assertAlmostEqual(motif.counts["A", 0], 40.513432405)
        self.assertAlmostEqual(motif.counts["A", 1], 10.877470983)
        self.assertAlmostEqual(motif.counts["A", 2], 21.716570782)
        self.assertAlmostEqual(motif.counts["A", 3], 2.546513251)
        self.assertAlmostEqual(motif.counts["A", 4], 0.000000000)
        self.assertAlmostEqual(motif.counts["A", 5], 3.441039751)
        self.assertAlmostEqual(motif.counts["A", 6], 0.000000000)
        self.assertAlmostEqual(motif.counts["A", 7], 0.000000000)
        self.assertAlmostEqual(motif.counts["A", 8], 43.079223333)
        self.assertAlmostEqual(motif.counts["T", 0], 38.773634853)
        self.assertAlmostEqual(motif.counts["T", 1], 96.547239851)
        self.assertAlmostEqual(motif.counts["T", 2], 67.652320196)
        self.assertAlmostEqual(motif.counts["T", 3], 4.231336967)
        self.assertAlmostEqual(motif.counts["T", 4], 2.107459242)
        self.assertAlmostEqual(motif.counts["T", 5], 0.351243207)
        self.assertAlmostEqual(motif.counts["T", 6], 149.815191211)
        self.assertAlmostEqual(motif.counts["T", 7], 0.000000000)
        self.assertAlmostEqual(motif.counts["T", 8], 27.844049228)
        self.assertAlmostEqual(motif.counts["C", 0], 18.259112548)
        self.assertAlmostEqual(motif.counts["C", 1], 11.870876720)
        self.assertAlmostEqual(motif.counts["C", 2], 43.883079838)
        self.assertAlmostEqual(motif.counts["C", 3], 1.317162026)
        self.assertAlmostEqual(motif.counts["C", 4], 150.358474505)
        self.assertAlmostEqual(motif.counts["C", 5], 0.790297216)
        self.assertAlmostEqual(motif.counts["C", 6], 3.441039751)
        self.assertAlmostEqual(motif.counts["C", 7], 0.000000000)
        self.assertAlmostEqual(motif.counts["C", 8], 66.875582269)
        self.assertEqual(motif.consensus, "GTTGCGTGC")
        self.assertEqual(motif.degenerate_consensus, "NTNGCGTGN")
        self.assertTrue(
            np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        0.09662409645348236,
                        0.5383413903068038,
                        0.17471270188228985,
                        1.6270151623731723,
                        1.8170663607301638,
                        1.7760800937680195,
                        1.803660112630464,
                        2.0,
                        0.17577786614573548,
                    ]
                ),
            )
        )
        self.assertEqual(motif[1:-2].consensus, "TTGCGT")
        motif = record[6]
        self.assertIsNone(motif.name)
        self.assertEqual(motif.length, 8)
        self.assertEqual(motif.alphabet, "GATC")
        self.assertAlmostEqual(motif.counts["G", 0], 0.098612000)
        self.assertAlmostEqual(motif.counts["G", 1], 0.025056000)
        self.assertAlmostEqual(motif.counts["G", 2], 0.918728000)
        self.assertAlmostEqual(motif.counts["G", 3], 0.029759000)
        self.assertAlmostEqual(motif.counts["G", 4], 0.104968000)
        self.assertAlmostEqual(motif.counts["G", 5], 0.006667000)
        self.assertAlmostEqual(motif.counts["G", 6], 0.026928000)
        self.assertAlmostEqual(motif.counts["G", 7], 0.005737000)
        self.assertAlmostEqual(motif.counts["A", 0], 0.772949000)
        self.assertAlmostEqual(motif.counts["A", 1], 0.026652000)
        self.assertAlmostEqual(motif.counts["A", 2], 0.017663000)
        self.assertAlmostEqual(motif.counts["A", 3], 0.919596000)
        self.assertAlmostEqual(motif.counts["A", 4], 0.060312000)
        self.assertAlmostEqual(motif.counts["A", 5], 0.037406000)
        self.assertAlmostEqual(motif.counts["A", 6], 0.047316000)
        self.assertAlmostEqual(motif.counts["A", 7], 0.948639000)
        self.assertAlmostEqual(motif.counts["T", 0], 0.038860000)
        self.assertAlmostEqual(motif.counts["T", 1], 0.943639000)
        self.assertAlmostEqual(motif.counts["T", 2], 0.040264000)
        self.assertAlmostEqual(motif.counts["T", 3], 0.025231000)
        self.assertAlmostEqual(motif.counts["T", 4], 0.062462000)
        self.assertAlmostEqual(motif.counts["T", 5], 0.935284000)
        self.assertAlmostEqual(motif.counts["T", 6], 0.026732000)
        self.assertAlmostEqual(motif.counts["T", 7], 0.026128000)
        self.assertAlmostEqual(motif.counts["C", 0], 0.089579000)
        self.assertAlmostEqual(motif.counts["C", 1], 0.004653000)
        self.assertAlmostEqual(motif.counts["C", 2], 0.023344000)
        self.assertAlmostEqual(motif.counts["C", 3], 0.025414000)
        self.assertAlmostEqual(motif.counts["C", 4], 0.772259000)
        self.assertAlmostEqual(motif.counts["C", 5], 0.020643000)
        self.assertAlmostEqual(motif.counts["C", 6], 0.899024000)
        self.assertAlmostEqual(motif.counts["C", 7], 0.019497000)
        self.assertEqual(motif.consensus, "ATGACTCA")
        self.assertEqual(motif.degenerate_consensus, "ATGACTCA")
        self.assertTrue(
            np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        0.889358068874075,
                        1.6123293058245811,
                        1.471654165929799,
                        1.4693092198124151,
                        0.8764628815119266,
                        1.5686388858173408,
                        1.37357038822754,
                        1.6369796776980579,
                    ]
                ),
            )
        )
        self.assertEqual(motif[1:-2].consensus, "TGACT")
        motif = record[7]
        self.assertIsNone(motif.name)
        self.assertEqual(motif.length, 11)
        self.assertEqual(motif.alphabet, "GATC")
        self.assertAlmostEqual(motif.counts["G", 0], 28.0)
        self.assertAlmostEqual(motif.counts["G", 1], 0.0)
        self.assertAlmostEqual(motif.counts["G", 2], 14.0)
        self.assertAlmostEqual(motif.counts["G", 3], 0.0)
        self.assertAlmostEqual(motif.counts["G", 4], 0.0)
        self.assertAlmostEqual(motif.counts["G", 5], 7.0)
        self.assertAlmostEqual(motif.counts["G", 6], 11.0)
        self.assertAlmostEqual(motif.counts["G", 7], 38.0)
        self.assertAlmostEqual(motif.counts["G", 8], 0.0)
        self.assertAlmostEqual(motif.counts["G", 9], 25.0)
        self.assertAlmostEqual(motif.counts["G", 10], 0.0)
        self.assertAlmostEqual(motif.counts["A", 0], 0.0)
        self.assertAlmostEqual(motif.counts["A", 1], 0.0)
        self.assertAlmostEqual(motif.counts["A", 2], 55.0)
        self.assertAlmostEqual(motif.counts["A", 3], 99.0)
        self.assertAlmostEqual(motif.counts["A", 4], 78.0)
        self.assertAlmostEqual(motif.counts["A", 5], 52.0)
        self.assertAlmostEqual(motif.counts["A", 6], 46.0)
        self.assertAlmostEqual(motif.counts["A", 7], 60.0)
        self.assertAlmostEqual(motif.counts["A", 8], 33.0)
        self.assertAlmostEqual(motif.counts["A", 9], 0.0)
        self.assertAlmostEqual(motif.counts["A", 10], 0.0)
        self.assertAlmostEqual(motif.counts["T", 0], 30.0)
        self.assertAlmostEqual(motif.counts["T", 1], 0.0)
        self.assertAlmostEqual(motif.counts["T", 2], 0.0)
        self.assertAlmostEqual(motif.counts["T", 3], 0.0)
        self.assertAlmostEqual(motif.counts["T", 4], 20.0)
        self.assertAlmostEqual(motif.counts["T", 5], 0.0)
        self.assertAlmostEqual(motif.counts["T", 6], 19.0)
        self.assertAlmostEqual(motif.counts["T", 7], 0.0)
        self.assertAlmostEqual(motif.counts["T", 8], 0.0)
        self.assertAlmostEqual(motif.counts["T", 9], 73.0)
        self.assertAlmostEqual(motif.counts["T", 10], 99.0)
        self.assertAlmostEqual(motif.counts["C", 0], 40.0)
        self.assertAlmostEqual(motif.counts["C", 1], 99.0)
        self.assertAlmostEqual(motif.counts["C", 2], 29.0)
        self.assertAlmostEqual(motif.counts["C", 3], 0.0)
        self.assertAlmostEqual(motif.counts["C", 4], 0.0)
        self.assertAlmostEqual(motif.counts["C", 5], 39.0)
        self.assertAlmostEqual(motif.counts["C", 6], 22.0)
        self.assertAlmostEqual(motif.counts["C", 7], 0.0)
        self.assertAlmostEqual(motif.counts["C", 8], 66.0)
        self.assertAlmostEqual(motif.counts["C", 9], 0.0)
        self.assertAlmostEqual(motif.counts["C", 10], 0.0)
        self.assertEqual(motif.consensus, "CCAAAAAACTT")
        self.assertEqual(motif.degenerate_consensus, "BCMAAMNRMTT")
        self.assertTrue(
            np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        0.43314504855176084,
                        2.0,
                        0.6114044621231828,
                        2.0,
                        1.2699833698542062,
                        0.7139129756130338,
                        0.1909607288346033,
                        1.0366644543273158,
                        1.0817041659455104,
                        1.180735028768561,
                        2.0,
                    ]
                ),
            )
        )
        self.assertEqual(motif[1:-2].consensus, "CAAAAAAC")

    def test_pfm_four_rows_parsing(self):
        """Test if Bio.motifs.pfm can parse motifs in position frequency matrix format (4 rows)."""
        with open("motifs/fourrows.pfm") as stream:
            record = motifs.parse(stream, "pfm-four-rows")
        self.assertEqual(len(record), 8)
        motif = record[0]
        self.assertEqual(motif.name, "")
        self.assertEqual(motif.length, 6)
        self.assertEqual(motif.alphabet, "GATC")
        self.assertAlmostEqual(motif.counts["G", 0], 5.0)
        self.assertAlmostEqual(motif.counts["G", 1], 0.0)
        self.assertAlmostEqual(motif.counts["G", 2], 0.0)
        self.assertAlmostEqual(motif.counts["G", 3], 0.0)
        self.assertAlmostEqual(motif.counts["G", 4], 3.0)
        self.assertAlmostEqual(motif.counts["G", 5], 0.0)
        self.assertAlmostEqual(motif.counts["A", 0], 0.0)
        self.assertAlmostEqual(motif.counts["A", 1], 5.0)
        self.assertAlmostEqual(motif.counts["A", 2], 6.0)
        self.assertAlmostEqual(motif.counts["A", 3], 5.0)
        self.assertAlmostEqual(motif.counts["A", 4], 1.0)
        self.assertAlmostEqual(motif.counts["A", 5], 0.0)
        self.assertAlmostEqual(motif.counts["T", 0], 0.0)
        self.assertAlmostEqual(motif.counts["T", 1], 0.0)
        self.assertAlmostEqual(motif.counts["T", 2], 0.0)
        self.assertAlmostEqual(motif.counts["T", 3], 1.0)
        self.assertAlmostEqual(motif.counts["T", 4], 2.0)
        self.assertAlmostEqual(motif.counts["T", 5], 2.0)
        self.assertAlmostEqual(motif.counts["C", 0], 1.0)
        self.assertAlmostEqual(motif.counts["C", 1], 1.0)
        self.assertAlmostEqual(motif.counts["C", 2], 0.0)
        self.assertAlmostEqual(motif.counts["C", 3], 0.0)
        self.assertAlmostEqual(motif.counts["C", 4], 0.0)
        self.assertAlmostEqual(motif.counts["C", 5], 4.0)
        self.assertEqual(motif.consensus, "GAAAGC")
        self.assertEqual(motif.degenerate_consensus, "GAAAKY")
        self.assertTrue(
            np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        1.349977578351646,
                        1.349977578351646,
                        2.0,
                        1.349977578351646,
                        0.5408520829727552,
                        1.0817041659455104,
                    ]
                ),
            )
        )
        self.assertEqual(motif[:-2].consensus, "GAAA")
        motif = record[1]
        self.assertEqual(motif.name, "")
        self.assertEqual(motif.length, 15)
        self.assertEqual(motif.alphabet, "GATC")
        self.assertAlmostEqual(motif.counts["G", 0], 0.000000000)
        self.assertAlmostEqual(motif.counts["G", 1], 1.000000000)
        self.assertAlmostEqual(motif.counts["G", 2], 0.000000000)
        self.assertAlmostEqual(motif.counts["G", 3], 0.250000000)
        self.assertAlmostEqual(motif.counts["G", 4], 0.250000000)
        self.assertAlmostEqual(motif.counts["G", 5], 0.250000000)
        self.assertAlmostEqual(motif.counts["G", 6], 0.250000000)
        self.assertAlmostEqual(motif.counts["G", 7], 0.250000000)
        self.assertAlmostEqual(motif.counts["G", 8], 0.250000000)
        self.assertAlmostEqual(motif.counts["G", 9], 0.250000000)
        self.assertAlmostEqual(motif.counts["G", 10], 0.250000000)
        self.assertAlmostEqual(motif.counts["G", 11], 0.250000000)
        self.assertAlmostEqual(motif.counts["G", 12], 0.000000000)
        self.assertAlmostEqual(motif.counts["G", 13], 1.000000000)
        self.assertAlmostEqual(motif.counts["G", 14], 0.250000000)
        self.assertAlmostEqual(motif.counts["A", 0], 0.500000000)
        self.assertAlmostEqual(motif.counts["A", 1], 0.000000000)
        self.assertAlmostEqual(motif.counts["A", 2], 0.000000000)
        self.assertAlmostEqual(motif.counts["A", 3], 0.250000000)
        self.assertAlmostEqual(motif.counts["A", 4], 0.250000000)
        self.assertAlmostEqual(motif.counts["A", 5], 0.250000000)
        self.assertAlmostEqual(motif.counts["A", 6], 0.250000000)
        self.assertAlmostEqual(motif.counts["A", 7], 0.250000000)
        self.assertAlmostEqual(motif.counts["A", 8], 0.250000000)
        self.assertAlmostEqual(motif.counts["A", 9], 0.250000000)
        self.assertAlmostEqual(motif.counts["A", 10], 0.250000000)
        self.assertAlmostEqual(motif.counts["A", 11], 0.250000000)
        self.assertAlmostEqual(motif.counts["A", 12], 0.500000000)
        self.assertAlmostEqual(motif.counts["A", 13], 0.000000000)
        self.assertAlmostEqual(motif.counts["A", 14], 0.083333333)
        self.assertAlmostEqual(motif.counts["T", 0], 0.000000000)
        self.assertAlmostEqual(motif.counts["T", 1], 0.000000000)
        self.assertAlmostEqual(motif.counts["T", 2], 0.000000000)
        self.assertAlmostEqual(motif.counts["T", 3], 0.250000000)
        self.assertAlmostEqual(motif.counts["T", 4], 0.250000000)
        self.assertAlmostEqual(motif.counts["T", 5], 0.250000000)
        self.assertAlmostEqual(motif.counts["T", 6], 0.250000000)
        self.assertAlmostEqual(motif.counts["T", 7], 0.250000000)
        self.assertAlmostEqual(motif.counts["T", 8], 0.250000000)
        self.assertAlmostEqual(motif.counts["T", 9], 0.250000000)
        self.assertAlmostEqual(motif.counts["T", 10], 0.250000000)
        self.assertAlmostEqual(motif.counts["T", 11], 0.250000000)
        self.assertAlmostEqual(motif.counts["T", 12], 0.000000000)
        self.assertAlmostEqual(motif.counts["T", 13], 0.000000000)
        self.assertAlmostEqual(motif.counts["T", 14], 0.083333333)
        self.assertAlmostEqual(motif.counts["C", 0], 0.500000000)
        self.assertAlmostEqual(motif.counts["C", 1], 0.000000000)
        self.assertAlmostEqual(motif.counts["C", 2], 1.000000000)
        self.assertAlmostEqual(motif.counts["C", 3], 0.250000000)
        self.assertAlmostEqual(motif.counts["C", 4], 0.250000000)
        self.assertAlmostEqual(motif.counts["C", 5], 0.250000000)
        self.assertAlmostEqual(motif.counts["C", 6], 0.250000000)
        self.assertAlmostEqual(motif.counts["C", 7], 0.250000000)
        self.assertAlmostEqual(motif.counts["C", 8], 0.250000000)
        self.assertAlmostEqual(motif.counts["C", 9], 0.250000000)
        self.assertAlmostEqual(motif.counts["C", 10], 0.250000000)
        self.assertAlmostEqual(motif.counts["C", 11], 0.250000000)
        self.assertAlmostEqual(motif.counts["C", 12], 0.500000000)
        self.assertAlmostEqual(motif.counts["C", 13], 0.000000000)
        self.assertAlmostEqual(motif.counts["C", 14], 0.583333333)
        self.assertEqual(motif.consensus, "AGCGGGGGGGGGAGC")
        self.assertEqual(motif.degenerate_consensus, "MGCNNNNNNNNNMGC")
        self.assertTrue(
            np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        1.0,
                        2.0,
                        2.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        1.0,
                        2.0,
                        0.44890182844369547,
                    ]
                ),
            )
        )
        self.assertEqual(motif[:-2].consensus, "AGCGGGGGGGGGA")
        motif = record[2]
        self.assertEqual(motif.name, "")
        self.assertEqual(motif.length, 15)
        self.assertEqual(motif.alphabet, "GATC")
        self.assertAlmostEqual(motif.counts["G", 0], 270.0)
        self.assertAlmostEqual(motif.counts["G", 1], 398.0)
        self.assertAlmostEqual(motif.counts["G", 2], 54.0)
        self.assertAlmostEqual(motif.counts["G", 3], 164.0)
        self.assertAlmostEqual(motif.counts["G", 4], 7.0)
        self.assertAlmostEqual(motif.counts["G", 5], 659.0)
        self.assertAlmostEqual(motif.counts["G", 6], 1.0)
        self.assertAlmostEqual(motif.counts["G", 7], 750.0)
        self.assertAlmostEqual(motif.counts["G", 8], 755.0)
        self.assertAlmostEqual(motif.counts["G", 9], 65.0)
        self.assertAlmostEqual(motif.counts["G", 10], 1.0)
        self.assertAlmostEqual(motif.counts["G", 11], 41.0)
        self.assertAlmostEqual(motif.counts["G", 12], 202.0)
        self.assertAlmostEqual(motif.counts["G", 13], 234.0)
        self.assertAlmostEqual(motif.counts["G", 14], 205.0)
        self.assertAlmostEqual(motif.counts["A", 0], 92.0)
        self.assertAlmostEqual(motif.counts["A", 1], 106.0)
        self.assertAlmostEqual(motif.counts["A", 2], 231.0)
        self.assertAlmostEqual(motif.counts["A", 3], 135.0)
        self.assertAlmostEqual(motif.counts["A", 4], 0.0)
        self.assertAlmostEqual(motif.counts["A", 5], 1.0)
        self.assertAlmostEqual(motif.counts["A", 6], 780.0)
        self.assertAlmostEqual(motif.counts["A", 7], 28.0)
        self.assertAlmostEqual(motif.counts["A", 8], 0.0)
        self.assertAlmostEqual(motif.counts["A", 9], 700.0)
        self.assertAlmostEqual(motif.counts["A", 10], 739.0)
        self.assertAlmostEqual(motif.counts["A", 11], 94.0)
        self.assertAlmostEqual(motif.counts["A", 12], 60.0)
        self.assertAlmostEqual(motif.counts["A", 13], 127.0)
        self.assertAlmostEqual(motif.counts["A", 14], 130.0)
        self.assertAlmostEqual(motif.counts["T", 0], 290.0)
        self.assertAlmostEqual(motif.counts["T", 1], 204.0)
        self.assertAlmostEqual(motif.counts["T", 2], 375.0)
        self.assertAlmostEqual(motif.counts["T", 3], 411.0)
        self.assertAlmostEqual(motif.counts["T", 4], 9.0)
        self.assertAlmostEqual(motif.counts["T", 5], 127.0)
        self.assertAlmostEqual(motif.counts["T", 6], 6.0)
        self.assertAlmostEqual(motif.counts["T", 7], 11.0)
        self.assertAlmostEqual(motif.counts["T", 8], 36.0)
        self.assertAlmostEqual(motif.counts["T", 9], 20.0)
        self.assertAlmostEqual(motif.counts["T", 10], 31.0)
        self.assertAlmostEqual(motif.counts["T", 11], 605.0)
        self.assertAlmostEqual(motif.counts["T", 12], 335.0)
        self.assertAlmostEqual(motif.counts["T", 13], 307.0)
        self.assertAlmostEqual(motif.counts["T", 14], 308.0)
        self.assertAlmostEqual(motif.counts["C", 0], 138.0)
        self.assertAlmostEqual(motif.counts["C", 1], 82.0)
        self.assertAlmostEqual(motif.counts["C", 2], 129.0)
        self.assertAlmostEqual(motif.counts["C", 3], 81.0)
        self.assertAlmostEqual(motif.counts["C", 4], 774.0)
        self.assertAlmostEqual(motif.counts["C", 5], 1.0)
        self.assertAlmostEqual(motif.counts["C", 6], 3.0)
        self.assertAlmostEqual(motif.counts["C", 7], 1.0)
        self.assertAlmostEqual(motif.counts["C", 8], 0.0)
        self.assertAlmostEqual(motif.counts["C", 9], 6.0)
        self.assertAlmostEqual(motif.counts["C", 10], 17.0)
        self.assertAlmostEqual(motif.counts["C", 11], 49.0)
        self.assertAlmostEqual(motif.counts["C", 12], 193.0)
        self.assertAlmostEqual(motif.counts["C", 13], 122.0)
        self.assertAlmostEqual(motif.counts["C", 14], 148.0)
        self.assertEqual(motif.consensus, "TGTTCGAGGAATTTT")
        self.assertEqual(motif.degenerate_consensus, "NKWTCGAGGAATNNN")
        self.assertTrue(
            np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        0.13892143832881046,
                        0.2692660952911542,
                        0.27915566353819243,
                        0.2665840150038887,
                        1.8371160692433293,
                        1.3354706334248059,
                        1.8856611660889357,
                        1.6600123906824402,
                        1.7329826640509962,
                        1.3601399752384014,
                        1.5978925123167893,
                        0.8698961051280728,
                        0.19290147849975406,
                        0.11003972948477392,
                        0.08469189143040626,
                    ]
                ),
            )
        )
        self.assertEqual(motif[:-2].consensus, "TGTTCGAGGAATT")
        motif = record[3]
        self.assertEqual(motif.name, "")
        self.assertEqual(motif.length, 6)
        self.assertEqual(motif.alphabet, "GATC")
        self.assertAlmostEqual(motif.counts["G", 0], 2.0)
        self.assertAlmostEqual(motif.counts["G", 1], 1.0)
        self.assertAlmostEqual(motif.counts["G", 2], 1.0)
        self.assertAlmostEqual(motif.counts["G", 3], 1.0)
        self.assertAlmostEqual(motif.counts["G", 4], 97.0)
        self.assertAlmostEqual(motif.counts["G", 5], 2.0)
        self.assertAlmostEqual(motif.counts["A", 0], 9.0)
        self.assertAlmostEqual(motif.counts["A", 1], 1.0)
        self.assertAlmostEqual(motif.counts["A", 2], 1.0)
        self.assertAlmostEqual(motif.counts["A", 3], 97.0)
        self.assertAlmostEqual(motif.counts["A", 4], 1.0)
        self.assertAlmostEqual(motif.counts["A", 5], 94.0)
        self.assertAlmostEqual(motif.counts["T", 0], 80.0)
        self.assertAlmostEqual(motif.counts["T", 1], 1.0)
        self.assertAlmostEqual(motif.counts["T", 2], 97.0)
        self.assertAlmostEqual(motif.counts["T", 3], 1.0)
        self.assertAlmostEqual(motif.counts["T", 4], 1.0)
        self.assertAlmostEqual(motif.counts["T", 5], 2.0)
        self.assertAlmostEqual(motif.counts["C", 0], 9.0)
        self.assertAlmostEqual(motif.counts["C", 1], 97.0)
        self.assertAlmostEqual(motif.counts["C", 2], 1.0)
        self.assertAlmostEqual(motif.counts["C", 3], 1.0)
        self.assertAlmostEqual(motif.counts["C", 4], 1.0)
        self.assertAlmostEqual(motif.counts["C", 5], 2.0)
        self.assertEqual(motif.consensus, "TCTAGA")
        self.assertEqual(motif.degenerate_consensus, "TCTAGA")
        self.assertTrue(
            np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        1.0042727863947818,
                        1.758059267146789,
                        1.7580592671467892,
                        1.7580592671467892,
                        1.7580592671467892,
                        1.5774573308022544,
                    ]
                ),
            )
        )
        self.assertEqual(motif[:-2].consensus, "TCTA")
        motif = record[4]
        self.assertEqual(motif.name, "")
        self.assertEqual(motif.length, 6)
        self.assertEqual(motif.alphabet, "GATC")
        self.assertAlmostEqual(motif.counts["G", 0], 0.02)
        self.assertAlmostEqual(motif.counts["G", 1], 0.01)
        self.assertAlmostEqual(motif.counts["G", 2], 0.01)
        self.assertAlmostEqual(motif.counts["G", 3], 0.01)
        self.assertAlmostEqual(motif.counts["G", 4], 0.97)
        self.assertAlmostEqual(motif.counts["G", 5], 0.02)
        self.assertAlmostEqual(motif.counts["A", 0], 0.09)
        self.assertAlmostEqual(motif.counts["A", 1], 0.01)
        self.assertAlmostEqual(motif.counts["A", 2], 0.01)
        self.assertAlmostEqual(motif.counts["A", 3], 0.97)
        self.assertAlmostEqual(motif.counts["A", 4], 0.01)
        self.assertAlmostEqual(motif.counts["A", 5], 0.94)
        self.assertAlmostEqual(motif.counts["T", 0], 0.80)
        self.assertAlmostEqual(motif.counts["T", 1], 0.01)
        self.assertAlmostEqual(motif.counts["T", 2], 0.97)
        self.assertAlmostEqual(motif.counts["T", 3], 0.01)
        self.assertAlmostEqual(motif.counts["T", 4], 0.01)
        self.assertAlmostEqual(motif.counts["T", 5], 0.02)
        self.assertAlmostEqual(motif.counts["C", 0], 0.09)
        self.assertAlmostEqual(motif.counts["C", 1], 0.97)
        self.assertAlmostEqual(motif.counts["C", 2], 0.01)
        self.assertAlmostEqual(motif.counts["C", 3], 0.01)
        self.assertAlmostEqual(motif.counts["C", 4], 0.01)
        self.assertAlmostEqual(motif.counts["C", 5], 0.02)
        self.assertEqual(motif.consensus, "TCTAGA")
        self.assertEqual(motif.degenerate_consensus, "TCTAGA")
        self.assertTrue(
            np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        1.0042727863947818,
                        1.758059267146789,
                        1.7580592671467892,
                        1.7580592671467892,
                        1.7580592671467892,
                        1.5774573308022544,
                    ]
                ),
            )
        )
        self.assertEqual(motif[:-2].consensus, "TCTA")
        motif = record[5]
        self.assertEqual(motif.name, "abd-A")
        self.assertEqual(motif.length, 8)
        self.assertEqual(motif.alphabet, "GATC")
        self.assertAlmostEqual(motif.counts["G", 0], 0.455991516)
        self.assertAlmostEqual(motif.counts["G", 1], 0.069194062)
        self.assertAlmostEqual(motif.counts["G", 2], 0.010869565)
        self.assertAlmostEqual(motif.counts["G", 3], 0.021739130)
        self.assertAlmostEqual(motif.counts["G", 4], 0.028499470)
        self.assertAlmostEqual(motif.counts["G", 5], 0.028499470)
        self.assertAlmostEqual(motif.counts["G", 6], 0.016304348)
        self.assertAlmostEqual(motif.counts["G", 7], 0.160127253)
        self.assertAlmostEqual(motif.counts["A", 0], 0.218451750)
        self.assertAlmostEqual(motif.counts["A", 1], 0.023064687)
        self.assertAlmostEqual(motif.counts["A", 2], 0.656680806)
        self.assertAlmostEqual(motif.counts["A", 3], 0.898197243)
        self.assertAlmostEqual(motif.counts["A", 4], 0.040694592)
        self.assertAlmostEqual(motif.counts["A", 5], 0.132953340)
        self.assertAlmostEqual(motif.counts["A", 6], 0.749072110)
        self.assertAlmostEqual(motif.counts["A", 7], 0.628313892)
        self.assertAlmostEqual(motif.counts["T", 0], 0.235949099)
        self.assertAlmostEqual(motif.counts["T", 1], 0.590402969)
        self.assertAlmostEqual(motif.counts["T", 2], 0.010869565)
        self.assertAlmostEqual(motif.counts["T", 3], 0.033934252)
        self.assertAlmostEqual(motif.counts["T", 4], 0.880567338)
        self.assertAlmostEqual(motif.counts["T", 5], 0.797852598)
        self.assertAlmostEqual(motif.counts["T", 6], 0.206124072)
        self.assertAlmostEqual(motif.counts["T", 7], 0.177624602)
        self.assertAlmostEqual(motif.counts["C", 0], 0.089607635)
        self.assertAlmostEqual(motif.counts["C", 1], 0.317338282)
        self.assertAlmostEqual(motif.counts["C", 2], 0.321580064)
        self.assertAlmostEqual(motif.counts["C", 3], 0.046129374)
        self.assertAlmostEqual(motif.counts["C", 4], 0.050238600)
        self.assertAlmostEqual(motif.counts["C", 5], 0.040694592)
        self.assertAlmostEqual(motif.counts["C", 6], 0.028499470)
        self.assertAlmostEqual(motif.counts["C", 7], 0.033934252)
        self.assertEqual(motif.consensus, "GTAATTAA")
        self.assertEqual(motif.degenerate_consensus, "NYAATTAA")
        self.assertTrue(
            np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        0.2005361303021225,
                        0.6336277209668335,
                        0.933405467206956,
                        1.3704286046679186,
                        1.2873833086962072,
                        1.0187720746919493,
                        0.975022432438911,
                        0.547109562258496,
                    ]
                ),
            )
        )
        self.assertEqual(motif[:-2].consensus, "GTAATT")
        motif = record[6]
        self.assertEqual(motif.name, "MA0001.1 AGL3")
        self.assertEqual(motif.length, 10)
        self.assertEqual(motif.alphabet, "GATC")
        self.assertAlmostEqual(motif.counts["G", 0], 1.0)
        self.assertAlmostEqual(motif.counts["G", 1], 0.0)
        self.assertAlmostEqual(motif.counts["G", 2], 3.0)
        self.assertAlmostEqual(motif.counts["G", 3], 4.0)
        self.assertAlmostEqual(motif.counts["G", 4], 1.0)
        self.assertAlmostEqual(motif.counts["G", 5], 0.0)
        self.assertAlmostEqual(motif.counts["G", 6], 5.0)
        self.assertAlmostEqual(motif.counts["G", 7], 3.0)
        self.assertAlmostEqual(motif.counts["G", 8], 28.0)
        self.assertAlmostEqual(motif.counts["G", 9], 88.0)
        self.assertAlmostEqual(motif.counts["A", 0], 0.0)
        self.assertAlmostEqual(motif.counts["A", 1], 3.0)
        self.assertAlmostEqual(motif.counts["A", 2], 79.0)
        self.assertAlmostEqual(motif.counts["A", 3], 40.0)
        self.assertAlmostEqual(motif.counts["A", 4], 66.0)
        self.assertAlmostEqual(motif.counts["A", 5], 48.0)
        self.assertAlmostEqual(motif.counts["A", 6], 65.0)
        self.assertAlmostEqual(motif.counts["A", 7], 11.0)
        self.assertAlmostEqual(motif.counts["A", 8], 65.0)
        self.assertAlmostEqual(motif.counts["A", 9], 0.0)
        self.assertAlmostEqual(motif.counts["T", 0], 2.0)
        self.assertAlmostEqual(motif.counts["T", 1], 19.0)
        self.assertAlmostEqual(motif.counts["T", 2], 11.0)
        self.assertAlmostEqual(motif.counts["T", 3], 50.0)
        self.assertAlmostEqual(motif.counts["T", 4], 29.0)
        self.assertAlmostEqual(motif.counts["T", 5], 47.0)
        self.assertAlmostEqual(motif.counts["T", 6], 22.0)
        self.assertAlmostEqual(motif.counts["T", 7], 81.0)
        self.assertAlmostEqual(motif.counts["T", 8], 1.0)
        self.assertAlmostEqual(motif.counts["T", 9], 6.0)
        self.assertAlmostEqual(motif.counts["C", 0], 94.0)
        self.assertAlmostEqual(motif.counts["C", 1], 75.0)
        self.assertAlmostEqual(motif.counts["C", 2], 4.0)
        self.assertAlmostEqual(motif.counts["C", 3], 3.0)
        self.assertAlmostEqual(motif.counts["C", 4], 1.0)
        self.assertAlmostEqual(motif.counts["C", 5], 2.0)
        self.assertAlmostEqual(motif.counts["C", 6], 5.0)
        self.assertAlmostEqual(motif.counts["C", 7], 2.0)
        self.assertAlmostEqual(motif.counts["C", 8], 3.0)
        self.assertAlmostEqual(motif.counts["C", 9], 3.0)
        self.assertEqual(motif.consensus, "CCATAAATAG")
        self.assertEqual(motif.degenerate_consensus, "CCAWAWATAG")
        self.assertTrue(
            np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        1.7725753233561499,
                        1.0972718180683638,
                        1.0578945228970464,
                        0.6353945886004412,
                        0.9651537633423314,
                        0.8757972203228152,
                        0.6864859661195083,
                        1.1561334005018244,
                        0.8724039945822116,
                        1.4691041160249607,
                    ]
                ),
            )
        )
        self.assertEqual(motif[:-2].consensus, "CCATAAAT")
        motif = record[7]
        self.assertEqual(motif.name, "MA0001.1 AGL3")
        self.assertEqual(motif.length, 10)
        self.assertEqual(motif.alphabet, "GATC")
        self.assertAlmostEqual(motif.counts["G", 0], 1.0)
        self.assertAlmostEqual(motif.counts["G", 1], 0.0)
        self.assertAlmostEqual(motif.counts["G", 2], 3.0)
        self.assertAlmostEqual(motif.counts["G", 3], 4.0)
        self.assertAlmostEqual(motif.counts["G", 4], 1.0)
        self.assertAlmostEqual(motif.counts["G", 5], 0.0)
        self.assertAlmostEqual(motif.counts["G", 6], 5.0)
        self.assertAlmostEqual(motif.counts["G", 7], 3.0)
        self.assertAlmostEqual(motif.counts["G", 8], 28.0)
        self.assertAlmostEqual(motif.counts["G", 9], 88.0)
        self.assertAlmostEqual(motif.counts["A", 0], 0.0)
        self.assertAlmostEqual(motif.counts["A", 1], 3.0)
        self.assertAlmostEqual(motif.counts["A", 2], 79.0)
        self.assertAlmostEqual(motif.counts["A", 3], 40.0)
        self.assertAlmostEqual(motif.counts["A", 4], 66.0)
        self.assertAlmostEqual(motif.counts["A", 5], 48.0)
        self.assertAlmostEqual(motif.counts["A", 6], 65.0)
        self.assertAlmostEqual(motif.counts["A", 7], 11.0)
        self.assertAlmostEqual(motif.counts["A", 8], 65.0)
        self.assertAlmostEqual(motif.counts["A", 9], 0.0)
        self.assertAlmostEqual(motif.counts["T", 0], 2.0)
        self.assertAlmostEqual(motif.counts["T", 1], 19.0)
        self.assertAlmostEqual(motif.counts["T", 2], 11.0)
        self.assertAlmostEqual(motif.counts["T", 3], 50.0)
        self.assertAlmostEqual(motif.counts["T", 4], 29.0)
        self.assertAlmostEqual(motif.counts["T", 5], 47.0)
        self.assertAlmostEqual(motif.counts["T", 6], 22.0)
        self.assertAlmostEqual(motif.counts["T", 7], 81.0)
        self.assertAlmostEqual(motif.counts["T", 8], 1.0)
        self.assertAlmostEqual(motif.counts["T", 9], 6.0)
        self.assertAlmostEqual(motif.counts["C", 0], 94.0)
        self.assertAlmostEqual(motif.counts["C", 1], 75.0)
        self.assertAlmostEqual(motif.counts["C", 2], 4.0)
        self.assertAlmostEqual(motif.counts["C", 3], 3.0)
        self.assertAlmostEqual(motif.counts["C", 4], 1.0)
        self.assertAlmostEqual(motif.counts["C", 5], 2.0)
        self.assertAlmostEqual(motif.counts["C", 6], 5.0)
        self.assertAlmostEqual(motif.counts["C", 7], 2.0)
        self.assertAlmostEqual(motif.counts["C", 8], 3.0)
        self.assertAlmostEqual(motif.counts["C", 9], 3.0)
        self.assertEqual(motif.consensus, "CCATAAATAG")
        self.assertEqual(motif.degenerate_consensus, "CCAWAWATAG")
        self.assertTrue(
            np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        1.7725753233561499,
                        1.0972718180683638,
                        1.0578945228970464,
                        0.6353945886004412,
                        0.9651537633423314,
                        0.8757972203228152,
                        0.6864859661195083,
                        1.1561334005018244,
                        0.8724039945822116,
                        1.4691041160249607,
                    ]
                ),
            )
        )
        self.assertEqual(motif[:-2].consensus, "CCATAAAT")

    def test_sites_parsing(self):
        """Test if Bio.motifs can parse JASPAR-style sites files."""
        with open("motifs/Arnt.sites") as stream:
            m = motifs.read(stream, "sites")
        self.assertEqual(m.length, 6)
        self.assertEqual(m.alignment.sequences[0], "CACGTG")
        self.assertEqual(m.alignment.sequences[1], "CACGTG")
        self.assertEqual(m.alignment.sequences[2], "CACGTG")
        self.assertEqual(m.alignment.sequences[3], "CACGTG")
        self.assertEqual(m.alignment.sequences[4], "CACGTG")
        self.assertEqual(m.alignment.sequences[5], "CACGTG")
        self.assertEqual(m.alignment.sequences[6], "CACGTG")
        self.assertEqual(m.alignment.sequences[7], "CACGTG")
        self.assertEqual(m.alignment.sequences[8], "CACGTG")
        self.assertEqual(m.alignment.sequences[9], "CACGTG")
        self.assertEqual(m.alignment.sequences[10], "CACGTG")
        self.assertEqual(m.alignment.sequences[11], "CACGTG")
        self.assertEqual(m.alignment.sequences[12], "CACGTG")
        self.assertEqual(m.alignment.sequences[13], "CACGTG")
        self.assertEqual(m.alignment.sequences[14], "CACGTG")
        self.assertEqual(m.alignment.sequences[15], "AACGTG")
        self.assertEqual(m.alignment.sequences[16], "AACGTG")
        self.assertEqual(m.alignment.sequences[17], "AACGTG")
        self.assertEqual(m.alignment.sequences[18], "AACGTG")
        self.assertEqual(m.alignment.sequences[19], "CGCGTG")
        self.assertAlmostEqual(m.counts["A", 0], 4)
        self.assertAlmostEqual(m.counts["A", 1], 19)
        self.assertAlmostEqual(m.counts["A", 2], 0)
        self.assertAlmostEqual(m.counts["A", 3], 0)
        self.assertAlmostEqual(m.counts["A", 4], 0)
        self.assertAlmostEqual(m.counts["A", 5], 0)
        self.assertAlmostEqual(m.counts["C", 0], 16)
        self.assertAlmostEqual(m.counts["C", 1], 0)
        self.assertAlmostEqual(m.counts["C", 2], 20)
        self.assertAlmostEqual(m.counts["C", 3], 0)
        self.assertAlmostEqual(m.counts["C", 4], 0)
        self.assertAlmostEqual(m.counts["C", 5], 0)
        self.assertAlmostEqual(m.counts["G", 0], 0)
        self.assertAlmostEqual(m.counts["G", 1], 1)
        self.assertAlmostEqual(m.counts["G", 2], 0)
        self.assertAlmostEqual(m.counts["G", 3], 20)
        self.assertAlmostEqual(m.counts["G", 4], 0)
        self.assertAlmostEqual(m.counts["G", 5], 20)
        self.assertAlmostEqual(m.counts["T", 0], 0)
        self.assertAlmostEqual(m.counts["T", 1], 0)
        self.assertAlmostEqual(m.counts["T", 2], 0)
        self.assertAlmostEqual(m.counts["T", 3], 0)
        self.assertAlmostEqual(m.counts["T", 4], 20)
        self.assertAlmostEqual(m.counts["T", 5], 0)
        self.assertEqual(m.consensus, "CACGTG")
        self.assertEqual(m.degenerate_consensus, "CACGTG")
        self.assertTrue(
            np.allclose(
                m.relative_entropy,
                np.array([1.278071905112638, 1.7136030428840439, 2.0, 2.0, 2.0, 2.0]),
            )
        )
        self.assertEqual(m[::2].consensus, "CCT")


class TestMEME(unittest.TestCase):
    def test_meme_parser_1(self):
        """Parse motifs/meme.INO_up800.classic.oops.xml file."""
        with open("motifs/meme.INO_up800.classic.oops.xml") as stream:
            record = motifs.parse(stream, "meme")
        self.assertEqual(record.version, "5.0.1")
        self.assertEqual(record.datafile, "common/INO_up800.s")
        self.assertEqual(record.alphabet, "ACGT")
        self.assertEqual(len(record.sequences), 7)
        self.assertEqual(record.sequences[0], "sequence_0")
        self.assertEqual(record.sequences[1], "sequence_1")
        self.assertEqual(record.sequences[2], "sequence_2")
        self.assertEqual(record.sequences[3], "sequence_3")
        self.assertEqual(record.sequences[4], "sequence_4")
        self.assertEqual(record.sequences[5], "sequence_5")
        self.assertEqual(record.sequences[6], "sequence_6")
        self.assertEqual(
            record.command,
            "meme common/INO_up800.s -oc results/meme10 -mod oops -dna -revcomp -bfile common/yeast.nc.6.freq -nmotifs 2 -objfun classic -minw 8 -nostatus ",
        )
        self.assertEqual(len(record), 2)
        motif = record[0]
        self.assertEqual(motif.name, "GSKGCATGTGAAA")
        self.assertEqual(record["GSKGCATGTGAAA"], motif)
        self.assertEqual(motif.num_occurrences, 7)
        self.assertAlmostEqual(motif.evalue, 0.19)
        self.assertEqual(motif.alphabet, "ACGT")
        # using the old instances property:
        with self.assertWarns(BiopythonDeprecationWarning):
            self.assertEqual(len(motif.instances), 7)
            self.assertAlmostEqual(motif.instances[0].pvalue, 1.21e-08)
            self.assertAlmostEqual(motif.instances[1].pvalue, 1.87e-08)
            self.assertAlmostEqual(motif.instances[2].pvalue, 6.62e-08)
            self.assertAlmostEqual(motif.instances[3].pvalue, 1.05e-07)
            self.assertAlmostEqual(motif.instances[4].pvalue, 1.69e-07)
            self.assertAlmostEqual(motif.instances[5].pvalue, 5.62e-07)
            self.assertAlmostEqual(motif.instances[6].pvalue, 1.08e-06)
            self.assertEqual(motif.instances[0].sequence_name, "INO1")
            self.assertEqual(motif.instances[1].sequence_name, "FAS1")
            self.assertEqual(motif.instances[2].sequence_name, "ACC1")
            self.assertEqual(motif.instances[3].sequence_name, "CHO2")
            self.assertEqual(motif.instances[4].sequence_name, "CHO1")
            self.assertEqual(motif.instances[5].sequence_name, "FAS2")
            self.assertEqual(motif.instances[6].sequence_name, "OPI3")
            self.assertEqual(motif.instances[0].sequence_id, "sequence_5")
            self.assertEqual(motif.instances[1].sequence_id, "sequence_2")
            self.assertEqual(motif.instances[2].sequence_id, "sequence_4")
            self.assertEqual(motif.instances[3].sequence_id, "sequence_1")
            self.assertEqual(motif.instances[4].sequence_id, "sequence_0")
            self.assertEqual(motif.instances[5].sequence_id, "sequence_3")
            self.assertEqual(motif.instances[6].sequence_id, "sequence_6")
            self.assertEqual(motif.instances[0].strand, "+")
            self.assertEqual(motif.instances[1].strand, "-")
            self.assertEqual(motif.instances[2].strand, "-")
            self.assertEqual(motif.instances[3].strand, "-")
            self.assertEqual(motif.instances[4].strand, "-")
            self.assertEqual(motif.instances[5].strand, "-")
            self.assertEqual(motif.instances[6].strand, "+")
            self.assertEqual(motif.instances[0].length, 13)
            self.assertEqual(motif.instances[1].length, 13)
            self.assertEqual(motif.instances[2].length, 13)
            self.assertEqual(motif.instances[3].length, 13)
            self.assertEqual(motif.instances[4].length, 13)
            self.assertEqual(motif.instances[5].length, 13)
            self.assertEqual(motif.instances[6].length, 13)
            self.assertEqual(motif.instances[0].start, 620)
            self.assertEqual(motif.instances[1].start, 94)
            self.assertEqual(motif.instances[2].start, 82)
            self.assertEqual(motif.instances[3].start, 353)
            self.assertEqual(motif.instances[4].start, 639)
            self.assertEqual(motif.instances[5].start, 566)
            self.assertEqual(motif.instances[6].start, 585)
            self.assertEqual(motif.instances[0], "GCGGCATGTGAAA")
            self.assertEqual(motif.instances[1], "GCGGCATGTGAAG")
            self.assertEqual(motif.instances[2], "GGGCCATGTGAAG")
            self.assertEqual(motif.instances[3], "GCGGCATGAGAAA")
            self.assertEqual(motif.instances[4], "GGTCCATGTGAAA")
            self.assertEqual(motif.instances[5], "GTAGCATGTGAAA")
            self.assertEqual(motif.instances[6], "AGTGCATGTGGAA")
        self.assertEqual(len(motif.alignment.sequences), 7)
        self.assertAlmostEqual(motif.alignment.sequences[0].pvalue, 1.21e-08)
        self.assertAlmostEqual(motif.alignment.sequences[1].pvalue, 1.87e-08)
        self.assertAlmostEqual(motif.alignment.sequences[2].pvalue, 6.62e-08)
        self.assertAlmostEqual(motif.alignment.sequences[3].pvalue, 1.05e-07)
        self.assertAlmostEqual(motif.alignment.sequences[4].pvalue, 1.69e-07)
        self.assertAlmostEqual(motif.alignment.sequences[5].pvalue, 5.62e-07)
        self.assertAlmostEqual(motif.alignment.sequences[6].pvalue, 1.08e-06)
        self.assertEqual(motif.alignment.sequences[0].sequence_name, "INO1")
        self.assertEqual(motif.alignment.sequences[1].sequence_name, "FAS1")
        self.assertEqual(motif.alignment.sequences[2].sequence_name, "ACC1")
        self.assertEqual(motif.alignment.sequences[3].sequence_name, "CHO2")
        self.assertEqual(motif.alignment.sequences[4].sequence_name, "CHO1")
        self.assertEqual(motif.alignment.sequences[5].sequence_name, "FAS2")
        self.assertEqual(motif.alignment.sequences[6].sequence_name, "OPI3")
        self.assertEqual(motif.alignment.sequences[0].sequence_id, "sequence_5")
        self.assertEqual(motif.alignment.sequences[1].sequence_id, "sequence_2")
        self.assertEqual(motif.alignment.sequences[2].sequence_id, "sequence_4")
        self.assertEqual(motif.alignment.sequences[3].sequence_id, "sequence_1")
        self.assertEqual(motif.alignment.sequences[4].sequence_id, "sequence_0")
        self.assertEqual(motif.alignment.sequences[5].sequence_id, "sequence_3")
        self.assertEqual(motif.alignment.sequences[6].sequence_id, "sequence_6")
        self.assertEqual(motif.alignment.sequences[0].strand, "+")
        self.assertEqual(motif.alignment.sequences[1].strand, "-")
        self.assertEqual(motif.alignment.sequences[2].strand, "-")
        self.assertEqual(motif.alignment.sequences[3].strand, "-")
        self.assertEqual(motif.alignment.sequences[4].strand, "-")
        self.assertEqual(motif.alignment.sequences[5].strand, "-")
        self.assertEqual(motif.alignment.sequences[6].strand, "+")
        self.assertEqual(motif.alignment.sequences[0].length, 13)
        self.assertEqual(motif.alignment.sequences[1].length, 13)
        self.assertEqual(motif.alignment.sequences[2].length, 13)
        self.assertEqual(motif.alignment.sequences[3].length, 13)
        self.assertEqual(motif.alignment.sequences[4].length, 13)
        self.assertEqual(motif.alignment.sequences[5].length, 13)
        self.assertEqual(motif.alignment.sequences[6].length, 13)
        self.assertEqual(motif.alignment.sequences[0].start, 620)
        self.assertEqual(motif.alignment.sequences[1].start, 94)
        self.assertEqual(motif.alignment.sequences[2].start, 82)
        self.assertEqual(motif.alignment.sequences[3].start, 353)
        self.assertEqual(motif.alignment.sequences[4].start, 639)
        self.assertEqual(motif.alignment.sequences[5].start, 566)
        self.assertEqual(motif.alignment.sequences[6].start, 585)
        self.assertEqual(motif.alignment.sequences[0], "GCGGCATGTGAAA")
        self.assertEqual(motif.alignment.sequences[1], "GCGGCATGTGAAG")
        self.assertEqual(motif.alignment.sequences[2], "GGGCCATGTGAAG")
        self.assertEqual(motif.alignment.sequences[3], "GCGGCATGAGAAA")
        self.assertEqual(motif.alignment.sequences[4], "GGTCCATGTGAAA")
        self.assertEqual(motif.alignment.sequences[5], "GTAGCATGTGAAA")
        self.assertEqual(motif.alignment.sequences[6], "AGTGCATGTGGAA")
        self.assertEqual(motif.consensus, "GCGGCATGTGAAA")
        self.assertEqual(motif.degenerate_consensus, "GSKGCATGTGAAA")
        self.assertTrue(
            np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        1.4083272214176723,
                        0.5511843642748154,
                        0.6212165065138244,
                        1.136879431433369,
                        2.0,
                        2.0,
                        2.0,
                        2.0,
                        1.4083272214176723,
                        2.0,
                        1.4083272214176723,
                        2.0,
                        1.136879431433369,
                    ]
                ),
            )
        )
        self.assertEqual(motif[1::2].consensus, "CGAGGA")
        motif = record[1]
        self.assertEqual(motif.name, "TTGACWCYTGCYCWG")
        self.assertEqual(record["TTGACWCYTGCYCWG"], motif)
        self.assertEqual(motif.num_occurrences, 7)
        self.assertAlmostEqual(motif.evalue, 54)
        self.assertEqual(motif.alphabet, "ACGT")
        # using the old instances property:
        with self.assertWarns(BiopythonDeprecationWarning):
            self.assertEqual(len(motif.instances), 7)
            self.assertAlmostEqual(motif.instances[0].pvalue, 7.2e-10)
            self.assertAlmostEqual(motif.instances[1].pvalue, 2.56e-08)
            self.assertAlmostEqual(motif.instances[2].pvalue, 1.59e-07)
            self.assertAlmostEqual(motif.instances[3].pvalue, 2.05e-07)
            self.assertAlmostEqual(motif.instances[4].pvalue, 3.85e-07)
            self.assertAlmostEqual(motif.instances[5].pvalue, 5.11e-07)
            self.assertAlmostEqual(motif.instances[6].pvalue, 8.01e-07)
            self.assertEqual(motif.instances[0].sequence_id, "sequence_1")
            self.assertEqual(motif.instances[1].sequence_id, "sequence_6")
            self.assertEqual(motif.instances[2].sequence_id, "sequence_4")
            self.assertEqual(motif.instances[3].sequence_id, "sequence_0")
            self.assertEqual(motif.instances[4].sequence_id, "sequence_2")
            self.assertEqual(motif.instances[5].sequence_id, "sequence_3")
            self.assertEqual(motif.instances[6].sequence_id, "sequence_5")
            self.assertEqual(motif.instances[0].strand, "+")
            self.assertEqual(motif.instances[1].strand, "-")
            self.assertEqual(motif.instances[2].strand, "-")
            self.assertEqual(motif.instances[3].strand, "+")
            self.assertEqual(motif.instances[4].strand, "+")
            self.assertEqual(motif.instances[5].strand, "-")
            self.assertEqual(motif.instances[6].strand, "+")
            self.assertEqual(motif.instances[0].length, 15)
            self.assertEqual(motif.instances[1].length, 15)
            self.assertEqual(motif.instances[2].length, 15)
            self.assertEqual(motif.instances[3].length, 15)
            self.assertEqual(motif.instances[4].length, 15)
            self.assertEqual(motif.instances[5].length, 15)
            self.assertEqual(motif.instances[6].length, 15)
            self.assertEqual(motif.instances[0].start, 104)
            self.assertEqual(motif.instances[1].start, 566)
            self.assertEqual(motif.instances[2].start, 585)
            self.assertEqual(motif.instances[3].start, 30)
            self.assertEqual(motif.instances[4].start, 54)
            self.assertEqual(motif.instances[5].start, 272)
            self.assertEqual(motif.instances[6].start, 214)
            self.assertEqual(motif.instances[0], "TTGACACCTGCCCAG")
            self.assertEqual(motif.instances[1], "TTGACACCTACCCTG")
            self.assertEqual(motif.instances[2], "TTGTCTCTTGCTCTG")
            self.assertEqual(motif.instances[3], "TTGACACTTGATCAG")
            self.assertEqual(motif.instances[4], "TTCACTACTCCCCTG")
            self.assertEqual(motif.instances[5], "TTGACAACGGCTGGG")
            self.assertEqual(motif.instances[6], "TTCACGCTTGCTACG")
        self.assertEqual(len(motif.alignment.sequences), 7)
        self.assertAlmostEqual(motif.alignment.sequences[0].pvalue, 7.2e-10)
        self.assertAlmostEqual(motif.alignment.sequences[1].pvalue, 2.56e-08)
        self.assertAlmostEqual(motif.alignment.sequences[2].pvalue, 1.59e-07)
        self.assertAlmostEqual(motif.alignment.sequences[3].pvalue, 2.05e-07)
        self.assertAlmostEqual(motif.alignment.sequences[4].pvalue, 3.85e-07)
        self.assertAlmostEqual(motif.alignment.sequences[5].pvalue, 5.11e-07)
        self.assertAlmostEqual(motif.alignment.sequences[6].pvalue, 8.01e-07)
        self.assertEqual(motif.alignment.sequences[0].sequence_id, "sequence_1")
        self.assertEqual(motif.alignment.sequences[1].sequence_id, "sequence_6")
        self.assertEqual(motif.alignment.sequences[2].sequence_id, "sequence_4")
        self.assertEqual(motif.alignment.sequences[3].sequence_id, "sequence_0")
        self.assertEqual(motif.alignment.sequences[4].sequence_id, "sequence_2")
        self.assertEqual(motif.alignment.sequences[5].sequence_id, "sequence_3")
        self.assertEqual(motif.alignment.sequences[6].sequence_id, "sequence_5")
        self.assertEqual(motif.alignment.sequences[0].strand, "+")
        self.assertEqual(motif.alignment.sequences[1].strand, "-")
        self.assertEqual(motif.alignment.sequences[2].strand, "-")
        self.assertEqual(motif.alignment.sequences[3].strand, "+")
        self.assertEqual(motif.alignment.sequences[4].strand, "+")
        self.assertEqual(motif.alignment.sequences[5].strand, "-")
        self.assertEqual(motif.alignment.sequences[6].strand, "+")
        self.assertEqual(motif.alignment.sequences[0].length, 15)
        self.assertEqual(motif.alignment.sequences[1].length, 15)
        self.assertEqual(motif.alignment.sequences[2].length, 15)
        self.assertEqual(motif.alignment.sequences[3].length, 15)
        self.assertEqual(motif.alignment.sequences[4].length, 15)
        self.assertEqual(motif.alignment.sequences[5].length, 15)
        self.assertEqual(motif.alignment.sequences[6].length, 15)
        self.assertEqual(motif.alignment.sequences[0].start, 104)
        self.assertEqual(motif.alignment.sequences[1].start, 566)
        self.assertEqual(motif.alignment.sequences[2].start, 585)
        self.assertEqual(motif.alignment.sequences[3].start, 30)
        self.assertEqual(motif.alignment.sequences[4].start, 54)
        self.assertEqual(motif.alignment.sequences[5].start, 272)
        self.assertEqual(motif.alignment.sequences[6].start, 214)
        self.assertEqual(motif.alignment.sequences[0], "TTGACACCTGCCCAG")
        self.assertEqual(motif.alignment.sequences[1], "TTGACACCTACCCTG")
        self.assertEqual(motif.alignment.sequences[2], "TTGTCTCTTGCTCTG")
        self.assertEqual(motif.alignment.sequences[3], "TTGACACTTGATCAG")
        self.assertEqual(motif.alignment.sequences[4], "TTCACTACTCCCCTG")
        self.assertEqual(motif.alignment.sequences[5], "TTGACAACGGCTGGG")
        self.assertEqual(motif.alignment.sequences[6], "TTCACGCTTGCTACG")
        self.assertEqual(motif.consensus, "TTGACACCTGCTCTG")
        self.assertEqual(motif.degenerate_consensus, "TTGACWCYTGCYCNG")
        self.assertTrue(
            np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        2.0,
                        2.0,
                        1.136879431433369,
                        1.4083272214176723,
                        2.0,
                        0.6212165065138244,
                        1.136879431433369,
                        1.0147718639657484,
                        1.4083272214176723,
                        0.8511651457190834,
                        1.4083272214176723,
                        1.0147718639657484,
                        0.8511651457190834,
                        0.15762900682289133,
                        2.0,
                    ]
                ),
            )
        )
        self.assertEqual(motif[1::2].consensus, "TAACGTT")

    def test_meme_parser_2(self):
        """Parsing motifs/meme.adh.classic.oops.xml file."""
        with open("motifs/meme.adh.classic.oops.xml") as stream:
            record = motifs.parse(stream, "meme")
        self.assertEqual(record.version, "5.0.1")
        self.assertEqual(record.datafile, "common/adh.s")
        self.assertEqual(record.alphabet, "ACDEFGHIKLMNPQRSTVWY")
        self.assertEqual(len(record.sequences), 33)
        self.assertEqual(record.sequences[0], "sequence_0")
        self.assertEqual(record.sequences[1], "sequence_1")
        self.assertEqual(record.sequences[2], "sequence_2")
        self.assertEqual(record.sequences[3], "sequence_3")
        self.assertEqual(record.sequences[4], "sequence_4")
        self.assertEqual(record.sequences[5], "sequence_5")
        self.assertEqual(record.sequences[6], "sequence_6")
        self.assertEqual(record.sequences[7], "sequence_7")
        self.assertEqual(record.sequences[8], "sequence_8")
        self.assertEqual(record.sequences[9], "sequence_9")
        self.assertEqual(record.sequences[10], "sequence_10")
        self.assertEqual(record.sequences[11], "sequence_11")
        self.assertEqual(record.sequences[12], "sequence_12")
        self.assertEqual(record.sequences[13], "sequence_13")
        self.assertEqual(record.sequences[14], "sequence_14")
        self.assertEqual(record.sequences[15], "sequence_15")
        self.assertEqual(record.sequences[16], "sequence_16")
        self.assertEqual(record.sequences[17], "sequence_17")
        self.assertEqual(record.sequences[18], "sequence_18")
        self.assertEqual(record.sequences[19], "sequence_19")
        self.assertEqual(record.sequences[20], "sequence_20")
        self.assertEqual(record.sequences[21], "sequence_21")
        self.assertEqual(record.sequences[22], "sequence_22")
        self.assertEqual(record.sequences[23], "sequence_23")
        self.assertEqual(record.sequences[24], "sequence_24")
        self.assertEqual(record.sequences[25], "sequence_25")
        self.assertEqual(record.sequences[26], "sequence_26")
        self.assertEqual(record.sequences[27], "sequence_27")
        self.assertEqual(record.sequences[28], "sequence_28")
        self.assertEqual(record.sequences[29], "sequence_29")
        self.assertEqual(record.sequences[30], "sequence_30")
        self.assertEqual(record.sequences[31], "sequence_31")
        self.assertEqual(record.sequences[32], "sequence_32")
        self.assertEqual(
            record.command,
            "meme common/adh.s -oc results/meme4 -mod oops -protein -nmotifs 2 -objfun classic -minw 8 -nostatus ",
        )
        self.assertEqual(len(record), 2)
        motif = record[0]
        self.assertEqual(motif.id, "motif_1")
        self.assertEqual(motif.name, "GKVALVTGAASGJGKATAKAL")
        self.assertEqual(motif.alt_id, "MEME-1")
        self.assertEqual(record["GKVALVTGAASGJGKATAKAL"], motif)
        self.assertEqual(motif.num_occurrences, 33)
        self.assertAlmostEqual(motif.evalue, 3.6e-165)
        self.assertEqual(motif.alphabet, "ACDEFGHIKLMNPQRSTVWY")
        # using the old instances property:
        with self.assertWarns(BiopythonDeprecationWarning):
            self.assertEqual(len(motif.instances), 33)
            self.assertAlmostEqual(motif.instances[0].pvalue, 8.78e-18)
            self.assertAlmostEqual(motif.instances[1].pvalue, 1.41e-17)
            self.assertAlmostEqual(motif.instances[2].pvalue, 1.42e-16)
            self.assertAlmostEqual(motif.instances[3].pvalue, 2.75e-16)
            self.assertAlmostEqual(motif.instances[4].pvalue, 3.55e-16)
            self.assertAlmostEqual(motif.instances[5].pvalue, 3.55e-16)
            self.assertAlmostEqual(motif.instances[6].pvalue, 1.74e-15)
            self.assertAlmostEqual(motif.instances[7].pvalue, 3.87e-15)
            self.assertAlmostEqual(motif.instances[8].pvalue, 4.84e-15)
            self.assertAlmostEqual(motif.instances[9].pvalue, 1.04e-14)
            self.assertAlmostEqual(motif.instances[10].pvalue, 1.58e-14)
            self.assertAlmostEqual(motif.instances[11].pvalue, 1.76e-14)
            self.assertAlmostEqual(motif.instances[12].pvalue, 2.16e-14)
            self.assertAlmostEqual(motif.instances[13].pvalue, 2.94e-14)
            self.assertAlmostEqual(motif.instances[14].pvalue, 3.25e-14)
            self.assertAlmostEqual(motif.instances[15].pvalue, 3.98e-14)
            self.assertAlmostEqual(motif.instances[16].pvalue, 4.39e-14)
            self.assertAlmostEqual(motif.instances[17].pvalue, 4.39e-14)
            self.assertAlmostEqual(motif.instances[18].pvalue, 4.85e-14)
            self.assertAlmostEqual(motif.instances[19].pvalue, 6.52e-14)
            self.assertAlmostEqual(motif.instances[20].pvalue, 1.41e-13)
            self.assertAlmostEqual(motif.instances[21].pvalue, 1.55e-13)
            self.assertAlmostEqual(motif.instances[22].pvalue, 3.07e-12)
            self.assertAlmostEqual(motif.instances[23].pvalue, 5.43e-12)
            self.assertAlmostEqual(motif.instances[24].pvalue, 6.91e-12)
            self.assertAlmostEqual(motif.instances[25].pvalue, 8.76e-12)
            self.assertAlmostEqual(motif.instances[26].pvalue, 9.48e-12)
            self.assertAlmostEqual(motif.instances[27].pvalue, 1.2e-11)
            self.assertAlmostEqual(motif.instances[28].pvalue, 1.19e-09)
            self.assertAlmostEqual(motif.instances[29].pvalue, 1.54e-09)
            self.assertAlmostEqual(motif.instances[30].pvalue, 1.99e-09)
            self.assertAlmostEqual(motif.instances[31].pvalue, 1.42e-06)
            self.assertAlmostEqual(motif.instances[32].pvalue, 3.43e-06)
            self.assertEqual(motif.instances[0].sequence_name, "BUDC_KLETE")
            self.assertEqual(motif.instances[1].sequence_name, "YINL_LISMO")
            self.assertEqual(motif.instances[2].sequence_name, "DHII_HUMAN")
            self.assertEqual(motif.instances[3].sequence_name, "HDE_CANTR")
            self.assertEqual(motif.instances[4].sequence_name, "YRTP_BACSU")
            self.assertEqual(motif.instances[5].sequence_name, "ENTA_ECOLI")
            self.assertEqual(motif.instances[6].sequence_name, "HDHA_ECOLI")
            self.assertEqual(motif.instances[7].sequence_name, "RIDH_KLEAE")
            self.assertEqual(motif.instances[8].sequence_name, "DHB2_HUMAN")
            self.assertEqual(motif.instances[9].sequence_name, "FIXR_BRAJA")
            self.assertEqual(motif.instances[10].sequence_name, "PCR_PEA")
            self.assertEqual(motif.instances[11].sequence_name, "DHCA_HUMAN")
            self.assertEqual(motif.instances[12].sequence_name, "BDH_HUMAN")
            self.assertEqual(motif.instances[13].sequence_name, "3BHD_COMTE")
            self.assertEqual(motif.instances[14].sequence_name, "DHGB_BACME")
            self.assertEqual(motif.instances[15].sequence_name, "DHMA_FLAS1")
            self.assertEqual(motif.instances[16].sequence_name, "FVT1_HUMAN")
            self.assertEqual(motif.instances[17].sequence_name, "BA72_EUBSP")
            self.assertEqual(motif.instances[18].sequence_name, "BPHB_PSEPS")
            self.assertEqual(motif.instances[19].sequence_name, "DHB3_HUMAN")
            self.assertEqual(motif.instances[20].sequence_name, "DHES_HUMAN")
            self.assertEqual(motif.instances[21].sequence_name, "AP27_MOUSE")
            self.assertEqual(motif.instances[22].sequence_name, "2BHD_STREX")
            self.assertEqual(motif.instances[23].sequence_name, "NODG_RHIME")
            self.assertEqual(motif.instances[24].sequence_name, "HMTR_LEIMA")
            self.assertEqual(motif.instances[25].sequence_name, "LIGD_PSEPA")
            self.assertEqual(motif.instances[26].sequence_name, "MAS1_AGRRA")
            self.assertEqual(motif.instances[27].sequence_name, "RFBB_NEIGO")
            self.assertEqual(motif.instances[28].sequence_name, "GUTD_ECOLI")
            self.assertEqual(motif.instances[29].sequence_name, "ADH_DROME")
            self.assertEqual(motif.instances[30].sequence_name, "FABI_ECOLI")
            self.assertEqual(motif.instances[31].sequence_name, "CSGA_MYXXA")
            self.assertEqual(motif.instances[32].sequence_name, "YURA_MYXXA")
            self.assertEqual(motif.instances[0].sequence_id, "sequence_7")
            self.assertEqual(motif.instances[1].sequence_id, "sequence_20")
            self.assertEqual(motif.instances[2].sequence_id, "sequence_10")
            self.assertEqual(motif.instances[3].sequence_id, "sequence_15")
            self.assertEqual(motif.instances[4].sequence_id, "sequence_21")
            self.assertEqual(motif.instances[5].sequence_id, "sequence_12")
            self.assertEqual(motif.instances[6].sequence_id, "sequence_16")
            self.assertEqual(motif.instances[7].sequence_id, "sequence_19")
            self.assertEqual(motif.instances[8].sequence_id, "sequence_23")
            self.assertEqual(motif.instances[9].sequence_id, "sequence_13")
            self.assertEqual(motif.instances[10].sequence_id, "sequence_30")
            self.assertEqual(motif.instances[11].sequence_id, "sequence_25")
            self.assertEqual(motif.instances[12].sequence_id, "sequence_5")
            self.assertEqual(motif.instances[13].sequence_id, "sequence_1")
            self.assertEqual(motif.instances[14].sequence_id, "sequence_9")
            self.assertEqual(motif.instances[15].sequence_id, "sequence_11")
            self.assertEqual(motif.instances[16].sequence_id, "sequence_27")
            self.assertEqual(motif.instances[17].sequence_id, "sequence_4")
            self.assertEqual(motif.instances[18].sequence_id, "sequence_6")
            self.assertEqual(motif.instances[19].sequence_id, "sequence_24")
            self.assertEqual(motif.instances[20].sequence_id, "sequence_8")
            self.assertEqual(motif.instances[21].sequence_id, "sequence_3")
            self.assertEqual(motif.instances[22].sequence_id, "sequence_0")
            self.assertEqual(motif.instances[23].sequence_id, "sequence_18")
            self.assertEqual(motif.instances[24].sequence_id, "sequence_28")
            self.assertEqual(motif.instances[25].sequence_id, "sequence_17")
            self.assertEqual(motif.instances[26].sequence_id, "sequence_29")
            self.assertEqual(motif.instances[27].sequence_id, "sequence_31")
            self.assertEqual(motif.instances[28].sequence_id, "sequence_14")
            self.assertEqual(motif.instances[29].sequence_id, "sequence_2")
            self.assertEqual(motif.instances[30].sequence_id, "sequence_26")
            self.assertEqual(motif.instances[31].sequence_id, "sequence_22")
            self.assertEqual(motif.instances[32].sequence_id, "sequence_32")
            self.assertEqual(motif.instances[0].strand, "+")
            self.assertEqual(motif.instances[1].strand, "+")
            self.assertEqual(motif.instances[2].strand, "+")
            self.assertEqual(motif.instances[3].strand, "+")
            self.assertEqual(motif.instances[4].strand, "+")
            self.assertEqual(motif.instances[5].strand, "+")
            self.assertEqual(motif.instances[6].strand, "+")
            self.assertEqual(motif.instances[7].strand, "+")
            self.assertEqual(motif.instances[8].strand, "+")
            self.assertEqual(motif.instances[9].strand, "+")
            self.assertEqual(motif.instances[10].strand, "+")
            self.assertEqual(motif.instances[11].strand, "+")
            self.assertEqual(motif.instances[12].strand, "+")
            self.assertEqual(motif.instances[13].strand, "+")
            self.assertEqual(motif.instances[14].strand, "+")
            self.assertEqual(motif.instances[15].strand, "+")
            self.assertEqual(motif.instances[16].strand, "+")
            self.assertEqual(motif.instances[17].strand, "+")
            self.assertEqual(motif.instances[18].strand, "+")
            self.assertEqual(motif.instances[19].strand, "+")
            self.assertEqual(motif.instances[20].strand, "+")
            self.assertEqual(motif.instances[21].strand, "+")
            self.assertEqual(motif.instances[22].strand, "+")
            self.assertEqual(motif.instances[23].strand, "+")
            self.assertEqual(motif.instances[24].strand, "+")
            self.assertEqual(motif.instances[25].strand, "+")
            self.assertEqual(motif.instances[26].strand, "+")
            self.assertEqual(motif.instances[27].strand, "+")
            self.assertEqual(motif.instances[28].strand, "+")
            self.assertEqual(motif.instances[29].strand, "+")
            self.assertEqual(motif.instances[30].strand, "+")
            self.assertEqual(motif.instances[31].strand, "+")
            self.assertEqual(motif.instances[32].strand, "+")
            self.assertEqual(motif.instances[0].length, 21)
            self.assertEqual(motif.instances[1].length, 21)
            self.assertEqual(motif.instances[2].length, 21)
            self.assertEqual(motif.instances[3].length, 21)
            self.assertEqual(motif.instances[4].length, 21)
            self.assertEqual(motif.instances[5].length, 21)
            self.assertEqual(motif.instances[6].length, 21)
            self.assertEqual(motif.instances[7].length, 21)
            self.assertEqual(motif.instances[8].length, 21)
            self.assertEqual(motif.instances[9].length, 21)
            self.assertEqual(motif.instances[10].length, 21)
            self.assertEqual(motif.instances[11].length, 21)
            self.assertEqual(motif.instances[12].length, 21)
            self.assertEqual(motif.instances[13].length, 21)
            self.assertEqual(motif.instances[14].length, 21)
            self.assertEqual(motif.instances[15].length, 21)
            self.assertEqual(motif.instances[16].length, 21)
            self.assertEqual(motif.instances[17].length, 21)
            self.assertEqual(motif.instances[18].length, 21)
            self.assertEqual(motif.instances[19].length, 21)
            self.assertEqual(motif.instances[20].length, 21)
            self.assertEqual(motif.instances[21].length, 21)
            self.assertEqual(motif.instances[22].length, 21)
            self.assertEqual(motif.instances[23].length, 21)
            self.assertEqual(motif.instances[24].length, 21)
            self.assertEqual(motif.instances[25].length, 21)
            self.assertEqual(motif.instances[26].length, 21)
            self.assertEqual(motif.instances[27].length, 21)
            self.assertEqual(motif.instances[28].length, 21)
            self.assertEqual(motif.instances[29].length, 21)
            self.assertEqual(motif.instances[30].length, 21)
            self.assertEqual(motif.instances[31].length, 21)
            self.assertEqual(motif.instances[32].length, 21)
            self.assertEqual(motif.instances[0].start, 2)
            self.assertEqual(motif.instances[1].start, 5)
            self.assertEqual(motif.instances[2].start, 34)
            self.assertEqual(motif.instances[3].start, 322)
            self.assertEqual(motif.instances[4].start, 6)
            self.assertEqual(motif.instances[5].start, 5)
            self.assertEqual(motif.instances[6].start, 11)
            self.assertEqual(motif.instances[7].start, 14)
            self.assertEqual(motif.instances[8].start, 82)
            self.assertEqual(motif.instances[9].start, 36)
            self.assertEqual(motif.instances[10].start, 86)
            self.assertEqual(motif.instances[11].start, 4)
            self.assertEqual(motif.instances[12].start, 55)
            self.assertEqual(motif.instances[13].start, 6)
            self.assertEqual(motif.instances[14].start, 7)
            self.assertEqual(motif.instances[15].start, 14)
            self.assertEqual(motif.instances[16].start, 32)
            self.assertEqual(motif.instances[17].start, 6)
            self.assertEqual(motif.instances[18].start, 5)
            self.assertEqual(motif.instances[19].start, 48)
            self.assertEqual(motif.instances[20].start, 2)
            self.assertEqual(motif.instances[21].start, 7)
            self.assertEqual(motif.instances[22].start, 6)
            self.assertEqual(motif.instances[23].start, 6)
            self.assertEqual(motif.instances[24].start, 6)
            self.assertEqual(motif.instances[25].start, 6)
            self.assertEqual(motif.instances[26].start, 245)
            self.assertEqual(motif.instances[27].start, 6)
            self.assertEqual(motif.instances[28].start, 2)
            self.assertEqual(motif.instances[29].start, 6)
            self.assertEqual(motif.instances[30].start, 6)
            self.assertEqual(motif.instances[31].start, 13)
            self.assertEqual(motif.instances[32].start, 116)
            self.assertEqual(motif.instances[0], "QKVALVTGAGQGIGKAIALRL")
            self.assertEqual(motif.instances[1], "NKVIIITGASSGIGKATALLL")
            self.assertEqual(motif.instances[2], "GKKVIVTGASKGIGREMAYHL")
            self.assertEqual(motif.instances[3], "DKVVLITGAGAGLGKEYAKWF")
            self.assertEqual(motif.instances[4], "HKTALITGGGRGIGRATALAL")
            self.assertEqual(motif.instances[5], "GKNVWVTGAGKGIGYATALAF")
            self.assertEqual(motif.instances[6], "GKCAIITGAGAGIGKEIAITF")
            self.assertEqual(motif.instances[7], "GKVAAITGAASGIGLECARTL")
            self.assertEqual(motif.instances[8], "QKAVLVTGGDCGLGHALCKYL")
            self.assertEqual(motif.instances[9], "PKVMLLTGASRGIGHATAKLF")
            self.assertEqual(motif.instances[10], "KGNVVITGASSGLGLATAKAL")
            self.assertEqual(motif.instances[11], "IHVALVTGGNKGIGLAIVRDL")
            self.assertEqual(motif.instances[12], "SKAVLVTGCDSGFGFSLAKHL")
            self.assertEqual(motif.instances[13], "GKVALVTGGASGVGLEVVKLL")
            self.assertEqual(motif.instances[14], "GKVVVITGSSTGLGKSMAIRF")
            self.assertEqual(motif.instances[15], "GKAAIVTGAAGGIGRATVEAY")
            self.assertEqual(motif.instances[16], "GAHVVVTGGSSGIGKCIAIEC")
            self.assertEqual(motif.instances[17], "DKVTIITGGTRGIGFAAAKIF")
            self.assertEqual(motif.instances[18], "GEAVLITGGASGLGRALVDRF")
            self.assertEqual(motif.instances[19], "GQWAVITGAGDGIGKAYSFEL")
            self.assertEqual(motif.instances[20], "RTVVLITGCSSGIGLHLAVRL")
            self.assertEqual(motif.instances[21], "GLRALVTGAGKGIGRDTVKAL")
            self.assertEqual(motif.instances[22], "GKTVIITGGARGLGAEAARQA")
            self.assertEqual(motif.instances[23], "GRKALVTGASGAIGGAIARVL")
            self.assertEqual(motif.instances[24], "VPVALVTGAAKRLGRSIAEGL")
            self.assertEqual(motif.instances[25], "DQVAFITGGASGAGFGQAKVF")
            self.assertEqual(motif.instances[26], "SPVILVSGSNRGVGKAIAEDL")
            self.assertEqual(motif.instances[27], "KKNILVTGGAGFIGSAVVRHI")
            self.assertEqual(motif.instances[28], "NQVAVVIGGGQTLGAFLCHGL")
            self.assertEqual(motif.instances[29], "NKNVIFVAGLGGIGLDTSKEL")
            self.assertEqual(motif.instances[30], "GKRILVTGVASKLSIAYGIAQ")
            self.assertEqual(motif.instances[31], "VDVLINNAGVSGLWCALGDVD")
            self.assertEqual(motif.instances[32], "IIDTNVTGAAATLSAVLPQMV")
        self.assertEqual(len(motif.alignment.sequences), 33)
        self.assertAlmostEqual(motif.alignment.sequences[0].pvalue, 8.78e-18)
        self.assertAlmostEqual(motif.alignment.sequences[1].pvalue, 1.41e-17)
        self.assertAlmostEqual(motif.alignment.sequences[2].pvalue, 1.42e-16)
        self.assertAlmostEqual(motif.alignment.sequences[3].pvalue, 2.75e-16)
        self.assertAlmostEqual(motif.alignment.sequences[4].pvalue, 3.55e-16)
        self.assertAlmostEqual(motif.alignment.sequences[5].pvalue, 3.55e-16)
        self.assertAlmostEqual(motif.alignment.sequences[6].pvalue, 1.74e-15)
        self.assertAlmostEqual(motif.alignment.sequences[7].pvalue, 3.87e-15)
        self.assertAlmostEqual(motif.alignment.sequences[8].pvalue, 4.84e-15)
        self.assertAlmostEqual(motif.alignment.sequences[9].pvalue, 1.04e-14)
        self.assertAlmostEqual(motif.alignment.sequences[10].pvalue, 1.58e-14)
        self.assertAlmostEqual(motif.alignment.sequences[11].pvalue, 1.76e-14)
        self.assertAlmostEqual(motif.alignment.sequences[12].pvalue, 2.16e-14)
        self.assertAlmostEqual(motif.alignment.sequences[13].pvalue, 2.94e-14)
        self.assertAlmostEqual(motif.alignment.sequences[14].pvalue, 3.25e-14)
        self.assertAlmostEqual(motif.alignment.sequences[15].pvalue, 3.98e-14)
        self.assertAlmostEqual(motif.alignment.sequences[16].pvalue, 4.39e-14)
        self.assertAlmostEqual(motif.alignment.sequences[17].pvalue, 4.39e-14)
        self.assertAlmostEqual(motif.alignment.sequences[18].pvalue, 4.85e-14)
        self.assertAlmostEqual(motif.alignment.sequences[19].pvalue, 6.52e-14)
        self.assertAlmostEqual(motif.alignment.sequences[20].pvalue, 1.41e-13)
        self.assertAlmostEqual(motif.alignment.sequences[21].pvalue, 1.55e-13)
        self.assertAlmostEqual(motif.alignment.sequences[22].pvalue, 3.07e-12)
        self.assertAlmostEqual(motif.alignment.sequences[23].pvalue, 5.43e-12)
        self.assertAlmostEqual(motif.alignment.sequences[24].pvalue, 6.91e-12)
        self.assertAlmostEqual(motif.alignment.sequences[25].pvalue, 8.76e-12)
        self.assertAlmostEqual(motif.alignment.sequences[26].pvalue, 9.48e-12)
        self.assertAlmostEqual(motif.alignment.sequences[27].pvalue, 1.2e-11)
        self.assertAlmostEqual(motif.alignment.sequences[28].pvalue, 1.19e-09)
        self.assertAlmostEqual(motif.alignment.sequences[29].pvalue, 1.54e-09)
        self.assertAlmostEqual(motif.alignment.sequences[30].pvalue, 1.99e-09)
        self.assertAlmostEqual(motif.alignment.sequences[31].pvalue, 1.42e-06)
        self.assertAlmostEqual(motif.alignment.sequences[32].pvalue, 3.43e-06)
        self.assertEqual(motif.alignment.sequences[0].sequence_name, "BUDC_KLETE")
        self.assertEqual(motif.alignment.sequences[1].sequence_name, "YINL_LISMO")
        self.assertEqual(motif.alignment.sequences[2].sequence_name, "DHII_HUMAN")
        self.assertEqual(motif.alignment.sequences[3].sequence_name, "HDE_CANTR")
        self.assertEqual(motif.alignment.sequences[4].sequence_name, "YRTP_BACSU")
        self.assertEqual(motif.alignment.sequences[5].sequence_name, "ENTA_ECOLI")
        self.assertEqual(motif.alignment.sequences[6].sequence_name, "HDHA_ECOLI")
        self.assertEqual(motif.alignment.sequences[7].sequence_name, "RIDH_KLEAE")
        self.assertEqual(motif.alignment.sequences[8].sequence_name, "DHB2_HUMAN")
        self.assertEqual(motif.alignment.sequences[9].sequence_name, "FIXR_BRAJA")
        self.assertEqual(motif.alignment.sequences[10].sequence_name, "PCR_PEA")
        self.assertEqual(motif.alignment.sequences[11].sequence_name, "DHCA_HUMAN")
        self.assertEqual(motif.alignment.sequences[12].sequence_name, "BDH_HUMAN")
        self.assertEqual(motif.alignment.sequences[13].sequence_name, "3BHD_COMTE")
        self.assertEqual(motif.alignment.sequences[14].sequence_name, "DHGB_BACME")
        self.assertEqual(motif.alignment.sequences[15].sequence_name, "DHMA_FLAS1")
        self.assertEqual(motif.alignment.sequences[16].sequence_name, "FVT1_HUMAN")
        self.assertEqual(motif.alignment.sequences[17].sequence_name, "BA72_EUBSP")
        self.assertEqual(motif.alignment.sequences[18].sequence_name, "BPHB_PSEPS")
        self.assertEqual(motif.alignment.sequences[19].sequence_name, "DHB3_HUMAN")
        self.assertEqual(motif.alignment.sequences[20].sequence_name, "DHES_HUMAN")
        self.assertEqual(motif.alignment.sequences[21].sequence_name, "AP27_MOUSE")
        self.assertEqual(motif.alignment.sequences[22].sequence_name, "2BHD_STREX")
        self.assertEqual(motif.alignment.sequences[23].sequence_name, "NODG_RHIME")
        self.assertEqual(motif.alignment.sequences[24].sequence_name, "HMTR_LEIMA")
        self.assertEqual(motif.alignment.sequences[25].sequence_name, "LIGD_PSEPA")
        self.assertEqual(motif.alignment.sequences[26].sequence_name, "MAS1_AGRRA")
        self.assertEqual(motif.alignment.sequences[27].sequence_name, "RFBB_NEIGO")
        self.assertEqual(motif.alignment.sequences[28].sequence_name, "GUTD_ECOLI")
        self.assertEqual(motif.alignment.sequences[29].sequence_name, "ADH_DROME")
        self.assertEqual(motif.alignment.sequences[30].sequence_name, "FABI_ECOLI")
        self.assertEqual(motif.alignment.sequences[31].sequence_name, "CSGA_MYXXA")
        self.assertEqual(motif.alignment.sequences[32].sequence_name, "YURA_MYXXA")
        self.assertEqual(motif.alignment.sequences[0].sequence_id, "sequence_7")
        self.assertEqual(motif.alignment.sequences[1].sequence_id, "sequence_20")
        self.assertEqual(motif.alignment.sequences[2].sequence_id, "sequence_10")
        self.assertEqual(motif.alignment.sequences[3].sequence_id, "sequence_15")
        self.assertEqual(motif.alignment.sequences[4].sequence_id, "sequence_21")
        self.assertEqual(motif.alignment.sequences[5].sequence_id, "sequence_12")
        self.assertEqual(motif.alignment.sequences[6].sequence_id, "sequence_16")
        self.assertEqual(motif.alignment.sequences[7].sequence_id, "sequence_19")
        self.assertEqual(motif.alignment.sequences[8].sequence_id, "sequence_23")
        self.assertEqual(motif.alignment.sequences[9].sequence_id, "sequence_13")
        self.assertEqual(motif.alignment.sequences[10].sequence_id, "sequence_30")
        self.assertEqual(motif.alignment.sequences[11].sequence_id, "sequence_25")
        self.assertEqual(motif.alignment.sequences[12].sequence_id, "sequence_5")
        self.assertEqual(motif.alignment.sequences[13].sequence_id, "sequence_1")
        self.assertEqual(motif.alignment.sequences[14].sequence_id, "sequence_9")
        self.assertEqual(motif.alignment.sequences[15].sequence_id, "sequence_11")
        self.assertEqual(motif.alignment.sequences[16].sequence_id, "sequence_27")
        self.assertEqual(motif.alignment.sequences[17].sequence_id, "sequence_4")
        self.assertEqual(motif.alignment.sequences[18].sequence_id, "sequence_6")
        self.assertEqual(motif.alignment.sequences[19].sequence_id, "sequence_24")
        self.assertEqual(motif.alignment.sequences[20].sequence_id, "sequence_8")
        self.assertEqual(motif.alignment.sequences[21].sequence_id, "sequence_3")
        self.assertEqual(motif.alignment.sequences[22].sequence_id, "sequence_0")
        self.assertEqual(motif.alignment.sequences[23].sequence_id, "sequence_18")
        self.assertEqual(motif.alignment.sequences[24].sequence_id, "sequence_28")
        self.assertEqual(motif.alignment.sequences[25].sequence_id, "sequence_17")
        self.assertEqual(motif.alignment.sequences[26].sequence_id, "sequence_29")
        self.assertEqual(motif.alignment.sequences[27].sequence_id, "sequence_31")
        self.assertEqual(motif.alignment.sequences[28].sequence_id, "sequence_14")
        self.assertEqual(motif.alignment.sequences[29].sequence_id, "sequence_2")
        self.assertEqual(motif.alignment.sequences[30].sequence_id, "sequence_26")
        self.assertEqual(motif.alignment.sequences[31].sequence_id, "sequence_22")
        self.assertEqual(motif.alignment.sequences[32].sequence_id, "sequence_32")
        self.assertEqual(motif.alignment.sequences[0].strand, "+")
        self.assertEqual(motif.alignment.sequences[1].strand, "+")
        self.assertEqual(motif.alignment.sequences[2].strand, "+")
        self.assertEqual(motif.alignment.sequences[3].strand, "+")
        self.assertEqual(motif.alignment.sequences[4].strand, "+")
        self.assertEqual(motif.alignment.sequences[5].strand, "+")
        self.assertEqual(motif.alignment.sequences[6].strand, "+")
        self.assertEqual(motif.alignment.sequences[7].strand, "+")
        self.assertEqual(motif.alignment.sequences[8].strand, "+")
        self.assertEqual(motif.alignment.sequences[9].strand, "+")
        self.assertEqual(motif.alignment.sequences[10].strand, "+")
        self.assertEqual(motif.alignment.sequences[11].strand, "+")
        self.assertEqual(motif.alignment.sequences[12].strand, "+")
        self.assertEqual(motif.alignment.sequences[13].strand, "+")
        self.assertEqual(motif.alignment.sequences[14].strand, "+")
        self.assertEqual(motif.alignment.sequences[15].strand, "+")
        self.assertEqual(motif.alignment.sequences[16].strand, "+")
        self.assertEqual(motif.alignment.sequences[17].strand, "+")
        self.assertEqual(motif.alignment.sequences[18].strand, "+")
        self.assertEqual(motif.alignment.sequences[19].strand, "+")
        self.assertEqual(motif.alignment.sequences[20].strand, "+")
        self.assertEqual(motif.alignment.sequences[21].strand, "+")
        self.assertEqual(motif.alignment.sequences[22].strand, "+")
        self.assertEqual(motif.alignment.sequences[23].strand, "+")
        self.assertEqual(motif.alignment.sequences[24].strand, "+")
        self.assertEqual(motif.alignment.sequences[25].strand, "+")
        self.assertEqual(motif.alignment.sequences[26].strand, "+")
        self.assertEqual(motif.alignment.sequences[27].strand, "+")
        self.assertEqual(motif.alignment.sequences[28].strand, "+")
        self.assertEqual(motif.alignment.sequences[29].strand, "+")
        self.assertEqual(motif.alignment.sequences[30].strand, "+")
        self.assertEqual(motif.alignment.sequences[31].strand, "+")
        self.assertEqual(motif.alignment.sequences[32].strand, "+")
        self.assertEqual(motif.alignment.sequences[0].length, 21)
        self.assertEqual(motif.alignment.sequences[1].length, 21)
        self.assertEqual(motif.alignment.sequences[2].length, 21)
        self.assertEqual(motif.alignment.sequences[3].length, 21)
        self.assertEqual(motif.alignment.sequences[4].length, 21)
        self.assertEqual(motif.alignment.sequences[5].length, 21)
        self.assertEqual(motif.alignment.sequences[6].length, 21)
        self.assertEqual(motif.alignment.sequences[7].length, 21)
        self.assertEqual(motif.alignment.sequences[8].length, 21)
        self.assertEqual(motif.alignment.sequences[9].length, 21)
        self.assertEqual(motif.alignment.sequences[10].length, 21)
        self.assertEqual(motif.alignment.sequences[11].length, 21)
        self.assertEqual(motif.alignment.sequences[12].length, 21)
        self.assertEqual(motif.alignment.sequences[13].length, 21)
        self.assertEqual(motif.alignment.sequences[14].length, 21)
        self.assertEqual(motif.alignment.sequences[15].length, 21)
        self.assertEqual(motif.alignment.sequences[16].length, 21)
        self.assertEqual(motif.alignment.sequences[17].length, 21)
        self.assertEqual(motif.alignment.sequences[18].length, 21)
        self.assertEqual(motif.alignment.sequences[19].length, 21)
        self.assertEqual(motif.alignment.sequences[20].length, 21)
        self.assertEqual(motif.alignment.sequences[21].length, 21)
        self.assertEqual(motif.alignment.sequences[22].length, 21)
        self.assertEqual(motif.alignment.sequences[23].length, 21)
        self.assertEqual(motif.alignment.sequences[24].length, 21)
        self.assertEqual(motif.alignment.sequences[25].length, 21)
        self.assertEqual(motif.alignment.sequences[26].length, 21)
        self.assertEqual(motif.alignment.sequences[27].length, 21)
        self.assertEqual(motif.alignment.sequences[28].length, 21)
        self.assertEqual(motif.alignment.sequences[29].length, 21)
        self.assertEqual(motif.alignment.sequences[30].length, 21)
        self.assertEqual(motif.alignment.sequences[31].length, 21)
        self.assertEqual(motif.alignment.sequences[32].length, 21)
        self.assertEqual(motif.alignment.sequences[0].start, 2)
        self.assertEqual(motif.alignment.sequences[1].start, 5)
        self.assertEqual(motif.alignment.sequences[2].start, 34)
        self.assertEqual(motif.alignment.sequences[3].start, 322)
        self.assertEqual(motif.alignment.sequences[4].start, 6)
        self.assertEqual(motif.alignment.sequences[5].start, 5)
        self.assertEqual(motif.alignment.sequences[6].start, 11)
        self.assertEqual(motif.alignment.sequences[7].start, 14)
        self.assertEqual(motif.alignment.sequences[8].start, 82)
        self.assertEqual(motif.alignment.sequences[9].start, 36)
        self.assertEqual(motif.alignment.sequences[10].start, 86)
        self.assertEqual(motif.alignment.sequences[11].start, 4)
        self.assertEqual(motif.alignment.sequences[12].start, 55)
        self.assertEqual(motif.alignment.sequences[13].start, 6)
        self.assertEqual(motif.alignment.sequences[14].start, 7)
        self.assertEqual(motif.alignment.sequences[15].start, 14)
        self.assertEqual(motif.alignment.sequences[16].start, 32)
        self.assertEqual(motif.alignment.sequences[17].start, 6)
        self.assertEqual(motif.alignment.sequences[18].start, 5)
        self.assertEqual(motif.alignment.sequences[19].start, 48)
        self.assertEqual(motif.alignment.sequences[20].start, 2)
        self.assertEqual(motif.alignment.sequences[21].start, 7)
        self.assertEqual(motif.alignment.sequences[22].start, 6)
        self.assertEqual(motif.alignment.sequences[23].start, 6)
        self.assertEqual(motif.alignment.sequences[24].start, 6)
        self.assertEqual(motif.alignment.sequences[25].start, 6)
        self.assertEqual(motif.alignment.sequences[26].start, 245)
        self.assertEqual(motif.alignment.sequences[27].start, 6)
        self.assertEqual(motif.alignment.sequences[28].start, 2)
        self.assertEqual(motif.alignment.sequences[29].start, 6)
        self.assertEqual(motif.alignment.sequences[30].start, 6)
        self.assertEqual(motif.alignment.sequences[31].start, 13)
        self.assertEqual(motif.alignment.sequences[32].start, 116)
        self.assertEqual(motif.alignment.sequences[0], "QKVALVTGAGQGIGKAIALRL")
        self.assertEqual(motif.alignment.sequences[1], "NKVIIITGASSGIGKATALLL")
        self.assertEqual(motif.alignment.sequences[2], "GKKVIVTGASKGIGREMAYHL")
        self.assertEqual(motif.alignment.sequences[3], "DKVVLITGAGAGLGKEYAKWF")
        self.assertEqual(motif.alignment.sequences[4], "HKTALITGGGRGIGRATALAL")
        self.assertEqual(motif.alignment.sequences[5], "GKNVWVTGAGKGIGYATALAF")
        self.assertEqual(motif.alignment.sequences[6], "GKCAIITGAGAGIGKEIAITF")
        self.assertEqual(motif.alignment.sequences[7], "GKVAAITGAASGIGLECARTL")
        self.assertEqual(motif.alignment.sequences[8], "QKAVLVTGGDCGLGHALCKYL")
        self.assertEqual(motif.alignment.sequences[9], "PKVMLLTGASRGIGHATAKLF")
        self.assertEqual(motif.alignment.sequences[10], "KGNVVITGASSGLGLATAKAL")
        self.assertEqual(motif.alignment.sequences[11], "IHVALVTGGNKGIGLAIVRDL")
        self.assertEqual(motif.alignment.sequences[12], "SKAVLVTGCDSGFGFSLAKHL")
        self.assertEqual(motif.alignment.sequences[13], "GKVALVTGGASGVGLEVVKLL")
        self.assertEqual(motif.alignment.sequences[14], "GKVVVITGSSTGLGKSMAIRF")
        self.assertEqual(motif.alignment.sequences[15], "GKAAIVTGAAGGIGRATVEAY")
        self.assertEqual(motif.alignment.sequences[16], "GAHVVVTGGSSGIGKCIAIEC")
        self.assertEqual(motif.alignment.sequences[17], "DKVTIITGGTRGIGFAAAKIF")
        self.assertEqual(motif.alignment.sequences[18], "GEAVLITGGASGLGRALVDRF")
        self.assertEqual(motif.alignment.sequences[19], "GQWAVITGAGDGIGKAYSFEL")
        self.assertEqual(motif.alignment.sequences[20], "RTVVLITGCSSGIGLHLAVRL")
        self.assertEqual(motif.alignment.sequences[21], "GLRALVTGAGKGIGRDTVKAL")
        self.assertEqual(motif.alignment.sequences[22], "GKTVIITGGARGLGAEAARQA")
        self.assertEqual(motif.alignment.sequences[23], "GRKALVTGASGAIGGAIARVL")
        self.assertEqual(motif.alignment.sequences[24], "VPVALVTGAAKRLGRSIAEGL")
        self.assertEqual(motif.alignment.sequences[25], "DQVAFITGGASGAGFGQAKVF")
        self.assertEqual(motif.alignment.sequences[26], "SPVILVSGSNRGVGKAIAEDL")
        self.assertEqual(motif.alignment.sequences[27], "KKNILVTGGAGFIGSAVVRHI")
        self.assertEqual(motif.alignment.sequences[28], "NQVAVVIGGGQTLGAFLCHGL")
        self.assertEqual(motif.alignment.sequences[29], "NKNVIFVAGLGGIGLDTSKEL")
        self.assertEqual(motif.alignment.sequences[30], "GKRILVTGVASKLSIAYGIAQ")
        self.assertEqual(motif.alignment.sequences[31], "VDVLINNAGVSGLWCALGDVD")
        self.assertEqual(motif.alignment.sequences[32], "IIDTNVTGAAATLSAVLPQMV")
        self.assertEqual(motif.consensus, "GKVALVTGAASGIGKATAKAL")
        self.assertEqual(motif[2:8].consensus, "VALVTG")
        motif = record[1]
        self.assertEqual(motif.name, "VGNPGASAYSASKAAVRGLTESLALELAP")
        self.assertEqual(motif.alt_id, "MEME-2")
        self.assertEqual(record["VGNPGASAYSASKAAVRGLTESLALELAP"], motif)
        self.assertEqual(motif.num_occurrences, 33)
        self.assertAlmostEqual(motif.evalue, 3.1e-130)
        self.assertEqual(motif.alphabet, "ACDEFGHIKLMNPQRSTVWY")
        # using the old instances property:
        with self.assertWarns(BiopythonDeprecationWarning):
            self.assertEqual(len(motif.instances), 33)
            self.assertAlmostEqual(motif.instances[0].pvalue, 2.09e-21)
            self.assertAlmostEqual(motif.instances[1].pvalue, 7.63e-20)
            self.assertAlmostEqual(motif.instances[2].pvalue, 6.49e-19)
            self.assertAlmostEqual(motif.instances[3].pvalue, 1.92e-18)
            self.assertAlmostEqual(motif.instances[4].pvalue, 5.46e-18)
            self.assertAlmostEqual(motif.instances[5].pvalue, 6.21e-18)
            self.assertAlmostEqual(motif.instances[6].pvalue, 4.52e-17)
            self.assertAlmostEqual(motif.instances[7].pvalue, 4.52e-17)
            self.assertAlmostEqual(motif.instances[8].pvalue, 9.21e-17)
            self.assertAlmostEqual(motif.instances[9].pvalue, 1.65e-16)
            self.assertAlmostEqual(motif.instances[10].pvalue, 2.07e-16)
            self.assertAlmostEqual(motif.instances[11].pvalue, 3.65e-16)
            self.assertAlmostEqual(motif.instances[12].pvalue, 5.7e-16)
            self.assertAlmostEqual(motif.instances[13].pvalue, 5.7e-16)
            self.assertAlmostEqual(motif.instances[14].pvalue, 7.93e-16)
            self.assertAlmostEqual(motif.instances[15].pvalue, 8.85e-16)
            self.assertAlmostEqual(motif.instances[16].pvalue, 1.1e-15)
            self.assertAlmostEqual(motif.instances[17].pvalue, 1.69e-15)
            self.assertAlmostEqual(motif.instances[18].pvalue, 3.54e-15)
            self.assertAlmostEqual(motif.instances[19].pvalue, 4.83e-15)
            self.assertAlmostEqual(motif.instances[20].pvalue, 7.27e-15)
            self.assertAlmostEqual(motif.instances[21].pvalue, 9.85e-15)
            self.assertAlmostEqual(motif.instances[22].pvalue, 2.41e-14)
            self.assertAlmostEqual(motif.instances[23].pvalue, 2.66e-14)
            self.assertAlmostEqual(motif.instances[24].pvalue, 1.22e-13)
            self.assertAlmostEqual(motif.instances[25].pvalue, 5.18e-13)
            self.assertAlmostEqual(motif.instances[26].pvalue, 1.24e-12)
            self.assertAlmostEqual(motif.instances[27].pvalue, 1.35e-12)
            self.assertAlmostEqual(motif.instances[28].pvalue, 5.59e-12)
            self.assertAlmostEqual(motif.instances[29].pvalue, 1.44e-10)
            self.assertAlmostEqual(motif.instances[30].pvalue, 1.61e-08)
            self.assertAlmostEqual(motif.instances[31].pvalue, 4.26e-08)
            self.assertAlmostEqual(motif.instances[32].pvalue, 1.16e-07)
            self.assertEqual(motif.instances[0].sequence_name, "BUDC_KLETE")
            self.assertEqual(motif.instances[1].sequence_name, "NODG_RHIME")
            self.assertEqual(motif.instances[2].sequence_name, "FVT1_HUMAN")
            self.assertEqual(motif.instances[3].sequence_name, "DHES_HUMAN")
            self.assertEqual(motif.instances[4].sequence_name, "DHB3_HUMAN")
            self.assertEqual(motif.instances[5].sequence_name, "YRTP_BACSU")
            self.assertEqual(motif.instances[6].sequence_name, "HMTR_LEIMA")
            self.assertEqual(motif.instances[7].sequence_name, "HDE_CANTR")
            self.assertEqual(motif.instances[8].sequence_name, "DHGB_BACME")
            self.assertEqual(motif.instances[9].sequence_name, "GUTD_ECOLI")
            self.assertEqual(motif.instances[10].sequence_name, "3BHD_COMTE")
            self.assertEqual(motif.instances[11].sequence_name, "DHII_HUMAN")
            self.assertEqual(motif.instances[12].sequence_name, "BPHB_PSEPS")
            self.assertEqual(motif.instances[13].sequence_name, "AP27_MOUSE")
            self.assertEqual(motif.instances[14].sequence_name, "BDH_HUMAN")
            self.assertEqual(motif.instances[15].sequence_name, "YINL_LISMO")
            self.assertEqual(motif.instances[16].sequence_name, "FIXR_BRAJA")
            self.assertEqual(motif.instances[17].sequence_name, "2BHD_STREX")
            self.assertEqual(motif.instances[18].sequence_name, "RFBB_NEIGO")
            self.assertEqual(motif.instances[19].sequence_name, "YURA_MYXXA")
            self.assertEqual(motif.instances[20].sequence_name, "RIDH_KLEAE")
            self.assertEqual(motif.instances[21].sequence_name, "DHMA_FLAS1")
            self.assertEqual(motif.instances[22].sequence_name, "DHB2_HUMAN")
            self.assertEqual(motif.instances[23].sequence_name, "HDHA_ECOLI")
            self.assertEqual(motif.instances[24].sequence_name, "ENTA_ECOLI")
            self.assertEqual(motif.instances[25].sequence_name, "LIGD_PSEPA")
            self.assertEqual(motif.instances[26].sequence_name, "CSGA_MYXXA")
            self.assertEqual(motif.instances[27].sequence_name, "BA72_EUBSP")
            self.assertEqual(motif.instances[28].sequence_name, "ADH_DROME")
            self.assertEqual(motif.instances[29].sequence_name, "MAS1_AGRRA")
            self.assertEqual(motif.instances[30].sequence_name, "PCR_PEA")
            self.assertEqual(motif.instances[31].sequence_name, "FABI_ECOLI")
            self.assertEqual(motif.instances[32].sequence_name, "DHCA_HUMAN")
            self.assertEqual(motif.instances[0].sequence_id, "sequence_7")
            self.assertEqual(motif.instances[1].sequence_id, "sequence_18")
            self.assertEqual(motif.instances[2].sequence_id, "sequence_27")
            self.assertEqual(motif.instances[3].sequence_id, "sequence_8")
            self.assertEqual(motif.instances[4].sequence_id, "sequence_24")
            self.assertEqual(motif.instances[5].sequence_id, "sequence_21")
            self.assertEqual(motif.instances[6].sequence_id, "sequence_28")
            self.assertEqual(motif.instances[7].sequence_id, "sequence_15")
            self.assertEqual(motif.instances[8].sequence_id, "sequence_9")
            self.assertEqual(motif.instances[9].sequence_id, "sequence_14")
            self.assertEqual(motif.instances[10].sequence_id, "sequence_1")
            self.assertEqual(motif.instances[11].sequence_id, "sequence_10")
            self.assertEqual(motif.instances[12].sequence_id, "sequence_6")
            self.assertEqual(motif.instances[13].sequence_id, "sequence_3")
            self.assertEqual(motif.instances[14].sequence_id, "sequence_5")
            self.assertEqual(motif.instances[15].sequence_id, "sequence_20")
            self.assertEqual(motif.instances[16].sequence_id, "sequence_13")
            self.assertEqual(motif.instances[17].sequence_id, "sequence_0")
            self.assertEqual(motif.instances[18].sequence_id, "sequence_31")
            self.assertEqual(motif.instances[19].sequence_id, "sequence_32")
            self.assertEqual(motif.instances[20].sequence_id, "sequence_19")
            self.assertEqual(motif.instances[21].sequence_id, "sequence_11")
            self.assertEqual(motif.instances[22].sequence_id, "sequence_23")
            self.assertEqual(motif.instances[23].sequence_id, "sequence_16")
            self.assertEqual(motif.instances[24].sequence_id, "sequence_12")
            self.assertEqual(motif.instances[25].sequence_id, "sequence_17")
            self.assertEqual(motif.instances[26].sequence_id, "sequence_22")
            self.assertEqual(motif.instances[27].sequence_id, "sequence_4")
            self.assertEqual(motif.instances[28].sequence_id, "sequence_2")
            self.assertEqual(motif.instances[29].sequence_id, "sequence_29")
            self.assertEqual(motif.instances[30].sequence_id, "sequence_30")
            self.assertEqual(motif.instances[31].sequence_id, "sequence_26")
            self.assertEqual(motif.instances[32].sequence_id, "sequence_25")
            self.assertEqual(motif.instances[0].start, 144)
            self.assertEqual(motif.instances[1].start, 144)
            self.assertEqual(motif.instances[2].start, 178)
            self.assertEqual(motif.instances[3].start, 147)
            self.assertEqual(motif.instances[4].start, 190)
            self.assertEqual(motif.instances[5].start, 147)
            self.assertEqual(motif.instances[6].start, 185)
            self.assertEqual(motif.instances[7].start, 459)
            self.assertEqual(motif.instances[8].start, 152)
            self.assertEqual(motif.instances[9].start, 146)
            self.assertEqual(motif.instances[10].start, 143)
            self.assertEqual(motif.instances[11].start, 175)
            self.assertEqual(motif.instances[12].start, 145)
            self.assertEqual(motif.instances[13].start, 141)
            self.assertEqual(motif.instances[14].start, 200)
            self.assertEqual(motif.instances[15].start, 146)
            self.assertEqual(motif.instances[16].start, 181)
            self.assertEqual(motif.instances[17].start, 144)
            self.assertEqual(motif.instances[18].start, 157)
            self.assertEqual(motif.instances[19].start, 152)
            self.assertEqual(motif.instances[20].start, 152)
            self.assertEqual(motif.instances[21].start, 157)
            self.assertEqual(motif.instances[22].start, 224)
            self.assertEqual(motif.instances[23].start, 151)
            self.assertEqual(motif.instances[24].start, 136)
            self.assertEqual(motif.instances[25].start, 149)
            self.assertEqual(motif.instances[26].start, 80)
            self.assertEqual(motif.instances[27].start, 149)
            self.assertEqual(motif.instances[28].start, 144)
            self.assertEqual(motif.instances[29].start, 384)
            self.assertEqual(motif.instances[30].start, 18)
            self.assertEqual(motif.instances[31].start, 177)
            self.assertEqual(motif.instances[32].start, 144)
            self.assertEqual(motif.instances[0], "VGNPELAVYSSSKFAVRGLTQTAARDLAP")
            self.assertEqual(motif.instances[1], "IGNPGQTNYCASKAGMIGFSKSLAQEIAT")
            self.assertEqual(motif.instances[2], "LGLFGFTAYSASKFAIRGLAEALQMEVKP")
            self.assertEqual(motif.instances[3], "MGLPFNDVYCASKFALEGLCESLAVLLLP")
            self.assertEqual(motif.instances[4], "FPWPLYSMYSASKAFVCAFSKALQEEYKA")
            self.assertEqual(motif.instances[5], "RGAAVTSAYSASKFAVLGLTESLMQEVRK")
            self.assertEqual(motif.instances[6], "QPLLGYTIYTMAKGALEGLTRSAALELAP")
            self.assertEqual(motif.instances[7], "YGNFGQANYSSSKAGILGLSKTMAIEGAK")
            self.assertEqual(motif.instances[8], "IPWPLFVHYAASKGGMKLMTETLALEYAP")
            self.assertEqual(motif.instances[9], "VGSKHNSGYSAAKFGGVGLTQSLALDLAE")
            self.assertEqual(motif.instances[10], "LPIEQYAGYSASKAAVSALTRAAALSCRK")
            self.assertEqual(motif.instances[11], "VAYPMVAAYSASKFALDGFFSSIRKEYSV")
            self.assertEqual(motif.instances[12], "YPNGGGPLYTAAKQAIVGLVRELAFELAP")
            self.assertEqual(motif.instances[13], "VTFPNLITYSSTKGAMTMLTKAMAMELGP")
            self.assertEqual(motif.instances[14], "MANPARSPYCITKFGVEAFSDCLRYEMYP")
            self.assertEqual(motif.instances[15], "KAYPGGAVYGATKWAVRDLMEVLRMESAQ")
            self.assertEqual(motif.instances[16], "VHPFAGSAYATSKAALASLTRELAHDYAP")
            self.assertEqual(motif.instances[17], "MGLALTSSYGASKWGVRGLSKLAAVELGT")
            self.assertEqual(motif.instances[18], "TPYAPSSPYSASKAAADHLVRAWQRTYRL")
            self.assertEqual(motif.instances[19], "FRGLPATRYSASKAFLSTFMESLRVDLRG")
            self.assertEqual(motif.instances[20], "VPVIWEPVYTASKFAVQAFVHTTRRQVAQ")
            self.assertEqual(motif.instances[21], "MAEPEAAAYVAAKGGVAMLTRAMAVDLAR")
            self.assertEqual(motif.instances[22], "APMERLASYGSSKAAVTMFSSVMRLELSK")
            self.assertEqual(motif.instances[23], "NKNINMTSYASSKAAASHLVRNMAFDLGE")
            self.assertEqual(motif.instances[24], "TPRIGMSAYGASKAALKSLALSVGLELAG")
            self.assertEqual(motif.instances[25], "MGSALAGPYSAAKAASINLMEGYRQGLEK")
            self.assertEqual(motif.instances[26], "NTDGGAYAYRMSKAALNMAVRSMSTDLRP")
            self.assertEqual(motif.instances[27], "FGSLSGVGYPASKASVIGLTHGLGREIIR")
            self.assertEqual(motif.instances[28], "NAIYQVPVYSGTKAAVVNFTSSLAKLAPI")
            self.assertEqual(motif.instances[29], "RVLNPLVGYNMTKHALGGLTKTTQHVGWD")
            self.assertEqual(motif.instances[30], "EGKIGASLKDSTLFGVSSLSDSLKGDFTS")
            self.assertEqual(motif.instances[31], "MGPEGVRVNAISAGPIRTLAASGIKDFRK")
            self.assertEqual(motif.instances[32], "RALKSCSPELQQKFRSETITEEELVGLMN")
        self.assertEqual(len(motif.alignment.sequences), 33)
        self.assertAlmostEqual(motif.alignment.sequences[0].pvalue, 2.09e-21)
        self.assertAlmostEqual(motif.alignment.sequences[1].pvalue, 7.63e-20)
        self.assertAlmostEqual(motif.alignment.sequences[2].pvalue, 6.49e-19)
        self.assertAlmostEqual(motif.alignment.sequences[3].pvalue, 1.92e-18)
        self.assertAlmostEqual(motif.alignment.sequences[4].pvalue, 5.46e-18)
        self.assertAlmostEqual(motif.alignment.sequences[5].pvalue, 6.21e-18)
        self.assertAlmostEqual(motif.alignment.sequences[6].pvalue, 4.52e-17)
        self.assertAlmostEqual(motif.alignment.sequences[7].pvalue, 4.52e-17)
        self.assertAlmostEqual(motif.alignment.sequences[8].pvalue, 9.21e-17)
        self.assertAlmostEqual(motif.alignment.sequences[9].pvalue, 1.65e-16)
        self.assertAlmostEqual(motif.alignment.sequences[10].pvalue, 2.07e-16)
        self.assertAlmostEqual(motif.alignment.sequences[11].pvalue, 3.65e-16)
        self.assertAlmostEqual(motif.alignment.sequences[12].pvalue, 5.7e-16)
        self.assertAlmostEqual(motif.alignment.sequences[13].pvalue, 5.7e-16)
        self.assertAlmostEqual(motif.alignment.sequences[14].pvalue, 7.93e-16)
        self.assertAlmostEqual(motif.alignment.sequences[15].pvalue, 8.85e-16)
        self.assertAlmostEqual(motif.alignment.sequences[16].pvalue, 1.1e-15)
        self.assertAlmostEqual(motif.alignment.sequences[17].pvalue, 1.69e-15)
        self.assertAlmostEqual(motif.alignment.sequences[18].pvalue, 3.54e-15)
        self.assertAlmostEqual(motif.alignment.sequences[19].pvalue, 4.83e-15)
        self.assertAlmostEqual(motif.alignment.sequences[20].pvalue, 7.27e-15)
        self.assertAlmostEqual(motif.alignment.sequences[21].pvalue, 9.85e-15)
        self.assertAlmostEqual(motif.alignment.sequences[22].pvalue, 2.41e-14)
        self.assertAlmostEqual(motif.alignment.sequences[23].pvalue, 2.66e-14)
        self.assertAlmostEqual(motif.alignment.sequences[24].pvalue, 1.22e-13)
        self.assertAlmostEqual(motif.alignment.sequences[25].pvalue, 5.18e-13)
        self.assertAlmostEqual(motif.alignment.sequences[26].pvalue, 1.24e-12)
        self.assertAlmostEqual(motif.alignment.sequences[27].pvalue, 1.35e-12)
        self.assertAlmostEqual(motif.alignment.sequences[28].pvalue, 5.59e-12)
        self.assertAlmostEqual(motif.alignment.sequences[29].pvalue, 1.44e-10)
        self.assertAlmostEqual(motif.alignment.sequences[30].pvalue, 1.61e-08)
        self.assertAlmostEqual(motif.alignment.sequences[31].pvalue, 4.26e-08)
        self.assertAlmostEqual(motif.alignment.sequences[32].pvalue, 1.16e-07)
        self.assertEqual(motif.alignment.sequences[0].sequence_name, "BUDC_KLETE")
        self.assertEqual(motif.alignment.sequences[1].sequence_name, "NODG_RHIME")
        self.assertEqual(motif.alignment.sequences[2].sequence_name, "FVT1_HUMAN")
        self.assertEqual(motif.alignment.sequences[3].sequence_name, "DHES_HUMAN")
        self.assertEqual(motif.alignment.sequences[4].sequence_name, "DHB3_HUMAN")
        self.assertEqual(motif.alignment.sequences[5].sequence_name, "YRTP_BACSU")
        self.assertEqual(motif.alignment.sequences[6].sequence_name, "HMTR_LEIMA")
        self.assertEqual(motif.alignment.sequences[7].sequence_name, "HDE_CANTR")
        self.assertEqual(motif.alignment.sequences[8].sequence_name, "DHGB_BACME")
        self.assertEqual(motif.alignment.sequences[9].sequence_name, "GUTD_ECOLI")
        self.assertEqual(motif.alignment.sequences[10].sequence_name, "3BHD_COMTE")
        self.assertEqual(motif.alignment.sequences[11].sequence_name, "DHII_HUMAN")
        self.assertEqual(motif.alignment.sequences[12].sequence_name, "BPHB_PSEPS")
        self.assertEqual(motif.alignment.sequences[13].sequence_name, "AP27_MOUSE")
        self.assertEqual(motif.alignment.sequences[14].sequence_name, "BDH_HUMAN")
        self.assertEqual(motif.alignment.sequences[15].sequence_name, "YINL_LISMO")
        self.assertEqual(motif.alignment.sequences[16].sequence_name, "FIXR_BRAJA")
        self.assertEqual(motif.alignment.sequences[17].sequence_name, "2BHD_STREX")
        self.assertEqual(motif.alignment.sequences[18].sequence_name, "RFBB_NEIGO")
        self.assertEqual(motif.alignment.sequences[19].sequence_name, "YURA_MYXXA")
        self.assertEqual(motif.alignment.sequences[20].sequence_name, "RIDH_KLEAE")
        self.assertEqual(motif.alignment.sequences[21].sequence_name, "DHMA_FLAS1")
        self.assertEqual(motif.alignment.sequences[22].sequence_name, "DHB2_HUMAN")
        self.assertEqual(motif.alignment.sequences[23].sequence_name, "HDHA_ECOLI")
        self.assertEqual(motif.alignment.sequences[24].sequence_name, "ENTA_ECOLI")
        self.assertEqual(motif.alignment.sequences[25].sequence_name, "LIGD_PSEPA")
        self.assertEqual(motif.alignment.sequences[26].sequence_name, "CSGA_MYXXA")
        self.assertEqual(motif.alignment.sequences[27].sequence_name, "BA72_EUBSP")
        self.assertEqual(motif.alignment.sequences[28].sequence_name, "ADH_DROME")
        self.assertEqual(motif.alignment.sequences[29].sequence_name, "MAS1_AGRRA")
        self.assertEqual(motif.alignment.sequences[30].sequence_name, "PCR_PEA")
        self.assertEqual(motif.alignment.sequences[31].sequence_name, "FABI_ECOLI")
        self.assertEqual(motif.alignment.sequences[32].sequence_name, "DHCA_HUMAN")
        self.assertEqual(motif.alignment.sequences[0].sequence_id, "sequence_7")
        self.assertEqual(motif.alignment.sequences[1].sequence_id, "sequence_18")
        self.assertEqual(motif.alignment.sequences[2].sequence_id, "sequence_27")
        self.assertEqual(motif.alignment.sequences[3].sequence_id, "sequence_8")
        self.assertEqual(motif.alignment.sequences[4].sequence_id, "sequence_24")
        self.assertEqual(motif.alignment.sequences[5].sequence_id, "sequence_21")
        self.assertEqual(motif.alignment.sequences[6].sequence_id, "sequence_28")
        self.assertEqual(motif.alignment.sequences[7].sequence_id, "sequence_15")
        self.assertEqual(motif.alignment.sequences[8].sequence_id, "sequence_9")
        self.assertEqual(motif.alignment.sequences[9].sequence_id, "sequence_14")
        self.assertEqual(motif.alignment.sequences[10].sequence_id, "sequence_1")
        self.assertEqual(motif.alignment.sequences[11].sequence_id, "sequence_10")
        self.assertEqual(motif.alignment.sequences[12].sequence_id, "sequence_6")
        self.assertEqual(motif.alignment.sequences[13].sequence_id, "sequence_3")
        self.assertEqual(motif.alignment.sequences[14].sequence_id, "sequence_5")
        self.assertEqual(motif.alignment.sequences[15].sequence_id, "sequence_20")
        self.assertEqual(motif.alignment.sequences[16].sequence_id, "sequence_13")
        self.assertEqual(motif.alignment.sequences[17].sequence_id, "sequence_0")
        self.assertEqual(motif.alignment.sequences[18].sequence_id, "sequence_31")
        self.assertEqual(motif.alignment.sequences[19].sequence_id, "sequence_32")
        self.assertEqual(motif.alignment.sequences[20].sequence_id, "sequence_19")
        self.assertEqual(motif.alignment.sequences[21].sequence_id, "sequence_11")
        self.assertEqual(motif.alignment.sequences[22].sequence_id, "sequence_23")
        self.assertEqual(motif.alignment.sequences[23].sequence_id, "sequence_16")
        self.assertEqual(motif.alignment.sequences[24].sequence_id, "sequence_12")
        self.assertEqual(motif.alignment.sequences[25].sequence_id, "sequence_17")
        self.assertEqual(motif.alignment.sequences[26].sequence_id, "sequence_22")
        self.assertEqual(motif.alignment.sequences[27].sequence_id, "sequence_4")
        self.assertEqual(motif.alignment.sequences[28].sequence_id, "sequence_2")
        self.assertEqual(motif.alignment.sequences[29].sequence_id, "sequence_29")
        self.assertEqual(motif.alignment.sequences[30].sequence_id, "sequence_30")
        self.assertEqual(motif.alignment.sequences[31].sequence_id, "sequence_26")
        self.assertEqual(motif.alignment.sequences[32].sequence_id, "sequence_25")
        self.assertEqual(motif.alignment.sequences[0].start, 144)
        self.assertEqual(motif.alignment.sequences[1].start, 144)
        self.assertEqual(motif.alignment.sequences[2].start, 178)
        self.assertEqual(motif.alignment.sequences[3].start, 147)
        self.assertEqual(motif.alignment.sequences[4].start, 190)
        self.assertEqual(motif.alignment.sequences[5].start, 147)
        self.assertEqual(motif.alignment.sequences[6].start, 185)
        self.assertEqual(motif.alignment.sequences[7].start, 459)
        self.assertEqual(motif.alignment.sequences[8].start, 152)
        self.assertEqual(motif.alignment.sequences[9].start, 146)
        self.assertEqual(motif.alignment.sequences[10].start, 143)
        self.assertEqual(motif.alignment.sequences[11].start, 175)
        self.assertEqual(motif.alignment.sequences[12].start, 145)
        self.assertEqual(motif.alignment.sequences[13].start, 141)
        self.assertEqual(motif.alignment.sequences[14].start, 200)
        self.assertEqual(motif.alignment.sequences[15].start, 146)
        self.assertEqual(motif.alignment.sequences[16].start, 181)
        self.assertEqual(motif.alignment.sequences[17].start, 144)
        self.assertEqual(motif.alignment.sequences[18].start, 157)
        self.assertEqual(motif.alignment.sequences[19].start, 152)
        self.assertEqual(motif.alignment.sequences[20].start, 152)
        self.assertEqual(motif.alignment.sequences[21].start, 157)
        self.assertEqual(motif.alignment.sequences[22].start, 224)
        self.assertEqual(motif.alignment.sequences[23].start, 151)
        self.assertEqual(motif.alignment.sequences[24].start, 136)
        self.assertEqual(motif.alignment.sequences[25].start, 149)
        self.assertEqual(motif.alignment.sequences[26].start, 80)
        self.assertEqual(motif.alignment.sequences[27].start, 149)
        self.assertEqual(motif.alignment.sequences[28].start, 144)
        self.assertEqual(motif.alignment.sequences[29].start, 384)
        self.assertEqual(motif.alignment.sequences[30].start, 18)
        self.assertEqual(motif.alignment.sequences[31].start, 177)
        self.assertEqual(motif.alignment.sequences[32].start, 144)
        self.assertEqual(motif.alignment.sequences[0], "VGNPELAVYSSSKFAVRGLTQTAARDLAP")
        self.assertEqual(motif.alignment.sequences[1], "IGNPGQTNYCASKAGMIGFSKSLAQEIAT")
        self.assertEqual(motif.alignment.sequences[2], "LGLFGFTAYSASKFAIRGLAEALQMEVKP")
        self.assertEqual(motif.alignment.sequences[3], "MGLPFNDVYCASKFALEGLCESLAVLLLP")
        self.assertEqual(motif.alignment.sequences[4], "FPWPLYSMYSASKAFVCAFSKALQEEYKA")
        self.assertEqual(motif.alignment.sequences[5], "RGAAVTSAYSASKFAVLGLTESLMQEVRK")
        self.assertEqual(motif.alignment.sequences[6], "QPLLGYTIYTMAKGALEGLTRSAALELAP")
        self.assertEqual(motif.alignment.sequences[7], "YGNFGQANYSSSKAGILGLSKTMAIEGAK")
        self.assertEqual(motif.alignment.sequences[8], "IPWPLFVHYAASKGGMKLMTETLALEYAP")
        self.assertEqual(motif.alignment.sequences[9], "VGSKHNSGYSAAKFGGVGLTQSLALDLAE")
        self.assertEqual(motif.alignment.sequences[10], "LPIEQYAGYSASKAAVSALTRAAALSCRK")
        self.assertEqual(motif.alignment.sequences[11], "VAYPMVAAYSASKFALDGFFSSIRKEYSV")
        self.assertEqual(motif.alignment.sequences[12], "YPNGGGPLYTAAKQAIVGLVRELAFELAP")
        self.assertEqual(motif.alignment.sequences[13], "VTFPNLITYSSTKGAMTMLTKAMAMELGP")
        self.assertEqual(motif.alignment.sequences[14], "MANPARSPYCITKFGVEAFSDCLRYEMYP")
        self.assertEqual(motif.alignment.sequences[15], "KAYPGGAVYGATKWAVRDLMEVLRMESAQ")
        self.assertEqual(motif.alignment.sequences[16], "VHPFAGSAYATSKAALASLTRELAHDYAP")
        self.assertEqual(motif.alignment.sequences[17], "MGLALTSSYGASKWGVRGLSKLAAVELGT")
        self.assertEqual(motif.alignment.sequences[18], "TPYAPSSPYSASKAAADHLVRAWQRTYRL")
        self.assertEqual(motif.alignment.sequences[19], "FRGLPATRYSASKAFLSTFMESLRVDLRG")
        self.assertEqual(motif.alignment.sequences[20], "VPVIWEPVYTASKFAVQAFVHTTRRQVAQ")
        self.assertEqual(motif.alignment.sequences[21], "MAEPEAAAYVAAKGGVAMLTRAMAVDLAR")
        self.assertEqual(motif.alignment.sequences[22], "APMERLASYGSSKAAVTMFSSVMRLELSK")
        self.assertEqual(motif.alignment.sequences[23], "NKNINMTSYASSKAAASHLVRNMAFDLGE")
        self.assertEqual(motif.alignment.sequences[24], "TPRIGMSAYGASKAALKSLALSVGLELAG")
        self.assertEqual(motif.alignment.sequences[25], "MGSALAGPYSAAKAASINLMEGYRQGLEK")
        self.assertEqual(motif.alignment.sequences[26], "NTDGGAYAYRMSKAALNMAVRSMSTDLRP")
        self.assertEqual(motif.alignment.sequences[27], "FGSLSGVGYPASKASVIGLTHGLGREIIR")
        self.assertEqual(motif.alignment.sequences[28], "NAIYQVPVYSGTKAAVVNFTSSLAKLAPI")
        self.assertEqual(motif.alignment.sequences[29], "RVLNPLVGYNMTKHALGGLTKTTQHVGWD")
        self.assertEqual(motif.alignment.sequences[30], "EGKIGASLKDSTLFGVSSLSDSLKGDFTS")
        self.assertEqual(motif.alignment.sequences[31], "MGPEGVRVNAISAGPIRTLAASGIKDFRK")
        self.assertEqual(motif.alignment.sequences[32], "RALKSCSPELQQKFRSETITEEELVGLMN")
        self.assertEqual(motif.consensus, "MGLPGASAYSASKAAVRGLTESLALELAP")
        self.assertEqual(motif[-8:-2].consensus, "SLALEL")

    def test_meme_parser_3(self):
        """Parse motifs/meme.farntrans5.classic.anr.xml file."""
        with open("motifs/meme.farntrans5.classic.anr.xml") as stream:
            record = motifs.parse(stream, "meme")
        self.assertEqual(record.version, "5.0.1")
        self.assertEqual(record.datafile, "common/farntrans5.s")
        self.assertEqual(record.alphabet, "ACDEFGHIKLMNPQRSTVWY")
        self.assertEqual(len(record.sequences), 5)
        self.assertEqual(record.sequences[0], "sequence_0")
        self.assertEqual(record.sequences[1], "sequence_1")
        self.assertEqual(record.sequences[2], "sequence_2")
        self.assertEqual(record.sequences[3], "sequence_3")
        self.assertEqual(record.sequences[4], "sequence_4")
        self.assertEqual(
            record.command,
            "meme common/farntrans5.s -oc results/meme15 -mod anr -protein -nmotifs 2 -objfun classic -minw 8 -nostatus ",
        )
        self.assertEqual(len(record), 2)
        motif = record[0]
        self.assertEqual(motif.name, "GGFGGRPGKEVDLCYTYCALAALAJLGSLD")
        self.assertEqual(record["GGFGGRPGKEVDLCYTYCALAALAJLGSLD"], motif)
        self.assertEqual(motif.num_occurrences, 24)
        self.assertAlmostEqual(motif.evalue, 2.2e-94)
        self.assertEqual(motif.alphabet, "ACDEFGHIKLMNPQRSTVWY")
        # using the old instances property:
        with self.assertWarns(BiopythonDeprecationWarning):
            self.assertEqual(len(motif.instances), 24)
            self.assertAlmostEqual(motif.instances[0].pvalue, 6.98e-22)
            self.assertAlmostEqual(motif.instances[1].pvalue, 4.67e-21)
            self.assertAlmostEqual(motif.instances[2].pvalue, 1.25e-19)
            self.assertAlmostEqual(motif.instances[3].pvalue, 1.56e-19)
            self.assertAlmostEqual(motif.instances[4].pvalue, 2.44e-19)
            self.assertAlmostEqual(motif.instances[5].pvalue, 6.47e-19)
            self.assertAlmostEqual(motif.instances[6].pvalue, 8.9e-19)
            self.assertAlmostEqual(motif.instances[7].pvalue, 2.53e-18)
            self.assertAlmostEqual(motif.instances[8].pvalue, 1.27e-17)
            self.assertAlmostEqual(motif.instances[9].pvalue, 2.77e-17)
            self.assertAlmostEqual(motif.instances[10].pvalue, 4.93e-17)
            self.assertAlmostEqual(motif.instances[11].pvalue, 7.19e-17)
            self.assertAlmostEqual(motif.instances[12].pvalue, 8.68e-17)
            self.assertAlmostEqual(motif.instances[13].pvalue, 2.62e-16)
            self.assertAlmostEqual(motif.instances[14].pvalue, 2.87e-16)
            self.assertAlmostEqual(motif.instances[15].pvalue, 7.66e-15)
            self.assertAlmostEqual(motif.instances[16].pvalue, 2.21e-14)
            self.assertAlmostEqual(motif.instances[17].pvalue, 3.29e-14)
            self.assertAlmostEqual(motif.instances[18].pvalue, 7.21e-14)
            self.assertAlmostEqual(motif.instances[19].pvalue, 1.14e-13)
            self.assertAlmostEqual(motif.instances[20].pvalue, 1.67e-13)
            self.assertAlmostEqual(motif.instances[21].pvalue, 4.42e-13)
            self.assertAlmostEqual(motif.instances[22].pvalue, 5.11e-13)
            self.assertAlmostEqual(motif.instances[23].pvalue, 2.82e-10)
            self.assertEqual(motif.instances[0].sequence_name, "BET2_YEAST")
            self.assertEqual(motif.instances[1].sequence_name, "RATRABGERB")
            self.assertEqual(motif.instances[2].sequence_name, "CAL1_YEAST")
            self.assertEqual(motif.instances[3].sequence_name, "PFTB_RAT")
            self.assertEqual(motif.instances[4].sequence_name, "PFTB_RAT")
            self.assertEqual(motif.instances[5].sequence_name, "RATRABGERB")
            self.assertEqual(motif.instances[6].sequence_name, "RATRABGERB")
            self.assertEqual(motif.instances[7].sequence_name, "BET2_YEAST")
            self.assertEqual(motif.instances[8].sequence_name, "RATRABGERB")
            self.assertEqual(motif.instances[9].sequence_name, "BET2_YEAST")
            self.assertEqual(motif.instances[10].sequence_name, "RAM1_YEAST")
            self.assertEqual(motif.instances[11].sequence_name, "BET2_YEAST")
            self.assertEqual(motif.instances[12].sequence_name, "RAM1_YEAST")
            self.assertEqual(motif.instances[13].sequence_name, "PFTB_RAT")
            self.assertEqual(motif.instances[14].sequence_name, "RAM1_YEAST")
            self.assertEqual(motif.instances[15].sequence_name, "PFTB_RAT")
            self.assertEqual(motif.instances[16].sequence_name, "RATRABGERB")
            self.assertEqual(motif.instances[17].sequence_name, "PFTB_RAT")
            self.assertEqual(motif.instances[18].sequence_name, "BET2_YEAST")
            self.assertEqual(motif.instances[19].sequence_name, "CAL1_YEAST")
            self.assertEqual(motif.instances[20].sequence_name, "RAM1_YEAST")
            self.assertEqual(motif.instances[21].sequence_name, "CAL1_YEAST")
            self.assertEqual(motif.instances[22].sequence_name, "RAM1_YEAST")
            self.assertEqual(motif.instances[23].sequence_name, "BET2_YEAST")
            self.assertEqual(motif.instances[0].sequence_id, "sequence_2")
            self.assertEqual(motif.instances[1].sequence_id, "sequence_3")
            self.assertEqual(motif.instances[2].sequence_id, "sequence_4")
            self.assertEqual(motif.instances[3].sequence_id, "sequence_1")
            self.assertEqual(motif.instances[4].sequence_id, "sequence_1")
            self.assertEqual(motif.instances[5].sequence_id, "sequence_3")
            self.assertEqual(motif.instances[6].sequence_id, "sequence_3")
            self.assertEqual(motif.instances[7].sequence_id, "sequence_2")
            self.assertEqual(motif.instances[8].sequence_id, "sequence_3")
            self.assertEqual(motif.instances[9].sequence_id, "sequence_2")
            self.assertEqual(motif.instances[10].sequence_id, "sequence_0")
            self.assertEqual(motif.instances[11].sequence_id, "sequence_2")
            self.assertEqual(motif.instances[12].sequence_id, "sequence_0")
            self.assertEqual(motif.instances[13].sequence_id, "sequence_1")
            self.assertEqual(motif.instances[14].sequence_id, "sequence_0")
            self.assertEqual(motif.instances[15].sequence_id, "sequence_1")
            self.assertEqual(motif.instances[16].sequence_id, "sequence_3")
            self.assertEqual(motif.instances[17].sequence_id, "sequence_1")
            self.assertEqual(motif.instances[18].sequence_id, "sequence_2")
            self.assertEqual(motif.instances[19].sequence_id, "sequence_4")
            self.assertEqual(motif.instances[20].sequence_id, "sequence_0")
            self.assertEqual(motif.instances[21].sequence_id, "sequence_4")
            self.assertEqual(motif.instances[22].sequence_id, "sequence_0")
            self.assertEqual(motif.instances[23].sequence_id, "sequence_2")
            self.assertEqual(motif.instances[0].strand, "+")
            self.assertEqual(motif.instances[1].strand, "+")
            self.assertEqual(motif.instances[2].strand, "+")
            self.assertEqual(motif.instances[3].strand, "+")
            self.assertEqual(motif.instances[4].strand, "+")
            self.assertEqual(motif.instances[5].strand, "+")
            self.assertEqual(motif.instances[6].strand, "+")
            self.assertEqual(motif.instances[7].strand, "+")
            self.assertEqual(motif.instances[8].strand, "+")
            self.assertEqual(motif.instances[9].strand, "+")
            self.assertEqual(motif.instances[10].strand, "+")
            self.assertEqual(motif.instances[11].strand, "+")
            self.assertEqual(motif.instances[12].strand, "+")
            self.assertEqual(motif.instances[13].strand, "+")
            self.assertEqual(motif.instances[14].strand, "+")
            self.assertEqual(motif.instances[15].strand, "+")
            self.assertEqual(motif.instances[16].strand, "+")
            self.assertEqual(motif.instances[17].strand, "+")
            self.assertEqual(motif.instances[18].strand, "+")
            self.assertEqual(motif.instances[19].strand, "+")
            self.assertEqual(motif.instances[20].strand, "+")
            self.assertEqual(motif.instances[21].strand, "+")
            self.assertEqual(motif.instances[22].strand, "+")
            self.assertEqual(motif.instances[23].strand, "+")
            self.assertEqual(motif.instances[0].length, 30)
            self.assertEqual(motif.instances[1].length, 30)
            self.assertEqual(motif.instances[2].length, 30)
            self.assertEqual(motif.instances[3].length, 30)
            self.assertEqual(motif.instances[4].length, 30)
            self.assertEqual(motif.instances[5].length, 30)
            self.assertEqual(motif.instances[6].length, 30)
            self.assertEqual(motif.instances[7].length, 30)
            self.assertEqual(motif.instances[8].length, 30)
            self.assertEqual(motif.instances[9].length, 30)
            self.assertEqual(motif.instances[10].length, 30)
            self.assertEqual(motif.instances[11].length, 30)
            self.assertEqual(motif.instances[12].length, 30)
            self.assertEqual(motif.instances[13].length, 30)
            self.assertEqual(motif.instances[14].length, 30)
            self.assertEqual(motif.instances[15].length, 30)
            self.assertEqual(motif.instances[16].length, 30)
            self.assertEqual(motif.instances[17].length, 30)
            self.assertEqual(motif.instances[18].length, 30)
            self.assertEqual(motif.instances[19].length, 30)
            self.assertEqual(motif.instances[20].length, 30)
            self.assertEqual(motif.instances[21].length, 30)
            self.assertEqual(motif.instances[22].length, 30)
            self.assertEqual(motif.instances[23].length, 30)
            self.assertEqual(motif.instances[0].start, 223)
            self.assertEqual(motif.instances[1].start, 227)
            self.assertEqual(motif.instances[2].start, 275)
            self.assertEqual(motif.instances[3].start, 237)
            self.assertEqual(motif.instances[4].start, 138)
            self.assertEqual(motif.instances[5].start, 179)
            self.assertEqual(motif.instances[6].start, 131)
            self.assertEqual(motif.instances[7].start, 172)
            self.assertEqual(motif.instances[8].start, 276)
            self.assertEqual(motif.instances[9].start, 124)
            self.assertEqual(motif.instances[10].start, 247)
            self.assertEqual(motif.instances[11].start, 272)
            self.assertEqual(motif.instances[12].start, 145)
            self.assertEqual(motif.instances[13].start, 286)
            self.assertEqual(motif.instances[14].start, 296)
            self.assertEqual(motif.instances[15].start, 348)
            self.assertEqual(motif.instances[16].start, 83)
            self.assertEqual(motif.instances[17].start, 189)
            self.assertEqual(motif.instances[18].start, 73)
            self.assertEqual(motif.instances[19].start, 205)
            self.assertEqual(motif.instances[20].start, 198)
            self.assertEqual(motif.instances[21].start, 327)
            self.assertEqual(motif.instances[22].start, 349)
            self.assertEqual(motif.instances[23].start, 24)
            self.assertEqual(motif.instances[0], "GGLNGRPSKLPDVCYSWWVLSSLAIIGRLD")
            self.assertEqual(motif.instances[1], "GGLNGRPEKLPDVCYSWWVLASLKIIGRLH")
            self.assertEqual(motif.instances[2], "GGFQGRENKFADTCYAFWCLNSLHLLTKDW")
            self.assertEqual(motif.instances[3], "GGIGGVPGMEAHGGYTFCGLAALVILKKER")
            self.assertEqual(motif.instances[4], "GGFGGGPGQYPHLAPTYAAVNALCIIGTEE")
            self.assertEqual(motif.instances[5], "GGFGCRPGSESHAGQIYCCTGFLAITSQLH")
            self.assertEqual(motif.instances[6], "GSFAGDIWGEIDTRFSFCAVATLALLGKLD")
            self.assertEqual(motif.instances[7], "GGFGLCPNAESHAAQAFTCLGALAIANKLD")
            self.assertEqual(motif.instances[8], "GGFADRPGDMVDPFHTLFGIAGLSLLGEEQ")
            self.assertEqual(motif.instances[9], "GSFQGDRFGEVDTRFVYTALSALSILGELT")
            self.assertEqual(motif.instances[10], "GFGSCPHVDEAHGGYTFCATASLAILRSMD")
            self.assertEqual(motif.instances[11], "GGISDRPENEVDVFHTVFGVAGLSLMGYDN")
            self.assertEqual(motif.instances[12], "GPFGGGPGQLSHLASTYAAINALSLCDNID")
            self.assertEqual(motif.instances[13], "GGFQGRCNKLVDGCYSFWQAGLLPLLHRAL")
            self.assertEqual(motif.instances[14], "RGFCGRSNKLVDGCYSFWVGGSAAILEAFG")
            self.assertEqual(motif.instances[15], "GGLLDKPGKSRDFYHTCYCLSGLSIAQHFG")
            self.assertEqual(motif.instances[16], "GGVSASIGHDPHLLYTLSAVQILTLYDSIH")
            self.assertEqual(motif.instances[17], "GSFLMHVGGEVDVRSAYCAASVASLTNIIT")
            self.assertEqual(motif.instances[18], "GAFAPFPRHDAHLLTTLSAVQILATYDALD")
            self.assertEqual(motif.instances[19], "YNGAFGAHNEPHSGYTSCALSTLALLSSLE")
            self.assertEqual(motif.instances[20], "GFKTCLEVGEVDTRGIYCALSIATLLNILT")
            self.assertEqual(motif.instances[21], "GGFSKNDEEDADLYHSCLGSAALALIEGKF")
            self.assertEqual(motif.instances[22], "PGLRDKPGAHSDFYHTNYCLLGLAVAESSY")
            self.assertEqual(motif.instances[23], "HNFEYWLTEHLRLNGIYWGLTALCVLDSPE")
        self.assertEqual(len(motif.alignment.sequences), 24)
        self.assertAlmostEqual(motif.alignment.sequences[0].pvalue, 6.98e-22)
        self.assertAlmostEqual(motif.alignment.sequences[1].pvalue, 4.67e-21)
        self.assertAlmostEqual(motif.alignment.sequences[2].pvalue, 1.25e-19)
        self.assertAlmostEqual(motif.alignment.sequences[3].pvalue, 1.56e-19)
        self.assertAlmostEqual(motif.alignment.sequences[4].pvalue, 2.44e-19)
        self.assertAlmostEqual(motif.alignment.sequences[5].pvalue, 6.47e-19)
        self.assertAlmostEqual(motif.alignment.sequences[6].pvalue, 8.9e-19)
        self.assertAlmostEqual(motif.alignment.sequences[7].pvalue, 2.53e-18)
        self.assertAlmostEqual(motif.alignment.sequences[8].pvalue, 1.27e-17)
        self.assertAlmostEqual(motif.alignment.sequences[9].pvalue, 2.77e-17)
        self.assertAlmostEqual(motif.alignment.sequences[10].pvalue, 4.93e-17)
        self.assertAlmostEqual(motif.alignment.sequences[11].pvalue, 7.19e-17)
        self.assertAlmostEqual(motif.alignment.sequences[12].pvalue, 8.68e-17)
        self.assertAlmostEqual(motif.alignment.sequences[13].pvalue, 2.62e-16)
        self.assertAlmostEqual(motif.alignment.sequences[14].pvalue, 2.87e-16)
        self.assertAlmostEqual(motif.alignment.sequences[15].pvalue, 7.66e-15)
        self.assertAlmostEqual(motif.alignment.sequences[16].pvalue, 2.21e-14)
        self.assertAlmostEqual(motif.alignment.sequences[17].pvalue, 3.29e-14)
        self.assertAlmostEqual(motif.alignment.sequences[18].pvalue, 7.21e-14)
        self.assertAlmostEqual(motif.alignment.sequences[19].pvalue, 1.14e-13)
        self.assertAlmostEqual(motif.alignment.sequences[20].pvalue, 1.67e-13)
        self.assertAlmostEqual(motif.alignment.sequences[21].pvalue, 4.42e-13)
        self.assertAlmostEqual(motif.alignment.sequences[22].pvalue, 5.11e-13)
        self.assertAlmostEqual(motif.alignment.sequences[23].pvalue, 2.82e-10)
        self.assertEqual(motif.alignment.sequences[0].sequence_name, "BET2_YEAST")
        self.assertEqual(motif.alignment.sequences[1].sequence_name, "RATRABGERB")
        self.assertEqual(motif.alignment.sequences[2].sequence_name, "CAL1_YEAST")
        self.assertEqual(motif.alignment.sequences[3].sequence_name, "PFTB_RAT")
        self.assertEqual(motif.alignment.sequences[4].sequence_name, "PFTB_RAT")
        self.assertEqual(motif.alignment.sequences[5].sequence_name, "RATRABGERB")
        self.assertEqual(motif.alignment.sequences[6].sequence_name, "RATRABGERB")
        self.assertEqual(motif.alignment.sequences[7].sequence_name, "BET2_YEAST")
        self.assertEqual(motif.alignment.sequences[8].sequence_name, "RATRABGERB")
        self.assertEqual(motif.alignment.sequences[9].sequence_name, "BET2_YEAST")
        self.assertEqual(motif.alignment.sequences[10].sequence_name, "RAM1_YEAST")
        self.assertEqual(motif.alignment.sequences[11].sequence_name, "BET2_YEAST")
        self.assertEqual(motif.alignment.sequences[12].sequence_name, "RAM1_YEAST")
        self.assertEqual(motif.alignment.sequences[13].sequence_name, "PFTB_RAT")
        self.assertEqual(motif.alignment.sequences[14].sequence_name, "RAM1_YEAST")
        self.assertEqual(motif.alignment.sequences[15].sequence_name, "PFTB_RAT")
        self.assertEqual(motif.alignment.sequences[16].sequence_name, "RATRABGERB")
        self.assertEqual(motif.alignment.sequences[17].sequence_name, "PFTB_RAT")
        self.assertEqual(motif.alignment.sequences[18].sequence_name, "BET2_YEAST")
        self.assertEqual(motif.alignment.sequences[19].sequence_name, "CAL1_YEAST")
        self.assertEqual(motif.alignment.sequences[20].sequence_name, "RAM1_YEAST")
        self.assertEqual(motif.alignment.sequences[21].sequence_name, "CAL1_YEAST")
        self.assertEqual(motif.alignment.sequences[22].sequence_name, "RAM1_YEAST")
        self.assertEqual(motif.alignment.sequences[23].sequence_name, "BET2_YEAST")
        self.assertEqual(motif.alignment.sequences[0].sequence_id, "sequence_2")
        self.assertEqual(motif.alignment.sequences[1].sequence_id, "sequence_3")
        self.assertEqual(motif.alignment.sequences[2].sequence_id, "sequence_4")
        self.assertEqual(motif.alignment.sequences[3].sequence_id, "sequence_1")
        self.assertEqual(motif.alignment.sequences[4].sequence_id, "sequence_1")
        self.assertEqual(motif.alignment.sequences[5].sequence_id, "sequence_3")
        self.assertEqual(motif.alignment.sequences[6].sequence_id, "sequence_3")
        self.assertEqual(motif.alignment.sequences[7].sequence_id, "sequence_2")
        self.assertEqual(motif.alignment.sequences[8].sequence_id, "sequence_3")
        self.assertEqual(motif.alignment.sequences[9].sequence_id, "sequence_2")
        self.assertEqual(motif.alignment.sequences[10].sequence_id, "sequence_0")
        self.assertEqual(motif.alignment.sequences[11].sequence_id, "sequence_2")
        self.assertEqual(motif.alignment.sequences[12].sequence_id, "sequence_0")
        self.assertEqual(motif.alignment.sequences[13].sequence_id, "sequence_1")
        self.assertEqual(motif.alignment.sequences[14].sequence_id, "sequence_0")
        self.assertEqual(motif.alignment.sequences[15].sequence_id, "sequence_1")
        self.assertEqual(motif.alignment.sequences[16].sequence_id, "sequence_3")
        self.assertEqual(motif.alignment.sequences[17].sequence_id, "sequence_1")
        self.assertEqual(motif.alignment.sequences[18].sequence_id, "sequence_2")
        self.assertEqual(motif.alignment.sequences[19].sequence_id, "sequence_4")
        self.assertEqual(motif.alignment.sequences[20].sequence_id, "sequence_0")
        self.assertEqual(motif.alignment.sequences[21].sequence_id, "sequence_4")
        self.assertEqual(motif.alignment.sequences[22].sequence_id, "sequence_0")
        self.assertEqual(motif.alignment.sequences[23].sequence_id, "sequence_2")
        self.assertEqual(motif.alignment.sequences[0].strand, "+")
        self.assertEqual(motif.alignment.sequences[1].strand, "+")
        self.assertEqual(motif.alignment.sequences[2].strand, "+")
        self.assertEqual(motif.alignment.sequences[3].strand, "+")
        self.assertEqual(motif.alignment.sequences[4].strand, "+")
        self.assertEqual(motif.alignment.sequences[5].strand, "+")
        self.assertEqual(motif.alignment.sequences[6].strand, "+")
        self.assertEqual(motif.alignment.sequences[7].strand, "+")
        self.assertEqual(motif.alignment.sequences[8].strand, "+")
        self.assertEqual(motif.alignment.sequences[9].strand, "+")
        self.assertEqual(motif.alignment.sequences[10].strand, "+")
        self.assertEqual(motif.alignment.sequences[11].strand, "+")
        self.assertEqual(motif.alignment.sequences[12].strand, "+")
        self.assertEqual(motif.alignment.sequences[13].strand, "+")
        self.assertEqual(motif.alignment.sequences[14].strand, "+")
        self.assertEqual(motif.alignment.sequences[15].strand, "+")
        self.assertEqual(motif.alignment.sequences[16].strand, "+")
        self.assertEqual(motif.alignment.sequences[17].strand, "+")
        self.assertEqual(motif.alignment.sequences[18].strand, "+")
        self.assertEqual(motif.alignment.sequences[19].strand, "+")
        self.assertEqual(motif.alignment.sequences[20].strand, "+")
        self.assertEqual(motif.alignment.sequences[21].strand, "+")
        self.assertEqual(motif.alignment.sequences[22].strand, "+")
        self.assertEqual(motif.alignment.sequences[23].strand, "+")
        self.assertEqual(motif.alignment.sequences[0].length, 30)
        self.assertEqual(motif.alignment.sequences[1].length, 30)
        self.assertEqual(motif.alignment.sequences[2].length, 30)
        self.assertEqual(motif.alignment.sequences[3].length, 30)
        self.assertEqual(motif.alignment.sequences[4].length, 30)
        self.assertEqual(motif.alignment.sequences[5].length, 30)
        self.assertEqual(motif.alignment.sequences[6].length, 30)
        self.assertEqual(motif.alignment.sequences[7].length, 30)
        self.assertEqual(motif.alignment.sequences[8].length, 30)
        self.assertEqual(motif.alignment.sequences[9].length, 30)
        self.assertEqual(motif.alignment.sequences[10].length, 30)
        self.assertEqual(motif.alignment.sequences[11].length, 30)
        self.assertEqual(motif.alignment.sequences[12].length, 30)
        self.assertEqual(motif.alignment.sequences[13].length, 30)
        self.assertEqual(motif.alignment.sequences[14].length, 30)
        self.assertEqual(motif.alignment.sequences[15].length, 30)
        self.assertEqual(motif.alignment.sequences[16].length, 30)
        self.assertEqual(motif.alignment.sequences[17].length, 30)
        self.assertEqual(motif.alignment.sequences[18].length, 30)
        self.assertEqual(motif.alignment.sequences[19].length, 30)
        self.assertEqual(motif.alignment.sequences[20].length, 30)
        self.assertEqual(motif.alignment.sequences[21].length, 30)
        self.assertEqual(motif.alignment.sequences[22].length, 30)
        self.assertEqual(motif.alignment.sequences[23].length, 30)
        self.assertEqual(motif.alignment.sequences[0].start, 223)
        self.assertEqual(motif.alignment.sequences[1].start, 227)
        self.assertEqual(motif.alignment.sequences[2].start, 275)
        self.assertEqual(motif.alignment.sequences[3].start, 237)
        self.assertEqual(motif.alignment.sequences[4].start, 138)
        self.assertEqual(motif.alignment.sequences[5].start, 179)
        self.assertEqual(motif.alignment.sequences[6].start, 131)
        self.assertEqual(motif.alignment.sequences[7].start, 172)
        self.assertEqual(motif.alignment.sequences[8].start, 276)
        self.assertEqual(motif.alignment.sequences[9].start, 124)
        self.assertEqual(motif.alignment.sequences[10].start, 247)
        self.assertEqual(motif.alignment.sequences[11].start, 272)
        self.assertEqual(motif.alignment.sequences[12].start, 145)
        self.assertEqual(motif.alignment.sequences[13].start, 286)
        self.assertEqual(motif.alignment.sequences[14].start, 296)
        self.assertEqual(motif.alignment.sequences[15].start, 348)
        self.assertEqual(motif.alignment.sequences[16].start, 83)
        self.assertEqual(motif.alignment.sequences[17].start, 189)
        self.assertEqual(motif.alignment.sequences[18].start, 73)
        self.assertEqual(motif.alignment.sequences[19].start, 205)
        self.assertEqual(motif.alignment.sequences[20].start, 198)
        self.assertEqual(motif.alignment.sequences[21].start, 327)
        self.assertEqual(motif.alignment.sequences[22].start, 349)
        self.assertEqual(motif.alignment.sequences[23].start, 24)
        self.assertEqual(motif.alignment.sequences[0], "GGLNGRPSKLPDVCYSWWVLSSLAIIGRLD")
        self.assertEqual(motif.alignment.sequences[1], "GGLNGRPEKLPDVCYSWWVLASLKIIGRLH")
        self.assertEqual(motif.alignment.sequences[2], "GGFQGRENKFADTCYAFWCLNSLHLLTKDW")
        self.assertEqual(motif.alignment.sequences[3], "GGIGGVPGMEAHGGYTFCGLAALVILKKER")
        self.assertEqual(motif.alignment.sequences[4], "GGFGGGPGQYPHLAPTYAAVNALCIIGTEE")
        self.assertEqual(motif.alignment.sequences[5], "GGFGCRPGSESHAGQIYCCTGFLAITSQLH")
        self.assertEqual(motif.alignment.sequences[6], "GSFAGDIWGEIDTRFSFCAVATLALLGKLD")
        self.assertEqual(motif.alignment.sequences[7], "GGFGLCPNAESHAAQAFTCLGALAIANKLD")
        self.assertEqual(motif.alignment.sequences[8], "GGFADRPGDMVDPFHTLFGIAGLSLLGEEQ")
        self.assertEqual(motif.alignment.sequences[9], "GSFQGDRFGEVDTRFVYTALSALSILGELT")
        self.assertEqual(
            motif.alignment.sequences[10], "GFGSCPHVDEAHGGYTFCATASLAILRSMD"
        )
        self.assertEqual(
            motif.alignment.sequences[11], "GGISDRPENEVDVFHTVFGVAGLSLMGYDN"
        )
        self.assertEqual(
            motif.alignment.sequences[12], "GPFGGGPGQLSHLASTYAAINALSLCDNID"
        )
        self.assertEqual(
            motif.alignment.sequences[13], "GGFQGRCNKLVDGCYSFWQAGLLPLLHRAL"
        )
        self.assertEqual(
            motif.alignment.sequences[14], "RGFCGRSNKLVDGCYSFWVGGSAAILEAFG"
        )
        self.assertEqual(
            motif.alignment.sequences[15], "GGLLDKPGKSRDFYHTCYCLSGLSIAQHFG"
        )
        self.assertEqual(
            motif.alignment.sequences[16], "GGVSASIGHDPHLLYTLSAVQILTLYDSIH"
        )
        self.assertEqual(
            motif.alignment.sequences[17], "GSFLMHVGGEVDVRSAYCAASVASLTNIIT"
        )
        self.assertEqual(
            motif.alignment.sequences[18], "GAFAPFPRHDAHLLTTLSAVQILATYDALD"
        )
        self.assertEqual(
            motif.alignment.sequences[19], "YNGAFGAHNEPHSGYTSCALSTLALLSSLE"
        )
        self.assertEqual(
            motif.alignment.sequences[20], "GFKTCLEVGEVDTRGIYCALSIATLLNILT"
        )
        self.assertEqual(
            motif.alignment.sequences[21], "GGFSKNDEEDADLYHSCLGSAALALIEGKF"
        )
        self.assertEqual(
            motif.alignment.sequences[22], "PGLRDKPGAHSDFYHTNYCLLGLAVAESSY"
        )
        self.assertEqual(
            motif.alignment.sequences[23], "HNFEYWLTEHLRLNGIYWGLTALCVLDSPE"
        )
        self.assertEqual(motif.consensus, "GGFGGRPGKEVDLCYTFCALAALALLGSLD")
        self.assertEqual(motif[3:-8].consensus, "GGRPGKEVDLCYTFCALAA")
        motif = record[1]
        self.assertEqual(motif.name, "JNKEKLLEYILSCQ")
        self.assertEqual(record["JNKEKLLEYILSCQ"], motif)
        self.assertEqual(motif.num_occurrences, 21)
        self.assertAlmostEqual(motif.evalue, 3.1e-19)
        self.assertEqual(motif.alphabet, "ACDEFGHIKLMNPQRSTVWY")
        self.assertEqual(len(motif.alignment.sequences), 21)
        self.assertAlmostEqual(motif.alignment.sequences[0].pvalue, 2.71e-12)
        self.assertAlmostEqual(motif.alignment.sequences[1].pvalue, 5.7e-12)
        self.assertAlmostEqual(motif.alignment.sequences[2].pvalue, 6.43e-12)
        self.assertAlmostEqual(motif.alignment.sequences[3].pvalue, 2.61e-11)
        self.assertAlmostEqual(motif.alignment.sequences[4].pvalue, 6.3e-11)
        self.assertAlmostEqual(motif.alignment.sequences[5].pvalue, 2.7e-10)
        self.assertAlmostEqual(motif.alignment.sequences[6].pvalue, 4.03e-10)
        self.assertAlmostEqual(motif.alignment.sequences[7].pvalue, 1.27e-09)
        self.assertAlmostEqual(motif.alignment.sequences[8].pvalue, 3.17e-09)
        self.assertAlmostEqual(motif.alignment.sequences[9].pvalue, 6.39e-09)
        self.assertAlmostEqual(motif.alignment.sequences[10].pvalue, 6.96e-09)
        self.assertAlmostEqual(motif.alignment.sequences[11].pvalue, 1.06e-08)
        self.assertAlmostEqual(motif.alignment.sequences[12].pvalue, 1.26e-08)
        self.assertAlmostEqual(motif.alignment.sequences[13].pvalue, 1.37e-08)
        self.assertAlmostEqual(motif.alignment.sequences[14].pvalue, 2.07e-08)
        self.assertAlmostEqual(motif.alignment.sequences[15].pvalue, 4.96e-08)
        self.assertAlmostEqual(motif.alignment.sequences[16].pvalue, 1.15e-07)
        self.assertAlmostEqual(motif.alignment.sequences[17].pvalue, 1.44e-07)
        self.assertAlmostEqual(motif.alignment.sequences[18].pvalue, 1.55e-07)
        self.assertAlmostEqual(motif.alignment.sequences[19].pvalue, 1.93e-07)
        self.assertAlmostEqual(motif.alignment.sequences[20].pvalue, 5.2e-07)
        self.assertEqual(motif.alignment.sequences[0].sequence_name, "RATRABGERB")
        self.assertEqual(motif.alignment.sequences[1].sequence_name, "BET2_YEAST")
        self.assertEqual(motif.alignment.sequences[2].sequence_name, "RATRABGERB")
        self.assertEqual(motif.alignment.sequences[3].sequence_name, "RATRABGERB")
        self.assertEqual(motif.alignment.sequences[4].sequence_name, "CAL1_YEAST")
        self.assertEqual(motif.alignment.sequences[5].sequence_name, "RAM1_YEAST")
        self.assertEqual(motif.alignment.sequences[6].sequence_name, "PFTB_RAT")
        self.assertEqual(motif.alignment.sequences[7].sequence_name, "RATRABGERB")
        self.assertEqual(motif.alignment.sequences[8].sequence_name, "BET2_YEAST")
        self.assertEqual(motif.alignment.sequences[9].sequence_name, "PFTB_RAT")
        self.assertEqual(motif.alignment.sequences[10].sequence_name, "RAM1_YEAST")
        self.assertEqual(motif.alignment.sequences[11].sequence_name, "CAL1_YEAST")
        self.assertEqual(motif.alignment.sequences[12].sequence_name, "PFTB_RAT")
        self.assertEqual(motif.alignment.sequences[13].sequence_name, "BET2_YEAST")
        self.assertEqual(motif.alignment.sequences[14].sequence_name, "RAM1_YEAST")
        self.assertEqual(motif.alignment.sequences[15].sequence_name, "RAM1_YEAST")
        self.assertEqual(motif.alignment.sequences[16].sequence_name, "RATRABGERB")
        self.assertEqual(motif.alignment.sequences[17].sequence_name, "RAM1_YEAST")
        self.assertEqual(motif.alignment.sequences[18].sequence_name, "PFTB_RAT")
        self.assertEqual(motif.alignment.sequences[19].sequence_name, "BET2_YEAST")
        self.assertEqual(motif.alignment.sequences[20].sequence_name, "CAL1_YEAST")
        self.assertEqual(motif.alignment.sequences[0].sequence_id, "sequence_3")
        self.assertEqual(motif.alignment.sequences[1].sequence_id, "sequence_2")
        self.assertEqual(motif.alignment.sequences[2].sequence_id, "sequence_3")
        self.assertEqual(motif.alignment.sequences[3].sequence_id, "sequence_3")
        self.assertEqual(motif.alignment.sequences[4].sequence_id, "sequence_4")
        self.assertEqual(motif.alignment.sequences[5].sequence_id, "sequence_0")
        self.assertEqual(motif.alignment.sequences[6].sequence_id, "sequence_1")
        self.assertEqual(motif.alignment.sequences[7].sequence_id, "sequence_3")
        self.assertEqual(motif.alignment.sequences[8].sequence_id, "sequence_2")
        self.assertEqual(motif.alignment.sequences[9].sequence_id, "sequence_1")
        self.assertEqual(motif.alignment.sequences[10].sequence_id, "sequence_0")
        self.assertEqual(motif.alignment.sequences[11].sequence_id, "sequence_4")
        self.assertEqual(motif.alignment.sequences[12].sequence_id, "sequence_1")
        self.assertEqual(motif.alignment.sequences[13].sequence_id, "sequence_2")
        self.assertEqual(motif.alignment.sequences[14].sequence_id, "sequence_0")
        self.assertEqual(motif.alignment.sequences[15].sequence_id, "sequence_0")
        self.assertEqual(motif.alignment.sequences[16].sequence_id, "sequence_3")
        self.assertEqual(motif.alignment.sequences[17].sequence_id, "sequence_0")
        self.assertEqual(motif.alignment.sequences[18].sequence_id, "sequence_1")
        self.assertEqual(motif.alignment.sequences[19].sequence_id, "sequence_2")
        self.assertEqual(motif.alignment.sequences[20].sequence_id, "sequence_4")
        self.assertEqual(motif.alignment.sequences[0].strand, "+")
        self.assertEqual(motif.alignment.sequences[1].strand, "+")
        self.assertEqual(motif.alignment.sequences[2].strand, "+")
        self.assertEqual(motif.alignment.sequences[3].strand, "+")
        self.assertEqual(motif.alignment.sequences[4].strand, "+")
        self.assertEqual(motif.alignment.sequences[5].strand, "+")
        self.assertEqual(motif.alignment.sequences[6].strand, "+")
        self.assertEqual(motif.alignment.sequences[7].strand, "+")
        self.assertEqual(motif.alignment.sequences[8].strand, "+")
        self.assertEqual(motif.alignment.sequences[9].strand, "+")
        self.assertEqual(motif.alignment.sequences[10].strand, "+")
        self.assertEqual(motif.alignment.sequences[11].strand, "+")
        self.assertEqual(motif.alignment.sequences[12].strand, "+")
        self.assertEqual(motif.alignment.sequences[13].strand, "+")
        self.assertEqual(motif.alignment.sequences[14].strand, "+")
        self.assertEqual(motif.alignment.sequences[15].strand, "+")
        self.assertEqual(motif.alignment.sequences[16].strand, "+")
        self.assertEqual(motif.alignment.sequences[17].strand, "+")
        self.assertEqual(motif.alignment.sequences[18].strand, "+")
        self.assertEqual(motif.alignment.sequences[19].strand, "+")
        self.assertEqual(motif.alignment.sequences[20].strand, "+")
        self.assertEqual(motif.alignment.sequences[0].length, 14)
        self.assertEqual(motif.alignment.sequences[1].length, 14)
        self.assertEqual(motif.alignment.sequences[2].length, 14)
        self.assertEqual(motif.alignment.sequences[3].length, 14)
        self.assertEqual(motif.alignment.sequences[4].length, 14)
        self.assertEqual(motif.alignment.sequences[5].length, 14)
        self.assertEqual(motif.alignment.sequences[6].length, 14)
        self.assertEqual(motif.alignment.sequences[7].length, 14)
        self.assertEqual(motif.alignment.sequences[8].length, 14)
        self.assertEqual(motif.alignment.sequences[9].length, 14)
        self.assertEqual(motif.alignment.sequences[10].length, 14)
        self.assertEqual(motif.alignment.sequences[11].length, 14)
        self.assertEqual(motif.alignment.sequences[12].length, 14)
        self.assertEqual(motif.alignment.sequences[13].length, 14)
        self.assertEqual(motif.alignment.sequences[14].length, 14)
        self.assertEqual(motif.alignment.sequences[15].length, 14)
        self.assertEqual(motif.alignment.sequences[16].length, 14)
        self.assertEqual(motif.alignment.sequences[17].length, 14)
        self.assertEqual(motif.alignment.sequences[18].length, 14)
        self.assertEqual(motif.alignment.sequences[19].length, 14)
        self.assertEqual(motif.alignment.sequences[20].length, 14)
        self.assertEqual(motif.alignment.sequences[0].start, 66)
        self.assertEqual(motif.alignment.sequences[1].start, 254)
        self.assertEqual(motif.alignment.sequences[2].start, 258)
        self.assertEqual(motif.alignment.sequences[3].start, 162)
        self.assertEqual(motif.alignment.sequences[4].start, 190)
        self.assertEqual(motif.alignment.sequences[5].start, 278)
        self.assertEqual(motif.alignment.sequences[6].start, 172)
        self.assertEqual(motif.alignment.sequences[7].start, 114)
        self.assertEqual(motif.alignment.sequences[8].start, 7)
        self.assertEqual(motif.alignment.sequences[9].start, 268)
        self.assertEqual(motif.alignment.sequences[10].start, 414)
        self.assertEqual(motif.alignment.sequences[11].start, 126)
        self.assertEqual(motif.alignment.sequences[12].start, 220)
        self.assertEqual(motif.alignment.sequences[13].start, 55)
        self.assertEqual(motif.alignment.sequences[14].start, 229)
        self.assertEqual(motif.alignment.sequences[15].start, 330)
        self.assertEqual(motif.alignment.sequences[16].start, 18)
        self.assertEqual(motif.alignment.sequences[17].start, 180)
        self.assertEqual(motif.alignment.sequences[18].start, 73)
        self.assertEqual(motif.alignment.sequences[19].start, 107)
        self.assertEqual(motif.alignment.sequences[20].start, 36)
        self.assertEqual(motif.alignment.sequences[0], "MNKEEILVFIKSCQ")
        self.assertEqual(motif.alignment.sequences[1], "INYEKLTEFILKCQ")
        self.assertEqual(motif.alignment.sequences[2], "IDREKLRSFILACQ")
        self.assertEqual(motif.alignment.sequences[3], "INVEKAIEFVLSCM")
        self.assertEqual(motif.alignment.sequences[4], "IDTEKLLGYIMSQQ")
        self.assertEqual(motif.alignment.sequences[5], "INVEKLLEWSSARQ")
        self.assertEqual(motif.alignment.sequences[6], "INREKLLQYLYSLK")
        self.assertEqual(motif.alignment.sequences[7], "INVDKVVAYVQSLQ")
        self.assertEqual(motif.alignment.sequences[8], "LLKEKHIRYIESLD")
        self.assertEqual(motif.alignment.sequences[9], "LNLKSLLQWVTSRQ")
        self.assertEqual(motif.alignment.sequences[10], "ENVRKIIHYFKSNL")
        self.assertEqual(motif.alignment.sequences[11], "LDKRSLARFVSKCQ")
        self.assertEqual(motif.alignment.sequences[12], "DLFEGTAEWIARCQ")
        self.assertEqual(motif.alignment.sequences[13], "FVKEEVISFVLSCW")
        self.assertEqual(motif.alignment.sequences[14], "ELTEGVLNYLKNCQ")
        self.assertEqual(motif.alignment.sequences[15], "FNKHALRDYILYCC")
        self.assertEqual(motif.alignment.sequences[16], "LLLEKHADYIASYG")
        self.assertEqual(motif.alignment.sequences[17], "IDRKGIYQWLISLK")
        self.assertEqual(motif.alignment.sequences[18], "LQREKHFHYLKRGL")
        self.assertEqual(motif.alignment.sequences[19], "DRKVRLISFIRGNQ")
        self.assertEqual(motif.alignment.sequences[20], "VNRMAIIFYSISGL")
        # using the old instances property:
        with self.assertWarns(BiopythonDeprecationWarning):
            self.assertEqual(len(motif.instances), 21)
            self.assertAlmostEqual(motif.instances[0].pvalue, 2.71e-12)
            self.assertAlmostEqual(motif.instances[1].pvalue, 5.7e-12)
            self.assertAlmostEqual(motif.instances[2].pvalue, 6.43e-12)
            self.assertAlmostEqual(motif.instances[3].pvalue, 2.61e-11)
            self.assertAlmostEqual(motif.instances[4].pvalue, 6.3e-11)
            self.assertAlmostEqual(motif.instances[5].pvalue, 2.7e-10)
            self.assertAlmostEqual(motif.instances[6].pvalue, 4.03e-10)
            self.assertAlmostEqual(motif.instances[7].pvalue, 1.27e-09)
            self.assertAlmostEqual(motif.instances[8].pvalue, 3.17e-09)
            self.assertAlmostEqual(motif.instances[9].pvalue, 6.39e-09)
            self.assertAlmostEqual(motif.instances[10].pvalue, 6.96e-09)
            self.assertAlmostEqual(motif.instances[11].pvalue, 1.06e-08)
            self.assertAlmostEqual(motif.instances[12].pvalue, 1.26e-08)
            self.assertAlmostEqual(motif.instances[13].pvalue, 1.37e-08)
            self.assertAlmostEqual(motif.instances[14].pvalue, 2.07e-08)
            self.assertAlmostEqual(motif.instances[15].pvalue, 4.96e-08)
            self.assertAlmostEqual(motif.instances[16].pvalue, 1.15e-07)
            self.assertAlmostEqual(motif.instances[17].pvalue, 1.44e-07)
            self.assertAlmostEqual(motif.instances[18].pvalue, 1.55e-07)
            self.assertAlmostEqual(motif.instances[19].pvalue, 1.93e-07)
            self.assertAlmostEqual(motif.instances[20].pvalue, 5.2e-07)
            self.assertEqual(motif.instances[0].sequence_name, "RATRABGERB")
            self.assertEqual(motif.instances[1].sequence_name, "BET2_YEAST")
            self.assertEqual(motif.instances[2].sequence_name, "RATRABGERB")
            self.assertEqual(motif.instances[3].sequence_name, "RATRABGERB")
            self.assertEqual(motif.instances[4].sequence_name, "CAL1_YEAST")
            self.assertEqual(motif.instances[5].sequence_name, "RAM1_YEAST")
            self.assertEqual(motif.instances[6].sequence_name, "PFTB_RAT")
            self.assertEqual(motif.instances[7].sequence_name, "RATRABGERB")
            self.assertEqual(motif.instances[8].sequence_name, "BET2_YEAST")
            self.assertEqual(motif.instances[9].sequence_name, "PFTB_RAT")
            self.assertEqual(motif.instances[10].sequence_name, "RAM1_YEAST")
            self.assertEqual(motif.instances[11].sequence_name, "CAL1_YEAST")
            self.assertEqual(motif.instances[12].sequence_name, "PFTB_RAT")
            self.assertEqual(motif.instances[13].sequence_name, "BET2_YEAST")
            self.assertEqual(motif.instances[14].sequence_name, "RAM1_YEAST")
            self.assertEqual(motif.instances[15].sequence_name, "RAM1_YEAST")
            self.assertEqual(motif.instances[16].sequence_name, "RATRABGERB")
            self.assertEqual(motif.instances[17].sequence_name, "RAM1_YEAST")
            self.assertEqual(motif.instances[18].sequence_name, "PFTB_RAT")
            self.assertEqual(motif.instances[19].sequence_name, "BET2_YEAST")
            self.assertEqual(motif.instances[20].sequence_name, "CAL1_YEAST")
            self.assertEqual(motif.instances[0].sequence_id, "sequence_3")
            self.assertEqual(motif.instances[1].sequence_id, "sequence_2")
            self.assertEqual(motif.instances[2].sequence_id, "sequence_3")
            self.assertEqual(motif.instances[3].sequence_id, "sequence_3")
            self.assertEqual(motif.instances[4].sequence_id, "sequence_4")
            self.assertEqual(motif.instances[5].sequence_id, "sequence_0")
            self.assertEqual(motif.instances[6].sequence_id, "sequence_1")
            self.assertEqual(motif.instances[7].sequence_id, "sequence_3")
            self.assertEqual(motif.instances[8].sequence_id, "sequence_2")
            self.assertEqual(motif.instances[9].sequence_id, "sequence_1")
            self.assertEqual(motif.instances[10].sequence_id, "sequence_0")
            self.assertEqual(motif.instances[11].sequence_id, "sequence_4")
            self.assertEqual(motif.instances[12].sequence_id, "sequence_1")
            self.assertEqual(motif.instances[13].sequence_id, "sequence_2")
            self.assertEqual(motif.instances[14].sequence_id, "sequence_0")
            self.assertEqual(motif.instances[15].sequence_id, "sequence_0")
            self.assertEqual(motif.instances[16].sequence_id, "sequence_3")
            self.assertEqual(motif.instances[17].sequence_id, "sequence_0")
            self.assertEqual(motif.instances[18].sequence_id, "sequence_1")
            self.assertEqual(motif.instances[19].sequence_id, "sequence_2")
            self.assertEqual(motif.instances[20].sequence_id, "sequence_4")
            self.assertEqual(motif.instances[0].strand, "+")
            self.assertEqual(motif.instances[1].strand, "+")
            self.assertEqual(motif.instances[2].strand, "+")
            self.assertEqual(motif.instances[3].strand, "+")
            self.assertEqual(motif.instances[4].strand, "+")
            self.assertEqual(motif.instances[5].strand, "+")
            self.assertEqual(motif.instances[6].strand, "+")
            self.assertEqual(motif.instances[7].strand, "+")
            self.assertEqual(motif.instances[8].strand, "+")
            self.assertEqual(motif.instances[9].strand, "+")
            self.assertEqual(motif.instances[10].strand, "+")
            self.assertEqual(motif.instances[11].strand, "+")
            self.assertEqual(motif.instances[12].strand, "+")
            self.assertEqual(motif.instances[13].strand, "+")
            self.assertEqual(motif.instances[14].strand, "+")
            self.assertEqual(motif.instances[15].strand, "+")
            self.assertEqual(motif.instances[16].strand, "+")
            self.assertEqual(motif.instances[17].strand, "+")
            self.assertEqual(motif.instances[18].strand, "+")
            self.assertEqual(motif.instances[19].strand, "+")
            self.assertEqual(motif.instances[20].strand, "+")
            self.assertEqual(motif.instances[0].length, 14)
            self.assertEqual(motif.instances[1].length, 14)
            self.assertEqual(motif.instances[2].length, 14)
            self.assertEqual(motif.instances[3].length, 14)
            self.assertEqual(motif.instances[4].length, 14)
            self.assertEqual(motif.instances[5].length, 14)
            self.assertEqual(motif.instances[6].length, 14)
            self.assertEqual(motif.instances[7].length, 14)
            self.assertEqual(motif.instances[8].length, 14)
            self.assertEqual(motif.instances[9].length, 14)
            self.assertEqual(motif.instances[10].length, 14)
            self.assertEqual(motif.instances[11].length, 14)
            self.assertEqual(motif.instances[12].length, 14)
            self.assertEqual(motif.instances[13].length, 14)
            self.assertEqual(motif.instances[14].length, 14)
            self.assertEqual(motif.instances[15].length, 14)
            self.assertEqual(motif.instances[16].length, 14)
            self.assertEqual(motif.instances[17].length, 14)
            self.assertEqual(motif.instances[18].length, 14)
            self.assertEqual(motif.instances[19].length, 14)
            self.assertEqual(motif.instances[20].length, 14)
            self.assertEqual(motif.instances[0].start, 66)
            self.assertEqual(motif.instances[1].start, 254)
            self.assertEqual(motif.instances[2].start, 258)
            self.assertEqual(motif.instances[3].start, 162)
            self.assertEqual(motif.instances[4].start, 190)
            self.assertEqual(motif.instances[5].start, 278)
            self.assertEqual(motif.instances[6].start, 172)
            self.assertEqual(motif.instances[7].start, 114)
            self.assertEqual(motif.instances[8].start, 7)
            self.assertEqual(motif.instances[9].start, 268)
            self.assertEqual(motif.instances[10].start, 414)
            self.assertEqual(motif.instances[11].start, 126)
            self.assertEqual(motif.instances[12].start, 220)
            self.assertEqual(motif.instances[13].start, 55)
            self.assertEqual(motif.instances[14].start, 229)
            self.assertEqual(motif.instances[15].start, 330)
            self.assertEqual(motif.instances[16].start, 18)
            self.assertEqual(motif.instances[17].start, 180)
            self.assertEqual(motif.instances[18].start, 73)
            self.assertEqual(motif.instances[19].start, 107)
            self.assertEqual(motif.instances[20].start, 36)
            self.assertEqual(motif.instances[0], "MNKEEILVFIKSCQ")
            self.assertEqual(motif.instances[1], "INYEKLTEFILKCQ")
            self.assertEqual(motif.instances[2], "IDREKLRSFILACQ")
            self.assertEqual(motif.instances[3], "INVEKAIEFVLSCM")
            self.assertEqual(motif.instances[4], "IDTEKLLGYIMSQQ")
            self.assertEqual(motif.instances[5], "INVEKLLEWSSARQ")
            self.assertEqual(motif.instances[6], "INREKLLQYLYSLK")
            self.assertEqual(motif.instances[7], "INVDKVVAYVQSLQ")
            self.assertEqual(motif.instances[8], "LLKEKHIRYIESLD")
            self.assertEqual(motif.instances[9], "LNLKSLLQWVTSRQ")
            self.assertEqual(motif.instances[10], "ENVRKIIHYFKSNL")
            self.assertEqual(motif.instances[11], "LDKRSLARFVSKCQ")
            self.assertEqual(motif.instances[12], "DLFEGTAEWIARCQ")
            self.assertEqual(motif.instances[13], "FVKEEVISFVLSCW")
            self.assertEqual(motif.instances[14], "ELTEGVLNYLKNCQ")
            self.assertEqual(motif.instances[15], "FNKHALRDYILYCC")
            self.assertEqual(motif.instances[16], "LLLEKHADYIASYG")
            self.assertEqual(motif.instances[17], "IDRKGIYQWLISLK")
            self.assertEqual(motif.instances[18], "LQREKHFHYLKRGL")
            self.assertEqual(motif.instances[19], "DRKVRLISFIRGNQ")
            self.assertEqual(motif.instances[20], "VNRMAIIFYSISGL")
        self.assertEqual(motif.consensus, "INKEKLIEYILSCQ")
        self.assertEqual(motif[3:-8].consensus, "EKL")

    def test_minimal_meme_parser(self):
        """Parse motifs/minimal_test.meme file."""
        with open("motifs/minimal_test.meme") as stream:
            record = motifs.parse(stream, "minimal")
        self.assertEqual(record.version, "4")
        self.assertEqual(record.alphabet, "ACGT")
        self.assertEqual(len(record.sequences), 0)
        self.assertEqual(record.command, "")
        self.assertEqual(len(record), 2)
        motif = record[0]
        self.assertEqual(motif.name, "KRP")
        self.assertEqual(record["KRP"], motif)
        self.assertEqual(motif.num_occurrences, 17)
        self.assertEqual(motif.length, 19)
        self.assertAlmostEqual(motif.background["A"], 0.30269730269730266)
        self.assertAlmostEqual(motif.background["C"], 0.1828171828171828)
        self.assertAlmostEqual(motif.background["G"], 0.20879120879120877)
        self.assertAlmostEqual(motif.background["T"], 0.30569430569430567)
        self.assertAlmostEqual(motif.evalue, 4.1e-09)
        self.assertEqual(motif.alphabet, "ACGT")
        self.assertIsNone(motif.alignment)
        # using the old instances property:
        with self.assertWarns(BiopythonDeprecationWarning):
            self.assertIsNone(motif.instances)
        self.assertEqual(motif.consensus, "TGTGATCGAGGTCACACTT")
        self.assertEqual(motif.degenerate_consensus, "TGTGANNNWGNTCACAYWW")
        self.assertTrue(
            np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        1.1684297174927525,
                        0.9432809925744818,
                        1.4307101633876265,
                        1.1549413780465179,
                        0.9308256303218774,
                        0.009164393966550805,
                        0.20124190687894253,
                        0.17618542656995528,
                        0.36777933103380855,
                        0.6635834532368525,
                        0.07729943368061855,
                        0.9838293592717438,
                        1.72489868427398,
                        0.8397561713453014,
                        1.72489868427398,
                        0.8455332015343343,
                        0.3106481207768122,
                        0.7382733641762232,
                        0.537435993300495,
                    ]
                ),
            )
        )
        self.assertEqual(motif[2:9].consensus, "TGATCGA")
        motif = record[1]
        self.assertEqual(motif.name, "IFXA")
        self.assertEqual(record["IFXA"], motif)
        self.assertEqual(motif.num_occurrences, 14)
        self.assertEqual(motif.length, 18)
        self.assertAlmostEqual(motif.background["A"], 0.30269730269730266)
        self.assertAlmostEqual(motif.background["C"], 0.1828171828171828)
        self.assertAlmostEqual(motif.background["G"], 0.20879120879120877)
        self.assertAlmostEqual(motif.background["T"], 0.30569430569430567)
        self.assertAlmostEqual(motif.evalue, 4.1e-09)
        self.assertEqual(motif.alphabet, "ACGT")
        self.assertIsNone(motif.alignment)
        self.assertEqual(motif.consensus, "TACTGTATATATATCCAG")
        self.assertEqual(motif.degenerate_consensus, "TACTGTATATAHAWMCAG")
        self.assertTrue(
            np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        0.9632889858595118,
                        1.02677956765017,
                        2.451526420551951,
                        1.7098384161433415,
                        2.2598671267551107,
                        1.7098384161433415,
                        1.02677956765017,
                        1.391583804103081,
                        1.02677956765017,
                        1.1201961888781142,
                        0.27822438781180836,
                        0.36915366971717867,
                        1.7240522753630425,
                        0.3802185945622609,
                        0.790937683007783,
                        2.451526420551951,
                        1.7240522753630425,
                        1.3924085743645374,
                    ]
                ),
            )
        )
        self.assertEqual(motif[2:9].consensus, "CTGTATA")
        # using the old instances property:
        with self.assertWarns(BiopythonDeprecationWarning):
            self.assertIsNone(motif.instances)

    def test_meme_parser_rna(self):
        """Test if Bio.motifs can parse MEME output files using RNA."""


class TestMAST(unittest.TestCase):
    """MAST format tests."""

    def test_mast_parser_1(self):
        """Parse motifs/mast.crp0.de.oops.txt.xml file."""
        with open("motifs/mast.crp0.de.oops.txt.xml") as stream:
            record = motifs.parse(stream, "MAST")
        self.assertEqual(record.version, "5.0.1")
        self.assertEqual(record.database, "common/crp0.s")
        self.assertEqual(record.alphabet, "DNA")
        self.assertEqual(len(record), 2)
        self.assertEqual(len(record.sequences), 18)
        self.assertEqual(record.sequences[0], "lac")
        self.assertEqual(record.sequences[1], "bglr1")
        self.assertEqual(record.sequences[2], "tdc")
        self.assertEqual(record.sequences[3], "deop2")
        self.assertEqual(record.sequences[4], "pbr322")
        self.assertEqual(record.sequences[5], "malk")
        self.assertEqual(record.sequences[6], "tnaa")
        self.assertEqual(record.sequences[7], "male")
        self.assertEqual(record.sequences[8], "ara")
        self.assertEqual(record.sequences[9], "cya")
        self.assertEqual(record.sequences[10], "ompa")
        self.assertEqual(record.sequences[11], "ilv")
        self.assertEqual(record.sequences[12], "gale")
        self.assertEqual(record.sequences[13], "malt")
        self.assertEqual(record.sequences[14], "crp")
        self.assertEqual(record.sequences[15], "ce1cg")
        self.assertEqual(record.sequences[16], "trn9cat")
        self.assertEqual(record.sequences[17], "uxu1")
        self.assertEqual(record.diagrams["lac"], "[+1]-2-[-2]-79")
        self.assertEqual(record.diagrams["bglr1"], "79-[+2]-14")
        self.assertEqual(record.diagrams["tdc"], "30-[+1]-39-[+2]-12")
        self.assertEqual(record.diagrams["deop2"], "19-[+1]-74")
        self.assertEqual(record.diagrams["pbr322"], "58-[-2]-35")
        self.assertEqual(record.diagrams["malk"], "32-[+2]-61")
        self.assertEqual(record.diagrams["tnaa"], "105")
        self.assertEqual(record.diagrams["male"], "105")
        self.assertEqual(record.diagrams["ara"], "105")
        self.assertEqual(record.diagrams["cya"], "105")
        self.assertEqual(record.diagrams["ompa"], "105")
        self.assertEqual(record.diagrams["ilv"], "105")
        self.assertEqual(record.diagrams["gale"], "105")
        self.assertEqual(record.diagrams["malt"], "105")
        self.assertEqual(record.diagrams["crp"], "105")
        self.assertEqual(record.diagrams["ce1cg"], "105")
        self.assertEqual(record.diagrams["trn9cat"], "105")
        self.assertEqual(record.diagrams["uxu1"], "105")
        motif = record[0]
        self.assertIs(record["1"], motif)
        self.assertEqual(motif.name, "1")
        self.assertEqual(motif.alphabet, "DNA")
        self.assertEqual(motif.length, 12)
        self.assertEqual(motif[1:-2].length, 9)
        motif = record[1]
        self.assertIs(record["2"], motif)
        self.assertEqual(motif.name, "2")
        self.assertEqual(motif.alphabet, "DNA")
        self.assertEqual(motif.length, 12)
        self.assertEqual(motif[1:40:3].length, 4)

    def test_mast_parser_2(self):
        """Parse motifs/mast.adh.de.oops.html.xml file."""
        with open("motifs/mast.adh.de.oops.html.xml") as stream:
            record = motifs.parse(stream, "MAST")
        self.assertEqual(record.version, "5.0.1")
        self.assertEqual(record.database, "common/adh.s")
        self.assertEqual(record.alphabet, "Protein")
        self.assertEqual(len(record.sequences), 33)
        self.assertEqual(record.sequences[0], "ENTA_ECOLI")
        self.assertEqual(record.sequences[1], "DHII_HUMAN")
        self.assertEqual(record.sequences[2], "YINL_LISMO")
        self.assertEqual(record.sequences[3], "FIXR_BRAJA")
        self.assertEqual(record.sequences[4], "HDHA_ECOLI")
        self.assertEqual(record.sequences[5], "BUDC_KLETE")
        self.assertEqual(record.sequences[6], "AP27_MOUSE")
        self.assertEqual(record.sequences[7], "FVT1_HUMAN")
        self.assertEqual(record.sequences[8], "YRTP_BACSU")
        self.assertEqual(record.sequences[9], "DHMA_FLAS1")
        self.assertEqual(record.sequences[10], "HDE_CANTR")
        self.assertEqual(record.sequences[11], "3BHD_COMTE")
        self.assertEqual(record.sequences[12], "BDH_HUMAN")
        self.assertEqual(record.sequences[13], "2BHD_STREX")
        self.assertEqual(record.sequences[14], "BA72_EUBSP")
        self.assertEqual(record.sequences[15], "RIDH_KLEAE")
        self.assertEqual(record.sequences[16], "DHGB_BACME")
        self.assertEqual(record.sequences[17], "PCR_PEA")
        self.assertEqual(record.sequences[18], "RFBB_NEIGO")
        self.assertEqual(record.sequences[19], "BPHB_PSEPS")
        self.assertEqual(record.sequences[20], "DHB2_HUMAN")
        self.assertEqual(record.sequences[21], "NODG_RHIME")
        self.assertEqual(record.sequences[22], "MAS1_AGRRA")
        self.assertEqual(record.sequences[23], "DHCA_HUMAN")
        self.assertEqual(record.sequences[24], "DHES_HUMAN")
        self.assertEqual(record.sequences[25], "DHB3_HUMAN")
        self.assertEqual(record.sequences[26], "HMTR_LEIMA")
        self.assertEqual(record.sequences[27], "ADH_DROME")
        self.assertEqual(record.sequences[28], "YURA_MYXXA")
        self.assertEqual(record.sequences[29], "LIGD_PSEPA")
        self.assertEqual(record.sequences[30], "FABI_ECOLI")
        self.assertEqual(record.sequences[31], "GUTD_ECOLI")
        self.assertEqual(record.sequences[32], "CSGA_MYXXA")
        self.assertEqual(record.diagrams["ENTA_ECOLI"], "[1]-[2]-224")
        self.assertEqual(record.diagrams["DHII_HUMAN"], "29-[1]-[2]-239")
        self.assertEqual(record.diagrams["YINL_LISMO"], "[1]-[2]-224")
        self.assertEqual(record.diagrams["FIXR_BRAJA"], "43-[2]-149-[1]-62")
        self.assertEqual(record.diagrams["HDHA_ECOLI"], "6-[1]-[2]-144-[1]-69")
        self.assertEqual(record.diagrams["BUDC_KLETE"], "9-[2]-53-[1]-81-[1]-62")
        self.assertEqual(record.diagrams["AP27_MOUSE"], "2-[1]-[2]-138-[1]-68")
        self.assertEqual(record.diagrams["FVT1_HUMAN"], "39-[2]-150-[1]-119")
        self.assertEqual(record.diagrams["YRTP_BACSU"], "1-[1]-[2]-145-[1]-56")
        self.assertEqual(record.diagrams["DHMA_FLAS1"], "9-[1]-[2]-147-[1]-78")
        self.assertEqual(record.diagrams["HDE_CANTR"], "3-[1]-[2]-290-[1]-[2]-565")
        self.assertEqual(record.diagrams["3BHD_COMTE"], "1-[1]-[2]-50-[1]-166")
        self.assertEqual(record.diagrams["BDH_HUMAN"], "50-[1]-[2]-269")
        self.assertEqual(record.diagrams["2BHD_STREX"], "1-[1]-[2]-142-[1]-76")
        self.assertEqual(record.diagrams["BA72_EUBSP"], "1-[1]-[2]-125-[2]-10-[1]-65")
        self.assertEqual(record.diagrams["RIDH_KLEAE"], "9-[1]-[2]-216")
        self.assertEqual(record.diagrams["DHGB_BACME"], "2-[1]-[2]-149-[1]-75")
        self.assertEqual(record.diagrams["PCR_PEA"], "81-[1]-[2]-108-[1]-174")
        self.assertEqual(record.diagrams["RFBB_NEIGO"], "1-[1]-[2]-321")
        self.assertEqual(record.diagrams["BPHB_PSEPS"], "[1]-[2]-251")
        self.assertEqual(record.diagrams["DHB2_HUMAN"], "77-[1]-[2]-286")
        self.assertEqual(record.diagrams["NODG_RHIME"], "1-[1]-[2]-142-[1]-66")
        self.assertEqual(record.diagrams["MAS1_AGRRA"], "252-[2]-36-[1]-164")
        self.assertEqual(record.diagrams["DHCA_HUMAN"], "11-[2]-54-[1]-101-[1]-74")
        self.assertEqual(record.diagrams["DHES_HUMAN"], "9-[2]-108-[1]-186")
        self.assertEqual(record.diagrams["DHB3_HUMAN"], "55-[2]-146-[1]-85")
        self.assertEqual(record.diagrams["HMTR_LEIMA"], "24-[2]-172-[1]-67")
        self.assertEqual(record.diagrams["ADH_DROME"], "13-[2]-217-[1]-1")
        self.assertEqual(record.diagrams["YURA_MYXXA"], "94-[2]-69-[1]-71")
        self.assertEqual(record.diagrams["LIGD_PSEPA"], "1-[1]-[2]-280")
        self.assertEqual(record.diagrams["FABI_ECOLI"], "1-[1]-161-[1]-76")
        self.assertEqual(record.diagrams["GUTD_ECOLI"], "147-[2]-10-[1]-78")
        self.assertEqual(record.diagrams["CSGA_MYXXA"], "12-[1]-53-[2]-77")
        self.assertEqual(len(record), 2)
        motif = record[0]
        self.assertIs(record["1"], motif)
        self.assertEqual(motif.alphabet, "Protein")
        self.assertEqual(motif.length, 12)
        self.assertEqual(motif.name, "1")
        self.assertEqual(motif[1:-2].length, 9)
        motif = record[1]
        self.assertIs(record["2"], motif)
        self.assertEqual(motif.alphabet, "Protein")
        self.assertEqual(motif.length, 12)
        self.assertEqual(motif.name, "2")
        self.assertEqual(motif[-20:-2].length, 10)

    def test_mast_parser_3(self):
        """Parse motifs/mast.Klf1-200.cd.oops.xml.xml file."""
        with open("motifs/mast.Klf1-200.cd.oops.xml.xml") as stream:
            record = motifs.parse(stream, "MAST")
        self.assertEqual(record.version, "5.0.1")
        self.assertEqual(record.database, "common/Klf1-200.fa")
        self.assertEqual(record.alphabet, "DNA")
        self.assertEqual(len(record.sequences), 113)
        self.assertEqual(record.sequences[0], "chr3:104843905-104844405")
        self.assertEqual(record.sequences[1], "chr12:114390660-114391160")
        self.assertEqual(record.sequences[2], "chr12:27135944-27136444")
        self.assertEqual(record.sequences[3], "chr10:59256089-59256589")
        self.assertEqual(record.sequences[4], "chr4:135733850-135734350")
        self.assertEqual(record.sequences[5], "chr1:137838164-137838664")
        self.assertEqual(record.sequences[6], "chr17:47735006-47735506")
        self.assertEqual(record.sequences[7], "chr6:72223026-72223526")
        self.assertEqual(record.sequences[8], "chr13:3866266-3866766")
        self.assertEqual(record.sequences[9], "chr1:133343883-133344383")
        self.assertEqual(record.sequences[10], "chr11:117187372-117187872")
        self.assertEqual(record.sequences[11], "chr13:76003199-76003699")
        self.assertEqual(record.sequences[12], "chr5:65202593-65203093")
        self.assertEqual(record.sequences[13], "chr14:79702844-79703344")
        self.assertEqual(record.sequences[14], "chr12:112796794-112797294")
        self.assertEqual(record.sequences[15], "chr13:112863645-112864145")
        self.assertEqual(record.sequences[16], "chr7:111007530-111008030")
        self.assertEqual(record.sequences[17], "chr1:43307690-43308190")
        self.assertEqual(record.sequences[18], "chr14:47973722-47974222")
        self.assertEqual(record.sequences[19], "chr9:120025371-120025871")
        self.assertEqual(record.sequences[20], "chr7:105490727-105491227")
        self.assertEqual(record.sequences[21], "chr5:37127175-37127675")
        self.assertEqual(record.sequences[22], "chr5:45951565-45952065")
        self.assertEqual(record.sequences[23], "chr7:91033422-91033922")
        self.assertEqual(record.sequences[24], "chr4:154285745-154286245")
        self.assertEqual(record.sequences[25], "chr13:100518008-100518508")
        self.assertEqual(record.sequences[26], "chr1:36977019-36977519")
        self.assertEqual(record.sequences[27], "chr7:151917814-151918314")
        self.assertEqual(record.sequences[28], "chr7:110976195-110976695")
        self.assertEqual(record.sequences[29], "chr15:58719281-58719781")
        self.assertEqual(record.sequences[30], "chr11:57590460-57590960")
        self.assertEqual(record.sequences[31], "chr8:83025150-83025650")
        self.assertEqual(record.sequences[32], "chr13:54345922-54346422")
        self.assertEqual(record.sequences[33], "chr12:82044358-82044858")
        self.assertEqual(record.sequences[34], "chr11:105013714-105014214")
        self.assertEqual(record.sequences[35], "chr10:93585404-93585904")
        self.assertEqual(record.sequences[36], "chr7:19832207-19832707")
        self.assertEqual(record.sequences[37], "chr8:97323995-97324495")
        self.assertEqual(record.sequences[38], "chr10:126642277-126642777")
        self.assertEqual(record.sequences[39], "chr1:156887119-156887619")
        self.assertEqual(record.sequences[40], "chr15:81700367-81700867")
        self.assertEqual(record.sequences[41], "chr6:121187425-121187925")
        self.assertEqual(record.sequences[42], "chr4:43977111-43977611")
        self.assertEqual(record.sequences[43], "chr11:102236405-102236905")
        self.assertEqual(record.sequences[44], "chr17:5112057-5112557")
        self.assertEqual(record.sequences[45], "chr10:110604369-110604869")
        self.assertEqual(record.sequences[46], "chr1:169314208-169314708")
        self.assertEqual(record.sequences[47], "chr9:57618594-57619094")
        self.assertEqual(record.sequences[48], "chr10:128184604-128185104")
        self.assertEqual(record.sequences[49], "chr4:109112541-109113041")
        self.assertEqual(record.sequences[50], "chr3:97461668-97462168")
        self.assertEqual(record.sequences[51], "chr9:102674395-102674895")
        self.assertEqual(record.sequences[52], "chr17:24289205-24289705")
        self.assertEqual(record.sequences[53], "chr17:28960252-28960752")
        self.assertEqual(record.sequences[54], "chr2:73323093-73323593")
        self.assertEqual(record.sequences[55], "chr11:32150818-32151318")
        self.assertEqual(record.sequences[56], "chr7:103853792-103854292")
        self.assertEqual(record.sequences[57], "chr16:49839621-49840121")
        self.assertEqual(record.sequences[58], "chr6:135115628-135116128")
        self.assertEqual(record.sequences[59], "chr3:88305500-88306000")
        self.assertEqual(record.sequences[60], "chr18:57137388-57137888")
        self.assertEqual(record.sequences[61], "chr5:97380648-97381148")
        self.assertEqual(record.sequences[62], "chr15:91082416-91082916")
        self.assertEqual(record.sequences[63], "chr14:61272713-61273213")
        self.assertEqual(record.sequences[64], "chr5:33616214-33616714")
        self.assertEqual(record.sequences[65], "chr18:23982470-23982970")
        self.assertEqual(record.sequences[66], "chr9:24715045-24715545")
        self.assertEqual(record.sequences[67], "chr10:116195445-116195945")
        self.assertEqual(record.sequences[68], "chr11:77795184-77795684")
        self.assertEqual(record.sequences[69], "chr16:32508975-32509475")
        self.assertEqual(record.sequences[70], "chr18:80416880-80417380")
        self.assertEqual(record.sequences[71], "chr10:57252236-57252736")
        self.assertEqual(record.sequences[72], "chr5:34915767-34916267")
        self.assertEqual(record.sequences[73], "chr9:98389943-98390443")
        self.assertEqual(record.sequences[74], "chr19:5845899-5846399")
        self.assertEqual(record.sequences[75], "chr3:151777796-151778296")
        self.assertEqual(record.sequences[76], "chr4:76585120-76585620")
        self.assertEqual(record.sequences[77], "chr7:104332488-104332988")
        self.assertEqual(record.sequences[78], "chr5:138127197-138127697")
        self.assertEqual(record.sequences[79], "chr11:60988820-60989320")
        self.assertEqual(record.sequences[80], "chr8:19984030-19984530")
        self.assertEqual(record.sequences[81], "chr11:31712262-31712762")
        self.assertEqual(record.sequences[82], "chr15:41338514-41339014")
        self.assertEqual(record.sequences[83], "chr9:21362671-21363171")
        self.assertEqual(record.sequences[84], "chr18:58822702-58823202")
        self.assertEqual(record.sequences[85], "chr1:173447614-173448114")
        self.assertEqual(record.sequences[86], "chr6:81915769-81916269")
        self.assertEqual(record.sequences[87], "chr1:169322898-169323398")
        self.assertEqual(record.sequences[88], "chr12:70860461-70860961")
        self.assertEqual(record.sequences[89], "chr9:59598186-59598686")
        self.assertEqual(record.sequences[90], "chr3:19550495-19550995")
        self.assertEqual(record.sequences[91], "chr7:36132953-36133453")
        self.assertEqual(record.sequences[92], "chr7:38970375-38970875")
        self.assertEqual(record.sequences[93], "chr15:78243390-78243890")
        self.assertEqual(record.sequences[94], "chr7:87847381-87847881")
        self.assertEqual(record.sequences[95], "chr1:33631214-33631714")
        self.assertEqual(record.sequences[96], "chr4:135407873-135408373")
        self.assertEqual(record.sequences[97], "chr7:101244829-101245329")
        self.assertEqual(record.sequences[98], "chr10:60612190-60612690")
        self.assertEqual(record.sequences[99], "chr19:56465963-56466463")
        self.assertEqual(record.sequences[100], "chr4:41334759-41335259")
        self.assertEqual(record.sequences[101], "chr8:92969521-92970021")
        self.assertEqual(record.sequences[102], "chr6:145703215-145703715")
        self.assertEqual(record.sequences[103], "chr13:57679178-57679678")
        self.assertEqual(record.sequences[104], "chr19:45121628-45122128")
        self.assertEqual(record.sequences[105], "chr15:79757891-79758391")
        self.assertEqual(record.sequences[106], "chr1:134264178-134264678")
        self.assertEqual(record.sequences[107], "chr13:81067500-81068000")
        self.assertEqual(record.sequences[108], "chr11:69714224-69714724")
        self.assertEqual(record.sequences[109], "chr2:103728071-103728571")
        self.assertEqual(record.sequences[110], "chr5:105994747-105995247")
        self.assertEqual(record.sequences[111], "chr17:84209565-84210065")
        self.assertEqual(record.sequences[112], "chr7:16507689-16508189")
        self.assertEqual(
            record.diagrams["chr3:104843905-104844405"], "115-[-1]-209-[-2]-126"
        )
        self.assertEqual(
            record.diagrams["chr12:114390660-114391160"],
            "3-[+2]-[+2]-3-[+1]-173-[+1]-3-[-2]-188",
        )
        self.assertEqual(
            record.diagrams["chr12:27135944-27136444"], "275-[-1]-89-[+2]-4-[+2]-52"
        )
        self.assertEqual(
            record.diagrams["chr10:59256089-59256589"], "247-[+2]-17-[-1]-186"
        )
        self.assertEqual(
            record.diagrams["chr4:135733850-135734350"], "183-[-1]-263-[+2]-4"
        )
        self.assertEqual(
            record.diagrams["chr1:137838164-137838664"], "192-[-2]-1-[+1]-44-[-1]-193"
        )
        self.assertEqual(
            record.diagrams["chr17:47735006-47735506"], "203-[+2]-15-[+1]-97-[-1]-115"
        )
        self.assertEqual(
            record.diagrams["chr6:72223026-72223526"],
            "52-[-2]-7-[+2]-162-[-1]-42-[-1]-137",
        )
        self.assertEqual(
            record.diagrams["chr13:3866266-3866766"], "241-[+1]-2-[-1]-217"
        )
        self.assertEqual(
            record.diagrams["chr1:133343883-133344383"], "190-[+2]-15-[+1]-245"
        )
        self.assertEqual(
            record.diagrams["chr11:117187372-117187872"], "242-[+1]-46-[-2]-71-[+1]-71"
        )
        self.assertEqual(
            record.diagrams["chr13:76003199-76003699"], "230-[+2]-15-[+2]-60-[-1]-115"
        )
        self.assertEqual(
            record.diagrams["chr5:65202593-65203093"],
            "24-[-2]-36-[+2]-193-[-1]-11-[+1]-10-[+1]-106",
        )
        self.assertEqual(
            record.diagrams["chr14:79702844-79703344"], "247-[-1]-46-[-2]-157"
        )
        self.assertEqual(
            record.diagrams["chr12:112796794-112797294"], "232-[+1]-41-[+1]-187"
        )
        self.assertEqual(
            record.diagrams["chr13:112863645-112864145"], "228-[+1]-20-[-1]-212"
        )
        self.assertEqual(
            record.diagrams["chr7:111007530-111008030"], "217-[+1]-83-[+2]-150"
        )
        self.assertEqual(
            record.diagrams["chr1:43307690-43308190"], "164-[-2]-52-[-2]-224"
        )
        self.assertEqual(
            record.diagrams["chr14:47973722-47974222"], "21-[+1]-181-[+1]-20-[-2]-208"
        )
        self.assertEqual(
            record.diagrams["chr9:120025371-120025871"], "110-[-2]-58-[+1]-282"
        )
        self.assertEqual(
            record.diagrams["chr7:105490727-105491227"], "100-[-2]-111-[-1]-239"
        )
        self.assertEqual(
            record.diagrams["chr5:37127175-37127675"], "234-[-2]-24-[+1]-192"
        )
        self.assertEqual(record.diagrams["chr5:45951565-45952065"], "261-[-1]-219")
        self.assertEqual(record.diagrams["chr7:91033422-91033922"], "465-[-1]-15")
        self.assertEqual(
            record.diagrams["chr4:154285745-154286245"], "235-[+1]-20-[-2]-195"
        )
        self.assertEqual(
            record.diagrams["chr13:100518008-100518508"], "226-[-2]-18-[-1]-206"
        )
        self.assertEqual(
            record.diagrams["chr1:36977019-36977519"], "88-[+1]-187-[+2]-60-[-1]-95"
        )
        self.assertEqual(
            record.diagrams["chr7:151917814-151918314"], "219-[+1]-80-[+2]-151"
        )
        self.assertEqual(
            record.diagrams["chr7:110976195-110976695"], "287-[+2]-12-[+1]-151"
        )
        self.assertEqual(record.diagrams["chr15:58719281-58719781"], "212-[-2]-258")
        self.assertEqual(
            record.diagrams["chr11:57590460-57590960"], "56-[-1]-271-[-1]-75-[+2]-28"
        )
        self.assertEqual(
            record.diagrams["chr8:83025150-83025650"], "219-[+1]-87-[+2]-144"
        )
        self.assertEqual(
            record.diagrams["chr13:54345922-54346422"], "283-[-2]-161-[+1]-6"
        )
        self.assertEqual(
            record.diagrams["chr12:82044358-82044858"], "50-[+2]-160-[+1]-39-[+2]-171"
        )
        self.assertEqual(
            record.diagrams["chr11:105013714-105014214"],
            "115-[-2]-160-[+1]-26-[-1]-129",
        )
        self.assertEqual(
            record.diagrams["chr10:93585404-93585904"], "141-[+2]-48-[+1]-261"
        )
        self.assertEqual(record.diagrams["chr7:19832207-19832707"], "229-[-1]-251")
        self.assertEqual(
            record.diagrams["chr8:97323995-97324495"], "177-[-1]-40-[-2]-139-[+1]-74"
        )
        self.assertEqual(
            record.diagrams["chr10:126642277-126642777"], "252-[-1]-92-[-2]-106"
        )
        self.assertEqual(
            record.diagrams["chr1:156887119-156887619"], "189-[-2]-78-[-1]-183"
        )
        self.assertEqual(
            record.diagrams["chr15:81700367-81700867"], "109-[-1]-99-[-1]-252"
        )
        self.assertEqual(
            record.diagrams["chr6:121187425-121187925"], "29-[+2]-313-[-1]-108"
        )
        self.assertEqual(
            record.diagrams["chr4:43977111-43977611"], "60-[+1]-148-[+1]-252"
        )
        self.assertEqual(
            record.diagrams["chr11:102236405-102236905"],
            "10-[+2]-145-[-1]-3-[-1]-6-[+2]-60-[+1]-156",
        )
        self.assertEqual(record.diagrams["chr17:5112057-5112557"], "249-[+1]-231")
        self.assertEqual(record.diagrams["chr10:110604369-110604869"], "232-[+1]-248")
        self.assertEqual(
            record.diagrams["chr1:169314208-169314708"], "192-[-1]-[-1]-11-[-2]-227"
        )
        self.assertEqual(
            record.diagrams["chr9:57618594-57619094"], "125-[+2]-151-[-1]-4-[-1]-150"
        )
        self.assertEqual(
            record.diagrams["chr10:128184604-128185104"], "30-[-2]-128-[+1]-292"
        )
        self.assertEqual(
            record.diagrams["chr4:109112541-109113041"], "21-[-1]-13-[+1]-94-[+2]-302"
        )
        self.assertEqual(
            record.diagrams["chr3:97461668-97462168"],
            "18-[+2]-256-[-1]-81-[+1]-21-[+1]-34",
        )
        self.assertEqual(record.diagrams["chr9:102674395-102674895"], "372-[+2]-98")
        self.assertEqual(record.diagrams["chr17:24289205-24289705"], "262-[-1]-218")
        self.assertEqual(
            record.diagrams["chr17:28960252-28960752"], "221-[+1]-81-[+1]-158"
        )
        self.assertEqual(record.diagrams["chr2:73323093-73323593"], "49-[-2]-421")
        self.assertEqual(
            record.diagrams["chr11:32150818-32151318"], "151-[-1]-27-[-1]-118-[-2]-134"
        )
        self.assertEqual(
            record.diagrams["chr7:103853792-103854292"], "212-[-2]-42-[+1]-196"
        )
        self.assertEqual(
            record.diagrams["chr16:49839621-49840121"], "192-[+2]-47-[-1]-17-[+2]-164"
        )
        self.assertEqual(record.diagrams["chr6:135115628-135116128"], "231-[-1]-249")
        self.assertEqual(record.diagrams["chr3:88305500-88306000"], "229-[+1]-251")
        self.assertEqual(record.diagrams["chr18:57137388-57137888"], "296-[+2]-174")
        self.assertEqual(record.diagrams["chr5:97380648-97381148"], "188-[-2]-282")
        self.assertEqual(
            record.diagrams["chr15:91082416-91082916"], "239-[-1]-104-[-1]-73-[+2]-14"
        )
        self.assertEqual(
            record.diagrams["chr14:61272713-61273213"], "216-[+2]-104-[+1]-130"
        )
        self.assertEqual(record.diagrams["chr5:33616214-33616714"], "247-[-1]-233")
        self.assertEqual(record.diagrams["chr18:23982470-23982970"], "285-[-1]-195")
        self.assertEqual(
            record.diagrams["chr9:24715045-24715545"], "214-[-1]-153-[+1]-93"
        )
        self.assertEqual(record.diagrams["chr10:116195445-116195945"], "400-[+2]-70")
        self.assertEqual(
            record.diagrams["chr11:77795184-77795684"], "247-[+1]-42-[-2]-67-[-2]-64"
        )
        self.assertEqual(
            record.diagrams["chr16:32508975-32509475"], "213-[+2]-29-[-1]-208"
        )
        self.assertEqual(record.diagrams["chr18:80416880-80417380"], "239-[-1]-241")
        self.assertEqual(
            record.diagrams["chr10:57252236-57252736"], "155-[+1]-158-[+2]-137"
        )
        self.assertEqual(
            record.diagrams["chr5:34915767-34916267"], "179-[+2]-29-[-1]-242"
        )
        self.assertEqual(record.diagrams["chr9:98389943-98390443"], "252-[-1]-228")
        self.assertEqual(
            record.diagrams["chr19:5845899-5846399"], "136-[+1]-193-[+1]-131"
        )
        self.assertEqual(
            record.diagrams["chr3:151777796-151778296"], "30-[-2]-58-[-1]-362"
        )
        self.assertEqual(record.diagrams["chr4:76585120-76585620"], "329-[+2]-141")
        self.assertEqual(
            record.diagrams["chr7:104332488-104332988"], "164-[+2]-23-[-1]-222-[+1]-21"
        )
        self.assertEqual(record.diagrams["chr5:138127197-138127697"], "238-[+1]-242")
        self.assertEqual(
            record.diagrams["chr11:60988820-60989320"], "115-[+1]-68-[+1]-47-[+1]-210"
        )
        self.assertEqual(
            record.diagrams["chr8:19984030-19984530"], "103-[-1]-81-[+2]-266"
        )
        self.assertEqual(
            record.diagrams["chr11:31712262-31712762"], "118-[+2]-53-[+2]-269"
        )
        self.assertEqual(
            record.diagrams["chr15:41338514-41339014"], "173-[+2]-75-[+2]-192"
        )
        self.assertEqual(
            record.diagrams["chr9:21362671-21363171"], "105-[+1]-131-[+1]-224"
        )
        self.assertEqual(record.diagrams["chr18:58822702-58823202"], "467-[-2]-3")
        self.assertEqual(record.diagrams["chr1:173447614-173448114"], "369-[-1]-111")
        self.assertEqual(record.diagrams["chr6:81915769-81916269"], "197-[+1]-283")
        self.assertEqual(record.diagrams["chr1:169322898-169323398"], "253-[-1]-227")
        self.assertEqual(
            record.diagrams["chr12:70860461-70860961"], "197-[+2]-22-[-1]-231"
        )
        self.assertEqual(
            record.diagrams["chr9:59598186-59598686"], "163-[-2]-10-[-1]-277"
        )
        self.assertEqual(record.diagrams["chr3:19550495-19550995"], "452-[-2]-18")
        self.assertEqual(record.diagrams["chr7:36132953-36133453"], "157-[-1]-323")
        self.assertEqual(
            record.diagrams["chr7:38970375-38970875"], "49-[+1]-114-[+1]-297"
        )
        self.assertEqual(record.diagrams["chr15:78243390-78243890"], "234-[+1]-246")
        self.assertEqual(
            record.diagrams["chr7:87847381-87847881"], "99-[+2]-2-[-1]-230-[-1]-99"
        )
        self.assertEqual(record.diagrams["chr1:33631214-33631714"], "358-[-1]-122")
        self.assertEqual(
            record.diagrams["chr4:135407873-135408373"], "116-[-1]-64-[+2]-270"
        )
        self.assertEqual(record.diagrams["chr7:101244829-101245329"], "311-[-2]-159")
        self.assertEqual(record.diagrams["chr10:60612190-60612690"], "215-[+1]-265")
        self.assertEqual(
            record.diagrams["chr19:56465963-56466463"], "306-[+1]-36-[+1]-18-[+1]-80"
        )
        self.assertEqual(record.diagrams["chr4:41334759-41335259"], "204-[+1]-276")
        self.assertEqual(record.diagrams["chr8:92969521-92970021"], "453-[+2]-17")
        self.assertEqual(
            record.diagrams["chr6:145703215-145703715"], "154-[-2]-58-[+2]-228"
        )
        self.assertEqual(record.diagrams["chr13:57679178-57679678"], "217-[-1]-263")
        self.assertEqual(record.diagrams["chr19:45121628-45122128"], "35-[-2]-435")
        self.assertEqual(record.diagrams["chr15:79757891-79758391"], "310-[+1]-170")
        self.assertEqual(record.diagrams["chr1:134264178-134264678"], "23-[+2]-447")
        self.assertEqual(record.diagrams["chr13:81067500-81068000"], "252-[+1]-228")
        self.assertEqual(record.diagrams["chr11:69714224-69714724"], "145-[+2]-325")
        self.assertEqual(record.diagrams["chr2:103728071-103728571"], "369-[+1]-111")
        self.assertEqual(
            record.diagrams["chr5:105994747-105995247"], "93-[+2]-153-[-2]-194"
        )
        self.assertEqual(record.diagrams["chr17:84209565-84210065"], "64-[-2]-406")
        self.assertEqual(record.diagrams["chr7:16507689-16508189"], "231-[+2]-239")
        self.assertEqual(len(record), 2)
        motif = record[0]
        self.assertIs(record["1"], motif)
        self.assertEqual(motif.alphabet, "DNA")
        self.assertEqual(motif.length, 20)
        self.assertEqual(motif.name, "1")
        self.assertEqual(motif[1:-2].length, 17)
        motif = record[1]
        self.assertIs(record["2"], motif)
        self.assertEqual(motif.alphabet, "DNA")
        self.assertEqual(motif.length, 30)
        self.assertEqual(motif.name, "2")
        self.assertEqual(motif[10:20].length, 10)


class TestTransfac(unittest.TestCase):
    """Transfac format tests."""

    def test_transfac_parser(self):
        """Parse motifs/transfac.dat file."""
        with open("motifs/transfac.dat") as stream:
            record = motifs.parse(stream, "TRANSFAC")
        motif = record[0]
        self.assertEqual(motif["ID"], "motif1")
        self.assertEqual(len(motif.counts), 4)
        self.assertEqual(motif.counts.length, 12)
        self.assertEqual(motif.counts["A", 0], 1)
        self.assertEqual(motif.counts["A", 1], 2)
        self.assertEqual(motif.counts["A", 2], 3)
        self.assertEqual(motif.counts["A", 3], 0)
        self.assertEqual(motif.counts["A", 4], 5)
        self.assertEqual(motif.counts["A", 5], 0)
        self.assertEqual(motif.counts["A", 6], 0)
        self.assertEqual(motif.counts["A", 7], 0)
        self.assertEqual(motif.counts["A", 8], 0)
        self.assertEqual(motif.counts["A", 9], 0)
        self.assertEqual(motif.counts["A", 10], 0)
        self.assertEqual(motif.counts["A", 11], 1)
        self.assertEqual(motif.counts["C", 0], 2)
        self.assertEqual(motif.counts["C", 1], 1)
        self.assertEqual(motif.counts["C", 2], 0)
        self.assertEqual(motif.counts["C", 3], 5)
        self.assertEqual(motif.counts["C", 4], 0)
        self.assertEqual(motif.counts["C", 5], 0)
        self.assertEqual(motif.counts["C", 6], 1)
        self.assertEqual(motif.counts["C", 7], 0)
        self.assertEqual(motif.counts["C", 8], 0)
        self.assertEqual(motif.counts["C", 9], 1)
        self.assertEqual(motif.counts["C", 10], 2)
        self.assertEqual(motif.counts["C", 11], 0)
        self.assertEqual(motif.counts["G", 0], 2)
        self.assertEqual(motif.counts["G", 1], 2)
        self.assertEqual(motif.counts["G", 2], 1)
        self.assertEqual(motif.counts["G", 3], 0)
        self.assertEqual(motif.counts["G", 4], 0)
        self.assertEqual(motif.counts["G", 5], 4)
        self.assertEqual(motif.counts["G", 6], 4)
        self.assertEqual(motif.counts["G", 7], 0)
        self.assertEqual(motif.counts["G", 8], 5)
        self.assertEqual(motif.counts["G", 9], 2)
        self.assertEqual(motif.counts["G", 10], 0)
        self.assertEqual(motif.counts["G", 11], 3)
        self.assertEqual(motif.counts["T", 0], 0)
        self.assertEqual(motif.counts["T", 1], 0)
        self.assertEqual(motif.counts["T", 2], 1)
        self.assertEqual(motif.counts["T", 3], 0)
        self.assertEqual(motif.counts["T", 4], 0)
        self.assertEqual(motif.counts["T", 5], 1)
        self.assertEqual(motif.counts["T", 6], 0)
        self.assertEqual(motif.counts["T", 7], 5)
        self.assertEqual(motif.counts["T", 8], 0)
        self.assertEqual(motif.counts["T", 9], 2)
        self.assertEqual(motif.counts["T", 10], 3)
        self.assertEqual(motif.counts["T", 11], 1)
        self.assertEqual(motif.degenerate_consensus, "SRACAGGTGKYG")
        self.assertTrue(
            np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        0.4780719051126377,
                        0.4780719051126377,
                        0.6290494055453314,
                        2.0,
                        2.0,
                        1.278071905112638,
                        1.278071905112638,
                        2.0,
                        2.0,
                        0.4780719051126377,
                        1.0290494055453312,
                        0.6290494055453314,
                    ]
                ),
            )
        )
        self.assertEqual(motif[1:-2].degenerate_consensus, "RACAGGTGK")
        self.assertTrue(
            np.allclose(
                motif[1:-2].relative_entropy,
                np.array(
                    [
                        0.4780719051126377,
                        0.6290494055453314,
                        2.0,
                        2.0,
                        1.278071905112638,
                        1.278071905112638,
                        2.0,
                        2.0,
                        0.4780719051126377,
                    ]
                ),
            )
        )
        motif = record[1]
        self.assertEqual(motif["ID"], "motif2")
        self.assertEqual(len(motif.counts), 4)
        self.assertEqual(motif.counts.length, 10)
        self.assertEqual(motif.counts["A", 0], 2)
        self.assertEqual(motif.counts["A", 1], 1)
        self.assertEqual(motif.counts["A", 2], 0)
        self.assertEqual(motif.counts["A", 3], 3)
        self.assertEqual(motif.counts["A", 4], 0)
        self.assertEqual(motif.counts["A", 5], 5)
        self.assertEqual(motif.counts["A", 6], 0)
        self.assertEqual(motif.counts["A", 7], 0)
        self.assertEqual(motif.counts["A", 8], 0)
        self.assertEqual(motif.counts["A", 9], 0)
        self.assertEqual(motif.counts["C", 0], 1)
        self.assertEqual(motif.counts["C", 1], 2)
        self.assertEqual(motif.counts["C", 2], 5)
        self.assertEqual(motif.counts["C", 3], 0)
        self.assertEqual(motif.counts["C", 4], 0)
        self.assertEqual(motif.counts["C", 5], 0)
        self.assertEqual(motif.counts["C", 6], 1)
        self.assertEqual(motif.counts["C", 7], 0)
        self.assertEqual(motif.counts["C", 8], 0)
        self.assertEqual(motif.counts["C", 9], 2)
        self.assertEqual(motif.counts["G", 0], 2)
        self.assertEqual(motif.counts["G", 1], 2)
        self.assertEqual(motif.counts["G", 2], 0)
        self.assertEqual(motif.counts["G", 3], 1)
        self.assertEqual(motif.counts["G", 4], 4)
        self.assertEqual(motif.counts["G", 5], 0)
        self.assertEqual(motif.counts["G", 6], 4)
        self.assertEqual(motif.counts["G", 7], 5)
        self.assertEqual(motif.counts["G", 8], 0)
        self.assertEqual(motif.counts["G", 9], 0)
        self.assertEqual(motif.counts["T", 0], 0)
        self.assertEqual(motif.counts["T", 1], 0)
        self.assertEqual(motif.counts["T", 2], 0)
        self.assertEqual(motif.counts["T", 3], 1)
        self.assertEqual(motif.counts["T", 4], 1)
        self.assertEqual(motif.counts["T", 5], 0)
        self.assertEqual(motif.counts["T", 6], 0)
        self.assertEqual(motif.counts["T", 7], 0)
        self.assertEqual(motif.counts["T", 8], 5)
        self.assertEqual(motif.counts["T", 9], 3)
        self.assertEqual(motif.degenerate_consensus, "RSCAGAGGTY")
        self.assertTrue(
            np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        0.4780719051126377,
                        0.4780719051126377,
                        2.0,
                        0.6290494055453314,
                        1.278071905112638,
                        2.0,
                        1.278071905112638,
                        2.0,
                        2.0,
                        1.0290494055453312,
                    ]
                ),
            )
        )
        self.assertEqual(motif[::2].degenerate_consensus, "RCGGT")
        self.assertTrue(
            np.allclose(
                motif[::2].relative_entropy,
                np.array(
                    [0.4780719051126377, 2.0, 1.278071905112638, 1.278071905112638, 2.0]
                ),
            )
        )

    def test_permissive_transfac_parser(self):
        """Parse the TRANSFAC-like file motifs/MA0056.1.transfac."""
        # The test file MA0056.1.transfac was obtained from the JASPAR database
        # in a TRANSFAC-like format.
        # Khan, A. et al. JASPAR 2018: update of the open-access database of
        # transcription factor binding profiles and its web framework.
        # Nucleic Acids Res. 2018; 46:D260-D266,
        path = "motifs/MA0056.1.transfac"
        with open(path) as stream:
            self.assertRaises(ValueError, motifs.parse, stream, "TRANSFAC")
        with open(path) as stream:
            records = motifs.parse(stream, "TRANSFAC", strict=False)
        motif = records[0]
        self.assertEqual(sorted(motif.keys()), ["AC", "DE", "ID"])
        self.assertEqual(motif["AC"], "MA0056.1")
        self.assertEqual(motif["DE"], "MA0056.1 MZF1 ; From JASPAR 2018")
        self.assertEqual(motif["ID"], "MZF1")
        self.assertEqual(motif.counts.length, 6)
        self.assertEqual(len(motif.counts), 4)
        self.assertEqual(motif.counts["A", 0], 3.0)
        self.assertEqual(motif.counts["A", 1], 0.0)
        self.assertEqual(motif.counts["A", 2], 2.0)
        self.assertEqual(motif.counts["A", 3], 0.0)
        self.assertEqual(motif.counts["A", 4], 0.0)
        self.assertEqual(motif.counts["A", 5], 18.0)
        self.assertEqual(motif.counts["C", 0], 5.0)
        self.assertEqual(motif.counts["C", 1], 0.0)
        self.assertEqual(motif.counts["C", 2], 0.0)
        self.assertEqual(motif.counts["C", 3], 0.0)
        self.assertEqual(motif.counts["C", 4], 0.0)
        self.assertEqual(motif.counts["C", 5], 0.0)
        self.assertEqual(motif.counts["G", 0], 4.0)
        self.assertEqual(motif.counts["G", 1], 19.0)
        self.assertEqual(motif.counts["G", 2], 18.0)
        self.assertEqual(motif.counts["G", 3], 19.0)
        self.assertEqual(motif.counts["G", 4], 20.0)
        self.assertEqual(motif.counts["G", 5], 2.0)
        self.assertEqual(motif.counts["T", 0], 8.0)
        self.assertEqual(motif.counts["T", 1], 1.0)
        self.assertEqual(motif.counts["T", 2], 0.0)
        self.assertEqual(motif.counts["T", 3], 1.0)
        self.assertEqual(motif.counts["T", 4], 0.0)
        self.assertEqual(motif.counts["T", 5], 0.0)
        self.assertEqual(motif.consensus, "TGGGGA")
        self.assertEqual(motif.degenerate_consensus, "NGGGGA")
        self.assertTrue(
            np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        0.09629830394265171,
                        1.7136030428840439,
                        1.5310044064107189,
                        1.7136030428840439,
                        2.0,
                        1.5310044064107189,
                    ]
                ),
            )
        )
        self.assertEqual(motif[1:-3].degenerate_consensus, "GG")
        self.assertTrue(
            np.allclose(
                motif[1:-3].relative_entropy,
                np.array([1.7136030428840439, 1.5310044064107189]),
            )
        )

    def test_TFoutput(self):
        """Ensure that we can write proper TransFac output files."""
        m = motifs.create([Seq("ATATA")])
        with tempfile.TemporaryFile("w") as stream:
            stream.write(format(m, "transfac"))


class MotifTestPWM(unittest.TestCase):
    """PWM motif tests."""

    with open("motifs/SRF.pfm") as stream:
        m = motifs.read(stream, "pfm")

    s = Seq("ACGTGTGCGTAGTGCGT")

    def test_getitem(self):
        counts = self.m.counts
        python_integers = range(13)
        numpy_integers = np.array(python_integers)
        integers = {"python": python_integers, "numpy": numpy_integers}
        for int_type in ("python", "numpy"):
            i0, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12 = integers[int_type]
            msg = f"using {int_type} integers as indices"
            # slice, slice
            d = counts[i1::i2, i2:i12:i3]
            self.assertIsInstance(d, dict, msg=msg)
            self.assertEqual(len(d), 2, msg=msg)
            self.assertEqual(len(d["C"]), 4, msg=msg)
            self.assertEqual(len(d["T"]), 4, msg=msg)
            self.assertAlmostEqual(d["C"][i0], 45.0, msg=msg)
            self.assertAlmostEqual(d["C"][i1], 1.0, msg=msg)
            self.assertAlmostEqual(d["C"][i2], 0.0, msg=msg)
            self.assertAlmostEqual(d["C"][i3], 1.0, msg=msg)
            self.assertAlmostEqual(d["T"][i0], 0.0, msg=msg)
            self.assertAlmostEqual(d["T"][i1], 42.0, msg=msg)
            self.assertAlmostEqual(d["T"][i2], 3.0, msg=msg)
            self.assertAlmostEqual(d["T"][i3], 0.0, msg=msg)
            # slice, int
            d = counts[i1::i2, i4]
            self.assertIsInstance(d, dict, msg=msg)
            self.assertEqual(len(d), 2, msg=msg)
            self.assertAlmostEqual(d["C"], 1.0, msg=msg)
            self.assertAlmostEqual(d["T"], 13.0, msg=msg)
            # int, slice
            t = counts[i2, i3:i12:i2]
            self.assertIsInstance(t, tuple, msg=msg)
            self.assertAlmostEqual(t[i0], 0.0, msg=msg)
            self.assertAlmostEqual(t[i1], 0.0, msg=msg)
            self.assertAlmostEqual(t[i2], 0.0, msg=msg)
            self.assertAlmostEqual(t[i3], 0.0, msg=msg)
            self.assertAlmostEqual(t[i4], 43.0, msg=msg)
            # int, int
            v = counts[i1, i5]
            self.assertAlmostEqual(v, 1.0, msg=msg)
            # tuple, slice
            d = counts[(i0, i3), i3:i12:i2]
            self.assertIsInstance(d, dict, msg=msg)
            self.assertEqual(len(d), 2, msg=msg)
            self.assertEqual(len(d["A"]), 5, msg=msg)
            self.assertEqual(len(d["T"]), 5, msg=msg)
            self.assertAlmostEqual(d["A"][i0], 1.0, msg=msg)
            self.assertAlmostEqual(d["A"][i1], 3.0, msg=msg)
            self.assertAlmostEqual(d["A"][i2], 1.0, msg=msg)
            self.assertAlmostEqual(d["A"][i3], 15.0, msg=msg)
            self.assertAlmostEqual(d["A"][i4], 2.0, msg=msg)
            self.assertAlmostEqual(d["T"][i0], 0.0, msg=msg)
            self.assertAlmostEqual(d["T"][i1], 42.0, msg=msg)
            self.assertAlmostEqual(d["T"][i2], 45.0, msg=msg)
            self.assertAlmostEqual(d["T"][i3], 30.0, msg=msg)
            self.assertAlmostEqual(d["T"][i4], 0.0, msg=msg)
            # tuple, int
            d = counts[(i0, i3), i5]
            self.assertIsInstance(d, dict, msg=msg)
            self.assertEqual(len(d), 2, msg=msg)
            self.assertAlmostEqual(d["A"], 3.0, msg=msg)
            self.assertAlmostEqual(d["T"], 42.0, msg=msg)
            # str, slice
            t = counts["C", i2:i12:i4]
            self.assertIsInstance(t, tuple, msg=msg)
            self.assertAlmostEqual(t[i0], 45.0, msg=msg)
            self.assertAlmostEqual(t[i1], 0.0, msg=msg)
            self.assertAlmostEqual(t[i2], 0.0, msg=msg)
            # str, int
            self.assertAlmostEqual(counts["T", i4], 13.0, msg=msg)

    def test_simple(self):
        """Test if Bio.motifs PWM scoring works."""
        counts = self.m.counts
        pwm = counts.normalize(pseudocounts=0.25)
        pssm = pwm.log_odds()
        result = pssm.calculate(self.s)
        self.assertEqual(6, len(result))
        # The fast C-code in Bio/motifs/_pwm.c stores all results as 32-bit
        # floats; the slower Python code in Bio/motifs/__init__.py uses 64-bit
        # doubles. The C-code and Python code results will therefore not be
        # exactly equal. Test the first 5 decimal places only to avoid either
        # the C-code or the Python code to inadvertently fail this test.
        self.assertAlmostEqual(result[0], -29.18363571, places=5)
        self.assertAlmostEqual(result[1], -38.3365097, places=5)
        self.assertAlmostEqual(result[2], -29.17756271, places=5)
        self.assertAlmostEqual(result[3], -38.04542542, places=5)
        self.assertAlmostEqual(result[4], -20.3014183, places=5)
        self.assertAlmostEqual(result[5], -25.18009186, places=5)

    def test_with_mixed_case(self):
        """Test if Bio.motifs PWM scoring works with mixed case."""
        counts = self.m.counts
        pwm = counts.normalize(pseudocounts=0.25)
        pssm = pwm.log_odds()
        result = pssm.calculate(Seq("AcGTgTGCGtaGTGCGT"))
        self.assertEqual(6, len(result))
        self.assertAlmostEqual(result[0], -29.18363571, places=5)
        self.assertAlmostEqual(result[1], -38.3365097, places=5)
        self.assertAlmostEqual(result[2], -29.17756271, places=5)
        self.assertAlmostEqual(result[3], -38.04542542, places=5)
        self.assertAlmostEqual(result[4], -20.3014183, places=5)
        self.assertAlmostEqual(result[5], -25.18009186, places=5)

    def test_with_bad_char(self):
        """Test if Bio.motifs PWM scoring works with unexpected letters like N."""
        counts = self.m.counts
        pwm = counts.normalize(pseudocounts=0.25)
        pssm = pwm.log_odds()
        result = pssm.calculate(Seq("ACGTGTGCGTAGTGCGTN"))
        self.assertEqual(7, len(result))
        self.assertAlmostEqual(result[0], -29.18363571, places=5)
        self.assertAlmostEqual(result[1], -38.3365097, places=5)
        self.assertAlmostEqual(result[2], -29.17756271, places=5)
        self.assertAlmostEqual(result[3], -38.04542542, places=5)
        self.assertAlmostEqual(result[4], -20.3014183, places=5)
        self.assertAlmostEqual(result[5], -25.18009186, places=5)
        self.assertTrue(math.isnan(result[6]), f"Expected nan, not {result[6]!r}")

    def test_calculate_pseudocounts(self):
        pseudocounts = motifs.jaspar.calculate_pseudocounts(self.m)
        self.assertAlmostEqual(pseudocounts["A"], 1.695582495781317, places=5)
        self.assertAlmostEqual(pseudocounts["C"], 1.695582495781317, places=5)
        self.assertAlmostEqual(pseudocounts["G"], 1.695582495781317, places=5)
        self.assertAlmostEqual(pseudocounts["T"], 1.695582495781317, places=5)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
