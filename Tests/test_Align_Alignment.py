# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


"""Tests for the Alignment class in Bio.Align."""

import os
import unittest

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install numpy if you want to use Bio.Align."
    ) from None

from Bio import Align, SeqIO
from Bio.Seq import Seq, reverse_complement
from Bio.SeqUtils import GC


class TestPairwiseAlignment(unittest.TestCase):
    def check_indexing_slicing(self, alignment, msg):
        self.assertEqual(
            repr(alignment),
            "<Bio.Align.Alignment object (2 rows x 12 columns) at 0x%x>"
            % id(alignment),
        )
        self.assertEqual(
            str(alignment),
            """\
AACCGGGA-CCG
|-|-||-|-|--
A-C-GG-AAC--
""",
            msg=msg,
        )
        self.assertAlmostEqual(alignment.score, 6.0)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 12))
        self.assertEqual(alignment[0], "AACCGGGA-CCG", msg=msg)
        self.assertEqual(alignment[1], "A-C-GG-AAC--", msg=msg)
        self.assertEqual(alignment[-2], "AACCGGGA-CCG", msg=msg)
        self.assertEqual(alignment[-1], "A-C-GG-AAC--", msg=msg)
        self.assertEqual(alignment[0, :], "AACCGGGA-CCG", msg=msg)
        self.assertEqual(alignment[1, :], "A-C-GG-AAC--", msg=msg)
        self.assertEqual(alignment[-2, :], "AACCGGGA-CCG", msg=msg)
        self.assertEqual(alignment[-1, :], "A-C-GG-AAC--", msg=msg)
        self.assertEqual(alignment[:, 0], "AA", msg=msg)
        self.assertEqual(alignment[:, 1], "A-", msg=msg)
        self.assertEqual(alignment[:, 2], "CC", msg=msg)
        self.assertEqual(alignment[:, 3], "C-", msg=msg)
        self.assertEqual(alignment[:, 4], "GG", msg=msg)
        self.assertEqual(alignment[:, 5], "GG", msg=msg)
        self.assertEqual(alignment[:, 6], "G-", msg=msg)
        self.assertEqual(alignment[:, 7], "AA", msg=msg)
        self.assertEqual(alignment[:, 8], "-A", msg=msg)
        self.assertEqual(alignment[:, 9], "CC", msg=msg)
        self.assertEqual(alignment[:, 10], "C-", msg=msg)
        self.assertEqual(alignment[:, 11], "G-", msg=msg)
        self.assertEqual(alignment[:, -12], "AA", msg=msg)
        self.assertEqual(alignment[:, -11], "A-", msg=msg)
        self.assertEqual(alignment[:, -10], "CC", msg=msg)
        self.assertEqual(alignment[:, -9], "C-", msg=msg)
        self.assertEqual(alignment[:, -8], "GG", msg=msg)
        self.assertEqual(alignment[:, -7], "GG", msg=msg)
        self.assertEqual(alignment[:, -6], "G-", msg=msg)
        self.assertEqual(alignment[:, -5], "AA", msg=msg)
        self.assertEqual(alignment[:, -4], "-A", msg=msg)
        self.assertEqual(alignment[:, -3], "CC", msg=msg)
        self.assertEqual(alignment[:, -2], "C-", msg=msg)
        self.assertEqual(alignment[:, -1], "G-", msg=msg)
        self.assertEqual(alignment[1, range(1, 12, 2)], "--GAC-", msg=msg)
        self.assertEqual(alignment[0, (1, 4, 9)], "AGC", msg=msg)
        self.assertEqual(alignment[1, (1, 4, 9)], "-GC", msg=msg)
        self.assertEqual(alignment[0, range(0, 12, 2)], "ACGG-C", msg=msg)
        self.assertAlmostEqual(alignment[:, :].score, 6.0, msg=msg)
        self.assertEqual(
            str(alignment[:, :]),
            """\
AACCGGGA-CCG
|-|-||-|-|--
A-C-GG-AAC--
""",
            msg=msg,
        )
        self.assertEqual(alignment[0, 0:12], "AACCGGGA-CCG", msg=msg)
        self.assertEqual(alignment[1, 0:12], "A-C-GG-AAC--", msg=msg)
        self.assertEqual(alignment[0, 0:], "AACCGGGA-CCG", msg=msg)
        self.assertEqual(alignment[1, 0:], "A-C-GG-AAC--", msg=msg)
        self.assertAlmostEqual(alignment[:, 0:].score, 6.0, msg=msg)
        self.assertEqual(
            str(alignment[:, 0:]),
            """\
AACCGGGA-CCG
|-|-||-|-|--
A-C-GG-AAC--
""",
            msg=msg,
        )
        self.assertEqual(alignment[0, :12], "AACCGGGA-CCG", msg=msg)
        self.assertEqual(alignment[1, :12], "A-C-GG-AAC--", msg=msg)
        self.assertEqual(alignment[0, 1:], "ACCGGGA-CCG", msg=msg)
        self.assertEqual(alignment[1, 1:], "-C-GG-AAC--", msg=msg)
        self.assertEqual(alignment[0, 2:], "CCGGGA-CCG", msg=msg)
        self.assertEqual(alignment[1, 2:], "C-GG-AAC--", msg=msg)
        self.assertEqual(alignment[0, 3:], "CGGGA-CCG", msg=msg)
        self.assertEqual(alignment[1, 3:], "-GG-AAC--", msg=msg)
        self.assertEqual(alignment[0, 4:], "GGGA-CCG", msg=msg)
        self.assertEqual(alignment[1, 4:], "GG-AAC--", msg=msg)
        self.assertEqual(alignment[0, 5:], "GGA-CCG", msg=msg)
        self.assertEqual(alignment[1, 5:], "G-AAC--", msg=msg)
        self.assertEqual(alignment[0, 6:], "GA-CCG", msg=msg)
        self.assertEqual(alignment[1, 6:], "-AAC--", msg=msg)
        self.assertEqual(alignment[0, 7:], "A-CCG", msg=msg)
        self.assertEqual(alignment[1, 7:], "AAC--", msg=msg)
        self.assertEqual(alignment[0, 8:], "-CCG", msg=msg)
        self.assertEqual(alignment[1, 8:], "AC--", msg=msg)
        self.assertEqual(alignment[0, 9:], "CCG", msg=msg)
        self.assertEqual(alignment[1, 9:], "C--", msg=msg)
        self.assertEqual(alignment[0, 10:], "CG", msg=msg)
        self.assertEqual(alignment[1, 10:], "--", msg=msg)
        self.assertEqual(alignment[0, 11:], "G", msg=msg)
        self.assertEqual(alignment[1, 11:], "-", msg=msg)
        self.assertEqual(alignment[0, 12:], "", msg=msg)
        self.assertEqual(alignment[1, 12:], "", msg=msg)
        self.assertEqual(alignment[0, :-1], "AACCGGGA-CC", msg=msg)
        self.assertEqual(alignment[1, :-1], "A-C-GG-AAC-", msg=msg)
        self.assertEqual(alignment[0, :-2], "AACCGGGA-C", msg=msg)
        self.assertEqual(alignment[1, :-2], "A-C-GG-AAC", msg=msg)
        self.assertEqual(alignment[0, :-3], "AACCGGGA-", msg=msg)
        self.assertEqual(alignment[1, :-3], "A-C-GG-AA", msg=msg)
        self.assertEqual(alignment[0, 1:-1], "ACCGGGA-CC", msg=msg)
        self.assertEqual(alignment[1, 1:-1], "-C-GG-AAC-", msg=msg)
        self.assertEqual(alignment[0, 1:-2], "ACCGGGA-C", msg=msg)
        self.assertEqual(alignment[1, 1:-2], "-C-GG-AAC", msg=msg)
        self.assertEqual(alignment[0, 2:-1], "CCGGGA-CC", msg=msg)
        self.assertEqual(alignment[1, 2:-1], "C-GG-AAC-", msg=msg)
        self.assertEqual(alignment[0, 2:-2], "CCGGGA-C", msg=msg)
        self.assertEqual(alignment[1, 2:-2], "C-GG-AAC", msg=msg)
        self.assertAlmostEqual(alignment[:, :12].score, 6.0, msg=msg)
        self.assertEqual(
            str(alignment[:, :12]),
            """\
AACCGGGA-CCG
|-|-||-|-|--
A-C-GG-AAC--
""",
            msg=msg,
        )
        self.assertAlmostEqual(alignment[:, 0:12].score, 6.0, msg=msg)
        self.assertEqual(
            str(alignment[:, 0:12]),
            """\
AACCGGGA-CCG
|-|-||-|-|--
A-C-GG-AAC--
""",
            msg=msg,
        )
        self.assertEqual(
            str(alignment[:, 1:]),
            """\
AACCGGGA-CCG
 -|-||-|-|--
A-C-GG-AAC--
""",
            msg=msg,
        )
        self.assertEqual(
            str(alignment[:, 2:]),
            """\
AACCGGGA-CCG
  |-||-|-|--
 AC-GG-AAC--
""",
            msg=msg,
        )
        self.assertEqual(
            str(alignment[:, 3:]),
            """\
AACCGGGA-CCG
   -||-|-|--
 AC-GG-AAC--
""",
            msg=msg,
        )
        self.assertEqual(
            str(alignment[:, 4:]),
            """\
AACCGGGA-CCG
    ||-|-|--
  ACGG-AAC--
""",
            msg=msg,
        )
        self.assertEqual(
            str(alignment[:, :-1]),
            """\
AACCGGGA-CCG
|-|-||-|-|-
A-C-GG-AAC-
""",
            msg=msg,
        )
        self.assertEqual(
            str(alignment[:, :-2]),
            """\
AACCGGGA-CCG
|-|-||-|-|
A-C-GG-AAC
""",
            msg=msg,
        )
        self.assertEqual(
            str(alignment[:, :-3]),
            """\
AACCGGGA-CCG
|-|-||-|-
A-C-GG-AAC
""",
            msg=msg,
        )
        self.assertEqual(
            str(alignment[:, 1:-1]),
            """\
AACCGGGA-CCG
 -|-||-|-|-
A-C-GG-AAC-
""",
            msg=msg,
        )
        self.assertEqual(
            str(alignment[:, 1:-2]),
            """\
AACCGGGA-CCG
 -|-||-|-|
A-C-GG-AAC
""",
            msg=msg,
        )
        self.assertEqual(
            str(alignment[:, 2:-1]),
            """\
AACCGGGA-CCG
  |-||-|-|-
 AC-GG-AAC-
""",
            msg=msg,
        )
        self.assertEqual(
            str(alignment[:, 2:-2]),
            """\
AACCGGGA-CCG
  |-||-|-|
 AC-GG-AAC
""",
            msg=msg,
        )
        self.assertEqual(
            str(alignment[:, ::2]),
            """\
ACGG-C
|||---
ACG-A-
""",
            msg=msg,
        )
        self.assertEqual(
            str(alignment[:, range(0, 12, 2)]),
            """\
ACGG-C
|||---
ACG-A-
""",
            msg=msg,
        )
        self.assertEqual(
            str(alignment[:, (1, 8, 5)]),
            """\
A-G
--|
-AG
""",
            msg=msg,
        )
        with self.assertRaises(NotImplementedError, msg=msg):
            alignment[:1]
        with self.assertRaises(NotImplementedError, msg=msg):
            alignment[:1, :]

    def test_indexing_slicing(self):
        target = "AACCGGGACCG"
        query = "ACGGAAC"
        sequences = (target, query)
        coordinates = numpy.array(
            [[0, 1, 2, 3, 4, 6, 7, 8, 8, 9, 11], [0, 1, 1, 2, 2, 4, 4, 5, 6, 7, 7]]
        )
        alignment = Align.Alignment(sequences, coordinates)
        alignment.score = 6.0
        msg = "forward strand"
        self.check_indexing_slicing(alignment, msg)
        query = reverse_complement(query)
        sequences = (target, query)
        coordinates = numpy.array(
            [[0, 1, 2, 3, 4, 6, 7, 8, 8, 9, 11], [7, 6, 6, 5, 5, 3, 3, 2, 1, 0, 0]]
        )
        alignment = Align.Alignment(sequences, coordinates)
        alignment.score = 6.0
        msg = "reverse strand"
        self.check_indexing_slicing(alignment, msg)

    def test_sort(self):
        target = Seq("ACTT")
        query = Seq("ACCT")
        sequences = (target, query)
        coordinates = numpy.array([[0, 4], [0, 4]])
        alignment = Align.Alignment(sequences, coordinates)
        self.assertEqual(
            str(alignment),
            """\
ACTT
||.|
ACCT
""",
        )
        alignment.sort()
        self.assertEqual(
            str(alignment),
            """\
ACCT
||.|
ACTT
""",
        )
        alignment.sort(reverse=True)
        self.assertEqual(
            str(alignment),
            """\
ACTT
||.|
ACCT
""",
        )
        target.id = "seq1"
        query.id = "seq2"
        alignment.sort()
        self.assertEqual(
            str(alignment),
            """\
ACTT
||.|
ACCT
""",
        )
        alignment.sort(reverse=True)
        self.assertEqual(
            str(alignment),
            """\
ACCT
||.|
ACTT
""",
        )
        alignment.sort(key=GC)
        self.assertEqual(
            str(alignment),
            """\
ACTT
||.|
ACCT
""",
        )
        alignment.sort(key=GC, reverse=True)
        self.assertEqual(
            str(alignment),
            """\
ACCT
||.|
ACTT
""",
        )

    def test_substitutions(self):
        path = os.path.join("Align", "ecoli.fa")
        record = SeqIO.read(path, "fasta")
        target = record.seq
        path = os.path.join("Align", "bsubtilis.fa")
        record = SeqIO.read(path, "fasta")
        query = record.seq
        coordinates = numpy.array(
            [
                [
                    503,
                    744,
                    744,
                    747,
                    748,
                    820,
                    820,
                    822,
                    822,
                    823,
                    823,
                    828,
                    828,
                    833,
                    833,
                    845,
                    848,
                    850,
                    851,
                    854,
                    857,
                    1003,
                    1004,
                    1011,
                    1011,
                    1017,
                    1017,
                    1020,
                    1021,
                    1116,
                    1116,
                    1119,
                    1120,
                    1132,
                    1133,
                    1242,
                    1243,
                    1246,
                    1246,
                    1289,
                    1289,
                    1292,
                    1293,
                    1413,
                ],
                [
                    512,
                    753,
                    754,
                    757,
                    757,
                    829,
                    831,
                    833,
                    834,
                    835,
                    838,
                    843,
                    844,
                    849,
                    850,
                    862,
                    862,
                    864,
                    864,
                    867,
                    867,
                    1013,
                    1013,
                    1020,
                    1021,
                    1027,
                    1028,
                    1031,
                    1031,
                    1126,
                    1127,
                    1130,
                    1130,
                    1142,
                    1142,
                    1251,
                    1251,
                    1254,
                    1255,
                    1298,
                    1299,
                    1302,
                    1302,
                    1422,
                ],
            ]
        )
        sequences = (target, query)
        alignment = Align.Alignment(sequences, coordinates)
        m = alignment.substitutions
        self.assertEqual(
            str(m),
            """\
      A     C     G     T
A 191.0   3.0  15.0  13.0
C   5.0 186.0   9.0  14.0
G  12.0  11.0 248.0   8.0
T  11.0  19.0   6.0 145.0
""",
        )
        self.assertAlmostEqual(m["T", "C"], 19.0)
        self.assertAlmostEqual(m["C", "T"], 14.0)
        m += m.transpose()
        m /= 2.0
        self.assertEqual(
            str(m),
            """\
      A     C     G     T
A 191.0   4.0  13.5  12.0
C   4.0 186.0  10.0  16.5
G  13.5  10.0 248.0   7.0
T  12.0  16.5   7.0 145.0
""",
        )
        self.assertAlmostEqual(m["C", "T"], 16.5)
        self.assertAlmostEqual(m["T", "C"], 16.5)

    def test_target_query_properties(self):
        target = "ABCD"
        query = "XYZ"
        sequences = [target, query]
        coordinates = numpy.array([[0, 3, 4], [0, 3, 3]])
        alignment = Align.Alignment(sequences, coordinates)
        self.assertEqual(alignment.sequences[0], target)
        self.assertEqual(alignment.sequences[1], query)
        self.assertEqual(alignment.target, target)
        self.assertEqual(alignment.query, query)
        target = "EFGH"
        query = "UVW"
        sequences = [target, query]
        alignment.sequences = sequences
        self.assertEqual(alignment.sequences[0], target)
        self.assertEqual(alignment.sequences[1], query)
        self.assertEqual(alignment.target, target)
        self.assertEqual(alignment.query, query)
        target = "IJKL"
        query = "RST"
        sequences = [target, query]
        alignment.sequences = sequences
        self.assertEqual(alignment.sequences[0], target)
        self.assertEqual(alignment.sequences[1], query)
        self.assertEqual(alignment.target, target)
        self.assertEqual(alignment.query, query)


class TestMultipleAlignment(unittest.TestCase):
    def setUp(self):
        from Bio.Align import clustal

        path = "Clustalw/opuntia.aln"
        with open(path) as stream:
            alignments = clustal.AlignmentIterator(stream)
            self.alignment = next(alignments)

    def tearDown(self):
        del self.alignment

    def test_target_query_properties(self):
        target = "ABCD"
        query = "XYZ"
        alignment = self.alignment
        with self.assertRaises(ValueError):
            alignment.target
        with self.assertRaises(ValueError):
            alignment.query
        with self.assertRaises(ValueError):
            alignment.target = target
        with self.assertRaises(ValueError):
            alignment.query = query

    def test_comparison(self):
        alignment = self.alignment
        self.assertEqual(alignment.shape, (7, 156))
        sequences = alignment.sequences
        coordinates = numpy.array(alignment.coordinates)
        other = Align.Alignment(sequences, coordinates)
        self.assertEqual(alignment, other)
        self.assertLessEqual(alignment, other)
        self.assertGreaterEqual(other, alignment)
        other = Align.Alignment(sequences, coordinates[:, 1:])
        self.assertNotEqual(alignment, other)
        self.assertLess(alignment, other)
        self.assertLessEqual(alignment, other)
        self.assertGreater(other, alignment)
        self.assertGreaterEqual(other, alignment)

    def check_indexing_slicing(self, alignment, msg):
        self.assertEqual(
            repr(alignment),
            "<Bio.Align.Alignment object (7 rows x 156 columns) at 0x%x>"
            % id(alignment),
        )
        # self.assertEqual(str(alignment), ..., msg=msg)  # FIXME
        self.assertEqual(len(alignment), 7)
        self.assertEqual(alignment.shape, (7, 156))
        self.assertEqual(
            alignment[0],
            "TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCCATTGATTTAGTGTACCAGA",
            msg=msg,
        )
        self.assertEqual(
            alignment[1],
            "TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA",
            msg=msg,
        )
        self.assertEqual(
            alignment[-2],
            "TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTATACCAGA",
            msg=msg,
        )
        self.assertEqual(
            alignment[-1],
            "TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATATATATAATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA",
            msg=msg,
        )
        self.assertEqual(
            alignment[0, :],
            "TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCCATTGATTTAGTGTACCAGA",
            msg=msg,
        )
        self.assertEqual(
            alignment[1, :],
            "TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA",
            msg=msg,
        )
        self.assertEqual(
            alignment[-2, :],
            "TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTATACCAGA",
            msg=msg,
        )
        self.assertEqual(
            alignment[-1, :],
            "TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATATATATAATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA",
            msg=msg,
        )
        self.assertEqual(alignment[:, 0], "TTTTTTT", msg=msg)
        self.assertEqual(alignment[:, 1], "AAAAAAA", msg=msg)
        self.assertEqual(alignment[:, 2], "TTTTTTT", msg=msg)
        self.assertEqual(alignment[:, 3], "AAAAAAA", msg=msg)
        self.assertEqual(alignment[:, 4], "CCCCCCC", msg=msg)
        self.assertEqual(alignment[:, 5], "AAAAAAA", msg=msg)
        self.assertEqual(alignment[:, 6], "TTTTTTT", msg=msg)
        self.assertEqual(alignment[:, 7], "TTTATTT", msg=msg)
        self.assertEqual(alignment[:, 8], "AAAAAAA", msg=msg)
        self.assertEqual(alignment[:, 9], "AAAAAAA", msg=msg)
        self.assertEqual(alignment[:, 10], "AAAAAAA", msg=msg)
        self.assertEqual(alignment[:, 11], "GGGGGGG", msg=msg)
        self.assertEqual(alignment[:, 12], "AAAAGGG", msg=msg)
        self.assertEqual(alignment[:, -156], "TTTTTTT", msg=msg)
        self.assertEqual(alignment[:, -155], "AAAAAAA", msg=msg)
        self.assertEqual(alignment[:, -154], "TTTTTTT", msg=msg)
        self.assertEqual(alignment[:, -9], "TTTTTTT", msg=msg)
        self.assertEqual(alignment[:, -8], "GGGGGAG", msg=msg)
        self.assertEqual(alignment[:, -7], "TTTTTTT", msg=msg)
        self.assertEqual(alignment[:, -6], "AAAAAAA", msg=msg)
        self.assertEqual(alignment[:, -5], "CCCCCCC", msg=msg)
        self.assertEqual(alignment[:, -4], "CCCCCCC", msg=msg)
        self.assertEqual(alignment[:, -3], "AAAAAAA", msg=msg)
        self.assertEqual(alignment[:, -2], "GGGGGGG", msg=msg)
        self.assertEqual(alignment[:, -1], "AAAAAAA", msg=msg)
        self.assertEqual(
            alignment[0, range(0, 156, 2)],
            "TTCTAAAGGGTCGTATGAGCAAAAATTT-----AAATCATTCTTTCCATTAATTTAAATGTATTAAATCTGTTGGACG",
            msg=msg,
        )
        self.assertEqual(
            alignment[1, range(1, 156, 2)],
            "AAATAGAGGAGGAAAGAAGGAGAGAAAAA----TTTTAATCTAAACAAAAAAACATATAAGAACAGACATATATTCAA",
            msg=msg,
        )
        self.assertEqual(alignment[0, (1, 4, 9)], "ACA", msg=msg)
        self.assertEqual(alignment[1, (1, 57, 58)], "AA-", msg=msg)
        # self.assertEqual(str(alignment[:, :]), ..., msg=msg)
        self.assertEqual(
            alignment[0, 0:156],
            "TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCCATTGATTTAGTGTACCAGA",
            msg=msg,
        )
        self.assertEqual(
            alignment[1, 0:156],
            "TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA",
            msg=msg,
        )
        self.assertEqual(
            alignment[0, 0:],
            "TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCCATTGATTTAGTGTACCAGA",
            msg=msg,
        )
        self.assertEqual(
            alignment[1, 0:],
            "TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA",
            msg=msg,
        )
        # self.assertEqual(str(alignment[:, 0:]), ..., msg=msg)
        self.assertEqual(
            alignment[0, :156],
            "TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCCATTGATTTAGTGTACCAGA",
            msg=msg,
        )
        self.assertEqual(
            alignment[1, :156],
            "TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA",
            msg=msg,
        )
        self.assertEqual(
            alignment[0, 1:],
            "ATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCCATTGATTTAGTGTACCAGA",
            msg=msg,
        )
        self.assertEqual(
            alignment[1, 1:],
            "ATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA",
            msg=msg,
        )
        self.assertEqual(
            alignment[0, 2:],
            "TACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCCATTGATTTAGTGTACCAGA",
            msg=msg,
        )
        self.assertEqual(
            alignment[1, 2:],
            "TACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA",
            msg=msg,
        )
        self.assertEqual(
            alignment[0, 60:],
            "------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCCATTGATTTAGTGTACCAGA",
            msg=msg,
        )
        self.assertEqual(
            alignment[1, 60:],
            "------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA",
            msg=msg,
        )
        self.assertEqual(alignment[0, 156:], "", msg=msg)
        self.assertEqual(alignment[1, 156:], "", msg=msg)
        self.assertEqual(
            alignment[0, :-1],
            "TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCCATTGATTTAGTGTACCAG",
            msg=msg,
        )
        self.assertEqual(
            alignment[1, :-1],
            "TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAG",
            msg=msg,
        )
        self.assertEqual(
            alignment[0, :-2],
            "TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCCATTGATTTAGTGTACCA",
            msg=msg,
        )
        self.assertEqual(
            alignment[1, :-2],
            "TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCA",
            msg=msg,
        )
        self.assertEqual(
            alignment[0, :-3],
            "TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCCATTGATTTAGTGTACC",
            msg=msg,
        )
        self.assertEqual(
            alignment[1, :-3],
            "TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACC",
            msg=msg,
        )
        self.assertEqual(
            alignment[0, 1:-1],
            "ATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCCATTGATTTAGTGTACCAG",
            msg=msg,
        )
        self.assertEqual(
            alignment[1, 1:-1],
            "ATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAG",
            msg=msg,
        )
        self.assertEqual(
            alignment[0, 1:-2],
            "ATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCCATTGATTTAGTGTACCA",
            msg=msg,
        )
        self.assertEqual(
            alignment[1, 1:-2],
            "ATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCA",
            msg=msg,
        )
        self.assertEqual(
            alignment[0, 2:-1],
            "TACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCCATTGATTTAGTGTACCAG",
            msg=msg,
        )
        self.assertEqual(
            alignment[1, 2:-1],
            "TACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAG",
            msg=msg,
        )
        self.assertEqual(
            alignment[0, 2:-2],
            "TACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCCATTGATTTAGTGTACCA",
            msg=msg,
        )
        self.assertEqual(
            alignment[1, 2:-2],
            "TACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCA",
            msg=msg,
        )
        # self.assertEqual(str(alignment[:, :156]), ..., msg=msg)
        # self.assertEqual(str(alignment[:, 0:156]), ..., msg=msg)
        # self.assertEqual(str(alignment[:, 1:]), ..., msg=msg)
        # self.assertEqual(str(alignment[:, 2:]), ..., msg=msg)
        # self.assertEqual(str(alignment[:, 3:]), ..., msg=msg)
        # self.assertEqual(str(alignment[:, 4:]), ..., msg=msg)
        # self.assertEqual(str(alignment[:, :-1]), ..., msg=msg)
        # self.assertEqual(str(alignment[:, :-2]), ..., msg=msg)
        # self.assertEqual(str(alignment[:, :-3]), ..., msg=msg)
        # self.assertEqual(str(alignment[:, 1:-1]), ..., msg=msg)
        # self.assertEqual(str(alignment[:, 1:-2]), ..., msg=msg)
        # self.assertEqual(str(alignment[:, 2:-1]), ..., msg=msg)
        # self.assertEqual(str(alignment[:, 2:-2]), ..., msg=msg)
        # self.assertEqual(str(alignment[:, ::2]), ..., msg=msg)
        # self.assertEqual(str(alignment[:, range(0, 156, 2)]), ..., msg=msg)
        # self.assertEqual(str(alignment[:, (1, 8, 5)]), ..., msg=msg)
        with self.assertRaises(NotImplementedError, msg=msg):
            alignment[:1]
        with self.assertRaises(NotImplementedError, msg=msg):
            alignment[:1, :]

    def test_indexing_slicing(self):
        alignment = self.alignment
        msg = "forward strand"
        self.check_indexing_slicing(alignment, msg)
        alignment.sequences[2] = alignment.sequences[2].reverse_complement()
        n = len(alignment.sequences[2])
        alignment.coordinates[2, :] = n - alignment.coordinates[2, :]
        msg = "reverse strand"
        self.check_indexing_slicing(alignment, msg)

    def test_sort(self):
        alignment = self.alignment[:, 40:100]
        self.assertEqual(
            "\n".join(row for row in alignment),  # str(alignment),
            """\
TATACATTAAAGAAGG----------GGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
TATACATTAAAGAAGGGG--------GATGCGGATAAATGGAAAGGCGAAAGAAAGAATA
TATACATTAAAGAAGG----------GGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
TATACATAAAAGAAGG----------GGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
TATACATTAAAGGAGGGGGA------TGCGGATAAATGGAAAGGCGAAAGAAAGAATATA
TATACATTAAAGGAGGGGGA------TGCGGATAAATGGAAAGGCGAAAGAAAGAATATA
TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA""",
        )
        self.assertEqual(
            tuple(sequence.id for sequence in alignment.sequences),
            (
                "gi|6273285|gb|AF191659.1|AF191",
                "gi|6273284|gb|AF191658.1|AF191",
                "gi|6273287|gb|AF191661.1|AF191",
                "gi|6273286|gb|AF191660.1|AF191",
                "gi|6273290|gb|AF191664.1|AF191",
                "gi|6273289|gb|AF191663.1|AF191",
                "gi|6273291|gb|AF191665.1|AF191",
            ),
        )
        alignment.sort()
        self.assertEqual(
            "\n".join(row for row in alignment),  # str(alignment),
            """\
TATACATTAAAGAAGGGG--------GATGCGGATAAATGGAAAGGCGAAAGAAAGAATA
TATACATTAAAGAAGG----------GGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
TATACATAAAAGAAGG----------GGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
TATACATTAAAGAAGG----------GGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
TATACATTAAAGGAGGGGGA------TGCGGATAAATGGAAAGGCGAAAGAAAGAATATA
TATACATTAAAGGAGGGGGA------TGCGGATAAATGGAAAGGCGAAAGAAAGAATATA
TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA""",
        )
        self.assertEqual(
            tuple(sequence.id for sequence in alignment.sequences),
            (
                "gi|6273284|gb|AF191658.1|AF191",
                "gi|6273285|gb|AF191659.1|AF191",
                "gi|6273286|gb|AF191660.1|AF191",
                "gi|6273287|gb|AF191661.1|AF191",
                "gi|6273289|gb|AF191663.1|AF191",
                "gi|6273290|gb|AF191664.1|AF191",
                "gi|6273291|gb|AF191665.1|AF191",
            ),
        )
        alignment.sort(reverse=True)
        self.assertEqual(
            "\n".join(row for row in alignment),  # str(alignment),
            """\
TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA
TATACATTAAAGGAGGGGGA------TGCGGATAAATGGAAAGGCGAAAGAAAGAATATA
TATACATTAAAGGAGGGGGA------TGCGGATAAATGGAAAGGCGAAAGAAAGAATATA
TATACATTAAAGAAGG----------GGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
TATACATAAAAGAAGG----------GGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
TATACATTAAAGAAGG----------GGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
TATACATTAAAGAAGGGG--------GATGCGGATAAATGGAAAGGCGAAAGAAAGAATA""",
        )
        self.assertEqual(
            tuple(sequence.id for sequence in alignment.sequences),
            (
                "gi|6273291|gb|AF191665.1|AF191",
                "gi|6273290|gb|AF191664.1|AF191",
                "gi|6273289|gb|AF191663.1|AF191",
                "gi|6273287|gb|AF191661.1|AF191",
                "gi|6273286|gb|AF191660.1|AF191",
                "gi|6273285|gb|AF191659.1|AF191",
                "gi|6273284|gb|AF191658.1|AF191",
            ),
        )
        for i, sequence in enumerate(alignment.sequences[::-1]):
            sequence.id = "seq%d" % (i + 1)
        self.assertEqual(
            tuple(sequence.id for sequence in alignment.sequences),
            ("seq7", "seq6", "seq5", "seq4", "seq3", "seq2", "seq1"),
        )
        alignment.sort()
        self.assertEqual(
            "\n".join(row for row in alignment),  # str(alignment),
            """\
TATACATTAAAGAAGGGG--------GATGCGGATAAATGGAAAGGCGAAAGAAAGAATA
TATACATTAAAGAAGG----------GGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
TATACATAAAAGAAGG----------GGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
TATACATTAAAGAAGG----------GGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
TATACATTAAAGGAGGGGGA------TGCGGATAAATGGAAAGGCGAAAGAAAGAATATA
TATACATTAAAGGAGGGGGA------TGCGGATAAATGGAAAGGCGAAAGAAAGAATATA
TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA""",
        )
        self.assertEqual(
            tuple(sequence.id for sequence in alignment.sequences),
            ("seq1", "seq2", "seq3", "seq4", "seq5", "seq6", "seq7"),
        )
        alignment.sort(reverse=True)
        self.assertEqual(
            "\n".join(row for row in alignment),  # str(alignment),
            """\
TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA
TATACATTAAAGGAGGGGGA------TGCGGATAAATGGAAAGGCGAAAGAAAGAATATA
TATACATTAAAGGAGGGGGA------TGCGGATAAATGGAAAGGCGAAAGAAAGAATATA
TATACATTAAAGAAGG----------GGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
TATACATAAAAGAAGG----------GGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
TATACATTAAAGAAGG----------GGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
TATACATTAAAGAAGGGG--------GATGCGGATAAATGGAAAGGCGAAAGAAAGAATA""",
        )
        self.assertEqual(
            tuple(sequence.id for sequence in alignment.sequences),
            ("seq7", "seq6", "seq5", "seq4", "seq3", "seq2", "seq1"),
        )
        alignment.sort(key=lambda record: GC(record.seq))
        self.assertEqual(
            "\n".join(row for row in alignment),  # str(alignment),
            """\
TATACATAAAAGAAGG----------GGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA
TATACATTAAAGAAGG----------GGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
TATACATTAAAGGAGGGGGA------TGCGGATAAATGGAAAGGCGAAAGAAAGAATATA
TATACATTAAAGAAGGGG--------GATGCGGATAAATGGAAAGGCGAAAGAAAGAATA
TATACATTAAAGGAGGGGGA------TGCGGATAAATGGAAAGGCGAAAGAAAGAATATA
TATACATTAAAGAAGG----------GGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA""",
        )
        self.assertEqual(
            tuple(sequence.id for sequence in alignment.sequences),
            ("seq3", "seq7", "seq4", "seq5", "seq1", "seq6", "seq2"),
        )
        alignment.sort(key=lambda record: GC(record.seq), reverse=True)
        self.assertEqual(
            "\n".join(row for row in alignment),  # str(alignment),
            """\
TATACATTAAAGAAGG----------GGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
TATACATTAAAGGAGGGGGA------TGCGGATAAATGGAAAGGCGAAAGAAAGAATATA
TATACATTAAAGAAGGGG--------GATGCGGATAAATGGAAAGGCGAAAGAAAGAATA
TATACATTAAAGGAGGGGGA------TGCGGATAAATGGAAAGGCGAAAGAAAGAATATA
TATACATTAAAGAAGG----------GGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA
TATACATAAAAGAAGG----------GGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA""",
        )
        self.assertEqual(
            tuple(sequence.id for sequence in alignment.sequences),
            ("seq2", "seq6", "seq1", "seq5", "seq4", "seq7", "seq3"),
        )

    def test_substitutions(self):
        alignment = self.alignment
        m = alignment.substitutions
        self.assertEqual(
            str(m),
            """\
       A     C     G     T
A 1395.0   3.0  13.0   6.0
C    3.0 271.0   0.0  16.0
G    5.0   0.0 480.0   0.0
T    6.0  12.0   0.0 874.0
""",
        )
        self.assertAlmostEqual(m["T", "C"], 12.0)
        self.assertAlmostEqual(m["C", "T"], 16.0)
        m += m.transpose()
        m /= 2.0
        self.assertEqual(
            str(m),
            """\
       A     C     G     T
A 1395.0   3.0   9.0   6.0
C    3.0 271.0   0.0  14.0
G    9.0   0.0 480.0   0.0
T    6.0  14.0   0.0 874.0
""",
        )
        self.assertAlmostEqual(m["C", "T"], 14.0)
        self.assertAlmostEqual(m["T", "C"], 14.0)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
