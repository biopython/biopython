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
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction


class TestPairwiseAlignment(unittest.TestCase):
    target = "AACCGGGACCG"
    query = "ACGGAAC"
    query_rc = reverse_complement(query)
    forward_coordinates = numpy.array(
        [[0, 1, 2, 3, 4, 6, 7, 8, 8, 9, 11], [0, 1, 1, 2, 2, 4, 4, 5, 6, 7, 7]]
    )
    reverse_coordinates = numpy.array(
        [[0, 1, 2, 3, 4, 6, 7, 8, 8, 9, 11], [7, 6, 6, 5, 5, 3, 3, 2, 1, 0, 0]]
    )

    def check_indexing_slicing(self, alignment, cls, strand):
        msg = "%s, %s strand" % (cls.__name__, strand)
        self.assertEqual(
            repr(alignment),
            "<Alignment object (2 rows x 12 columns) at 0x%x>" % id(alignment),
        )
        if strand == "forward":
            self.assertEqual(
                str(alignment),
                """\
target            0 AACCGGGA-CCG 11
                  0 |-|-||-|-|-- 12
query             0 A-C-GG-AAC--  7
""",
                msg=msg,
            )
        if strand == "reverse":
            self.assertEqual(
                str(alignment),
                """\
target            0 AACCGGGA-CCG 11
                  0 |-|-||-|-|-- 12
query             7 A-C-GG-AAC--  0
""",
                msg=msg,
            )
        self.assertAlmostEqual(alignment.score, 6.0)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 12))
        self.assertIsInstance(alignment.sequences[0], cls)
        self.assertIsInstance(alignment.sequences[1], cls)
        self.assertEqual(alignment[0], "AACCGGGA-CCG", msg=msg)
        self.assertEqual(alignment[1], "A-C-GG-AAC--", msg=msg)
        self.assertEqual(alignment[-2], "AACCGGGA-CCG", msg=msg)
        self.assertEqual(alignment[-1], "A-C-GG-AAC--", msg=msg)
        self.assertEqual(alignment[0, 0], "A", msg=msg)
        self.assertEqual(alignment[0, 1], "A", msg=msg)
        self.assertEqual(alignment[0, 2], "C", msg=msg)
        self.assertEqual(alignment[0, 3], "C", msg=msg)
        self.assertEqual(alignment[0, 4], "G", msg=msg)
        self.assertEqual(alignment[0, 5], "G", msg=msg)
        self.assertEqual(alignment[0, 6], "G", msg=msg)
        self.assertEqual(alignment[0, 7], "A", msg=msg)
        self.assertEqual(alignment[0, 8], "-", msg=msg)
        self.assertEqual(alignment[0, 9], "C", msg=msg)
        self.assertEqual(alignment[0, 10], "C", msg=msg)
        self.assertEqual(alignment[0, 11], "G", msg=msg)
        self.assertEqual(alignment[1, 0], "A", msg=msg)
        self.assertEqual(alignment[1, 1], "-", msg=msg)
        self.assertEqual(alignment[1, 2], "C", msg=msg)
        self.assertEqual(alignment[1, 3], "-", msg=msg)
        self.assertEqual(alignment[1, 4], "G", msg=msg)
        self.assertEqual(alignment[1, 5], "G", msg=msg)
        self.assertEqual(alignment[1, 6], "-", msg=msg)
        self.assertEqual(alignment[1, 7], "A", msg=msg)
        self.assertEqual(alignment[1, 8], "A", msg=msg)
        self.assertEqual(alignment[1, 9], "C", msg=msg)
        self.assertEqual(alignment[1, 10], "-", msg=msg)
        self.assertEqual(alignment[1, 11], "-", msg=msg)
        self.assertEqual(alignment[0, :], "AACCGGGA-CCG", msg=msg)
        self.assertEqual(alignment[1, :], "A-C-GG-AAC--", msg=msg)
        self.assertEqual(alignment[-2, :], "AACCGGGA-CCG", msg=msg)
        self.assertEqual(alignment[-1, :], "A-C-GG-AAC--", msg=msg)
        self.assertEqual(alignment[0, 1:2], "A", msg=msg)
        self.assertEqual(alignment[1, 1:2], "-", msg=msg)
        self.assertEqual(alignment[0, 4:5], "G", msg=msg)
        self.assertEqual(alignment[1, 4:5], "G", msg=msg)
        self.assertEqual(alignment[0, 10:11], "C", msg=msg)
        self.assertEqual(alignment[1, 10:11], "-", msg=msg)
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
        subalignment = alignment[:, :]
        self.assertAlmostEqual(subalignment.score, 6.0, msg=msg)
        if strand == "forward":
            self.assertEqual(
                str(alignment[:, :]),
                """\
target            0 AACCGGGA-CCG 11
                  0 |-|-||-|-|-- 12
query             0 A-C-GG-AAC--  7
""",
                msg=msg,
            )
        if strand == "reverse":
            self.assertEqual(
                str(alignment[:, :]),
                """\
target            0 AACCGGGA-CCG 11
                  0 |-|-||-|-|-- 12
query             7 A-C-GG-AAC--  0
""",
                msg=msg,
            )
        self.assertIsInstance(subalignment.sequences[0], cls)
        self.assertIsInstance(subalignment.sequences[1], cls)
        self.assertEqual(alignment[0, 0:12], "AACCGGGA-CCG", msg=msg)
        self.assertEqual(alignment[1, 0:12], "A-C-GG-AAC--", msg=msg)
        self.assertEqual(alignment[0, 0:], "AACCGGGA-CCG", msg=msg)
        self.assertEqual(alignment[1, 0:], "A-C-GG-AAC--", msg=msg)
        subalignment = alignment[:, 0:]
        self.assertAlmostEqual(subalignment.score, 6.0, msg=msg)
        if strand == "forward":
            self.assertEqual(
                str(subalignment),
                """\
target            0 AACCGGGA-CCG 11
                  0 |-|-||-|-|-- 12
query             0 A-C-GG-AAC--  7
""",
                msg=msg,
            )
        if strand == "reverse":
            self.assertEqual(
                str(subalignment),
                """\
target            0 AACCGGGA-CCG 11
                  0 |-|-||-|-|-- 12
query             7 A-C-GG-AAC--  0
""",
                msg=msg,
            )
        self.assertIsInstance(subalignment.sequences[0], cls)
        self.assertIsInstance(subalignment.sequences[1], cls)
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
        subalignment = alignment[:, :12]
        self.assertAlmostEqual(subalignment.score, 6.0, msg=msg)
        if strand == "forward":
            self.assertEqual(
                str(subalignment),
                """\
target            0 AACCGGGA-CCG 11
                  0 |-|-||-|-|-- 12
query             0 A-C-GG-AAC--  7
""",
                msg=msg,
            )
        if strand == "reverse":
            self.assertEqual(
                str(subalignment),
                """\
target            0 AACCGGGA-CCG 11
                  0 |-|-||-|-|-- 12
query             7 A-C-GG-AAC--  0
""",
                msg=msg,
            )
        self.assertIsInstance(subalignment.sequences[0], cls)
        self.assertIsInstance(subalignment.sequences[1], cls)
        subalignment = alignment[:, 0:12]
        self.assertAlmostEqual(alignment.score, 6.0, msg=msg)
        if strand == "forward":
            self.assertEqual(
                str(subalignment),
                """\
target            0 AACCGGGA-CCG 11
                  0 |-|-||-|-|-- 12
query             0 A-C-GG-AAC--  7
""",
                msg=msg,
            )
        if strand == "reverse":
            self.assertEqual(
                str(subalignment),
                """\
target            0 AACCGGGA-CCG 11
                  0 |-|-||-|-|-- 12
query             7 A-C-GG-AAC--  0
""",
                msg=msg,
            )
        self.assertIsInstance(subalignment.sequences[0], cls)
        self.assertIsInstance(subalignment.sequences[1], cls)
        subalignment = alignment[:, 1:]
        if strand == "forward":
            self.assertEqual(
                str(subalignment),
                """\
target            1 ACCGGGA-CCG 11
                  0 -|-||-|-|-- 11
query             1 -C-GG-AAC--  7
""",
                msg=msg,
            )
        if strand == "reverse":
            self.assertEqual(
                str(subalignment),
                """\
target            1 ACCGGGA-CCG 11
                  0 -|-||-|-|-- 11
query             6 -C-GG-AAC--  0
""",
                msg=msg,
            )
        self.assertIsInstance(subalignment.sequences[0], cls)
        self.assertIsInstance(subalignment.sequences[1], cls)
        subalignment = alignment[:, 2:]
        if strand == "forward":
            self.assertEqual(
                str(subalignment),
                """\
target            2 CCGGGA-CCG 11
                  0 |-||-|-|-- 10
query             1 C-GG-AAC--  7
""",
                msg=msg,
            )
        if strand == "reverse":
            self.assertEqual(
                str(subalignment),
                """\
target            2 CCGGGA-CCG 11
                  0 |-||-|-|-- 10
query             6 C-GG-AAC--  0
""",
                msg=msg,
            )
        self.assertIsInstance(subalignment.sequences[0], cls)
        self.assertIsInstance(subalignment.sequences[1], cls)
        subalignment = alignment[:, 3:]
        if strand == "forward":
            self.assertEqual(
                str(subalignment),
                """\
target            3 CGGGA-CCG 11
                  0 -||-|-|--  9
query             2 -GG-AAC--  7
""",
                msg=msg,
            )
        if strand == "reverse":
            self.assertEqual(
                str(subalignment),
                """\
target            3 CGGGA-CCG 11
                  0 -||-|-|--  9
query             5 -GG-AAC--  0
""",
                msg=msg,
            )
        self.assertIsInstance(subalignment.sequences[0], cls)
        self.assertIsInstance(subalignment.sequences[1], cls)
        subalignment = alignment[:, 4:]
        if strand == "forward":
            self.assertEqual(
                str(subalignment),
                """\
target            4 GGGA-CCG 11
                  0 ||-|-|--  8
query             2 GG-AAC--  7
""",
                msg=msg,
            )
        if strand == "reverse":
            self.assertEqual(
                str(subalignment),
                """\
target            4 GGGA-CCG 11
                  0 ||-|-|--  8
query             5 GG-AAC--  0
""",
                msg=msg,
            )
        self.assertIsInstance(subalignment.sequences[0], cls)
        self.assertIsInstance(subalignment.sequences[1], cls)
        subalignment = alignment[:, :-1]
        if strand == "forward":
            self.assertEqual(
                str(subalignment),
                """\
target            0 AACCGGGA-CC 10
                  0 |-|-||-|-|- 11
query             0 A-C-GG-AAC-  7
""",
                msg=msg,
            )
        if strand == "reverse":
            self.assertEqual(
                str(subalignment),
                """\
target            0 AACCGGGA-CC 10
                  0 |-|-||-|-|- 11
query             7 A-C-GG-AAC-  0
""",
                msg=msg,
            )
        self.assertIsInstance(subalignment.sequences[0], cls)
        self.assertIsInstance(subalignment.sequences[1], cls)
        subalignment = alignment[:, :-2]
        if strand == "forward":
            self.assertEqual(
                str(subalignment),
                """\
target            0 AACCGGGA-C  9
                  0 |-|-||-|-| 10
query             0 A-C-GG-AAC  7
""",
                msg=msg,
            )
        if strand == "reverse":
            self.assertEqual(
                str(subalignment),
                """\
target            0 AACCGGGA-C  9
                  0 |-|-||-|-| 10
query             7 A-C-GG-AAC  0
""",
                msg=msg,
            )
        self.assertIsInstance(subalignment.sequences[0], cls)
        self.assertIsInstance(subalignment.sequences[1], cls)
        subalignment = alignment[:, :-3]
        if strand == "forward":
            self.assertEqual(
                str(subalignment),
                """\
target            0 AACCGGGA- 8
                  0 |-|-||-|- 9
query             0 A-C-GG-AA 6
""",
                msg=msg,
            )
        if strand == "reverse":
            self.assertEqual(
                str(subalignment),
                """\
target            0 AACCGGGA- 8
                  0 |-|-||-|- 9
query             7 A-C-GG-AA 1
""",
                msg=msg,
            )
        self.assertIsInstance(subalignment.sequences[0], cls)
        self.assertIsInstance(subalignment.sequences[1], cls)
        subalignment = alignment[:, 1:-1]
        if strand == "forward":
            self.assertEqual(
                str(subalignment),
                """\
target            1 ACCGGGA-CC 10
                  0 -|-||-|-|- 10
query             1 -C-GG-AAC-  7
""",
                msg=msg,
            )
        if strand == "reverse":
            self.assertEqual(
                str(subalignment),
                """\
target            1 ACCGGGA-CC 10
                  0 -|-||-|-|- 10
query             6 -C-GG-AAC-  0
""",
                msg=msg,
            )
        self.assertIsInstance(subalignment.sequences[0], cls)
        self.assertIsInstance(subalignment.sequences[1], cls)
        subalignment = alignment[:, 1:-2]
        if strand == "forward":
            self.assertEqual(
                str(subalignment),
                """\
target            1 ACCGGGA-C 9
                  0 -|-||-|-| 9
query             1 -C-GG-AAC 7
""",
                msg=msg,
            )
        if strand == "reverse":
            self.assertEqual(
                str(subalignment),
                """\
target            1 ACCGGGA-C 9
                  0 -|-||-|-| 9
query             6 -C-GG-AAC 0
""",
                msg=msg,
            )
        self.assertIsInstance(subalignment.sequences[0], cls)
        self.assertIsInstance(subalignment.sequences[1], cls)
        subalignment = alignment[:, 2:-1]
        if strand == "forward":
            self.assertEqual(
                str(subalignment),
                """\
target            2 CCGGGA-CC 10
                  0 |-||-|-|-  9
query             1 C-GG-AAC-  7
""",
                msg=msg,
            )
        if strand == "reverse":
            self.assertEqual(
                str(subalignment),
                """\
target            2 CCGGGA-CC 10
                  0 |-||-|-|-  9
query             6 C-GG-AAC-  0
""",
                msg=msg,
            )
        self.assertIsInstance(subalignment.sequences[0], cls)
        self.assertIsInstance(subalignment.sequences[1], cls)
        subalignment = alignment[:, 2:-2]
        if strand == "forward":
            self.assertEqual(
                str(subalignment),
                """\
target            2 CCGGGA-C 9
                  0 |-||-|-| 8
query             1 C-GG-AAC 7
""",
                msg=msg,
            )
        if strand == "reverse":
            self.assertEqual(
                str(subalignment),
                """\
target            2 CCGGGA-C 9
                  0 |-||-|-| 8
query             6 C-GG-AAC 0
""",
                msg=msg,
            )
        self.assertIsInstance(subalignment.sequences[0], cls)
        self.assertIsInstance(subalignment.sequences[1], cls)
        subalignment = alignment[:, ::2]
        self.assertEqual(
            str(subalignment),
            """\
target            0 ACGG-C 5
                  0 |||--- 6
query             0 ACG-A- 4
""",
            msg=msg,
        )
        self.assertIsInstance(subalignment.sequences[0], cls)
        self.assertIsInstance(subalignment.sequences[1], cls)
        subalignment = alignment[:, range(0, 12, 2)]
        self.assertEqual(
            str(subalignment),
            """\
target            0 ACGG-C 5
                  0 |||--- 6
query             0 ACG-A- 4
""",
            msg=msg,
        )
        self.assertIsInstance(subalignment.sequences[0], cls)
        self.assertIsInstance(subalignment.sequences[1], cls)
        subalignment = alignment[:, (1, 8, 5)]
        self.assertEqual(
            str(subalignment),
            """\
target            0 A-G 2
                  0 --| 3
query             0 -AG 2
""",
            msg=msg,
        )
        subalignment = alignment[:1]
        self.assertEqual(len(subalignment.sequences), 1)
        sequence = subalignment.sequences[0]
        self.assertIsInstance(sequence, cls)
        try:
            sequence = sequence.seq
        except AttributeError:
            pass
        self.assertEqual(sequence, "AACCGGGACCG")
        self.assertTrue(
            numpy.array_equal(
                subalignment.coordinates,
                numpy.array([[0, 1, 2, 3, 4, 6, 7, 8, 8, 9, 11]]),
            )
        )
        subalignment = alignment[:1, :]
        self.assertEqual(len(subalignment.sequences), 1)
        sequence = subalignment.sequences[0]
        self.assertIsInstance(sequence, cls)
        try:
            sequence = sequence.seq
        except AttributeError:
            pass
        self.assertEqual(sequence, "AACCGGGACCG")
        self.assertTrue(
            numpy.array_equal(
                subalignment.coordinates,
                numpy.array([[0, 1, 2, 3, 4, 6, 7, 8, 8, 9, 11]]),
            )
        )
        subalignment = alignment[:]
        self.assertEqual(alignment, subalignment)
        self.assertIsInstance(subalignment.sequences[0], cls)
        self.assertIsInstance(subalignment.sequences[1], cls)

    def test_indexing_slicing(self):
        sequences = (self.target, self.query)
        alignment = Align.Alignment(sequences, self.forward_coordinates)
        alignment.score = 6.0
        self.check_indexing_slicing(alignment, str, "forward")
        sequences = (self.target, self.query_rc)
        alignment = Align.Alignment(sequences, self.reverse_coordinates)
        alignment.score = 6.0
        self.check_indexing_slicing(alignment, str, "reverse")
        target = Seq(self.target)
        query = Seq(self.query)
        query_rc = Seq(self.query_rc)
        sequences = (target, query)
        alignment = Align.Alignment(sequences, self.forward_coordinates)
        alignment.score = 6.0
        self.check_indexing_slicing(alignment, Seq, "forward")
        sequences = (target, query_rc)
        alignment = Align.Alignment(sequences, self.reverse_coordinates)
        alignment.score = 6.0
        self.check_indexing_slicing(alignment, Seq, "reverse")
        target = SeqRecord(target, id=None)
        query = SeqRecord(query, id=None)
        query_rc = SeqRecord(query_rc, id=None)
        sequences = (target, query)
        alignment = Align.Alignment(sequences, self.forward_coordinates)
        alignment.score = 6.0
        self.check_indexing_slicing(alignment, SeqRecord, "forward")
        sequences = (target, query_rc)
        alignment = Align.Alignment(sequences, self.reverse_coordinates)
        alignment.score = 6.0
        self.check_indexing_slicing(alignment, SeqRecord, "reverse")

    def test_aligned_indices(self):
        sequences = (self.target, self.query)
        alignment = Align.Alignment(sequences, self.forward_coordinates)
        self.assertEqual(
            str(alignment),
            """\
target            0 AACCGGGA-CCG 11
                  0 |-|-||-|-|-- 12
query             0 A-C-GG-AAC--  7
""",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.aligned,
                # fmt: off
# flake8: noqa
                numpy.array([[[0, 1],
                              [2, 3],
                              [4, 6],
                              [7, 8],
                              [8, 9]],

                             [[0, 1],
                              [1, 2],
                              [2, 4],
                              [4, 5],
                              [6, 7]]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.indices,
                numpy.array(
                    [
                        # fmt: off
# flake8: noqa
                        [0,  1, 2,  3, 4, 5,  6, 7, -1, 8,  9, 10],
                        [0, -1, 1, -1, 2, 3, -1, 4,  5, 6, -1, -1],
                        # fmt: on
                    ]
                ),
            )
        )
        inverse_indices = alignment.inverse_indices
        self.assertEqual(len(inverse_indices), 2)
        self.assertTrue(
            # fmt: off
# flake8: noqa
            numpy.array_equal(
                inverse_indices[0],
                numpy.array([0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11])
            )
            # fmt: on
        )
        self.assertTrue(
            # fmt: off
# flake8: noqa
            numpy.array_equal(
                inverse_indices[1],
                numpy.array([0, 2, 4, 5, 7, 8, 9])
            )
            # fmt: on
        )
        alignment = Align.Alignment(sequences, self.forward_coordinates[:, 1:])
        self.assertEqual(
            str(alignment),
            """\
target            1 ACCGGGA-CCG 11
                  0 -|-||-|-|-- 11
query             1 -C-GG-AAC--  7
""",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.aligned,
                # fmt: off
# flake8: noqa
                numpy.array([[[2, 3],
                              [4, 6],
                              [7, 8],
                              [8, 9]],

                             [[1, 2],
                              [2, 4],
                              [4, 5],
                              [6, 7]]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.indices,
                numpy.array(
                    [
                        # fmt: off
# flake8: noqa
                        [ 1, 2,  3, 4, 5,  6, 7, -1, 8,  9, 10],
                        [-1, 1, -1, 2, 3, -1, 4,  5, 6, -1, -1],
                        # fmt: on
                    ]
                ),
            )
        )
        inverse_indices = alignment.inverse_indices
        self.assertEqual(len(inverse_indices), 2)
        self.assertTrue(
            # fmt: off
# flake8: noqa
            numpy.array_equal(
                inverse_indices[0],
                numpy.array([-1, 0, 1, 2, 3, 4, 5, 6, 8, 9, 10])
            )
            # fmt: on
        )
        self.assertTrue(
            # fmt: off
# flake8: noqa
            numpy.array_equal(
                inverse_indices[1],
                numpy.array([-1, 1, 3, 4, 6, 7, 8])
            )
            # fmt: on
        )
        alignment = Align.Alignment(sequences, self.forward_coordinates[:, :-1])
        self.assertEqual(
            str(alignment),
            """\
target            0 AACCGGGA-C  9
                  0 |-|-||-|-| 10
query             0 A-C-GG-AAC  7
""",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.aligned,
                # fmt: off
# flake8: noqa
                numpy.array([[[0, 1],
                              [2, 3],
                              [4, 6],
                              [7, 8],
                              [8, 9]],

                             [[0, 1],
                              [1, 2],
                              [2, 4],
                              [4, 5],
                              [6, 7]]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.indices,
                numpy.array(
                    [
                        # fmt: off
# flake8: noqa
                        [0,  1, 2,  3, 4, 5,  6, 7, -1, 8],
                        [0, -1, 1, -1, 2, 3, -1, 4,  5, 6],
                        # fmt: on
                    ]
                ),
            )
        )
        inverse_indices = alignment.inverse_indices
        self.assertEqual(len(inverse_indices), 2)
        self.assertTrue(
            # fmt: off
# flake8: noqa
            numpy.array_equal(
                inverse_indices[0],
                numpy.array([0, 1, 2, 3, 4, 5, 6, 7, 9, -1, -1])
            )
            # fmt: on
        )
        self.assertTrue(
            # fmt: off
# flake8: noqa
            numpy.array_equal(
                inverse_indices[1],
                numpy.array([0, 2, 4, 5, 7, 8, 9])
            )
            # fmt: on
        )
        alignment = Align.Alignment(sequences, self.forward_coordinates[:, 1:-1])
        self.assertEqual(
            str(alignment),
            """\
target            1 ACCGGGA-C 9
                  0 -|-||-|-| 9
query             1 -C-GG-AAC 7
""",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.aligned,
                # fmt: off
# flake8: noqa
                numpy.array([[[2, 3],
                              [4, 6],
                              [7, 8],
                              [8, 9]],

                             [[1, 2],
                              [2, 4],
                              [4, 5],
                              [6, 7]]])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.indices,
                numpy.array(
                    [
                        # fmt: off
# flake8: noqa
                        [ 1, 2,  3, 4, 5,  6, 7, -1, 8],
                        [-1, 1, -1, 2, 3, -1, 4,  5, 6],
                        # fmt: on
                    ]
                ),
            )
        )
        inverse_indices = alignment.inverse_indices
        self.assertEqual(len(inverse_indices), 2)
        self.assertTrue(
            # fmt: off
# flake8: noqa
            numpy.array_equal(
                inverse_indices[0],
                numpy.array([-1, 0, 1, 2, 3, 4, 5, 6, 8, -1, -1])
            )
            # fmt: on
        )
        self.assertTrue(
            # fmt: off
# flake8: noqa
            numpy.array_equal(
                inverse_indices[1],
                numpy.array([-1, 1, 3, 4, 6, 7, 8])
            )
            # fmt: on
        )
        sequences = (self.target, self.query_rc)
        alignment = Align.Alignment(sequences, self.reverse_coordinates)
        self.assertTrue(
            numpy.array_equal(
                alignment.aligned,
                # fmt: off
# flake8: noqa
                numpy.array([[[0, 1],
                              [2, 3],
                              [4, 6],
                              [7, 8],
                              [8, 9]],

                             [[7, 6],
                              [6, 5],
                              [5, 3],
                              [3, 2],
                              [1, 0]]])
                # fmt: on
            )
        )
        self.assertEqual(
            str(alignment),
            """\
target            0 AACCGGGA-CCG 11
                  0 |-|-||-|-|-- 12
query             7 A-C-GG-AAC--  0
""",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.indices,
                numpy.array(
                    [
                        # fmt: off
# flake8: noqa
                        [0,  1, 2,  3, 4, 5,  6, 7, -1, 8,  9, 10],
                        [6, -1, 5, -1, 4, 3, -1, 2,  1, 0, -1, -1],
                        # fmt: on
                    ]
                ),
            )
        )
        inverse_indices = alignment.inverse_indices
        self.assertEqual(len(inverse_indices), 2)
        self.assertTrue(
            # fmt: off
# flake8: noqa
            numpy.array_equal(
                inverse_indices[0],
                numpy.array([0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11])
            )
            # fmt: on
        )
        self.assertTrue(
            # fmt: off
# flake8: noqa
            numpy.array_equal(
                inverse_indices[1],
                numpy.array([9, 8, 7, 5, 4, 2, 0])
            )
            # fmt: on
        )
        alignment = Align.Alignment(sequences, self.reverse_coordinates[:, 1:])
        self.assertTrue(
            numpy.array_equal(
                alignment.aligned,
                # fmt: off
# flake8: noqa
                numpy.array([[[2, 3],
                              [4, 6],
                              [7, 8],
                              [8, 9]],

                             [[6, 5],
                              [5, 3],
                              [3, 2],
                              [1, 0]]])
                # fmt: on
            )
        )
        self.assertEqual(
            str(alignment),
            """\
target            1 ACCGGGA-CCG 11
                  0 -|-||-|-|-- 11
query             6 -C-GG-AAC--  0
""",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.indices,
                numpy.array(
                    [
                        # fmt: off
# flake8: noqa
                        [ 1, 2,  3, 4, 5,  6, 7, -1, 8,  9, 10],
                        [-1, 5, -1, 4, 3, -1, 2,  1, 0, -1, -1],
                        # fmt: on
                    ]
                ),
            )
        )
        inverse_indices = alignment.inverse_indices
        self.assertEqual(len(inverse_indices), 2)
        self.assertTrue(
            # fmt: off
# flake8: noqa
            numpy.array_equal(
                inverse_indices[0],
                numpy.array([-1, 0, 1, 2, 3, 4, 5, 6, 8, 9, 10])
            )
            # fmt: on
        )
        self.assertTrue(
            # fmt: off
# flake8: noqa
            numpy.array_equal(
                inverse_indices[1],
                numpy.array([8, 7, 6, 4, 3, 1, -1])
            )
            # fmt: on
        )
        alignment = Align.Alignment(sequences, self.reverse_coordinates[:, :-1])
        self.assertTrue(
            numpy.array_equal(
                alignment.aligned,
                # fmt: off
# flake8: noqa
                numpy.array([[[0, 1],
                              [2, 3],
                              [4, 6],
                              [7, 8],
                              [8, 9]],

                             [[7, 6],
                              [6, 5],
                              [5, 3],
                              [3, 2],
                              [1, 0]]])
                # fmt: on
            )
        )
        self.assertEqual(
            str(alignment),
            """\
target            0 AACCGGGA-C  9
                  0 |-|-||-|-| 10
query             7 A-C-GG-AAC  0
""",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.indices,
                numpy.array(
                    [
                        # fmt: off
# flake8: noqa
                        [0,  1, 2,  3, 4, 5,  6, 7, -1, 8],
                        [6, -1, 5, -1, 4, 3, -1, 2,  1, 0],
                        # fmt: on
                    ]
                ),
            )
        )
        inverse_indices = alignment.inverse_indices
        self.assertEqual(len(inverse_indices), 2)
        self.assertTrue(
            # fmt: off
# flake8: noqa
            numpy.array_equal(
                inverse_indices[0],
                numpy.array([0, 1, 2, 3, 4, 5, 6, 7, 9, -1, -1])
            )
            # fmt: on
        )
        self.assertTrue(
            # fmt: off
# flake8: noqa
            numpy.array_equal(
                inverse_indices[1],
                numpy.array([9, 8, 7, 5, 4, 2, 0])
            )
            # fmt: on
        )
        alignment = Align.Alignment(sequences, self.reverse_coordinates[:, 1:-1])
        self.assertTrue(
            numpy.array_equal(
                alignment.aligned,
                # fmt: off
# flake8: noqa
                numpy.array([[[2, 3],
                              [4, 6],
                              [7, 8],
                              [8, 9]],

                             [[6, 5],
                              [5, 3],
                              [3, 2],
                              [1, 0]]])
                # fmt: on
            )
        )
        self.assertEqual(
            str(alignment),
            """\
target            1 ACCGGGA-C 9
                  0 -|-||-|-| 9
query             6 -C-GG-AAC 0
""",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.indices,
                numpy.array(
                    [
                        # fmt: off
# flake8: noqa
                        [ 1, 2,  3, 4, 5,  6, 7, -1, 8],
                        [-1, 5, -1, 4, 3, -1, 2,  1, 0],
                        # fmt: on
                    ]
                ),
            )
        )
        inverse_indices = alignment.inverse_indices
        self.assertEqual(len(inverse_indices), 2)
        self.assertTrue(
            # fmt: off
# flake8: noqa
            numpy.array_equal(
                inverse_indices[0],
                numpy.array([-1, 0, 1, 2, 3, 4, 5, 6, 8, -1, -1])
            )
            # fmt: on
        )
        self.assertTrue(
            # fmt: off
# flake8: noqa
            numpy.array_equal(
                inverse_indices[1],
                numpy.array([8, 7, 6, 4, 3, 1, -1])
            )
            # fmt: on
        )

    def test_sort(self):
        target = Seq("ACTT")
        query = Seq("ACCT")
        sequences = (target, query)
        coordinates = numpy.array([[0, 4], [0, 4]])
        alignment = Align.Alignment(sequences, coordinates)
        self.assertEqual(
            str(alignment),
            """\
target            0 ACTT 4
                  0 ||.| 4
query             0 ACCT 4
""",
        )
        alignment.sort()
        self.assertEqual(
            str(alignment),
            """\
target            0 ACCT 4
                  0 ||.| 4
query             0 ACTT 4
""",
        )
        alignment.sort(reverse=True)
        self.assertEqual(
            str(alignment),
            """\
target            0 ACTT 4
                  0 ||.| 4
query             0 ACCT 4
""",
        )
        target.id = "seq1"
        query.id = "seq2"
        alignment.sort()
        self.assertEqual(
            str(alignment),
            """\
seq1              0 ACTT 4
                  0 ||.| 4
seq2              0 ACCT 4
""",
        )
        alignment.sort(reverse=True)
        self.assertEqual(
            str(alignment),
            """\
seq2              0 ACCT 4
                  0 ||.| 4
seq1              0 ACTT 4
""",
        )
        alignment.sort(key=gc_fraction)
        self.assertEqual(
            str(alignment),
            """\
seq1              0 ACTT 4
                  0 ||.| 4
seq2              0 ACCT 4
""",
        )
        alignment.sort(key=gc_fraction, reverse=True)
        self.assertEqual(
            str(alignment),
            """\
seq2              0 ACCT 4
                  0 ||.| 4
seq1              0 ACTT 4
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
            # fmt: off
# flake8: noqa
            [
                [
                     503,  744,  744,  747,  748,  820,  820,  822,  822,  823,
                     823,  828,  828,  833,  833,  845,  848,  850,  851,  854,
                     857, 1003, 1004, 1011, 1011, 1017, 1017, 1020, 1021, 1116,
                    1116, 1119, 1120, 1132, 1133, 1242, 1243, 1246, 1246, 1289,
                    1289, 1292, 1293, 1413,
                ],
                [
                     512,  753,  754,  757,  757,  829,  831,  833,  834,  835,
                     838,  843,  844,  849,  850,  862,  862,  864,  864,  867,
                     867, 1013, 1013, 1020, 1021, 1027, 1028, 1031, 1031, 1126,
                    1127, 1130, 1130, 1142, 1142, 1251, 1251, 1254, 1255, 1298,
                    1299, 1302, 1302, 1422,
                ],
            ]
            # fmt: on
        )
        sequences = (target, query)
        forward_alignment = Align.Alignment(sequences, coordinates)
        sequences = (target, query.reverse_complement())
        coordinates = coordinates.copy()
        coordinates[1, :] = len(query) - coordinates[1, :]
        reverse_alignment = Align.Alignment(sequences, coordinates)
        for alignment in (forward_alignment, reverse_alignment):
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
        path = "Clustalw/opuntia.aln"
        with open(path) as stream:
            alignments = Align.parse(stream, "clustal")
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

    def check_indexing_slicing(self, alignment, strand):
        msg = "%s strand" % strand
        self.assertEqual(
            repr(alignment),
            "<Alignment object (7 rows x 156 columns) at 0x%x>" % id(alignment),
        )
        if strand == "forward":
            self.assertEqual(
                str(alignment),
                """\
gi|627328         0 TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----
gi|627328         0 TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--
gi|627328         0 TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----
gi|627328         0 TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----
gi|627329         0 TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA
gi|627328         0 TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA
gi|627329         0 TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA

gi|627328        56 ------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGA
gi|627328        58 ------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGA
gi|627328        56 ------ATATATTTCAAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
gi|627328        56 ------ATATATTTATAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
gi|627329        60 ------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
gi|627328        60 ------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
gi|627329        60 TATATAATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA

gi|627328       110 TGAATATCAAAGAATCCATTGATTTAGTGTACCAGA 146
gi|627328       112 TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA 148
gi|627328       110 TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA 146
gi|627328       110 TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA 146
gi|627329       114 TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA 150
gi|627328       114 TGAATATCAAAGAATCTATTGATTTAGTATACCAGA 150
gi|627329       120 TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA 156
""",
                msg=msg,
            )
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
        if strand == "reverse":
            self.assertEqual(
                str(alignment),
                """\
gi|627328         0 TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----
gi|627328         0 TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--
gi|627328       146 TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----
gi|627328         0 TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----
gi|627329         0 TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA
gi|627328         0 TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA
gi|627329         0 TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA

gi|627328        56 ------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGA
gi|627328        58 ------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGA
gi|627328        90 ------ATATATTTCAAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
gi|627328        56 ------ATATATTTATAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
gi|627329        60 ------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
gi|627328        60 ------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
gi|627329        60 TATATAATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA

gi|627328       110 TGAATATCAAAGAATCCATTGATTTAGTGTACCAGA 146
gi|627328       112 TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA 148
gi|627328        36 TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA   0
gi|627328       110 TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA 146
gi|627329       114 TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA 150
gi|627328       114 TGAATATCAAAGAATCTATTGATTTAGTATACCAGA 150
gi|627329       120 TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA 156
""",
                msg=msg,
            )
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
        if strand == "forward":
            self.assertEqual(
                str(alignment),
                """\
gi|627328         0 TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----
gi|627328         0 TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--
gi|627328         0 TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----
gi|627328         0 TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----
gi|627329         0 TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA
gi|627328         0 TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA
gi|627329         0 TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA

gi|627328        56 ------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGA
gi|627328        58 ------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGA
gi|627328        56 ------ATATATTTCAAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
gi|627328        56 ------ATATATTTATAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
gi|627329        60 ------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
gi|627328        60 ------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
gi|627329        60 TATATAATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA

gi|627328       110 TGAATATCAAAGAATCCATTGATTTAGTGTACCAGA 146
gi|627328       112 TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA 148
gi|627328       110 TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA 146
gi|627328       110 TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA 146
gi|627329       114 TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA 150
gi|627328       114 TGAATATCAAAGAATCTATTGATTTAGTATACCAGA 150
gi|627329       120 TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA 156
""",
                msg=msg,
            )
        if strand == "reverse":
            self.assertEqual(
                str(alignment),
                """\
gi|627328         0 TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----
gi|627328         0 TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--
gi|627328       146 TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----
gi|627328         0 TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----
gi|627329         0 TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA
gi|627328         0 TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA
gi|627329         0 TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA

gi|627328        56 ------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGA
gi|627328        58 ------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGA
gi|627328        90 ------ATATATTTCAAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
gi|627328        56 ------ATATATTTATAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
gi|627329        60 ------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
gi|627328        60 ------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
gi|627329        60 TATATAATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA

gi|627328       110 TGAATATCAAAGAATCCATTGATTTAGTGTACCAGA 146
gi|627328       112 TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA 148
gi|627328        36 TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA   0
gi|627328       110 TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA 146
gi|627329       114 TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA 150
gi|627328       114 TGAATATCAAAGAATCTATTGATTTAGTATACCAGA 150
gi|627329       120 TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA 156
""",
                msg=msg,
            )
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
        if strand == "forward":
            self.assertEqual(
                str(alignment),
                """\
gi|627328         0 TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----
gi|627328         0 TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--
gi|627328         0 TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----
gi|627328         0 TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----
gi|627329         0 TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA
gi|627328         0 TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA
gi|627329         0 TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA

gi|627328        56 ------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGA
gi|627328        58 ------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGA
gi|627328        56 ------ATATATTTCAAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
gi|627328        56 ------ATATATTTATAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
gi|627329        60 ------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
gi|627328        60 ------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
gi|627329        60 TATATAATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA

gi|627328       110 TGAATATCAAAGAATCCATTGATTTAGTGTACCAGA 146
gi|627328       112 TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA 148
gi|627328       110 TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA 146
gi|627328       110 TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA 146
gi|627329       114 TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA 150
gi|627328       114 TGAATATCAAAGAATCTATTGATTTAGTATACCAGA 150
gi|627329       120 TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA 156
""",
                msg=msg,
            )
        if strand == "reverse":
            self.assertEqual(
                str(alignment),
                """\
gi|627328         0 TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----
gi|627328         0 TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--
gi|627328       146 TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----
gi|627328         0 TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----
gi|627329         0 TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA
gi|627328         0 TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA
gi|627329         0 TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA

gi|627328        56 ------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGA
gi|627328        58 ------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGA
gi|627328        90 ------ATATATTTCAAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
gi|627328        56 ------ATATATTTATAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
gi|627329        60 ------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
gi|627328        60 ------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
gi|627329        60 TATATAATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA

gi|627328       110 TGAATATCAAAGAATCCATTGATTTAGTGTACCAGA 146
gi|627328       112 TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA 148
gi|627328        36 TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA   0
gi|627328       110 TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA 146
gi|627329       114 TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA 150
gi|627328       114 TGAATATCAAAGAATCTATTGATTTAGTATACCAGA 150
gi|627329       120 TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA 156
""",
                msg=msg,
            )
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
        subalignment = alignment[:, :156]
        self.assertEqual(alignment, subalignment, msg=msg)
        self.assertEqual(
            alignment.column_annotations, subalignment.column_annotations, msg=msg
        )
        subalignment = alignment[:, 0:156]
        self.assertEqual(alignment, subalignment, msg=msg)
        self.assertEqual(
            alignment.column_annotations, subalignment.column_annotations, msg=msg
        )
        self.assertEqual(len(subalignment.column_annotations), 1)
        self.assertEqual(
            subalignment.column_annotations["clustal_consensus"],
            "******* **** *******************************************          ********  **** ********* ********************************************* *********** *******",
            msg=msg,
        )
        subalignment = alignment[:, 60:]
        if strand == "forward":
            self.assertEqual(
                str(subalignment),
                """\
gi|627328        56 ------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGA
gi|627328        58 ------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGA
gi|627328        56 ------ATATATTTCAAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
gi|627328        56 ------ATATATTTATAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
gi|627329        60 ------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
gi|627328        60 ------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
gi|627329        60 TATATAATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA

gi|627328       110 TGAATATCAAAGAATCCATTGATTTAGTGTACCAGA 146
gi|627328       112 TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA 148
gi|627328       110 TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA 146
gi|627328       110 TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA 146
gi|627329       114 TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA 150
gi|627328       114 TGAATATCAAAGAATCTATTGATTTAGTATACCAGA 150
gi|627329       120 TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA 156
""",
                msg=msg,
            )
        if strand == "reverse":
            self.assertEqual(
                str(subalignment),
                """\
gi|627328        56 ------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGA
gi|627328        58 ------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGA
gi|627328        90 ------ATATATTTCAAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
gi|627328        56 ------ATATATTTATAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
gi|627329        60 ------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
gi|627328        60 ------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA
gi|627329        60 TATATAATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGA

gi|627328       110 TGAATATCAAAGAATCCATTGATTTAGTGTACCAGA 146
gi|627328       112 TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA 148
gi|627328        36 TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA   0
gi|627328       110 TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA 146
gi|627329       114 TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA 150
gi|627328       114 TGAATATCAAAGAATCTATTGATTTAGTATACCAGA 150
gi|627329       120 TGAATATCAAAGAATCTATTGATTTAGTGTACCAGA 156
""",
                msg=msg,
            )
        self.assertEqual(len(subalignment.column_annotations), 1)
        self.assertEqual(
            subalignment.column_annotations["clustal_consensus"],
            "      ********  **** ********* ********************************************* *********** *******",
            msg=msg,
        )
        subalignment = alignment[:, :-60]
        if strand == "forward":
            self.assertEqual(
                str(subalignment),
                """\
gi|627328         0 TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----
gi|627328         0 TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--
gi|627328         0 TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----
gi|627328         0 TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----
gi|627329         0 TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA
gi|627328         0 TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA
gi|627329         0 TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA

gi|627328        56 ------ATATATTTCAAATTTCCTTATATACCCAAA 86
gi|627328        58 ------ATATATTTCAAATTTCCTTATATACCCAAA 88
gi|627328        56 ------ATATATTTCAAATTTCCTTATATATCCAAA 86
gi|627328        56 ------ATATATTTATAATTTCCTTATATATCCAAA 86
gi|627329        60 ------ATATATTTCAAATTCCCTTATATATCCAAA 90
gi|627328        60 ------ATATATTTCAAATTCCCTTATATATCCAAA 90
gi|627329        60 TATATAATATATTTCAAATTCCCTTATATATCCAAA 96
""",
                msg=msg,
            )
        if strand == "reverse":
            self.assertEqual(
                str(subalignment),
                """\
gi|627328         0 TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----
gi|627328         0 TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--
gi|627328       146 TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----
gi|627328         0 TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----
gi|627329         0 TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA
gi|627328         0 TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA
gi|627329         0 TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA

gi|627328        56 ------ATATATTTCAAATTTCCTTATATACCCAAA 86
gi|627328        58 ------ATATATTTCAAATTTCCTTATATACCCAAA 88
gi|627328        90 ------ATATATTTCAAATTTCCTTATATATCCAAA 60
gi|627328        56 ------ATATATTTATAATTTCCTTATATATCCAAA 86
gi|627329        60 ------ATATATTTCAAATTCCCTTATATATCCAAA 90
gi|627328        60 ------ATATATTTCAAATTCCCTTATATATCCAAA 90
gi|627329        60 TATATAATATATTTCAAATTCCCTTATATATCCAAA 96
""",
                msg=msg,
            )
        self.assertEqual(len(subalignment.column_annotations), 1)
        self.assertEqual(
            subalignment.column_annotations["clustal_consensus"],
            "******* **** *******************************************          ********  **** ********* *****",
            msg=msg,
        )
        subalignment = alignment[:, 20:-60]
        if strand == "forward":
            self.assertEqual(
                str(subalignment),
                """\
gi|627328        20 TGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTCAAATT
gi|627328        20 TGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--------ATATATTTCAAATT
gi|627328        20 TGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTCAAATT
gi|627328        20 TGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTATAATT
gi|627329        20 TGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA------ATATATTTCAAATT
gi|627328        20 TGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA------ATATATTTCAAATT
gi|627329        20 TGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATATATATAATATATTTCAAATT

gi|627328        70 TCCTTATATACCCAAA 86
gi|627328        72 TCCTTATATACCCAAA 88
gi|627328        70 TCCTTATATATCCAAA 86
gi|627328        70 TCCTTATATATCCAAA 86
gi|627329        74 CCCTTATATATCCAAA 90
gi|627328        74 CCCTTATATATCCAAA 90
gi|627329        80 CCCTTATATATCCAAA 96
""",
                msg=msg,
            )
        if strand == "reverse":
            self.assertEqual(
                str(subalignment),
                """\
gi|627328        20 TGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTCAAATT
gi|627328        20 TGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--------ATATATTTCAAATT
gi|627328       126 TGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTCAAATT
gi|627328        20 TGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTATAATT
gi|627329        20 TGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA------ATATATTTCAAATT
gi|627328        20 TGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA------ATATATTTCAAATT
gi|627329        20 TGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATATATATAATATATTTCAAATT

gi|627328        70 TCCTTATATACCCAAA 86
gi|627328        72 TCCTTATATACCCAAA 88
gi|627328        76 TCCTTATATATCCAAA 60
gi|627328        70 TCCTTATATATCCAAA 86
gi|627329        74 CCCTTATATATCCAAA 90
gi|627328        74 CCCTTATATATCCAAA 90
gi|627329        80 CCCTTATATATCCAAA 96
""",
                msg=msg,
            )
        self.assertEqual(len(subalignment.column_annotations), 1)
        self.assertEqual(
            subalignment.column_annotations["clustal_consensus"],
            "************************************          ********  **** ********* *****",
            msg=msg,
        )
        subalignment = alignment[:, ::2]
        self.assertEqual(
            str(subalignment),
            """\
gi|627328         0 TTCTAAAGGGTCGTATGAGCAAAAATTT-----AAATCATTCTTTCCATTAATTTAAATG
gi|627328         0 TTCTAAAGGGTCGTATGAGCAAAAATTTT----AAATCATTCTTTCCATTAATTTAAATG
gi|627328         0 TTCTAAAGGGTCGTATGAGCAAAAATTT-----AAATCATTCTTTTCATTAATTTAAATG
gi|627328         0 TTCTAAAGGGTCGTATGAGCAAAAATTT-----AAATAATTCTTTTCATTAATTTAAATG
gi|627329         0 TTCTAAGGGGTCGTATGAGCAAAAATTTTT---AAATCATCCTTTTCATTAATTTAAATG
gi|627328         0 TTCTAAGGGGTCGTATGAGCAAAAATTTTT---AAATCATCCTTTTCATTAATTTAAATG
gi|627329         0 TTCTAAGGGGTCGTATGAGCAAAAATTTTTTTTAAATCATCCTTTTCATTAATTTAAATG

gi|627328        55 TATTAAATCTGTTGGACG 73
gi|627328        56 TATTAAATTTGTTGGACG 74
gi|627328        55 TATTAAATTTGTTGGACG 73
gi|627328        55 TATTAAATTTGTTGGACG 73
gi|627329        57 TATTAAATTTGTTGGACG 75
gi|627328        57 TATTAAATTTGTTGAACG 75
gi|627329        60 TATTAAATTTGTTGGACG 78
""",
            msg=msg,
        )
        self.assertEqual(len(subalignment.column_annotations), 1)
        self.assertEqual(
            subalignment.column_annotations["clustal_consensus"],
            "****** *********************     **** ** **** ********************** ***** ***",
            msg=msg,
        )
        subalignment = alignment[:, range(0, 156, 2)]
        self.assertEqual(len(subalignment.column_annotations), 1)
        self.assertEqual(
            subalignment.column_annotations["clustal_consensus"],
            "****** *********************     **** ** **** ********************** ***** ***",
            msg=msg,
        )
        subalignment = alignment[:, (1, 7, 5)]
        self.assertEqual(len(subalignment.column_annotations), 1)
        self.assertEqual(
            subalignment.column_annotations["clustal_consensus"], "* *", msg=msg
        )
        self.assertEqual(
            str(alignment[1::3]),
            """\
gi|627328         0 TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--
                  0 ||||||||||||.|||||||||||||||||||||||||||||||||||||||||||||--
gi|627329         0 TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA

gi|627328        58 ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATA
                 60 ||||||||||||||.|||||||||.|||||||||||||||||||||||||||||||||||
gi|627329        60 ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATA

gi|627328       118 TCAAAGAATCTATTGATTTAGTGTACCAGA 148
                120 |||||||||||||||||||||||||||||| 150
gi|627329       120 TCAAAGAATCTATTGATTTAGTGTACCAGA 150
""",
            msg=msg,
        )
        self.assertEqual(
            str(alignment[1::3, :]),
            """\
gi|627328         0 TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--
                  0 ||||||||||||.|||||||||||||||||||||||||||||||||||||||||||||--
gi|627329         0 TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA

gi|627328        58 ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATA
                 60 ||||||||||||||.|||||||||.|||||||||||||||||||||||||||||||||||
gi|627329        60 ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATA

gi|627328       118 TCAAAGAATCTATTGATTTAGTGTACCAGA 148
                120 |||||||||||||||||||||||||||||| 150
gi|627329       120 TCAAAGAATCTATTGATTTAGTGTACCAGA 150
""",
            msg=msg,
        )
        self.assertEqual(alignment, alignment[:])

    def test_indexing_slicing(self):
        alignment = self.alignment
        strand = "forward"
        self.check_indexing_slicing(alignment, strand)
        name = alignment.sequences[2].id
        alignment.sequences[2] = alignment.sequences[2].reverse_complement()
        alignment.sequences[2].id = name
        n = len(alignment.sequences[2])
        alignment.coordinates[2, :] = n - alignment.coordinates[2, :]
        strand = "reverse"
        self.check_indexing_slicing(alignment, strand)

    def test_sort(self):
        alignment = self.alignment[:, 40:100]
        self.assertEqual(
            str(alignment),
            """\
gi|627328        40 AAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATA
gi|627328        40 AAAGAAAGAATATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATA
gi|627328        40 AAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATATCCAAATATA
gi|627328        40 AAAGAAAGAATATATA----------ATATATTTATAATTTCCTTATATATCCAAATATA
gi|627329        40 AAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATA
gi|627328        40 AAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATA
gi|627329        40 AAAGAAAGAATATATATATATATATAATATATTTCAAATTCCCTTATATATCCAAATATA

gi|627328        90
gi|627328        92
gi|627328        90
gi|627328        90
gi|627329        94
gi|627328        94
gi|627329       100
""",
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
            str(alignment),
            """\
gi|627328        40 AAAGAAAGAATATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATA
gi|627328        40 AAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATA
gi|627328        40 AAAGAAAGAATATATA----------ATATATTTATAATTTCCTTATATATCCAAATATA
gi|627328        40 AAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATATCCAAATATA
gi|627328        40 AAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATA
gi|627329        40 AAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATA
gi|627329        40 AAAGAAAGAATATATATATATATATAATATATTTCAAATTCCCTTATATATCCAAATATA

gi|627328        92
gi|627328        90
gi|627328        90
gi|627328        90
gi|627328        94
gi|627329        94
gi|627329       100
""",
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
            str(alignment),
            """\
gi|627329        40 AAAGAAAGAATATATATATATATATAATATATTTCAAATTCCCTTATATATCCAAATATA
gi|627329        40 AAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATA
gi|627328        40 AAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATA
gi|627328        40 AAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATATCCAAATATA
gi|627328        40 AAAGAAAGAATATATA----------ATATATTTATAATTTCCTTATATATCCAAATATA
gi|627328        40 AAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATA
gi|627328        40 AAAGAAAGAATATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATA

gi|627329       100
gi|627329        94
gi|627328        94
gi|627328        90
gi|627328        90
gi|627328        90
gi|627328        92
""",
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
            str(alignment),
            """\
seq1             40 AAAGAAAGAATATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATA
seq2             40 AAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATA
seq3             40 AAAGAAAGAATATATA----------ATATATTTATAATTTCCTTATATATCCAAATATA
seq4             40 AAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATATCCAAATATA
seq5             40 AAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATA
seq6             40 AAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATA
seq7             40 AAAGAAAGAATATATATATATATATAATATATTTCAAATTCCCTTATATATCCAAATATA

seq1             92
seq2             90
seq3             90
seq4             90
seq5             94
seq6             94
seq7            100
""",
        )
        self.assertEqual(
            tuple(sequence.id for sequence in alignment.sequences),
            ("seq1", "seq2", "seq3", "seq4", "seq5", "seq6", "seq7"),
        )
        alignment.sort(reverse=True)
        self.assertEqual(
            str(alignment),
            """\
seq7             40 AAAGAAAGAATATATATATATATATAATATATTTCAAATTCCCTTATATATCCAAATATA
seq6             40 AAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATA
seq5             40 AAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATA
seq4             40 AAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATATCCAAATATA
seq3             40 AAAGAAAGAATATATA----------ATATATTTATAATTTCCTTATATATCCAAATATA
seq2             40 AAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATA
seq1             40 AAAGAAAGAATATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATA

seq7            100
seq6             94
seq5             94
seq4             90
seq3             90
seq2             90
seq1             92
""",
        )
        self.assertEqual(
            tuple(sequence.id for sequence in alignment.sequences),
            ("seq7", "seq6", "seq5", "seq4", "seq3", "seq2", "seq1"),
        )
        alignment.sort(key=lambda record: gc_fraction(record.seq))
        self.assertEqual(
            str(alignment),
            """\
seq3             40 AAAGAAAGAATATATA----------ATATATTTATAATTTCCTTATATATCCAAATATA
seq7             40 AAAGAAAGAATATATATATATATATAATATATTTCAAATTCCCTTATATATCCAAATATA
seq4             40 AAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATATCCAAATATA
seq5             40 AAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATA
seq1             40 AAAGAAAGAATATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATA
seq6             40 AAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATA
seq2             40 AAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATA

seq3             90
seq7            100
seq4             90
seq5             94
seq1             92
seq6             94
seq2             90
""",
        )
        self.assertEqual(
            tuple(sequence.id for sequence in alignment.sequences),
            ("seq3", "seq7", "seq4", "seq5", "seq1", "seq6", "seq2"),
        )
        alignment.sort(key=lambda record: gc_fraction(record.seq), reverse=True)
        self.assertEqual(
            str(alignment),
            """\
seq2             40 AAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATA
seq6             40 AAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATA
seq1             40 AAAGAAAGAATATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATA
seq5             40 AAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATA
seq4             40 AAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATATCCAAATATA
seq7             40 AAAGAAAGAATATATATATATATATAATATATTTCAAATTCCCTTATATATCCAAATATA
seq3             40 AAAGAAAGAATATATA----------ATATATTTATAATTTCCTTATATATCCAAATATA

seq2             90
seq6             94
seq1             92
seq5             94
seq4             90
seq7            100
seq3             90
""",
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


class TestAlignment_format(unittest.TestCase):
    def setUp(self):
        path = "Clustalw/muscle.a2m"
        with open(path) as stream:
            alignments = Align.parse(stream, "a2m")
            alignment = next(alignments)
        self.alignment = alignment[:2, :]

    def test_a2m(self):
        self.assertEqual(
            self.alignment.format("a2m"),
            """\
>Test1seq <unknown description>
.................................................................AGTTACAATAACTGACGAAGCTAAGTAGGCTACTAATTAACGTCATCAACCTAATACATAGCACTTAGAAAAAAGTGAAGTAAGAAAATATAAAATAATAAAAGGGTGGGTTATCAATTGATAGTGTAAATCATCGTATTCCGGTGATATACCCTACCACAAAAACTCAAACCGACTTGATTCAAATCATCTCAATAAATTAGCGCCAAAATAATGAAAAAAATAATAACAAACAAAAACAAACCAAAATAAGAAAAAACATTACGCAAAACATAATAATTTACTCTTCGTTATTGTATTAACAAATCAAAGAGCTGAATTTTGATCACCTGCTAATACTACTTTCTGTATTGATCCTATATCAACGTAAACAAAGATACTAATAATTAACTAAAAGTACGTTCATCGATCGTGTTCGTTGACGAAGAAGAGCTCTATCTCCGGCGGAGCAAAGAAAACGATCTGTCTCCGTCGTAACACACGGTCGCTAGAGAAACTTTGCTTCTTCGGCGCCGGTGGACACGTCAGCATCTCCGGTATCCTAGACTTCTTGGCTTTCGGGGTACAACAACCGCGTGGTGACGTCAGCACCGCTGCTGGGGATGGAGAGGGAACAGAGTT.
>AT3G20900.1-SEQ <unknown description>
atgaacaaagtagcgaggaagaacaaaacatcaggtgaacaaaaaaaaaactcaatccacatcaaAGTTACAATAACTGACGAAGCTAAGTAGGCTAGAAATTAAAGTCATCAACCTAATACATAGCACTTAGAAAAAAGTGAAGCAAGAAAATATAAAATAATAAAAGGGTGGGTTATCAATTGATAGTGTAAATCATAGTTGATTTTTGATATACCCTACCACAAAAACTCAAACCGACTTGATTCAAATCATCTCAAAAAACAAGCGCCAAAATAATGAAAAAAATAATAACAAAAACAAACAAACCAAAATAAGAAAAAACATTACGCAAAACATAATAATTTACTCTTCGTTATTGTATTAACAAATCAAAGAGATGAATTTTGATCACCTGCTAATACTACTTTCTGTATTGATCCTATATCAAAAAAAAAAAAGATACTAATAATTAACTAAAAGTACGTTCATCGATCGTGTGCGTTGACGAAGAAGAGCTCTATCTCCGGCGGAGCAAAGAAAACGATCTGTCTCCGTCGTAACACACAGTTTTTCGAGACCCTTTGCTTCTTCGGCGCCGGTGGACACGTCAGCATCTCCGGTATCCTAGACTTCTTGGCTTTCGGGGTACAACAACCGCCTGGTGACGTCAGCACCGCTGCTGGGGATGGAGAGGGAACAGAGTAg
""",
        )

    def test_bed(self):
        self.alignment.score = 100
        self.assertEqual(
            self.alignment.format("bed"),
            """\
Test1seq	0	621	AT3G20900.1-SEQ	100	+	0	621	0	6	213,23,30,9,172,174,	0,213,236,266,275,447,
""",
        )

    def test_clustal(self):
        self.assertEqual(
            self.alignment.format("clustal"),
            """\
Test1seq                            --------------------------------------------------
AT3G20900.1-SEQ                     ATGAACAAAGTAGCGAGGAAGAACAAAACATCAGGTGAACAAAAAAAAAA

Test1seq                            ---------------AGTTACAATAACTGACGAAGCTAAGTAGGCTACTA
AT3G20900.1-SEQ                     CTCAATCCACATCAAAGTTACAATAACTGACGAAGCTAAGTAGGCTAGAA

Test1seq                            ATTAACGTCATCAACCTAATACATAGCACTTAGAAAAAAGTGAAGTAAGA
AT3G20900.1-SEQ                     ATTAAAGTCATCAACCTAATACATAGCACTTAGAAAAAAGTGAAGCAAGA

Test1seq                            AAATATAAAATAATAAAAGGGTGGGTTATCAATTGATAGTGTAAATCATC
AT3G20900.1-SEQ                     AAATATAAAATAATAAAAGGGTGGGTTATCAATTGATAGTGTAAATCATA

Test1seq                            GTATTCCGGTGATATACCCTACCACAAAAACTCAAACCGACTTGATTCAA
AT3G20900.1-SEQ                     GTTGATTTTTGATATACCCTACCACAAAAACTCAAACCGACTTGATTCAA

Test1seq                            ATCATCTCAATAAATTAGCGCCAAAATAATGAAAAAAATAATAACAAACA
AT3G20900.1-SEQ                     ATCATCTCAAAAAACAAGCGCCAAAATAATGAAAAAAATAATAACAAAAA

Test1seq                            AAAACAAACCAAAATAAGAAAAAACATTACGCAAAACATAATAATTTACT
AT3G20900.1-SEQ                     CAAACAAACCAAAATAAGAAAAAACATTACGCAAAACATAATAATTTACT

Test1seq                            CTTCGTTATTGTATTAACAAATCAAAGAGCTGAATTTTGATCACCTGCTA
AT3G20900.1-SEQ                     CTTCGTTATTGTATTAACAAATCAAAGAGATGAATTTTGATCACCTGCTA

Test1seq                            ATACTACTTTCTGTATTGATCCTATATCAACGTAAACAAAGATACTAATA
AT3G20900.1-SEQ                     ATACTACTTTCTGTATTGATCCTATATCAAAAAAAAAAAAGATACTAATA

Test1seq                            ATTAACTAAAAGTACGTTCATCGATCGTGTTCGTTGACGAAGAAGAGCTC
AT3G20900.1-SEQ                     ATTAACTAAAAGTACGTTCATCGATCGTGTGCGTTGACGAAGAAGAGCTC

Test1seq                            TATCTCCGGCGGAGCAAAGAAAACGATCTGTCTCCGTCGTAACACACGGT
AT3G20900.1-SEQ                     TATCTCCGGCGGAGCAAAGAAAACGATCTGTCTCCGTCGTAACACACAGT

Test1seq                            CGCTAGAGAAACTTTGCTTCTTCGGCGCCGGTGGACACGTCAGCATCTCC
AT3G20900.1-SEQ                     TTTTCGAGACCCTTTGCTTCTTCGGCGCCGGTGGACACGTCAGCATCTCC

Test1seq                            GGTATCCTAGACTTCTTGGCTTTCGGGGTACAACAACCGCGTGGTGACGT
AT3G20900.1-SEQ                     GGTATCCTAGACTTCTTGGCTTTCGGGGTACAACAACCGCCTGGTGACGT

Test1seq                            CAGCACCGCTGCTGGGGATGGAGAGGGAACAGAGTT-
AT3G20900.1-SEQ                     CAGCACCGCTGCTGGGGATGGAGAGGGAACAGAGTAG


""",
        )

    def test_bigbed(self):
        self.assertRaisesRegex(
            ValueError,
            "Formatting alignments has not yet been implemented for the bigbed format",
            self.alignment.format,
            "bigbed",
        )

    def test_bigmaf(self):
        self.assertRaisesRegex(
            ValueError,
            "Formatting alignments has not yet been implemented for the bigmaf format",
            self.alignment.format,
            "bigmaf",
        )

    def test_bigpsl(self):
        self.assertRaisesRegex(
            ValueError,
            "Formatting alignments has not yet been implemented for the bigpsl format",
            self.alignment.format,
            "bigpsl",
        )

    def test_emboss(self):
        self.assertRaisesRegex(
            ValueError,
            "Formatting alignments has not yet been implemented for the emboss format",
            self.alignment.format,
            "emboss",
        )

    def test_exonerate(self):
        self.alignment.score = 100
        self.assertEqual(
            self.alignment.format("exonerate"),
            """\
vulgar: AT3G20900.1-SEQ 0 687 + Test1seq 0 621 + 100 G 65 0 M 213 213 M 23 23 M 30 30 M 9 9 M 172 172 M 174 174 G 1 0
""",
        )
        self.assertEqual(
            self.alignment.format("exonerate", "vulgar"),
            """\
vulgar: AT3G20900.1-SEQ 0 687 + Test1seq 0 621 + 100 G 65 0 M 213 213 M 23 23 M 30 30 M 9 9 M 172 172 M 174 174 G 1 0
""",
        )
        self.assertEqual(
            self.alignment.format("exonerate", "cigar"),
            """\
cigar: AT3G20900.1-SEQ 0 687 + Test1seq 0 621 + 100 I 65 M 213 M 23 M 30 M 9 M 172 M 174 I 1
""",
        )

    def test_fasta(self):
        self.assertEqual(
            self.alignment.format("fasta"),
            """\
>Test1seq <unknown description>
-----------------------------------------------------------------AGTTACAATAACTGACGAAGCTAAGTAGGCTACTAATTAACGTCATCAACCTAATACATAGCACTTAGAAAAAAGTGAAGTAAGAAAATATAAAATAATAAAAGGGTGGGTTATCAATTGATAGTGTAAATCATCGTATTCCGGTGATATACCCTACCACAAAAACTCAAACCGACTTGATTCAAATCATCTCAATAAATTAGCGCCAAAATAATGAAAAAAATAATAACAAACAAAAACAAACCAAAATAAGAAAAAACATTACGCAAAACATAATAATTTACTCTTCGTTATTGTATTAACAAATCAAAGAGCTGAATTTTGATCACCTGCTAATACTACTTTCTGTATTGATCCTATATCAACGTAAACAAAGATACTAATAATTAACTAAAAGTACGTTCATCGATCGTGTTCGTTGACGAAGAAGAGCTCTATCTCCGGCGGAGCAAAGAAAACGATCTGTCTCCGTCGTAACACACGGTCGCTAGAGAAACTTTGCTTCTTCGGCGCCGGTGGACACGTCAGCATCTCCGGTATCCTAGACTTCTTGGCTTTCGGGGTACAACAACCGCGTGGTGACGTCAGCACCGCTGCTGGGGATGGAGAGGGAACAGAGTT-
>AT3G20900.1-SEQ <unknown description>
ATGAACAAAGTAGCGAGGAAGAACAAAACATCAGGTGAACAAAAAAAAAACTCAATCCACATCAAAGTTACAATAACTGACGAAGCTAAGTAGGCTAGAAATTAAAGTCATCAACCTAATACATAGCACTTAGAAAAAAGTGAAGCAAGAAAATATAAAATAATAAAAGGGTGGGTTATCAATTGATAGTGTAAATCATAGTTGATTTTTGATATACCCTACCACAAAAACTCAAACCGACTTGATTCAAATCATCTCAAAAAACAAGCGCCAAAATAATGAAAAAAATAATAACAAAAACAAACAAACCAAAATAAGAAAAAACATTACGCAAAACATAATAATTTACTCTTCGTTATTGTATTAACAAATCAAAGAGATGAATTTTGATCACCTGCTAATACTACTTTCTGTATTGATCCTATATCAAAAAAAAAAAAGATACTAATAATTAACTAAAAGTACGTTCATCGATCGTGTGCGTTGACGAAGAAGAGCTCTATCTCCGGCGGAGCAAAGAAAACGATCTGTCTCCGTCGTAACACACAGTTTTTCGAGACCCTTTGCTTCTTCGGCGCCGGTGGACACGTCAGCATCTCCGGTATCCTAGACTTCTTGGCTTTCGGGGTACAACAACCGCCTGGTGACGTCAGCACCGCTGCTGGGGATGGAGAGGGAACAGAGTAG
""",
        )

    def test_hhr(self):
        self.assertRaisesRegex(
            ValueError,
            "Formatting alignments has not yet been implemented for the hhr format",
            self.alignment.format,
            "hhr",
        )

    def test_maf(self):
        self.alignment.score = 100
        self.assertEqual(
            self.alignment.format("maf"),
            """\
a score=100.000000
s Test1seq        0 621 + 621 -----------------------------------------------------------------AGTTACAATAACTGACGAAGCTAAGTAGGCTACTAATTAACGTCATCAACCTAATACATAGCACTTAGAAAAAAGTGAAGTAAGAAAATATAAAATAATAAAAGGGTGGGTTATCAATTGATAGTGTAAATCATCGTATTCCGGTGATATACCCTACCACAAAAACTCAAACCGACTTGATTCAAATCATCTCAATAAATTAGCGCCAAAATAATGAAAAAAATAATAACAAACAAAAACAAACCAAAATAAGAAAAAACATTACGCAAAACATAATAATTTACTCTTCGTTATTGTATTAACAAATCAAAGAGCTGAATTTTGATCACCTGCTAATACTACTTTCTGTATTGATCCTATATCAACGTAAACAAAGATACTAATAATTAACTAAAAGTACGTTCATCGATCGTGTTCGTTGACGAAGAAGAGCTCTATCTCCGGCGGAGCAAAGAAAACGATCTGTCTCCGTCGTAACACACGGTCGCTAGAGAAACTTTGCTTCTTCGGCGCCGGTGGACACGTCAGCATCTCCGGTATCCTAGACTTCTTGGCTTTCGGGGTACAACAACCGCGTGGTGACGTCAGCACCGCTGCTGGGGATGGAGAGGGAACAGAGTT-
s AT3G20900.1-SEQ 0 687 + 687 ATGAACAAAGTAGCGAGGAAGAACAAAACATCAGGTGAACAAAAAAAAAACTCAATCCACATCAAAGTTACAATAACTGACGAAGCTAAGTAGGCTAGAAATTAAAGTCATCAACCTAATACATAGCACTTAGAAAAAAGTGAAGCAAGAAAATATAAAATAATAAAAGGGTGGGTTATCAATTGATAGTGTAAATCATAGTTGATTTTTGATATACCCTACCACAAAAACTCAAACCGACTTGATTCAAATCATCTCAAAAAACAAGCGCCAAAATAATGAAAAAAATAATAACAAAAACAAACAAACCAAAATAAGAAAAAACATTACGCAAAACATAATAATTTACTCTTCGTTATTGTATTAACAAATCAAAGAGATGAATTTTGATCACCTGCTAATACTACTTTCTGTATTGATCCTATATCAAAAAAAAAAAAGATACTAATAATTAACTAAAAGTACGTTCATCGATCGTGTGCGTTGACGAAGAAGAGCTCTATCTCCGGCGGAGCAAAGAAAACGATCTGTCTCCGTCGTAACACACAGTTTTTCGAGACCCTTTGCTTCTTCGGCGCCGGTGGACACGTCAGCATCTCCGGTATCCTAGACTTCTTGGCTTTCGGGGTACAACAACCGCCTGGTGACGTCAGCACCGCTGCTGGGGATGGAGAGGGAACAGAGTAG

""",
        )

    def test_mauve(self):
        alignment = self.alignment
        metadata = {"File": "testfile.fa"}
        for index, record in enumerate(alignment.sequences):
            record.id = str(index + 1)
        self.assertEqual(
            alignment.format("mauve", metadata=metadata),
            """\
> 2:1-621 + testfile.fa
-----------------------------------------------------------------AGTTACAATAACTGACGAAGCTAAGTAGGCTACTAATTAACGTCATCAACCTAATACATAGCACTTAGAAAAAAGTGAAGTAAGAAAATATAAAATAATAAAAGGGTGGGTTATCAATTGATAGTGTAAATCATCGTATTCCGGTGATATACCCTACCACAAAAACTCAAACCGACTTGATTCAAATCATCTCAATAAATTAGCGCCAAAATAATGAAAAAAATAATAACAAACAAAAACAAACCAAAATAAGAAAAAACATTACGCAAAACATAATAATTTACTCTTCGTTATTGTATTAACAAATCAAAGAGCTGAATTTTGATCACCTGCTAATACTACTTTCTGTATTGATCCTATATCAACGTAAACAAAGATACTAATAATTAACTAAAAGTACGTTCATCGATCGTGTTCGTTGACGAAGAAGAGCTCTATCTCCGGCGGAGCAAAGAAAACGATCTGTCTCCGTCGTAACACACGGTCGCTAGAGAAACTTTGCTTCTTCGGCGCCGGTGGACACGTCAGCATCTCCGGTATCCTAGACTTCTTGGCTTTCGGGGTACAACAACCGCGTGGTGACGTCAGCACCGCTGCTGGGGATGGAGAGGGAACAGAGTT-
> 3:1-687 + testfile.fa
ATGAACAAAGTAGCGAGGAAGAACAAAACATCAGGTGAACAAAAAAAAAACTCAATCCACATCAAAGTTACAATAACTGACGAAGCTAAGTAGGCTAGAAATTAAAGTCATCAACCTAATACATAGCACTTAGAAAAAAGTGAAGCAAGAAAATATAAAATAATAAAAGGGTGGGTTATCAATTGATAGTGTAAATCATAGTTGATTTTTGATATACCCTACCACAAAAACTCAAACCGACTTGATTCAAATCATCTCAAAAAACAAGCGCCAAAATAATGAAAAAAATAATAACAAAAACAAACAAACCAAAATAAGAAAAAACATTACGCAAAACATAATAATTTACTCTTCGTTATTGTATTAACAAATCAAAGAGATGAATTTTGATCACCTGCTAATACTACTTTCTGTATTGATCCTATATCAAAAAAAAAAAAGATACTAATAATTAACTAAAAGTACGTTCATCGATCGTGTGCGTTGACGAAGAAGAGCTCTATCTCCGGCGGAGCAAAGAAAACGATCTGTCTCCGTCGTAACACACAGTTTTTCGAGACCCTTTGCTTCTTCGGCGCCGGTGGACACGTCAGCATCTCCGGTATCCTAGACTTCTTGGCTTTCGGGGTACAACAACCGCCTGGTGACGTCAGCACCGCTGCTGGGGATGGAGAGGGAACAGAGTAG
=
""",
        )

    def test_msf(self):
        self.assertRaisesRegex(
            ValueError,
            "Formatting alignments has not yet been implemented for the msf format",
            self.alignment.format,
            "msf",
        )

    def test_nexus(self):
        alignment = self.alignment
        for record in alignment.sequences:
            record.annotations["molecule_type"] = "DNA"
        self.assertEqual(
            self.alignment.format("nexus"),
            """\
#NEXUS
begin data;
dimensions ntax=2 nchar=687;
format datatype=dna missing=? gap=-;
matrix
Test1seq          -----------------------------------------------------------------AGTTACAATAACTGACGAAGCTAAGTAGGCTACTAATTAACGTCATCAACCTAATACATAGCACTTAGAAAAAAGTGAAGTAAGAAAATATAAAATAATAAAAGGGTGGGTTATCAATTGATAGTGTAAATCATCGTATTCCGGTGATATACCCTACCACAAAAACTCAAACCGACTTGATTCAAATCATCTCAATAAATTAGCGCCAAAATAATGAAAAAAATAATAACAAACAAAAACAAACCAAAATAAGAAAAAACATTACGCAAAACATAATAATTTACTCTTCGTTATTGTATTAACAAATCAAAGAGCTGAATTTTGATCACCTGCTAATACTACTTTCTGTATTGATCCTATATCAACGTAAACAAAGATACTAATAATTAACTAAAAGTACGTTCATCGATCGTGTTCGTTGACGAAGAAGAGCTCTATCTCCGGCGGAGCAAAGAAAACGATCTGTCTCCGTCGTAACACACGGTCGCTAGAGAAACTTTGCTTCTTCGGCGCCGGTGGACACGTCAGCATCTCCGGTATCCTAGACTTCTTGGCTTTCGGGGTACAACAACCGCGTGGTGACGTCAGCACCGCTGCTGGGGATGGAGAGGGAACAGAGTT-
'AT3G20900.1-SEQ' ATGAACAAAGTAGCGAGGAAGAACAAAACATCAGGTGAACAAAAAAAAAACTCAATCCACATCAAAGTTACAATAACTGACGAAGCTAAGTAGGCTAGAAATTAAAGTCATCAACCTAATACATAGCACTTAGAAAAAAGTGAAGCAAGAAAATATAAAATAATAAAAGGGTGGGTTATCAATTGATAGTGTAAATCATAGTTGATTTTTGATATACCCTACCACAAAAACTCAAACCGACTTGATTCAAATCATCTCAAAAAACAAGCGCCAAAATAATGAAAAAAATAATAACAAAAACAAACAAACCAAAATAAGAAAAAACATTACGCAAAACATAATAATTTACTCTTCGTTATTGTATTAACAAATCAAAGAGATGAATTTTGATCACCTGCTAATACTACTTTCTGTATTGATCCTATATCAAAAAAAAAAAAGATACTAATAATTAACTAAAAGTACGTTCATCGATCGTGTGCGTTGACGAAGAAGAGCTCTATCTCCGGCGGAGCAAAGAAAACGATCTGTCTCCGTCGTAACACACAGTTTTTCGAGACCCTTTGCTTCTTCGGCGCCGGTGGACACGTCAGCATCTCCGGTATCCTAGACTTCTTGGCTTTCGGGGTACAACAACCGCCTGGTGACGTCAGCACCGCTGCTGGGGATGGAGAGGGAACAGAGTAG
;
end;
""",
        )

    def test_phylip(self):
        self.assertEqual(
            self.alignment.format("phylip"),
            """\
2 687
Test1seq  -----------------------------------------------------------------AGTTACAATAACTGACGAAGCTAAGTAGGCTACTAATTAACGTCATCAACCTAATACATAGCACTTAGAAAAAAGTGAAGTAAGAAAATATAAAATAATAAAAGGGTGGGTTATCAATTGATAGTGTAAATCATCGTATTCCGGTGATATACCCTACCACAAAAACTCAAACCGACTTGATTCAAATCATCTCAATAAATTAGCGCCAAAATAATGAAAAAAATAATAACAAACAAAAACAAACCAAAATAAGAAAAAACATTACGCAAAACATAATAATTTACTCTTCGTTATTGTATTAACAAATCAAAGAGCTGAATTTTGATCACCTGCTAATACTACTTTCTGTATTGATCCTATATCAACGTAAACAAAGATACTAATAATTAACTAAAAGTACGTTCATCGATCGTGTTCGTTGACGAAGAAGAGCTCTATCTCCGGCGGAGCAAAGAAAACGATCTGTCTCCGTCGTAACACACGGTCGCTAGAGAAACTTTGCTTCTTCGGCGCCGGTGGACACGTCAGCATCTCCGGTATCCTAGACTTCTTGGCTTTCGGGGTACAACAACCGCGTGGTGACGTCAGCACCGCTGCTGGGGATGGAGAGGGAACAGAGTT-
AT3G20900.ATGAACAAAGTAGCGAGGAAGAACAAAACATCAGGTGAACAAAAAAAAAACTCAATCCACATCAAAGTTACAATAACTGACGAAGCTAAGTAGGCTAGAAATTAAAGTCATCAACCTAATACATAGCACTTAGAAAAAAGTGAAGCAAGAAAATATAAAATAATAAAAGGGTGGGTTATCAATTGATAGTGTAAATCATAGTTGATTTTTGATATACCCTACCACAAAAACTCAAACCGACTTGATTCAAATCATCTCAAAAAACAAGCGCCAAAATAATGAAAAAAATAATAACAAAAACAAACAAACCAAAATAAGAAAAAACATTACGCAAAACATAATAATTTACTCTTCGTTATTGTATTAACAAATCAAAGAGATGAATTTTGATCACCTGCTAATACTACTTTCTGTATTGATCCTATATCAAAAAAAAAAAAGATACTAATAATTAACTAAAAGTACGTTCATCGATCGTGTGCGTTGACGAAGAAGAGCTCTATCTCCGGCGGAGCAAAGAAAACGATCTGTCTCCGTCGTAACACACAGTTTTTCGAGACCCTTTGCTTCTTCGGCGCCGGTGGACACGTCAGCATCTCCGGTATCCTAGACTTCTTGGCTTTCGGGGTACAACAACCGCCTGGTGACGTCAGCACCGCTGCTGGGGATGGAGAGGGAACAGAGTAG
""",
        )

    def test_psl(self):
        self.assertEqual(
            self.alignment.format("psl"),
            """\
589	32	0	0	0	0	0	0	+	AT3G20900.1-SEQ	687	65	686	Test1seq	621	0	621	6	213,23,30,9,172,174,	65,278,301,331,340,512,	0,213,236,266,275,447,
""",
        )
        self.assertEqual(
            self.alignment.format("psl", mask="upper"),
            """\
0	32	589	0	0	0	0	0	+	AT3G20900.1-SEQ	687	65	686	Test1seq	621	0	621	6	213,23,30,9,172,174,	65,278,301,331,340,512,	0,213,236,266,275,447,
""",
        )
        self.assertEqual(
            self.alignment.format("psl", wildcard="A"),
            """\
362	13	0	246	0	0	0	0	+	AT3G20900.1-SEQ	687	65	686	Test1seq	621	0	621	6	213,23,30,9,172,174,	65,278,301,331,340,512,	0,213,236,266,275,447,
""",
        )

    def test_sam(self):
        self.alignment.score = 100
        self.assertEqual(
            self.alignment.format("sam"),
            """\
AT3G20900.1-SEQ	0	Test1seq	1	255	65I213M23M30M9M172M174M1I	*	0	0	ATGAACAAAGTAGCGAGGAAGAACAAAACATCAGGTGAACAAAAAAAAAACTCAATCCACATCAAAGTTACAATAACTGACGAAGCTAAGTAGGCTAGAAATTAAAGTCATCAACCTAATACATAGCACTTAGAAAAAAGTGAAGCAAGAAAATATAAAATAATAAAAGGGTGGGTTATCAATTGATAGTGTAAATCATAGTTGATTTTTGATATACCCTACCACAAAAACTCAAACCGACTTGATTCAAATCATCTCAAAAAACAAGCGCCAAAATAATGAAAAAAATAATAACAAAAACAAACAAACCAAAATAAGAAAAAACATTACGCAAAACATAATAATTTACTCTTCGTTATTGTATTAACAAATCAAAGAGATGAATTTTGATCACCTGCTAATACTACTTTCTGTATTGATCCTATATCAAAAAAAAAAAAGATACTAATAATTAACTAAAAGTACGTTCATCGATCGTGTGCGTTGACGAAGAAGAGCTCTATCTCCGGCGGAGCAAAGAAAACGATCTGTCTCCGTCGTAACACACAGTTTTTCGAGACCCTTTGCTTCTTCGGCGCCGGTGGACACGTCAGCATCTCCGGTATCCTAGACTTCTTGGCTTTCGGGGTACAACAACCGCCTGGTGACGTCAGCACCGCTGCTGGGGATGGAGAGGGAACAGAGTAG	*	AS:i:100
""",
        )
        self.assertEqual(
            self.alignment.format("sam", md=True),
            """\
AT3G20900.1-SEQ	0	Test1seq	1	255	65I213M23M30M9M172M174M1I	*	0	0	ATGAACAAAGTAGCGAGGAAGAACAAAACATCAGGTGAACAAAAAAAAAACTCAATCCACATCAAAGTTACAATAACTGACGAAGCTAAGTAGGCTAGAAATTAAAGTCATCAACCTAATACATAGCACTTAGAAAAAAGTGAAGCAAGAAAATATAAAATAATAAAAGGGTGGGTTATCAATTGATAGTGTAAATCATAGTTGATTTTTGATATACCCTACCACAAAAACTCAAACCGACTTGATTCAAATCATCTCAAAAAACAAGCGCCAAAATAATGAAAAAAATAATAACAAAAACAAACAAACCAAAATAAGAAAAAACATTACGCAAAACATAATAATTTACTCTTCGTTATTGTATTAACAAATCAAAGAGATGAATTTTGATCACCTGCTAATACTACTTTCTGTATTGATCCTATATCAAAAAAAAAAAAGATACTAATAATTAACTAAAAGTACGTTCATCGATCGTGTGCGTTGACGAAGAAGAGCTCTATCTCCGGCGGAGCAAAGAAAACGATCTGTCTCCGTCGTAACACACAGTTTTTCGAGACCCTTTGCTTCTTCGGCGCCGGTGGACACGTCAGCATCTCCGGTATCCTAGACTTCTTGGCTTTCGGGGTACAACAACCGCCTGGTGACGTCAGCACCGCTGCTGGGGATGGAGAGGGAACAGAGTAG	*	MD:Z:32C0T6C39T53C2A0T0T0C0C0G0G51T3T0T32C1A78C50C0G0T3C43T66G2C0G0C1A4A0A79G44T	AS:i:100
""",
        )

    def test_stockholm(self):
        alignment = self.alignment
        del alignment.column_annotations["state"]
        self.assertEqual(
            self.alignment.format("stockholm"),
            """\
# STOCKHOLM 1.0
#=GF SQ   2
Test1seq                        -----------------------------------------------------------------AGTTACAATAACTGACGAAGCTAAGTAGGCTACTAATTAACGTCATCAACCTAATACATAGCACTTAGAAAAAAGTGAAGTAAGAAAATATAAAATAATAAAAGGGTGGGTTATCAATTGATAGTGTAAATCATCGTATTCCGGTGATATACCCTACCACAAAAACTCAAACCGACTTGATTCAAATCATCTCAATAAATTAGCGCCAAAATAATGAAAAAAATAATAACAAACAAAAACAAACCAAAATAAGAAAAAACATTACGCAAAACATAATAATTTACTCTTCGTTATTGTATTAACAAATCAAAGAGCTGAATTTTGATCACCTGCTAATACTACTTTCTGTATTGATCCTATATCAACGTAAACAAAGATACTAATAATTAACTAAAAGTACGTTCATCGATCGTGTTCGTTGACGAAGAAGAGCTCTATCTCCGGCGGAGCAAAGAAAACGATCTGTCTCCGTCGTAACACACGGTCGCTAGAGAAACTTTGCTTCTTCGGCGCCGGTGGACACGTCAGCATCTCCGGTATCCTAGACTTCTTGGCTTTCGGGGTACAACAACCGCGTGGTGACGTCAGCACCGCTGCTGGGGATGGAGAGGGAACAGAGTT-
AT3G20900.1-SEQ                 ATGAACAAAGTAGCGAGGAAGAACAAAACATCAGGTGAACAAAAAAAAAACTCAATCCACATCAAAGTTACAATAACTGACGAAGCTAAGTAGGCTAGAAATTAAAGTCATCAACCTAATACATAGCACTTAGAAAAAAGTGAAGCAAGAAAATATAAAATAATAAAAGGGTGGGTTATCAATTGATAGTGTAAATCATAGTTGATTTTTGATATACCCTACCACAAAAACTCAAACCGACTTGATTCAAATCATCTCAAAAAACAAGCGCCAAAATAATGAAAAAAATAATAACAAAAACAAACAAACCAAAATAAGAAAAAACATTACGCAAAACATAATAATTTACTCTTCGTTATTGTATTAACAAATCAAAGAGATGAATTTTGATCACCTGCTAATACTACTTTCTGTATTGATCCTATATCAAAAAAAAAAAAGATACTAATAATTAACTAAAAGTACGTTCATCGATCGTGTGCGTTGACGAAGAAGAGCTCTATCTCCGGCGGAGCAAAGAAAACGATCTGTCTCCGTCGTAACACACAGTTTTTCGAGACCCTTTGCTTCTTCGGCGCCGGTGGACACGTCAGCATCTCCGGTATCCTAGACTTCTTGGCTTTCGGGGTACAACAACCGCCTGGTGACGTCAGCACCGCTGCTGGGGATGGAGAGGGAACAGAGTAG
//
""",
        )

    def test_tabular(self):
        self.assertRaisesRegex(
            ValueError,
            "Formatting alignments has not yet been implemented for the tabular format",
            self.alignment.format,
            "tabular",
        )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
