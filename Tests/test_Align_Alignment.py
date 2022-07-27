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
from Bio.SeqUtils import gc_content


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
        self.assertEqual(
            str(alignment[:, :]),
            """\
AACCGGGA-CCG
|-|-||-|-|--
A-C-GG-AAC--
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
        self.assertEqual(
            str(subalignment),
            """\
AACCGGGA-CCG
|-|-||-|-|--
A-C-GG-AAC--
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
        self.assertEqual(
            str(subalignment),
            """\
AACCGGGA-CCG
|-|-||-|-|--
A-C-GG-AAC--
""",
            msg=msg,
        )
        self.assertIsInstance(subalignment.sequences[0], cls)
        self.assertIsInstance(subalignment.sequences[1], cls)
        subalignment = alignment[:, 0:12]
        self.assertAlmostEqual(alignment.score, 6.0, msg=msg)
        self.assertEqual(
            str(subalignment),
            """\
AACCGGGA-CCG
|-|-||-|-|--
A-C-GG-AAC--
""",
            msg=msg,
        )
        self.assertIsInstance(subalignment.sequences[0], cls)
        self.assertIsInstance(subalignment.sequences[1], cls)
        subalignment = alignment[:, 1:]
        self.assertEqual(
            str(subalignment),
            """\
AACCGGGA-CCG
 -|-||-|-|--
A-C-GG-AAC--
""",
            msg=msg,
        )
        self.assertIsInstance(subalignment.sequences[0], cls)
        self.assertIsInstance(subalignment.sequences[1], cls)
        subalignment = alignment[:, 2:]
        self.assertEqual(
            str(subalignment),
            """\
AACCGGGA-CCG
  |-||-|-|--
 AC-GG-AAC--
""",
            msg=msg,
        )
        self.assertIsInstance(subalignment.sequences[0], cls)
        self.assertIsInstance(subalignment.sequences[1], cls)
        subalignment = alignment[:, 3:]
        self.assertEqual(
            str(subalignment),
            """\
AACCGGGA-CCG
   -||-|-|--
 AC-GG-AAC--
""",
            msg=msg,
        )
        self.assertIsInstance(subalignment.sequences[0], cls)
        self.assertIsInstance(subalignment.sequences[1], cls)
        subalignment = alignment[:, 4:]
        self.assertEqual(
            str(subalignment),
            """\
AACCGGGA-CCG
    ||-|-|--
  ACGG-AAC--
""",
            msg=msg,
        )
        self.assertIsInstance(subalignment.sequences[0], cls)
        self.assertIsInstance(subalignment.sequences[1], cls)
        subalignment = alignment[:, :-1]
        self.assertEqual(
            str(subalignment),
            """\
AACCGGGA-CCG
|-|-||-|-|-
A-C-GG-AAC-
""",
            msg=msg,
        )
        self.assertIsInstance(subalignment.sequences[0], cls)
        self.assertIsInstance(subalignment.sequences[1], cls)
        subalignment = alignment[:, :-2]
        self.assertEqual(
            str(subalignment),
            """\
AACCGGGA-CCG
|-|-||-|-|
A-C-GG-AAC
""",
            msg=msg,
        )
        self.assertIsInstance(subalignment.sequences[0], cls)
        self.assertIsInstance(subalignment.sequences[1], cls)
        subalignment = alignment[:, :-3]
        self.assertEqual(
            str(subalignment),
            """\
AACCGGGA-CCG
|-|-||-|-
A-C-GG-AAC
""",
            msg=msg,
        )
        self.assertIsInstance(subalignment.sequences[0], cls)
        self.assertIsInstance(subalignment.sequences[1], cls)
        subalignment = alignment[:, 1:-1]
        self.assertEqual(
            str(subalignment),
            """\
AACCGGGA-CCG
 -|-||-|-|-
A-C-GG-AAC-
""",
            msg=msg,
        )
        self.assertIsInstance(subalignment.sequences[0], cls)
        self.assertIsInstance(subalignment.sequences[1], cls)
        subalignment = alignment[:, 1:-2]
        self.assertEqual(
            str(subalignment),
            """\
AACCGGGA-CCG
 -|-||-|-|
A-C-GG-AAC
""",
            msg=msg,
        )
        self.assertIsInstance(subalignment.sequences[0], cls)
        self.assertIsInstance(subalignment.sequences[1], cls)
        subalignment = alignment[:, 2:-1]
        self.assertEqual(
            str(subalignment),
            """\
AACCGGGA-CCG
  |-||-|-|-
 AC-GG-AAC-
""",
            msg=msg,
        )
        self.assertIsInstance(subalignment.sequences[0], cls)
        self.assertIsInstance(subalignment.sequences[1], cls)
        subalignment = alignment[:, 2:-2]
        self.assertEqual(
            str(subalignment),
            """\
AACCGGGA-CCG
  |-||-|-|
 AC-GG-AAC
""",
            msg=msg,
        )
        self.assertIsInstance(subalignment.sequences[0], cls)
        self.assertIsInstance(subalignment.sequences[1], cls)
        subalignment = alignment[:, ::2]
        self.assertEqual(
            str(subalignment),
            """\
ACGG-C
|||---
ACG-A-
""",
            msg=msg,
        )
        self.assertIsInstance(subalignment.sequences[0], cls)
        self.assertIsInstance(subalignment.sequences[1], cls)
        subalignment = alignment[:, range(0, 12, 2)]
        self.assertEqual(
            str(subalignment),
            """\
ACGG-C
|||---
ACG-A-
""",
            msg=msg,
        )
        self.assertIsInstance(subalignment.sequences[0], cls)
        self.assertIsInstance(subalignment.sequences[1], cls)
        subalignment = alignment[:, (1, 8, 5)]
        self.assertEqual(
            str(subalignment),
            """\
A-G
--|
-AG
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
        target = SeqRecord(target)
        query = SeqRecord(query)
        query_rc = SeqRecord(query_rc)
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
AACCGGGA-CCG
|-|-||-|-|--
A-C-GG-AAC--
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
AACCGGGA-CCG
 -|-||-|-|--
A-C-GG-AAC--
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
AACCGGGA-CCG
|-|-||-|-|
A-C-GG-AAC
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
AACCGGGA-CCG
 -|-||-|-|
A-C-GG-AAC
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
AACCGGGA-CCG
|-|-||-|-|--
A-C-GG-AAC--
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
AACCGGGA-CCG
 -|-||-|-|--
A-C-GG-AAC--
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
AACCGGGA-CCG
|-|-||-|-|
A-C-GG-AAC
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
AACCGGGA-CCG
 -|-||-|-|
A-C-GG-AAC
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
        alignment.sort(key=gc_content)
        self.assertEqual(
            str(alignment),
            """\
ACTT
||.|
ACCT
""",
        )
        alignment.sort(key=gc_content, reverse=True)
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
        import warnings
        from Bio import BiopythonExperimentalWarning

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonExperimentalWarning)
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
            "<Alignment object (7 rows x 156 columns) at 0x%x>" % id(alignment),
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
        self.assertEqual(
            "\n".join(row for row in subalignment),  # str(subalignment),
            """\
------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCCATTGATTTAGTGTACCAGA
------ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA
------ATATATTTCAAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA
------ATATATTTATAATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA
------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA
------ATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTATACCAGA
TATATAATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA""",
            msg=msg,
        )
        self.assertEqual(len(subalignment.column_annotations), 1)
        self.assertEqual(
            subalignment.column_annotations["clustal_consensus"],
            "      ********  **** ********* ********************************************* *********** *******",
            msg=msg,
        )
        subalignment = alignment[:, :-60]
        self.assertEqual(
            "\n".join(row for row in subalignment),  # str(subalignment),
            """\
TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATACCCAAA
TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--------ATATATTTCAAATTTCCTTATATACCCAAA
TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATATCCAAA
TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTATAATTTCCTTATATATCCAAA
TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAA
TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAA
TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATATATATAATATATTTCAAATTCCCTTATATATCCAAA""",
            msg=msg,
        )
        self.assertEqual(len(subalignment.column_annotations), 1)
        self.assertEqual(
            subalignment.column_annotations["clustal_consensus"],
            "******* **** *******************************************          ********  **** ********* *****",
            msg=msg,
        )
        subalignment = alignment[:, 20:-60]
        self.assertEqual(
            "\n".join(row for row in subalignment),  # str(subalignment),
            """\
TGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATACCCAAA
TGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--------ATATATTTCAAATTTCCTTATATACCCAAA
TGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATATCCAAA
TGCGGATAAATGGAAAGGCGAAAGAAAGAATATATA----------ATATATTTATAATTTCCTTATATATCCAAA
TGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAA
TGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAA
TGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATATATATAATATATTTCAAATTCCCTTATATATCCAAA""",
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
            "\n".join(row for row in subalignment),  # str(subalignment),
            """\
TTCTAAAGGGTCGTATGAGCAAAAATTT-----AAATCATTCTTTCCATTAATTTAAATGTATTAAATCTGTTGGACG
TTCTAAAGGGTCGTATGAGCAAAAATTTT----AAATCATTCTTTCCATTAATTTAAATGTATTAAATTTGTTGGACG
TTCTAAAGGGTCGTATGAGCAAAAATTT-----AAATCATTCTTTTCATTAATTTAAATGTATTAAATTTGTTGGACG
TTCTAAAGGGTCGTATGAGCAAAAATTT-----AAATAATTCTTTTCATTAATTTAAATGTATTAAATTTGTTGGACG
TTCTAAGGGGTCGTATGAGCAAAAATTTTT---AAATCATCCTTTTCATTAATTTAAATGTATTAAATTTGTTGGACG
TTCTAAGGGGTCGTATGAGCAAAAATTTTT---AAATCATCCTTTTCATTAATTTAAATGTATTAAATTTGTTGAACG
TTCTAAGGGGTCGTATGAGCAAAAATTTTTTTTAAATCATCCTTTTCATTAATTTAAATGTATTAAATTTGTTGGACG""",
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
TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA
||||||||||||.|||||||||||||||||||||||||||||||||||||||||||||--||||||||||||||.|||||||||.|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATAATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA
""",
            msg=msg,
        )
        self.assertEqual(
            str(alignment[1::3, :]),
            """\
TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATA--ATATATTTCAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA
||||||||||||.|||||||||||||||||||||||||||||||||||||||||||||--||||||||||||||.|||||||||.|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATAATATATTTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA
""",
            msg=msg,
        )
        self.assertEqual(alignment, alignment[:])

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
AAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATA
AAAGAAAGAATATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATA
AAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATATCCAAATATA
AAAGAAAGAATATATA----------ATATATTTATAATTTCCTTATATATCCAAATATA
AAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATA
AAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATA
AAAGAAAGAATATATATATATATATAATATATTTCAAATTCCCTTATATATCCAAATATA""",
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
AAAGAAAGAATATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATA
AAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATA
AAAGAAAGAATATATA----------ATATATTTATAATTTCCTTATATATCCAAATATA
AAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATATCCAAATATA
AAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATA
AAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATA
AAAGAAAGAATATATATATATATATAATATATTTCAAATTCCCTTATATATCCAAATATA""",
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
AAAGAAAGAATATATATATATATATAATATATTTCAAATTCCCTTATATATCCAAATATA
AAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATA
AAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATA
AAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATATCCAAATATA
AAAGAAAGAATATATA----------ATATATTTATAATTTCCTTATATATCCAAATATA
AAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATA
AAAGAAAGAATATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATA""",
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
AAAGAAAGAATATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATA
AAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATA
AAAGAAAGAATATATA----------ATATATTTATAATTTCCTTATATATCCAAATATA
AAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATATCCAAATATA
AAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATA
AAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATA
AAAGAAAGAATATATATATATATATAATATATTTCAAATTCCCTTATATATCCAAATATA""",
        )
        self.assertEqual(
            tuple(sequence.id for sequence in alignment.sequences),
            ("seq1", "seq2", "seq3", "seq4", "seq5", "seq6", "seq7"),
        )
        alignment.sort(reverse=True)
        self.assertEqual(
            "\n".join(row for row in alignment),  # str(alignment),
            """\
AAAGAAAGAATATATATATATATATAATATATTTCAAATTCCCTTATATATCCAAATATA
AAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATA
AAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATA
AAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATATCCAAATATA
AAAGAAAGAATATATA----------ATATATTTATAATTTCCTTATATATCCAAATATA
AAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATA
AAAGAAAGAATATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATA""",
        )
        self.assertEqual(
            tuple(sequence.id for sequence in alignment.sequences),
            ("seq7", "seq6", "seq5", "seq4", "seq3", "seq2", "seq1"),
        )
        alignment.sort(key=lambda record: gc_content(record.seq))
        self.assertEqual(
            "\n".join(row for row in alignment),  # str(alignment),
            """\
AAAGAAAGAATATATA----------ATATATTTATAATTTCCTTATATATCCAAATATA
AAAGAAAGAATATATATATATATATAATATATTTCAAATTCCCTTATATATCCAAATATA
AAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATATCCAAATATA
AAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATA
AAAGAAAGAATATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATA
AAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATA
AAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATA""",
        )
        self.assertEqual(
            tuple(sequence.id for sequence in alignment.sequences),
            ("seq3", "seq7", "seq4", "seq5", "seq1", "seq6", "seq2"),
        )
        alignment.sort(key=lambda record: gc_content(record.seq), reverse=True)
        self.assertEqual(
            "\n".join(row for row in alignment),  # str(alignment),
            """\
AAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATA
AAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATA
AAAGAAAGAATATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATA
AAAGAAAGAATATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATA
AAAGAAAGAATATATA----------ATATATTTCAAATTTCCTTATATATCCAAATATA
AAAGAAAGAATATATATATATATATAATATATTTCAAATTCCCTTATATATCCAAATATA
AAAGAAAGAATATATA----------ATATATTTATAATTTCCTTATATATCCAAATATA""",
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
