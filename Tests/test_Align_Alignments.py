# Copyright 2022 by Michiel de Hoon.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for the Alignments class in Bio.Align.interfaces."""
import unittest
import warnings

from Bio import BiopythonExperimentalWarning

with warnings.catch_warnings():
    warnings.simplefilter("ignore", BiopythonExperimentalWarning)
    from Bio.Align.emboss import AlignmentIterator


from Bio.Align import Alignments, PairwiseAligner
from Bio.Align._pairwise import PairwiseAlignments


try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install numpy if you want to use Bio.Align.interfaces."
    ) from None


class TestParsedAlignments(unittest.TestCase):

    path = "Emboss/needle.txt"
    same_repr = (
        "<Alignment object (2 rows x 124 columns) at ",
        "<Alignment object (2 rows x 119 columns) at ",
        "<Alignment object (2 rows x 120 columns) at ",
        "<Alignment object (2 rows x 118 columns) at ",
        "<Alignment object (2 rows x 125 columns) at ",
    )
    other_repr = ("<Alignment object (2 rows x 131 columns) at ",)

    def setUp(self):
        self.alignments = AlignmentIterator(self.path)

    def test_repr(self):
        representation = "<Alignments object at %s>" % hex(id(self.alignments))
        self.assertIs(type(self.alignments), AlignmentIterator)
        self.assertEqual(repr(self.alignments), representation)
        alignments = self.alignments[:]
        self.assertIs(type(self.alignments), Alignments)
        self.assertEqual(repr(self.alignments), representation)
        representation = "<Alignments object at %s>" % hex(id(alignments))
        self.assertEqual(repr(alignments), representation)

    def test_lt_false(self):
        same_alignments = AlignmentIterator(self.path)
        self.assertIs(type(self.alignments), AlignmentIterator)
        self.assertIs(type(same_alignments), AlignmentIterator)
        result = self.alignments < same_alignments
        self.assertIs(type(self.alignments), Alignments)
        self.assertIs(type(same_alignments), Alignments)
        self.assertFalse(result)

    def test_lt_true(self):
        other_alignments = AlignmentIterator("Emboss/water.txt")
        self.assertIs(type(self.alignments), AlignmentIterator)
        self.assertIs(type(other_alignments), AlignmentIterator)
        result = self.alignments < other_alignments
        self.assertIs(type(self.alignments), Alignments)
        self.assertIs(type(other_alignments), Alignments)
        self.assertTrue(result)

    def test_le_false(self):
        other_alignments = AlignmentIterator("Emboss/water.txt")
        self.assertIs(type(self.alignments), AlignmentIterator)
        self.assertIs(type(other_alignments), AlignmentIterator)
        result = other_alignments <= self.alignments
        self.assertIs(type(self.alignments), Alignments)
        self.assertIs(type(other_alignments), Alignments)
        self.assertFalse(result)

    def test_le_true(self):
        same_alignments = AlignmentIterator(self.path)
        self.assertIs(type(self.alignments), AlignmentIterator)
        self.assertIs(type(same_alignments), AlignmentIterator)
        result = self.alignments <= same_alignments
        self.assertIs(type(self.alignments), Alignments)
        self.assertIs(type(same_alignments), Alignments)
        self.assertTrue(result)

    def test_eq_false(self):
        other_alignments = AlignmentIterator("Emboss/water.txt")
        self.assertIs(type(self.alignments), AlignmentIterator)
        self.assertIs(type(other_alignments), AlignmentIterator)
        result = other_alignments == self.alignments
        self.assertIs(type(self.alignments), Alignments)
        self.assertIs(type(other_alignments), Alignments)
        self.assertFalse(result)

    def test_eq_true(self):
        same_alignments = AlignmentIterator(self.path)
        self.assertIs(type(self.alignments), AlignmentIterator)
        self.assertIs(type(same_alignments), AlignmentIterator)
        result = self.alignments == same_alignments
        self.assertIs(type(self.alignments), Alignments)
        self.assertIs(type(same_alignments), Alignments)
        self.assertTrue(result)

    def test_gt_false(self):
        same_alignments = AlignmentIterator(self.path)
        self.assertIs(type(self.alignments), AlignmentIterator)
        self.assertIs(type(same_alignments), AlignmentIterator)
        result = same_alignments > self.alignments
        self.assertIs(type(self.alignments), Alignments)
        self.assertIs(type(same_alignments), Alignments)
        self.assertFalse(result)

    def test_gt_true(self):
        other_alignments = AlignmentIterator("Emboss/water.txt")
        self.assertIs(type(self.alignments), AlignmentIterator)
        self.assertIs(type(other_alignments), AlignmentIterator)
        result = other_alignments > self.alignments
        self.assertIs(type(self.alignments), Alignments)
        self.assertIs(type(other_alignments), Alignments)
        self.assertTrue(result)

    def test_ge_false(self):
        other_alignments = AlignmentIterator("Emboss/water.txt")
        self.assertIs(type(self.alignments), AlignmentIterator)
        self.assertIs(type(other_alignments), AlignmentIterator)
        result = self.alignments >= other_alignments
        self.assertIs(type(self.alignments), Alignments)
        self.assertIs(type(other_alignments), Alignments)
        self.assertFalse(result)

    def test_ge_true(self):
        same_alignments = AlignmentIterator(self.path)
        self.assertIs(type(self.alignments), AlignmentIterator)
        self.assertIs(type(same_alignments), AlignmentIterator)
        result = self.alignments >= same_alignments
        self.assertIs(type(self.alignments), Alignments)
        self.assertIs(type(same_alignments), Alignments)
        self.assertTrue(result)

    def test_contains_false(self):
        other_alignments = AlignmentIterator("Emboss/water.txt")
        alignment = next(other_alignments)
        self.assertIs(type(self.alignments), AlignmentIterator)
        result = alignment in self.alignments
        self.assertIs(type(self.alignments), Alignments)
        self.assertFalse(result)

    def test_contains_true(self):
        same_alignments = AlignmentIterator(self.path)
        alignment = next(same_alignments)
        self.assertIs(type(self.alignments), AlignmentIterator)
        result = alignment in self.alignments
        self.assertIs(type(self.alignments), Alignments)
        self.assertTrue(result)

    def test_len(self):
        self.assertIs(type(self.alignments), AlignmentIterator)
        length = len(self.alignments)
        self.assertIs(type(self.alignments), Alignments)
        self.assertEqual(length, 5)

    def test_getitem(self):
        self.assertIs(type(self.alignments), AlignmentIterator)
        alignment = self.alignments[3]
        self.assertIs(type(self.alignments), Alignments)
        same_alignments = AlignmentIterator(self.path)
        zeroth_alignment = next(same_alignments)
        self.assertNotEqual(zeroth_alignment, alignment)
        first_alignment = next(same_alignments)
        self.assertNotEqual(first_alignment, alignment)
        second_alignment = next(same_alignments)
        self.assertNotEqual(second_alignment, alignment)
        third_alignment = next(same_alignments)
        self.assertEqual(third_alignment, alignment)
        fourth_alignment = next(same_alignments)
        self.assertNotEqual(fourth_alignment, alignment)
        with self.assertRaises(StopIteration):
            next(same_alignments)
        self.assertIs(type(same_alignments), AlignmentIterator)

    def test_setitem(self):
        same_alignments = AlignmentIterator(self.path)
        zeroth_alignment = next(same_alignments)
        first_alignment = next(same_alignments)
        second_alignment = next(same_alignments)
        third_alignment = next(same_alignments)
        fourth_alignment = next(same_alignments)
        self.assertIs(type(same_alignments), AlignmentIterator)
        self.assertIs(type(self.alignments), AlignmentIterator)
        self.alignments[3] = first_alignment
        self.assertIs(type(self.alignments), Alignments)
        alignment = next(self.alignments)
        self.assertEqual(zeroth_alignment, alignment)
        alignment = next(self.alignments)
        self.assertEqual(first_alignment, alignment)
        alignment = next(self.alignments)
        self.assertEqual(second_alignment, alignment)
        alignment = next(self.alignments)
        self.assertEqual(first_alignment, alignment)
        alignment = next(self.alignments)
        self.assertEqual(fourth_alignment, alignment)
        with self.assertRaises(StopIteration):
            next(self.alignments)

    def test_delitem(self):
        same_alignments = AlignmentIterator(self.path)
        zeroth_alignment = next(same_alignments)
        first_alignment = next(same_alignments)
        second_alignment = next(same_alignments)
        third_alignment = next(same_alignments)
        fourth_alignment = next(same_alignments)
        self.assertIs(type(same_alignments), AlignmentIterator)
        self.assertIs(type(self.alignments), AlignmentIterator)
        del self.alignments[3]
        self.assertIs(type(self.alignments), Alignments)
        alignment = next(self.alignments)
        self.assertEqual(zeroth_alignment, alignment)
        alignment = next(self.alignments)
        self.assertEqual(first_alignment, alignment)
        alignment = next(self.alignments)
        self.assertEqual(second_alignment, alignment)
        alignment = next(self.alignments)
        self.assertEqual(fourth_alignment, alignment)
        with self.assertRaises(StopIteration):
            next(self.alignments)

    def test_add(self):
        same_alignments = AlignmentIterator(self.path)
        self.assertEqual(self.alignments.metadata["Align_format"], "srspair")
        self.assertEqual(self.alignments.metadata["Program"], "needle")
        self.assertEqual(
            self.alignments.metadata["Rundate"], "Sun 27 Apr 2007 17:20:35"
        )
        self.assertEqual(
            self.alignments.metadata["Command line"],
            "needle [-asequence] Spo0F.faa [-bsequence] paired_r.faa -sformat2 pearson",
        )
        self.assertEqual(self.alignments.metadata["Report_file"], "ref_rec .needle")
        self.assertEqual(len(self.alignments.metadata), 5)
        self.assertIs(type(self.alignments), AlignmentIterator)
        self.assertIs(type(same_alignments), AlignmentIterator)
        alignments = self.alignments + same_alignments
        self.assertIs(type(self.alignments), Alignments)
        self.assertIs(type(same_alignments), Alignments)
        for alignment, representation in zip(alignments, 2 * self.same_repr):
            self.assertTrue(repr(alignment).startswith(representation))
        self.assertEqual(len(alignments), 10)
        with self.assertRaises(AttributeError):
            alignments.metadata
        other_alignments = AlignmentIterator("Emboss/water.txt")
        self.assertEqual(other_alignments.metadata["Align_format"], "srspair")
        self.assertEqual(other_alignments.metadata["Program"], "water")
        self.assertEqual(
            other_alignments.metadata["Rundate"], "Wed Jan 16 17:23:19 2002"
        )
        self.assertEqual(other_alignments.metadata["Report_file"], "stdout")
        self.assertEqual(len(other_alignments.metadata), 4)

        self.assertIs(type(other_alignments), AlignmentIterator)
        alignments = self.alignments + other_alignments
        self.assertIs(type(other_alignments), Alignments)
        for alignment, representation in zip(
            alignments, self.same_repr + self.other_repr
        ):
            self.assertTrue(repr(alignment).startswith(representation))
        self.assertEqual(len(alignments), 6)
        with self.assertRaises(AttributeError):
            alignments.metadata

    def test_radd(self):
        self.assertEqual(self.alignments.metadata["Align_format"], "srspair")
        self.assertEqual(self.alignments.metadata["Program"], "needle")
        self.assertEqual(
            self.alignments.metadata["Rundate"], "Sun 27 Apr 2007 17:20:35"
        )
        self.assertEqual(
            self.alignments.metadata["Command line"],
            "needle [-asequence] Spo0F.faa [-bsequence] paired_r.faa -sformat2 pearson",
        )
        self.assertEqual(self.alignments.metadata["Report_file"], "ref_rec .needle")
        self.assertEqual(len(self.alignments.metadata), 5)
        other_alignments = AlignmentIterator("Emboss/water.txt")
        self.assertEqual(other_alignments.metadata["Align_format"], "srspair")
        self.assertEqual(other_alignments.metadata["Program"], "water")
        self.assertEqual(
            other_alignments.metadata["Rundate"], "Wed Jan 16 17:23:19 2002"
        )
        self.assertEqual(other_alignments.metadata["Report_file"], "stdout")
        self.assertEqual(len(other_alignments.metadata), 4)
        self.assertIs(type(self.alignments), AlignmentIterator)
        self.assertIs(type(other_alignments), AlignmentIterator)
        alignments = list(other_alignments) + self.alignments
        self.assertIs(type(self.alignments), Alignments)
        self.assertIs(type(other_alignments), Alignments)
        for alignment, representation in zip(
            alignments, self.other_repr + self.same_repr
        ):
            self.assertTrue(repr(alignment).startswith(representation))
        self.assertEqual(len(alignments), 6)
        with self.assertRaises(AttributeError):
            alignments.metadata

    def test_iadd(self):
        same_alignments = AlignmentIterator(self.path)
        self.assertEqual(self.alignments.metadata["Align_format"], "srspair")
        self.assertEqual(self.alignments.metadata["Program"], "needle")
        self.assertEqual(
            self.alignments.metadata["Rundate"], "Sun 27 Apr 2007 17:20:35"
        )
        self.assertEqual(
            self.alignments.metadata["Command line"],
            "needle [-asequence] Spo0F.faa [-bsequence] paired_r.faa -sformat2 pearson",
        )
        self.assertEqual(self.alignments.metadata["Report_file"], "ref_rec .needle")
        self.assertEqual(len(self.alignments.metadata), 5)
        self.assertIs(type(self.alignments), AlignmentIterator)
        self.assertIs(type(same_alignments), AlignmentIterator)
        self.alignments += same_alignments
        self.assertIs(type(self.alignments), Alignments)
        self.assertIs(type(same_alignments), Alignments)
        for alignment, representation in zip(self.alignments, 2 * self.same_repr):
            self.assertTrue(repr(alignment).startswith(representation))
        self.assertEqual(len(self.alignments), 10)
        self.assertEqual(self.alignments.metadata["Align_format"], "srspair")
        self.assertEqual(self.alignments.metadata["Program"], "needle")
        self.assertEqual(
            self.alignments.metadata["Rundate"], "Sun 27 Apr 2007 17:20:35"
        )
        self.assertEqual(
            self.alignments.metadata["Command line"],
            "needle [-asequence] Spo0F.faa [-bsequence] paired_r.faa -sformat2 pearson",
        )
        self.assertEqual(self.alignments.metadata["Report_file"], "ref_rec .needle")
        other_alignments = AlignmentIterator("Emboss/water.txt")
        self.assertEqual(other_alignments.metadata["Align_format"], "srspair")
        self.assertEqual(other_alignments.metadata["Program"], "water")
        self.assertEqual(
            other_alignments.metadata["Rundate"], "Wed Jan 16 17:23:19 2002"
        )
        self.assertEqual(other_alignments.metadata["Report_file"], "stdout")
        self.assertEqual(len(other_alignments.metadata), 4)

        self.assertIs(type(other_alignments), AlignmentIterator)
        self.alignments += other_alignments
        self.assertIs(type(other_alignments), Alignments)
        for alignment, representation in zip(
            self.alignments, 2 * self.same_repr + self.other_repr
        ):
            self.assertTrue(repr(alignment).startswith(representation))
        self.assertEqual(len(self.alignments), 11)
        self.assertEqual(self.alignments.metadata["Align_format"], "srspair")
        self.assertEqual(self.alignments.metadata["Program"], "needle")
        self.assertEqual(
            self.alignments.metadata["Rundate"], "Sun 27 Apr 2007 17:20:35"
        )
        self.assertEqual(
            self.alignments.metadata["Command line"],
            "needle [-asequence] Spo0F.faa [-bsequence] paired_r.faa -sformat2 pearson",
        )
        self.assertEqual(self.alignments.metadata["Report_file"], "ref_rec .needle")

    def test_mul(self):
        self.assertEqual(self.alignments.metadata["Align_format"], "srspair")
        self.assertEqual(self.alignments.metadata["Program"], "needle")
        self.assertEqual(
            self.alignments.metadata["Rundate"], "Sun 27 Apr 2007 17:20:35"
        )
        self.assertEqual(
            self.alignments.metadata["Command line"],
            "needle [-asequence] Spo0F.faa [-bsequence] paired_r.faa -sformat2 pearson",
        )
        self.assertEqual(self.alignments.metadata["Report_file"], "ref_rec .needle")
        self.assertEqual(len(self.alignments.metadata), 5)
        self.assertIs(type(self.alignments), AlignmentIterator)
        alignments = 3 * self.alignments
        self.assertIs(type(self.alignments), Alignments)
        for alignment, representation in zip(alignments, 3 * self.same_repr):
            self.assertTrue(repr(alignment).startswith(representation))
        self.assertEqual(len(alignments), 15)
        with self.assertRaises(AttributeError):
            alignments.metadata

    def test_rmul(self):
        self.assertEqual(self.alignments.metadata["Align_format"], "srspair")
        self.assertEqual(self.alignments.metadata["Program"], "needle")
        self.assertEqual(
            self.alignments.metadata["Rundate"], "Sun 27 Apr 2007 17:20:35"
        )
        self.assertEqual(
            self.alignments.metadata["Command line"],
            "needle [-asequence] Spo0F.faa [-bsequence] paired_r.faa -sformat2 pearson",
        )
        self.assertEqual(self.alignments.metadata["Report_file"], "ref_rec .needle")
        self.assertEqual(len(self.alignments.metadata), 5)
        self.assertIs(type(self.alignments), AlignmentIterator)
        alignments = self.alignments * 3
        self.assertIs(type(self.alignments), Alignments)
        for alignment, representation in zip(alignments, 3 * self.same_repr):
            self.assertTrue(repr(alignment).startswith(representation))
        self.assertEqual(len(alignments), 15)
        with self.assertRaises(AttributeError):
            alignments.metadata

    def test_imul(self):
        self.assertEqual(self.alignments.metadata["Align_format"], "srspair")
        self.assertEqual(self.alignments.metadata["Program"], "needle")
        self.assertEqual(
            self.alignments.metadata["Rundate"], "Sun 27 Apr 2007 17:20:35"
        )
        self.assertEqual(
            self.alignments.metadata["Command line"],
            "needle [-asequence] Spo0F.faa [-bsequence] paired_r.faa -sformat2 pearson",
        )
        self.assertEqual(self.alignments.metadata["Report_file"], "ref_rec .needle")
        self.assertEqual(len(self.alignments.metadata), 5)
        self.assertIs(type(self.alignments), AlignmentIterator)
        self.alignments *= 3
        self.assertIs(type(self.alignments), Alignments)
        for alignment, representation in zip(self.alignments, 3 * self.same_repr):
            self.assertTrue(repr(alignment).startswith(representation))
        self.assertEqual(len(self.alignments), 15)
        self.assertEqual(self.alignments.metadata["Align_format"], "srspair")
        self.assertEqual(self.alignments.metadata["Program"], "needle")
        self.assertEqual(
            self.alignments.metadata["Rundate"], "Sun 27 Apr 2007 17:20:35"
        )
        self.assertEqual(
            self.alignments.metadata["Command line"],
            "needle [-asequence] Spo0F.faa [-bsequence] paired_r.faa -sformat2 pearson",
        )
        self.assertEqual(self.alignments.metadata["Report_file"], "ref_rec .needle")
        self.assertEqual(len(self.alignments.metadata), 5)

    def test_append(self):
        same_alignments = AlignmentIterator(self.path)
        zeroth_alignment = next(same_alignments)
        first_alignment = next(same_alignments)
        second_alignment = next(same_alignments)
        third_alignment = next(same_alignments)
        fourth_alignment = next(same_alignments)
        self.assertIs(type(same_alignments), AlignmentIterator)
        self.assertIs(type(self.alignments), AlignmentIterator)
        self.alignments.append(second_alignment)
        self.assertIs(type(self.alignments), Alignments)
        alignment = next(self.alignments)
        self.assertEqual(zeroth_alignment, alignment)
        alignment = next(self.alignments)
        self.assertEqual(first_alignment, alignment)
        alignment = next(self.alignments)
        self.assertEqual(second_alignment, alignment)
        alignment = next(self.alignments)
        self.assertEqual(third_alignment, alignment)
        alignment = next(self.alignments)
        self.assertEqual(fourth_alignment, alignment)
        alignment = next(self.alignments)
        self.assertEqual(second_alignment, alignment)
        with self.assertRaises(StopIteration):
            next(self.alignments)

    def test_insert(self):
        same_alignments = AlignmentIterator(self.path)
        zeroth_alignment = next(same_alignments)
        first_alignment = next(same_alignments)
        second_alignment = next(same_alignments)
        third_alignment = next(same_alignments)
        fourth_alignment = next(same_alignments)
        self.assertIs(type(same_alignments), AlignmentIterator)
        self.assertIs(type(self.alignments), AlignmentIterator)
        self.alignments.insert(1, fourth_alignment)
        self.assertIs(type(self.alignments), Alignments)
        alignment = next(self.alignments)
        self.assertEqual(zeroth_alignment, alignment)
        alignment = next(self.alignments)
        self.assertEqual(fourth_alignment, alignment)
        alignment = next(self.alignments)
        self.assertEqual(first_alignment, alignment)
        alignment = next(self.alignments)
        self.assertEqual(second_alignment, alignment)
        alignment = next(self.alignments)
        self.assertEqual(third_alignment, alignment)
        alignment = next(self.alignments)
        self.assertEqual(fourth_alignment, alignment)
        with self.assertRaises(StopIteration):
            next(self.alignments)

    def test_pop(self):
        self.assertIs(type(self.alignments), AlignmentIterator)
        alignment = self.alignments.pop(1)
        self.assertIs(type(self.alignments), Alignments)
        same_alignments = AlignmentIterator(self.path)
        zeroth_alignment = next(same_alignments)
        self.assertEqual(zeroth_alignment, next(self.alignments))
        first_alignment = next(same_alignments)
        self.assertEqual(first_alignment, alignment)
        second_alignment = next(same_alignments)
        self.assertEqual(second_alignment, next(self.alignments))
        third_alignment = next(same_alignments)
        self.assertEqual(third_alignment, next(self.alignments))
        fourth_alignment = next(same_alignments)
        self.assertEqual(fourth_alignment, next(self.alignments))
        with self.assertRaises(StopIteration):
            next(self.alignments)
        with self.assertRaises(StopIteration):
            next(same_alignments)
        self.assertIs(type(same_alignments), AlignmentIterator)

    def test_remove(self):
        same_alignments = AlignmentIterator(self.path)
        zeroth_alignment = next(same_alignments)
        first_alignment = next(same_alignments)
        second_alignment = next(same_alignments)
        third_alignment = next(same_alignments)
        fourth_alignment = next(same_alignments)
        with self.assertRaises(StopIteration):
            next(same_alignments)
        self.assertIs(type(same_alignments), AlignmentIterator)
        self.assertIs(type(self.alignments), AlignmentIterator)
        result = self.alignments.remove(first_alignment)
        self.assertIs(type(self.alignments), Alignments)
        self.assertIsNone(result)
        self.assertEqual(zeroth_alignment, next(self.alignments))
        self.assertEqual(second_alignment, next(self.alignments))
        self.assertEqual(third_alignment, next(self.alignments))
        self.assertEqual(fourth_alignment, next(self.alignments))
        with self.assertRaises(StopIteration):
            next(self.alignments)

    def test_clear(self):
        self.assertIs(type(self.alignments), AlignmentIterator)
        result = self.alignments.clear()
        self.assertIsNone(result)
        self.assertIs(type(self.alignments), AlignmentIterator)
        with self.assertRaises(AttributeError):
            self.alignments._stream
        with self.assertRaises(StopIteration):
            next(self.alignments)
        self.assertEqual(len(self.alignments), 0)

    def test_copy(self):
        self.assertIs(type(self.alignments), AlignmentIterator)
        copied_alignments = self.alignments.copy()
        self.assertIs(type(self.alignments), Alignments)
        self.assertIs(type(copied_alignments), Alignments)
        self.assertEqual(len(self.alignments), 5)
        for alignment, copied_alignment in zip(self.alignments, copied_alignments):
            self.assertEqual(id(alignment), id(copied_alignment))
        with self.assertRaises(StopIteration):
            next(self.alignments)
        with self.assertRaises(StopIteration):
            next(copied_alignments)
        self.assertEqual(self.alignments.metadata, copied_alignments.metadata)

    def test_count(self):
        same_alignments = AlignmentIterator(self.path)
        zeroth_alignment = next(same_alignments)
        first_alignment = next(same_alignments)
        second_alignment = next(same_alignments)
        third_alignment = next(same_alignments)
        fourth_alignment = next(same_alignments)
        with self.assertRaises(StopIteration):
            next(same_alignments)
        self.assertIs(type(same_alignments), AlignmentIterator)
        self.assertIs(type(self.alignments), AlignmentIterator)
        count = self.alignments.count(first_alignment)
        self.assertIs(type(self.alignments), Alignments)
        self.assertEqual(count, 1)
        count = self.alignments.count(second_alignment)
        self.assertEqual(count, 1)
        self.alignments.insert(3, first_alignment)
        count = self.alignments.count(first_alignment)
        self.assertEqual(count, 2)
        same_alignments.insert(3, first_alignment)
        self.assertIs(type(same_alignments), Alignments)
        count = self.alignments.count(first_alignment)
        self.assertEqual(count, 2)

    def test_index(self):
        same_alignments = AlignmentIterator(self.path)
        zeroth_alignment = next(same_alignments)
        first_alignment = next(same_alignments)
        second_alignment = next(same_alignments)
        third_alignment = next(same_alignments)
        fourth_alignment = next(same_alignments)
        with self.assertRaises(StopIteration):
            next(same_alignments)
        self.assertIs(type(same_alignments), AlignmentIterator)
        self.assertIs(type(self.alignments), AlignmentIterator)
        index = self.alignments.index(third_alignment)
        self.assertIs(type(self.alignments), Alignments)
        self.assertEqual(index, 3)

    def test_reverse(self):
        same_alignments = AlignmentIterator(self.path)
        zeroth_alignment = next(same_alignments)
        first_alignment = next(same_alignments)
        second_alignment = next(same_alignments)
        third_alignment = next(same_alignments)
        fourth_alignment = next(same_alignments)
        with self.assertRaises(StopIteration):
            next(same_alignments)
        self.assertIs(type(same_alignments), AlignmentIterator)
        self.assertIs(type(self.alignments), AlignmentIterator)
        self.alignments.reverse()
        self.assertIs(type(self.alignments), Alignments)
        self.assertEqual(next(self.alignments), fourth_alignment)
        self.assertEqual(next(self.alignments), third_alignment)
        self.assertEqual(next(self.alignments), second_alignment)
        self.assertEqual(next(self.alignments), first_alignment)
        self.assertEqual(next(self.alignments), zeroth_alignment)
        with self.assertRaises(StopIteration):
            next(self.alignments)
        alignments = AlignmentIterator(self.path)
        self.assertIs(type(alignments), AlignmentIterator)
        alignments.reverse()
        self.assertIs(type(alignments), Alignments)
        self.assertEqual(alignments[0], fourth_alignment)
        self.assertEqual(alignments[1], third_alignment)
        self.assertEqual(alignments[2], second_alignment)
        self.assertEqual(alignments[3], first_alignment)
        self.assertEqual(alignments[4], zeroth_alignment)
        self.assertEqual(len(alignments), 5)
        self.assertEqual(next(alignments), fourth_alignment)
        self.assertEqual(next(alignments), third_alignment)
        self.assertEqual(next(alignments), second_alignment)
        self.assertEqual(next(alignments), first_alignment)
        self.assertEqual(next(alignments), zeroth_alignment)
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_sort(self):
        same_alignments = AlignmentIterator(self.path)
        zeroth_alignment = next(same_alignments)
        first_alignment = next(same_alignments)
        second_alignment = next(same_alignments)
        third_alignment = next(same_alignments)
        fourth_alignment = next(same_alignments)
        with self.assertRaises(StopIteration):
            next(same_alignments)
        self.assertIs(type(self.alignments), AlignmentIterator)
        self.alignments.sort(key=lambda alignment: alignment.score, reverse=True)
        self.assertIs(type(self.alignments), Alignments)
        alignment = next(self.alignments)
        self.assertEqual(alignment, fourth_alignment)
        self.assertAlmostEqual(alignment.score, 156.5)
        alignment = next(self.alignments)
        self.assertEqual(alignment, first_alignment)
        self.assertAlmostEqual(alignment.score, 154.0)
        alignment = next(self.alignments)
        self.assertAlmostEqual(alignment.score, 126.0)
        self.assertEqual(alignment, third_alignment)
        alignment = next(self.alignments)
        self.assertAlmostEqual(alignment.score, 121.0)
        self.assertEqual(alignment, second_alignment)
        alignment = next(self.alignments)
        self.assertAlmostEqual(alignment.score, 112.0)
        self.assertEqual(alignment, zeroth_alignment)
        with self.assertRaises(StopIteration):
            next(self.alignments)

    def test_extend(self):
        same_alignments = AlignmentIterator(self.path)
        self.assertEqual(self.alignments.metadata["Align_format"], "srspair")
        self.assertEqual(self.alignments.metadata["Program"], "needle")
        self.assertEqual(
            self.alignments.metadata["Rundate"], "Sun 27 Apr 2007 17:20:35"
        )
        self.assertEqual(
            self.alignments.metadata["Command line"],
            "needle [-asequence] Spo0F.faa [-bsequence] paired_r.faa -sformat2 pearson",
        )
        self.assertEqual(self.alignments.metadata["Report_file"], "ref_rec .needle")
        self.assertEqual(len(self.alignments.metadata), 5)
        self.assertIs(type(self.alignments), AlignmentIterator)
        self.assertIs(type(same_alignments), AlignmentIterator)
        self.alignments.extend(same_alignments)
        self.assertIs(type(self.alignments), Alignments)
        self.assertIs(type(same_alignments), Alignments)
        for alignment, representation in zip(self.alignments, 2 * self.same_repr):
            self.assertTrue(repr(alignment).startswith(representation))
        self.assertEqual(len(self.alignments), 10)
        self.assertEqual(self.alignments.metadata["Align_format"], "srspair")
        self.assertEqual(self.alignments.metadata["Program"], "needle")
        self.assertEqual(
            self.alignments.metadata["Rundate"], "Sun 27 Apr 2007 17:20:35"
        )
        self.assertEqual(
            self.alignments.metadata["Command line"],
            "needle [-asequence] Spo0F.faa [-bsequence] paired_r.faa -sformat2 pearson",
        )
        self.assertEqual(self.alignments.metadata["Report_file"], "ref_rec .needle")
        other_alignments = AlignmentIterator("Emboss/water.txt")
        self.assertEqual(other_alignments.metadata["Align_format"], "srspair")
        self.assertEqual(other_alignments.metadata["Program"], "water")
        self.assertEqual(
            other_alignments.metadata["Rundate"], "Wed Jan 16 17:23:19 2002"
        )
        self.assertEqual(other_alignments.metadata["Report_file"], "stdout")
        self.assertEqual(len(other_alignments.metadata), 4)

        self.assertIs(type(other_alignments), AlignmentIterator)
        self.alignments.extend(other_alignments)
        self.assertIs(type(other_alignments), Alignments)
        for alignment, representation in zip(
            self.alignments, 2 * self.same_repr + self.other_repr
        ):
            self.assertTrue(repr(alignment).startswith(representation))
        self.assertEqual(len(self.alignments), 11)
        self.assertEqual(self.alignments.metadata["Align_format"], "srspair")
        self.assertEqual(self.alignments.metadata["Program"], "needle")
        self.assertEqual(
            self.alignments.metadata["Rundate"], "Sun 27 Apr 2007 17:20:35"
        )
        self.assertEqual(
            self.alignments.metadata["Command line"],
            "needle [-asequence] Spo0F.faa [-bsequence] paired_r.faa -sformat2 pearson",
        )
        self.assertEqual(self.alignments.metadata["Report_file"], "ref_rec .needle")


class TestPairwiseAlignments(unittest.TestCase):
    target = "AACCGG"
    query = "ACCG"
    other_query = "CCGG"
    same_str = [
        """\
AACCGG
|-|||-
A-CCG-
""",
        """\
AACCGG
-||||-
-ACCG-
""",
        """\
AACCGG
|-||-|
A-CC-G
""",
        """\
AACCGG
-|||-|
-ACC-G
""",
    ]
    other_str = [
        """\
AACCGG
--||||
--CCGG
"""
    ]

    def setUp(self):
        self.aligner = PairwiseAligner()
        self.alignments = self.aligner.align(self.target, self.query)

    def test_repr(self):
        representation = "<Alignments object at %s>" % hex(id(self.alignments))
        self.assertIs(type(self.alignments), PairwiseAlignments)
        self.assertEqual(repr(self.alignments), representation)
        alignments = self.alignments[:]
        self.assertIs(type(self.alignments), Alignments)
        self.assertEqual(repr(self.alignments), representation)
        representation = "<Alignments object at %s>" % hex(id(alignments))
        self.assertEqual(repr(alignments), representation)

    def test_lt_false(self):
        same_alignments = self.aligner.align(self.target, self.query)
        self.assertIs(type(self.alignments), PairwiseAlignments)
        self.assertIs(type(same_alignments), PairwiseAlignments)
        result = self.alignments < same_alignments
        self.assertIs(type(self.alignments), Alignments)
        self.assertIs(type(same_alignments), Alignments)
        self.assertFalse(result)

    def test_lt_true(self):
        other_alignments = self.aligner.align(self.target, self.other_query)
        self.assertIs(type(self.alignments), PairwiseAlignments)
        self.assertIs(type(other_alignments), PairwiseAlignments)
        result = self.alignments < other_alignments
        self.assertIs(type(self.alignments), Alignments)
        self.assertIs(type(other_alignments), Alignments)
        self.assertTrue(result)

    def test_le_false(self):
        other_alignments = self.aligner.align(self.target, self.other_query)
        self.assertIs(type(self.alignments), PairwiseAlignments)
        self.assertIs(type(other_alignments), PairwiseAlignments)
        result = other_alignments <= self.alignments
        self.assertIs(type(self.alignments), Alignments)
        self.assertIs(type(other_alignments), Alignments)
        self.assertFalse(result)

    def test_le_true(self):
        same_alignments = self.aligner.align(self.target, self.query)
        self.assertIs(type(self.alignments), PairwiseAlignments)
        self.assertIs(type(same_alignments), PairwiseAlignments)
        result = self.alignments <= same_alignments
        self.assertIs(type(self.alignments), Alignments)
        self.assertIs(type(same_alignments), Alignments)
        self.assertTrue(result)

    def test_eq_false(self):
        other_alignments = self.aligner.align(self.target, self.other_query)
        self.assertIs(type(self.alignments), PairwiseAlignments)
        self.assertIs(type(other_alignments), PairwiseAlignments)
        result = other_alignments == self.alignments
        self.assertIs(type(self.alignments), Alignments)
        self.assertIs(type(other_alignments), Alignments)
        self.assertFalse(result)

    def test_eq_true(self):
        same_alignments = self.aligner.align(self.target, self.query)
        self.assertIs(type(self.alignments), PairwiseAlignments)
        self.assertIs(type(same_alignments), PairwiseAlignments)
        result = self.alignments == same_alignments
        self.assertIs(type(self.alignments), Alignments)
        self.assertIs(type(same_alignments), Alignments)
        self.assertTrue(result)

    def test_gt_false(self):
        same_alignments = self.aligner.align(self.target, self.query)
        self.assertIs(type(self.alignments), PairwiseAlignments)
        self.assertIs(type(same_alignments), PairwiseAlignments)
        result = same_alignments > self.alignments
        self.assertIs(type(self.alignments), Alignments)
        self.assertIs(type(same_alignments), Alignments)
        self.assertFalse(result)

    def test_gt_true(self):
        other_alignments = self.aligner.align(self.target, self.other_query)
        self.assertIs(type(self.alignments), PairwiseAlignments)
        self.assertIs(type(other_alignments), PairwiseAlignments)
        result = other_alignments > self.alignments
        self.assertIs(type(self.alignments), Alignments)
        self.assertIs(type(other_alignments), Alignments)
        self.assertTrue(result)

    def test_ge_false(self):
        other_alignments = self.aligner.align(self.target, self.other_query)
        self.assertIs(type(self.alignments), PairwiseAlignments)
        self.assertIs(type(other_alignments), PairwiseAlignments)
        result = self.alignments >= other_alignments
        self.assertIs(type(self.alignments), Alignments)
        self.assertIs(type(other_alignments), Alignments)
        self.assertFalse(result)

    def test_ge_true(self):
        same_alignments = self.aligner.align(self.target, self.query)
        self.assertIs(type(self.alignments), PairwiseAlignments)
        self.assertIs(type(same_alignments), PairwiseAlignments)
        result = self.alignments >= same_alignments
        self.assertIs(type(self.alignments), Alignments)
        self.assertIs(type(same_alignments), Alignments)
        self.assertTrue(result)

    def test_contains_false(self):
        other_alignments = self.aligner.align(self.target, self.other_query)
        alignment = next(other_alignments)
        self.assertIs(type(self.alignments), PairwiseAlignments)
        result = alignment in self.alignments
        self.assertIs(type(self.alignments), Alignments)
        self.assertFalse(result)

    def test_contains_true(self):
        same_alignments = self.aligner.align(self.target, self.query)
        alignment = next(same_alignments)
        self.assertIs(type(self.alignments), PairwiseAlignments)
        result = alignment in self.alignments
        self.assertIs(type(self.alignments), Alignments)
        self.assertTrue(result)

    def test_len(self):
        self.assertIs(type(self.alignments), PairwiseAlignments)
        length = len(self.alignments)
        self.assertIs(type(self.alignments), PairwiseAlignments)
        self.assertEqual(length, 4)

    def test_getitem(self):
        self.assertIs(type(self.alignments), PairwiseAlignments)
        alignment = self.alignments[2]
        self.assertIs(type(self.alignments), PairwiseAlignments)
        same_alignments = self.aligner.align(self.target, self.query)
        zeroth_alignment = next(same_alignments)
        self.assertNotEqual(zeroth_alignment, alignment)
        first_alignment = next(same_alignments)
        self.assertNotEqual(first_alignment, alignment)
        second_alignment = next(same_alignments)
        self.assertEqual(second_alignment, alignment)
        third_alignment = next(same_alignments)
        self.assertNotEqual(third_alignment, alignment)
        with self.assertRaises(StopIteration):
            next(same_alignments)
        self.assertIs(type(same_alignments), PairwiseAlignments)

    def test_setitem(self):
        same_alignments = self.aligner.align(self.target, self.query)
        zeroth_alignment = next(same_alignments)
        first_alignment = next(same_alignments)
        second_alignment = next(same_alignments)
        third_alignment = next(same_alignments)
        self.assertIs(type(same_alignments), PairwiseAlignments)
        self.assertIs(type(self.alignments), PairwiseAlignments)
        self.alignments[2] = first_alignment
        self.assertIs(type(self.alignments), Alignments)
        alignment = next(self.alignments)
        self.assertEqual(zeroth_alignment, alignment)
        alignment = next(self.alignments)
        self.assertEqual(first_alignment, alignment)
        alignment = next(self.alignments)
        self.assertEqual(first_alignment, alignment)
        alignment = next(self.alignments)
        self.assertEqual(third_alignment, alignment)
        with self.assertRaises(StopIteration):
            next(self.alignments)

    def test_delitem(self):
        same_alignments = self.aligner.align(self.target, self.query)
        zeroth_alignment = next(same_alignments)
        first_alignment = next(same_alignments)
        second_alignment = next(same_alignments)
        third_alignment = next(same_alignments)
        self.assertIs(type(same_alignments), PairwiseAlignments)
        self.assertIs(type(self.alignments), PairwiseAlignments)
        del self.alignments[2]
        self.assertIs(type(self.alignments), Alignments)
        alignment = next(self.alignments)
        self.assertEqual(zeroth_alignment, alignment)
        alignment = next(self.alignments)
        self.assertEqual(first_alignment, alignment)
        alignment = next(self.alignments)
        self.assertEqual(third_alignment, alignment)
        with self.assertRaises(StopIteration):
            next(self.alignments)

    def test_add(self):
        same_alignments = self.aligner.align(self.target, self.query)
        self.assertAlmostEqual(self.alignments.score, 4.0)
        self.assertIs(type(self.alignments), PairwiseAlignments)
        self.assertIs(type(same_alignments), PairwiseAlignments)
        alignments = self.alignments + same_alignments
        self.assertIs(type(self.alignments), PairwiseAlignments)
        self.assertIs(type(same_alignments), PairwiseAlignments)
        for alignment, text in zip(alignments, 2 * self.same_str):
            self.assertEqual(str(alignment), text)
        self.assertEqual(len(alignments), 8)
        with self.assertRaises(AttributeError):
            alignments.score
        other_alignments = self.aligner.align(self.target, self.other_query)
        self.assertAlmostEqual(self.alignments.score, 4.0)
        self.assertIs(type(other_alignments), PairwiseAlignments)
        alignments = self.alignments + other_alignments
        self.assertIs(type(other_alignments), PairwiseAlignments)
        for alignment, text in zip(alignments, self.same_str + self.other_str):
            self.assertEqual(str(alignment), text)
        self.assertEqual(len(alignments), 5)
        with self.assertRaises(AttributeError):
            alignments.metadata

    def test_radd(self):
        self.assertAlmostEqual(self.alignments.score, 4.0)
        other_alignments = self.aligner.align(self.target, self.other_query)
        self.assertAlmostEqual(other_alignments.score, 4.0)
        self.assertIs(type(self.alignments), PairwiseAlignments)
        self.assertIs(type(other_alignments), PairwiseAlignments)
        alignments = list(other_alignments) + self.alignments
        self.assertIs(type(self.alignments), PairwiseAlignments)
        self.assertIs(type(other_alignments), PairwiseAlignments)
        self.assertIs(type(alignments), Alignments)
        for alignment, text in zip(alignments, self.other_str + self.same_str):
            self.assertEqual(str(alignment), text)
        self.assertEqual(len(alignments), 5)
        with self.assertRaises(AttributeError):
            alignments.score

    def test_iadd(self):
        same_alignments = self.aligner.align(self.target, self.query)
        self.assertAlmostEqual(self.alignments.score, 4.0)
        self.assertIs(type(self.alignments), PairwiseAlignments)
        self.assertIs(type(same_alignments), PairwiseAlignments)
        self.assertEqual(len(self.alignments), 4)
        self.alignments += same_alignments
        self.assertIs(type(self.alignments), Alignments)
        self.assertIs(type(same_alignments), PairwiseAlignments)
        for alignment, text in zip(self.alignments, 2 * self.same_str):
            self.assertEqual(str(alignment), text)
        self.assertEqual(len(self.alignments), 8)
        self.assertAlmostEqual(self.alignments.score, 4.0)
        other_alignments = self.aligner.align(self.target, self.other_query)
        self.assertAlmostEqual(other_alignments.score, 4.0)
        self.assertIs(type(other_alignments), PairwiseAlignments)
        self.alignments += other_alignments
        self.assertIs(type(other_alignments), PairwiseAlignments)
        for alignment, text in zip(self.alignments, 2 * self.same_str + self.other_str):
            self.assertEqual(str(alignment), text)
        self.assertEqual(len(self.alignments), 9)
        self.assertAlmostEqual(self.alignments.score, 4.0)

    def test_mul(self):
        self.assertAlmostEqual(self.alignments.score, 4.0)
        self.assertIs(type(self.alignments), PairwiseAlignments)
        alignments = 3 * self.alignments
        self.assertIs(type(self.alignments), PairwiseAlignments)
        self.assertIs(type(alignments), Alignments)
        for alignment, text in zip(alignments, 3 * self.same_str):
            self.assertEqual(str(alignment), text)
        self.assertEqual(len(alignments), 12)
        with self.assertRaises(AttributeError):
            alignments.score

    def test_rmul(self):
        self.assertAlmostEqual(self.alignments.score, 4.0)
        self.assertIs(type(self.alignments), PairwiseAlignments)
        alignments = self.alignments * 3
        self.assertIs(type(self.alignments), PairwiseAlignments)
        self.assertIs(type(alignments), Alignments)
        for alignment, text in zip(alignments, 3 * self.same_str):
            self.assertEqual(str(alignment), text)
        self.assertEqual(len(alignments), 12)
        with self.assertRaises(AttributeError):
            alignments.score

    def test_imul(self):
        self.assertAlmostEqual(self.alignments.score, 4.0)
        self.assertIs(type(self.alignments), PairwiseAlignments)
        self.alignments *= 3
        self.assertIs(type(self.alignments), Alignments)
        for alignment, text in zip(self.alignments, 3 * self.same_str):
            self.assertEqual(str(alignment), text)
        self.assertEqual(len(self.alignments), 12)
        self.assertAlmostEqual(self.alignments.score, 4.0)

    def test_append(self):
        same_alignments = self.aligner.align(self.target, self.query)
        zeroth_alignment = next(same_alignments)
        first_alignment = next(same_alignments)
        second_alignment = next(same_alignments)
        third_alignment = next(same_alignments)
        with self.assertRaises(StopIteration):
            next(same_alignments)
        self.assertIs(type(same_alignments), PairwiseAlignments)
        self.assertIs(type(self.alignments), PairwiseAlignments)
        self.alignments.append(second_alignment)
        self.assertIs(type(self.alignments), Alignments)
        alignment = next(self.alignments)
        self.assertEqual(zeroth_alignment, alignment)
        alignment = next(self.alignments)
        self.assertEqual(first_alignment, alignment)
        alignment = next(self.alignments)
        self.assertEqual(second_alignment, alignment)
        alignment = next(self.alignments)
        self.assertEqual(third_alignment, alignment)
        alignment = next(self.alignments)
        self.assertEqual(second_alignment, alignment)
        with self.assertRaises(StopIteration):
            next(self.alignments)

    def test_insert(self):
        same_alignments = self.aligner.align(self.target, self.query)
        zeroth_alignment = next(same_alignments)
        first_alignment = next(same_alignments)
        second_alignment = next(same_alignments)
        third_alignment = next(same_alignments)
        with self.assertRaises(StopIteration):
            next(same_alignments)
        self.assertIs(type(same_alignments), PairwiseAlignments)
        self.assertIs(type(self.alignments), PairwiseAlignments)
        self.alignments.insert(1, third_alignment)
        self.assertIs(type(self.alignments), Alignments)
        alignment = next(self.alignments)
        self.assertEqual(zeroth_alignment, alignment)
        alignment = next(self.alignments)
        self.assertEqual(third_alignment, alignment)
        alignment = next(self.alignments)
        self.assertEqual(first_alignment, alignment)
        alignment = next(self.alignments)
        self.assertEqual(second_alignment, alignment)
        alignment = next(self.alignments)
        self.assertEqual(third_alignment, alignment)
        with self.assertRaises(StopIteration):
            next(self.alignments)

    def test_pop(self):
        self.assertIs(type(self.alignments), PairwiseAlignments)
        alignment = self.alignments.pop(1)
        self.assertIs(type(self.alignments), Alignments)
        same_alignments = self.aligner.align(self.target, self.query)
        zeroth_alignment = next(same_alignments)
        self.assertEqual(zeroth_alignment, next(self.alignments))
        first_alignment = next(same_alignments)
        self.assertEqual(first_alignment, alignment)
        second_alignment = next(same_alignments)
        self.assertEqual(second_alignment, next(self.alignments))
        third_alignment = next(same_alignments)
        self.assertEqual(third_alignment, next(self.alignments))
        with self.assertRaises(StopIteration):
            next(self.alignments)
        with self.assertRaises(StopIteration):
            next(same_alignments)
        self.assertIs(type(same_alignments), PairwiseAlignments)

    def test_remove(self):
        same_alignments = self.aligner.align(self.target, self.query)
        zeroth_alignment = next(same_alignments)
        first_alignment = next(same_alignments)
        second_alignment = next(same_alignments)
        third_alignment = next(same_alignments)
        with self.assertRaises(StopIteration):
            next(same_alignments)
        self.assertIs(type(self.alignments), PairwiseAlignments)
        self.assertIs(type(same_alignments), PairwiseAlignments)
        result = self.alignments.remove(first_alignment)
        self.assertIs(type(self.alignments), Alignments)
        self.assertIsNone(result)
        self.assertEqual(zeroth_alignment, next(self.alignments))
        self.assertEqual(second_alignment, next(self.alignments))
        self.assertEqual(third_alignment, next(self.alignments))
        with self.assertRaises(StopIteration):
            next(self.alignments)

    def test_clear(self):
        self.assertIs(type(self.alignments), PairwiseAlignments)
        result = self.alignments.clear()
        self.assertIsNone(result)
        self.assertIs(type(self.alignments), Alignments)
        with self.assertRaises(AttributeError):
            self.alignments.paths
        with self.assertRaises(StopIteration):
            next(self.alignments)
        self.assertEqual(len(self.alignments), 0)

    def test_copy(self):
        self.assertIs(type(self.alignments), PairwiseAlignments)
        copied_alignments = self.alignments.copy()
        self.assertIs(type(self.alignments), Alignments)
        self.assertIs(type(copied_alignments), Alignments)
        self.assertEqual(len(self.alignments), 4)
        for alignment, copied_alignment in zip(self.alignments, copied_alignments):
            self.assertEqual(id(alignment), id(copied_alignment))
        with self.assertRaises(StopIteration):
            next(self.alignments)
        with self.assertRaises(StopIteration):
            next(copied_alignments)
        self.assertAlmostEqual(self.alignments.score, copied_alignments.score)

    def test_count(self):
        same_alignments = self.aligner.align(self.target, self.query)
        zeroth_alignment = next(same_alignments)
        first_alignment = next(same_alignments)
        second_alignment = next(same_alignments)
        third_alignment = next(same_alignments)
        with self.assertRaises(StopIteration):
            next(same_alignments)
        self.assertIs(type(same_alignments), PairwiseAlignments)
        self.assertIs(type(self.alignments), PairwiseAlignments)
        count = self.alignments.count(first_alignment)
        self.assertIs(type(self.alignments), Alignments)
        self.assertEqual(count, 1)
        count = self.alignments.count(second_alignment)
        self.assertEqual(count, 1)
        self.alignments.insert(3, first_alignment)
        count = self.alignments.count(first_alignment)
        self.assertEqual(count, 2)
        same_alignments.insert(3, first_alignment)
        self.assertIs(type(same_alignments), Alignments)
        count = self.alignments.count(first_alignment)
        self.assertEqual(count, 2)

    def test_index(self):
        same_alignments = self.aligner.align(self.target, self.query)
        zeroth_alignment = next(same_alignments)
        first_alignment = next(same_alignments)
        second_alignment = next(same_alignments)
        third_alignment = next(same_alignments)
        with self.assertRaises(StopIteration):
            next(same_alignments)
        self.assertIs(type(same_alignments), PairwiseAlignments)
        self.assertIs(type(self.alignments), PairwiseAlignments)
        index = self.alignments.index(second_alignment)
        self.assertIs(type(self.alignments), Alignments)
        self.assertEqual(index, 2)

    def test_reverse(self):
        same_alignments = self.aligner.align(self.target, self.query)
        zeroth_alignment = next(same_alignments)
        first_alignment = next(same_alignments)
        second_alignment = next(same_alignments)
        third_alignment = next(same_alignments)
        with self.assertRaises(StopIteration):
            next(same_alignments)
        self.assertIs(type(same_alignments), PairwiseAlignments)
        self.assertIs(type(self.alignments), PairwiseAlignments)
        self.alignments.reverse()
        self.assertIs(type(self.alignments), Alignments)
        self.assertEqual(next(self.alignments), third_alignment)
        self.assertEqual(next(self.alignments), second_alignment)
        self.assertEqual(next(self.alignments), first_alignment)
        self.assertEqual(next(self.alignments), zeroth_alignment)
        with self.assertRaises(StopIteration):
            next(self.alignments)
        alignments = self.aligner.align(self.target, self.query)
        self.assertIs(type(alignments), PairwiseAlignments)
        alignments.reverse()
        self.assertIs(type(alignments), Alignments)
        self.assertEqual(alignments[0], third_alignment)
        self.assertEqual(alignments[1], second_alignment)
        self.assertEqual(alignments[2], first_alignment)
        self.assertEqual(alignments[3], zeroth_alignment)
        self.assertEqual(len(alignments), 4)
        self.assertEqual(next(alignments), third_alignment)
        self.assertEqual(next(alignments), second_alignment)
        self.assertEqual(next(alignments), first_alignment)
        self.assertEqual(next(alignments), zeroth_alignment)
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_sort(self):
        same_alignments = self.aligner.align(self.target, self.query)
        zeroth_alignment = next(same_alignments)
        first_alignment = next(same_alignments)
        second_alignment = next(same_alignments)
        third_alignment = next(same_alignments)
        with self.assertRaises(StopIteration):
            next(same_alignments)
        self.assertIs(type(self.alignments), PairwiseAlignments)
        self.alignments.sort(key=lambda alignment: str(alignment), reverse=True)
        self.assertIs(type(self.alignments), Alignments)
        alignment = next(self.alignments)
        self.assertEqual(alignment, zeroth_alignment)
        alignment = next(self.alignments)
        self.assertEqual(alignment, second_alignment)
        alignment = next(self.alignments)
        self.assertEqual(alignment, first_alignment)
        alignment = next(self.alignments)
        self.assertEqual(alignment, third_alignment)
        with self.assertRaises(StopIteration):
            next(self.alignments)

    def test_extend(self):
        self.assertAlmostEqual(self.alignments.score, 4.0)
        same_alignments = self.aligner.align(self.target, self.query)
        self.assertIs(type(self.alignments), PairwiseAlignments)
        self.assertIs(type(same_alignments), PairwiseAlignments)
        self.alignments.extend(same_alignments)
        self.assertIs(type(self.alignments), Alignments)
        self.assertIs(type(same_alignments), PairwiseAlignments)
        for alignment, text in zip(self.alignments, 2 * self.same_str):
            self.assertEqual(str(alignment), text)
        self.assertEqual(len(self.alignments), 8)
        self.assertAlmostEqual(self.alignments.score, 4.0)
        other_alignments = self.aligner.align(self.target, self.other_query)
        self.assertAlmostEqual(other_alignments.score, 4.0)
        self.assertIs(type(other_alignments), PairwiseAlignments)
        self.alignments.extend(other_alignments)
        self.assertIs(type(other_alignments), PairwiseAlignments)
        for alignment, text in zip(self.alignments, 2 * self.same_str + self.other_str):
            self.assertEqual(str(alignment), text)
        self.assertEqual(len(self.alignments), 9)
        self.assertAlmostEqual(self.alignments.score, 4.0)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
