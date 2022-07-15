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


try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install numpy if you want to use Bio.Align.interfaces."
    ) from None


class TestAlignments(unittest.TestCase):
    def setUp(self):
        self.path = "Emboss/needle.txt"
        self.alignments = AlignmentIterator(self.path)

    def test_repr(self):
        representation = "<Bio.Align.emboss.AlignmentIterator object at %s>" % hex(id(self.alignments))
        self.assertFalse(self.alignments._loaded)
        self.assertEqual(repr(self.alignments), representation)
        alignments = self.alignments[:]
        self.assertTrue(self.alignments._loaded)
        self.assertEqual(repr(self.alignments), representation)
        representation = "<Bio.Align.interfaces.AlignmentIterator object at %s>" % hex(id(alignments))
        self.assertEqual(repr(alignments), representation)

    def test_lt_false(self):
        same_alignments = AlignmentIterator(self.path)
        self.assertFalse(self.alignments._loaded)
        self.assertFalse(same_alignments._loaded)
        result = (self.alignments < same_alignments)
        self.assertTrue(self.alignments._loaded)
        self.assertTrue(same_alignments._loaded)
        self.assertFalse(result)

    def test_lt_true(self):
        other_alignments = AlignmentIterator("Emboss/water.txt")
        self.assertFalse(self.alignments._loaded)
        self.assertFalse(other_alignments._loaded)
        result = (self.alignments < other_alignments)
        self.assertTrue(self.alignments._loaded)
        self.assertTrue(other_alignments._loaded)
        self.assertTrue(result)

    def test_le_false(self):
        other_alignments = AlignmentIterator("Emboss/water.txt")
        self.assertFalse(self.alignments._loaded)
        self.assertFalse(other_alignments._loaded)
        result = (other_alignments <= self.alignments)
        self.assertTrue(self.alignments._loaded)
        self.assertTrue(other_alignments._loaded)
        self.assertFalse(result)

    def test_le_true(self):
        same_alignments = AlignmentIterator(self.path)
        self.assertFalse(self.alignments._loaded)
        self.assertFalse(same_alignments._loaded)
        result = (self.alignments <= same_alignments)
        self.assertTrue(self.alignments._loaded)
        self.assertTrue(same_alignments._loaded)
        self.assertTrue(result)

    def test_eq_false(self):
        other_alignments = AlignmentIterator("Emboss/water.txt")
        self.assertFalse(self.alignments._loaded)
        self.assertFalse(other_alignments._loaded)
        result = (other_alignments == self.alignments)
        self.assertTrue(self.alignments._loaded)
        self.assertTrue(other_alignments._loaded)
        self.assertFalse(result)

    def test_eq_true(self):
        same_alignments = AlignmentIterator(self.path)
        self.assertFalse(self.alignments._loaded)
        self.assertFalse(same_alignments._loaded)
        result = (self.alignments == same_alignments)
        self.assertTrue(self.alignments._loaded)
        self.assertTrue(same_alignments._loaded)
        self.assertTrue(result)

    def test_gt_false(self):
        same_alignments = AlignmentIterator(self.path)
        self.assertFalse(self.alignments._loaded)
        self.assertFalse(same_alignments._loaded)
        result = (same_alignments > self.alignments)
        self.assertTrue(self.alignments._loaded)
        self.assertTrue(same_alignments._loaded)
        self.assertFalse(result)

    def test_gt_true(self):
        other_alignments = AlignmentIterator("Emboss/water.txt")
        self.assertFalse(self.alignments._loaded)
        self.assertFalse(other_alignments._loaded)
        result = (other_alignments > self.alignments)
        self.assertTrue(self.alignments._loaded)
        self.assertTrue(other_alignments._loaded)
        self.assertTrue(result)

    def test_ge_false(self):
        other_alignments = AlignmentIterator("Emboss/water.txt")
        self.assertFalse(self.alignments._loaded)
        self.assertFalse(other_alignments._loaded)
        result = (self.alignments >= other_alignments)
        self.assertTrue(self.alignments._loaded)
        self.assertTrue(other_alignments._loaded)
        self.assertFalse(result)

    def test_ge_true(self):
        same_alignments = AlignmentIterator(self.path)
        self.assertFalse(self.alignments._loaded)
        self.assertFalse(same_alignments._loaded)
        result = (self.alignments >= same_alignments)
        self.assertTrue(self.alignments._loaded)
        self.assertTrue(same_alignments._loaded)
        self.assertTrue(result)

    def test_contains_false(self):
        other_alignments = AlignmentIterator("Emboss/water.txt")
        alignment = next(other_alignments)
        self.assertFalse(self.alignments._loaded)
        result = (alignment in self.alignments)
        self.assertTrue(self.alignments._loaded)
        self.assertFalse(result)

    def test_contains_true(self):
        same_alignments = AlignmentIterator(self.path)
        alignment = next(same_alignments)
        self.assertFalse(self.alignments._loaded)
        result = (alignment in self.alignments)
        self.assertTrue(self.alignments._loaded)
        self.assertTrue(result)

    def test_len(self):
        self.assertFalse(self.alignments._loaded)
        length = len(self.alignments)
        self.assertTrue(self.alignments._loaded)
        self.assertEqual(length, 5)

    def test_getitem(self):
        self.assertFalse(self.alignments._loaded)
        alignment = self.alignments[3]
        self.assertTrue(self.alignments._loaded)
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
        self.assertFalse(same_alignments._loaded)

    def test_setitem(self):
        same_alignments = AlignmentIterator(self.path)
        zeroth_alignment = next(same_alignments)
        first_alignment = next(same_alignments)
        second_alignment = next(same_alignments)
        third_alignment = next(same_alignments)
        fourth_alignment = next(same_alignments)
        self.assertFalse(same_alignments._loaded)
        self.assertFalse(self.alignments._loaded)
        self.alignments[3] = first_alignment
        self.assertTrue(self.alignments._loaded)
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
        self.assertFalse(same_alignments._loaded)
        self.assertFalse(self.alignments._loaded)
        del self.alignments[3]
        self.assertTrue(self.alignments._loaded)
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
        same_repr = ("<Bio.Align.Alignment object (2 rows x 124 columns) at ",
                     "<Bio.Align.Alignment object (2 rows x 119 columns) at ",
                     "<Bio.Align.Alignment object (2 rows x 120 columns) at ",
                     "<Bio.Align.Alignment object (2 rows x 118 columns) at ",
                     "<Bio.Align.Alignment object (2 rows x 125 columns) at ",
                    )
        other_repr = ("<Bio.Align.Alignment object (2 rows x 131 columns) at ",
                     )
        same_alignments = AlignmentIterator(self.path)
        self.assertEqual(self.alignments.metadata['Align_format'], 'srspair')
        self.assertEqual(self.alignments.metadata['Program'], 'needle')
        self.assertEqual(self.alignments.metadata['Rundate'], 'Sun 27 Apr 2007 17:20:35')
        self.assertEqual(self.alignments.metadata['Command line'], 'needle [-asequence] Spo0F.faa [-bsequence] paired_r.faa -sformat2 pearson')
        self.assertEqual(self.alignments.metadata['Report_file'], 'ref_rec .needle')
        self.assertEqual(len(self.alignments.metadata), 5)
        self.assertFalse(self.alignments._loaded)
        self.assertFalse(same_alignments._loaded)
        alignments = self.alignments + same_alignments
        self.assertTrue(self.alignments._loaded)
        self.assertTrue(same_alignments._loaded)
        self.assertTrue(alignments._loaded)
        for alignment, representation in zip(alignments, 2 * same_repr):
            self.assertTrue(repr(alignment).startswith(representation))
        self.assertEqual(len(alignments), 10)
        with self.assertRaises(AttributeError):
            alignments.metadata
        other_alignments = AlignmentIterator("Emboss/water.txt")
        self.assertEqual(other_alignments.metadata['Align_format'], 'srspair')
        self.assertEqual(other_alignments.metadata['Program'], 'water')
        self.assertEqual(other_alignments.metadata['Rundate'], 'Wed Jan 16 17:23:19 2002')
        self.assertEqual(other_alignments.metadata['Report_file'], 'stdout')
        self.assertEqual(len(other_alignments.metadata), 4)

        self.assertFalse(other_alignments._loaded)
        alignments = self.alignments + other_alignments
        self.assertTrue(other_alignments._loaded)
        self.assertTrue(alignments._loaded)
        for alignment, representation in zip(alignments, same_repr + other_repr):
            self.assertTrue(repr(alignment).startswith(representation))
        self.assertEqual(len(alignments), 6)
        with self.assertRaises(AttributeError):
            alignments.metadata

    def test_radd(self):
        same_repr = ("<Bio.Align.Alignment object (2 rows x 124 columns) at ",
                     "<Bio.Align.Alignment object (2 rows x 119 columns) at ",
                     "<Bio.Align.Alignment object (2 rows x 120 columns) at ",
                     "<Bio.Align.Alignment object (2 rows x 118 columns) at ",
                     "<Bio.Align.Alignment object (2 rows x 125 columns) at ",
                    )
        other_repr = ("<Bio.Align.Alignment object (2 rows x 131 columns) at ",
                     )
        self.assertEqual(self.alignments.metadata['Align_format'], 'srspair')
        self.assertEqual(self.alignments.metadata['Program'], 'needle')
        self.assertEqual(self.alignments.metadata['Rundate'], 'Sun 27 Apr 2007 17:20:35')
        self.assertEqual(self.alignments.metadata['Command line'], 'needle [-asequence] Spo0F.faa [-bsequence] paired_r.faa -sformat2 pearson')
        self.assertEqual(self.alignments.metadata['Report_file'], 'ref_rec .needle')
        self.assertEqual(len(self.alignments.metadata), 5)
        other_alignments = AlignmentIterator("Emboss/water.txt")
        self.assertEqual(other_alignments.metadata['Align_format'], 'srspair')
        self.assertEqual(other_alignments.metadata['Program'], 'water')
        self.assertEqual(other_alignments.metadata['Rundate'], 'Wed Jan 16 17:23:19 2002')
        self.assertEqual(other_alignments.metadata['Report_file'], 'stdout')
        self.assertEqual(len(other_alignments.metadata), 4)
        self.assertFalse(self.alignments._loaded)
        self.assertFalse(other_alignments._loaded)
        alignments = list(other_alignments) + self.alignments
        self.assertTrue(self.alignments._loaded)
        self.assertTrue(other_alignments._loaded)
        self.assertTrue(alignments._loaded)
        for alignment, representation in zip(alignments, other_repr + same_repr):
            self.assertTrue(repr(alignment).startswith(representation))
        self.assertEqual(len(alignments), 6)
        with self.assertRaises(AttributeError):
            alignments.metadata

    def test_iadd(self):
        same_repr = ("<Bio.Align.Alignment object (2 rows x 124 columns) at ",
                     "<Bio.Align.Alignment object (2 rows x 119 columns) at ",
                     "<Bio.Align.Alignment object (2 rows x 120 columns) at ",
                     "<Bio.Align.Alignment object (2 rows x 118 columns) at ",
                     "<Bio.Align.Alignment object (2 rows x 125 columns) at ",
                    )
        other_repr = ("<Bio.Align.Alignment object (2 rows x 131 columns) at ",
                     )
        same_alignments = AlignmentIterator(self.path)
        self.assertEqual(self.alignments.metadata['Align_format'], 'srspair')
        self.assertEqual(self.alignments.metadata['Program'], 'needle')
        self.assertEqual(self.alignments.metadata['Rundate'], 'Sun 27 Apr 2007 17:20:35')
        self.assertEqual(self.alignments.metadata['Command line'], 'needle [-asequence] Spo0F.faa [-bsequence] paired_r.faa -sformat2 pearson')
        self.assertEqual(self.alignments.metadata['Report_file'], 'ref_rec .needle')
        self.assertEqual(len(self.alignments.metadata), 5)
        self.assertFalse(self.alignments._loaded)
        self.assertFalse(same_alignments._loaded)
        self.alignments += same_alignments
        self.assertTrue(self.alignments._loaded)
        self.assertFalse(same_alignments._loaded)
        for alignment, representation in zip(self.alignments, 2 * same_repr):
            self.assertTrue(repr(alignment).startswith(representation))
        self.assertEqual(len(self.alignments), 10)
        self.assertEqual(self.alignments.metadata['Align_format'], 'srspair')
        self.assertEqual(self.alignments.metadata['Program'], 'needle')
        self.assertEqual(self.alignments.metadata['Rundate'], 'Sun 27 Apr 2007 17:20:35')
        self.assertEqual(self.alignments.metadata['Command line'], 'needle [-asequence] Spo0F.faa [-bsequence] paired_r.faa -sformat2 pearson')
        self.assertEqual(self.alignments.metadata['Report_file'], 'ref_rec .needle')
        other_alignments = AlignmentIterator("Emboss/water.txt")
        self.assertEqual(other_alignments.metadata['Align_format'], 'srspair')
        self.assertEqual(other_alignments.metadata['Program'], 'water')
        self.assertEqual(other_alignments.metadata['Rundate'], 'Wed Jan 16 17:23:19 2002')
        self.assertEqual(other_alignments.metadata['Report_file'], 'stdout')
        self.assertEqual(len(other_alignments.metadata), 4)

        self.assertFalse(other_alignments._loaded)
        self.alignments += other_alignments
        self.assertFalse(other_alignments._loaded)
        for alignment, representation in zip(self.alignments, 2 * same_repr + other_repr):
            self.assertTrue(repr(alignment).startswith(representation))
        self.assertEqual(len(self.alignments), 11)
        self.assertEqual(self.alignments.metadata['Align_format'], 'srspair')
        self.assertEqual(self.alignments.metadata['Program'], 'needle')
        self.assertEqual(self.alignments.metadata['Rundate'], 'Sun 27 Apr 2007 17:20:35')
        self.assertEqual(self.alignments.metadata['Command line'], 'needle [-asequence] Spo0F.faa [-bsequence] paired_r.faa -sformat2 pearson')
        self.assertEqual(self.alignments.metadata['Report_file'], 'ref_rec .needle')

    def test_mul(self):
        same_repr = ("<Bio.Align.Alignment object (2 rows x 124 columns) at ",
                     "<Bio.Align.Alignment object (2 rows x 119 columns) at ",
                     "<Bio.Align.Alignment object (2 rows x 120 columns) at ",
                     "<Bio.Align.Alignment object (2 rows x 118 columns) at ",
                     "<Bio.Align.Alignment object (2 rows x 125 columns) at ",
                    )
        self.assertEqual(self.alignments.metadata['Align_format'], 'srspair')
        self.assertEqual(self.alignments.metadata['Program'], 'needle')
        self.assertEqual(self.alignments.metadata['Rundate'], 'Sun 27 Apr 2007 17:20:35')
        self.assertEqual(self.alignments.metadata['Command line'], 'needle [-asequence] Spo0F.faa [-bsequence] paired_r.faa -sformat2 pearson')
        self.assertEqual(self.alignments.metadata['Report_file'], 'ref_rec .needle')
        self.assertEqual(len(self.alignments.metadata), 5)
        self.assertFalse(self.alignments._loaded)
        alignments = 3 * self.alignments
        self.assertTrue(self.alignments._loaded)
        self.assertTrue(alignments._loaded)
        for alignment, representation in zip(alignments, 3 * same_repr):
            self.assertTrue(repr(alignment).startswith(representation))
        self.assertEqual(len(alignments), 15)
        with self.assertRaises(AttributeError):
            alignments.metadata

    def test_rmul(self):
        same_repr = ("<Bio.Align.Alignment object (2 rows x 124 columns) at ",
                     "<Bio.Align.Alignment object (2 rows x 119 columns) at ",
                     "<Bio.Align.Alignment object (2 rows x 120 columns) at ",
                     "<Bio.Align.Alignment object (2 rows x 118 columns) at ",
                     "<Bio.Align.Alignment object (2 rows x 125 columns) at ",
                    )
        self.assertEqual(self.alignments.metadata['Align_format'], 'srspair')
        self.assertEqual(self.alignments.metadata['Program'], 'needle')
        self.assertEqual(self.alignments.metadata['Rundate'], 'Sun 27 Apr 2007 17:20:35')
        self.assertEqual(self.alignments.metadata['Command line'], 'needle [-asequence] Spo0F.faa [-bsequence] paired_r.faa -sformat2 pearson')
        self.assertEqual(self.alignments.metadata['Report_file'], 'ref_rec .needle')
        self.assertEqual(len(self.alignments.metadata), 5)
        self.assertFalse(self.alignments._loaded)
        alignments = self.alignments * 3
        self.assertTrue(self.alignments._loaded)
        self.assertTrue(alignments._loaded)
        for alignment, representation in zip(alignments, 3 * same_repr):
            self.assertTrue(repr(alignment).startswith(representation))
        self.assertEqual(len(alignments), 15)
        with self.assertRaises(AttributeError):
            alignments.metadata

    def test_imul(self):
        same_repr = ("<Bio.Align.Alignment object (2 rows x 124 columns) at ",
                     "<Bio.Align.Alignment object (2 rows x 119 columns) at ",
                     "<Bio.Align.Alignment object (2 rows x 120 columns) at ",
                     "<Bio.Align.Alignment object (2 rows x 118 columns) at ",
                     "<Bio.Align.Alignment object (2 rows x 125 columns) at ",
                    )
        self.assertEqual(self.alignments.metadata['Align_format'], 'srspair')
        self.assertEqual(self.alignments.metadata['Program'], 'needle')
        self.assertEqual(self.alignments.metadata['Rundate'], 'Sun 27 Apr 2007 17:20:35')
        self.assertEqual(self.alignments.metadata['Command line'], 'needle [-asequence] Spo0F.faa [-bsequence] paired_r.faa -sformat2 pearson')
        self.assertEqual(self.alignments.metadata['Report_file'], 'ref_rec .needle')
        self.assertEqual(len(self.alignments.metadata), 5)
        self.assertFalse(self.alignments._loaded)
        self.alignments *= 3
        self.assertTrue(self.alignments._loaded)
        for alignment, representation in zip(self.alignments, 3 * same_repr):
            self.assertTrue(repr(alignment).startswith(representation))
        self.assertEqual(len(self.alignments), 15)
        self.assertEqual(self.alignments.metadata['Align_format'], 'srspair')
        self.assertEqual(self.alignments.metadata['Program'], 'needle')
        self.assertEqual(self.alignments.metadata['Rundate'], 'Sun 27 Apr 2007 17:20:35')
        self.assertEqual(self.alignments.metadata['Command line'], 'needle [-asequence] Spo0F.faa [-bsequence] paired_r.faa -sformat2 pearson')
        self.assertEqual(self.alignments.metadata['Report_file'], 'ref_rec .needle')
        self.assertEqual(len(self.alignments.metadata), 5)

    def test_append(self):
        same_alignments = AlignmentIterator(self.path)
        zeroth_alignment = next(same_alignments)
        first_alignment = next(same_alignments)
        second_alignment = next(same_alignments)
        third_alignment = next(same_alignments)
        fourth_alignment = next(same_alignments)
        self.assertFalse(same_alignments._loaded)
        self.assertFalse(self.alignments._loaded)
        self.alignments.append(second_alignment)
        self.assertTrue(self.alignments._loaded)
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
        self.assertFalse(same_alignments._loaded)
        self.assertFalse(self.alignments._loaded)
        self.alignments.insert(1, fourth_alignment)
        self.assertTrue(self.alignments._loaded)
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
        self.assertFalse(self.alignments._loaded)
        alignment = self.alignments.pop(1)
        self.assertTrue(self.alignments._loaded)
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
        self.assertFalse(same_alignments._loaded)

    def test_remove(self):
        same_alignments = AlignmentIterator(self.path)
        zeroth_alignment = next(same_alignments)
        first_alignment = next(same_alignments)
        second_alignment = next(same_alignments)
        third_alignment = next(same_alignments)
        fourth_alignment = next(same_alignments)
        with self.assertRaises(StopIteration):
            next(same_alignments)
        self.assertFalse(same_alignments._loaded)
        self.assertFalse(self.alignments._loaded)
        result = self.alignments.remove(first_alignment)
        self.assertTrue(self.alignments._loaded)
        self.assertIsNone(result)
        self.assertEqual(zeroth_alignment, next(self.alignments))
        self.assertEqual(second_alignment, next(self.alignments))
        self.assertEqual(third_alignment, next(self.alignments))
        self.assertEqual(fourth_alignment, next(self.alignments))
        with self.assertRaises(StopIteration):
            next(self.alignments)

    def test_clear(self):
        self.assertFalse(self.alignments._loaded)
        result = self.alignments.clear()
        self.assertIsNone(result)
        self.assertFalse(self.alignments._loaded)
        with self.assertRaises(AttributeError):
            self.alignments._stream
        with self.assertRaises(StopIteration):
            next(self.alignments)
        self.assertEqual(len(self.alignments), 0)

    def test_copy(self):
        self.assertFalse(self.alignments._loaded)
        alignments = self.alignments.copy()
        self.assertTrue(self.alignments._loaded)
        self.assertTrue(alignments._loaded)
        self.assertEqual(len(self.alignments), 5)
        for alignment in self.alignments:
            self.assertEqual(alignment, next(alignments))
        with self.assertRaises(StopIteration):
            next(self.alignments)
        with self.assertRaises(StopIteration):
            next(alignments)
        self.assertEqual(self.alignments.metadata, alignments.metadata)

    def test_count(self):
        same_alignments = AlignmentIterator(self.path)
        zeroth_alignment = next(same_alignments)
        first_alignment = next(same_alignments)
        second_alignment = next(same_alignments)
        third_alignment = next(same_alignments)
        fourth_alignment = next(same_alignments)
        with self.assertRaises(StopIteration):
            next(same_alignments)
        self.assertFalse(same_alignments._loaded)
        self.assertFalse(self.alignments._loaded)
        count = self.alignments.count(first_alignment)
        self.assertTrue(self.alignments._loaded)
        self.assertEqual(count, 1)
        count = self.alignments.count(second_alignment)
        self.assertEqual(count, 1)
        self.alignments.insert(3, first_alignment)
        count = self.alignments.count(first_alignment)
        self.assertEqual(count, 2)
        same_alignments.insert(3, first_alignment)
        self.assertTrue(same_alignments._loaded)
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
        self.assertFalse(same_alignments._loaded)
        self.assertFalse(self.alignments._loaded)
        index = self.alignments.index(third_alignment)
        self.assertTrue(self.alignments._loaded)
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
        self.assertFalse(same_alignments._loaded)
        self.assertFalse(self.alignments._loaded)
        self.alignments.reverse()
        self.assertTrue(self.alignments._loaded)
        self.assertEqual(next(self.alignments), fourth_alignment)
        self.assertEqual(next(self.alignments), third_alignment)
        self.assertEqual(next(self.alignments), second_alignment)
        self.assertEqual(next(self.alignments), first_alignment)
        self.assertEqual(next(self.alignments), zeroth_alignment)
        with self.assertRaises(StopIteration):
            next(self.alignments)
        alignments = AlignmentIterator(self.path)
        self.assertFalse(alignments._loaded)
        alignments.reverse()
        self.assertTrue(alignments._loaded)
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
        self.assertFalse(self.alignments._loaded)
        self.alignments.sort(key=lambda alignment: alignment.score, reverse=True)
        self.assertTrue(self.alignments._loaded)
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
        same_repr = ("<Bio.Align.Alignment object (2 rows x 124 columns) at ",
                     "<Bio.Align.Alignment object (2 rows x 119 columns) at ",
                     "<Bio.Align.Alignment object (2 rows x 120 columns) at ",
                     "<Bio.Align.Alignment object (2 rows x 118 columns) at ",
                     "<Bio.Align.Alignment object (2 rows x 125 columns) at ",
                    )
        other_repr = ("<Bio.Align.Alignment object (2 rows x 131 columns) at ",
                     )
        same_alignments = AlignmentIterator(self.path)
        self.assertEqual(self.alignments.metadata['Align_format'], 'srspair')
        self.assertEqual(self.alignments.metadata['Program'], 'needle')
        self.assertEqual(self.alignments.metadata['Rundate'], 'Sun 27 Apr 2007 17:20:35')
        self.assertEqual(self.alignments.metadata['Command line'], 'needle [-asequence] Spo0F.faa [-bsequence] paired_r.faa -sformat2 pearson')
        self.assertEqual(self.alignments.metadata['Report_file'], 'ref_rec .needle')
        self.assertEqual(len(self.alignments.metadata), 5)
        self.assertFalse(self.alignments._loaded)
        self.assertFalse(same_alignments._loaded)
        self.alignments.extend(same_alignments)
        self.assertTrue(self.alignments._loaded)
        self.assertFalse(same_alignments._loaded)
        for alignment, representation in zip(self.alignments, 2 * same_repr):
            self.assertTrue(repr(alignment).startswith(representation))
        self.assertEqual(len(self.alignments), 10)
        self.assertEqual(self.alignments.metadata['Align_format'], 'srspair')
        self.assertEqual(self.alignments.metadata['Program'], 'needle')
        self.assertEqual(self.alignments.metadata['Rundate'], 'Sun 27 Apr 2007 17:20:35')
        self.assertEqual(self.alignments.metadata['Command line'], 'needle [-asequence] Spo0F.faa [-bsequence] paired_r.faa -sformat2 pearson')
        self.assertEqual(self.alignments.metadata['Report_file'], 'ref_rec .needle')
        other_alignments = AlignmentIterator("Emboss/water.txt")
        self.assertEqual(other_alignments.metadata['Align_format'], 'srspair')
        self.assertEqual(other_alignments.metadata['Program'], 'water')
        self.assertEqual(other_alignments.metadata['Rundate'], 'Wed Jan 16 17:23:19 2002')
        self.assertEqual(other_alignments.metadata['Report_file'], 'stdout')
        self.assertEqual(len(other_alignments.metadata), 4)

        self.assertFalse(other_alignments._loaded)
        self.alignments.extend(other_alignments)
        self.assertFalse(other_alignments._loaded)
        for alignment, representation in zip(self.alignments, 2 * same_repr + other_repr):
            self.assertTrue(repr(alignment).startswith(representation))
        self.assertEqual(len(self.alignments), 11)
        self.assertEqual(self.alignments.metadata['Align_format'], 'srspair')
        self.assertEqual(self.alignments.metadata['Program'], 'needle')
        self.assertEqual(self.alignments.metadata['Rundate'], 'Sun 27 Apr 2007 17:20:35')
        self.assertEqual(self.alignments.metadata['Command line'], 'needle [-asequence] Spo0F.faa [-bsequence] paired_r.faa -sformat2 pearson')
        self.assertEqual(self.alignments.metadata['Report_file'], 'ref_rec .needle')


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
