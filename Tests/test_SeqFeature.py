# Copyright 2015-2017 by Kai Blin.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests Bio.SeqFeature.
"""
import unittest
from os import path
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, AfterPosition, BeforePosition
from Bio.SeqFeature import CompoundLocation, UnknownPosition


class TestReference(unittest.TestCase):
    """Tests for the SeqFeature.Reference class"""

    def test_eq_identical(self):
        """Test two identical references eq() to True"""
        testfile = path.join('GenBank', 'origin_line.gb')
        rec1 = SeqIO.read(testfile, 'genbank')
        rec2 = SeqIO.read(testfile, 'genbank')

        self.assertEqual(rec1.annotations['references'][0], rec1.annotations['references'][0])
        cmp1, cmp2 = rec1.annotations['references'][0], rec2.annotations['references'][0]
        self.assertEqual(cmp1, cmp2)
        self.assertNotEqual(rec1.annotations['references'][0], rec1.annotations['references'][1])
        self.assertNotEqual(rec1.annotations['references'][0], rec2.annotations['references'][1])
        self.assertEqual(rec1.annotations['references'][1], rec1.annotations['references'][1])
        self.assertEqual(rec1.annotations['references'][1], rec2.annotations['references'][1])


class TestFeatureLocation(unittest.TestCase):
    """Tests for the SeqFeature.FeatureLocation class"""

    def test_eq_identical(self):
        """Test two identical locations are equal"""
        loc1 = FeatureLocation(23, 42, 1)
        loc2 = FeatureLocation(23, 42, 1)
        self.assertEqual(loc1, loc2)

        loc1 = FeatureLocation(23, 42, -1)
        loc2 = FeatureLocation(23, 42, -1)
        self.assertEqual(loc1, loc2)

        loc1 = FeatureLocation(BeforePosition(23), AfterPosition(42), 1)
        loc2 = FeatureLocation(23, 42, 1)
        self.assertEqual(loc1, loc2)

        loc1 = FeatureLocation(23, 42, 1, 'foo', 'bar')
        loc2 = FeatureLocation(23, 42, 1, 'foo', 'bar')
        self.assertEqual(loc1, loc2)

    def test_eq_not_identical(self):
        """Test two different locations are not equal"""
        loc1 = FeatureLocation(22, 42, 1)
        loc2 = FeatureLocation(23, 42, 1)
        self.assertNotEqual(loc1, loc2)

        loc1 = FeatureLocation(23, 42, 1)
        loc2 = FeatureLocation(23, 43, 1)
        self.assertNotEqual(loc1, loc2)

        loc1 = FeatureLocation(23, 42, 1)
        loc2 = FeatureLocation(23, 42, -1)
        self.assertNotEqual(loc1, loc2)

        loc1 = FeatureLocation(23, 42, 1)
        loc2 = (23, 42, 1)
        self.assertNotEqual(loc1, loc2)

        loc1 = FeatureLocation(23, 42, 1, 'foo')
        loc2 = FeatureLocation(23, 42, 1, 'bar')
        self.assertNotEqual(loc1, loc2)

        loc1 = FeatureLocation(23, 42, 1, 'foo', 'bar')
        loc2 = FeatureLocation(23, 42, 1, 'foo', 'baz')
        self.assertNotEqual(loc1, loc2)

    def test_start_before_end(self):
        expected = "must be greater than or equal to start location"
        with self.assertRaises(ValueError) as err:
            FeatureLocation(42, 23, 1)
        assert expected in str(err.exception)

        with self.assertRaises(ValueError) as err:
            FeatureLocation(42, 0, 1)
        assert expected in str(err.exception)

        with self.assertRaises(ValueError) as err:
            FeatureLocation(BeforePosition(42), AfterPosition(23), -1)
        assert expected in str(err.exception)

        with self.assertRaises(ValueError) as err:
            FeatureLocation(42, AfterPosition(0), 1)
        assert expected in str(err.exception)

        # Features with UnknownPositions should pass check
        FeatureLocation(42, UnknownPosition())
        FeatureLocation(UnknownPosition(), 42)

        # Same start and end should pass check
        FeatureLocation(42, 42)


class TestCompoundLocation(unittest.TestCase):
    """Tests for the SeqFeature.CompoundLocation class"""

    def test_eq_identical(self):
        """Test two identical locations are equal"""
        loc1 = FeatureLocation(12, 17, 1) + FeatureLocation(23, 42, 1)
        loc2 = FeatureLocation(12, 17, 1) + FeatureLocation(23, 42, 1)
        self.assertEqual(loc1, loc2)

        loc1 = FeatureLocation(12, 17, 1) + FeatureLocation(23, 42, 1)
        loc2 = CompoundLocation([FeatureLocation(12, 17, 1), FeatureLocation(23, 42, 1)])
        self.assertEqual(loc1, loc2)

    def test_eq_not_identical(self):
        """Test two different locations are not equal"""

        loc1 = FeatureLocation(12, 17, 1) + FeatureLocation(23, 42, 1)
        loc2 = FeatureLocation(12, 17, 1) + FeatureLocation(23, 42, 1) + FeatureLocation(50, 60, 1)
        self.assertNotEqual(loc1, loc2)

        loc1 = FeatureLocation(12, 17, 1) + FeatureLocation(23, 42, 1)
        loc2 = FeatureLocation(12, 17, -1) + FeatureLocation(23, 42, -1)
        self.assertNotEqual(loc1, loc2)

        loc1 = CompoundLocation([FeatureLocation(12, 17, 1), FeatureLocation(23, 42, 1)])
        loc2 = CompoundLocation([FeatureLocation(12, 17, 1), FeatureLocation(23, 42, 1)], 'order')
        self.assertNotEqual(loc1, loc2)

        loc1 = FeatureLocation(12, 17, 1) + FeatureLocation(23, 42, 1)
        loc2 = 5
        self.assertNotEqual(loc1, loc2)
