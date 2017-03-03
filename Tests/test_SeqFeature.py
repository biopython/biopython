# Copyright 2015 by Kai Blin.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests Bio.SeqFeature.
"""
import unittest
from os import path
from Bio import SeqIO


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
