# Copyright 2008 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Align.stockholm module."""
import unittest

from Bio.Align import stockholm


class TestAlignIO_reading(unittest.TestCase):

    def test_reading_alignments_stockholm1(self):
        path = "Stockholm/simple.sth"
        with open(path) as stream:
            alignments = stockholm.AlignmentIterator(path)
            alignment = next(alignments)
            self.assertRaises(StopIteration, next, alignments)

    def test_reading_alignments_stockholm2(self):
        path = "Stockholm/funny.sth"
        with open(path) as stream:
            alignments = stockholm.AlignmentIterator(path)
            with self.assertRaises(ValueError) as cm:
                alignment = next(alignments)
            self.assertEqual(str(cm.exception), "Start and end of sequence O83071 are not consistent with sequence length")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
