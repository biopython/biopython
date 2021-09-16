# Copyright 2008 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for AlignIO module."""
import unittest
import warnings


from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Align import MultipleSeqAlignment


class TestAlignIO_reading(unittest.TestCase):

    def test_reading_alignments_msf1(self):
        path = "msf/DOA_prot.msf"
        with self.assertRaisesRegex(
            ValueError,
            "GCG MSF header said alignment length 62, "
            "but 11 of 12 sequences said Len: 250",
        ):
            AlignIO.read(path, "msf")

    def test_reading_alignments_msf2(self):
        path = "msf/W_prot.msf"
        with warnings.catch_warnings(record=True) as w:
            with open(path) as handle:
                alignments = list(AlignIO.parse(handle, format='msf'))
            self.assertEqual(len(alignments), 1)
            self.assertEqual(len(alignments[0]), 11)
            # Try using the iterator with a for loop and a filename not handle
            counter = 0
            for record in AlignIO.parse(path, format='msf'):
                counter += 1
            self.assertEqual(counter, 1)
            # Try using the iterator with the next() method
            counter = 0
            alignments = AlignIO.parse(path, format="msf")
            while True:
                try:
                    alignment = next(alignments)
                except StopIteration:
                    break
                self.assertIsNotNone(alignment)
                counter += 1
            self.assertEqual(counter, 1)
            # Try a mixture of next() and list
            counter = 0
            alignments = AlignIO.parse(path, format='msf')
            alignment = next(alignments)
            counter = 1
            counter += len(list(alignments))
            self.assertEqual(counter, 1)
            # Try a mixture of next() and for loop
            alignments = AlignIO.parse(path, format='msf')
            alignment = next(alignments)
            counter = 1
            for alignment in alignments:
                counter += 1
            self.assertEqual(counter, 1)
            # Check Bio.AlignIO.read(...)
            with open(path) as handle:
                alignment = AlignIO.read(handle, format='msf')
            self.assertIsInstance(alignment, MultipleSeqAlignment)
            self.assertEqual(len(alignment), 11)
            self.assertEqual(alignment.get_alignment_length(), 99)
        warning_msgs = {str(_.message) for _ in w}
        self.assertIn(
            "One of more alignment sequences were truncated and have been gap padded",
            warning_msgs,
        )
        # Compare each sequence column
        self.assertEqual(alignment[:, 0], "GGGGGGGGGGG")
        self.assertEqual(alignment[:, 1], "LLLLLLLLLLL")
        self.assertEqual(alignment[:, 2], "TTTTTTTTTTT")
        self.assertEqual(alignment[:, 3], "PPPPPPPPPPP")
        self.assertEqual(alignment[:, 4], "FFFFFFSSSSS")
        self.assertEqual(alignment[:, -1], "LLLLLL----L")
        summary = AlignInfo.SummaryInfo(alignment)
        dumb_consensus = summary.dumb_consensus()


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
