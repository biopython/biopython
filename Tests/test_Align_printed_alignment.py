# Copyright 2026 by Nicola Trinca. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for PrintedAlignmentParser and Alignment.parse_printed_alignment."""

import unittest

import numpy as np

from Bio.Align import Alignment, _aligncore


class TestPrintedAlignmentParser(unittest.TestCase):
    def test_parse_printed_alignment_valid_input(self):
        lines = [b"TAGGCATACGTG", b"AACG--TACGT-", b"-ACGCATACTTG"]
        sequences, coordinates = Alignment.parse_printed_alignment(lines)
        self.assertEqual(sequences, [b"TAGGCATACGTG", b"AACGTACGT", b"ACGCATACTTG"])
        expected = np.array(
            [[0, 1, 4, 6, 11, 12], [0, 1, 4, 4, 9, 9], [0, 0, 3, 5, 10, 11]],
            np.intp,
        )
        self.assertEqual(coordinates.dtype, np.intp)
        np.testing.assert_array_equal(coordinates, expected)

    def test_feed_valid_offset(self):
        parser = _aligncore.PrintedAlignmentParser()
        nbytes, sequence = parser.feed(b"..AC-GT\njunk", 2)
        self.assertEqual(nbytes, 5)
        self.assertEqual(sequence, b"ACGT")

    def test_feed_rejects_large_offset(self):
        parser = _aligncore.PrintedAlignmentParser()
        with self.assertRaises(ValueError):
            parser.feed(b"ACGT", 100000)

    def test_feed_rejects_offset_equal_to_line_length(self):
        parser = _aligncore.PrintedAlignmentParser()
        with self.assertRaises(ValueError):
            parser.feed(b"ACGT", 4)

    def test_feed_rejects_negative_offset(self):
        parser = _aligncore.PrintedAlignmentParser()
        for offset in (-1, -(1 << 40)):
            with self.subTest(offset=offset), self.assertRaises(ValueError):
                parser.feed(b"ACGT", offset)

    def test_feed_rejects_very_large_offset(self):
        # 1 << 40 segfaulted the interpreter before the fix (issue #5272)
        parser = _aligncore.PrintedAlignmentParser()
        with self.assertRaises(ValueError):
            parser.feed(b"ACGT", 1 << 40)

    def test_feed_rejects_ragged_rows(self):
        parser = _aligncore.PrintedAlignmentParser()
        parser.feed(b"ACGT--ACGTACGT--AC")  # 18 columns
        with self.assertRaises(ValueError):
            parser.feed(b"ACGTTTACGT")  # 10 columns

    def test_parse_printed_alignment_rejects_ragged_input(self):
        lines = [b"ACGT--ACGTACGT--AC", b"ACGTTTACGT"]  # 18 and 10 columns
        with self.assertRaises(ValueError):
            Alignment.parse_printed_alignment(lines)


if __name__ == "__main__":
    unittest.main(verbosity=2)
