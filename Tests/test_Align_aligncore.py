"""Unit tests for PrintedAlignmentParser.feed() offset handling.

These tests verify normal parsing, valid offsets, and rejection of offsets
outside the bounds of the input bytes object.
"""

import unittest

from Bio.Align._aligncore import PrintedAlignmentParser


class TestPrintedAlignmentParser(unittest.TestCase):
    def setUp(self):
        self.parser = PrintedAlignmentParser(b";")

    def test_feed_without_offset(self):
        result = self.parser.feed(b"AC-GT;")
        self.assertEqual(result, (5, b"ACGT"))

    def test_valid_offset_parses_remaining_bytes(self):
        result = self.parser.feed(b"XXAC-GT;", 2)
        self.assertEqual(result, (5, b"ACGT"))

    def test_offset_equal_to_length_is_accepted(self):
        result = self.parser.feed(b"ACGT", 4)
        self.assertEqual(result, (0, b""))

    def test_negative_offset_rejected(self):
        with self.assertRaisesRegex(ValueError, "offset -1 is out of range"):
            self.parser.feed(b"ACGT", -1)

    def test_large_offset_rejected(self):
        with self.assertRaisesRegex(ValueError, "offset 100 is out of range"):
            self.parser.feed(b"ACGT", 100)

    def test_offset_one_past_end_rejected(self):
        with self.assertRaisesRegex(ValueError, "offset 5 is out of range"):
            self.parser.feed(b"ACGT", 5)


if __name__ == "__main__":
    unittest.main()
