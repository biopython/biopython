# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for dna-2 encoding scheme."""

import importlib
import itertools
import unittest

_nucleotides = set("ATCG")


class TestCodec(unittest.TestCase):
    def setUp(self):
        importlib.import_module("Bio.codecs")

    def test_encode_and_decode(self):
        for length in range(1, 10):
            for sequence in itertools.product(_nucleotides, repeat=length):
                sequence = "".join(sequence)
                decoded = sequence.encode("dna-2").decode("dna-2")
                self.assertEqual(sequence, decoded)
                decoded = sequence.lower().encode("dna-2").decode("dna-2")
                self.assertEqual(sequence, decoded)

    def test_strict_errors(self):
        def encode_bad_sequence() -> None:
            """Call encode with a bad DNA sequence."""
            "actgatcgatcgatgctagctdacatgcttagct".encode("dna-2")

        self.assertRaises(ValueError, encode_bad_sequence)

    def test_ignore_errors(self):
        sequence = "actgatcgatcgatgctagctdacatgcttagct"
        expected = "actgatcgatcgatgctagctacatgcttagct".upper()
        decoded = sequence.encode("dna-2", errors="ignore").decode("dna-2")

        self.assertEqual(expected, decoded)
