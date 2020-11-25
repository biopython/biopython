"""Tests for SeqIO TwoBitIO module."""


import os
import random
import unittest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.SeqIO import TwoBitIO


class Parsing(unittest.TestCase):
    """Test parsing 2bit files."""

    def setUp(self):
        path = "TwoBit/sequence.fa"
        records = SeqIO.parse(path, "fasta")
        self.records = list(records)

    def test_littleendian(self):
        path = "TwoBit/sequence.littleendian.2bit"
        with open(path, "rb") as handle:
            records = TwoBitIO.TwoBitIterator(handle)
            self.assertEqual(records.byteorder, "little")
            self.assertEqual(len(self.records), len(records))
            for record1, record2 in zip(self.records, records):
                self.assertEqual(record1.id, record2.id)
                seq1 = record1.seq
                seq2 = record2.seq
                self.assertEqual(seq1, seq2)
                n = len(seq1)
                for i in range(n):
                    for j in range(i, n):
                        self.assertEqual(seq1[i:j], seq2[i:j])
                        self.assertEqual(repr(seq1[i:j]), repr(seq2[i:j]))

    def test_bigendian(self):
        path = "TwoBit/sequence.bigendian.2bit"
        with open(path, "rb") as handle:
            records = TwoBitIO.TwoBitIterator(handle)
            self.assertEqual(len(records), 6)
            self.assertEqual(records.byteorder, "big")
            for record1, record2 in zip(self.records, records):
                self.assertEqual(record1.id, record2.id)
                seq1 = record1.seq
                seq2 = record2.seq
                self.assertEqual(seq1, seq2)
                n = len(seq1)
                for i in range(n):
                    for j in range(i, n):
                        self.assertEqual(seq1[i:j], seq2[i:j])
                        self.assertEqual(repr(seq1[i:j]), repr(seq2[i:j]))

    def test_sequence_long(self):
        path = "TwoBit/sequence.long.2bit"
        with open(path, "rb") as handle:
            with self.assertRaises(ValueError) as cm:
                TwoBitIO.TwoBitIterator(handle)
            self.assertEqual(
                str(cm.exception),
                "version-1 twoBit files with 64-bit offsets for index are currently not supported",
            )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
