"""Tests for SeqIO TwoBitIO module."""


import os
import random
import unittest

from Bio.Seq import Seq, MutableSeq, UndefinedSequenceError
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
        with open(path, "rb") as stream:
            records = TwoBitIO.TwoBitIterator(stream)
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
        with open(path, "rb") as stream:
            records = TwoBitIO.TwoBitIterator(stream)
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
        with open(path, "rb") as stream:
            with self.assertRaises(ValueError) as cm:
                TwoBitIO.TwoBitIterator(stream)
            self.assertEqual(
                str(cm.exception),
                "version-1 twoBit files with 64-bit offsets for index are currently not supported",
            )


class Comparisons(unittest.TestCase):
    """Test comparisons of sequences read from 2bit files to Seq and other objects."""

    def setUp(self):
        path = "TwoBit/sequence.bigendian.2bit"
        self.stream = open(path, "rb")
        records = TwoBitIO.TwoBitIterator(self.stream)
        record1 = next(records)
        record2 = next(records)
        self.seq1 = record1.seq
        self.seq2 = record2.seq

    def tearDown(self):
        self.stream.close()

    def test_eq(self):
        seq1 = self.seq1
        seq2 = self.seq2
        self.assertEqual(seq1[2:4], seq2[3:5])
        self.assertEqual(seq1[2:4], "AT")
        self.assertEqual(seq2[3:5], "AT")
        self.assertEqual(seq1[2:4], b"AT")
        self.assertEqual(seq2[3:5], b"AT")
        self.assertEqual(seq1[2:4], Seq("AT"))
        self.assertEqual(seq2[3:5], Seq("AT"))
        self.assertEqual(seq1[2:4], MutableSeq("AT"))
        self.assertEqual(seq2[3:5], MutableSeq("AT"))
        self.assertEqual(seq2[3:5], seq1[2:4])
        self.assertEqual("AT", seq1[2:4])
        self.assertEqual("AT", seq2[3:5])
        self.assertEqual(b"AT", seq1[2:4])
        self.assertEqual(b"AT", seq2[3:5])
        self.assertEqual(Seq("AT"), seq1[2:4])
        self.assertEqual(Seq("AT"), seq2[3:5])
        self.assertEqual(MutableSeq("AT"), seq1[2:4])
        self.assertEqual(MutableSeq("AT"), seq2[3:5])
        with self.assertRaises(UndefinedSequenceError):
            seq1 == Seq(None, len(seq1))
        with self.assertRaises(UndefinedSequenceError):
            seq2 == Seq(None, len(seq2))
        with self.assertRaises(UndefinedSequenceError):
            seq1 == Seq(None, 10)
        with self.assertRaises(UndefinedSequenceError):
            seq2 == Seq(None, 10)
        with self.assertRaises(UndefinedSequenceError):
            Seq(None, len(seq1)) == seq1
        with self.assertRaises(UndefinedSequenceError):
            Seq(None, len(seq2)) == seq2
        with self.assertRaises(UndefinedSequenceError):
            Seq(None, 10) == seq1
        with self.assertRaises(UndefinedSequenceError):
            Seq(None, 10) == seq2

    def test_ne(self):
        seq1 = self.seq1
        seq2 = self.seq2
        self.assertNotEqual(seq1, seq2)
        self.assertNotEqual(seq1, "AT")
        self.assertNotEqual(seq2, "AT")
        self.assertNotEqual(seq1, b"AT")
        self.assertNotEqual(seq2, b"AT")
        self.assertNotEqual(seq1, Seq("AT"))
        self.assertNotEqual(seq2, Seq("AT"))
        self.assertNotEqual(seq1, MutableSeq("AT"))
        self.assertNotEqual(seq2, MutableSeq("AT"))
        self.assertNotEqual(seq1[2:4], "CG")
        self.assertNotEqual(seq2[3:5], "CG")
        self.assertNotEqual(seq1[2:4], b"CG")
        self.assertNotEqual(seq2[3:5], b"CG")
        self.assertNotEqual(seq1[2:4], Seq("CG"))
        self.assertNotEqual(seq2[3:5], Seq("CG"))
        self.assertNotEqual(seq1[2:4], MutableSeq("CG"))
        self.assertNotEqual(seq2[3:5], MutableSeq("CG"))
        self.assertNotEqual(seq2, seq1)
        self.assertNotEqual("AT", seq1)
        self.assertNotEqual("AT", seq2)
        self.assertNotEqual(b"AT", seq1)
        self.assertNotEqual(b"AT", seq2)
        self.assertNotEqual(Seq("AT"), seq1)
        self.assertNotEqual(Seq("AT"), seq2)
        self.assertNotEqual(MutableSeq("AT"), seq1)
        self.assertNotEqual(MutableSeq("AT"), seq2)
        self.assertNotEqual("CG", seq1[2:4])
        self.assertNotEqual("CG", seq2[3:5])
        self.assertNotEqual(b"CG", seq1[2:4])
        self.assertNotEqual(b"CG", seq2[3:5])
        self.assertNotEqual(Seq("CG"), seq1[2:4])
        self.assertNotEqual(Seq("CG"), seq2[3:5])
        self.assertNotEqual(MutableSeq("CG"), seq1[2:4])
        self.assertNotEqual(MutableSeq("CG"), seq2[3:5])
        with self.assertRaises(UndefinedSequenceError):
            seq1 != Seq(None, len(seq1))
        with self.assertRaises(UndefinedSequenceError):
            seq2 != Seq(None, len(seq2))
        with self.assertRaises(UndefinedSequenceError):
            seq1 != Seq(None, 10)
        with self.assertRaises(UndefinedSequenceError):
            seq2 != Seq(None, 10)
        with self.assertRaises(UndefinedSequenceError):
            Seq(None, len(seq1)) != seq1
        with self.assertRaises(UndefinedSequenceError):
            Seq(None, len(seq2)) != seq2
        with self.assertRaises(UndefinedSequenceError):
            Seq(None, 10) != seq1
        with self.assertRaises(UndefinedSequenceError):
            Seq(None, 10) != seq2

    def test_lt(self):
        seq1 = self.seq1
        seq2 = self.seq2
        self.assertLess(seq1, seq2)
        self.assertLess("AA", seq1)
        self.assertLess(seq1, "TT")
        self.assertLess("AA", seq2)
        self.assertLess(seq2, "TTT")
        self.assertLess(b"AA", seq1)
        self.assertLess(seq1, b"TT")
        self.assertLess(b"AA", seq2)
        self.assertLess(seq2, b"TTT")
        self.assertLess(Seq("AA"), seq1)
        self.assertLess(seq1, Seq("TT"))
        self.assertLess(Seq("AA"), seq2)
        self.assertLess(seq2, Seq("TTT"))
        self.assertLess(MutableSeq("AA"), seq1)
        self.assertLess(seq1, MutableSeq("TT"))
        self.assertLess(MutableSeq("AA"), seq2)
        self.assertLess(seq2, MutableSeq("TTT"))
        with self.assertRaises(UndefinedSequenceError):
            seq1 < Seq(None, len(seq1))
        with self.assertRaises(UndefinedSequenceError):
            seq2 < Seq(None, len(seq2))
        with self.assertRaises(UndefinedSequenceError):
            seq1 < Seq(None, 10)
        with self.assertRaises(UndefinedSequenceError):
            seq2 < Seq(None, 10)
        self.assertLess("AA", seq1[2:4])
        self.assertLess("AA", seq2[3:5])
        self.assertLess(b"AA", seq1[2:4])
        self.assertLess(b"AA", seq2[3:5])
        self.assertLess(seq2[3:5], seq1[2:6])
        self.assertLess(Seq("AA"), seq1[2:4])
        self.assertLess(Seq("AA"), seq2[3:5])
        self.assertLess(MutableSeq("AA"), seq1[2:4])
        self.assertLess(MutableSeq("AA"), seq2[3:5])
        self.assertLess(seq1[2:4], "TT")
        self.assertLess(seq2[3:5], "TT")
        self.assertLess(seq1[2:4], b"TT")
        self.assertLess(seq2[3:5], b"TT")
        self.assertLess(seq1[2:4], Seq("TT"))
        self.assertLess(seq2[3:5], Seq("TT"))
        self.assertLess(seq1[2:4], MutableSeq("TT"))
        self.assertLess(seq2[3:5], MutableSeq("TT"))

    def test_le(self):
        seq1 = self.seq1
        seq2 = self.seq2
        self.assertLessEqual(seq1, seq2)
        self.assertLessEqual(seq1, "TT")
        self.assertLessEqual("TT", seq2)
        self.assertLessEqual(seq1, b"TT")
        self.assertLessEqual("TT", seq2)
        self.assertLessEqual(seq1, Seq("TT"))
        self.assertLessEqual("TT", seq2)
        self.assertLessEqual(seq1, MutableSeq("TT"))
        self.assertLessEqual(MutableSeq("TT"), seq2)
        self.assertLessEqual("AA", seq1)
        self.assertLessEqual("AA", seq2)
        self.assertLessEqual(b"AA", seq1)
        self.assertLessEqual(b"AA", seq2)
        self.assertLessEqual(Seq("AA"), seq1)
        self.assertLessEqual(Seq("AA"), seq2)
        self.assertLessEqual(MutableSeq("AA"), seq1)
        self.assertLessEqual(MutableSeq("AA"), seq2)
        self.assertLessEqual("GC", seq1)
        self.assertLessEqual("GC", seq2)
        self.assertLessEqual(b"GC", seq1)
        self.assertLessEqual(b"GC", seq2)
        self.assertLessEqual(Seq("GC"), seq1)
        self.assertLessEqual(Seq("GC"), seq2)
        self.assertLessEqual(MutableSeq("GC"), seq1)
        self.assertLessEqual(MutableSeq("GC"), seq2)
        with self.assertRaises(UndefinedSequenceError):
            seq1 <= Seq(None, len(seq1))
        with self.assertRaises(UndefinedSequenceError):
            seq2 <= Seq(None, len(seq2))
        with self.assertRaises(UndefinedSequenceError):
            seq1 <= Seq(None, 10)
        with self.assertRaises(UndefinedSequenceError):
            seq2 <= Seq(None, 10)
        self.assertLessEqual("AA", seq1[2:4])
        self.assertLessEqual("AA", seq2[3:5])
        self.assertLessEqual(b"AA", seq1[2:4])
        self.assertLessEqual(b"AA", seq2[3:5])
        self.assertLessEqual(seq1[2:4], seq2[3:5])
        self.assertLessEqual(Seq("AA"), seq1[2:4])
        self.assertLessEqual(Seq("AA"), seq2[3:5])
        self.assertLessEqual(MutableSeq("AA"), seq1[2:4])
        self.assertLessEqual(MutableSeq("AA"), seq2[3:5])
        self.assertLessEqual(seq1[2:4], "TT")
        self.assertLessEqual(seq2[3:5], "TT")
        self.assertLessEqual(seq1[2:4], b"TT")
        self.assertLessEqual(seq2[3:5], b"TT")
        self.assertLessEqual(seq1[2:4], seq2[3:5])
        self.assertLessEqual(seq1[2:4], Seq("TT"))
        self.assertLessEqual(seq2[3:5], Seq("TT"))
        self.assertLessEqual(seq1[2:4], MutableSeq("TT"))
        self.assertLessEqual(seq2[3:5], MutableSeq("TT"))

    def test_gt(self):
        seq1 = self.seq1
        seq2 = self.seq2
        self.assertGreater(seq2, seq1)
        self.assertGreater("TT", seq1)
        self.assertGreater(seq2, "TT")
        self.assertGreater(b"TT", seq1)
        self.assertGreater(seq2, b"TT")
        self.assertGreater(Seq("TT"), seq1)
        self.assertGreater(seq2, Seq("TT"))
        self.assertGreater(MutableSeq("TT"), seq1)
        self.assertGreater(seq2, MutableSeq("TT"))
        self.assertGreater(seq1, "AA")
        self.assertGreater(seq2, "AA")
        self.assertGreater(seq1, b"AA")
        self.assertGreater(seq2, b"AA")
        self.assertGreater(seq1, Seq("AA"))
        self.assertGreater(seq2, Seq("AA"))
        self.assertGreater(seq1, MutableSeq("AA"))
        self.assertGreater(seq2, MutableSeq("AA"))
        self.assertGreater(seq1, "GC")
        self.assertGreater(seq2, "GC")
        self.assertGreater(seq1, b"GC")
        self.assertGreater(seq2, b"GC")
        self.assertGreater(seq1, Seq("GC"))
        self.assertGreater(seq2, Seq("GC"))
        self.assertGreater(seq1, MutableSeq("GC"))
        self.assertGreater(seq2, MutableSeq("GC"))
        with self.assertRaises(UndefinedSequenceError):
            seq1 > Seq(None, len(seq1))
        with self.assertRaises(UndefinedSequenceError):
            seq2 > Seq(None, len(seq2))
        with self.assertRaises(UndefinedSequenceError):
            seq1 > Seq(None, 10)
        with self.assertRaises(UndefinedSequenceError):
            seq2 > Seq(None, 10)
        self.assertGreater(seq1[2:4], "AA")
        self.assertGreater(seq2[3:5], "AA")
        self.assertGreater(seq1[2:4], b"AA")
        self.assertGreater(seq2[3:5], b"AA")
        self.assertGreater(seq1[2:6], seq2[3:5])
        self.assertGreater(seq1[2:4], Seq("AA"))
        self.assertGreater(seq2[3:5], Seq("AA"))
        self.assertGreater(seq1[2:4], MutableSeq("AA"))
        self.assertGreater(seq2[3:5], MutableSeq("AA"))
        self.assertGreater("TT", seq1[2:4])
        self.assertGreater("TT", seq2[3:5])
        self.assertGreater(b"TT", seq1[2:4])
        self.assertGreater(b"TT", seq2[3:5])
        self.assertGreater(Seq("TT"), seq1[2:4])
        self.assertGreater(Seq("TT"), seq2[3:5])
        self.assertGreater(MutableSeq("TT"), seq1[2:4])
        self.assertGreater(MutableSeq("TT"), seq2[3:5])

    def test_ge(self):
        seq1 = self.seq1
        seq2 = self.seq2
        self.assertGreaterEqual(seq2, seq1)
        self.assertGreaterEqual("TT", seq1)
        self.assertGreaterEqual(seq2, "TT")
        self.assertGreaterEqual(b"TT", seq1)
        self.assertGreaterEqual(seq2, b"TT")
        self.assertGreaterEqual(Seq("TT"), seq1)
        self.assertGreaterEqual(seq2, Seq("TT"))
        self.assertGreaterEqual(MutableSeq("TT"), seq1)
        self.assertGreaterEqual(seq2, MutableSeq("TT"))
        self.assertGreaterEqual(seq1, "AA")
        self.assertGreaterEqual(seq2, "AA")
        self.assertGreaterEqual(seq1, b"AA")
        self.assertGreaterEqual(seq2, b"AA")
        self.assertGreaterEqual(seq1, Seq("AA"))
        self.assertGreaterEqual(seq2, Seq("AA"))
        self.assertGreaterEqual(seq1, MutableSeq("AA"))
        self.assertGreaterEqual(seq2, MutableSeq("AA"))
        self.assertGreaterEqual(seq1, "GC")
        self.assertGreaterEqual(seq2, "GC")
        self.assertGreaterEqual(seq1, b"GC")
        self.assertGreaterEqual(seq2, b"GC")
        self.assertGreaterEqual(seq1, Seq("GC"))
        self.assertGreaterEqual(seq2, Seq("GC"))
        self.assertGreaterEqual(seq1, MutableSeq("GC"))
        self.assertGreaterEqual(seq2, MutableSeq("GC"))
        with self.assertRaises(UndefinedSequenceError):
            seq1 >= Seq(None, len(seq1))
        with self.assertRaises(UndefinedSequenceError):
            seq2 >= Seq(None, len(seq2))
        with self.assertRaises(UndefinedSequenceError):
            seq1 >= Seq(None, 10)
        with self.assertRaises(UndefinedSequenceError):
            seq2 >= Seq(None, 10)
        self.assertGreaterEqual(seq1[2:4], "AA")
        self.assertGreaterEqual(seq2[3:5], "AA")
        self.assertGreaterEqual(seq1[2:4], b"AA")
        self.assertGreaterEqual(seq2[3:5], b"AA")
        self.assertGreaterEqual(seq1[2:4], seq2[3:5])
        self.assertGreaterEqual(seq1[2:4], Seq("AA"))
        self.assertGreaterEqual(seq2[3:5], Seq("AA"))
        self.assertGreaterEqual(seq1[2:4], MutableSeq("AA"))
        self.assertGreaterEqual(seq2[3:5], MutableSeq("AA"))
        self.assertGreaterEqual("TT", seq1[2:4])
        self.assertGreaterEqual("TT", seq2[3:5])
        self.assertGreaterEqual(b"TT", seq1[2:4])
        self.assertGreaterEqual(b"TT", seq2[3:5])
        self.assertGreaterEqual(seq1[2:4], seq2[3:5])
        self.assertGreaterEqual(Seq("TT"), seq1[2:4])
        self.assertGreaterEqual(Seq("TT"), seq2[3:5])
        self.assertGreaterEqual(MutableSeq("TT"), seq1[2:4])
        self.assertGreaterEqual(MutableSeq("TT"), seq2[3:5])


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
