"""Tests for SeqIO TwoBitIO module."""
import os
import random
import unittest

from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.Seq import Seq
from Bio.Seq import UndefinedSequenceError
from Bio.SeqRecord import SeqRecord


class Parsing(unittest.TestCase):
    """Test parsing 2bit files."""

    def setUp(self):
        path = "TwoBit/sequence.fa"
        records = SeqIO.parse(path, "fasta")
        self.records = list(records)

    def test_littleendian(self, step=5):
        path = "TwoBit/sequence.littleendian.2bit"
        with open(path, "rb") as stream:
            records = SeqIO.parse(stream, "twobit")
            self.assertEqual(records.byteorder, "little")
            self.assertEqual(len(self.records), len(records))
            for record1, record2 in zip(self.records, records):
                self.assertEqual(record1.id, record2.id)
                seq1 = record1.seq
                seq2 = record2.seq
                self.assertEqual(seq1, seq2)
                n = len(seq1)
                for i in range(0, n, step):
                    for j in range(i, n, step):
                        self.assertEqual(seq1[i:j], seq2[i:j])
                        self.assertEqual(repr(seq1[i:j]), repr(seq2[i:j]))

    def test_bigendian(self, step=5):
        path = "TwoBit/sequence.bigendian.2bit"
        with open(path, "rb") as stream:
            records = SeqIO.parse(stream, "twobit")
            self.assertEqual(len(records), 6)
            self.assertEqual(records.byteorder, "big")
            for record1, record2 in zip(self.records, records):
                self.assertEqual(record1.id, record2.id)
                seq1 = record1.seq
                seq2 = record2.seq
                self.assertEqual(seq1, seq2)
                n = len(seq1)
                for i in range(0, n, step):
                    for j in range(i, n, step):
                        self.assertEqual(seq1[i:j], seq2[i:j])
                        self.assertEqual(repr(seq1[i:j]), repr(seq2[i:j]))

    def test_sequence_long(self):
        path = "TwoBit/sequence.long.2bit"
        with open(path, "rb") as stream:
            with self.assertRaises(ValueError) as cm:
                SeqIO.parse(stream, "twobit")
            self.assertEqual(
                str(cm.exception),
                "version-1 twoBit files with 64-bit offsets for index are currently not supported",
            )


class TestComparisons(unittest.TestCase):
    """Test comparisons of sequences read from 2bit files to Seq and other objects."""

    def setUp(self):
        path = "TwoBit/sequence.bigendian.2bit"
        self.stream = open(path, "rb")
        records = SeqIO.parse(self.stream, "twobit")
        record1 = next(records)
        record2 = next(records)
        self.seq1a = record1.seq
        self.seq2a = record2.seq
        path = "TwoBit/sequence.fa"
        records = SeqIO.parse(path, "fasta")
        record1 = next(records)
        record2 = next(records)
        self.seq1b = record1.seq
        self.seq2b = record2.seq

    def tearDown(self):
        self.stream.close()

    def test_eq(self):
        seq1a = self.seq1a
        seq2a = self.seq2a
        seq1b = self.seq1b
        seq2b = self.seq2b
        self.assertEqual(seq1a, seq1b)
        self.assertEqual(seq2a, seq2b)
        self.assertEqual(seq1a, seq1a)
        self.assertEqual(seq2a, seq2a)
        with self.assertRaises(UndefinedSequenceError):
            seq1a == Seq(None, len(seq1a))
        with self.assertRaises(UndefinedSequenceError):
            seq2a == Seq(None, len(seq2a))
        with self.assertRaises(UndefinedSequenceError):
            seq1a == Seq(None, 10)
        with self.assertRaises(UndefinedSequenceError):
            seq2a == Seq(None, 10)
        with self.assertRaises(UndefinedSequenceError):
            Seq(None, len(seq1a)) == seq1a
        with self.assertRaises(UndefinedSequenceError):
            Seq(None, len(seq2a)) == seq2a
        with self.assertRaises(UndefinedSequenceError):
            Seq(None, 10) == seq1a
        with self.assertRaises(UndefinedSequenceError):
            Seq(None, 10) == seq2a

    def test_ne(self):
        seq1a = self.seq1a
        seq2a = self.seq2a
        seq1b = self.seq1b
        seq2b = self.seq2b
        self.assertNotEqual(seq1a, seq2a)
        self.assertNotEqual(seq1a, seq2b)
        with self.assertRaises(UndefinedSequenceError):
            seq1a != Seq(None, len(seq1a))
        with self.assertRaises(UndefinedSequenceError):
            seq2a != Seq(None, len(seq2a))
        with self.assertRaises(UndefinedSequenceError):
            seq1a != Seq(None, 10)
        with self.assertRaises(UndefinedSequenceError):
            seq2a != Seq(None, 10)
        with self.assertRaises(UndefinedSequenceError):
            Seq(None, len(seq1a)) != seq1a
        with self.assertRaises(UndefinedSequenceError):
            Seq(None, len(seq2a)) != seq2a
        with self.assertRaises(UndefinedSequenceError):
            Seq(None, 10) != seq1a
        with self.assertRaises(UndefinedSequenceError):
            Seq(None, 10) != seq2a

    def test_lt(self):
        seq1 = self.seq1a
        seq2 = self.seq2a
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

    def test_le(self):
        seq1 = self.seq1a
        seq2 = self.seq2a
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

    def test_gt(self):
        seq1 = self.seq1a
        seq2 = self.seq2a
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

    def test_ge(self):
        seq1 = self.seq1a
        seq2 = self.seq2a
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


class TestBaseClassMethods(unittest.TestCase):
    """Test if methods from the base class are called correctly."""

    def setUp(self):
        path = "TwoBit/sequence.bigendian.2bit"
        self.stream = open(path, "rb")
        records = SeqIO.parse(self.stream, "twobit")
        self.record1_twobit = next(records)
        self.seq1_twobit = self.record1_twobit.seq
        self.record2_twobit = next(records)
        self.seq2_twobit = self.record2_twobit.seq
        path = "TwoBit/sequence.fa"
        records = SeqIO.parse(path, "fasta")
        self.record1_fasta = next(records)
        self.seq1_fasta = self.record1_fasta.seq
        self.record2_fasta = next(records)
        self.seq2_fasta = self.record2_fasta.seq

    def tearDown(self):
        self.stream.close()

    def test_getitem(self):
        self.assertEqual(self.seq1_twobit, self.seq1_fasta)
        self.assertEqual(self.seq2_twobit, self.seq2_fasta)
        self.assertEqual(self.seq1_twobit[:], self.seq1_fasta[:])
        self.assertEqual(self.seq2_twobit[:], self.seq2_fasta[:])
        self.assertEqual(self.seq1_twobit[30], self.seq1_fasta[30])
        self.assertEqual(self.seq2_twobit[30], self.seq2_fasta[30])
        self.assertEqual(self.seq1_twobit[-30], self.seq1_fasta[-30])
        self.assertEqual(self.seq2_twobit[-30], self.seq2_fasta[-30])
        self.assertEqual(self.record1_twobit.seq, self.record1_fasta.seq)
        self.assertEqual(self.record2_twobit.seq, self.record2_fasta.seq)
        self.assertEqual(self.record1_twobit[:].seq, self.record1_fasta[:].seq)
        self.assertEqual(self.record2_twobit[:].seq, self.record2_fasta[:].seq)
        self.assertEqual(self.record1_twobit[30], self.record1_fasta[30])
        self.assertEqual(self.record2_twobit[30], self.record2_fasta[30])
        self.assertEqual(self.record1_twobit[-30], self.record1_fasta[-30])
        self.assertEqual(self.record2_twobit[-30], self.record2_fasta[-30])

    def test_bytes(self):
        b = bytes(self.seq1_twobit)
        self.assertIsInstance(b, bytes)
        self.assertEqual(len(b), 480)
        self.assertEqual(b, bytes(self.seq1_fasta))
        b = bytes(self.seq1_twobit[:10])
        self.assertEqual(len(b), 10)
        self.assertIsInstance(b, bytes)
        self.assertEqual(b, b"GTATACCCCT")

    def test_hash(self):
        self.assertEqual(hash(self.seq1_twobit), hash(self.seq1_fasta))

    def test_add(self):
        self.assertIsInstance(self.seq1_twobit + "ABCD", Seq)
        self.assertEqual(self.seq1_twobit + "ABCD", self.seq1_fasta + "ABCD")
        self.assertIsInstance(self.record1_twobit + "ABCD", SeqRecord)
        record1_twobit = self.record1_twobit + "ABCD"
        record1_fasta = self.record1_fasta + "ABCD"
        self.assertEqual(self.record1_twobit.seq, self.record1_fasta.seq)

    def test_radd(self):
        self.assertIsInstance("ABCD" + self.seq1_twobit, Seq)
        self.assertEqual("ABCD" + self.seq1_twobit, "ABCD" + self.seq1_fasta)
        self.assertIsInstance("ABCD" + self.record1_twobit, SeqRecord)
        record1_twobit = "ABCD" + self.record1_twobit
        record1_fasta = "ABCD" + self.record1_fasta
        self.assertEqual(self.record1_twobit.seq, self.record1_fasta.seq)

    def test_mul(self):
        self.assertIsInstance(2 * self.seq1_twobit, Seq)
        self.assertEqual(2 * self.seq1_twobit, 2 * self.seq1_fasta)
        self.assertIsInstance(self.seq1_twobit * 2, Seq)
        self.assertEqual(self.seq1_twobit * 2, self.seq1_fasta * 2)

    def test_contains(self):
        for seq in (
            self.seq1_twobit,
            self.seq1_fasta,
            self.record1_twobit,
            self.record1_fasta,
        ):
            self.assertIn("ACCCCT", seq)
            self.assertNotIn("ACGTACGT", seq)

    def test_repr(self):
        self.assertIsInstance(repr(self.seq1_twobit), str)
        self.assertEqual(repr(self.seq1_twobit), repr(self.seq1_fasta))

    def test_str(self):
        self.assertIsInstance(str(self.seq1_twobit), str)
        self.assertEqual(str(self.seq1_twobit), str(self.seq1_fasta))

    def test_count(self):
        self.assertEqual(self.seq1_twobit.count("CT"), self.seq1_fasta.count("CT"))
        self.assertEqual(
            self.seq1_twobit.count("CT", 75), self.seq1_fasta.count("CT", 75)
        )
        self.assertEqual(
            self.seq1_twobit.count("CT", 125, 250),
            self.seq1_fasta.count("CT", 125, 250),
        )
        self.assertEqual(
            self.record1_twobit.count("CT"), self.record1_fasta.count("CT")
        )
        self.assertEqual(
            self.record1_twobit.count("CT", 75), self.record1_fasta.count("CT", 75)
        )
        self.assertEqual(
            self.record1_twobit.count("CT", 125, 250),
            self.record1_fasta.count("CT", 125, 250),
        )

    def test_find(self):
        self.assertEqual(self.seq1_twobit.find("CT"), self.seq1_fasta.find("CT"))
        self.assertEqual(
            self.seq1_twobit.find("CT", 75), self.seq1_fasta.find("CT", 75)
        )
        self.assertEqual(
            self.seq1_twobit.find("CT", 75, 100), self.seq1_fasta.find("CT", 75, 100)
        )
        self.assertEqual(
            self.seq1_twobit.find("CT", None, 100),
            self.seq1_fasta.find("CT", None, 100),
        )

    def test_rfind(self):
        self.assertEqual(self.seq1_twobit.rfind("CT"), self.seq1_fasta.rfind("CT"))
        self.assertEqual(
            self.seq1_twobit.rfind("CT", 450), self.seq1_fasta.rfind("CT", 450)
        )
        self.assertEqual(
            self.seq1_twobit.rfind("CT", None, 100),
            self.seq1_fasta.rfind("CT", None, 100),
        )
        self.assertEqual(
            self.seq1_twobit.rfind("CT", 75, 100), self.seq1_fasta.rfind("CT", 75, 100)
        )

    def test_index(self):
        self.assertEqual(self.seq1_twobit.index("CT"), self.seq1_fasta.index("CT"))
        self.assertEqual(
            self.seq1_twobit.index("CT", 75), self.seq1_fasta.index("CT", 75)
        )
        self.assertEqual(
            self.seq1_twobit.index("CT", None, 100),
            self.seq1_fasta.index("CT", None, 100),
        )
        for seq in (self.seq1_twobit, self.seq1_fasta):
            self.assertRaises(ValueError, seq.index, "CT", 75, 100)

    def test_rindex(self):
        self.assertEqual(self.seq1_twobit.rindex("CT"), self.seq1_fasta.rindex("CT"))
        self.assertEqual(
            self.seq1_twobit.rindex("CT", None, 100),
            self.seq1_fasta.rindex("CT", None, 100),
        )
        for seq in (self.seq1_twobit, self.seq1_fasta):
            self.assertRaises(ValueError, seq.rindex, "CT", 450)
            self.assertRaises(ValueError, seq.rindex, "CT", 75, 100)

    def test_startswith(self):
        for seq in (self.seq1_twobit, self.seq1_fasta):
            self.assertTrue(seq.startswith("GTAT"))
            self.assertTrue(seq.startswith("TGGG", start=10))
            self.assertTrue(seq.startswith("TGGG", start=10, end=14))
            self.assertFalse(seq.startswith("TGGG", start=10, end=12))

    def test_endswith(self):
        for seq in (self.seq1_twobit, self.seq1_fasta):
            self.assertTrue(seq.endswith("ACCG"))
            self.assertTrue(seq.endswith("ACCG", 476))
            self.assertTrue(seq.endswith("GCAC", 472, 478))
            self.assertFalse(seq.endswith("GCAC", 476, 478))

    def test_split(self):
        self.assertEqual(self.seq1_twobit.split(), self.seq1_fasta.split())
        self.assertEqual(self.seq1_twobit.split("C"), self.seq1_fasta.split("C"))
        self.assertEqual(self.seq1_twobit.split("C", 1), self.seq1_fasta.split("C", 1))

    def test_rsplit(self):
        self.assertEqual(self.seq1_twobit.rsplit(), self.seq1_fasta.rsplit())
        self.assertEqual(self.seq1_twobit.rsplit("C"), self.seq1_fasta.rsplit("C"))
        self.assertEqual(
            self.seq1_twobit.rsplit("C", 1), self.seq1_fasta.rsplit("C", 1)
        )

    def test_strip(self):
        self.assertEqual(self.seq1_twobit.strip("G"), self.seq1_fasta.strip("G"))

    def test_lstrip(self, chars=None):
        self.assertEqual(self.seq1_twobit.lstrip("G"), self.seq1_fasta.lstrip("G"))

    def test_rstrip(self, chars=None):
        self.assertEqual(self.seq1_twobit.rstrip("G"), self.seq1_fasta.rstrip("G"))

    def test_upper(self):
        seq1_twobit_upper = self.seq1_twobit.upper()
        seq1_fasta_upper = self.seq1_fasta.upper()
        self.assertEqual(seq1_twobit_upper, seq1_fasta_upper)
        self.assertEqual(seq1_twobit_upper[140:210], seq1_fasta_upper[140:210])
        seq2_twobit_upper = self.seq2_twobit.upper()
        seq2_fasta_upper = self.seq2_fasta.upper()
        self.assertEqual(seq2_twobit_upper, seq2_fasta_upper)
        self.assertEqual(seq2_twobit_upper[140:210], seq2_fasta_upper[140:210])

    def test_lower(self):
        seq1_twobit_lower = self.seq1_twobit.lower()
        seq1_fasta_lower = self.seq1_fasta.lower()
        self.assertEqual(seq1_twobit_lower, seq1_fasta_lower)
        self.assertEqual(seq1_twobit_lower[140:210], seq1_fasta_lower[140:210])
        seq2_twobit_lower = self.seq2_twobit.lower()
        seq2_fasta_lower = self.seq2_fasta.lower()
        self.assertEqual(seq2_twobit_lower, seq2_fasta_lower)
        self.assertEqual(seq2_twobit_lower[140:210], seq2_fasta_lower[140:210])

    def test_isupper(self):
        self.assertEqual(self.seq1_twobit.isupper(), self.seq1_fasta.isupper())
        self.assertEqual(self.seq2_twobit.isupper(), self.seq2_fasta.isupper())
        self.assertEqual(self.record1_twobit.isupper(), self.record1_fasta.isupper())
        self.assertEqual(self.record2_twobit.isupper(), self.record2_fasta.isupper())

    def test_islower(self):
        self.assertEqual(self.seq1_twobit.islower(), self.seq1_fasta.islower())
        self.assertEqual(self.seq2_twobit.islower(), self.seq2_fasta.islower())
        self.assertEqual(self.record1_twobit.islower(), self.record1_fasta.islower())
        self.assertEqual(self.record2_twobit.islower(), self.record2_fasta.islower())

    def test_replace(self):
        # seq.transcribe uses seq._data.replace
        self.assertEqual(self.seq1_twobit.transcribe(), self.seq1_fasta.transcribe())

    def test_translate(self):
        # seq.reverse_complement uses seq._data.translate
        self.assertEqual(
            self.seq1_twobit.reverse_complement(), self.seq1_fasta.reverse_complement()
        )
        record1_twobit = self.record1_twobit.reverse_complement()
        record1_fasta = self.record1_fasta.reverse_complement()
        self.assertEqual(record1_twobit.seq, record1_fasta.seq)

    def test_defined(self):
        self.assertTrue(self.seq1_twobit.defined)
        self.assertTrue(self.seq2_twobit.defined)
        self.assertEqual(self.seq1_twobit.defined_ranges, ((0, len(self.seq1_twobit)),))
        self.assertEqual(self.seq2_twobit.defined_ranges, ((0, len(self.seq2_twobit)),))


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
