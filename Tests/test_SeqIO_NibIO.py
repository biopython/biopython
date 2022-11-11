"""Tests for SeqIO NibIO module."""
import unittest

from io import BytesIO

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class TestNibReaderWriter(unittest.TestCase):
    def test_read_even(self):
        with open("Nib/test_even.fa") as handle:
            record = SeqIO.read(handle, "fasta")
        sequence = record.seq
        with open("Nib/test_even_bigendian.nib", "rb") as handle:
            record = SeqIO.read(handle, "nib")
        self.assertEqual(sequence, record.seq)
        with open("Nib/test_even_littleendian.nib", "rb") as handle:
            record = SeqIO.read(handle, "nib")
        self.assertEqual(sequence, record.seq)

    def test_read_odd(self):
        with open("Nib/test_odd.fa") as handle:
            record = SeqIO.read(handle, "fasta")
        sequence = record.seq
        with open("Nib/test_odd_bigendian.nib", "rb") as handle:
            record = SeqIO.read(handle, "nib")
        self.assertEqual(sequence, record.seq)
        with open("Nib/test_odd_littleendian.nib", "rb") as handle:
            record = SeqIO.read(handle, "nib")
        self.assertEqual(sequence, record.seq)

    def test_write_even(self):
        with open("Nib/test_even.fa") as handle:
            record = SeqIO.read(handle, "fasta")
        sequence = record.seq
        handle = BytesIO()
        n = SeqIO.write(record, handle, "nib")
        self.assertEqual(n, 1)
        handle.flush()
        handle.seek(0)
        record = SeqIO.read(handle, "nib")
        handle.close()
        self.assertEqual(sequence, record.seq)

    def test_write_odd(self):
        with open("Nib/test_odd.fa") as handle:
            record = SeqIO.read(handle, "fasta")
        sequence = record.seq
        handle = BytesIO()
        n = SeqIO.write(record, handle, "nib")
        self.assertEqual(n, 1)
        handle.flush()
        handle.seek(0)
        record = SeqIO.read(handle, "nib")
        handle.close()
        self.assertEqual(sequence, record.seq)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
