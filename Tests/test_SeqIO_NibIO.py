"""Tests for SeqIO NibIO module."""

import unittest

from io import BytesIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class TestNibReaderWriter(unittest.TestCase):
    def test_read_even(self):
        handle = open("Nib/test_even.fa")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        sequence = str(record.seq)
        handle = open("Nib/test_even_bigendian.nib", "rb")
        record = SeqIO.read(handle, "nib")
        handle.close()
        self.assertEqual(sequence, str(record.seq))
        handle = open("Nib/test_even_littleendian.nib", "rb")
        record = SeqIO.read(handle, "nib")
        handle.close()
        self.assertEqual(sequence, str(record.seq))

    def test_read_odd(self):
        handle = open("Nib/test_odd.fa")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        sequence = str(record.seq)
        handle = open("Nib/test_odd_bigendian.nib", "rb")
        record = SeqIO.read(handle, "nib")
        handle.close()
        self.assertEqual(sequence, str(record.seq))
        handle = open("Nib/test_odd_littleendian.nib", "rb")
        record = SeqIO.read(handle, "nib")
        handle.close()
        self.assertEqual(sequence, str(record.seq))

    def test_write_even(self):
        handle = open("Nib/test_even.fa")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        sequence = str(record.seq)
        handle = BytesIO()
        n = SeqIO.write(record, handle, "nib")
        self.assertEqual(n, 1)
        handle.flush()
        handle.seek(0)
        record = SeqIO.read(handle, "nib")
        handle.close()
        self.assertEqual(sequence, str(record.seq))

    def test_write_odd(self):
        handle = open("Nib/test_odd.fa")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        sequence = str(record.seq)
        handle = BytesIO()
        n = SeqIO.write(record, handle, "nib")
        self.assertEqual(n, 1)
        handle.flush()
        handle.seek(0)
        record = SeqIO.read(handle, "nib")
        handle.close()
        self.assertEqual(sequence, str(record.seq))


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
