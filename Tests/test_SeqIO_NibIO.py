"""Tests for SeqIO NibIO module."""

import unittest

from io import BytesIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class TestNibReaderWriter(unittest.TestCase):

    nucleotides = 'ACGTAAACCGTACCCGTANANCANNNNACNANNANCN'

    def test_read_bigendian(self):
        handle = open('Nib/test_bigendian.nib', 'rb')
        records = SeqIO.parse(handle, 'nib')
        record = next(records)
        handle.close()
        self.assertEqual(str(record.seq), self.nucleotides)

    def test_read_littleendian(self):
        handle = open('Nib/test_littleendian.nib', 'rb')
        records = SeqIO.parse(handle, 'nib')
        record = next(records)
        handle.close()
        self.assertEqual(str(record.seq), self.nucleotides)

    def test_write_and_read(self):
        handle = BytesIO()
        sequence = Seq(self.nucleotides)
        record = SeqRecord(sequence)
        n = SeqIO.write(record, handle, 'nib')
        self.assertEqual(n, 1)
        handle.flush()
        handle.seek(0)
        record = SeqIO.read(handle, 'nib')
        handle.close()
        sequence = record.seq
        self.assertEqual(str(sequence), self.nucleotides)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
