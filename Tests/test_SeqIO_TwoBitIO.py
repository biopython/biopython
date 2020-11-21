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
        records = SeqIO.parse(path, 'fasta')
        self.records = list(records)

    def test_littleendian(self):
        path = "TwoBit/sequence.littleendian.2bit"
        handle = open(path)
        records = TwoBitIO.TwoBitIterator(handle)
        self.assertEqual(len(sequences), 5)
        self.assertFalse(records.isByteSwapped)
        for record1, record2 in zip(self.records, records):
            self.assertEqual(record1.seq, record2.seq)
            self.assertEqual(record1.id, record2.id)
        handle.close()

    def test_bigendian(self):
        path = "TwoBit/sequence.bigendian.2bit"
        handle = open(path)
        records = TwoBitIO.TwoBitIterator(handle)
        self.assertEqual(len(sequences), 5)
        self.assertFalse(records.isByteSwapped)
        for record1, record2 in zip(self.records, records):
            self.assertEqual(record1.seq, record2.seq)
            self.assertEqual(record1.id, record2.id)
        handle.close()

if False:
  for length in range(1,21):
    records = []
    for i in range(10):
        nucleotides = ['ACGTNacgtn'[random.randint(0,9)] for i in range(length)]
        sequence = ''.join(nucleotides)
        seq = Seq(sequence)
        record = SeqRecord(seq, id='name_%d' % i)
        records.append(record)
    handle = open("test.fa", 'w')
    SeqIO.write(records, handle, 'fasta')
    handle.close()
    os.system("faToTwoBit test.fa test.2bit")
    handle = open("test.2bit")
    sequences = TwoBitIO.TwoBitIterator(handle)
    for sequence, record in zip(sequences, records):
        seq1 = sequence
        seq2 = str(record.seq)
        assert str(seq1) == str(seq2)
        for start in range(length):
            for end in range(start+1, length+1):
                print("Testing sequence length %d start %d end %d" % (length, start, end))
                for i in range(10):
                    assert str(seq1[start:end]) == str(seq2[start:end])
                    for step in range(1, end-start+1):
                        assert str(seq1[start:end:step]) == str(seq2[start:end:step])
                        assert str(seq1[end:start:-step]) == str(seq2[end:start:-step])



if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
