"""Tests for SeqIO GFA module."""

import unittest

from Bio import BiopythonWarning
from Bio import SeqIO


class TestRead(unittest.TestCase):
    def test_read_GFA1(self):
        """Test parsing valid GFA 1.x files."""
        records = list(SeqIO.parse("GFA/seq.gfa", "gfa1"))
        self.assertEqual(len(records), 8)
        self.assertEqual(records[6].id, "MTh13014")
        self.assertEqual(
            records[6].seq,
            "TTAGGTCTCCACCCCTGACTCCCCTCAGCCATAGAAGGCCCCACCCCAGTCTCAGCCCTACTCCACTCAAGCACTATAGTTGTAGCAGGAATCTTCTTACTCATCCGCTTCCACCCCCTAGCAGAAAATAGCCCACTAATCCAAACTCTAACACTATGCTTAGGCGCTATCACCACTCTGTTCGCAGCAGTCTGCGCCCTTACACAAAATGACATCAAAAAAATCGTAGCCTTCTCCACTTCAAGTCAACTAGGACTCATAATAGTTACAATCGGCATCAACCAACCACACCTAGCATTCCTGCACATCTGTACCCACGCCTTCTTCAAAGCCATACTATTTATGTGCTCCGGGTCCATCATCCACAACCTTAACAATGAACAAGATATTCGAAAAATAGGAGGACTACTCAAAACCATACCTCTCACTTCAACCTCCCTCACCATTGGCAGCCTAGCATTAGCAGGAATACCTTTCCTCACAGGTTTCTACTCCAAAGACC",
        )
        self.assertEqual(records[0].annotations["SN"], ("Z", "MT_human"))
        self.assertEqual(records[0].annotations["SO"], ("i", "0"))

        records = list(SeqIO.parse("GFA/seq_with_len.gfa", "gfa1"))
        self.assertEqual(len(records), 9)
        self.assertEqual(
            records[8].seq,
            "GAAAAATTGCCCTTGGTTTTCGCTTCGCTCAAACTCTATTGAACTTCGCTTTCGCTCAGTTCGTCGGGGCAATTTTTTGGTTAATACTT",
        )

        records = list(SeqIO.parse("GFA/fake_with_checksum.gfa", "gfa1"))
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].seq, "AAA")

        records = list(SeqIO.parse("GFA/no_seq.gfa", "gfa1"))
        self.assertEqual(len(records), 9)
        self.assertEqual(len(records[0]), 528)

    def test_read_GFA2(self):
        """Test parsing valid GFA 2.0 files."""
        records = list(SeqIO.parse("GFA/fake_gfa2.gfa", "gfa2"))
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].seq, "AAA")


class TestCorrupt(unittest.TestCase):
    def test_corrupt_gfa2(self):
        """Check a GFA 1.x file does not parse in GFA 2."""
        with self.assertRaises(ValueError):
            list(SeqIO.parse("GFA/seq.gfa", "gfa2"))

    def test_corrupt_segment_fields(self):
        """Check a GFA file with invalid fields on a segment line."""
        with self.assertRaises(ValueError):
            list(SeqIO.parse("GFA/corrupt_segment_fields.gfa", "gfa1"))

    def test_corrupt_len(self):
        """Check a GFA file with an incorrect length."""
        with self.assertWarns(BiopythonWarning):
            list(SeqIO.parse("GFA/corrupt_len.gfa", "gfa1"))

    def test_corrupt_checksum(self):
        """Check a GFA file with an incorrect checksum."""
        with self.assertWarns(BiopythonWarning):
            list(SeqIO.parse("GFA/corrupt_checksum.gfa", "gfa1"))

    def test_corrupt_tag_name(self):
        """Check a GFA file with an invalid tag name."""
        with self.assertWarns(BiopythonWarning):
            list(SeqIO.parse("GFA/corrupt_tag_name.gfa", "gfa1"))

    def test_corrupt_tag_type(self):
        """Check a GFA file with an incorrect tag type."""
        with self.assertWarns(BiopythonWarning):
            list(SeqIO.parse("GFA/corrupt_tag_type.gfa", "gfa1"))


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
