# Copyright 2001 by Iddo Friedberg.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for CAPS module."""

import unittest

from Bio import CAPS
from Bio.Restriction import EcoRI, AluI
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import Alignment, MultipleSeqAlignment


def createAlignment(sequences):
    """Create an Alignment object from a list of sequences."""
    return Alignment(
        [
            SeqRecord(Seq(s), id="sequence%i" % (i + 1))
            for (i, s) in enumerate(sequences)
        ]
    )


def createMultipleSeqAlignment(sequences):
    """Create a MultipleSeqAlignment object from a list of sequences."""
    return MultipleSeqAlignment(
        SeqRecord(Seq(s), id="sequence%i" % (i + 1)) for (i, s) in enumerate(sequences)
    )


class TestCAPS(unittest.TestCase):
    def test_trivial(self):
        enzymes = [EcoRI]
        alignment = ["gaattc", "gaactc"]
        align = createAlignment(alignment)
        capsmap = CAPS.CAPSMap(align, enzymes)

        self.assertEqual(len(capsmap.dcuts), 1)
        self.assertEqual(capsmap.dcuts[0].enzyme, EcoRI)
        self.assertEqual(capsmap.dcuts[0].start, 1)
        self.assertEqual(capsmap.dcuts[0].cuts_in, [0])
        self.assertEqual(capsmap.dcuts[0].blocked_in, [1])

    def test_trivial_msa(self):
        enzymes = [EcoRI]
        alignment = ["gaattc", "gaactc"]
        align = createMultipleSeqAlignment(alignment)
        capsmap = CAPS.CAPSMap(align, enzymes)

        self.assertEqual(len(capsmap.dcuts), 1)
        self.assertEqual(capsmap.dcuts[0].enzyme, EcoRI)
        self.assertEqual(capsmap.dcuts[0].start, 1)
        self.assertEqual(capsmap.dcuts[0].cuts_in, [0])
        self.assertEqual(capsmap.dcuts[0].blocked_in, [1])

    def test(self):
        alignment = [
            "AAAagaattcTAGATATACCAAACCAGAGAAAACAAATACATAATCGGAGAAATACAGAT"
            "AGAGAGCGAGAGAGATCGACGGCGAAGCTCTTTACCCGGAAACCATTGAAATCGGACGGT"
            "TTAGTGAAAATGGAGGATCAAGTagAtTTTGGGTTCCGTCCGAACGACGAGGAGCTCGTT"
            "GGTCACTATCTCCGTAACAAAATCGAAGGAAACACTAGCCGCGACGTTGAAGTAGCCATC"
            "AGCGAGGTCAACATCTGTAGCTACGATCCTTGGAACTTGCGCTGTAAGTTCCGAATTTTC",
            "AAAagaTttcTAGATATACCAAACCAGAGAAAACAAATACATAATCGGAGAAATACAGAT"
            "AGAGAGCGAGAGAGATCGACGGCGAAGCTCTTTACCCGGAAACCATTGAAATCGGACGGT"
            "TTAGTGAAAATGGAGGATCAAGTagctTTTGGGTTCCGTCCGAACGACGAGGAGCTCGTT"
            "GGTCACTATCTCCGTAACAAAATCGAAGGAAACACTAGCCGCGACGTTGAAGTAGCCATC"
            "AGCGAGGTCAACATCTGTAGCTACGATCCTTGGAACTTGCGCTGTAAGTTCCGAATTTTC",
            "AAAagaTttcTAGATATACCAAACCAGAGAAAACAAATACATAATCGGAGAAATACAGAT"
            "AGAGAGCGAGAGAGATCGACGGCGAAGCTCTTTACCCGGAAACCATTGAAATCGGACGGT"
            "TTAGTGAAAATGGAGGATCAAGTagctTTTGGGTTCCGTCCGAACGACGAGGAGCTCGTT"
            "GGTCACTATCTCCGTAACAAAATCGAAGGAAACACTAGCCGCGACGTTGAAGTAGCCATC"
            "AGCGAGGTCAACATCTGTAGCTACGATCCTTGGAACTTGCGCTGTAAGTTCCGAATTTTC",
        ]
        self.assertEqual(len(alignment), 3)
        enzymes = [EcoRI, AluI]
        align = createAlignment(alignment)
        capsmap = CAPS.CAPSMap(align, enzymes)

        self.assertEqual(len(capsmap.dcuts), 2)
        self.assertEqual(capsmap.dcuts[0].enzyme, EcoRI)
        self.assertEqual(capsmap.dcuts[0].start, 5)
        self.assertEqual(capsmap.dcuts[0].cuts_in, [0])
        self.assertEqual(capsmap.dcuts[0].blocked_in, [1, 2])
        self.assertEqual(capsmap.dcuts[1].enzyme, AluI)
        self.assertEqual(capsmap.dcuts[1].start, 144)
        self.assertEqual(capsmap.dcuts[1].cuts_in, [1, 2])
        self.assertEqual(capsmap.dcuts[1].blocked_in, [0])

    def test_msa(self):
        alignment = [
            "AAAagaattcTAGATATACCAAACCAGAGAAAACAAATACATAATCGGAGAAATACAGAT"
            "AGAGAGCGAGAGAGATCGACGGCGAAGCTCTTTACCCGGAAACCATTGAAATCGGACGGT"
            "TTAGTGAAAATGGAGGATCAAGTagAtTTTGGGTTCCGTCCGAACGACGAGGAGCTCGTT"
            "GGTCACTATCTCCGTAACAAAATCGAAGGAAACACTAGCCGCGACGTTGAAGTAGCCATC"
            "AGCGAGGTCAACATCTGTAGCTACGATCCTTGGAACTTGCGCTGTAAGTTCCGAATTTTC",
            "AAAagaTttcTAGATATACCAAACCAGAGAAAACAAATACATAATCGGAGAAATACAGAT"
            "AGAGAGCGAGAGAGATCGACGGCGAAGCTCTTTACCCGGAAACCATTGAAATCGGACGGT"
            "TTAGTGAAAATGGAGGATCAAGTagctTTTGGGTTCCGTCCGAACGACGAGGAGCTCGTT"
            "GGTCACTATCTCCGTAACAAAATCGAAGGAAACACTAGCCGCGACGTTGAAGTAGCCATC"
            "AGCGAGGTCAACATCTGTAGCTACGATCCTTGGAACTTGCGCTGTAAGTTCCGAATTTTC",
            "AAAagaTttcTAGATATACCAAACCAGAGAAAACAAATACATAATCGGAGAAATACAGAT"
            "AGAGAGCGAGAGAGATCGACGGCGAAGCTCTTTACCCGGAAACCATTGAAATCGGACGGT"
            "TTAGTGAAAATGGAGGATCAAGTagctTTTGGGTTCCGTCCGAACGACGAGGAGCTCGTT"
            "GGTCACTATCTCCGTAACAAAATCGAAGGAAACACTAGCCGCGACGTTGAAGTAGCCATC"
            "AGCGAGGTCAACATCTGTAGCTACGATCCTTGGAACTTGCGCTGTAAGTTCCGAATTTTC",
        ]
        self.assertEqual(len(alignment), 3)
        enzymes = [EcoRI, AluI]
        align = createMultipleSeqAlignment(alignment)
        capsmap = CAPS.CAPSMap(align, enzymes)

        self.assertEqual(len(capsmap.dcuts), 2)
        self.assertEqual(capsmap.dcuts[0].enzyme, EcoRI)
        self.assertEqual(capsmap.dcuts[0].start, 5)
        self.assertEqual(capsmap.dcuts[0].cuts_in, [0])
        self.assertEqual(capsmap.dcuts[0].blocked_in, [1, 2])
        self.assertEqual(capsmap.dcuts[1].enzyme, AluI)
        self.assertEqual(capsmap.dcuts[1].start, 144)
        self.assertEqual(capsmap.dcuts[1].cuts_in, [1, 2])
        self.assertEqual(capsmap.dcuts[1].blocked_in, [0])

    def testNoCAPS(self):
        alignment = ["aaaaaaaaaaaaaaaaaaaa", "aaaaaaaaaaaaaaaaaaaa"]
        enzymes = []
        align = createAlignment(alignment)
        capsmap = CAPS.CAPSMap(align, enzymes)
        self.assertEqual(capsmap.dcuts, [])

    def testNoCAPS_msa(self):
        alignment = ["aaaaaaaaaaaaaaaaaaaa", "aaaaaaaaaaaaaaaaaaaa"]
        enzymes = []
        align = createMultipleSeqAlignment(alignment)
        capsmap = CAPS.CAPSMap(align, enzymes)
        self.assertEqual(capsmap.dcuts, [])

    def test_uneven(self):
        alignment = [
            "aaaaaaaaaaaaaa",
            "aaaaaaaaaaaaaa",  # we'll change this below
            "aaaaaaaaaaaaaa",
        ]
        align = createMultipleSeqAlignment(alignment)
        align[1].seq = align[1].seq[:8]  # evil
        self.assertRaises(CAPS.AlignmentHasDifferentLengthsError, CAPS.CAPSMap, align)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
