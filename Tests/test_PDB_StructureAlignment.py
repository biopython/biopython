# Copyright 2017 by Francesco Gastaldello. All rights reserved.
#
# Converted by Francesco Gastaldello from an older unit test copyright 2004
# by Thomas Hamelryck.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Unit tests for the Bio.PDB.StructureAlignment module."""

import unittest

from Bio import Align
from Bio import AlignIO
from Bio.PDB import PDBParser
from Bio.PDB import StructureAlignment


class StructureAlignTests(unittest.TestCase):
    """Test module StructureAlignment."""

    def test_StructAlign(self):
        """Tests on module to align two proteins according to a FASTA file."""
        p = PDBParser(QUIET=1)

        al_file = "PDB/alignment_file.fa"

        # Using Bio.AlignIO, which returns a MultipleSeqAlignment object:
        with open(al_file) as handle:
            records = AlignIO.read(handle, "fasta")

        # Using Bio.Align, which returns an Alignment object:
        with open(al_file) as handle:
            alignment = Align.read(handle, "fasta")

        s1 = p.get_structure("1", "PDB/2XHE.pdb")
        s2 = p.get_structure("2", "PDB/1A8O.pdb")
        m1 = s1[0]
        m2 = s2[0]

        for argument in (records, alignment):
            al = StructureAlignment(argument, m1, m2)
            self.assertNotEqual(al.map12, al.map21)
            self.assertTrue(len(al.map12), 566)
            self.assertTrue(len(al.map21), 70)
            chain1_A = m1["A"]
            chain2_A = m2["A"]
            self.assertEqual(chain1_A[202].get_resname(), "ILE")
            self.assertEqual(chain2_A[202].get_resname(), "LEU")
            self.assertNotEqual(
                chain1_A[291].get_resname(), chain2_A[181].get_resname()
            )

    def test_StructAlign_auto_generate(self):
        """Test auto-generation when fasta_align is None."""

        p = PDBParser(QUIET=1)

        s1 = p.get_structure("1", "PDB/2XHE.pdb")
        s2 = p.get_structure("2", "PDB/1A8O.pdb")
        m1 = s1[0]
        m2 = s2[0]

        al = StructureAlignment(fasta_align=None, m1=m1, m2=m2)

        # Basic sanity checks
        self.assertIsNotNone(al.map12)
        self.assertIsNotNone(al.map21)
        self.assertIsNotNone(al.duos)

        self.assertGreater(len(al.duos), 0)

        # Test that the iterator works
        duo_count = 0
        gap_positions = 0
        for duo in al.get_iterator():
            duo_count += 1
            self.assertIsInstance(duo, tuple)
            self.assertEqual(len(duo), 2)
            # Count positions where either sequence has a gap (None)
            if duo[0] is None or duo[1] is None:
                gap_positions += 1

        self.assertEqual(duo_count, len(al.duos))

        non_gap_positions = duo_count - gap_positions
        self.assertGreater(non_gap_positions, 0)

    def test_StructAlign_model_validation(self):
        """Test that ValueError is raised when models are None."""
        p = PDBParser(QUIET=1)

        s1 = p.get_structure("1", "PDB/2XHE.pdb")
        m1 = s1[0]

        # Test with m1=None
        with self.assertRaises(ValueError) as context:
            StructureAlignment(fasta_align=None, m1=None, m2=m1)
        self.assertIn(
            "Both m1 and m2 models must be provided and cannot be None",
            str(context.exception),
        )

        # Test with m2=None
        with self.assertRaises(ValueError) as context:
            StructureAlignment(fasta_align=None, m1=m1, m2=None)
        self.assertIn(
            "Both m1 and m2 models must be provided and cannot be None",
            str(context.exception),
        )

        # Test with both None
        with self.assertRaises(ValueError) as context:
            StructureAlignment(fasta_align=None, m1=None, m2=None)
        self.assertIn(
            "Both m1 and m2 models must be provided and cannot be None",
            str(context.exception),
        )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
