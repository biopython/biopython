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

    def test_StructAlign_msa_vs_alignment_objects(self):
        """Test that MultipleSeqAlignment and Alignment objects produce identical results.

        This test verifies that when given the same sequence alignment data,
        StructureAlignment produces identical mappings whether the input is:
        1. A Bio.Align.MultipleSeqAlignment object (from AlignIO.read)
        2. A Bio.Align.Alignment object (from Align.read)

        Both should map the same residues between structures in the same way.
        """
        p = PDBParser(QUIET=1)

        s1 = p.get_structure("1", "PDB/2XHE.pdb")
        s2 = p.get_structure("2", "PDB/1A8O.pdb")
        m1 = s1[0]
        m2 = s2[0]

        # Read the same alignment file that works in the first test
        al_file = "PDB/alignment_file.fa"

        # Test 1: Using Bio.AlignIO to get MultipleSeqAlignment
        with open(al_file) as handle:
            msa_records = AlignIO.read(handle, "fasta")

        # Test 2: Using Bio.Align to get Alignment object
        with open(al_file) as handle:
            alignment_obj = Align.read(handle, "fasta")

        # Create StructureAlignment with both types
        al_msa = StructureAlignment(msa_records, m1, m2)
        al_align = StructureAlignment(alignment_obj, m1, m2)

        # Results should be identical
        self.assertEqual(len(al_msa.duos), len(al_align.duos))
        self.assertEqual(len(al_msa.map12), len(al_align.map12))
        self.assertEqual(len(al_msa.map21), len(al_align.map21))

        # Compare mappings
        for residue in al_msa.map12:
            self.assertIn(residue, al_align.map12)
            self.assertEqual(al_msa.map12[residue], al_align.map12[residue])

        for residue in al_msa.map21:
            self.assertIn(residue, al_align.map21)
            self.assertEqual(al_msa.map21[residue], al_align.map21[residue])

        # Compare duos
        for i, (duo_msa, duo_align) in enumerate(zip(al_msa.duos, al_align.duos)):
            self.assertEqual(
                duo_msa,
                duo_align,
                f"Duo mismatch at position {i}: {duo_msa} vs {duo_align}",
            )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
