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

import os
import unittest
import warnings

from Bio import Align
from Bio import AlignIO
from Bio.Align import PairwiseAligner
from Bio.Align import MultipleSeqAlignment
from Bio.Data import PDBData
from Bio.PDB import PDBParser
from Bio.PDB import StructureAlignment
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.PDB.Selection import unfold_entities
from Bio.PDB.Polypeptide import is_aa
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
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

    def test_msa_vs_alignment_objects(self):
        """Test that MultipleSeqAlignment and Alignment objects produce identical results.

        This test verifies that when given the same sequence alignment data,
        StructureAlignment produces identical mappings whether the input is:
        1. A Bio.Align.MultipleSeqAlignment object (from AlignIO.read)
        2. A Bio.Align.Alignment object (from Align.read)
        """
        p = PDBParser(QUIET=1)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
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

    def test_custom_vs_automatic_alignment(self):
        """Test that custom PairwiseAligner alignment vs automatic alignment work identically.

        This test verifies that providing a custom alignment created with PairwiseAligner
        produces the same results as letting StructureAlignment generate the alignment
        automatically.
        """

        p = PDBParser(QUIET=1)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            s1 = p.get_structure("1", "PDB/2XHE.pdb")
            s2 = p.get_structure("2", "PDB/1A8O.pdb")
        m1 = s1[0]
        m2 = s2[0]

        # Extract sequences from models (same logic as in StructureAlignment)
        def extract_sequence(model, seq_id):
            residues = unfold_entities(model, "R")
            sequence = ""
            for residue in residues:
                if is_aa(residue):
                    resname = residue.get_resname()
                    if resname in PDBData.protein_letters_3to1_extended:
                        aa_code = PDBData.protein_letters_3to1_extended[resname]
                        sequence += aa_code
                    else:
                        sequence += "X"
            return SeqRecord(Seq(sequence), id=seq_id)

        seq1_record = extract_sequence(m1, "structure1")
        seq2_record = extract_sequence(m2, "structure2")

        # Create custom alignment using PairwiseAligner
        aligner = PairwiseAligner("blastp")
        alignments = aligner.align(seq1_record.seq, seq2_record.seq)
        best_alignment = alignments[0]

        # Convert to MultipleSeqAlignment
        aligned_seq1_str = str(best_alignment[0])
        aligned_seq2_str = str(best_alignment[1])
        aligned_seq1 = SeqRecord(Seq(aligned_seq1_str), id="structure1")
        aligned_seq2 = SeqRecord(Seq(aligned_seq2_str), id="structure2")
        custom_alignment = MultipleSeqAlignment([aligned_seq1, aligned_seq2])

        # Create StructureAlignment with custom alignment
        al_custom = StructureAlignment(custom_alignment, m1, m2)

        # Create StructureAlignment with automatic alignment
        al_auto = StructureAlignment(m1=m1, m2=m2)

        # Results should be identical since they use the same alignment algorithm
        self.assertEqual(len(al_custom.duos), len(al_auto.duos))
        self.assertEqual(len(al_custom.map12), len(al_auto.map12))
        self.assertEqual(len(al_custom.map21), len(al_auto.map21))

        # Compare mappings
        for residue in al_custom.map12:
            self.assertIn(residue, al_auto.map12)
            self.assertEqual(al_custom.map12[residue], al_auto.map12[residue])

        for residue in al_custom.map21:
            self.assertIn(residue, al_auto.map21)
            self.assertEqual(al_custom.map21[residue], al_auto.map21[residue])

        # Compare duos
        for i, (duo_custom, duo_auto) in enumerate(zip(al_custom.duos, al_auto.duos)):
            self.assertEqual(
                duo_custom,
                duo_auto,
                f"Duo mismatch at position {i}: {duo_custom} vs {duo_auto}",
            )

    def test_automatic_alignment_generation(self):
        """Test that automatic alignment generation works correctly.

        This test verifies that StructureAlignment can generate alignments
        automatically when no alignment is provided.
        """
        p = PDBParser(QUIET=1)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            s1 = p.get_structure("1", "PDB/2XHE.pdb")
            s2 = p.get_structure("2", "PDB/1A8O.pdb")
        m1 = s1[0]
        m2 = s2[0]

        # Create StructureAlignment with automatic alignment generation
        al_auto = StructureAlignment(m1=m1, m2=m2)

        # Verify that mappings were created
        self.assertIsNotNone(al_auto.map12)
        self.assertIsNotNone(al_auto.map21)
        self.assertIsNotNone(al_auto.duos)
        self.assertGreater(len(al_auto.duos), 0)

        # Verify that some residues were mapped
        self.assertGreater(len(al_auto.map12), 0)
        self.assertGreater(len(al_auto.map21), 0)


if __name__ == "__main__":
    os.chdir("Tests")
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
