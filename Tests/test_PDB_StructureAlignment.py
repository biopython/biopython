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

from Bio.PDB import StructureAlignment
from Bio.PDB import PDBParser
from Bio import AlignIO


class StructureAlignTests(unittest.TestCase):
    """Test module StructureAlignment."""

    def test_StructAlign(self):
        """Tests on module to align two proteins according to a FASTA alignment file."""
        al_file = "PDB/alignment_file.fa"
        pdb2 = "PDB/1A8O.pdb"
        pdb1 = "PDB/2XHE.pdb"
        with open(al_file, 'r') as handle:
            records = AlignIO.read(handle, "fasta")
        p = PDBParser()
        s1 = p.get_structure('1', pdb1)
        p = PDBParser()
        s2 = p.get_structure('2', pdb2)
        m1 = s1[0]
        m2 = s2[0]
        al = StructureAlignment(records, m1, m2)
        self.assertFalse(al.map12 == al.map21)
        self.assertTrue(len(al.map12), 566)
        self.assertTrue(len(al.map21), 70)
        chain1_A = m1["A"]
        chain2_A = m2["A"]
        self.assertEqual(chain1_A[202].get_resname(), 'ILE')
        self.assertEqual(chain2_A[202].get_resname(), 'LEU')
        self.assertEqual(chain1_A[291].get_resname(), chain2_A[180].get_resname())
        self.assertNotEqual(chain1_A[291].get_resname(), chain2_A[181].get_resname())


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
