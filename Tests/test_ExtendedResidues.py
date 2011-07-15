# Copyright 2011 by Mikael Trellet.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for the Bio.PDB module."""
import os
import unittest
import warnings
from StringIO import StringIO

from Bio.PDB import Residue
from Bio.PDB import PDBParser

class ExtendedResidues(unittest.TestCase):
    " Test ExtendedResidues Module from Residue.py"

    def setUp(self):
        "Setup"

        P = PDBParser()
        structure = P.get_structure('test', 'PDB/1MOT.pdb')[0]
        self.structure = structure

    def test_ExtendResidue(self):
        "Test ExtendedResidue conversion from Residue"

        s=self.structure
        res = list(s	.get_residues())[0]
        res.extend()

        self.assertEqual(res.__repr__(), '<ExtendedResidue ALA het=  resseq=249 icode= >')
        self.assertEqual(res.id, (' ', 249, ' '))
        self.assertEqual(res.parent.id, 'A')
        self.assertEqual(res.extended, 1)
        self.assertEqual(res.hydrophobicity, 0.62)
        self.assertEqual(res.charge, 0)
        self.assertEqual("%2.3f" %res.mass, "%2.3f" %(73.094))

if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
