# Copyright 2017 by Francesco Gastaldello. All rights reserved.
#
# Converted by Francesco Gastaldello from an older unit test copyright 2004
# by Thomas Hamelryck.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Unit tests for the Bio.PDB.Superimposer module."""

import unittest

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.PDB.")

from Bio.PDB import Superimposer, Selection
from Bio.PDB import PDBParser


class SuperimposerTests(unittest.TestCase):
    """Test Superimposer module."""

    def test_Superimposer(self):
        """Test on module that superimpose two protein structures."""
        pdb1 = "PDB/1A8O.pdb"
        p = PDBParser()
        s1 = p.get_structure("FIXED", pdb1)
        fixed = Selection.unfold_entities(s1, "A")
        s2 = p.get_structure("MOVING", pdb1)
        moving = Selection.unfold_entities(s2, "A")
        rot = numpy.identity(3).astype('f')
        tran = numpy.array((1.0, 2.0, 3.0), 'f')
        for atom in moving:
            atom.transform(rot, tran)
        sup = Superimposer()
        sup.set_atoms(fixed, moving)
        self.assertTrue(numpy.allclose(sup.rotran[0], numpy.identity(3)))
        self.assertTrue(numpy.allclose(sup.rotran[1], numpy.array([-1.0, -2.0, -3.0])))
        self.assertAlmostEqual(sup.rms, 0.0, places=3)
        atom_list = ['N', 'C', 'C', 'O', 'C', 'C', 'SE', 'C', 'N', 'C', 'C',
                     'O', 'C', 'C', 'O', 'O', 'N', 'C', 'C', 'O', 'C', 'C',
                     'C', 'C', 'N', 'C', 'C', 'O', 'C', 'C', 'C', 'N', 'C',
                     'N', 'N', 'N', 'C', 'C', 'O', 'C', 'C', 'C', 'O', 'N',
                     'N', 'C', 'C', 'O', 'N', 'C', 'C', 'O', 'C', 'C', 'C',
                     'N', 'C', 'C', 'O', 'C', 'C', 'C', 'C', 'N', 'N', 'C',
                     'C', 'O', 'C', 'C', 'C', 'O', 'O', 'N', 'C', 'C', 'O',
                     'C', 'C', 'C', 'N', 'C', 'C', 'O', 'C', 'C', 'C', 'C',
                     'C', 'C', 'C', 'N', 'C', 'C', 'O', 'C', 'C', 'C', 'N',
                     'C', 'N', 'N', 'N', 'C', 'C', 'O', 'C', 'C', 'O', 'O',
                     'N', 'C', 'C', 'O', 'C', 'C', 'C', 'C', 'C', 'C', 'C',
                     'O', 'N', 'C', 'C', 'O', 'C', 'C', 'C', 'N', 'C', 'C',
                     'O', 'C', 'C', 'O', 'O', 'N', 'C', 'C', 'O', 'C', 'C',
                     'C', 'N', 'C', 'N', 'N', 'N', 'C', 'C', 'O', 'C', 'C',
                     'C', 'C', 'C', 'C', 'C', 'N', 'C', 'C', 'O', 'C', 'C',
                     'C', 'C', 'C', 'C', 'C', 'O', 'N', 'C', 'C', 'O', 'C',
                     'C', 'C', 'C', 'N', 'N', 'C', 'C', 'O', 'C', 'O', 'C',
                     'N', 'C', 'C', 'O', 'C', 'C', 'C', 'C', 'N', 'C', 'C',
                     'O', 'C', 'C', 'C', 'N', 'C', 'N', 'N', 'N', 'C', 'C',
                     'O', 'C', 'N', 'C', 'C', 'O', 'C', 'C', 'C', 'O', 'O',
                     'N', 'C', 'C', 'O', 'C', 'C', 'C', 'O', 'N', 'N', 'C',
                     'C', 'O', 'C', 'N', 'C', 'C', 'O', 'C', 'O', 'N', 'C',
                     'C', 'O', 'C', 'C', 'C', 'O', 'N', 'N', 'C', 'C', 'O',
                     'C', 'C', 'C', 'O', 'O', 'N', 'C', 'C', 'O', 'C', 'C',
                     'C', 'N', 'C', 'C', 'O', 'C', 'C', 'C', 'C', 'N', 'N',
                     'C', 'C', 'O', 'C', 'C', 'O', 'N', 'N', 'C', 'C', 'O',
                     'C', 'C', 'C', 'C', 'N', 'C', 'C', 'C', 'C', 'C', 'N',
                     'C', 'C', 'O', 'C', 'C', 'SE', 'C', 'N', 'C', 'C', 'O',
                     'C', 'O', 'C', 'N', 'C', 'C', 'O', 'C', 'C', 'C', 'O',
                     'O', 'N', 'C', 'C', 'O', 'C', 'O', 'C', 'N', 'C', 'C',
                     'O', 'C', 'C', 'C', 'C', 'N', 'C', 'C', 'O', 'C', 'C',
                     'C', 'C', 'N', 'C', 'C', 'O', 'C', 'C', 'C', 'N', 'C',
                     'C', 'O', 'C', 'C', 'C', 'O', 'N', 'N', 'C', 'C', 'O',
                     'C', 'C', 'O', 'N', 'N', 'C', 'C', 'O', 'C', 'N', 'C',
                     'C', 'O', 'C', 'C', 'O', 'N', 'N', 'C', 'C', 'O', 'C',
                     'C', 'C', 'N', 'C', 'C', 'O', 'C', 'C', 'O', 'O', 'N',
                     'C', 'C', 'O', 'C', 'S', 'N', 'C', 'C', 'O', 'C', 'C',
                     'C', 'C', 'N', 'N', 'C', 'C', 'O', 'C', 'O', 'C', 'N',
                     'C', 'C', 'O', 'C', 'C', 'C', 'C', 'N', 'C', 'C', 'O',
                     'C', 'C', 'C', 'C', 'N', 'C', 'C', 'O', 'C', 'C', 'C',
                     'C', 'N', 'N', 'C', 'C', 'O', 'C', 'N', 'C', 'C', 'O',
                     'C', 'C', 'C', 'C', 'N', 'C', 'C', 'O', 'N', 'C', 'C',
                     'O', 'C', 'C', 'C', 'N', 'C', 'C', 'O', 'N', 'C', 'C',
                     'O', 'C', 'N', 'C', 'C', 'O', 'C', 'O', 'C', 'N', 'C',
                     'C', 'O', 'C', 'C', 'C', 'C', 'N', 'C', 'C', 'O', 'C',
                     'C', 'C', 'O', 'O', 'N', 'C', 'C', 'O', 'C', 'C', 'C',
                     'O', 'O', 'N', 'C', 'C', 'O', 'C', 'C', 'SE', 'C', 'N',
                     'C', 'C', 'O', 'C', 'C', 'SE', 'C', 'N', 'C', 'C', 'O',
                     'C', 'O', 'C', 'N', 'C', 'C', 'O', 'C', 'N', 'C', 'C',
                     'O', 'C', 'S', 'N', 'C', 'C', 'O', 'C', 'C', 'C', 'O',
                     'N', 'N', 'C', 'C', 'O', 'O', 'O', 'O', 'O', 'O', 'O',
                     'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O',
                     'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O',
                     'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O',
                     'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O',
                     'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O',
                     'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O',
                     'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O',
                     'O', 'O', 'O', 'O', 'O', 'O']
        sup.apply(moving)
        atom_moved = []
        for aa in moving:
            atom_moved.append(aa.element)
        self.assertEqual(atom_moved, atom_list)


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
