# Revisions copyright 2020 Joao Rodrigues. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Unit tests for disordered atoms in the Bio.PDB module."""

import os
import tempfile
import unittest

from Bio.PDB import PDBParser, PDBIO


class TestDisordered(unittest.TestCase):
    """Tests for operations on DisorderedEntities."""

    def unpack_all_atoms(self, structure):
        """Return a list of all atoms in the structure."""
        return [a for r in structure.get_residues() for a in r.get_unpacked_list()]

    def test_copy_disordered_atom(self):
        """Copies disordered atoms and all their children."""
        parser = PDBParser(QUIET=1)
        s = parser.get_structure("x", "PDB/disordered.pdb")

        resi27 = s[0]["A"][27]
        resi27_copy = resi27.copy()

        self.assertNotEqual(id(resi27), id(resi27_copy))  # did we really copy

        resi27_atoms = resi27.get_unpacked_list()
        resi27_copy_atoms = resi27.get_unpacked_list()
        self.assertEqual(len(resi27_atoms), len(resi27_copy_atoms))

        for ai, aj in zip(resi27_atoms, resi27_copy_atoms):
            self.assertEqual(ai.name, aj.name)

    def test_copy_entire_chain(self):
        """Copy propagates throughout SMCRA object."""
        parser = PDBParser(QUIET=1)
        s = parser.get_structure("x", "PDB/disordered.pdb")

        s_copy = s.copy()

        self.assertNotEqual(id(s), id(s_copy))  # did we really copy

        atoms = self.unpack_all_atoms(s)
        copy_atoms = self.unpack_all_atoms(s_copy)
        self.assertEqual(len(atoms), len(copy_atoms))

        for ai, aj in zip(atoms, copy_atoms):
            self.assertEqual(ai.name, aj.name)

    def test_transform_disordered(self):
        """Transform propagates through disordered atoms."""
        # This test relates to issue #455 where applying a transformation
        # to a copied structure did not work for disordered atoms.
        parser = PDBParser(QUIET=1)
        s = parser.get_structure("x", "PDB/disordered.pdb")

        s_copy = s.copy()

        mtx = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
        tr_vec = (20.0, 0.0, 0.0)

        s_copy.transform(mtx, tr_vec)  # transform copy

        atoms = self.unpack_all_atoms(s)
        copy_atoms = self.unpack_all_atoms(s_copy)
        self.assertEqual(len(atoms), len(copy_atoms))
        for ai, aj in zip(atoms, copy_atoms):
            self.assertEqual(ai - aj, 20.0)  # check distance == 20.0

    def test_copy_and_write_disordered(self):
        """Extract, save, and parse again disordered atoms."""
        parser = PDBParser(QUIET=1)
        writer = PDBIO()

        s = parser.get_structure("x", "PDB/disordered.pdb")

        # Extract the chain object
        chain = s[0]["A"]

        writer.set_structure(chain)

        filenumber, filename = tempfile.mkstemp()  # save to temp file
        os.close(filenumber)
        try:
            writer.save(filename)

            # Parse again
            s2 = parser.get_structure("x_copy", filename)

            # Do we have the same stuff?
            atoms1 = self.unpack_all_atoms(s)
            atoms2 = self.unpack_all_atoms(s2)
            self.assertEqual(len(atoms1), len(atoms2))
            for ai, aj in zip(atoms1, atoms2):
                self.assertEqual(ai.name, aj.name)

        finally:
            os.remove(filename)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
