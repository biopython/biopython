"""
Run Bio.PDB.Selection tests.

Currently only tests unfold_entities.
"""

import unittest

from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBException
from Bio.PDB.Selection import unfold_entities
from Bio.PDB.Residue import Residue


def res_full_id(res: Residue):
    """Return full residue identifier for thoroughly comparing residues."""
    return (res.get_resname(), *res.get_id())


class UnfoldEntitiesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.parser = PDBParser(PERMISSIVE=True)
        cls.structure = cls.parser.get_structure("scr", "PDB/1A8O.pdb")

    def test_from_structure_level(self):
        """Unfold from highest level to all levels."""
        struct_unfold = unfold_entities(self.structure, "S")[0]
        for res1, res2 in zip(
            self.structure.get_residues(), struct_unfold.get_residues()
        ):
            assert res_full_id(res1) == res_full_id(res2)

        model_unfold = unfold_entities(self.structure, "M")[0]
        for res1, res2 in zip(
            self.structure.get_residues(), model_unfold.get_residues()
        ):
            assert res_full_id(res1) == res_full_id(res2)

        residue_unfold = unfold_entities(self.structure, "R")
        for res1, res2 in zip(self.structure.get_residues(), residue_unfold):
            assert res_full_id(res1) == res_full_id(res2)

        atom_unfold = unfold_entities(self.structure, "A")
        for at1, at2 in zip(self.structure.get_atoms(), atom_unfold):
            assert at1 is at2

    def test_from_model_level(self):
        """Unfold from model to all levels."""
        structure_models = list(self.structure.get_models())

        struct_unfold = unfold_entities(structure_models, "S")[0]
        for res1, res2 in zip(
            self.structure.get_residues(), struct_unfold.get_residues()
        ):
            assert res_full_id(res1) == res_full_id(res2)

        model_unfold = unfold_entities(structure_models, "M")[0]
        for res1, res2 in zip(
            self.structure.get_residues(), model_unfold.get_residues()
        ):
            assert res_full_id(res1) == res_full_id(res2)

        residue_unfold = unfold_entities(structure_models, "R")
        for res1, res2 in zip(self.structure.get_residues(), residue_unfold):
            assert res_full_id(res1) == res_full_id(res2)

        atom_unfold = unfold_entities(structure_models, "A")
        for at1, at2 in zip(self.structure.get_atoms(), atom_unfold):
            assert at1 is at2

    def test_from_chain_level(self):
        """Unfold from chain level to all levels."""
        structure_chains = list(self.structure.get_chains())

        struct_unfold = unfold_entities(structure_chains, "S")[0]
        for res1, res2 in zip(
            self.structure.get_residues(), struct_unfold.get_residues()
        ):
            assert res_full_id(res1) == res_full_id(res2)

        model_unfold = unfold_entities(structure_chains, "M")[0]
        for res1, res2 in zip(
            self.structure.get_residues(), model_unfold.get_residues()
        ):
            assert res_full_id(res1) == res_full_id(res2)

        residue_unfold = unfold_entities(structure_chains, "R")
        for res1, res2 in zip(self.structure.get_residues(), residue_unfold):
            assert res_full_id(res1) == res_full_id(res2)

        atom_unfold = unfold_entities(structure_chains, "A")
        for at1, at2 in zip(self.structure.get_atoms(), atom_unfold):
            assert at1 is at2

    def test_from_residue_level(self):
        """Unfold from chain level to all levels."""
        structure_residues = list(self.structure.get_residues())

        struct_unfold = unfold_entities(structure_residues, "S")[0]
        for res1, res2 in zip(structure_residues, struct_unfold.get_residues()):
            assert res_full_id(res1) == res_full_id(res2)

        model_unfold = unfold_entities(structure_residues, "M")[0]
        for res1, res2 in zip(structure_residues, model_unfold.get_residues()):
            assert res_full_id(res1) == res_full_id(res2)

        residue_unfold = unfold_entities(structure_residues, "R")
        for res1, res2 in zip(structure_residues, residue_unfold):
            assert res_full_id(res1) == res_full_id(res2)

        atom_unfold = unfold_entities(structure_residues, "A")
        for at1, at2 in zip(self.structure.get_atoms(), atom_unfold):
            assert at1 is at2

    def test_from_atom_level(self):
        """Unfold from lowest level to all levels."""
        structure_atoms = list(self.structure.get_atoms())

        struct_unfold = unfold_entities(structure_atoms, "S")[0]
        for res1, res2 in zip(
            self.structure.get_residues(), struct_unfold.get_residues()
        ):
            assert res_full_id(res1) == res_full_id(res2)

        model_unfold = unfold_entities(structure_atoms, "M")[0]
        for res1, res2 in zip(
            self.structure.get_residues(), model_unfold.get_residues()
        ):
            assert res_full_id(res1) == res_full_id(res2)

        residue_unfold = unfold_entities(structure_atoms, "R")
        for res1, res2 in zip(self.structure.get_residues(), residue_unfold):
            assert res_full_id(res1) == res_full_id(res2)

        atom_unfold = unfold_entities(structure_atoms, "A")
        for at1, at2 in zip(structure_atoms, atom_unfold):
            assert at1 is at2

    def test_invalid_level(self):
        with self.assertRaises(PDBException):
            unfold_entities(self.structure, "Z")

    def test_entities_not_homogenous(self):
        structure_atom = next(self.structure.get_atoms())
        structure_chain = next(self.structure.get_chains())

        with self.assertRaises(PDBException):
            unfold_entities([structure_atom, structure_chain], "A")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
