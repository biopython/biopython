# Copyright 2009-2011 by Eric Talevich.  All rights reserved.
# Revisions copyright 2009-2013 by Peter Cock.  All rights reserved.
# Revisions copyright 2013 Lenna X. Peterson. All rights reserved.
# Revisions copyright 2020 Joao Rodrigues. All rights reserved.
#
# Converted by Eric Talevich from an older unit test copyright 2002
# by Thomas Hamelryck.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Generic unit tests for the SMCRA classes of the Bio.PDB module."""

from copy import deepcopy
import unittest
import warnings

try:
    import numpy
    from numpy import dot  # Missing on old PyPy's micronumpy

    del dot
    from numpy.linalg import svd, det  # Missing in PyPy 2.0 numpypy

    del svd, det
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.PDB."
    ) from None

from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.PDB import rotmat, Vector
from Bio.PDB import Atom


class Atom_Element(unittest.TestCase):
    """induces Atom Element from Atom Name."""

    def test_atom_element_assignment(self):
        """Atom Element."""
        parser = PDBParser(PERMISSIVE=True, QUIET=True)
        structure = parser.get_structure("X", "PDB/a_structure.pdb")
        residue = structure[0]["A"][("H_PCA", 1, " ")]

        atoms = residue.child_list
        self.assertEqual("N", atoms[0].element)  # N
        self.assertEqual("C", atoms[1].element)  # Alpha Carbon
        self.assertEqual("D", atoms[4].element)  # Deuterium
        self.assertEqual("CA", atoms[8].element)  # Calcium

    def test_assign_unknown_element(self):
        """Unknown element is assigned 'X'."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            a = Atom.Atom(
                "XE1", None, None, None, None, " XE1", None  # serial 5170 - 4CP4
            )
        self.assertEqual(a.element, "X")

    def test_ions(self):
        """Element for magnesium is assigned correctly."""
        parser = PDBParser(PERMISSIVE=True)
        structure = parser.get_structure("X", "PDB/ions.pdb")
        # check magnesium atom
        atoms = structure[0]["A"][("H_MG", 1, " ")].child_list
        self.assertEqual("MG", atoms[0].element)

    def test_hydrogens(self):
        def quick_assign(fullname):
            return Atom.Atom(
                fullname.strip(), None, None, None, None, fullname, None
            ).element

        pdb_elements = {
            "H": (
                " H  ",
                " HA ",
                " HB ",
                " HD1",
                " HD2",
                " HE ",
                " HE1",
                " HE2",
                " HE3",
                " HG ",
                " HG1",
                " HH ",
                " HH2",
                " HZ ",
                " HZ2",
                " HZ3",
                "1H  ",
                "1HA ",
                "1HB ",
                "1HD ",
                "1HD1",
                "1HD2",
                "1HE ",
                "1HE2",
                "1HG ",
                "1HG1",
                "1HG2",
                "1HH1",
                "1HH2",
                "1HZ ",
                "2H  ",
                "2HA ",
                "2HB ",
                "2HD ",
                "2HD1",
                "2HD2",
                "2HE ",
                "2HE2",
                "2HG ",
                "2HG1",
                "2HG2",
                "2HH1",
                "2HH2",
                "2HZ ",
                "3H  ",
                "3HB ",
                "3HD1",
                "3HD2",
                "3HE ",
                "3HG1",
                "3HG2",
                "3HZ ",
                "HE21",
            ),
            "O": (" OH ",),  # noqa: E741
            "C": (" CH2",),
            "N": (" NH1", " NH2"),
        }

        for element, atom_names in pdb_elements.items():
            for fullname in atom_names:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", PDBConstructionWarning)
                    e = quick_assign(fullname)
                self.assertEqual(e, element)


class SortingTests(unittest.TestCase):
    """Tests for sorting elements of the SMCRA representation."""

    def test_residue_sort(self):
        """Test atoms are sorted correctly in residues."""
        parser = PDBParser()
        structure = parser.get_structure("example", "PDB/1A8O.pdb")

        for residue in structure.get_residues():
            old = [a.name for a in residue]
            new = [a.name for a in sorted(residue)]
            special = []
            for a in ["N", "CA", "C", "O"]:
                if a in old:
                    special.append(a)
            special_len = len(special)
            self.assertEqual(
                new[0:special_len],
                special,
                f"Sorted residue did not place N, CA, C, O first: {new}",
            )
            self.assertEqual(
                new[special_len:],
                sorted(new[special_len:]),
                f"After N, CA, C, O should be alphabet: {new}",
            )

    # Tests for sorting methods
    def test_comparison_entities(self):
        """Test comparing and sorting the several SMCRA objects."""
        parser = PDBParser(QUIET=True)
        struct = parser.get_structure("example", "PDB/a_structure.pdb")

        # Test deepcopy of a structure with disordered atoms
        struct2 = deepcopy(struct)

        # Sorting (<, >, <=, <=)
        # Chains (same code as models)
        model = struct[1]
        chains = [c.id for c in sorted(model)]
        self.assertEqual(chains, ["A", "B", "C", " "])
        # Residues
        residues = [r.id[1] for r in sorted(struct[1]["C"])]
        self.assertEqual(residues, [1, 2, 3, 4, 0])
        # Atoms
        for residue in struct.get_residues():
            old = [a.name for a in residue]
            new = [a.name for a in sorted(residue)]

            special = [a for a in ("N", "CA", "C", "O") if a in old]
            len_special = len(special)
            # Placed N, CA, C, O first?
            self.assertEqual(
                new[:len_special],
                special,
                f"Sorted residue did not place N, CA, C, O first: {new}",
            )
            # Placed everyone else alphabetically?
            self.assertEqual(
                new[len_special:],
                sorted(new[len_special:]),
                f"After N, CA, C, O order Should be alphabetical: {new}",
            )
        # DisorderedResidue
        residues = [r.id[1] for r in sorted(struct[1]["A"])][79:81]
        self.assertEqual(residues, [80, 81])
        # Insertion code + hetflag + chain
        residues = list(struct[1]["B"]) + [struct[1]["A"][44]]
        self.assertEqual(
            [("{}" * 4).format(r.parent.id, *r.id) for r in sorted(residues)],
            [
                "A 44 ",
                "B 44 ",
                "B 46 ",
                "B 47 ",
                "B 48 ",
                "B 49 ",
                "B 50 ",
                "B 51 ",
                "B 51A",
                "B 52 ",
                "BH_SEP45 ",
                "BW0 ",
            ],
        )
        # DisorderedAtom
        atoms = [a.altloc for a in sorted(struct[1]["A"][74]["OD1"])]
        self.assertEqual(atoms, ["A", "B"])

        # Comparisons
        # Structure
        self.assertEqual(struct, struct2)
        self.assertLessEqual(struct, struct2)
        self.assertGreaterEqual(struct, struct2)
        struct2.id = "new_id"
        self.assertNotEqual(struct, struct2)
        self.assertLess(struct, struct2)
        self.assertLessEqual(struct, struct2)
        self.assertGreater(struct2, struct)
        self.assertGreaterEqual(struct2, struct)

        # Model

        self.assertEqual(model, model)  # __eq__ same type
        self.assertNotEqual(struct[0], struct[1])

        self.assertNotEqual(struct[0], [])  # __eq__ diff. types
        self.assertNotEqual(struct, model)

        # residues with same ID string should not be equal if the parent is not equal
        res1, res2, res3 = residues[0], residues[-1], struct2[1]["A"][44]
        self.assertEqual(res1.id, res2.id)
        self.assertEqual(
            res2, res3
        )  # Equality of identical residues with different structure ID
        self.assertNotEqual(res1, res2)
        self.assertGreater(res1, res2)
        self.assertGreaterEqual(res1, res2)
        self.assertLess(res2, res1)
        self.assertLessEqual(res2, res1)

        # atom should not be equal if the parent is not equal
        atom1, atom2, atom3 = res1["CA"], res2["CA"], res3["CA"]
        self.assertEqual(
            atom2, atom3
        )  # Equality of identical atoms with different structure ID
        self.assertGreater(atom1, atom2)
        self.assertGreaterEqual(atom1, atom2)
        self.assertGreaterEqual(atom2, atom3)
        self.assertNotEqual(atom1, atom2)
        self.assertLess(atom2, atom1)
        self.assertLessEqual(atom2, atom1)
        self.assertLessEqual(atom2, atom3)


class IterationTests(unittest.TestCase):
    """Tests iterating over the SMCRA hierarchy."""

    @classmethod
    def setUpClass(cls):
        parser = PDBParser(PERMISSIVE=True, QUIET=True)
        cls.structure = parser.get_structure("X", "PDB/a_structure.pdb")

    def test_get_chains(self):
        """Yields chains from different models separately."""
        chains = [chain.id for chain in self.structure.get_chains()]
        self.assertEqual(chains, ["A", "A", "B", "C", " "])

    def test_get_residues(self):
        """Yields all residues from all models."""
        residues = [resi.id for resi in self.structure.get_residues()]
        self.assertEqual(len(residues), 179)

    def test_get_atoms(self):
        """Yields all atoms from the structure, excluding duplicates and ALTLOCs which are not parsed."""
        atoms = [
            "%12s" % str((atom.id, atom.altloc)) for atom in self.structure.get_atoms()
        ]
        self.assertEqual(len(atoms), 835)


class ChangingIdTests(unittest.TestCase):
    """Tests changing properties of SMCRA objects."""

    def setUp(self):
        parser = PDBParser(PERMISSIVE=True, QUIET=True)
        self.structure = parser.get_structure("X", "PDB/a_structure.pdb")

    def test_change_model_id(self):
        """Change the id of a model."""
        for model in self.structure:
            break  # Get first model in structure
        model.id = 2
        self.assertEqual(model.id, 2)
        self.assertIn(2, self.structure)
        self.assertNotIn(0, self.structure)

    def test_change_model_id_raises(self):
        """Cannot change id to a value already in use by another child."""
        model = next(iter(self.structure))
        with self.assertRaises(ValueError):
            model.id = 1
        # Make sure nothing was changed
        self.assertEqual(model.id, 0)
        self.assertIn(0, self.structure)
        self.assertIn(1, self.structure)

    def test_change_chain_id(self):
        """Change the id of a model."""
        chain = next(iter(self.structure.get_chains()))
        chain.id = "R"
        self.assertEqual(chain.id, "R")
        model = next(iter(self.structure))
        self.assertIn("R", model)

    def test_change_id_to_self(self):
        """Changing the id to itself does nothing (does not raise)."""
        chain = next(iter(self.structure.get_chains()))
        chain_id = chain.id
        chain.id = chain_id
        self.assertEqual(chain.id, chain_id)

    def test_change_residue_id(self):
        """Change the id of a residue."""
        chain = next(iter(self.structure.get_chains()))
        res = chain[("H_PCA", 1, " ")]
        res.id = (" ", 1, " ")

        self.assertEqual(res.id, (" ", 1, " "))
        self.assertIn((" ", 1, " "), chain)
        self.assertNotIn(("H_PCA", 1, " "), chain)
        self.assertEqual(chain[(" ", 1, " ")], res)

    def test_full_id_is_updated_residue(self):
        """Invalidate cached full_ids if an id is changed."""
        atom = next(iter(self.structure.get_atoms()))

        # Generate the original full id.
        original_id = atom.get_full_id()
        self.assertEqual(original_id, ("X", 0, "A", ("H_PCA", 1, " "), ("N", " ")))
        residue = next(iter(self.structure.get_residues()))

        # Make sure the full id was in fact cached,
        # so we need to invalidate it later.
        self.assertEqual(residue.full_id, ("X", 0, "A", ("H_PCA", 1, " ")))

        # Changing the residue's id should lead to an updated full id.
        residue.id = (" ", 1, " ")
        new_id = atom.get_full_id()
        self.assertNotEqual(original_id, new_id)
        self.assertEqual(new_id, ("X", 0, "A", (" ", 1, " "), ("N", " ")))

    def test_full_id_is_updated_chain(self):
        """Invalidate cached full_ids if an id is changed."""
        atom = next(iter(self.structure.get_atoms()))

        # Generate the original full id.
        original_id = atom.get_full_id()
        self.assertEqual(original_id, ("X", 0, "A", ("H_PCA", 1, " "), ("N", " ")))
        residue = next(iter(self.structure.get_residues()))

        # Make sure the full id was in fact cached,
        # so we need to invalidate it later.
        self.assertEqual(residue.full_id, ("X", 0, "A", ("H_PCA", 1, " ")))
        chain = next(iter(self.structure.get_chains()))

        # Changing the chain's id should lead to an updated full id.
        chain.id = "Q"
        new_id = atom.get_full_id()
        self.assertNotEqual(original_id, new_id)
        self.assertEqual(new_id, ("X", 0, "Q", ("H_PCA", 1, " "), ("N", " ")))


class TransformTests(unittest.TestCase):
    """Tests transforming the coordinates of Atoms in a structure."""

    def setUp(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            self.s = PDBParser(PERMISSIVE=True).get_structure(
                "X", "PDB/a_structure.pdb"
            )
        self.m = self.s.get_list()[0]
        self.c = self.m.get_list()[0]
        self.r = self.c.get_list()[0]
        self.a = self.r.get_list()[0]

    def get_total_pos(self, o):
        """Sum of positions of atoms in an entity along with the number of atoms."""
        if hasattr(o, "get_coord"):
            return o.get_coord(), 1
        total_pos = numpy.array((0.0, 0.0, 0.0))
        total_count = 0
        for p in o.get_list():
            pos, count = self.get_total_pos(p)
            total_pos += pos
            total_count += count
        return total_pos, total_count

    def get_pos(self, o):
        """Average atom position in an entity."""
        pos, count = self.get_total_pos(o)
        return pos / count

    def test_transform(self):
        """Transform entities (rotation and translation)."""
        for o in (self.s, self.m, self.c, self.r, self.a):
            rotation = rotmat(Vector(1, 3, 5), Vector(1, 0, 0))
            translation = numpy.array((2.4, 0, 1), "f")
            oldpos = self.get_pos(o)
            o.transform(rotation, translation)
            newpos = self.get_pos(o)
            newpos_check = numpy.dot(oldpos, rotation) + translation
            for i in range(0, 3):
                self.assertAlmostEqual(newpos[i], newpos_check[i])


class CopyTests(unittest.TestCase):
    """Tests copying SMCRA objects."""

    def setUp(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            self.s = PDBParser(PERMISSIVE=True).get_structure(
                "X", "PDB/a_structure.pdb"
            )
        self.m = self.s.get_list()[0]
        self.c = self.m.get_list()[0]
        self.r = self.c.get_list()[0]
        self.a = self.r.get_list()[0]

    def test_atom_copy(self):
        aa = self.a.copy()
        self.assertIsNot(self.a, aa)
        self.assertIsNot(self.a.get_coord(), aa.get_coord())

    def test_entity_copy(self):
        """Make a copy of a residue."""
        for e in (self.s, self.m, self.c, self.r):
            ee = e.copy()
            self.assertIsNot(e, ee)
            self.assertIsNot(e.get_list()[0], ee.get_list()[0])


class CenterOfMassTests(unittest.TestCase):
    """Tests calculating centers of mass/geometry."""

    @classmethod
    def setUpClass(cls):
        cls.parser = parser = PDBParser()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            cls.structure = parser.get_structure("a", "PDB/1LCD.pdb")

    def test_structure_com(self):
        """Calculate Structure center of mass."""
        com = self.structure.center_of_mass()

        self.assertTrue(numpy.allclose(com, [19.870, 25.455, 28.753], atol=1e-3))

    def test_structure_cog(self):
        """Calculate Structure center of geometry."""
        cog = self.structure.center_of_mass(geometric=True)

        self.assertTrue(numpy.allclose(cog, [19.882, 25.842, 28.333], atol=1e-3))

    def test_chain_cog(self):
        """Calculate center of geometry of individual chains."""
        expected = {
            "A": [20.271, 30.191, 23.563],
            "B": [19.272, 21.163, 33.711],
            "C": [19.610, 20.599, 32.708],
        }

        for chain in self.structure[0].get_chains():  # one model only
            cog = chain.center_of_mass(geometric=True)
            self.assertTrue(numpy.allclose(cog, expected[chain.id], atol=1e-3))

    def test_com_empty_structure(self):
        """Center of mass of empty structure raises ValueError."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            s = self.parser.get_structure("b", "PDB/disordered.pdb")  # smaller

        for child in list(s):
            s.detach_child(child.id)

        with self.assertRaises(ValueError):
            s.center_of_mass()


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
