"""
Run Bio.PDB.Selection tests.

Currently only tests unfold_entities.
"""

import itertools
import unittest
from typing import Union

from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.Atom import Atom
from Bio.PDB.Entity import Entity, DisorderedEntityWrapper
from Bio.PDB.PDBExceptions import PDBException
from Bio.PDB.Residue import Residue
from Bio.PDB.Selection import unfold_entities, _SelectParser, _AtomIndicator, select


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


class SelectParserTests(unittest.TestCase):
    def setUp(self):
        self.parser = _SelectParser()

    def test_identifiers(self):
        comparators = ("==", "!=", "<=", "<", ">=", ">")
        tests = [
            ("model 0", ("model", "0")),
            *(
                (f"model {comparator} 0", ("model", comparator, "0"))
                for comparator in comparators
            ),
            ("chain A", ("chain", "A")),
            *(
                (f"chain {comparator} A", ("chain", comparator, "A"))
                for comparator in comparators
            ),
            ("resn ALA", ("resn", "ALA")),
            *(
                (f"resn {comparator} ALA", ("resn", comparator, "ALA"))
                for comparator in ("==", "!=")
            ),
            ("resi 10", ("resi", "10")),
            *(
                (f"resi {comparator} 10", ("resi", comparator, "10"))
                for comparator in comparators
            ),
            ("name C", ("name", "C")),
            *(
                (f"name {comparator} C", ("name", comparator, "C"))
                for comparator in ("==", "!=")
            ),
        ]

        for query, expected_result in tests:
            result = tuple(self.parser(query)[0])
            self.assertEqual(result, expected_result)

    def test_boolean_operators(self):
        tests = (
            (
                "model 0 and chain B",
                [["model", "0"], "and", ["chain", "B"]],
            ),
            (
                "model 0 and (chain B or chain E)",
                [
                    ["model", "0"],
                    "and",
                    [["chain", "B"], "or", ["chain", "E"]],
                ],
            ),
            (
                "model 0 or chain B and not resn ALA",
                [
                    ["model", "0"],
                    "or",
                    [["chain", "B"], "and", ["not", ["resn", "ALA"]]],
                ],
            ),
            (
                "(model 0 or chain B) and resn ALA",
                [
                    [["model", "0"], "or", ["chain", "B"]],
                    "and",
                    ["resn", "ALA"],
                ],
            ),
            (
                "not chain A and chain B",
                [["not", ["chain", "A"]], "and", ["chain", "B"]],
            ),
            (
                "not (chain A and chain B)",
                ["not", [["chain", "A"], "and", ["chain", "B"]]],
            ),
        )

        for query, expected_result in tests:
            result = self.parser(query)[0].as_list()
            self.assertEqual(result, expected_result)

    def test_bad_input(self):
        try:
            import pyparsing as pp
        except ImportError:
            from Bio import MissingPythonDependencyError

            raise MissingPythonDependencyError(
                "Install pyparsing to use Bio.PDB.Selection (e.g. pip install pyparsing)"
            ) from None

        parser = self.parser

        with self.assertRaises(pp.ParseException):
            parser("bad input")
        with self.assertRaises(pp.ParseException):
            parser("model <=")
        with self.assertRaises(pp.ParseException):
            parser("chain chain A")


class AtomIndicatorTests(unittest.TestCase):
    def setUp(self):
        parser = _SelectParser()

        def create_indicator(query: str):
            return _AtomIndicator(parse_results=parser(query)[0])

        self.create_indicator = create_indicator
        structure_parser = MMCIFParser(QUIET=True)
        self.structure = structure_parser.get_structure("7CFN", "PDB/7CFN.cif")

    def test_query(self):
        structure = self.structure
        model = structure[0]
        indicator = self.create_indicator("not (chain A or chain B) and resn != CYS")

        for residue in itertools.chain(model["A"], model["B"]):
            for atom in residue:
                self.assertFalse(indicator(atom))

        for residue in itertools.chain(model["G"], model["N"], model["R"]):
            for atom in residue:
                self.assertEqual(indicator(atom), residue.resname != "CYS")

    def test_query2(self):
        structure = self.structure
        model = structure[0]
        indicator = self.create_indicator(
            "model 0 and not chain B and resi >= 10 and resi < 100 and name != C"
        )

        for chain in model:
            for residue in chain:
                for atom in residue:
                    self.assertEqual(
                        indicator(atom),
                        chain.id != "B"
                        and residue.id[1] in range(10, 100)
                        and atom.name != "C",
                    )


def get_atoms(entity: Union[Entity, DisorderedEntityWrapper]):
    """
    This helper is necessary because structure.get_atoms() does not yield
    all of the entities in a disordered entity.
    """
    if isinstance(entity, Atom):
        yield entity
    else:
        for child in entity.child_dict.values():
            yield from get_atoms(child)


class End2EndSelectTests(unittest.TestCase):
    def setUp(self):
        structure_parser = MMCIFParser(QUIET=True)
        # Note that this structure contains disordered residues and atoms.
        self.structure = structure_parser.get_structure("3JQH", "PDB/3JQH.cif")

    def test_query(self):
        structure = self.structure
        selection = select(structure, "resn PRO or resn GLU and name CA")
        atoms = get_atoms(selection)
        serial_numbers = {atom.serial_number for atom in atoms}
        expected = {1, 2, 3, 4, 5, 6, 7, 15, 71, 111, 155, 202}

        self.assertEqual(serial_numbers, expected)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
