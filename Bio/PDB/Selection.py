# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Selection of atoms, residues, etc."""

import itertools
from copy import copy

import numpy as np

from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Entity import Entity
from Bio.PDB.Model import Model
from Bio.PDB.PDBExceptions import PDBException
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure

try:
    import pyparsing as pp
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install pyparsing to use Bio.PDB.Selection (e.g. pip install pyparsing)"
    ) from None

entity_levels = ["A", "R", "C", "M", "S"]


def uniqueify(items):
    """Return a list of the unique items in the given iterable.

    Order is NOT preserved.
    """
    return list(set(items))


def get_unique_parents(entity_list):
    """Translate a list of entities to a list of their (unique) parents."""
    unique_parents = {entity.get_parent() for entity in entity_list}
    return list(unique_parents)


def unfold_entities(entity_list, target_level):
    """Unfold entities list to a child level (e.g. residues in chain).

    Unfold a list of entities to a list of entities of another
    level.  E.g.:

    list of atoms -> list of residues
    list of modules -> list of atoms
    list of residues -> list of chains

    - entity_list - list of entities or a single entity
    - target_level - char (A, R, C, M, S)

    Note that if entity_list is an empty list, you get an empty list back:

    >>> unfold_entities([], "A")
    []

    """
    if target_level not in entity_levels:
        raise PDBException(f"{target_level}: Not an entity level.")
    if entity_list == []:
        return []
    if isinstance(entity_list, (Entity, Atom)):
        entity_list = [entity_list]

    level = entity_list[0].get_level()
    if not all(entity.get_level() == level for entity in entity_list):
        raise PDBException("Entity list is not homogeneous.")

    target_index = entity_levels.index(target_level)
    level_index = entity_levels.index(level)

    if level_index == target_index:  # already right level
        return entity_list

    entities = entity_list

    if level_index > target_index:  # we're going down, e.g. S->A
        for i in range(target_index, level_index):
            entities = itertools.chain.from_iterable(entities)
    else:  # we're going up, e.g. A->S
        for i in range(level_index, target_index):
            # get unique parents by removing duplicates while preserving order
            entities = {entity.get_parent(): None for entity in entities}

    return list(entities)


class _SelectParser:
    def __init__(self):
        model_pattern = pp.Keyword("model") + pp.Word(pp.alphanums)
        chain_pattern = pp.Keyword("chain") + pp.Word(pp.alphanums)
        resn_pattern = pp.Keyword("resn") + pp.Word(pp.alphanums)
        resi_pattern = pp.Keyword("resi") + pp.Word(pp.alphanums)
        name_pattern = pp.Keyword("name") + pp.Word(pp.alphanums)
        identifier_pattern = pp.Group(
            model_pattern | chain_pattern | resn_pattern | resi_pattern | name_pattern
        ).set_results_name("identifier")

        parentheses_pattern = pp.Forward()
        not_pattern = pp.Forward()
        and_pattern = pp.Forward()
        or_pattern = pp.Forward()

        parentheses_pattern <<= (
            pp.Group("(" + or_pattern + ")").set_results_name("parentheses")
            | identifier_pattern
        )
        not_pattern <<= (
            pp.Group(pp.Keyword("not") + not_pattern).set_results_name("not")
            | parentheses_pattern
        )
        and_pattern <<= (
            pp.Group(not_pattern + pp.Keyword("and") + and_pattern).set_results_name(
                "and"
            )
            | not_pattern
        )
        or_pattern <<= (
            pp.Group(and_pattern + pp.Keyword("or") + or_pattern).set_results_name("or")
            | and_pattern
        )

        self._parser = or_pattern.parse_string

    def __call__(self, *args, **kwargs):
        assert len(args) == 1
        assert not kwargs
        return self.parse(args[0])

    def parse(self, statement: str):
        return self._parser(statement)


class _SelectEvaluator:
    def __init__(self, *, structure: Structure):
        self._structure = structure
        self._atoms = list(structure.get_atoms())
        self._evaluators = {
            "identifier": self._evaluate_identifier,
            "parentheses": self._evaluate_parentheses,
            "not": self._evaluate_not,
            "and": self._evaluate_and,
            "or": self._evaluate_or,
        }

    def __call__(self, *args, **kwargs):
        assert len(args) == 1
        assert not kwargs
        return self.evaluate(args[0])

    def _evaluate_identifier(self, parse_result: pp.ParseResults) -> np.ndarray:
        assert len(parse_result) == 2

        id_type = parse_result[0]
        id_value = parse_result[1]

        if id_type == "model":
            atom_models = [
                atom.get_parent().get_parent().get_parent() for atom in self._atoms
            ]
            match_indicators = [str(model.id) == id_value for model in atom_models]
        elif id_type == "chain":
            atom_chains = [atom.get_parent().get_parent() for atom in self._atoms]
            match_indicators = [chain.id == id_value for chain in atom_chains]
        elif id_type == "resn":
            atom_residues = [atom.get_parent() for atom in self._atoms]
            match_indicators = [
                residue.get_resname() == id_value for residue in atom_residues
            ]
        elif id_type == "resi":
            # TODO: Residues have weird IDs
            atom_residues = [atom.get_parent() for atom in self._atoms]
            match_indicators = [
                residue.get_id() == id_value for residue in atom_residues
            ]
        elif id_type == "name":
            match_indicators = [atom.get_name() == id_value for atom in self._atoms]
        else:
            raise ValueError(f"Unexpected identifier type: {id_type}")

        return np.array(match_indicators)

    def _evaluate_parentheses(self, parse_result: pp.ParseResults) -> np.ndarray:
        assert len(parse_result) == 3
        assert parse_result[0] == "("
        assert parse_result[2] == ")"

        return self.evaluate(parse_result[1])

    def _evaluate_not(self, parse_result: pp.ParseResults) -> np.ndarray:
        assert len(parse_result) == 2
        assert parse_result[0] == "not"

        operand = self.evaluate(parse_result[1])

        return np.logical_not(operand)

    def _evaluate_and(self, parse_result: pp.ParseResults) -> np.ndarray:
        assert len(parse_result) == 3
        assert parse_result[1] == "and"

        left_operand = self.evaluate(parse_result[0])
        right_operaand = self.evaluate(parse_result[2])

        return np.logical_and(left_operand, right_operaand)

    def _evaluate_or(self, parse_result: pp.ParseResults) -> np.ndarray:
        assert len(parse_result) == 3
        assert parse_result[1] == "or"

        left_operand = self.evaluate(parse_result[0])
        right_operaand = self.evaluate(parse_result[2])

        return np.logical_or(left_operand, right_operaand)

    def evaluate(self, parse_result: pp.ParseResults) -> np.ndarray:
        operation_name = parse_result.get_name()
        evaluator = self._evaluators[operation_name]
        return evaluator(parse_result)


def select(structure: Structure, statement: str) -> Structure:
    """
    TODO: Fill out docstring
    """
    parser = _SelectParser()
    evaluator = _SelectEvaluator(structure=structure)
    parse_result = parser(statement)[0]
    evaluate_result = evaluator(parse_result)

    result_structure = Structure(structure.id)
    result_model = None
    result_chain = None
    result_residue = None

    for match_indicator, atom in zip(evaluate_result, structure.get_atoms()):
        if not match_indicator:
            continue

        residue = atom.get_parent()
        chain = residue.get_parent()
        model = chain.get_parent()

        assert isinstance(residue, Residue)
        assert isinstance(chain, Chain)
        assert isinstance(model, Model)
        # TODO: Handle disordered entities

        if not result_model or result_model.id != model.id:
            result_model = Model(model.id, model.serial_num)
            result_structure.add(result_model)
            result_chain = None
            result_residue = None
        if not result_chain or result_chain.id != chain.id:
            if result_model.has_id(chain.id):
                result_chain = result_model[chain.id]
            else:
                result_chain = Chain(chain.id)
                result_model.add(result_chain)
            result_residue = None
        if not result_residue or result_residue.id != residue.id:
            if result_chain.has_id(residue.id):
                result_residue = result_chain[residue.id]
            else:
                result_residue = Residue(residue.id, residue.resname, residue.segid)
                result_chain.add(result_residue)

        atom_copy = copy(atom)
        assert atom_copy is not atom
        atom_copy.set_parent(result_residue)
        result_residue.add(atom_copy)

    return result_structure


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
