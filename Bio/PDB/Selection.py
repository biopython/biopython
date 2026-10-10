# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Selection of atoms, residues, etc."""
import functools
import itertools
import operator
from collections.abc import Callable
from typing import Union

from Bio.PDB.Atom import Atom, DisorderedAtom
from Bio.PDB.Entity import Entity, DisorderedEntityWrapper
from Bio.PDB.PDBExceptions import PDBException
from Bio.PDB.Residue import Residue, DisorderedResidue
from Bio.PDB.Structure import Structure

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
    """
    A parser for atom selection statements.

    The parser accepts a query string as input
    and returns a parse tree that represents the query.
    """

    def __init__(self):
        try:
            import pyparsing as pp
        except ImportError:
            from Bio import MissingPythonDependencyError

            raise MissingPythonDependencyError(
                "Install pyparsing to use Bio.PDB.Selection (e.g. pip install pyparsing)"
            ) from None

        equality_pattern = pp.Keyword("==") | pp.Keyword("!=")
        subordinate_pattern = pp.Keyword("<=") | pp.Keyword("<")
        superordinate_pattern = pp.Keyword(">=") | pp.Keyword(">")
        comparator_pattern = (
            equality_pattern | subordinate_pattern | superordinate_pattern
        )

        model_pattern = (
            pp.Keyword("model") + pp.Opt(comparator_pattern) + pp.Word(pp.alphanums)
        )
        chain_pattern = (
            pp.Keyword("chain") + pp.Opt(comparator_pattern) + pp.Word(pp.alphanums)
        )
        resn_pattern = (
            pp.Keyword("resn") + pp.Opt(equality_pattern) + pp.Word(pp.alphanums)
        )
        resi_pattern = (
            pp.Keyword("resi") + pp.Opt(comparator_pattern) + pp.Word(pp.alphanums)
        )
        name_pattern = (
            pp.Keyword("name") + pp.Opt(equality_pattern) + pp.Word(pp.alphanums)
        )
        identifier_pattern = pp.Group(
            model_pattern | chain_pattern | resn_pattern | resi_pattern | name_pattern
        ).set_results_name("identifier")

        parentheses_pattern = pp.Forward()
        not_pattern = pp.Forward()
        and_pattern = pp.Forward()
        or_pattern = pp.Forward()

        parentheses_pattern <<= (
            pp.Suppress("(") + or_pattern + pp.Suppress(")")
        ).set_results_name("parentheses") | identifier_pattern
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

        self._parser = functools.partial(or_pattern.parse_string, parse_all=True)

    def __call__(self, *args, **kwargs):
        assert len(args) == 1
        assert not kwargs
        return self.parse(args[0])

    def parse(self, statement: str):
        return self._parser(statement)


class _AtomIndicator:
    """
    An atom indicator is a function that accepts an atom argument
    and returns a boolean result, indicating whether the query selects the atom.

    The atom indicator initializer accepts a parse tree representing the query.
    """

    def __init__(self, *, parse_results):
        self._indicator_composers = {
            "identifier": self._compose_identifier_indicator,
            "parentheses": self._compose_parentheses_indicator,
            "not": self._compose_not_indicator,
            "and": self._compose_and_indicator,
            "or": self._compose_or_indicator,
        }
        self._comparators = {
            "==": operator.eq,
            "!=": operator.ne,
            "<=": operator.le,
            "<": operator.lt,
            ">=": operator.ge,
            ">": operator.gt,
        }
        self._indicator = self._compose_indicator(parse_results)

    def _compose_identifier_indicator(self, parse_results) -> Callable[[Atom], bool]:
        assert len(parse_results) in {2, 3}

        id_type = parse_results[0]
        id_value = parse_results[-1]
        comparator = parse_results[1] if len(parse_results) == 3 else "=="
        comparator = self._comparators[comparator]

        if id_type == "model":

            def indicator(atom: Atom) -> bool:
                model = atom.get_parent().get_parent().get_parent()
                return comparator(model.id, int(id_value))

        elif id_type == "chain":

            def indicator(atom: Atom) -> bool:
                chain = atom.get_parent().get_parent()
                return comparator(chain.id, id_value)

        elif id_type == "resn":

            def indicator(atom: Atom) -> bool:
                residue = atom.get_parent()
                return comparator(residue.resname, id_value)

        elif id_type == "resi":

            def indicator(atom: Atom) -> bool:
                residue = atom.get_parent()
                return comparator(residue.id[1], int(id_value))

        elif id_type == "name":

            def indicator(atom: Atom) -> bool:
                return comparator(atom.name, id_value)

        else:
            raise ValueError(f"Unexpected identifier type: {id_type}")

        return indicator

    def _compose_parentheses_indicator(self, parse_results) -> Callable[[Atom], bool]:
        assert len(parse_results) == 3
        assert parse_results[0] == "("
        assert parse_results[2] == ")"

        return self._compose_indicator(parse_results[1])

    def _compose_not_indicator(self, parse_results) -> Callable[[Atom], bool]:
        assert len(parse_results) == 2
        assert parse_results[0] == "not"

        operand_indicator = self._compose_indicator(parse_results[1])

        def indicator(atom: Atom) -> bool:
            return not operand_indicator(atom)

        return indicator

    def _compose_and_indicator(self, parse_results) -> Callable[[Atom], bool]:
        assert len(parse_results) == 3
        assert parse_results[1] == "and"

        left_indicator = self._compose_indicator(parse_results[0])
        right_indicator = self._compose_indicator(parse_results[2])

        def indicator(atom: Atom) -> bool:
            return left_indicator(atom) and right_indicator(atom)

        return indicator

    def _compose_or_indicator(self, parse_results) -> Callable[[Atom], bool]:
        assert len(parse_results) == 3
        assert parse_results[1] == "or"

        left_indicator = self._compose_indicator(parse_results[0])
        right_indicator = self._compose_indicator(parse_results[2])

        def indicator(atom: Atom) -> bool:
            return left_indicator(atom) or right_indicator(atom)

        return indicator

    def _compose_indicator(self, parse_results) -> Callable[[Atom], bool]:
        operation_name = parse_results.get_name()
        composer = self._indicator_composers[operation_name]
        return composer(parse_results)

    def __call__(self, *args, **kwargs):
        assert len(args) == 1
        assert not kwargs

        atom = args[0]

        return self._indicator(atom)


class _StructurePruner:
    """
    A structure pruner uses an atom indicator to recursively remove
    entities that are not selected by the query from the structure.
    """

    def __init__(self, indicator: Callable[[Atom], bool]):
        self._atom_indicator = indicator

    def __call__(self, *args, **kwargs):
        assert len(args) == 1
        assert not kwargs

        self.prune(args[0])

    def prune(self, entity: Union[Entity, DisorderedEntityWrapper]) -> bool:
        indicator = self._atom_indicator

        if isinstance(entity, Atom):
            return not indicator(entity)

        pruned_all_children = True
        children = tuple(entity.child_dict.items())

        for child_id, child in children:
            if self.prune(child):
                if isinstance(entity, DisorderedResidue):
                    assert isinstance(child, Residue)
                    entity.disordered_remove(child.resname)
                elif isinstance(entity, DisorderedAtom):
                    assert isinstance(child, Atom)
                    entity.disordered_remove(child.altloc)
                else:
                    entity.detach_child(child_id)
            else:
                pruned_all_children = False

        return pruned_all_children


def select(structure: Structure, query: str) -> Structure:
    """
    Select atoms in a structure with a query,
    returning a new structure with only the selected atoms.

    :param Structure structure: the PDB structure
    :param str query: the query
    :return: the new PDB structure
    :rtype: Structure
    """
    parser = _SelectParser()
    parse_results = parser(query)[0]
    indicator = _AtomIndicator(parse_results=parse_results)
    pruner = _StructurePruner(indicator)
    result_structure = structure.copy()

    pruner(result_structure)

    return result_structure


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
