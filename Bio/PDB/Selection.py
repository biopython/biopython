# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Selection of atoms, residues, etc."""


import itertools

from Bio.PDB.Atom import Atom
from Bio.PDB.Entity import Entity
from Bio.PDB.PDBExceptions import PDBException


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
        raise PDBException("%s: Not an entity level." % target_level)
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

    if level_index > target_index:  # we're going down, e.g. S->A
        for i in range(target_index, level_index):
            entity_list = itertools.chain.from_iterable(entity_list)
    else:  # we're going up, e.g. A->S
        for i in range(level_index, target_index):
            # find unique parents
            entity_list = {entity.get_parent() for entity in entity_list}
    return list(entity_list)


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
