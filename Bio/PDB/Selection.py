# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from types import ListType

from PDBExceptions import PDBException

__doc__="Selection of atoms, residues, etc."


entity_levels=["A", "R", "C", "M", "S"]

def uniqueify(l):
    "Return unique items in list l."
    d={}
    for i in l:
        if not d.has_key(i):
            d[i]=None
    return d.keys()

def get_unique_parents(entity_list):
    """
    Translate a list of entities to a list of their
    (unique) parents.
    """ 
    l=[]
    for entity in entity_list:
        parent=entity.get_parent()
        l.append(parent)
    return uniqueify(l)

def unfold_entities(entity_list, target_level):
    """
    Unfold a list of entities to a list of entities of another 
    level.  E.g.:

    list of atoms -> list of residues
    list of modules -> list of atoms
    list of residues -> list of chains

    o entity_list - list of entities or a single entity
    o target_level - char (A, R, C, M, S)
    """
    if not target_level in entity_levels:
        raise PDBException("%s: Not an entity level." % target_level)
    if type(entity_list)!=ListType:
        # single entity
        entity_list=[entity_list]
    # level of entity list
    level=entity_list[0].get_level()
    for entity in entity_list:
        if not (entity.get_level()==level):
            raise PDBException("Entity list is not homogeneous.")
    target_index=entity_levels.index(target_level)
    level_index=entity_levels.index(level)
    if level_index==target_index:
        # already right level
        return entity_list
    if level_index>target_index:
        # we're going down, e.g. S->A
        for i in range(target_index, level_index):
            new_entity_list=[]
            for entity in entity_list:
                new_entity_list=new_entity_list+entity.get_list()
            entity_list=new_entity_list
    else:
        # we're going up, e.g. A->S
        for i in range(level_index, target_index):
            new_entity_list=[]  
            for entity in entity_list:
                parent=entity.get_parent()
                new_entity_list.append(parent)
            # find unique parents
            entity_list=uniqueify(new_entity_list)
    return entity_list

