# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Fast atom neighbor lookup using a KD tree (implemented in C++)."""

import numpy

from Bio.KDTree import KDTree

from Bio.PDB.PDBExceptions import PDBException
from Bio.PDB.Selection import unfold_entities, entity_levels, uniqueify


class NeighborSearch(object):
    """
    This class can be used for two related purposes:

    1. To find all atoms/residues/chains/models/structures within radius 
    of a given query position. 

    2. To find all atoms/residues/chains/models/structures that are within 
    a fixed radius of each other.

    NeighborSearch makes use of the Bio.KDTree C++ module, so it's fast.
    """
    def __init__(self, atom_list, bucket_size=10):
        """
        o atom_list - list of atoms. This list is used in the queries.
        It can contain atoms from different structures.
        o bucket_size - bucket size of KD tree. You can play around 
        with this to optimize speed if you feel like it.
        """
        self.atom_list=atom_list
        # get the coordinates
        coord_list = [a.get_coord() for a in atom_list]
        # to Nx3 array of type float
        self.coords=numpy.array(coord_list).astype("f")
        assert(bucket_size>1)
        assert(self.coords.shape[1]==3)
        self.kdt=KDTree(3, bucket_size)
        self.kdt.set_coords(self.coords)

    # Private

    def _get_unique_parent_pairs(self, pair_list):
        # translate a list of (entity, entity) tuples to 
        # a list of (parent entity, parent entity) tuples,
        # thereby removing duplicate (parent entity, parent entity)
        # pairs.
        # o pair_list - a list of (entity, entity) tuples
        parent_pair_list=[]
        for (e1, e2) in pair_list:
            p1=e1.get_parent()
            p2=e2.get_parent()
            if p1==p2:
                continue
            elif p1<p2:
                parent_pair_list.append((p1, p2))
            else:
                parent_pair_list.append((p2, p1))
        return uniqueify(parent_pair_list)

    # Public

    def search(self, center, radius, level="A"):
        """Neighbor search.

        Return all atoms/residues/chains/models/structures
        that have at least one atom within radius of center.
        What entitity level is returned (e.g. atoms or residues)
        is determined by level (A=atoms, R=residues, C=chains,
        M=models, S=structures).

        o center - Numeric array 
        o radius - float
        o level - char (A, R, C, M, S)
        """
        if not level in entity_levels:
            raise PDBException("%s: Unknown level" % level)
        self.kdt.search(center, radius)
        indices=self.kdt.get_indices()
        n_atom_list=[]
        atom_list=self.atom_list
        for i in indices:
            a=atom_list[i]
            n_atom_list.append(a)
        if level=="A":
            return n_atom_list
        else:
            return unfold_entities(n_atom_list, level)
            
    def search_all(self, radius, level="A"):
        """All neighbor search.

        Search all entities that have atoms pairs within
        radius. 

        o radius - float
        o level - char (A, R, C, M, S)
        """
        if not level in entity_levels:
            raise PDBException("%s: Unknown level" % level)
        self.kdt.all_search(radius)
        indices=self.kdt.all_get_indices()
        atom_list=self.atom_list
        atom_pair_list=[]
        for i1, i2 in indices:
            a1=atom_list[i1]
            a2=atom_list[i2]
            atom_pair_list.append((a1, a2))
        if level=="A":
            # return atoms
            return atom_pair_list
        next_level_pair_list=atom_pair_list
        for l in ["R", "C", "M", "S"]:
            next_level_pair_list=self._get_unique_parent_pairs(next_level_pair_list)
            if level==l:
                return next_level_pair_list 

if __name__=="__main__":

    from numpy.random import random

    class Atom(object):
        def __init__(self):
            self.coord=(100*random(3))

        def get_coord(self):
            return self.coord

    for i in range(0, 20):
        #Make a list of 100 atoms
        al = [Atom() for j in range(100)]

        ns=NeighborSearch(al)

        print "Found ", len(ns.search_all(5.0))

