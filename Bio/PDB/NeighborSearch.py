# Copyright (C) 2002, 2004 Thomas Hamelryck (thamelry@binf.ku.dk)
# All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Fast atom neighbor lookup using a KD tree (implemented in C)."""


import numpy

from Bio.PDB.PDBExceptions import PDBException
from Bio.PDB.Selection import unfold_entities, entity_levels, uniqueify


class NeighborSearch:
    """Class for neighbor searching.

    This class can be used for two related purposes:

     1. To find all atoms/residues/chains/models/structures within radius
        of a given query position.
     2. To find all atoms/residues/chains/models/structures that are within
        a fixed radius of each other.

    NeighborSearch makes use of the KDTree class implemented in C for speed.
    """

    def __init__(self, atom_list, bucket_size=10):
        """Create the object.

        Arguments:
         - atom_list - list of atoms. This list is used in the queries.
           It can contain atoms from different structures.
         - bucket_size - bucket size of KD tree. You can play around
           with this to optimize speed if you feel like it.

        """
        from Bio.PDB.kdtrees import KDTree

        self.atom_list = atom_list
        # get the coordinates
        coord_list = [a.get_coord() for a in atom_list]
        # to Nx3 array of type float
        self.coords = numpy.array(coord_list, dtype="d")
        assert bucket_size > 1
        assert self.coords.shape[1] == 3
        self.kdt = KDTree(self.coords, bucket_size)

    # Private

    def _get_unique_parent_pairs(self, pair_list):
        # translate a list of (entity, entity) tuples to
        # a list of (parent entity, parent entity) tuples,
        # thereby removing duplicate (parent entity, parent entity)
        # pairs.
        # o pair_list - a list of (entity, entity) tuples
        parent_pair_list = []
        for (e1, e2) in pair_list:
            p1 = e1.get_parent()
            p2 = e2.get_parent()
            if p1 == p2:
                continue
            elif p1 < p2:
                parent_pair_list.append((p1, p2))
            else:
                parent_pair_list.append((p2, p1))
        return uniqueify(parent_pair_list)

    # Public

    def search(self, center, radius, level="A"):
        """Neighbor search.

        Return all atoms/residues/chains/models/structures
        that have at least one atom within radius of center.
        What entity level is returned (e.g. atoms or residues)
        is determined by level (A=atoms, R=residues, C=chains,
        M=models, S=structures).

        Arguments:
         - center - Numeric array
         - radius - float
         - level - char (A, R, C, M, S)

        """
        if level not in entity_levels:
            raise PDBException(f"{level}: Unknown level")
        center = numpy.require(center, dtype="d", requirements="C")
        if center.shape != (3,):
            raise Exception("Expected a 3-dimensional NumPy array")
        points = self.kdt.search(center, radius)
        atom_list = [self.atom_list[point.index] for point in points]
        if level == "A":
            return atom_list
        else:
            return unfold_entities(atom_list, level)

    def search_all(self, radius, level="A"):
        """All neighbor search.

        Search all entities that have atoms pairs within
        radius.

        Arguments:
         - radius - float
         - level - char (A, R, C, M, S)

        """
        if level not in entity_levels:
            raise PDBException(f"{level}: Unknown level")
        neighbors = self.kdt.neighbor_search(radius)
        atom_list = self.atom_list
        atom_pair_list = []
        for neighbor in neighbors:
            i1 = neighbor.index1
            i2 = neighbor.index2
            a1 = atom_list[i1]
            a2 = atom_list[i2]
            atom_pair_list.append((a1, a2))
        if level == "A":
            # return atoms
            return atom_pair_list
        next_level_pair_list = atom_pair_list
        for next_level in ["R", "C", "M", "S"]:
            next_level_pair_list = self._get_unique_parent_pairs(next_level_pair_list)
            if level == next_level:
                return next_level_pair_list
