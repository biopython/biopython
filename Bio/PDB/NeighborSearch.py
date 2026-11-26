# Copyright (C) 2002, 2004 Thomas Hamelryck (thamelry@binf.ku.dk)
# All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Fast atom neighbor lookup using a KD tree (implemented in C)."""

import numpy as np

from Bio.PDB.PDBExceptions import PDBException
from Bio.PDB.Selection import entity_levels
from Bio.PDB.Selection import unfold_entities


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
        self.coords = np.array(coord_list, dtype="d")
        assert bucket_size > 1
        assert self.coords.shape[1] == 3
        self.kdt = KDTree(self.coords, bucket_size)

    # Public

    def search(self, center, radius, level="A"):
        """Neighbor search.

        Return all atoms/residues/chains/models/structures
        that have at least one atom within radius of center.
        What entity level is returned (e.g. atoms or residues)
        is determined by level (A=atoms, R=residues, C=chains,
        M=models, S=structures).

        Arguments:
         - center - NumPy array
         - radius - float
         - level - char (A, R, C, M, S)

        """
        if level not in entity_levels:
            raise PDBException(f"{level}: Unknown level")
        center = np.require(center, dtype="d", requirements="C")
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

        if level == "A":
            # return atoms
            return [(atom_list[n.index1], atom_list[n.index2]) for n in neighbors]
        # Atom indexes as pairs (N,2) for vectorised mapping
        pairs = np.fromiter(
            (i for n in neighbors for i in (n.index1, n.index2)),
            dtype=np.intp,
            count=2 * len(neighbors),
        ).reshape(-1, 2)

        depth = "ARCMS".index(level)
        entities, seen = [], {}
        parent_idx = np.empty(len(self.atom_list), dtype=np.intp)
        for k, atom in enumerate(self.atom_list):
            parent = atom
            # Map each atom to its parent at the requested level
            # Computed once per atom instead of once per pair
            for _ in range(depth):
                parent = parent.get_parent()
            # Assign a new index the first time the parent is seen
            key = id(parent)
            if key not in seen:
                seen[key] = len(entities)
                # List of all atom parents
                entities.append(parent)
            # Mapping each atom to its parents index
            parent_idx[k] = seen[key]

        # Translate atom pairs into parent pairs
        p = parent_idx[pairs]
        # Drop pairs where both atoms belong to the same parent
        p = p[p[:, 1] != p[:, 0]]
        # Order each pair so (a, b) and (b, a) become identical, then deduplicate
        p.sort(axis=1)
        unique_pairs = np.unique(p, axis=0)

        return [(entities[i], entities[j]) for i, j in unique_pairs.tolist()]
