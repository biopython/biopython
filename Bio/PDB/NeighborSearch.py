# Copyright (C) 2002, 2004 Thomas Hamelryck (thamelry@binf.ku.dk)
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


"""Fast atom neighbor lookup using a KD tree (implemented in C)."""

from __future__ import print_function

import numpy

from Bio.PDB.PDBExceptions import PDBException
from Bio.PDB.Selection import unfold_entities, entity_levels, uniqueify

from . import _kdtrees


class KDTree(_kdtrees.KDTree):
    """KD tree data structure for searching N-dimensional vectors.

    The KD tree data structure can be used for all kinds of searches that
    involve N-dimensional vectors. For example, neighbor searches (find all
    points within a radius of a given point) or finding all point pairs in a
    set that are within a certain radius of each other.

    Reference:

    Computational Geometry: Algorithms and Applications
    Second Edition
    Mark de Berg, Marc van Kreveld, Mark Overmars, Otfried Schwarzkopf
    published by Springer-Verlag
    2nd rev. ed. 2000.
    ISBN: 3-540-65620-0

    The KD tree data structure is described in chapter 5, pg. 99.

    The following article made clear to me that the nodes should
    contain more than one point (this leads to dramatic speed
    improvements for the "all fixed radius neighbor search", see
    below):

    JL Bentley, "K-d trees for semidynamic point sets," in Sixth Annual ACM
    Symposium on Computational Geometry, vol. 91. San Francisco, 1990

    This KD implementation also performs a "all fixed radius neighbor search",
    i.e. it can find all point pairs in a set that are within a certain radius
    of each other. As far as I know the algorithm has not been published.
    """

    def __init__(self, dim, bucket_size=1):
        """Initialize KDTree class."""
        self.dim = dim
        _kdtrees.KDTree.__init__(self, dim, bucket_size)
        self.built = 0

    # Set data

    def set_coords(self, coords):
        """Add the coordinates of the points.

        Arguments:
         - coords: two dimensional NumPy array. E.g. if the points
           have dimensionality D and there are N points, the coords
           array should be NxD dimensional.

        """
        if coords.min() <= -1e6 or coords.max() >= 1e6:
            raise Exception("Points should lie between -1e6 and 1e6")
        if len(coords.shape) != 2 or coords.shape[1] != self.dim:
            raise Exception("Expected a Nx%i NumPy array" % self.dim)
        self.set_data(coords)
        self.built = 1

    # Fixed radius search for a point

    def search(self, center, radius):
        """Search all points within radius of center.

        Arguments:
         - center: one dimensional NumPy array. E.g. if the points have
           dimensionality D, the center array should be D dimensional.
         - radius: float>0

        """
        if not self.built:
            raise Exception("No point set specified")
        if center.shape != (self.dim,):
            raise Exception("Expected a %i-dimensional NumPy array"
                            % self.dim)
        self.search_center_radius(center, radius)

    def get_radii(self):
        """Return radii.

        Return the list of distances from center after
        a neighbor search.
        """
        n = self.get_count()
        if n == 0:
            return []
        radii = numpy.empty(n, int)
        _kdtrees.KDTree.get_radii(self, radii)
        return radii

    def get_indices(self):
        """Return the list of indices.

        Return the list of indices after a neighbor search.
        The indices refer to the original coords NumPy array. The
        coordinates with these indices were within radius of center.

        For an index pair, the first index<second index.
        """
        n = self.get_count()
        if n == 0:
            return []
        indices = numpy.empty(n, int)
        _kdtrees.KDTree.get_indices(self, indices)
        return indices

    def all_search(self, radius):
        """All fixed neighbor search.

        Search all point pairs that are within radius.

        Arguments:
         - radius: float (>0)

        """
        # Fixed radius search for all points
        if not self.built:
            raise Exception("No point set specified")
        self.neighbors = self.neighbor_search(radius)

    def all_get_indices(self):
        """Return All Fixed Neighbor Search results.

        Return a Nx2 dim NumPy array containing
        the indices of the point pairs, where N
        is the number of neighbor pairs.
        """
        a = numpy.array([[neighbor.index1, neighbor.index2] for neighbor in self.neighbors])
        return a

    def all_get_radii(self):
        """Return All Fixed Neighbor Search results.

        Return an N-dim array containing the distances
        of all the point pairs, where N is the number
        of neighbor pairs..
        """
        return [neighbor.radius for neighbor in self.neighbors]


class NeighborSearch(object):
    """Class for neighbor searching.

    This class can be used for two related purposes:

     1. To find all atoms/residues/chains/models/structures within radius
        of a given query position.
     2. To find all atoms/residues/chains/models/structures that are within
        a fixed radius of each other.

    NeighborSearch makes use of the Bio.KDTree C++ module, so it's fast.
    """

    def __init__(self, atom_list, bucket_size=10):
        """Create the object.

        Arguments:
         - atom_list - list of atoms. This list is used in the queries.
           It can contain atoms from different structures.
         - bucket_size - bucket size of KD tree. You can play around
           with this to optimize speed if you feel like it.

        """
        self.atom_list = atom_list
        # get the coordinates
        coord_list = [a.get_coord() for a in atom_list]
        # to Nx3 array of type float
        self.coords = numpy.array(coord_list).astype("f")
        assert(bucket_size > 1)
        assert(self.coords.shape[1] == 3)
        self.kdt = KDTree(3, bucket_size)
        self.kdt.set_coords(self.coords)

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
            raise PDBException("%s: Unknown level" % level)
        self.kdt.search(center, radius)
        indices = self.kdt.get_indices()
        n_atom_list = []
        atom_list = self.atom_list
        for i in indices:
            a = atom_list[i]
            n_atom_list.append(a)
        if level == "A":
            return n_atom_list
        else:
            return unfold_entities(n_atom_list, level)

    def search_all(self, radius, level="A"):
        """All neighbor search.

        Search all entities that have atoms pairs within
        radius.

        Arguments:
         - radius - float
         - level - char (A, R, C, M, S)

        """
        if level not in entity_levels:
            raise PDBException("%s: Unknown level" % level)
        self.kdt.all_search(radius)
        indices = self.kdt.all_get_indices()
        atom_list = self.atom_list
        atom_pair_list = []
        for i1, i2 in indices:
            a1 = atom_list[i1]
            a2 = atom_list[i2]
            atom_pair_list.append((a1, a2))
        if level == "A":
            # return atoms
            return atom_pair_list
        next_level_pair_list = atom_pair_list
        for l in ["R", "C", "M", "S"]:
            next_level_pair_list = self._get_unique_parent_pairs(next_level_pair_list)
            if level == l:
                return next_level_pair_list


if __name__ == "__main__":

    from numpy.random import random

    class Atom(object):
        def __init__(self):
            """Initialize the class."""
            self.coord = (100 * random(3))

        def get_coord(self):
            return self.coord

    for i in range(0, 20):
        # Make a list of 100 atoms
        al = [Atom() for j in range(100)]
        ns = NeighborSearch(al)
        print("Found %i" % len(ns.search_all(5.0)))
