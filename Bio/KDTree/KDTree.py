# Copyright 2004 by Thomas Hamelryck.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""KD tree data structure for searching N-dimensional vectors.

The KD tree data structure can be used for all kinds of searches that
involve N-dimensional vectors, e.g.  neighbor searches (find all points
within a radius of a given point) or finding all point pairs in a set
that are within a certain radius of each other. See "Computational Geometry:
Algorithms and Applications" (Mark de Berg, Marc van Kreveld, Mark Overmars,
Otfried Schwarzkopf). Author: Thomas Hamelryck.
"""

# from __future__ import print_function
from numpy import array, empty
from Bio.KDTree import _CKDTree


class KDTree(object):
    """KD tree implementation in C++, SWIG python wrapper.

    The KD tree data structure can be used for all kinds of searches that
    involve N-dimensional vectors, e.g.  neighbor searches (find all points
    within a radius of a given point) or finding all point pairs in a set
    that are within a certain radius of each other.

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

    JL Bentley, "Kd trees for semidynamic point sets," in Sixth Annual ACM
    Symposium on Computational Geometry, vol. 91. San Francisco, 1990

    This KD implementation also performs a "all fixed radius neighbor search",
    i.e. it can find all point pairs in a set that are within a certain radius
    of each other. As far as I know the algorithm has not been published.
    """

    def __init__(self, dim, bucket_size=1):
        """Initialize KDTree class."""
        self.dim = dim
        self.kdt = _CKDTree.KDTree(dim, bucket_size)
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
        self.kdt.set_data(coords)
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
        self.kdt.search_center_radius(center, radius)

    def get_radii(self):
        """Return radii.

        Return the list of distances from center after
        a neighbor search.
        """
        n = self.kdt.get_count()
        if n == 0:
            return []
        radii = empty(n, int)
        self.kdt.get_radii(radii)
        return radii

    def get_indices(self):
        """Return the list of indices.

        Return the list of indices after a neighbor search.
        The indices refer to the original coords NumPy array. The
        coordinates with these indices were within radius of center.

        For an index pair, the first index<second index.
        """
        n = self.kdt.get_count()
        if n == 0:
            return []
        indices = empty(n, int)
        self.kdt.get_indices(indices)
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
        self.neighbors = self.kdt.neighbor_search(radius)

    def all_get_indices(self):
        """Return All Fixed Neighbor Search results.

        Return a Nx2 dim NumPy array containing
        the indices of the point pairs, where N
        is the number of neighbor pairs.
        """
        a = array([[neighbor.index1, neighbor.index2] for neighbor in self.neighbors])
        return a

    def all_get_radii(self):
        """Return All Fixed Neighbor Search results.

        Return an N-dim array containing the distances
        of all the point pairs, where N is the number
        of neighbor pairs..
        """
        return [neighbor.radius for neighbor in self.neighbors]
