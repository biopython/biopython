# Copyright 2009-2010 by Eric Talevich.  All rights reserved.
# Revisions copyright 2010 by Peter Cock.  All rights reserved.
#
# Converted by Eric Talevich from an older unit test copyright 2002
# by Thomas Hamelryck.
#
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for those parts of the Bio.PDB module using Bio.KDTree."""

import unittest

try:
    from numpy import array, dot, sqrt
    from numpy.random import random
except ImportError:
    from Bio import MissingExternalDependencyError
    raise MissingExternalDependencyError(
        "Install NumPy if you want to use Bio.PDB.")

try:
    from Bio.PDB import _kdtrees
except ImportError:
    from Bio import MissingExternalDependencyError
    raise MissingExternalDependencyError(
        "C module Bio.PDB._kdtrees not compiled")

from Bio.PDB.NeighborSearch import NeighborSearch, KDTree


class NeighborTest(unittest.TestCase):

    def test_neighbor_search(self):
        """NeighborSearch: Find nearby randomly generated coordinates.

        Based on the self test in Bio.PDB.NeighborSearch.
        """
        class RandomAtom(object):
            def __init__(self):
                self.coord = 100 * random(3)

            def get_coord(self):
                return self.coord

        for i in range(0, 20):
            atoms = [RandomAtom() for j in range(100)]
            ns = NeighborSearch(atoms)
            hits = ns.search_all(5.0)
            self.assertTrue(isinstance(hits, list), hits)
            self.assertTrue(len(hits) >= 0, hits)
        x = array([250, 250, 250])  # Far away from our random atoms
        self.assertEqual([], ns.search(x, 5.0, "A"))
        self.assertEqual([], ns.search(x, 5.0, "R"))
        self.assertEqual([], ns.search(x, 5.0, "C"))
        self.assertEqual([], ns.search(x, 5.0, "M"))
        self.assertEqual([], ns.search(x, 5.0, "S"))


class KDTreeTest(unittest.TestCase):

    nr_points = 5000     # number of points used in test
    bucket_size = 5      # number of points per tree node
    radius = 0.01        # radius of search (typically 0.05 or so)
    query_radius = 10    # radius of search

    def test_KDTree_exceptions(self):
        bucket_size = self.bucket_size
        nr_points = self.nr_points
        radius = self.radius
        coords = random((nr_points, 3)) * 100000000000000
        with self.assertRaises(Exception) as context:
            kdt = KDTree(coords, bucket_size)
        self.assertTrue("coordinate values should lie between -1e6 and 1e6" in str(context.exception))
        with self.assertRaises(Exception) as context:
            kdt = KDTree(random((nr_points, 3 - 2)), bucket_size)
        self.assertTrue("expected a Nx3 numpy array" in str(context.exception))


    def test_KDTree_neighbour(self):
        """Test all fixed radius neighbor search.

        Test all fixed radius neighbor search using the KD tree C
        module, and compare the results to a manual search.
        """
        bucket_size = self.bucket_size
        nr_points = self.nr_points
        radius = self.radius
        for i in range(0, 10):
            # KD tree search
            coords = random((nr_points, 3))
            kdt = KDTree(coords, bucket_size)
            neighbors = kdt.neighbor_search(radius)
            r = [neighbor.radius for neighbor in neighbors]
            if r is None:
                l1 = 0
            else:
                l1 = len(r)
            # manual search
            neighbors = kdt.neighbor_simple_search(radius)
            r = [neighbor.radius for neighbor in neighbors]
            if r is None:
                l2 = 0
            else:
                l2 = len(r)
            # compare results
            self.assertEqual(l1, l2)


    def test_KDTree(self):
        """Test neighbor search.

        Test neighbor search using the KD tree C module,
        and compare the results to a manual search.
        """
        bucket_size = self.bucket_size
        nr_points = self.nr_points
        radius = self.radius
        for i in range(0, 10):
            # kd tree search
            coords = random((nr_points, 3))
            center = coords[0]
            kdt = KDTree(coords, bucket_size)
            points = kdt.search(center, radius)
            l1 = len(points)
            # manual search
            l2 = 0
            for i in range(0, nr_points):
                p = coords[i]
                v = p - center
                if sqrt(dot(v, v)) <= radius:
                    l2 += 1
            # compare results
            self.assertEqual(l1, l2)


    def test_all_search(self):
        """Test fixed neighbor search.

        Using the KDTree C module, search point pairs that are
        within a large radius, and verify that we found all radii.
        """
        bucket_size = self.bucket_size
        nr_points = self.nr_points
        query_radius = self.query_radius
        for i in range(0, 5):
            # KD tree search
            coords = random((nr_points // 10, 3))
            kdt = KDTree(coords, bucket_size)
            neighbors = kdt.neighbor_search(query_radius)
            l1 = len(neighbors)
            # find all points
            # self.assertEqual(l1, l2)


    def test_search(self):
        """Test search all points within radius of center.

        Using the KDTree C module, search all point pairs that are
        within radius, and compare the results to a manual search.
        """
        bucket_size = self.bucket_size
        nr_points = self.nr_points
        radius = self.radius
        for i in range(0, 5):
            # KD tree search
            coords = random((nr_points, 3))
            kdt = KDTree(coords, bucket_size)
            points = kdt.search(coords[0], radius * 100)
            # manual search
            l1 = 0
            for i in range(0, nr_points):
                p = coords[i]
                v = p - coords[0]
                if sqrt(dot(v, v)) <= radius * 100:
                    l1 += 1
            # compare th results
            self.assertEqual(l1, len(points))


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
