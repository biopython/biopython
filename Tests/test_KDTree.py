# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import unittest

try:
    import numpy
    del numpy
except ImportError:
    from Bio import MissingExternalDependencyError
    raise MissingExternalDependencyError(
        "Install NumPy if you want to use Bio.KDTree.")

try:
    from Bio.KDTree import _CKDTree
    del _CKDTree
except ImportError:
    from Bio import MissingExternalDependencyError
    raise MissingExternalDependencyError(
            "C module in Bio.KDTree not compiled")

from Bio.KDTree.KDTree import KDTree
from numpy import sum, sqrt, array
from numpy import random

nr_points = 5000
dim = 3
bucket_size = 5
radius = 0.01
query_radius = 10


def _dist(p, q):
    diff = p - q
    return sqrt(sum(diff * diff))


def neighbor_test(nr_points, dim, bucket_size, radius):
    """Test all fixed radius neighbor search.

    Test all fixed radius neighbor search using the
    KD tree C module.

    Arguments:
     - nr_points: number of points used in test
     - dim: dimension of coords
     - bucket_size: nr of points per tree node
     - radius: radius of search (typically 0.05 or so)

    Returns true if the test passes.
    """
    # KD tree search
    kdt = KDTree(dim, bucket_size)
    coords = random.random((nr_points, dim))
    kdt.kdt.set_data(coords)
    neighbors = kdt.kdt.neighbor_search(radius)
    r = [neighbor.radius for neighbor in neighbors]
    if r is None:
        l1 = 0
    else:
        l1 = len(r)
    # now do a slow search to compare results
    neighbors = kdt.kdt.neighbor_simple_search(radius)
    r = [neighbor.radius for neighbor in neighbors]
    if r is None:
        l2 = 0
    else:
        l2 = len(r)
    if l1 == l2:
        # print("Passed.")
        return True
    else:
        print("Not passed: %i != %i." % (l1, l2))
        return False


def test(nr_points, dim, bucket_size, radius):
    """Test neighbor search.

    Test neighbor search using the KD tree C module.

    Arguments:
     - nr_points: number of points used in test
     - dim: dimension of coords
     - bucket_size: nr of points per tree node
     - radius: radius of search (typically 0.05 or so)

    Returns true if the test passes.
    """
    # kd tree search
    kdt = KDTree(dim, bucket_size)
    coords = random.random((nr_points, dim))
    center = coords[0]
    kdt.kdt.set_data(coords)
    kdt.kdt.search_center_radius(center, radius)
    r = kdt.get_indices()
    if r is None:
        l1 = 0
    else:
        l1 = len(r)
    l2 = 0
    # now do a manual search to compare results
    for i in range(0, nr_points):
        p = coords[i]
        if _dist(p, center) <= radius:
            l2 = l2 + 1
    if l1 == l2:
        # print("Passed.")
        return True
    else:
        print("Not passed: %i != %i." % (l1, l2))
        return False


def test_search(nr_points, dim, bucket_size, radius):
    """Test search all points within radius of center.

    Search all point pairs that are within radius.

    Arguments:
     - nr_points: number of points used in test
     - dim: dimension of coords
     - bucket_size: nr of points per tree node
     - radius: radius of search

    Returns true if the test passes.
    """
    kdt = KDTree(dim, bucket_size)
    coords = random.random((nr_points, dim))
    kdt.set_coords(coords)
    kdt.search(coords[0], radius * 100)
    radii = kdt.get_radii()
    l1 = 0
    for i in range(0, nr_points):
        p = coords[i]
        if _dist(p, coords[0]) <= radius * 100:
            l1 = l1 + 1
    if l1 == len(radii):
        return True
    else:
        return False


def test_all_search(nr_points, dim, bucket_size, query_radius):
    """Test fixed neighbor search.

    Search all point pairs that are within radius.

    Arguments:
     - nr_points: number of points used in test
     - dim: dimension of coords
     - bucket_size: nr of points per tree node
     - query_radius: radius of search

    Returns true if the test passes.
    """
    kdt = KDTree(dim, bucket_size)
    coords = random.random((nr_points, dim))
    kdt.set_coords(coords)
    kdt.all_search(query_radius)
    indices = kdt.all_get_indices()
    if indices is None:
        l1 = 0
    else:
        l1 = len(indices)
    radii = kdt.all_get_radii()
    if radii is None:
        l2 = 0
    else:
        l2 = len(radii)
    if l1 == l2:
        return True
    else:
        return False


class KDTreeTest(unittest.TestCase):

    def test_KDTree_exceptions(self):
        kdt = KDTree(dim, bucket_size)
        with self.assertRaises(Exception) as context:
            kdt.set_coords(random.random((nr_points, dim)) * 100000000000000)
        self.assertTrue("Points should lie between -1e6 and 1e6" in str(context.exception))
        with self.assertRaises(Exception) as context:
            kdt.set_coords(random.random((nr_points, dim - 2)))
        self.assertTrue("Expected a Nx%i NumPy array" % dim in str(context.exception))
        with self.assertRaises(Exception) as context:
            kdt.search(array([0, 0, 0]), radius)
        self.assertTrue("No point set specified" in str(context.exception))

    def test_KDTree_neighbour(self):
        for i in range(0, 10):
            self.assertTrue(neighbor_test(nr_points, dim, bucket_size, radius))

    def test_KDTree(self):
        for i in range(0, 10):
            self.assertTrue(test(nr_points, dim, bucket_size, radius))

    def test_all_search(self):
        for i in range(0, 5):
            self.assertTrue(test_all_search((nr_points // 10), dim, bucket_size, query_radius))

    def test_search(self):
        for i in range(0, 5):
            self.assertTrue(test_search(nr_points, dim, bucket_size, radius))


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
