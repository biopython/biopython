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

from Bio.KDTree.KDTree import _neighbor_test, _test

nr_points = 5000
dim = 3
bucket_size = 5
radius = 0.01


class KDTreeTest(unittest.TestCase):

    def test_KDTree_neighbour(self):
        for i in range(0, 10):
            self.assertTrue(_neighbor_test(nr_points, dim, bucket_size, radius))

    def test_KDTree(self):
        for i in range(0, 10):
            self.assertTrue(_test(nr_points, dim, bucket_size, radius))


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
