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
    from numpy import array
    from numpy.random import random
except ImportError:
    from Bio import MissingExternalDependencyError
    raise MissingExternalDependencyError(\
        "Install NumPy if you want to use Bio.PDB.")

try:
    from Bio.KDTree import _CKDTree
except ImportError:
    from Bio import MissingExternalDependencyError
    raise MissingExternalDependencyError(\
        "C module in Bio.KDTree not compiled")

from Bio.PDB.NeighborSearch import NeighborSearch

class NeighborTest(unittest.TestCase):
    def test_neighbor_search(self):
        """NeighborSearch: Find nearby randomly generated coordinates.
         
        Based on the self test in Bio.PDB.NeighborSearch.
        """
        class RandomAtom:
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
        x = array([250,250,250]) #Far away from our random atoms
        self.assertEqual([], ns.search(x, 5.0, "A"))
        self.assertEqual([], ns.search(x, 5.0, "R"))
        self.assertEqual([], ns.search(x, 5.0, "C"))
        self.assertEqual([], ns.search(x, 5.0, "M"))
        self.assertEqual([], ns.search(x, 5.0, "S"))


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
