# Copyright 2017 by Peter Cock.  All rights reserved.
# Based on code Copyright 2002 by Thomas Hamelryck.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

from Bio.PDB import NeighborSearch
from numpy.random import random
import unittest


class Atom(object):
    def __init__(self):
        self.coord = (100 * random(3))

    def get_coord(self):
        return self.coord


class NeighborSearchTest(unittest.TestCase):

    def test_NeighborSearch(self):
        for i in range(0, 20):
            al = [Atom() for j in range(100)]
            ns = NeighborSearch(al)
            for i in [0, 1, 2, 3, 4, 5, 6]:
                self.assertEqual(i, len(ns.search_all(5.0)))
