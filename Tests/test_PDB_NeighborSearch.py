from Bio.PDB import NeighborSearch

from numpy.random import random

import unittest


class Atom(object):
    def __init__(self):
        self.coord = (100 * random(3))

    def get_coord(self):
        return self.coord


class KDTreeTest(unittest.TestCase):
    
    def test_NeighborSearch():
        for i in range(0, 20):
            al = [Atom() for j in range(100)]
            ns = NeighborSearch(al)
            for i in range(0, 6):
                self.assertEqual(i, len(ns.search_all(5.0)))
