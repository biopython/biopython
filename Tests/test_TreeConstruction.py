# Copyright (C) 2013 by Yanbo Ye (yeyanbo289@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for the Bio.TreeConstruction module."""

import unittest

from Bio import AlignIO
from Bio.Phylo import TreeConstruction
from Bio.Phylo.TreeConstruction import DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceCaluculator

class DistanceMatrixTest(unittest.TestCase):
    """Test for DistanceMatrix construction and manipulation"""

    def test_good_construction(self):
        names = ['Alpha', 'Beta', 'Gamma', 'Delta']
        matrix = [[0], [1, 0], [2, 3, 0], [4, 5, 6, 0]]
        dm = DistanceMatrix(names, matrix)
        self.assertTrue(isinstance(dm, TreeConstruction.DistanceMatrix))
        self.assertEqual(dm.names[0], 'Alpha')
        self.assertEqual(dm.matrix[2][1], 3)
        self.assertEqual(len(dm), 4)

    def test_bad_construction(self):
        self.assertRaises(TypeError, DistanceMatrix, ['Alpha', 100, 'Gamma', 'Delta'], [[0], [0.1, 0], [0.2, 0.3, 0], [0.4, 0.5, 0.6, 0]])
        self.assertRaises(TypeError, DistanceMatrix, ['Alpha', 'Beta', 'Gamma', 'Delta'], [[0], ['a'], [0.2, 0.3], [0.4, 0.5, 0.6]])
        self.assertRaises(ValueError, DistanceMatrix, ['Alpha', 'Alpha', 'Gamma', 'Delta'], [[0], [0.1], [0.2, 0.3], [0.4, 0.5, 0.6]])
        self.assertRaises(ValueError, DistanceMatrix, ['Alpha', 'Beta', 'Gamma', 'Delta'], [[0], [0.2, 0], [0.4, 0.5, 0.6]])
        self.assertRaises(ValueError, DistanceMatrix, ['Alpha', 'Beta', 'Gamma', 'Delta'], [[0], [0.1], [0.2, 0.3, 0.4], [0.4, 0.5, 0.6]])

    def test_good_manipulation(self):
        names = ['Alpha', 'Beta', 'Gamma', 'Delta']
        matrix = [[0], [1, 0], [2, 3, 0], [4, 5, 6, 0]]
        dm = DistanceMatrix(names, matrix)
        # getitem
        self.assertEqual(dm[1], [1, 0, 3, 5])
        self.assertEqual(dm[2, 1], 3)
        self.assertEqual(dm['Alpha'], [0, 1, 2, 4])
        self.assertEqual(dm['Gamma', 'Delta'], 6)
        # setitem
        dm['Alpha'] = [0, 10, 20, 40]
        self.assertEqual(dm['Alpha'], [0, 10, 20, 40])
        # delitem insert item
        del dm[1]
        self.assertEqual(dm.names, ['Alpha', 'Gamma', 'Delta'])
        self.assertEqual(dm.matrix, [[0], [20, 0], [40, 6, 0]])
        dm.insert('Beta', [1, 0, 3, 5], 1)
        self.assertEqual(dm.names, names)
        self.assertEqual(dm.matrix, [[0], [1, 0], [20, 3, 0], [40, 5, 6, 0]])
        del dm['Alpha']
        self.assertEqual(dm.names, ['Beta', 'Gamma', 'Delta'])
        self.assertEqual(dm.matrix, [[0], [3, 0], [5, 6, 0]])
        dm.insert('Alpha', [1, 2, 4, 0])
        self.assertEqual(dm.names, ['Beta', 'Gamma', 'Delta', 'Alpha'])
        self.assertEqual(dm.matrix, [[0], [3, 0], [5, 6, 0], [1, 2, 4, 0]])

    def test_bad_manipulation(self):
        names = ['Alpha', 'Beta', 'Gamma', 'Delta']
        matrix = [[0], [1, 0], [2, 3, 0], [4, 5, 6, 0]]
        dm = DistanceMatrix(names, matrix)
        #getitem
        self.assertRaises(ValueError, dm.__getitem__, 'A')
        self.assertRaises(ValueError, dm.__getitem__, ('Alpha', 'A'))
        self.assertRaises(TypeError, dm.__getitem__, (1, 'A'))
        self.assertRaises(TypeError, dm.__getitem__, (1, 1.2))
        self.assertRaises(IndexError, dm.__getitem__, 6)
        self.assertRaises(IndexError, dm.__getitem__, (10, 10))
        #setitem: item or index test
        self.assertRaises(ValueError, dm.__setitem__, 'A', [1, 3, 4])
        self.assertRaises(ValueError, dm.__setitem__, ('Alpha', 'A'), 4)
        self.assertRaises(TypeError, dm.__setitem__, (1, 'A'), 3)
        self.assertRaises(TypeError, dm.__setitem__, (1, 1.2), 2)
        self.assertRaises(IndexError, dm.__setitem__, 6, [1, 3, 4])
        self.assertRaises(IndexError, dm.__setitem__, (10,10), 1)
        #setitem: value test
        self.assertRaises(ValueError, dm.__setitem__, 0, [1, 2])
        self.assertRaises(TypeError, dm.__setitem__, ('Alpha', 'Beta'), 'a')
        self.assertRaises(TypeError, dm.__setitem__, 'Alpha', ['a', 'b', 'c'])

class DistanceCaluculatorTest(unittest.TestCase):
    """Test DistanceCaluculator"""

    def test_distance_calculator(self):
        aln = AlignIO.read(open('TreeConstruction/msa.phy'), 'phylip')

        calculator = DistanceCaluculator(aln, 'identity')
        dm = calculator.get_distance()
        self.assertEqual(dm['Alpha', 'Beta'], 1 - (10 * 1.0 / 13))

        calculator = DistanceCaluculator(aln, 'blastn')
        dm = calculator.get_distance()
        self.assertEqual(dm['Alpha', 'Beta'], 1 - (38 * 1.0 / 65))

        calculator = DistanceCaluculator(aln, 'trans')
        dm = calculator.get_distance()
        self.assertEqual(dm['Alpha', 'Beta'], 1 - (49 * 1.0 / 78))

        calculator = DistanceCaluculator(aln, 'blosum62')
        dm = calculator.get_distance()
        self.assertEqual(dm['Alpha', 'Beta'], 1 - (53 * 1.0 / 84))

if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)