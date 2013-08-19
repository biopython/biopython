# Copyright (C) 2013 by Yanbo Ye (yeyanbo289@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for the Bio.Phylo.Consensus module."""
import unittest
import StringIO
from Bio import Phylo
from Bio.Phylo import Consensus
from Bio.Phylo.Consensus import *


class BitStringTest(unittest.TestCase):
    """Test for BitString class"""
    def test_bitstring(self):
        bitstr1 = BitString('0011')
        bitstr2 = BitString('0101')
        bitstr3 = BitString('0001')
        bitstr4 = BitString('0010')
        self.assertRaises(TypeError, BitString, '10O1')
        self.assertEqual(bitstr1 & bitstr2, BitString('0001'))
        self.assertEqual(bitstr1 | bitstr2, BitString('0111'))
        self.assertEqual(bitstr1 ^ bitstr2, BitString('0110'))
        self.assertFalse(bitstr1.contains(bitstr2))
        self.assertTrue(bitstr1.contains(bitstr1))
        self.assertTrue(bitstr1.contains(bitstr3))
        self.assertTrue(bitstr1.contains(bitstr4))
        self.assertFalse(bitstr1.independent(bitstr2))
        self.assertFalse(bitstr1.independent(bitstr4))
        self.assertTrue(bitstr2.independent(bitstr4))
        self.assertTrue(bitstr3.independent(bitstr4))
        self.assertFalse(bitstr1.iscompatible(bitstr2))
        self.assertTrue(bitstr1.iscompatible(bitstr3))
        self.assertTrue(bitstr1.iscompatible(bitstr4))
        self.assertTrue(bitstr2.iscompatible(bitstr4))
        self.assertTrue(bitstr3.iscompatible(bitstr4))

class ConsensusTest(unittest.TestCase):
    """Test for consensus methods"""

    def setUp(self):
        self.trees = list(Phylo.parse('./TreeConstruction/trees.tre', 'newick'))

    def test_count_clades(self):
        bitstr_counts = Consensus._count_clades(self.trees)
        self.assertEqual(len(bitstr_counts), 6)
        self.assertEqual(bitstr_counts[BitString('11111')][0], 3)
        self.assertEqual(bitstr_counts[BitString('11000')][0], 2)
        self.assertEqual(bitstr_counts[BitString('00111')][0], 3)
        self.assertEqual(bitstr_counts[BitString('00110')][0], 2)
        self.assertEqual(bitstr_counts[BitString('00011')][0], 1)
        self.assertEqual(bitstr_counts[BitString('01111')][0], 1)

    def test_strict_consensus(self):
        ref_trees = open('./TreeConstruction/strict_refs.tre')
        # three trees
        consensus_tree = strict_consensus(self.trees)
        tree_file = StringIO.StringIO()
        Phylo.write(consensus_tree, tree_file, 'newick')
        self.assertEqual(tree_file.getvalue(), ref_trees.readline())
        # tree 1 and tree 2
        consensus_tree = strict_consensus(self.trees[:2])
        tree_file = StringIO.StringIO()
        Phylo.write(consensus_tree, tree_file, 'newick')
        self.assertEqual(tree_file.getvalue(), ref_trees.readline())
        # tree 1 and tree 3
        consensus_tree = strict_consensus(self.trees[::2])
        tree_file = StringIO.StringIO()
        Phylo.write(consensus_tree, tree_file, 'newick')
        self.assertEqual(tree_file.getvalue(), ref_trees.readline())
        ref_trees.close()
        tree_file.close()

    def test_majority_consensus(self):
        # three trees
        ref_tree = open('./TreeConstruction/majority_ref.tre')
        consensus_tree = majority_consensus(self.trees)
        tree_file = StringIO.StringIO()
        Phylo.write(consensus_tree, tree_file, 'newick')
        self.assertEqual(tree_file.getvalue(), ref_tree.readline())
        consensus_tree = majority_consensus(self.trees, 1)
        tree_file = StringIO.StringIO()
        Phylo.write(consensus_tree, tree_file, 'newick')
        self.assertEqual(tree_file.getvalue(), ref_tree.readline())

    def test_adam_consensus(self):
        ref_trees = open('./TreeConstruction/adam_refs.tre')
        # three trees
        consensus_tree = adam_consensus(self.trees)
        tree_file = StringIO.StringIO()
        Phylo.write(consensus_tree, tree_file, 'newick')
        self.assertEqual(tree_file.getvalue(), ref_trees.readline())
        # tree 1 and tree 2
        consensus_tree = adam_consensus(self.trees[:2])
        tree_file = StringIO.StringIO()
        Phylo.write(consensus_tree, tree_file, 'newick')
        self.assertEqual(tree_file.getvalue(), ref_trees.readline())
        # tree 1 and tree 3
        consensus_tree = adam_consensus(self.trees[::2])
        tree_file = StringIO.StringIO()
        Phylo.write(consensus_tree, tree_file, 'newick')
        self.assertEqual(tree_file.getvalue(), ref_trees.readline())
        ref_trees.close()
        tree_file.close()


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
