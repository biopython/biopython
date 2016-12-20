# Copyright (C) 2013 by Yanbo Ye (yeyanbo289@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for the Bio.Phylo.Consensus module."""

import os
import unittest
import tempfile

# from Bio._py3k import StringIO
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo import BaseTree
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo import Consensus
from Bio.Phylo.Consensus import _BitString


temp_dir = tempfile.mkdtemp()


class BitStringTest(unittest.TestCase):
    """Test for _BitString class"""
    def test_bitstring(self):
        bitstr1 = _BitString('0011')
        bitstr2 = _BitString('0101')
        bitstr3 = _BitString('0001')
        bitstr4 = _BitString('0010')
        self.assertRaises(TypeError, _BitString, '10O1')
        self.assertEqual(bitstr1 & bitstr2, _BitString('0001'))
        self.assertEqual(bitstr1 | bitstr2, _BitString('0111'))
        self.assertEqual(bitstr1 ^ bitstr2, _BitString('0110'))
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
        bitstr_counts, len_trees = Consensus._count_clades(self.trees)
        self.assertEqual(len_trees, len(self.trees))
        self.assertEqual(len(bitstr_counts), 6)
        self.assertEqual(bitstr_counts[_BitString('11111')][0], 3)
        self.assertEqual(bitstr_counts[_BitString('11000')][0], 2)
        self.assertEqual(bitstr_counts[_BitString('00111')][0], 3)
        self.assertEqual(bitstr_counts[_BitString('00110')][0], 2)
        self.assertEqual(bitstr_counts[_BitString('00011')][0], 1)
        self.assertEqual(bitstr_counts[_BitString('01111')][0], 1)

    def test_strict_consensus(self):
        ref_trees = list(Phylo.parse('./TreeConstruction/strict_refs.tre', 'newick'))
        # three trees
        consensus_tree = Consensus.strict_consensus(self.trees)
        # tree_file = StringIO()
        # Phylo.write(consensus_tree, tree_file, 'newick')
        self.assertTrue(Consensus._equal_topology(consensus_tree, ref_trees[0]))
        # tree 1 and tree 2
        consensus_tree = Consensus.strict_consensus(self.trees[:2])
        # tree_file = StringIO()
        # Phylo.write(consensus_tree, tree_file, 'newick')
        self.assertTrue(Consensus._equal_topology(consensus_tree, ref_trees[1]))
        # tree 1 and tree 3
        consensus_tree = Consensus.strict_consensus(self.trees[::2])
        # tree_file = StringIO()
        # Phylo.write(consensus_tree, tree_file, 'newick')
        self.assertTrue(Consensus._equal_topology(consensus_tree, ref_trees[2]))
        # tree_file.close()

    def test_majority_consensus(self):
        ref_trees = Phylo.parse('./TreeConstruction/majority_ref.tre', 'newick')
        ref_tree = next(ref_trees)
        consensus_tree = Consensus.majority_consensus(self.trees)
        self.assertTrue(Consensus._equal_topology(consensus_tree, ref_tree))
        ref_tree = next(ref_trees)
        consensus_tree = Consensus.majority_consensus(self.trees, 1)
        self.assertTrue(Consensus._equal_topology(consensus_tree, ref_tree))

    def test_adam_consensus(self):
        ref_trees = list(Phylo.parse('./TreeConstruction/adam_refs.tre', 'newick'))
        # three trees
        consensus_tree = Consensus.adam_consensus(self.trees)
        # tree_file = '/home/yeyanbo/adam.tres'
        # tree_file = StringIO()
        # Phylo.write(consensus_tree, tree_file, 'newick')
        self.assertTrue(Consensus._equal_topology(consensus_tree, ref_trees[0]))
        # tree 1 and tree 2
        consensus_tree = Consensus.adam_consensus(self.trees[:2])
        # tree_file = StringIO()
        # Phylo.write(consensus_tree, tree_file, 'newick')
        self.assertTrue(Consensus._equal_topology(consensus_tree, ref_trees[1]))
        # tree 1 and tree 3
        consensus_tree = Consensus.adam_consensus(self.trees[::2])
        # tree_file = StringIO()
        # Phylo.write(consensus_tree, tree_file, 'newick')
        self.assertTrue(Consensus._equal_topology(consensus_tree, ref_trees[2]))
        # tree_file.close()

    def test_get_support(self):
        support_tree = Consensus.get_support(self.trees[0], self.trees)
        clade = support_tree.common_ancestor([support_tree.find_any(name="Beta"), support_tree.find_any(name="Gamma")])
        self.assertEqual(clade.confidence, 2 * 100.0 / 3)
        clade = support_tree.common_ancestor([support_tree.find_any(name="Alpha"), support_tree.find_any(name="Beta")])
        self.assertEqual(clade.confidence, 3 * 100.0 / 3)
        clade = support_tree.common_ancestor([support_tree.find_any(name="Delta"), support_tree.find_any(name="Epsilon")])
        self.assertEqual(clade.confidence, 2 * 100.0 / 3)


class BootstrapTest(unittest.TestCase):
    """Test for bootstrap methods"""

    def setUp(self):
        self.msa = AlignIO.read('TreeConstruction/msa.phy', 'phylip')

    def test_bootstrap(self):
        msa_list = list(Consensus.bootstrap(self.msa, 100))
        self.assertEqual(len(msa_list), 100)
        self.assertEqual(len(msa_list[0]), len(self.msa))
        self.assertEqual(len(msa_list[0][0]), len(self.msa[0]))

    def test_bootstrap_trees(self):
        calculator = DistanceCalculator('blosum62')
        constructor = DistanceTreeConstructor(calculator)
        trees = list(Consensus.bootstrap_trees(self.msa, 100, constructor))
        self.assertEqual(len(trees), 100)
        self.assertTrue(isinstance(trees[0], BaseTree.Tree))

    def test_bootstrap_consensus(self):
        calculator = DistanceCalculator('blosum62')
        constructor = DistanceTreeConstructor(calculator, 'nj')
        tree = Consensus.bootstrap_consensus(self.msa, 100, constructor, Consensus.majority_consensus)
        self.assertTrue(isinstance(tree, BaseTree.Tree))
        Phylo.write(tree, os.path.join(temp_dir, 'bootstrap_consensus.tre'), 'newick')


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
