# Copyright (C) 2013 by Yanbo Ye (yeyanbo289@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for the Bio.Phylo.Consensus module."""
import unittest
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


class ConsensusTest(unittest.TestCase):
    """Test for consensus methods"""
    def test_consensus(self):
        trees = list(Phylo.parse('./TreeConstruction/trees.tre', 'newick'))
        # test _count_clades
        bitstr_counts = Consensus._count_clades(trees)
        self.assertEqual(len(bitstr_counts), 6)
        self.assertEqual(bitstr_counts[BitString('11111')], 3)
        self.assertEqual(bitstr_counts[BitString('11000')], 2)
        self.assertEqual(bitstr_counts[BitString('00111')], 3)
        self.assertEqual(bitstr_counts[BitString('00110')], 2)
        self.assertEqual(bitstr_counts[BitString('00011')], 1)
        self.assertEqual(bitstr_counts[BitString('01111')], 1)
        # test stric_consensus
        consensus_tree = strict_consensus(trees)
        Phylo.write(consensus_tree, './TreeConstruction/consensus.tre', 'newick')
        consensus_tree = strict_consensus(trees[:2])
        Phylo.write(consensus_tree, './TreeConstruction/consensus1.tre', 'newick')
        consensus_tree = strict_consensus(trees[::2])
        Phylo.write(consensus_tree, './TreeConstruction/consensus2.tre', 'newick')

if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
