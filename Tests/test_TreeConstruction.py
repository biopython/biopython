# Copyright (C) 2013 by Yanbo Ye (yeyanbo289@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for the Bio.Phylo.TreeConstruction module."""

import os
import unittest
import tempfile

from io import StringIO
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo import BaseTree
from Bio.Phylo import TreeConstruction
from Bio.Phylo import Consensus
from Bio.Phylo.TreeConstruction import _Matrix
from Bio.Phylo.TreeConstruction import DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import ParsimonyScorer
from Bio.Phylo.TreeConstruction import NNITreeSearcher
from Bio.Phylo.TreeConstruction import ParsimonyTreeConstructor


temp_dir = tempfile.mkdtemp()


class DistanceMatrixTest(unittest.TestCase):
    """Test for DistanceMatrix construction and manipulation."""

    def setUp(self):
        self.names = ["Alpha", "Beta", "Gamma", "Delta"]
        self.matrix = [[0], [1, 0], [2, 3, 0], [4, 5, 6, 0]]

    def test_good_construction(self):
        dm = DistanceMatrix(self.names, self.matrix)
        self.assertIsInstance(dm, TreeConstruction.DistanceMatrix)
        self.assertEqual(dm.names[0], "Alpha")
        self.assertEqual(dm.matrix[2][1], 3)
        self.assertEqual(len(dm), 4)
        self.assertEqual(
            repr(dm),
            "DistanceMatrix(names=['Alpha', 'Beta', 'Gamma', 'Delta'], "
            "matrix=[[0], [1, 0], [2, 3, 0], [4, 5, 6, 0]])",
        )

    def test_bad_construction(self):
        self.assertRaises(
            TypeError,
            DistanceMatrix,
            ["Alpha", 100, "Gamma", "Delta"],
            [[0], [0.1, 0], [0.2, 0.3, 0], [0.4, 0.5, 0.6, 0]],
        )
        self.assertRaises(
            TypeError,
            DistanceMatrix,
            ["Alpha", "Beta", "Gamma", "Delta"],
            [[0], ["a"], [0.2, 0.3], [0.4, 0.5, 0.6]],
        )
        self.assertRaises(
            ValueError,
            DistanceMatrix,
            ["Alpha", "Alpha", "Gamma", "Delta"],
            [[0], [0.1], [0.2, 0.3], [0.4, 0.5, 0.6]],
        )
        self.assertRaises(
            ValueError,
            DistanceMatrix,
            ["Alpha", "Beta", "Gamma", "Delta"],
            [[0], [0.2, 0], [0.4, 0.5, 0.6]],
        )
        self.assertRaises(
            ValueError,
            DistanceMatrix,
            ["Alpha", "Beta", "Gamma", "Delta"],
            [[0], [0.1], [0.2, 0.3, 0.4], [0.4, 0.5, 0.6]],
        )

    def test_good_manipulation(self):
        dm = DistanceMatrix(self.names, self.matrix)
        # getitem
        self.assertEqual(dm[1], [1, 0, 3, 5])
        self.assertEqual(dm[2, 1], 3)
        self.assertEqual(dm[2][1], 3)
        self.assertEqual(dm[1, 2], 3)
        self.assertEqual(dm[1][2], 3)
        self.assertEqual(dm["Alpha"], [0, 1, 2, 4])
        self.assertEqual(dm["Gamma", "Delta"], 6)
        # setitem
        dm["Alpha"] = [0, 10, 20, 40]
        self.assertEqual(dm["Alpha"], [0, 10, 20, 40])
        # delitem insert item
        del dm[1]
        self.assertEqual(dm.names, ["Alpha", "Gamma", "Delta"])
        self.assertEqual(dm.matrix, [[0], [20, 0], [40, 6, 0]])
        dm.insert("Beta", [1, 0, 3, 5], 1)
        self.assertEqual(dm.names, self.names)
        self.assertEqual(dm.matrix, [[0], [1, 0], [20, 3, 0], [40, 5, 6, 0]])
        del dm["Alpha"]
        self.assertEqual(dm.names, ["Beta", "Gamma", "Delta"])
        self.assertEqual(dm.matrix, [[0], [3, 0], [5, 6, 0]])
        dm.insert("Alpha", [1, 2, 4, 0])
        self.assertEqual(dm.names, ["Beta", "Gamma", "Delta", "Alpha"])
        self.assertEqual(dm.matrix, [[0], [3, 0], [5, 6, 0], [1, 2, 4, 0]])

    def test_bad_manipulation(self):
        dm = DistanceMatrix(self.names, self.matrix)
        # getitem
        self.assertRaises(ValueError, dm.__getitem__, "A")
        self.assertRaises(ValueError, dm.__getitem__, ("Alpha", "A"))
        self.assertRaises(TypeError, dm.__getitem__, (1, "A"))
        self.assertRaises(TypeError, dm.__getitem__, (1, 1.2))
        self.assertRaises(IndexError, dm.__getitem__, 6)
        self.assertRaises(IndexError, dm.__getitem__, (10, 10))
        # setitem: item or index test
        self.assertRaises(ValueError, dm.__setitem__, "A", [1, 3, 4])
        self.assertRaises(ValueError, dm.__setitem__, ("Alpha", "A"), 4)
        self.assertRaises(TypeError, dm.__setitem__, (1, "A"), 3)
        self.assertRaises(TypeError, dm.__setitem__, (1, 1.2), 2)
        self.assertRaises(IndexError, dm.__setitem__, 6, [1, 3, 4])
        self.assertRaises(IndexError, dm.__setitem__, (10, 10), 1)
        # setitem: value test
        self.assertRaises(ValueError, dm.__setitem__, 0, [1, 2])
        self.assertRaises(TypeError, dm.__setitem__, ("Alpha", "Beta"), "a")
        self.assertRaises(TypeError, dm.__setitem__, "Alpha", ["a", "b", "c"])

    def test_format_phylip(self):
        dm = DistanceMatrix(self.names, self.matrix)
        handle = StringIO()
        dm.format_phylip(handle)
        lines = handle.getvalue().splitlines()
        self.assertEqual(len(lines), len(dm) + 1)
        self.assertTrue(lines[0].endswith(str(len(dm))))
        for name, line in zip(self.names, lines[1:]):
            self.assertTrue(line.startswith(name))


class DistanceCalculatorTest(unittest.TestCase):
    """Test DistanceCalculator."""

    def test_known_matrices(self):
        aln = AlignIO.read("TreeConstruction/msa.phy", "phylip")

        calculator = DistanceCalculator("identity")
        dm = calculator.get_distance(aln)
        self.assertEqual(dm["Alpha", "Beta"], 1 - 10 / 13)

        calculator = DistanceCalculator("blastn")
        dm = calculator.get_distance(aln)
        self.assertEqual(dm["Alpha", "Beta"], 1 - 38 / 65)

        calculator = DistanceCalculator("trans")
        dm = calculator.get_distance(aln)
        self.assertEqual(dm["Alpha", "Beta"], 1 - 54 / 65)

        calculator = DistanceCalculator("blosum62")
        dm = calculator.get_distance(aln)
        self.assertEqual(dm["Alpha", "Beta"], 1 - 53 / 84)

    def test_nonmatching_seqs(self):
        aln = AlignIO.read(StringIO(">Alpha\nA-A--\n>Gamma\n-Y-Y-"), "fasta")
        # With a proper scoring matrix -- no matches
        dmat = DistanceCalculator("blosum62").get_distance(aln)
        self.assertEqual(dmat["Alpha", "Alpha"], 0.0)
        self.assertEqual(dmat["Alpha", "Gamma"], 1.0)
        # Comparing characters only -- 4 misses, 1 match
        dmat = DistanceCalculator().get_distance(aln)
        self.assertEqual(dmat["Alpha", "Alpha"], 0.0)
        self.assertAlmostEqual(dmat["Alpha", "Gamma"], 4.0 / 5.0)


class DistanceTreeConstructorTest(unittest.TestCase):
    """Test DistanceTreeConstructor."""

    def setUp(self):
        self.aln = AlignIO.read("TreeConstruction/msa.phy", "phylip")
        calculator = DistanceCalculator("blosum62")
        self.dm = calculator.get_distance(self.aln)
        self.constructor = DistanceTreeConstructor(calculator)

    def test_upgma(self):
        tree = self.constructor.upgma(self.dm)
        self.assertIsInstance(tree, BaseTree.Tree)
        # tree_file = StringIO()
        # Phylo.write(tree, tree_file, 'newick')
        ref_tree = Phylo.read("./TreeConstruction/upgma.tre", "newick")
        self.assertTrue(Consensus._equal_topology(tree, ref_tree))
        # ref_tree.close()

    def test_nj(self):
        tree = self.constructor.nj(self.dm)
        self.assertIsInstance(tree, BaseTree.Tree)
        # tree_file = StringIO()
        # Phylo.write(tree, tree_file, 'newick')
        ref_tree = Phylo.read("./TreeConstruction/nj.tre", "newick")
        self.assertTrue(Consensus._equal_topology(tree, ref_tree))
        # ref_tree.close()

        # create a matrix of length 2
        calculator = DistanceCalculator("blosum62")
        self.min_dm = calculator.get_distance(self.aln)
        for i in range(len(self.min_dm) - 2):
            del self.min_dm[len(self.min_dm) - 1]

        min_tree = self.constructor.nj(self.min_dm)
        self.assertIsInstance(min_tree, BaseTree.Tree)

        ref_min_tree = Phylo.read("./TreeConstruction/nj_min.tre", "newick")
        self.assertTrue(Consensus._equal_topology(min_tree, ref_min_tree))

    def test_built_tree(self):
        tree = self.constructor.build_tree(self.aln)
        self.assertIsInstance(tree, BaseTree.Tree)
        # tree_file = StringIO()
        # Phylo.write(tree, tree_file, 'newick')
        ref_tree = Phylo.read("./TreeConstruction/nj.tre", "newick")
        self.assertTrue(Consensus._equal_topology(tree, ref_tree))
        # ref_tree.close()


class ParsimonyScorerTest(unittest.TestCase):
    """Test ParsimonyScorer."""

    def test_get_score(self):
        aln = AlignIO.read("TreeConstruction/msa.phy", "phylip")
        tree = Phylo.read("./TreeConstruction/upgma.tre", "newick")
        scorer = ParsimonyScorer()
        score = scorer.get_score(tree, aln)
        self.assertEqual(score, 2 + 1 + 2 + 2 + 1 + 1 + 1 + 3)

        alphabet = ["A", "T", "C", "G"]
        step_matrix = [[0], [2.5, 0], [2.5, 1, 0], [1, 2.5, 2.5, 0]]
        matrix = _Matrix(alphabet, step_matrix)
        scorer = ParsimonyScorer(matrix)
        score = scorer.get_score(tree, aln)
        self.assertEqual(score, 3.5 + 2.5 + 3.5 + 3.5 + 2.5 + 1 + 2.5 + 4.5)

        alphabet = [
            "A",
            "C",
            "D",
            "E",
            "F",
            "G",
            "H",
            "I",
            "K",
            "L",
            "M",
            "N",
            "P",
            "Q",
            "R",
            "1",
            "2",
            "T",
            "V",
            "W",
            "Y",
            "*",
            "-",
        ]
        step_matrix = [
            [0],
            [2, 0],
            [1, 2, 0],
            [1, 2, 1, 0],
            [2, 1, 2, 2, 0],
            [1, 1, 1, 1, 2, 0],
            [2, 2, 1, 2, 2, 2, 0],
            [2, 2, 2, 2, 1, 2, 2, 0],
            [2, 2, 2, 1, 2, 2, 2, 1, 0],
            [2, 2, 2, 2, 1, 2, 1, 1, 2, 0],
            [2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 0],
            [2, 2, 1, 2, 2, 2, 1, 1, 1, 2, 2, 0],
            [1, 2, 2, 2, 2, 2, 1, 2, 2, 1, 2, 2, 0],
            [2, 2, 2, 1, 2, 2, 1, 2, 1, 1, 2, 2, 1, 0],
            [2, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 0],
            [1, 1, 2, 2, 1, 2, 2, 2, 2, 1, 2, 2, 1, 2, 2, 0],
            [2, 1, 2, 2, 2, 1, 2, 1, 2, 2, 2, 1, 2, 2, 1, 2, 0],
            [1, 2, 2, 2, 2, 2, 2, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 0],
            [1, 2, 1, 1, 1, 1, 2, 1, 2, 1, 1, 2, 2, 2, 2, 2, 2, 2, 0],
            [2, 1, 2, 2, 2, 1, 2, 2, 2, 1, 2, 3, 2, 2, 1, 1, 2, 2, 2, 0],
            [2, 1, 1, 2, 1, 2, 1, 2, 2, 2, 3, 1, 2, 2, 2, 1, 2, 2, 2, 2, 0],
            [2, 1, 2, 1, 2, 1, 2, 2, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 1, 1, 0],
            [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0],
        ]

        matrix = _Matrix(alphabet, step_matrix)
        scorer = ParsimonyScorer(matrix)
        score = scorer.get_score(tree, aln)
        self.assertEqual(score, 3 + 1 + 3 + 3 + 2 + 1 + 2 + 5)


class NNITreeSearcherTest(unittest.TestCase):
    """Test NNITreeSearcher."""

    def test_get_neighbors(self):
        tree = Phylo.read("./TreeConstruction/upgma.tre", "newick")
        alphabet = ["A", "T", "C", "G"]
        step_matrix = [[0], [2.5, 0], [2.5, 1, 0], [1, 2.5, 2.5, 0]]
        matrix = _Matrix(alphabet, step_matrix)
        scorer = ParsimonyScorer(matrix)
        searcher = NNITreeSearcher(scorer)
        trees = searcher._get_neighbors(tree)
        self.assertEqual(len(trees), 2 * (5 - 3))
        Phylo.write(trees, os.path.join(temp_dir, "neighbor_trees.tre"), "newick")


class ParsimonyTreeConstructorTest(unittest.TestCase):
    """Test ParsimonyTreeConstructor."""

    def test_build_tree(self):
        aln = AlignIO.read("TreeConstruction/msa.phy", "phylip")
        tree1 = Phylo.read("./TreeConstruction/upgma.tre", "newick")
        tree2 = Phylo.read("./TreeConstruction/nj.tre", "newick")
        alphabet = ["A", "T", "C", "G"]
        step_matrix = [[0], [2.5, 0], [2.5, 1, 0], [1, 2.5, 2.5, 0]]
        matrix = _Matrix(alphabet, step_matrix)
        scorer = ParsimonyScorer(matrix)
        searcher = NNITreeSearcher(scorer)
        constructor = ParsimonyTreeConstructor(searcher, tree1)
        best_tree = constructor.build_tree(aln)
        Phylo.write(best_tree, os.path.join(temp_dir, "pars1.tre"), "newick")
        constructor.starting_tree = tree2
        best_tree = constructor.build_tree(aln)
        Phylo.write(best_tree, os.path.join(temp_dir, "pars2.tre"), "newick")
        constructor.starting_tree = None
        best_tree = constructor.build_tree(aln)
        Phylo.write(best_tree, os.path.join(temp_dir, "pars3.tre"), "newick")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
