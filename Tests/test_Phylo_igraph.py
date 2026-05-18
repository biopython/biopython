# Copyright (C) 2023 by Fabio Zanini (fabio.zanini@unsw.edu.au)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for Bio.Phylo functions with external dependencies."""

import unittest
from io import StringIO

# Check for any missing dependencies at the top level so we can skip
from Bio import MissingExternalDependencyError
from Bio import Phylo

try:
    import igraph as ig
except ImportError:
    raise MissingExternalDependencyError(
        "Install igraph if you wish to use it with Bio.Phylo"
    ) from None

# Example PhyloXML file
EX_DOLLO = "PhyloXML/o_tol_332_d_dollo.xml"
EX_APAF = "PhyloXML/apaf.xml"


class UtilTests(unittest.TestCase):
    """Tests for various utility functions."""

    def test_to_igraph_simple(self):
        """Tree to Graph conversion, simple tree."""
        tree = Phylo.read(StringIO("()()"), format="newick")
        for i, elem in enumerate(tree.find_clades()):
            elem.name = f"elem_{i + 1}"
            elem.width = f"width_{i}"
        graph = Phylo.to_igraph(tree, vertex_attributes=["name"])
        self.assertEqual(graph.vs["name"], ["elem_1", "elem_2", "elem_3"])
        self.assertEqual(graph.es["width"], ["width_1", "width_2"])

    def test_to_igraph_exd_ollo(self):
        """Tree to Graph conversion, if networkx is available."""
        tree = Phylo.read(EX_DOLLO, "phyloxml")
        graph = Phylo.to_igraph(tree)
        self.assertEqual(graph.vcount(), 659)
        self.assertEqual(graph.ecount(), graph.vcount() - 1)

    def test_to_igraph_ex_dollo_with_attributes(self):
        tree = Phylo.read(EX_DOLLO, "phyloxml")
        graph = Phylo.to_igraph(tree, vertex_attributes=["name"])
        self.assertEqual(graph.vcount(), 659)
        self.assertEqual(graph.ecount(), graph.vcount() - 1)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
