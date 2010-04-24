# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for Bio.Phylo functions with external dependencies."""

# Check for any missing dependencies at the top level so we can skip
from Bio import MissingExternalDependencyError

try:
    import networkx
except ImportError:
    raise MissingExternalDependencyError(
            "Install NetworkX if you want to use Bio.Phylo._utils.")

try:
    import pygraphviz
except ImportError:
    try:
        import pydot
    except ImportError:
        raise MissingExternalDependencyError(
                "Install PyGraphviz or Pydot if you want to use "
                "Bio.Phylo._utils.")

# OK, we can go ahead
import unittest

from Bio import Phylo

# Example PhyloXML file
EX_DOLLO = 'PhyloXML/o_tol_332_d_dollo.xml'

class UtilTests(unittest.TestCase):
    """Tests for various utility functions."""

    def test_to_networkx(self):
        """Tree to Graph conversion, if networkx is available."""
        tree = Phylo.read(EX_DOLLO, 'phyloxml')
        G = Phylo.to_networkx(tree)
        self.assertEqual(len(G.nodes()), 659)


# ---------------------------------------------------------

if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
