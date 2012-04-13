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

try:
    from matplotlib import pyplot
except ImportError:
    raise MissingExternalDependencyError(
            "Install matplotlib if you want to use Bio.Phylo._utils.")

# OK, we can go ahead
import unittest
from cStringIO import StringIO

from Bio import Phylo

# Example PhyloXML file
EX_DOLLO = 'PhyloXML/o_tol_332_d_dollo.xml'
EX_APAF = 'PhyloXML/apaf.xml'

class UtilTests(unittest.TestCase):
    """Tests for various utility functions."""

    def test_draw(self):
        """Run the tree layout algorithm, but don't display it."""
        pyplot.ioff()   # Turn off interactive display
        dollo = Phylo.read(EX_DOLLO, 'phyloxml')
        apaf = Phylo.read(EX_APAF, 'phyloxml')
        Phylo.draw(dollo, do_show=False)
        Phylo.draw(apaf, do_show=False)
        # Fancier options
        Phylo.draw(apaf, do_show=False, branch_labels={apaf.root: 'Root'})
        Phylo.draw(apaf, do_show=False, branch_labels=lambda c: c.branch_length)

    def test_draw_ascii(self):
        """Tree to Graph conversion, if networkx is available."""
        handle = StringIO()
        tree = Phylo.read(EX_APAF, 'phyloxml')
        Phylo.draw_ascii(tree, file=handle)
        Phylo.draw_ascii(tree, file=handle, column_width=120)
        handle.close()

    def test_to_networkx(self):
        """Tree to Graph conversion, if networkx is available."""
        tree = Phylo.read(EX_DOLLO, 'phyloxml')
        G = Phylo.to_networkx(tree)
        self.assertEqual(len(G.nodes()), 659)


# ---------------------------------------------------------

if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
