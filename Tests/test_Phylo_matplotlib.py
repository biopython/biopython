# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for Bio.Phylo functions with external dependencies."""

import unittest

from Bio._py3k import StringIO
from Bio import Phylo

# Check for any missing dependencies at the top level so we can skip
from Bio import MissingExternalDependencyError

try:
    import matplotlib
except ImportError:
    raise MissingExternalDependencyError(
            "Install matplotlib if you want to use Bio.Phylo._utils.")

# Don't use the Wx backend for matplotlib, use the simpler postscript
# backend -- we're not going to display or save the plot anyway, so it
# doesn't matter much, as long as it's not Wx.  See:
# http://lists.open-bio.org/pipermail/biopython-dev/2012-April/009559.html
matplotlib.use('ps')
try:
    from matplotlib import pyplot
except ImportError:
    # Can fail here with font problems
    raise MissingExternalDependencyError(
            "Install matplotlib if you want to use Bio.Phylo._utils.")


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

    def test_draw_with_label_colors_dict(self):
        """Layout tree with label colors as dict.

        Run the tree layout algorithm with a label_colors argument passed in
        as a dictionary. Don't display tree.
        """
        pyplot.ioff()   # Turn off interactive display
        dollo = Phylo.read(EX_DOLLO, 'phyloxml')
        apaf = Phylo.read(EX_APAF, 'phyloxml')
        label_colors_dollo = {
            'f_50': 'red',
            'f_34': 'blue',
        }
        label_colors_apaf = {
            '22_MOUSE': 'red',
            '18_NEMVE': 'blue',
        }
        Phylo.draw(dollo, label_colors=label_colors_dollo, do_show=False)
        Phylo.draw(apaf, label_colors=label_colors_apaf, do_show=False)

    def test_draw_with_label_colors_callable(self):
        """Layout tree with label colors as callable.

        Run the tree layout algorithm with a label_colors argument passed in
        as a callable. Don't display tree.
        """
        pyplot.ioff()   # Turn off interactive display
        dollo = Phylo.read(EX_DOLLO, 'phyloxml')
        apaf = Phylo.read(EX_APAF, 'phyloxml')

        label_colors_dollo = lambda label: 'r' if label == 'f_50' else 'k'
        label_colors_apaf = lambda label: 'r'

        Phylo.draw(dollo, label_colors=label_colors_dollo, do_show=False)
        Phylo.draw(apaf, label_colors=label_colors_apaf, do_show=False)

    def test_draw_ascii(self):
        """Tree to Graph conversion."""
        handle = StringIO()
        tree = Phylo.read(EX_APAF, 'phyloxml')
        Phylo.draw_ascii(tree, file=handle)
        Phylo.draw_ascii(tree, file=handle, column_width=120)
        handle.close()


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
