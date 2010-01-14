# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for the Bio.Phylo module."""

import unittest
import zipfile
from itertools import izip
from cStringIO import StringIO

from Bio import Phylo
from Bio.Phylo import PhyloXML


# Example PhyloXML files
EX_APAF = 'PhyloXML/apaf.xml'
EX_BCL2 = 'PhyloXML/bcl_2.xml'
EX_MADE = 'PhyloXML/made_up.xml'
EX_PHYLO = 'PhyloXML/phyloxml_examples.xml'
EX_MOLLUSCA = 'PhyloXML/ncbi_taxonomy_mollusca.xml.zip'


def unzip(fname):
    """Extract a single file from a Zip archive and return a handle to it."""
    assert zipfile.is_zipfile(fname)
    z = zipfile.ZipFile(fname)
    return StringIO(z.read(z.filelist[0].filename))


class UtilTests(unittest.TestCase):
    """Tests for various utility functions."""
    def test_pretty_print(self):
        """Check pretty_print by counting lines of output for each example.

        The line counts are liable to change whenever the object constructors
        change.
        """
        for source, count in izip(
                (EX_APAF, EX_BCL2, unzip(EX_MOLLUSCA),
                    # unzip(EX_METAZOA), unzip(EX_NCBI),
                    ),
                (386, 747, 16207, 214911, 648553)):
            tree = Phylo.read(source, 'phyloxml')
            output = StringIO()
            Phylo.pretty_print(tree, output)
            output.seek(0)
            self.assertEquals(len(output.readlines()), count)
            output = StringIO()
            Phylo.pretty_print(tree, output, show_all=True)
            output.seek(0)
            self.assertEquals(len(output.readlines()), count)


class TreeTests(unittest.TestCase):
    """Tests for methods on BaseTree.Tree objects."""
    # TODO: magic: iter, len, getitem

    def setUp(self):
        self.phylogenies = list(Phylo.parse(EX_PHYLO, 'phyloxml'))

    # Traversal methods

    def test_find_all(self):
        """TreeMixin: find_all() method."""
        # From the docstring example
        tree = self.phylogenies[5]
        matches = list(tree.find_all(PhyloXML.Taxonomy, code='OCTVU'))
        self.assertEqual(len(matches), 1)
        self.assert_(isinstance(matches[0], PhyloXML.Taxonomy))
        self.assertEqual(matches[0].code, 'OCTVU')
        self.assertEqual(matches[0].scientific_name, 'Octopus vulgaris')
        # Iteration and regexps
        tree = self.phylogenies[10]
        for point, alt in izip(tree.find_all(geodetic_datum=r'WGS\d{2}'),
                               (472, 10, 452)):
            self.assert_(isinstance(point, PhyloXML.Point))
            self.assertEqual(point.geodetic_datum, 'WGS84')
            self.assertAlmostEqual(point.alt, alt)
        # class filter
        tree = self.phylogenies[4]
        events = list(tree.find_all(PhyloXML.Events))
        self.assertEqual(len(events), 2)
        self.assertEqual(events[0].speciations, 1)
        self.assertEqual(events[1].duplications, 1)
        # integer filter
        tree = Phylo.read(EX_APAF, 'phyloxml')
        domains = list(tree.find_all(start=5))
        self.assertEqual(len(domains), 8)
        for dom in domains:
            self.assertEqual(dom.start, 5)
            self.assertEqual(dom.value, 'CARD')

    def test_find_clades(self):
        """TreeMixin: find_clades() method."""
        # boolean filter
        for clade, name in izip(self.phylogenies[10].find_clades(name=True),
                                list('ABCD')):
            self.assert_(isinstance(clade, PhyloXML.Clade))
            self.assertEqual(clade.name, name)
        # finding deeper attributes
        octo = list(self.phylogenies[5].find_clades(code='OCTVU'))
        self.assertEqual(len(octo), 1)
        self.assert_(isinstance(octo[0], PhyloXML.Clade))
        self.assertEqual(octo[0].taxonomies[0].code, 'OCTVU')

    def test_find_terminal(self):
        """TreeMixin: find_all() with terminal argument."""
        def iter_len(it, count=0):
            for elem in it: count += 1
            return count
        for tree, total, extern, intern in izip(
                self.phylogenies,
                (6, 6, 7, 18, 21, 27, 7, 9, 9, 19, 15, 9, 6),
                (3, 3, 3, 3,  3,  3,  3, 3, 3, 3,  4,  3, 3),
                (3, 3, 3, 3,  3,  3,  3, 3, 3, 3,  3,  3, 3),
                ):
            self.assertEqual(iter_len(tree.find_all()), total)
            self.assertEqual(iter_len(tree.find_all(terminal=True)), extern)
            self.assertEqual(iter_len(tree.find_all(terminal=False)), intern)

    def test_get_path(self):
        """TreeMixin: get_path() method."""
        path = self.phylogenies[1].get_path({'name': 'B'})
        self.assertEqual(len(path), 2)
        self.assertAlmostEqual(path[0].branch_length, 0.06)
        self.assertAlmostEqual(path[1].branch_length, 0.23)
        self.assertEqual(path[1].name, 'B')

    def test_trace(self):
        """TreeMixin: trace() method."""
        tree = self.phylogenies[1]
        path = tree.trace({'name': 'A'}, {'name': 'C'})
        self.assertEqual(len(path), 3)
        self.assertAlmostEqual(path[0].branch_length, 0.06)
        self.assertAlmostEqual(path[2].branch_length, 0.4)
        self.assertEqual(path[2].name, 'C')

    # Information methods

    def test_common_ancestor(self):
        """TreeMixin: common_ancestor() method."""
        tree = self.phylogenies[1]
        lca = tree.common_ancestor({'name': 'A'}, {'name': 'B'})
        self.assertEqual(lca, tree.clade[0])
        lca = tree.common_ancestor({'name': 'A'}, {'name': 'C'})
        self.assertEqual(lca, tree.clade)
        tree = self.phylogenies[10]
        lca = tree.common_ancestor({'name': 'A'}, {'name': 'B'}, {'name': 'C'})
        self.assertEqual(lca, tree.clade[0])

    def test_depths(self):
        """TreeMixin: depths() method."""
        tree = self.phylogenies[1]
        depths = tree.depths()
        self.assertEqual(len(depths), 5)
        for found, expect in zip(sorted(depths.itervalues()),
                                 [0, 0.060, 0.162, 0.290, 0.400]):
            self.assertAlmostEqual(found, expect)

    def test_distance(self):
        """TreeMixin: distance() method."""
        t = self.phylogenies[1]
        self.assertAlmostEqual(t.distance({'name': 'A'}), 0.162)
        self.assertAlmostEqual(t.distance({'name': 'B'}), 0.29)
        self.assertAlmostEqual(t.distance({'name': 'C'}), 0.4)
        self.assertAlmostEqual(t.distance({'name': 'A'}, {'name': 'B'}), 0.332)
        self.assertAlmostEqual(t.distance({'name': 'A'}, {'name': 'C'}), 0.562)
        self.assertAlmostEqual(t.distance({'name': 'B'}, {'name': 'C'}), 0.69)

    def test_total_branch_length(self):
        """TreeMixin: total_branch_length() method."""
        tree = self.phylogenies[1]
        self.assertAlmostEqual(tree.total_branch_length(), 0.792)
        self.assertAlmostEqual(tree.clade[0].total_branch_length(), 0.392)

    # Tree manipulation methods

    def test_collapse(self):
        """TreeMixin: collapse() method."""
        tree = self.phylogenies[1]
        parent = tree.collapse(tree.clade[0])
        self.assertEqual(len(parent), 3)
        for clade, name, blength in zip(parent,
                ('C', 'A', 'B'),
                (0.4, 0.162, 0.29)):
            self.assertEqual(clade.name, name)
            self.assertAlmostEqual(clade.branch_length, blength)

    def test_ladderize(self):
        """TreeMixin: ladderize() method."""
        def ordered_names(tree):
            return [n.name for n in tree.get_terminals()]
        tree = self.phylogenies[10]
        self.assertEqual(ordered_names(tree), list('ABCD'))
        tree.ladderize()
        self.assertEqual(ordered_names(tree), list('DABC'))
        tree.ladderize(reverse=True)
        self.assertEqual(ordered_names(tree), list('ABCD'))


# ---------------------------------------------------------

if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
