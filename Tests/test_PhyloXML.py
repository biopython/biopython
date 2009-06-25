# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for the Bio.PhyloXML module.

"""

import os.path
import unittest
import warnings
import zipfile
from itertools import izip
from cStringIO import StringIO

from Bio import PhyloXML
from Bio.PhyloXML import Tree


# Example PhyloXML files
EX_APAF = 'PhyloXML/apaf.xml'
EX_BCL2 = 'PhyloXML/bcl_2.xml'
EX_PHYLO = 'PhyloXML/phyloxml_examples.xml'
EX_MOLLUSCA = 'PhyloXML/ncbi_taxonomy_mollusca.xml.zip'
# Big files - not checked into git yet
EX_METAZOA = 'PhyloXML/ncbi_taxonomy_metazoa.xml.zip'
EX_NCBI = 'PhyloXML/ncbi_taxonomy.xml.zip'


def unzip(fname):
    """Extract a single file from a Zip archive and return a handle to it."""
    assert zipfile.is_zipfile(fname)
    z = zipfile.ZipFile(fname)
    return z.open(z.filelist[0].filename)


# ---------------------------------------------------------
# Utility tests

class UtilTests(unittest.TestCase):
    """Tests for various PhyloXML utility functions."""
    def test_dump_tags(self):
        """Count and confirm the number of tags in each example XML file."""
        for source, count in izip(
                (EX_APAF, EX_BCL2, EX_PHYLO, unzip(EX_MOLLUSCA),
                    # unzip(EX_METAZOA), unzip(EX_NCBI),
                    ),
                (509, 1496, 287, 24311, 322367, 972830),
                ):
            output = StringIO()
            PhyloXML.dump_tags(source, output)
            output.reset()
            self.assertEquals(len(output.readlines()), count)

    def test_pretty_print(self):
        """Check pretty_print by counting lines of output for each example.

        The line counts are liable to change whenever the object constructors
        change.
        """
        for source, count, count_all in izip(
                (EX_APAF, EX_BCL2, EX_PHYLO, unzip(EX_MOLLUSCA),
                    # unzip(EX_METAZOA), unzip(EX_NCBI),
                    ),
                (387, 748, 136, 16208, 214912, 648554),
                (773, 1495, 267, 32415, 429823, 1297107),
                ):
            handle = PhyloXML.read(source)
            output = StringIO()
            PhyloXML.pretty_print(handle, output=output)
            output.reset()
            self.assertEquals(len(output.readlines()), count)
            # print "Obj:", source, len(output.readlines()), 'eq?', count
            output = StringIO()
            PhyloXML.pretty_print(handle, show_all=True, output=output)
            output.reset()
            self.assertEquals(len(output.readlines()), count_all)
            # print "All:", source, len(output.readlines()), 'eq?', count_all


# ---------------------------------------------------------
# Parser tests

def _test_read_factory(source, count):
    """Generate a test method for read()ing the given source."""
    old_doc = """Read each example file to produce a phyloXML object.

        Test for existence of the root node, and count the number of phylogenies
        within each file.
        """
    fname = os.path.basename(source)
    if zipfile.is_zipfile(source):
        source = unzip(source)
    def test_read(self):
        phx = PhyloXML.read(source)
        self.assert_(phx)
        self.assertEquals(len(phx), count[0])
        self.assertEquals(len(phx.other), count[1])
    test_read.__doc__ = "Read %s to produce a phyloXML object." % fname
    return test_read


def _test_parse_factory(source, count):
    """Generate a test method for parse()ing the given source."""
    old_doc = "Extract each phylogenetic tree using the parse() function."
    fname = os.path.basename(source)
    if zipfile.is_zipfile(source):
        source = unzip(source)
    def test_parse(self):
        trees = PhyloXML.parse(source)
        self.assertEquals(len(list(trees)), count)
    test_parse.__doc__ = "Parse the phylogenies in %s." % fname
    return test_parse


def _test_shape_factory(source, shapes):
    """Generate a test method for checking tree shapes.

    Counts the branches at each level of branching in a phylogenetic tree, 3
    clades deep.
    """
    old_doc = "Check simple tree topologies, three clades deep."
    fname = os.path.basename(source)
    if zipfile.is_zipfile(source):
        source = unzip(source)
    def test_shape(self):
        trees = PhyloXML.parse(source)
        # self.assertEquals(len(list(trees)), len(shapes))
        for tree, shape_expect in izip(trees, shapes):
            # print "%s :: %s" % (source, shape_expect))
            self.assertEquals(len(tree.clade), len(shape_expect))
            for clade, sub_expect in izip(tree.clade, shape_expect):
                self.assertEquals(len(clade), sub_expect[0])
                for subclade, len_expect in izip(clade, sub_expect[1]):
                    # print "\tcompare %d == %d" % (len(subclade), len_expect)
                    self.assertEquals(len(subclade), len_expect)
    test_shape.__doc__ = "Check the branching structure of %s." % fname
    return test_shape


class ParseTests(unittest.TestCase):
    """Tests for proper parsing of example phyloXML files."""

    test_read_apaf = _test_read_factory(EX_APAF, (1, 0))
    test_read_bcl2 = _test_read_factory(EX_BCL2, (1, 0))
    test_read_phylo = _test_read_factory(EX_PHYLO, (13, 1))
    test_read_mollusca = _test_read_factory(EX_MOLLUSCA, (1, 0))
    # test_read_metazoa = _test_read_factory(EX_METAZOA, (1, 0))
    # test_read_ncbi = _test_read_factory(EX_NCBI, (1, 0))

    test_parse_apaf = _test_parse_factory(EX_APAF, 1)
    test_parse_bcl2 = _test_parse_factory(EX_BCL2, 1)
    test_parse_phylo = _test_parse_factory(EX_PHYLO, 13)
    test_parse_mollusca = _test_parse_factory(EX_MOLLUSCA, 1)
    # test_parse_metazoa = _test_parse_factory(EX_METAZOA, 1)
    # test_parse_ncbi = _test_parse_factory(EX_NCBI, 1)

    test_shape_apaf = _test_shape_factory(EX_APAF,
            # lvl-2 clades, sub-clade counts, lvl-3 clades
                ( ( (2, (2, 2)),
                    (2, (2, 2)),
                  ),
                ),
            )
    test_shape_bcl2 = _test_shape_factory(EX_BCL2,
                ( ( (2, (2, 2)),
                    (2, (2, 2)),
                  ),
                ),
            )
    test_shape_phylo = _test_shape_factory(EX_PHYLO,
                ( ( (2, (0, 0)),
                    (0, ()),
                  ),
                  (
                    (2, (0, 0)),
                    (0, ()),
                  ),
                  ( (2, (0, 0)),
                    (0, ()),
                  ),
                  (
                    (2, (0, 0)),
                    (0, ()),
                  ),
                  ( (2, (0, 0)),
                    (0, ()),
                  ),
                  ( (2, (0, 0)),
                    (0, ()),
                  ),
                  ( (2, (0, 0)),
                    (0, ()),
                  ),
                  ( (2, (0, 0)),
                    (0, ()),
                  ),
                  ( (2, (0, 0)),
                    (0, ()),
                  ),
                  ( (0, ()),
                    (2, (0, 0)),
                  ),
                  ( (3, (0, 0, 0)),
                    (0, ()),
                  ),
                  ( (2, (0, 0)),
                    (0, ()),
                  ),
                  ( (2, (0, 0)),
                    (0, ()),
                  ),
                ),
            )
    test_shape_zip = _test_shape_factory(EX_MOLLUSCA,
                ( ( (3, (5, 1, 4)),
                    (5, (6, 2, 2, 2, 1)),
                    (2, (1, 1)),
                    (1, (2,)),
                    (2, (4, 4)),
                    (2, (2, 2)),
                    (1, (1,)),
                  ),
                ),
            )


class TreeTests(unittest.TestCase):
    """Tests for instantiation and attributes of each complex type."""
    # NB: also test check_str() regexps wherever they're used
    def test_Phyloxml(self):
        """Test instantiation of Phyloxml objects."""
        phx = PhyloXML.read(EX_PHYLO)
        self.assert_(isinstance(phx, Tree.Phyloxml))
        self.assert_('schemaLocation' in phx.attributes)
        for tree in phx:
            self.assert_(isinstance(tree, Tree.Phylogeny))
        for otr in phx.other:
            self.assert_(isinstance(otr, Tree.Other))

    def test_Other(self):
        """Test instantiation of Other objects."""
        phx = PhyloXML.read(EX_PHYLO)
        otr = phx.other[0]
        self.assert_(isinstance(otr, Tree.Other))
        self.assertEquals(otr.tag, 'alignment')
        self.assertEquals(otr.namespace, 'http://example.org/align')
        self.assertEquals(len(otr.children), 3)
        for child, name, value in izip(otr, ('A', 'B', 'C'), (
            'acgtcgcggcccgtggaagtcctctcct', 'aggtcgcggcctgtggaagtcctctcct',
            'taaatcgc--cccgtgg-agtccc-cct')):
            self.assertEquals(child.tag, 'seq')
            self.assertEquals(child.attributes['name'], name)
            self.assertEquals(child.value, value)

    def test_Phylogeny(self):
        """Test instantiation of Phylogeny objects."""
        trees = list(PhyloXML.parse(EX_PHYLO))
        # Monitor lizards
        self.assertEquals(trees[9].name, 'monitor lizards')
        self.assertEquals(trees[9].description,
                'a pylogeny of some monitor lizards')
        self.assertEquals(trees[9].rooted, True)
        # Network (unrooted)
        tree6 = trees[6]
        self.assertEquals(trees[6].name,
                'network, node B is connected to TWO nodes: AB and C')
        self.assertEquals(trees[6].rooted, False)

    def test_Clade(self):
        """Test instantiation of Clade objects."""
        # ENH: check node_id, width (float) -- need an example
        tree = list(PhyloXML.parse(EX_PHYLO))[6]
        clade_ab, clade_c = tree.clade.clades
        clade_a, clade_b = clade_ab.clades
        for clade, id_source, name, blen in izip(
                (clade_ab, clade_a, clade_b, clade_c),
                ('ab', 'a', 'b', 'c'),
                ('AB', 'A', 'B', 'C'),
                (0.06, 0.102, 0.23, 0.4)):
            self.assert_(isinstance(clade, Tree.Clade))
            self.assertEqual(clade.id_source, id_source)
            self.assertEqual(clade.name, name)
            self.assertAlmostEqual(clade.branch_length, blen)

    # BinaryCharacterList -- not implemented
    # BinaryCharacters -- not implemented
    # BranchColor -- no example

    def test_CladeRelation(self):
        """Test instantiation of CladeRelation objects."""
        tree = list(PhyloXML.parse(EX_PHYLO))[6]
        crel = tree.clade_relations[0]
        self.assert_(isinstance(crel, Tree.CladeRelation))
        self.assertEqual(crel.id_ref_0, 'b')
        self.assertEqual(crel.id_ref_1, 'c')
        self.assertEqual(crel.type, 'network_connection')

    def test_Confidence(self):
        """Test instantiation of Confidence objects."""
        tree = PhyloXML.parse(EX_BCL2).next()
        conf = tree.clade.clades[0].confidences[0]
        self.assert_(isinstance(conf, Tree.Confidence))
        self.assertEqual(conf.type, 'bootstrap')
        self.assertAlmostEqual(conf.value, 33.0)

    def test_Date(self):
        """Test instantiation of Date objects."""
        tree = list(PhyloXML.parse(EX_PHYLO))[11]
        silurian = tree.clade.clades[0].clades[0].date
        devonian = tree.clade.clades[0].clades[1].date
        ediacaran = tree.clade.clades[1].date
        for date, rang, desc, val in izip(
                (silurian, devonian, ediacaran),
                (10, 20, 30),
                ('Silurian', 'Devonian', 'Ediacaran'),
                (425, 320, 600)):
            self.assert_(isinstance(date, Tree.Date))
            self.assertEqual(date.unit, 'mya')
            self.assertAlmostEqual(date.range, rang)
            self.assertEqual(date.desc, desc)
            self.assertAlmostEqual(date.value, val)

    def test_Distribution(self):
        """Test instantiation of Distribution objects.

        Also checks Point type and safe Unicode handling (?).
        """
        tree = list(PhyloXML.parse(EX_PHYLO))[10]
        hirschweg = tree.clade.clades[0].clades[0].distributions[0]
        nagoya = tree.clade.clades[0].clades[1].distributions[0]
        eth_zurich = tree.clade.clades[0].clades[2].distributions[0]
        san_diego = tree.clade.clades[1].distributions[0]
        for dist, desc, lat, long, alt in izip(
                (hirschweg, nagoya, eth_zurich, san_diego),
                ('Hirschweg, Winterthur, Switzerland',
                    'Nagoya, Aichi, Japan',
                    u'ETH Z\xfcrich',
                    'San Diego'),
                (47.481277, 35.155904, 47.376334, 32.880933),
                (8.769303, 136.915863, 8.548108, -117.217543),
                (472, 10, 452, 104)):
            self.assert_(isinstance(dist, Tree.Distribution))
            self.assertEqual(dist.desc, desc)
            point = dist.points[0]
            self.assert_(isinstance(point, Tree.Point))
            self.assertEqual(point.geodetic_datum, 'WGS84')
            self.assertEqual(point.lat, lat)
            self.assertEqual(point.long, long)
            self.assertEqual(point.alt, alt)

    def test_DomainArchitecture(self):
        """Test instantiation of DomainArchitecture objects.

        Also checks ProteinDomain type.
        """
        tree = PhyloXML.parse(EX_APAF).next()
        clade = tree.clade.clades[0].clades[0] \
                .clades[0].clades[0].clades[0].clades[0] \
                .clades[0].clades[0].clades[0].clades[0]
        darch = clade.sequences[0].domain_architecture
        self.assert_(isinstance(darch, Tree.DomainArchitecture))
        self.assertEqual(darch.length, 1249)
        for domain, start, end, conf, value in izip(darch.domains,
                (6, 109, 605, 647, 689, 733, 872, 993, 1075, 1117, 1168),
                (90, 414, 643, 685, 729, 771, 910, 1031, 1113, 1155, 1204),
                (7.0e-26, 7.2e-117, 2.4e-6, 1.1e-12, 2.4e-7, 4.7e-14, 2.5e-8,
                    4.6e-6, 6.3e-7, 1.4e-7, 0.3),
                ('CARD', 'NB-ARC', 'WD40', 'WD40', 'WD40', 'WD40', 'WD40',
                    'WD40', 'WD40', 'WD40', 'WD40')):
            self.assert_(isinstance(domain, Tree.ProteinDomain))
            self.assertEqual(domain.start, start)
            self.assertEqual(domain.end, end)
            self.assertAlmostEqual(domain.confidence, conf)
            self.assertEqual(domain.value, value)

    def test_Events(self):
        """Test instantiation of Events objects."""
        pass

    def test_Point(self):
        """Test instantiation of Point objects."""
        pass

    # Polygon -- not implemented

    def test_Property(self):
        """Test instantiation of Property objects."""
        pass

    def test_ProteinDomain(self):
        """Test instantiation of ProteinDomain objects."""
        pass

    # Reference -- not implemented

    def test_Sequence(self):
        """Test instantiation of Sequence objects."""
        pass

    def test_SequenceRelation(self):
        """Test instantiation of SequenceRelation objects."""
        pass

    def test_Taxonomy(self):
        """Test instantiation of Taxonomy objects."""
        pass

    def test_Uri(self):
        """Test instantiation of Uri objects."""
        tree = list(PhyloXML.parse(EX_PHYLO))[9]
        uri = tree.clade.taxonomies[0].uri
        self.assert_(isinstance(uri, Tree.Uri))
        self.assertEqual(uri.desc, 'EMBL REPTILE DATABASE')
        self.assertEqual(uri.value,
                'http://www.embl-heidelberg.de/~uetz/families/Varanidae.html')


# ---------------------------------------------------------

if __name__ == '__main__':
    warnings.simplefilter('ignore')
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
