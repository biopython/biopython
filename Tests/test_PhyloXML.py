# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for the Bio.PhyloXML module.

"""

import os
import unittest
import warnings
import zipfile
from itertools import izip
from cStringIO import StringIO

from Bio import PhyloXML
from Bio.Tree import PhyloXMLTree as Tree
from Bio.Tree import Utils


# Example PhyloXML files
EX_APAF = 'PhyloXML/apaf.xml'
EX_BCL2 = 'PhyloXML/bcl_2.xml'
EX_PHYLO = 'PhyloXML/phyloxml_examples.xml'
EX_DOLLO = 'PhyloXML/o_tol_332_d_dollo.xml'
EX_MOLLUSCA = 'PhyloXML/ncbi_taxonomy_mollusca.xml.zip'
# Big files - not checked into git yet
EX_TOL = 'PhyloXML/tol_life_on_earth_1.xml.zip'
EX_METAZOA = 'PhyloXML/ncbi_taxonomy_metazoa.xml.zip'
EX_NCBI = 'PhyloXML/ncbi_taxonomy.xml.zip'


def unzip(fname):
    """Extract a single file from a Zip archive and return a handle to it."""
    assert zipfile.is_zipfile(fname)
    z = zipfile.ZipFile(fname)
    return StringIO(z.read(z.filelist[0].filename))


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
            output.seek(0)
            self.assertEquals(len(output.readlines()), count)

    def test_pretty_print(self):
        """Check pretty_print by counting lines of output for each example.

        The line counts are liable to change whenever the object constructors
        change.
        """
        for source, count in izip(
                (EX_APAF, EX_BCL2, EX_PHYLO, unzip(EX_MOLLUSCA),
                    # unzip(EX_METAZOA), unzip(EX_NCBI),
                    ),
                (387, 748, 164, 16208, 214912, 648554)):
            phx = PhyloXML.read(source)
            output = StringIO()
            Utils.pretty_print(phx, output)
            output.seek(0)
            self.assertEquals(len(output.readlines()), count)
            output = StringIO()
            Utils.pretty_print(phx, output, show_all=True)
            output.seek(0)
            self.assertEquals(len(output.readlines()), count)


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
    test_read_dollo = _test_read_factory(EX_DOLLO, (1, 0))
    test_read_mollusca = _test_read_factory(EX_MOLLUSCA, (1, 0))
    # test_read_metazoa = _test_read_factory(EX_METAZOA, (1, 0))
    # test_read_ncbi = _test_read_factory(EX_NCBI, (1, 0))

    test_parse_apaf = _test_parse_factory(EX_APAF, 1)
    test_parse_bcl2 = _test_parse_factory(EX_BCL2, 1)
    test_parse_phylo = _test_parse_factory(EX_PHYLO, 13)
    test_parse_dollo = _test_parse_factory(EX_DOLLO, 1)
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
        """Instantiation of Phyloxml objects."""
        phx = PhyloXML.read(EX_PHYLO)
        self.assert_(isinstance(phx, Tree.Phyloxml))
        self.assert_('schemaLocation' in phx.attributes)
        for tree in phx:
            self.assert_(isinstance(tree, Tree.Phylogeny))
        for otr in phx.other:
            self.assert_(isinstance(otr, Tree.Other))

    def test_Other(self):
        """Instantiation of Other objects."""
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
        """Instantiation of Phylogeny objects."""
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
        """Instantiation of Clade objects."""
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

    def test_Annotation(self):
        """Instantiation of Annotation objects."""
        tree = list(PhyloXML.parse(EX_PHYLO))[3]
        ann = tree.clade[1].sequences[0].annotations[0]
        self.assert_(isinstance(ann, Tree.Annotation))
        self.assertEqual(ann.desc, 'alcohol dehydrogenase')
        self.assertAlmostEqual(ann.confidence.value, 0.67)
        self.assertEqual(ann.confidence.type, 'probability')

    def test_BinaryCharacters(self):
        """Instantiation of BinaryCharacters objects."""
        tree = PhyloXML.parse(EX_DOLLO).next()
        bchars = tree.clade[0,0].binary_characters
        self.assert_(isinstance(bchars, Tree.BinaryCharacters))
        self.assertEqual(bchars.type, 'parsimony inferred')
        for name, count, value in (
                ('gained',  2, ['Cofilin_ADF', 'Gelsolin']),
                ('lost',    0, []),
                ('present', 2, ['Cofilin_ADF', 'Gelsolin']),
                ('absent',  None, []),
                ):
            self.assertEqual(getattr(bchars, name+'_count'), count)
            self.assertEqual(getattr(bchars, name), value)

    # BranchColor -- no example

    def test_CladeRelation(self):
        """Instantiation of CladeRelation objects."""
        tree = list(PhyloXML.parse(EX_PHYLO))[6]
        crel = tree.clade_relations[0]
        self.assert_(isinstance(crel, Tree.CladeRelation))
        self.assertEqual(crel.id_ref_0, 'b')
        self.assertEqual(crel.id_ref_1, 'c')
        self.assertEqual(crel.type, 'network_connection')

    def test_Confidence(self):
        """Instantiation of Confidence objects."""
        tree = PhyloXML.parse(EX_BCL2).next()
        conf = tree.clade[0].confidences[0]
        self.assert_(isinstance(conf, Tree.Confidence))
        self.assertEqual(conf.type, 'bootstrap')
        self.assertAlmostEqual(conf.value, 33.0)

    def test_Date(self):
        """Instantiation of Date objects."""
        tree = list(PhyloXML.parse(EX_PHYLO))[11]
        silurian = tree.clade[0,0].date
        devonian = tree.clade[0,1].date
        ediacaran = tree.clade[1].date
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
        """Instantiation of Distribution objects.

        Also checks Point type and safe Unicode handling (?).
        """
        tree = list(PhyloXML.parse(EX_PHYLO))[10]
        hirschweg = tree.clade[0,0].distributions[0]
        nagoya = tree.clade[0,1].distributions[0]
        eth_zurich = tree.clade[0,2].distributions[0]
        san_diego = tree.clade[1].distributions[0]
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
        """Instantiation of DomainArchitecture objects.

        Also checks ProteinDomain type.
        """
        tree = PhyloXML.parse(EX_APAF).next()
        clade = tree.clade[0,0,0,0,0,0,0,0,0,0]
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
            self.assertEqual(domain.start + 1, start)
            self.assertEqual(domain.end, end)
            self.assertAlmostEqual(domain.confidence, conf)
            self.assertEqual(domain.value, value)

    def test_Events(self):
        """Instantiation of Events objects."""
        tree = list(PhyloXML.parse(EX_PHYLO))[4]
        event_s = tree.clade.events
        self.assert_(isinstance(event_s, Tree.Events))
        self.assertEqual(event_s.speciations, 1)
        event_d = tree.clade[0].events
        self.assert_(isinstance(event_d, Tree.Events))
        self.assertEqual(event_d.duplications, 1)

    # Polygon -- not implemented

    def test_Property(self):
        """Instantiation of Property objects."""
        tree = list(PhyloXML.parse(EX_PHYLO))[8]
        for prop, id_ref, value in izip(
                tree.properties,
                ('id_a', 'id_b', 'id_c'),
                ('1200', '2300', '200')):
            self.assert_(isinstance(prop, Tree.Property))
            self.assertEqual(prop.id_ref, id_ref)
            self.assertEqual(prop.datatype, "xsd:integer")
            self.assertEqual(prop.ref, "NOAA:depth")
            self.assertEqual(prop.applies_to, "node")
            self.assertEqual(prop.unit, "METRIC:m" )
            self.assertEqual(prop.value, value)

    # Reference -- not implemented
    def test_Reference(self):
        """Instantiation of Reference objects."""
        tree = PhyloXML.parse(EX_DOLLO).next()
        reference = tree.clade[0,0,0,0,0,0].references[0]
        self.assert_(isinstance(reference, Tree.Reference))
        self.assertEqual(reference.doi, '10.1038/nature06614')
        self.assertEqual(reference.desc, None)

    def test_Sequence(self):
        """Instantiation of Sequence objects.

        Also checks Accession and Annotation types.
        """
        trees = list(PhyloXML.parse(EX_PHYLO))
        # Simple element with id_source
        seq0 = trees[4].clade[1].sequences[0]
        self.assert_(isinstance(seq0, Tree.Sequence))
        self.assertEqual(seq0.id_source, 'z')
        self.assertEqual(seq0.symbol, 'ADHX')
        self.assertEqual(seq0.accession.source, 'ncbi')
        self.assertEqual(seq0.accession.value, 'Q17335')
        self.assertEqual(seq0.name, 'alcohol dehydrogenase')
        self.assertEqual(seq0.annotations[0].ref, 'InterPro:IPR002085')
        # More complete elements
        seq1 = trees[5].clade[0,0].sequences[0]
        seq2 = trees[5].clade[0,1].sequences[0]
        seq3 = trees[5].clade[1].sequences[0]
        for seq, sym, acc, name, mol_seq, ann_refs in izip(
                (seq1, seq2, seq3),
                ('ADHX', 'RT4I1', 'ADHB'),
                ('P81431', 'Q54II4', 'Q04945'),
                ('Alcohol dehydrogenase class-3',
                 'Reticulon-4-interacting protein 1 homolog, ' \
                         'mitochondrial precursor',
                 'NADH-dependent butanol dehydrogenase B'),
                ('TDATGKPIKCMAAIAWEAKKPLSIEEVEVAPPKSGEVRIKILHSGVCHTD',
                 'MKGILLNGYGESLDLLEYKTDLPVPKPIKSQVLIKIHSTSINPLDNVMRK',
                 'MVDFEYSIPTRIFFGKDKINVLGRELKKYGSKVLIVYGGGSIKRNGIYDK'),
                (("EC:1.1.1.1", "GO:0004022"),
                 ("GO:0008270", "GO:0016491"),
                 ("GO:0046872", "KEGG:Tetrachloroethene degradation")),
                ):
            self.assert_(isinstance(seq, Tree.Sequence))
            self.assertEqual(seq.symbol, sym)
            self.assertEqual(seq.accession.source, 'UniProtKB')
            self.assertEqual(seq.accession.value, acc)
            self.assertEqual(seq.name, name)
            self.assertEqual(seq.mol_seq, mol_seq)
            self.assertEqual(seq.annotations[0].ref, ann_refs[0])
            self.assertEqual(seq.annotations[1].ref, ann_refs[1])

    def test_SequenceRelation(self):
        """Instantiation of SequenceRelation objects."""
        tree = list(PhyloXML.parse(EX_PHYLO))[4]
        for seqrel, id_ref_0, id_ref_1, type in izip(
                tree.sequence_relations,
                ('x', 'x', 'y'), ('y', 'z', 'z'),
                ('paralogy', 'orthology', 'orthology')):
            self.assert_(isinstance(seqrel, Tree.SequenceRelation))
            self.assertEqual(seqrel.id_ref_0, id_ref_0)
            self.assertEqual(seqrel.id_ref_1, id_ref_1)
            self.assertEqual(seqrel.type, type)

    def test_Taxonomy(self):
        """Instantiation of Taxonomy objects.

        Also checks Id type.
        """
        trees = list(PhyloXML.parse(EX_PHYLO))
        # Octopus
        tax5 = trees[5].clade[0,0].taxonomies[0]
        self.assert_(isinstance(tax5, Tree.Taxonomy))
        self.assertEqual(tax5.id.value, '6645')
        self.assertEqual(tax5.id.type, 'NCBI')
        self.assertEqual(tax5.code, 'OCTVU')
        self.assertEqual(tax5.scientific_name, 'Octopus vulgaris')
        # Nile monitor
        tax9 = trees[9].clade[0].taxonomies[0]
        self.assert_(isinstance(tax9, Tree.Taxonomy))
        self.assertEqual(tax9.id.value, '62046')
        self.assertEqual(tax9.id.type, 'NCBI')
        self.assertEqual(tax9.scientific_name, 'Varanus niloticus')
        self.assertEqual(tax9.common_names[0], 'Nile monitor')
        self.assertEqual(tax9.rank, 'species')

    def test_Uri(self):
        """Instantiation of Uri objects."""
        tree = list(PhyloXML.parse(EX_PHYLO))[9]
        uri = tree.clade.taxonomies[0].uri
        self.assert_(isinstance(uri, Tree.Uri))
        self.assertEqual(uri.desc, 'EMBL REPTILE DATABASE')
        self.assertEqual(uri.value,
                'http://www.embl-heidelberg.de/~uetz/families/Varanidae.html')


# ---------------------------------------------------------
# Serialization tests

class WriterTests(unittest.TestCase):
    """Tests for serialization of objects to phyloXML format."""
    def _stash_rewrite_and_call(self, fname, test_cases):
        """Safely run a series of tests on a parsed and rewritten file.

        Specifically: Parse a file, rename the source file to a backup, rewrite
        the file from the parsed object, check the rewritten file with the
        given series of test functions, then restore the original by renaming
        the backup copy.

        Python 2.4 support: This would make more sense as a context manager
        that simply handles renaming and finally restoring the original.
        """
        phx = PhyloXML.read(fname)
        os.rename(fname, fname + '~')
        try:
            PhyloXML.write(phx, fname)
            for cls, tests in test_cases:
                inst = cls('setUp')
                for test in tests:
                    getattr(inst, test)()
        finally:
            os.rename(fname + '~', fname)

    def test_apaf(self):
        """Round-trip parsing and serialization of apaf.xml."""
        self._stash_rewrite_and_call(EX_APAF, (
            (ParseTests, [
                'test_read_apaf', 'test_parse_apaf', 'test_shape_apaf']),
            (TreeTests, ['test_DomainArchitecture']),
            ))

    def test_bcl2(self):
        """Round-trip parsing and serialization of bcl_2.xml."""
        self._stash_rewrite_and_call(EX_BCL2, (
            (ParseTests, [
                'test_read_bcl2', 'test_parse_bcl2', 'test_shape_bcl2']),
            (TreeTests, ['test_Confidence']),
            ))

    def test_phylo(self):
        """Round-trip parsing and serialization of phyloxml_examples.xml."""
        self._stash_rewrite_and_call(EX_PHYLO, (
            (ParseTests, [
                'test_read_phylo', 'test_parse_phylo', 'test_shape_phylo']),
            (TreeTests, [
                'test_Phyloxml',   'test_Other',
                'test_Phylogeny',  'test_Clade',
                'test_Annotation', 'test_CladeRelation',
                'test_Date',       'test_Distribution',
                'test_Events',     'test_Property',
                'test_Sequence',   'test_SequenceRelation',
                'test_Taxonomy',   'test_Uri',
                ]),
            ))

    def test_dollo(self):
        """Round-trip parsing and serialization of o_tol_332_d_dollo.xml."""
        self._stash_rewrite_and_call(EX_DOLLO, (
            (ParseTests, ['test_read_dollo', 'test_parse_dollo']),
            (TreeTests, ['test_BinaryCharacters']),
            ))


# ---------------------------------------------------------
# Method tests

class MethodTests(unittest.TestCase):
    """Tests for methods on specific classes/objects."""
    def setUp(self):
        self.phyloxml = PhyloXML.read(EX_PHYLO)

    # Type conversions

    def test_clade_to_phylogeny(self):
        clade = self.phyloxml.phylogenies[0].clade[0]
        tree = clade.to_phylogeny(rooted=True)
        self.assert_(isinstance(tree, Tree.Phylogeny))

    def test_phylogeny_to_phyloxml(self):
        tree = self.phyloxml.phylogenies[0]
        doc = tree.to_phyloxml()
        self.assert_(isinstance(doc, Tree.Phyloxml))

    def test_sequence_conversion(self):
        pass

    # Syntax sugar

    def test_clade_getitem(self):
        """Clade.__getitem__: get sub-clades by extended indexing."""
        tree = self.phyloxml.phylogenies[3]
        self.assertEqual(tree.clade[0,0], tree.clade.clades[0].clades[0])
        self.assertEqual(tree.clade[0,1], tree.clade.clades[0].clades[1])
        self.assertEqual(tree.clade[1], tree.clade.clades[1])
        self.assertEqual(len(tree.clade[:]), len(tree.clade.clades))
        self.assertEqual(len(tree.clade[0,:]),
                         len(tree.clade.clades[0].clades))

    def test_phyloxml_getitem(self):
        """Phyloxml.__getitem__: get phylogenies by name or index."""
        self.assert_(self.phyloxml.phylogenies[9] is self.phyloxml[9])
        self.assert_(self.phyloxml['monitor lizards'] is self.phyloxml[9])
        self.assertEqual(len(self.phyloxml[:]), len(self.phyloxml))

    def test_events(self):
        """Events: Mapping-type behavior."""
        evts = self.phyloxml.phylogenies[4].clade.events
        # Container behavior: __len__, __contains__
        self.assertEquals(len(evts), 1)
        self.assertEqual('speciations' in evts, True)
        self.assertEqual('duplications' in evts, False)
        # Attribute access: __get/set/delitem__
        self.assertEqual(evts['speciations'], 1)
        self.assertRaises(KeyError, lambda k: evts[k], 'duplications')
        evts['duplications'] = 3
        self.assertEqual(evts.duplications, 3)
        self.assertEqual(len(evts), 2)
        del evts['speciations']
        self.assertEqual(evts.speciations, None)
        self.assertEquals(len(evts), 1)
        # Iteration: __iter__, keys, values, items
        self.assertEqual(list(iter(evts)), ['duplications'])
        self.assertEqual(evts.keys(), ['duplications'])
        self.assertEqual(evts.values(), [3])
        self.assertEqual(evts.items(), [('duplications', 3)])

    def test_singlular(self):
        """Clade, Phylogeny: Singular properties for plural attributes."""
        conf = Tree.Confidence(0.9, 'bootstrap')
        taxo = Tree.Taxonomy(rank='genus')
        # Clade.taxonomy, Clade.confidence
        clade = Tree.Clade(confidences=[conf], taxonomies=[taxo])
        self.assertEqual(clade.confidence.type, 'bootstrap')
        self.assertEqual(clade.taxonomy.rank, 'genus')
        # raise if len > 1
        clade.confidences.append(conf)
        self.assertRaises(RuntimeError, getattr, clade, 'confidence')
        clade.taxonomies.append(taxo)
        self.assertRaises(RuntimeError, getattr, clade, 'taxonomy')
        # raise if []
        clade.confidences = []
        self.assertRaises(RuntimeError, getattr, clade, 'confidence')
        clade.taxonomies = []
        self.assertRaises(RuntimeError, getattr, clade, 'taxonomy')
        # Phylogeny.confidence
        tree = Tree.Phylogeny(True, confidences=[conf])
        self.assertEqual(tree.confidence.type, 'bootstrap')
        tree.confidences.append(conf)
        self.assertRaises(RuntimeError, getattr, tree, 'confidence')
        tree.confidences = []
        self.assertRaises(RuntimeError, getattr, tree, 'confidence')

    # Other methods

    def test_find(self):
        """Clade, Phylogeny: find() method."""
        # From the docstring example
        tree = self.phyloxml.phylogenies[5]
        matches = list(tree.find(Tree.Taxonomy, code='OCTVU'))
        self.assertEqual(len(matches), 1)
        self.assert_(isinstance(matches[0], Tree.Taxonomy))
        self.assertEqual(matches[0].code, 'OCTVU')
        self.assertEqual(matches[0].scientific_name, 'Octopus vulgaris')
        # Iteration and regexps
        tree = self.phyloxml.phylogenies[10]
        for point, alt in izip(tree.find(geodetic_datum=r'WGS\d{2}'),
                              (472, 10, 452)):
            self.assert_(isinstance(point, Tree.Point))
            self.assertEqual(point.geodetic_datum, 'WGS84')
            self.assertAlmostEqual(point.alt, alt)
        # boolean filter
        for clade, name in izip(tree.find(name=True), list('ABCD')):
            self.assert_(isinstance(clade, Tree.Clade))
            self.assertEqual(clade.name, name)
        # class filter
        tree = self.phyloxml.phylogenies[4]
        events = list(tree.find(Tree.Events))
        self.assertEqual(len(events), 2)
        self.assertEqual(events[0].speciations, 1)
        self.assertEqual(events[1].duplications, 1)
        # integer filter
        tree = PhyloXML.parse(EX_APAF).next()
        domains = list(tree.find(start=5))
        self.assertEqual(len(domains), 8)
        for dom in domains:
            self.assertEqual(dom.start, 5)
            self.assertEqual(dom.value, 'CARD')

    def test_color_rgb(self):
        """BranchColor: to_rgb() method."""
        black = Tree.BranchColor(0, 0, 0)
        self.assertEqual(black.to_rgb(), '000000')
        white = Tree.BranchColor(255, 255, 255)
        self.assertEqual(white.to_rgb(), 'ffffff')
        green = Tree.BranchColor(14, 192, 113)
        self.assertEqual(green.to_rgb(), '0ec071')

# ---------------------------------------------------------

if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
