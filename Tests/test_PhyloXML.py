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

from Bio import PhyloXML
from Bio.PhyloXML import Parser


# Example PhyloXML files
EX_APAF = 'PhyloXML/apaf.xml'
EX_BCL2 = 'PhyloXML/bcl_2.xml'
EX_PHYLO = 'PhyloXML/phyloxml_examples.xml'
EX_ZIP = 'PhyloXML/ncbi_taxonomy_mollusca.xml.zip'


def unzip(fname):
    """Extract a single file from a Zip archive and return a handle to it."""
    assert zipfile.is_zipfile(fname)
    z = zipfile.ZipFile(fname)
    return z.open(z.filelist[0].filename)


# ---------------------------------------------------------

class ParseNoOp(unittest.TestCase):
    """Tests for basic availability of library functions needed for parsing."""
    def test_noop(self):
        """Parse an XML document and dump its tags to standard output."""
        for source in (EX_APAF, EX_BCL2, EX_PHYLO):
            Parser._dump_tags(source)

    def test_zip(self):
        """Parse a Zip-compressed XML file and dump tags to standard output."""
        Parser._dump_tags(unzip(EX_ZIP))


# ---------------------------------------------------------

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
        phylo = PhyloXML.read(source)
        self.assert_(phylo)
        # warnings.warn("# phylogenies in %s: %s (expect %d)"
        #         "\n\trepr: %s" % (source, len(phylo.phylogenies), count,
        #             phylo.phylogenies))
        self.assertEquals(len(phylo), count)
        # TODO: check existence/length of 'other'
    test_read.__doc__ = "Read %s to produce a phyloXML object." % fname
    return test_read


def _test_parse_factory(source, count):
    """Generate a test method for parse()ing the given source."""
    old_doc = "Extract each phylogenetic tree using the parse() function."
    fname = os.path.basename(source)
    if zipfile.is_zipfile(source):
        source = unzip(source)
    def test_parse(self):
        phylo = PhyloXML.parse(source)
        self.assertEquals(len(list(phylo)), count)
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
        phylo = PhyloXML.parse(source)
        # self.assertEquals(len(phylo), len(shapes))
        for tree, shape_expect in izip(phylo, shapes):
            # warnings.warn("%s :: %s" % (source, shape_expect))
            self.assertEquals(len(tree.clade), len(shape_expect))
            for clade, sub_expect in izip(tree.clade, shape_expect):
                self.assertEquals(len(clade), sub_expect[0])
                for subclade, len_expect in izip(clade, sub_expect[1]):
                    # warnings.warn("\tcompare %d == %d"
                    #               % (len(subclade), len_expect))
                    self.assertEquals(len(subclade), len_expect)
    test_shape.__doc__ = "Check the branching structure of %s." % fname
    return test_shape


class ParsePhylo(unittest.TestCase):
    """Tests for proper parsing of example phyloXML files."""

    test_read_apaf = _test_read_factory(EX_APAF, 1)
    test_read_bcl2 = _test_read_factory(EX_BCL2, 1)
    test_read_phylo = _test_read_factory(EX_PHYLO, 13)
    test_read_zip = _test_read_factory(EX_ZIP, 1)

    test_parse_apaf = _test_parse_factory(EX_APAF, 1)
    test_parse_bcl2 = _test_parse_factory(EX_BCL2, 1)
    test_parse_phylo = _test_parse_factory(EX_PHYLO, 13)
    test_parse_zip = _test_parse_factory(EX_ZIP, 1)

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
    test_shape_zip = _test_shape_factory(EX_ZIP,
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


# ---------------------------------------------------------

if __name__ == '__main__':
    warnings.simplefilter('ignore')
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
