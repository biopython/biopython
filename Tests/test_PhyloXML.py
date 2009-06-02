# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for the Bio.PhyloXML module.

"""

import unittest
import warnings
import zipfile

from Bio import PhyloXML
from Bio.PhyloXML import Parser

# Example PhyloXML files
example_apaf = 'PhyloXML/apaf.xml'
example_bcl2 = 'PhyloXML/bcl_2.xml'
example_zip = 'PhyloXML/ncbi_taxonomy_mollusca.xml.zip'


def unzip(fname):
    """Extract a single file from a Zip archive and return a handle to it."""
    assert zipfile.is_zipfile(fname)
    z = zipfile.ZipFile(fname)
    return z.open(z.filelist[0].filename)


class ParseNoOp(unittest.TestCase):
    """Tests for basic availability of library functions needed for parsing."""
    def test_noop(self):
        """Parse an XML document and dump its tags to standard output."""
        for source in (example_apaf, example_bcl2):
            Parser._dump_tags(source)

    def test_zip(self):
        """Parse a Zip-compressed XML file and dump tags to standard output."""
        Parser._dump_tags(unzip(example_zip))


class ParsePhylo(unittest.TestCase):
    """Tests for proper parsing of example phyloXML files."""
    def test_root(self):
        """Read small example files to produce a tree object."""
        for source in (example_apaf, example_bcl2):
            tree = PhyloXML.read(source)
            self.assert_(tree)

    def test_zip(self):
        """Read a large Zip-compressed file to produce a tree object."""
        tree = PhyloXML.read(unzip(example_zip))
        self.assert_(tree)

    def test_core(self):
        """Verify the presence of core elements within the tree."""
        for source in (example_apaf, example_bcl2, unzip(example_zip)):
            tree = PhyloXML.read(source)
            self.assert_(len(tree.clades))


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
