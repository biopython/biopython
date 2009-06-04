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
example_phylo = 'PhyloXML/phyloxml_examples.xml'
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
        for source in (example_apaf, example_bcl2, example_phylo):
            Parser._dump_tags(source)

    def test_zip(self):
        """Parse a Zip-compressed XML file and dump tags to standard output."""
        Parser._dump_tags(unzip(example_zip))


class ParsePhylo(unittest.TestCase):
    """Tests for proper parsing of example phyloXML files."""
    def setUp(self):
        self.all_examples = (example_apaf, example_bcl2,
                             example_phylo, unzip(example_zip))

    def test_create(self):
        """Read each example file to produce a phyloXML object."""
        for source in self.all_examples:
            phylo = PhyloXML.read(source)
            self.assert_(phylo)

    def test_count_phylogenies(self):
        """Count the number of phylogenies within each file."""
        for source, count in zip(self.all_examples, (1, 1, 13, 1)):
            phylo = PhyloXML.read(source)
            self.assertEquals(len(phylo.phylogenies), count)


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
