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

    def test_topology_3_clades(self):
        """Check simple tree topologies, three clades deep."""
        for source, counts in zip(self.all_examples, (
            # num phylogenies, lvl-2 clades, sub-clade counts, lvl-3 clades
            (1, ( ( (2, (2, 2)),
                    (2, (2, 2)),
                  ),
                ),
            ),  # apaf
            (1, ( ( (2, (2, 2)),
                    (2, (2, 2)),
                  ),
                ),
            ),  # bcl_2
            (13, (( (2, (0, 0)),
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
            ), # phylo_example
            (1, ( ( (3, (5, 1, 4)),
                    (5, (6, 2, 2, 2, 1)),
                    (2, (1, 1)),
                    (1, (2,)),
                    (2, (4, 4)),
                    (2, (2, 2)),
                    (1, (1,)),
                  ),
                ),
            ),   # ncbi zip
            )):
            phylo = PhyloXML.read(source)
            self.assertEquals(len(phylo), counts[0])
            self.assertEquals(len(phylo), len(counts[1]))
            for tree, count1 in zip(phylo, counts[1]):
                # warnings.warn("%s :: %s" % (source, count1))
                self.assertEquals(len(tree), 1) # True for all examples...
                self.assertEquals(len(tree.clades[0]), len(count1))
                for clade, subcounts in zip(tree.clades[0], count1):
                    self.assertEquals(len(clade), subcounts[0])
                    for subclade, size in zip(clade, subcounts[1]):
                        # warnings.warn("\tcompare %d == %d"
                        #               % (len(subclade), size))
                        self.assertEquals(len(subclade), size)


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
