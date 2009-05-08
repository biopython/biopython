# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for the Bio.PhyloXML module.

"""

import unittest
import warnings

try:
    from xml.etree import cElementTree as ElementTree
except ImportError:
    try:
        from xml.etree import ElementTree as ElementTree
    except ImportError:
        # Python 2.4 -- check for 3rd-party implementations
        try:
            from lxml.etree import ElementTree
        except ImportError:
            try:
                import cElementTree as ElementTree
            except ImportError:
                try:
                    from elementtree import ElementTree
                except ImportError:
                    from Bio import MissingExternalDependencyError
                    raise MissingExternalDependencyError(
                            "No ElementTree module was found. " \
                            "Use Python 2.5+, lxml or elementtree if you " \
                            "want to use Bio.PhyloXML.")

from Bio import PhyloXML

# Example PhyloXML files
example_apaf = 'PhyloXML/apaf.xml'
example_bcl2 = 'PhyloXML/bcl_2.xml'
example_zip = 'PhyloXML/ncbi_taxonomy_mollusca.xml.zip'


class ParseNoOp(unittest.TestCase):
    """Tests for basic availability of library functions needed for parsing."""
    def test_noop(self):
        """Parse an XML document and dump its tags to standard output."""
        self._dump_tags(example_apaf)
        self._dump_tags(example_bcl2)

    def test_zip(self):
        """Parse a Zip-compressed XML file and dump tags to standard output."""
        import zipfile
        z = zipfile.ZipFile(example_zip)
        self._dump_tags(z.open(z.filelist[0].filename))
        z.close()

    @staticmethod
    def _dump_tags(source):
        events = ('start', 'end')
        for event, elem in ElementTree.iterparse(source, events=events):
            if event == 'start':
                print elem.tag
            else:
                elem.clear()


class ParsePhylo(unittest.TestCase):
    """Tests for proper parsing of example phyloXML files."""
    def test_root(self):
        """Read small example files to produce a tree object."""
        for source in (example_apaf, example_bcl2):
            tree = PhyloXML.read(source)
            self.assert_(tree)

    def test_zip(self):
        """Read a large Zip-compressed file to produce a tree object."""
        tree = PhyloXML.read(example_zip)
        self.assert_(tree)

    def test_core(self):
        """Verify the presence of core elements within the tree."""
        for source in (example_apaf, example_bcl2):
            tree = PhyloXML.read(source)
            self.assert_(len(tree.clades))


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
