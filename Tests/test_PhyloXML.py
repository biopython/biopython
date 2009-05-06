# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for the Bio.PhyloXML module.

"""

import unittest
import warnings

try:
    from xml.etree import cElementTree as ETree
except ImportError:
    try:
        from xml.etree import ElementTree as ETree
    except ImportError:
        raise MissingExternalDependencyError(
                "Use Python 2.5+ if you want to use Bio.PhyloXML.")

# from Bio import PhyloXML


class ParseNoOp(unittest.TestCase):
    def test_noop(self):
        self._dump_tags('PhyloXML/apaf.xml')
        self._dump_tags('PhyloXML/bcl_2.xml')

    def test_zip(self):
        source = 'PhyloXML/ncbi_taxonomy_mollusca.xml.zip'
        # TODO: unzip this file & deal with it

    @staticmethod
    def _dump_tags(source):
        events = ('start', 'end')
        for event, elem in ETree.iterparse(source, events=events):
            if event == 'start':
                print elem.tag
            else:
                elem.clear()


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
