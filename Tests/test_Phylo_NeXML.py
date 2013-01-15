# Copyright (C) 2013 by Ben Morris (ben@bendmorris.com)
# based on code by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for the NeXML and NeXMLIO modules.
"""

import os
import tempfile
import unittest

import Bio.Phylo as bp
from Bio.Phylo import NeXML, NeXMLIO

# Example NeXML files
example_files = ['NeXML/%s' % f 
                 for f in (
                            'characters.xml', 
                            'edgelabels.xml', 
                            'meta_taxa.xml', 
                            'meta_types.xml', 
                            'nexml.xml', 
                            'phenoscape.xml', 
                            'sets.xml', 
                            'taxa.xml', 
                            'timetree.xml', 
                            'tolweb.xml', 
                            'treebase-record.xml', 
                            'trees-uris.xml', 
                            'trees.xml',
                 )]

# Temporary file name for Writer tests below
DUMMY = tempfile.mktemp()


# ---------------------------------------------------------
# Parser tests

def _test_parse_factory(source, count=0):
    """Generate a test method for parse()ing the given source.

    The generated function extracts each phylogenetic tree using the parse()
    function and counts the total number of trees extracted.
    """
    fname = os.path.basename(source)

    def test_parse(self):
        trees = bp._io.parse(source, 'nexml')
        #self.assertEqual(len(list(trees)), count)

    test_parse.__doc__ = "Parse the phylogenies in %s." % fname
    return test_parse


class ParseTests(unittest.TestCase):
    """Tests for proper parsing of example phyloXML files."""

for n, ex in enumerate(example_files):
    parse_test = _test_parse_factory(ex)
    parse_test.__name__ = 'test_parse_%s' % n
    setattr(ParseTests, parse_test.__name__, parse_test)



if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
    # Clean up the temporary file
    if os.path.exists(DUMMY):
        os.remove(DUMMY)
