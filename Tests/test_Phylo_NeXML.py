# Copyright (C) 2013 by Ben Morris (ben@bendmorris.com)
# based on code by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for the NeXML and NeXMLIO modules."""

import os
import tempfile
import unittest

from io import BytesIO

from Bio import Phylo

# Example NeXML files
nexml_files = (
    "characters.xml",
    "edgelabels.xml",
    "meta_taxa.xml",
    "meta_types.xml",
    "nexml.xml",
    "phenoscape.xml",
    "sets.xml",
    "taxa.xml",
    "timetree.xml",
    "tolweb.xml",
    "treebase-record.xml",
    "trees-uris.xml",
    "trees.xml",
)

tree_counts = {
    "taxa.xml": 0,
    "timetree.xml": 38,
    "phenoscape.xml": 0,
    "nexml.xml": 0,
    "meta_types.xml": 0,
    "meta_taxa.xml": 0,
    "trees.xml": 2,
    "characters.xml": 0,
}


class ParseTests(unittest.TestCase):
    """Tests for proper parsing of example NeXML files."""

    def test_parse(self):
        """Extract and count phylogenetic trees using Phylo.parse."""
        for filename in nexml_files:
            count = tree_counts.get(filename, 1)
            path = os.path.join("NeXML", filename)
            msg = f"Failed parser test for {path}"
            trees = list(Phylo.parse(path, "nexml"))
            self.assertEqual(len(trees), count, msg=msg)


class WriterTests(unittest.TestCase):
    """NeXML writer tests."""

    def check(self, path):
        """Parse, rewrite and retest an example phylogeny file."""
        msg = f"Failed NeXMLIO writer test for {path}"

        # Using Phylo.NeXMLIO directly, binary mode handles
        with open(path, "rb") as stream:
            t1 = next(Phylo.NeXMLIO.Parser(stream).parse())
        stream = BytesIO()
        Phylo.NeXMLIO.write([t1], stream)
        stream.seek(0)
        t2 = next(Phylo.NeXMLIO.Parser(stream).parse())
        self.compare(t1, t2, msg)

        # Using Phylo.parse/write, using filenames
        msg = f"Failed Phylo API writer test for {path}"
        t1 = next(Phylo.parse(path, "nexml"))
        tmp = tempfile.NamedTemporaryFile().name
        Phylo.write([t1], tmp, "nexml")
        t2 = next(Phylo.parse(tmp, "nexml"))
        self.compare(t1, t2, msg)
        os.remove(tmp)

    def compare(self, t1, t2, msg=None):
        """Compare two trees."""
        for prop_name in ("name", "branch_length", "confidence"):
            p1 = sorted(
                getattr(n, prop_name)
                for n in t1.get_terminals()
                if getattr(n, prop_name)
            )
            p2 = sorted(
                getattr(n, prop_name)
                for n in t2.get_terminals()
                if getattr(n, prop_name)
            )
            self.assertEqual(p1, p2, msg=msg)

    def test_write(self):
        """Test for serialization of objects to NeXML format."""
        for filename in nexml_files:
            count = tree_counts.get(filename, 1)
            if count > 0:
                path = os.path.join("NeXML", filename)
                self.check(path)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
