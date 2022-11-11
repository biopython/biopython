#
# Copyright 2013 by Leighton Pritchard.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for general functionality of the KGML parser and pathway model."""

# Builtins
import os
import unittest
import tempfile

# Biopython Bio.KEGG.KGML (?)
from Bio.KEGG.KGML.KGML_parser import read


class PathwayData:
    """Convenience structure for testing pathway data."""

    def __init__(
        self,
        infilename,
        outfilename,
        element_counts,
        pathway_image,
        show_pathway_image=False,
    ):
        """Initialize the class."""
        self.infilename = infilename
        self.outfilename = outfilename
        self.element_counts = element_counts
        self.pathway_image = pathway_image
        self.show_pathway_image = show_pathway_image


class KGMLPathwayTest(unittest.TestCase):
    """KGML checks using ko01100 metabolic map.

    Import the ko01100 metabolic map from a local .xml KGML file, and from
    the KEGG site, and write valid KGML output for each
    """

    def setUp(self):
        # Does our output director exist?  If not, create it
        if not os.path.isdir("KEGG"):
            os.mkdir("KEGG")
        # Define some data to work with as a list of tuples:
        # (infilename, outfilename, (entry_count, ortholog_count,
        # compound_count, map_counts), pathway_image,
        # show_image_map)
        self.data = [
            PathwayData(
                os.path.join("KEGG", "ko01100.xml"),
                tempfile.gettempprefix() + ".ko01100.kgml",
                (3628, 1726, 1746, 149),
                os.path.join("KEGG", "map01100.png"),
            ),
            PathwayData(
                os.path.join("KEGG", "ko03070.xml"),
                tempfile.gettempprefix() + ".ko03070.kgml",
                (81, 72, 8, 1),
                os.path.join("KEGG", "map03070.png"),
                True,
            ),
        ]
        # A list of KO IDs that we're going to use to modify pathway
        # appearance. These are KO IDs for reactions that take part in ko00020,
        # the TCA cycle
        # Turn black code style off
        # fmt: off
        self.ko_ids = {
            "ko:K00239", "ko:K00240", "ko:K00241", "ko:K00242", "ko:K00244",
            "ko:K00245", "ko:K00246", "ko:K00247", "ko:K00174", "ko:K00175",
            "ko:K00177", "ko:K00176", "ko:K00382", "ko:K00164", "ko:K00164",
            "ko:K00658", "ko:K01902", "ko:K01903", "ko:K01899", "ko:K01900",
            "ko:K01899", "ko:K01900", "ko:K00031", "ko:K00030", "ko:K00031",
            "ko:K01648", "ko:K00234", "ko:K00235", "ko:K00236", "ko:K00237",
            "ko:K01676", "ko:K01677", "ko:K01678", "ko:K01679", "ko:K01681",
            "ko:K01682", "ko:K01681", "ko:K01682", "ko:K01647", "ko:K00025",
            "ko:K00026", "ko:K00024", "ko:K01958", "ko:K01959", "ko:K01960",
            "ko:K00163", "ko:K00161", "ko:K00162", "ko:K00163", "ko:K00161",
            "ko:K00162", "ko:K00382", "ko:K00627", "ko:K00169", "ko:K00170",
            "ko:K00172", "ko:K00171", "ko:K01643", "ko:K01644", "ko:K01646",
            "ko:K01610", "ko:K01596"
        }
        # Turn black code style on
        # fmt: on

    def tearDown(self):
        for p in self.data:
            if os.path.isfile(p.outfilename):
                os.remove(p.outfilename)

    def test_read_and_write_KGML_files(self):
        """Read KGML from, and write KGML to, local files.

        Check we read/write the correct number of elements.
        """
        for p in self.data:
            # Test opening file
            with open(p.infilename) as f:
                pathway = read(f)
                # Do we have the correct number of elements of each type
                self.assertEqual(
                    (
                        len(pathway.entries),
                        len(pathway.orthologs),
                        len(pathway.compounds),
                        len(pathway.maps),
                    ),
                    p.element_counts,
                )
            # Test writing file
            with open(p.outfilename, "w") as f:
                f.write(pathway.get_KGML())
            # Can we read the file we wrote?
            with open(p.outfilename) as f:
                pathway = read(f)
                # Do we have the correct number of elements of each type
                self.assertEqual(
                    (
                        len(pathway.entries),
                        len(pathway.orthologs),
                        len(pathway.compounds),
                        len(pathway.maps),
                    ),
                    p.element_counts,
                )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
