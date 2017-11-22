#!/usr/bin/env python
#
# Copyright 2014 by Leighton Pritchard.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for online functionality of the KGML modules."""

# Builtins
from __future__ import with_statement
import os
import unittest

import requires_internet
requires_internet.check()

# Biopython
from Bio.Graphics.ColorSpiral import ColorSpiral

# Do we have ReportLab?  Raise error if not present.
from Bio import MissingExternalDependencyError
try:
    from reportlab.pdfgen.canvas import Canvas
    from reportlab.lib.pagesizes import A4
except ImportError:
    raise MissingExternalDependencyError(
        "Install reportlab if you want to use Bio.Graphics.")

# Do we have PIL?
try:
    from PIL import Image
except ImportError:
    raise MissingExternalDependencyError(
        "Install Pillow or its predecessor PIL (Python Imaging Library) "
        "if you want to use bitmaps from KGML.")


# Biopython Bio.KEGG.KGML
from Bio.KEGG.KGML.KGML_parser import read
from Bio.Graphics.KGML_vis import KGMLCanvas

# test_KGML_graphics module
from test_KGML_graphics import PathwayData


class KGMLPathwayOnlineTest(unittest.TestCase):
    """Import XML file and write KGML - online tests.

    Import metabolic maps from a local .xml KGML file, and from
    the KEGG site, and write valid KGML output for each
    """
    def setUp(self):
        # Does our output directory exist?  If not, create it
        if not os.path.isdir('KEGG'):
            os.mkdir('KEGG')
        # Define some data to work with as a list of tuples:
        # (infilename, outfilename, (entry_count, ortholog_count,
        # compound_count, map_counts), pathway_image,
        # show_image_map)
        self.data = [
            PathwayData("01100", (3628, 1726, 1746, 149)),
            PathwayData("03070", (81, 72, 8, 1), True),
            ]

    def test_render_KGML_import_map(self):
        """Basic rendering of KGML: use imported imagemap

        Uses the URL indicated in the .xml file.

        This test may fail if the imagemap is not available (e.g. if
        there is not a web connection), and may look odd if the remote
        imagemap has changed since the local KGML file was downloaded.
        """
        # We test rendering of the original KEGG KGML using imported files
        for p in self.data:
            with open(p.infilename, 'rU') as f:
                pathway = read(f)
                kgml_map = KGMLCanvas(pathway, import_imagemap=True)
                kgml_map.draw(p.output_stem + '_importmap.pdf')


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
