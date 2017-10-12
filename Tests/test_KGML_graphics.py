#!/usr/bin/env python
#
# Copyright 2013, 2014 by Leighton Pritchard.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for general functionality of the KGML modules."""

# Builtins
from __future__ import with_statement
import os
import unittest

# Biopython
from Bio.Graphics.ColorSpiral import ColorSpiral

# Do we have ReportLab?  Raise error if not present.
from Bio import MissingExternalDependencyError
try:
    # Not actually using these imports directly:
    from reportlab.pdfgen.canvas import Canvas
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.colors import HexColor
except ImportError:
    raise MissingExternalDependencyError(
        "Install reportlab if you want to use Bio.Graphics.")

try:
    c = HexColor('#8080F780')
except TypeError:
    # Known to fail under ReportLab 2.6 with:
    # unsupported operand type(s) for &: 'int' and 'float'
    # ReportLab 2.7+ also offers hasAlpha=True rather than alpha=True
    raise MissingExternalDependencyError(
        "Install at least reportlab 2.7 for transparency support.")

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


# The PathwayData class is also imported by the online test module in
# test_KGML_graphics_online.py
class PathwayData(object):
    """Convenience structure for testing pathway data."""
    def __init__(self, name, element_counts, show_pathway_image=False):
        """Initialize."""
        self.infilename = os.path.join("KEGG", "ko%s.xml" % name)
        self.outfilename = os.path.join("KEGG", "ko%s.kgml" % name)
        self.element_counts = element_counts
        self.pathway_image = os.path.join("KEGG", "map%s.png" % name)
        self.show_pathway_image = show_pathway_image
        self.output_stem = "Graphics/map%s" % name


class KGMLPathwayTest(unittest.TestCase):
    """Import XML file and write KGML.

    Import the ko01100 metabolic map from a local .xml KGML file,
    and write valid KGML output for each.
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
        # A list of KO IDs that we're going to use to modify pathway
        # appearance. These are KO IDs for reactions that take part in ko00020,
        # the TCA cycle
        self.ko_ids = \
            set(['ko:K00239', 'ko:K00240', 'ko:K00241', 'ko:K00242',
                 'ko:K00244', 'ko:K00245', 'ko:K00246', 'ko:K00247',
                 'ko:K00174', 'ko:K00175', 'ko:K00177', 'ko:K00176',
                 'ko:K00382', 'ko:K00164', 'ko:K00164', 'ko:K00658',
                 'ko:K01902', 'ko:K01903', 'ko:K01899', 'ko:K01900',
                 'ko:K01899', 'ko:K01900', 'ko:K00031', 'ko:K00030',
                 'ko:K00031', 'ko:K01648', 'ko:K00234', 'ko:K00235',
                 'ko:K00236', 'ko:K00237', 'ko:K01676', 'ko:K01677',
                 'ko:K01678', 'ko:K01679', 'ko:K01681', 'ko:K01682',
                 'ko:K01681', 'ko:K01682', 'ko:K01647', 'ko:K00025',
                 'ko:K00026', 'ko:K00024', 'ko:K01958', 'ko:K01959',
                 'ko:K01960', 'ko:K00163', 'ko:K00161', 'ko:K00162',
                 'ko:K00163', 'ko:K00161', 'ko:K00162', 'ko:K00382',
                 'ko:K00627', 'ko:K00169', 'ko:K00170', 'ko:K00172',
                 'ko:K00171', 'ko:K01643', 'ko:K01644', 'ko:K01646',
                 'ko:K01610', 'ko:K01596'])

    def test_render_KGML_basic(self):
        """Basic rendering of KGML: write to PDF without modification."""
        # We test rendering of the original KEGG KGML using only local
        # files.
        for p in self.data:
            with open(p.infilename, 'rU') as f:
                pathway = read(f)
                pathway.image = p.pathway_image
                kgml_map = KGMLCanvas(pathway)
                kgml_map.import_imagemap = p.show_pathway_image
                kgml_map.draw(p.output_stem + '_original.pdf')

    def test_render_KGML_modify(self):
        """Rendering of KGML to PDF, with modification."""
        # We test rendering of the original KGML for KO01100,
        # modifying line width for the lipid pathway
        p = self.data
        with open(p[0].infilename) as f:
            pathway = read(f)
            mod_rs = [e for e in pathway.orthologs if
                      len(set(e.name.split()).intersection(self.ko_ids))]
            for r in mod_rs:
                for g in r.graphics:
                    g.width = 10
            kgml_map = KGMLCanvas(pathway)
            kgml_map.draw(p[0].output_stem + '_widths.pdf')
        # We test rendering of the original KGML for KO3070,
        # modifying the reaction colours for each ortholog entry
        with open(p[1].infilename) as f:
            pathway = read(f)
            orthologs = [e for e in pathway.orthologs]
            # Use Biopython's ColorSpiral to generate colours
            cs = ColorSpiral(a=2, b=0.2, v_init=0.85, v_final=0.5,
                             jitter=0.03)
            colors = cs.get_colors(len(orthologs))
            for o, c in zip(orthologs, colors):
                for g in o.graphics:
                    g.bgcolor = c
            kgml_map = KGMLCanvas(pathway)
            pathway.image = p[1].pathway_image
            kgml_map.import_imagemap = p[1].show_pathway_image
            kgml_map.draw(p[1].output_stem + '_colors.pdf')

    def test_render_KGML_transparency(self):
        """Rendering of KGML to PDF, with color alpha channel."""
        # We test rendering of the original KGML for KO01100,
        # modifying alpha channel for the lipid pathway
        p = self.data
        with open(p[0].infilename) as f:
            pathway = read(f)
            mod_rs = [e for e in pathway.orthologs if
                      len(set(e.name.split()).intersection(self.ko_ids))]
            for r in mod_rs:
                for g in r.graphics:
                    # Modify hex colour directly by appending alpha channel
                    # to hex string
                    g.fgcolor = g.fgcolor + "77"
                    g.width = 20
            kgml_map = KGMLCanvas(pathway)
            kgml_map.draw(p[0].output_stem + '_transparency.pdf')
        # We test rendering of the original KGML for KO3070,
        # modifying the alpha channel for each ortholog entry
        with open(p[1].infilename) as f:
            pathway = read(f)
            orthologs = [e for e in pathway.orthologs]
            # Use Biopython's ColorSpiral to generate colours
            cs = ColorSpiral(a=2, b=0.2, v_init=0.85, v_final=0.5,
                             jitter=0.03)
            colors = cs.get_colors(len(orthologs))
            for o, c in zip(orthologs, colors):
                # Modify color tuples to add alpha channel
                c = c + (0.5, )
                for g in o.graphics:
                    g.bgcolor = c
            kgml_map = KGMLCanvas(pathway)
            pathway.image = p[1].pathway_image
            kgml_map.import_imagemap = p[1].show_pathway_image
            kgml_map.draw(p[1].output_stem + '_transparency.pdf')


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
