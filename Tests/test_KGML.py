#!/usr/bin/env python
#
# Copyright 2013 by Leighton Pritchard.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

""" Tests for general functionality of the KGML parser, pathway and 
    visualisation modules
"""

# Builtins
import os
import unittest

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

# Do we have internet (code from Biopython's requires_internet.py)
from Bio import MissingExternalDependencyError

def check_internet():
    try:
        check.available
    except AttributeError:
        # I'm going to check for internet availability
        RELIABLE_DOMAIN = "biopython.org"
        import socket
        try:
            socket.getaddrinfo(RELIABLE_DOMAIN,
                               80,
                               socket.AF_UNSPEC,
                               socket.SOCK_STREAM)
        except socket.gaierror, x:
            check.available = False
        else:
            check.available = True
    if not check.available:
        raise MissingExternalDependencyError("internet not available")

# Biopython Bio.KEGG.KGML (?)
from Bio.KGML.KGML_parser import read
from Bio.Graphics.KGML_vis import KGMLCanvas
from Bio.KGML.KGML_scrape import retrieve_kgml_to_file, retrieve_KEGG_pathway

class PathwayData(object):
    """ Convenience structure for testing pathway data
    """
    def __init__(self, infilename, outfilename, element_counts,
                 pathway_image, show_pathway_image=False):
        self.infilename = infilename
        self.outfilename = outfilename
        self.element_counts = element_counts
        self.pathway_image = pathway_image
        self.show_pathway_image = show_pathway_image

class KGMLPathwayTest(unittest.TestCase):
    """ Import the ko01100 metabolic map from a local .xml KGML file, and from
        the KEGG site, and write valid KGML output for each
    """
    def setUp(self):
        # Does our output director exist?  If not, create it
        if not os.path.isdir('KEGG'):
            os.mkdir('KEGG')
        # Define some data to work with as a list of tuples:
        # (infilename, outfilename, (entry_count, ortholog_count,
        # compound_count, map_counts), pathway_image,
        # show_image_map)
        self.data = [PathwayData(os.path.join("KEGG", "ko01100.xml"),
                                 os.path.join("KEGG", "ko01100.kgml"),
                                 (3628, 1726, 1746, 149),
                                 os.path.join("KEGG", "map01100.png")),
                     PathwayData(os.path.join("KEGG", "ko03070.xml"),
                                 os.path.join("KEGG", "ko03070.kgml"),
                                 (81, 72, 8, 1),
                                 os.path.join("KEGG", "map03070.png"),
                                 True)]
        # A list of KO IDs that we're going to use to modify pathway 
        # appearance. These are KO IDs for reactions that take part in ko00020, 
        # the TCA cycle
        self.ko_ids = \
            set(['ko:K00239','ko:K00240','ko:K00241','ko:K00242','ko:K00244',
                 'ko:K00245','ko:K00246','ko:K00247','ko:K00174','ko:K00175',
                 'ko:K00177','ko:K00176','ko:K00382','ko:K00164','ko:K00164',
                 'ko:K00658','ko:K01902','ko:K01903','ko:K01899','ko:K01900',
                 'ko:K01899','ko:K01900','ko:K00031','ko:K00030','ko:K00031',
                 'ko:K01648','ko:K00234','ko:K00235','ko:K00236','ko:K00237',
                 'ko:K01676','ko:K01677','ko:K01678','ko:K01679','ko:K01681',
                 'ko:K01682','ko:K01681','ko:K01682','ko:K01647','ko:K00025',
                 'ko:K00026','ko:K00024','ko:K01958','ko:K01959','ko:K01960',
                 'ko:K00163','ko:K00161','ko:K00162','ko:K00163','ko:K00161',
                 'ko:K00162','ko:K00382','ko:K00627','ko:K00169','ko:K00170',
                 'ko:K00172','ko:K00171','ko:K01643','ko:K01644','ko:K01646',
                 'ko:K01610','ko:K01596'])


    def test_read_and_write_KGML_files(self):
        """ Read KGML from, and write KGML to, local files.
            Check we read/write the correct number of elements.
        """
        for p in self.data:
            # Test opening file
            with open(p.infilename, 'rU') as f:
                pathway = read(f)
                # Do we have the correct number of elements of each type
                self.assertEqual((len(pathway.entries), 
                                  len(pathway.orthologs),
                                  len(pathway.compounds),
                                  len(pathway.maps)),
                                 p.element_counts)
            # Test writing file
            with open(p.outfilename, 'w') as f:
                f.write(pathway.get_KGML())
            # Can we read the file we wrote?
            with open(p.outfilename, 'rU') as f:
                pathway = read(f)
                # Do we have the correct number of elements of each type
                self.assertEqual((len(pathway.entries), 
                                  len(pathway.orthologs),
                                  len(pathway.compounds),
                                  len(pathway.maps)),
                                 p.element_counts)


    def test_render_KGML_basic(self):
        """ Basic rendering of KGML: write to PDF without modification.
        """
        # We test rendering of both the original KEGG KGML and the 
        # roundtrip KGML written by this module, using only local files.
        for p in self.data:
            for filename, ext in [(p.infilename, '_original.pdf'),
                                  (p.outfilename, '_roundtrip.pdf')]:
                with open(filename, 'rU') as f:
                    pathway = read(f)
                    pathway.image = p.pathway_image
                    kgml_map = KGMLCanvas(pathway)
                    kgml_map.import_imagemap = p.show_pathway_image
                    kgml_map.draw(os.path.splitext(filename)[0] + ext)
        

    def test_render_KGML_modify(self):
        """ Rendering of KGML to PDF, with modification.
        """
        # We test rendering of the original and roundtrip KGML for KO01100,
        # modifying line width for the lipid pathway
        p = self.data[0]
        for filename, ext in [(p.infilename, '_mod_original.pdf'),
                              (p.outfilename, '_mod_roundtrip.pdf')]:
            with open(filename) as f:
                pathway = read(f)
                mod_rs = [e for e in pathway.orthologs if \
                        len(set(e.name.split()).intersection(self.ko_ids))]
                for r in mod_rs:
                    for g in r.graphics:
                        g.width = 10
                kgml_map = KGMLCanvas(pathway)
                kgml_map.draw(os.path.splitext(filename)[0] + ext)
        # We test rendering of the original and roundtrip KGML for KO3070,
        # modifying the reaction colours for each ortholog entry
        p = self.data[1]
        for filename, ext in [(p.infilename, '_mod_original.pdf'),
                              (p.outfilename, '_mod_roundtrip.pdf')]:
            with open(filename) as f:
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
                pathway.image = p.pathway_image
                kgml_map.import_imagemap = p.show_pathway_image
                kgml_map.draw(os.path.splitext(filename)[0] + ext)

    def test_KEGG_download(self):
        """ Download a KEGG pathway from the KEGG server and write KGML.
        """
        # Download the KEGG ko01110 pathway, and render it
        pathway = retrieve_KEGG_pathway('ko01110')
        kgml_map = KGMLCanvas(pathway)
        kgml_map.import_imagemap = True
        kgml_map.draw(os.path.join('KEGG', 'ko01110.pdf'))
        # Download the KEGG ko01120 pathway and write to file as KGML
        retrieve_kgml_to_file("ko01120", os.path.join("KEGG", "ko01120.xml"))
        


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner = runner)
