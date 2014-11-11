""" This module provides classes and functions to visualise a KGML Pathway Map

The KGML definition is as of release KGML v0.7.1
(http://www.kegg.jp/kegg/xml/docs/)

Classes:
"""

from __future__ import print_function

import os
import tempfile
from itertools import chain
from io import BytesIO

from reportlab.lib import pagesizes
from reportlab.lib import colors
from reportlab.lib.utils import ImageReader
from reportlab.graphics.shapes import *
from reportlab.pdfgen import canvas

from PIL import Image

from Bio._py3k import urlopen as _urlopen

from Bio.KEGG.KGML.KGML_pathway import Pathway


def hexdarken(hexcolor, factor=0.7):
    """Returns darkened hex color as a ReportLab RGB color.

    Take a passed hex color and return an RGB color that is
    slightly darker (if possible).
    """
    c = colors.HexColor(hexcolor)
    for a in ['red', 'green', 'blue']:
        setattr(c, a, factor * getattr(c, a))
    return c


def get_temp_imagefilename(url):
    """Returns filename of temporary file containing downloaded image.

    Create a new temporary file to hold the image file at the passed URL
    and return the filename.
    """
    img = _urlopen(url).read()
    im = Image.open(BtyesIO(img))
    # im.transpose(Image.FLIP_TOP_BOTTOM)
    f = tempfile.NamedTemporaryFile(delete=False, suffix='.png')
    fname = f.name
    f.close()
    im.save(fname, 'PNG')
    return fname


class KGMLCanvas(object):
    """Reportlab Canvas-based representation of a KGML pathway map."""

    def __init__(self, pathway, import_imagemap=False, label_compounds=True,
                 label_orthologs=True, label_reaction_entries=True,
                 label_maps=True, show_maps=False, fontname='Helvetica',
                 fontsize=6, draw_relations=True, show_orthologs=True,
                 show_compounds=True, show_genes=True,
                 margins=(0.02, 0.02)):
        self.pathway = pathway
        self.show_maps = show_maps
        self.show_orthologs = show_orthologs
        self.show_compounds = show_compounds
        self.show_genes = show_genes
        self.label_compounds = label_compounds
        self.label_orthologs = label_orthologs
        self.label_reaction_entries = label_reaction_entries
        self.label_maps = label_maps
        self.fontname = fontname
        self.fontsize = fontsize
        self.draw_relations = draw_relations
        self.non_reactant_transparency = 0.3
        self.import_imagemap = import_imagemap  # Import the map .png from URL
        # percentage of canvas that will be margin in on either side in the
        # X and Y directions
        self.margins = margins

    def draw(self, filename):
        """Add the map elements to the drawing."""
        # Instantiate the drawing, first
        # size x_max, y_max for now - we can add margins, later
        if self.import_imagemap:
            # We're drawing directly on the image, so we set the canvas to the
            # same size as the image
            if os.path.isfile(self.pathway.image):
                imfilename = self.pathway.image
            else:
                imfilename = get_temp_imagefilename(self.pathway.image)
            im = Image.open(imfilename)
            cwidth, cheight = im.size
        else:
            # No image, so we set the canvas size to accommodate visible
            # elements
            cwidth, cheight = (self.pathway.bounds[1][0],
                               self.pathway.bounds[1][1])
        # Instantiate canvas
        self.drawing = \
            canvas.Canvas(filename, bottomup=0,
                          pagesize=(cwidth *
                                        (1 + 2 * self.margins[0]),
                                    cheight *
                                        (1 + 2 * self.margins[1])))
        self.drawing.setFont(self.fontname, self.fontsize)
        # Transform the canvas to add the margins
        self.drawing.translate(self.margins[0] * self.pathway.bounds[1][0],
                               self.margins[1] * self.pathway.bounds[1][1])
        # Add the map image, if required
        if self.import_imagemap:
            self.drawing.saveState()
            self.drawing.scale(1, -1)
            self.drawing.translate(0, -cheight)
            self.drawing.drawImage(imfilename, 0, 0)
            self.drawing.restoreState()
        # Add the reactions, compounds and maps
        # Maps go on first, to be overlaid by more information.
        # By default, they're slightly transparent.
        if self.show_maps:
            self.__add_maps()
        self.__add_reaction_entries()
        if self.show_orthologs:
            self.__add_orthologs()
        if self.show_compounds:
            self.__add_compounds()
        if self.show_genes:
            self.__add_genes()
        # TODO: complete draw_relations code
        # if self.draw_relations:
        #    self.__add_relations()
        # Write the pathway map to PDF
        self.drawing.save()

    def __add_maps(self):
        """Adds maps to the drawing of the map.

        We do this first, as they're regional labels to be overlaid by
        information.  Also, we want to set the color to something subtle.

        We're using Hex colors because that's what KGML uses, and
        Reportlab doesn't mind.
        """
        for m in self.pathway.maps:
            for g in m.graphics:
                self.drawing.setStrokeColor('#888888')
                self.drawing.setFillColor('#DDDDDD')
                self.__add_graphics(g)
                if self.label_maps:
                    self.drawing.setFillColor('#888888')
                    self.__add_labels(g)

    def __add_graphics(self, graphics):
        """Adds the passed graphics object to the map.

        Add text, add after the graphics object, for sane Z-ordering.
        """
        if graphics.type == 'line':
            p = self.drawing.beginPath()
            x, y = graphics.coords[0]
            # There are optional settings for lines that aren't necessarily
            # part of the KGML DTD
            if graphics.width is not None:
                self.drawing.setLineWidth(graphics.width)
            else:
                self.drawing.setLineWidth(1)
            p.moveTo(x, y)
            for (x, y) in graphics.coords:
                p.lineTo(x, y)
            self.drawing.drawPath(p)
            self.drawing.setLineWidth(1)        # Return to default
        # KGML defines the (x, y) coordinates as the centre of the circle/
        # rectangle/roundrectangle, but Reportlab uses the co-ordinates of the
        # lower-left corner for rectangle/elif.
        if graphics.type == 'circle':
            self.drawing.circle(graphics.x, graphics.y, graphics.width*0.5,
                                stroke=1, fill=1)
        elif graphics.type == 'roundrectangle':
            self.drawing.roundRect(graphics.x - graphics.width * 0.5,
                                   graphics.y - graphics.height * 0.5,
                                   graphics.width, graphics.height,
                                   min(graphics.width, graphics.height) * 0.1,
                                   stroke=1, fill=1)
        elif graphics.type == 'rectangle':
            self.drawing.rect(graphics.x - graphics.width * 0.5,
                              graphics.y - graphics.height * 0.5,
                              graphics.width, graphics.height,
                              stroke=1, fill=1)

    def __add_labels(self, graphics):
        """Adds labels for the passed graphics objects to the map (PRIVATE).

        We don't check that the labels fit inside objects such as circles/
        rectangles/roundrectangles.
        """
        if graphics.type == 'line':
            # We use the midpoint of the line - sort of - we take the median
            # line segment (list-wise, not in terms of length), and use the
            # midpoint of that line.  We could have other options here,
            # maybe even parameterising it to a proportion of the total line
            # length.
            mid_idx = len(graphics.coords) * 0.5
            if not int(mid_idx) == mid_idx:
                idx1, idx2 = int(mid_idx - 0.5), int(mid_idx + 0.5)
            else:
                idx1, idx2 = int(mid_idx - 1), int(mid_idx)
            x1, y1 = graphics.coords[idx1]
            x2, y2 = graphics.coords[idx2]
            x, y = 0.5 * (x1 + x2), 0.5 * (y1 + y2)
        elif graphics.type == 'circle':
            x, y = graphics.x, graphics.y
        elif graphics.type in ('rectangle', 'roundrectangle'):
            x, y = graphics.x, graphics.y
        # How big so we want the text, and how many characters?
        if graphics._parent.type == 'map':
            text = graphics.name
            self.drawing.setFont(self.fontname, self.fontsize + 2)
        elif len(graphics.name) < 15:
            text = graphics.name
        else:
            text = graphics.name[:12] + '...'
        self.drawing.drawCentredString(x, y, text)
        self.drawing.setFont(self.fontname, self.fontsize)

    def __add_orthologs(self):
        """Adds 'ortholog' Entry elements to the drawing of the map (PRIVATE).

        In KGML, these are typically line objects, so we render them
        before the compound circles to cover the unsightly ends/junctions.
        """
        for ortholog in self.pathway.orthologs:
            for g in ortholog.graphics:
                self.drawing.setStrokeColor(g.fgcolor)
                self.drawing.setFillColor(g.bgcolor)
                self.__add_graphics(g)
                if self.label_orthologs:
                    # We want the label color to be slightly darker
                    # (where possible), so it can be read
                    self.drawing.setFillColor(hexdarken(g.fgcolor))
                    self.__add_labels(g)

    def __add_reaction_entries(self):
        """Adds Entry elements corresponding to Reactions to the map drawing (PRIVATE).

        In KGML, these are typically line objects, so we render them
        before the compound circles to cover the unsightly ends/junctions
        """
        for reaction in self.pathway.reaction_entries:
            for g in reaction.graphics:
                self.drawing.setStrokeColor(g.fgcolor)
                self.drawing.setFillColor(g.bgcolor)
                self.__add_graphics(g)
                if self.label_reaction_entries:
                    # We want the label color to be slightly darker
                    # (where possible), so it can be read
                    self.drawing.setFillColor(hexdarken(g.fgcolor))
                    self.__add_labels(g)

    def __add_compounds(self):
        """Adds compound elements to the drawing of the map (PRIVATE)."""
        for compound in self.pathway.compounds:
            for g in compound.graphics:
                # Modify transparency of compounds that don't participate
                # in reactions
                fillcolor = colors.HexColor(g.bgcolor)
                if not compound.is_reactant:
                    fillcolor.alpha *= self.non_reactant_transparency
                self.drawing.setStrokeColor(g.fgcolor)
                self.drawing.setFillColor(fillcolor)
                self.__add_graphics(g)
                if self.label_compounds:
                    if not compound.is_reactant:
                        t = 0.3
                    else:
                        t = 1
                    self.drawing.setFillColor(colors.Color(0.2, 0.2, 0.2, t))
                    self.__add_labels(g)

    def __add_genes(self):
        """Adds gene elements to the drawing of the map (PRIVATE)."""
        for gene in self.pathway.genes:
            for g in gene.graphics:
                fillcolor = colors.HexColor(g.bgcolor)
                self.drawing.setStrokeColor(g.fgcolor)
                self.drawing.setFillColor(fillcolor)
                self.__add_graphics(g)
                if self.label_compounds:
                    self.drawing.setFillColor(hexdarken(g.fgcolor))
                    self.__add_labels(g)

    def __add_relations(self):
        """Adds relations to the map (PRIVATE).

        This is tricky. There is no defined graphic in KGML for a
        relation, and the corresponding entries are typically defined
        as objects 'to be connected somehow'.  KEGG uses KegSketch, which
        is not public, and most third-party software draws straight line
        arrows, with heads to indicate the appropriate direction
        (at both ends for reversible reactions), using solid lines for
        ECrel relation types, and dashed lines for maplink relation types.

        The relation has:
        - entry1: 'from' node
        - entry2: 'to' node
        - subtype: what the relation refers to

        Typically we have entry1 = map/ortholog; entry2 = map/ortholog,
        subtype = compound.
        """
        # Dashed lines for maplinks, solid for everything else
        for relation in list(self.pathway.relations):
            if relation.type == 'maplink':
                self.drawing.setDash(6, 3)
            else:
                self.drawing.setDash()
            for s in relation.subtypes:
                subtype = self.pathway.entries[s[1]]
                # Our aim is to draw an arrow from the entry1 object to the
                # entry2 object, via the subtype object.
                # 1) Entry 1 to subtype
                self.__draw_arrow(relation.entry1, subtype)
                # 2) subtype to Entry 2
                self.__draw_arrow(subtype, relation.entry2)

    def __draw_arrow(self, g_from, g_to):
        """Draw an arrow between given Entry objects (PRIVATE).

        Draws an arrow from the g_from Entry object to the g_to
        Entry object; both must have Graphics objects.
        """
        # Centre and bound co-ordinates for the from and two objects
        bounds_from, bounds_to = g_from.bounds, g_to.bounds
        centre_from = (0.5 * (bounds_from[0][0] + bounds_from[1][0]),
                       0.5 * (bounds_from[0][1] + bounds_from[1][1]))
        centre_to = (0.5 * (bounds_to[0][0] + bounds_to[1][0]),
                     0.5 * (bounds_to[0][1] + bounds_to[1][1]))
        p = self.drawing.beginPath()
        # print(True, g_from.name, g_to.name, bounds_to, bounds_from)
        # If the 'from' and 'to' graphics are vertically-aligned, draw a line
        # from the 'from' to the 'to' entity
        if bounds_to[0][0] < centre_from[0] < bounds_to[1][0]:
            # print(True, g_from.name, g_to.name, bounds_to, bounds_from)
            if centre_to[1] > centre_from[1]:  # to above from
                p.moveTo(centre_from[0], bounds_from[1][1])
                p.lineTo(centre_from[0], bounds_to[0][1])
                # Draw arrow point - TODO
            else:                             # to below from
                p.moveTo(centre_from[0], bounds_from[0][1])
                p.lineTo(centre_from[0], bounds_to[1][1])
                # Draw arrow point - TODO
        elif bounds_from[0][0] < centre_to[0] < bounds_from[1][0]:
            # print(True, g_from.name, g_to.name, bounds_to, bounds_from)
            if centre_to[1] > centre_from[1]:  # to above from
                p.moveTo(centre_to[0], bounds_from[1][1])
                p.lineTo(centre_to[0], bounds_to[0][1])
                # Draw arrow point - TODO
            else:                             # to below from
                p.moveTo(centre_to[0], bounds_from[0][1])
                p.lineTo(centre_to[0], bounds_to[1][1])
                # Draw arrow point - TODO
        self.drawing.drawPath(p)    # Draw arrow shaft
        # print(g_from)
        # print(bounds_from)
        # print(g_to)
        # print(bounds_to)


if __name__ == '__main__':
    # Test production of Reportlab Canvas PDF visualisation
    # Try a default KO metabolic map with ortholog lines given
    pathway = KGML_parser.read(open('ko01100.xml', 'rU'))
    print(pathway)
    kgml_map = KGMLCanvas(pathway)
    kgml_map.draw('KGML_canvas_test.pdf')

    # Try a Dickeya metabolic map with ortholog lines, modifying reaction
    # graphics
    pathway = KGML_parser.read(open('ddc01100.xml', 'rU'))
    print(pathway)
    kgml_map = KGMLCanvas(pathway)
    # Set all reaction linewidths to 3 units
    for r in pathway.reaction_entries:
        for g in r.graphics:
            g.width = 3
    kgml_map.draw('KGML_canvas_ddc_test.pdf')

    # Try a KO metabolic map with no ortholog lines, using a local .png
    pathway = KGML_parser.read(open('ko_metabolic/ko00910.xml', 'rU'))
    pathway.image = 'map/map00910.png'
    print(pathway)
    kgml_map = KGMLCanvas(pathway)
    kgml_map.import_imagemap = True
    kgml_map.show_maps = False
    kgml_map.draw('KGML_canvas_map_local_test.pdf')

    # Try a KO metabolic map with no ortholog lines, using the KEGG .png
    pathway = KGML_parser.read(open('ko_metabolic/ko00253.xml', 'rU'))
    print(pathway)
    kgml_map = KGMLCanvas(pathway)
    kgml_map.import_imagemap = True
    kgml_map.show_maps = False
    kgml_map.draw('KGML_canvas_map_test.pdf')

    # Try a KO metabolic map with no ortholog lines using the KEGG .png,
    # but this time using a Dickeya XML file
    pathway = KGML_parser.read(open('dda00190.xml', 'rU'))
    print(pathway)
    kgml_map = KGMLCanvas(pathway)
    kgml_map.import_imagemap = True
    kgml_map.show_maps = True
    kgml_map.draw('KGML_canvas_dda_map_test.pdf')

    # Try a KO metabolic map with no ortholog lines using the KEGG .png,
    # but this time using a Dickeya XML file
    pathway = KGML_parser.read(open('test_retrieve_ddc00190.kgml', 'rU'))
    print(pathway)
    kgml_map = KGMLCanvas(pathway)
    kgml_map.import_imagemap = True
    kgml_map.show_maps = True
    kgml_map.draw('KGML_canvas_ddc_map_test.pdf')
