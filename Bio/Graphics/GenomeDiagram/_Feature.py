# Copyright 2003-2008 by Leighton Pritchard.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
#
# Contact:       Leighton Pritchard, The James Hutton Institute,
#                Invergowrie, Dundee, Scotland, DD2 5DA, UK
#                Leighton.Pritchard@hutton.ac.uk
################################################################################

"""Feature module.

Provides:
 - Feature - class to wrap Bio.SeqFeature objects with drawing information

For drawing capabilities, this module uses reportlab to define colors:
http://www.reportlab.com
"""

# ReportLab imports
from reportlab.lib import colors

# GenomeDiagram imports
from ._Colors import ColorTranslator


class Feature:
    """Class to wrap Bio.SeqFeature objects for GenomeDiagram.

    Attributes:
     - parent    FeatureSet, container for the object
     - id        Unique id
     - color    color.Color, color to draw the feature
     - hide      Boolean for whether the feature will be drawn or not
     - sigil     String denoting the type of sigil to use for the feature.
       Currently either "BOX" or "ARROW" are supported.
     - arrowhead_length  Float denoting length of the arrow head to be drawn,
       relative to the bounding box height.  The arrow shaft
       takes up the remainder of the bounding box's length.
     - arrowshaft_height  Float denoting length of the representative arrow
       shaft to be drawn, relative to the bounding box height.
       The arrow head takes the full height of the bound box.
     - name_qualifiers   List of Strings, describes the qualifiers that may
       contain feature names in the wrapped Bio.SeqFeature object
     - label     Boolean, 1 if the label should be shown
     - label_font    String describing the font to use for the feature label
     - label_size    Int describing the feature label font size
     - label_color  color.Color describing the feature label color
     - label_angle   Float describing the angle through which to rotate the
       feature label in degrees (default = 45, linear only)
     - label_position    String, 'start', 'end' or 'middle' denoting where
       to place the feature label. Leave as None for the default
       which is 'start' for linear diagrams, and at the bottom of
       the feature as drawn on circular diagrams.
     - label_strand  Integer -1 or +1 to explicitly place the label on the
       forward or reverse strand. Default (None) follows th
       feature's strand. Use -1 to put labels under (linear) or
       inside (circular) the track, +1 to put them above (linear)
       or outside (circular) the track.
     - locations     List of tuples of (start, end) ints describing where the
       feature and any subfeatures start and end
     - type      String denoting the feature type
     - name      String denoting the feature name
     - strand    Int describing the strand on which the feature is found

    """

    def __init__(
        self,
        parent=None,
        feature_id=None,
        feature=None,
        color=colors.lightgreen,
        label=0,
        border=None,
        colour=None,
    ):
        """Initialize.

        Arguments:
         - parent    FeatureSet containing the feature
         - feature_id    Unique id for the feature
         - feature   Bio.SeqFeature object to be wrapped
         - color    color.Color Color to draw the feature (overridden
           by backwards compatible argument with UK spelling, colour).
           Either argument is overridden if 'color' is found in feature
           qualifiers
         - border   color.Color Color to draw the feature border, use
           None for the same as the fill color, False for no border.
         - label     Boolean, 1 if the label should be shown

        """
        # Let the UK spelling (colour) override the USA spelling (color)
        if colour is not None:
            color = colour

        self._colortranslator = ColorTranslator()

        # Initialize attributes
        self.parent = parent
        self.id = feature_id
        self.color = color  # default color to draw the feature
        self.border = border
        self._feature = None  # Bio.SeqFeature object to wrap
        self.hide = 0  # show by default
        self.sigil = "BOX"
        self.arrowhead_length = 0.5  # 50% of the box height
        self.arrowshaft_height = 0.4  # 40% of the box height
        self.name_qualifiers = ["gene", "label", "name", "locus_tag", "product"]
        self.label = label
        self.label_font = "Helvetica"
        self.label_size = 6
        self.label_color = colors.black
        self.label_angle = 45
        self.label_position = None  # Expect 'start', 'middle', or 'end' (plus aliases)
        self.label_strand = None  # Expect +1 or -1 if overriding this

        if feature is not None:
            self.set_feature(feature)

    def set_feature(self, feature):
        """Define the Bio.SeqFeature object to be wrapped."""
        self._feature = feature
        self.__process_feature()

    def __process_feature(self):
        """Examine wrapped feature and set some properties accordingly (PRIVATE)."""
        self.locations = []
        bounds = []
        # This will be a list of length one for simple FeatureLocation:
        for location in self._feature.location.parts:
            start = location.nofuzzy_start
            end = location.nofuzzy_end
            # if start > end and self.strand == -1:
            #    start, end = end, start
            self.locations.append((start, end))
            bounds += [start, end]
        self.type = str(self._feature.type)  # Feature type
        # TODO - Strand can vary with subfeatures (e.g. mixed strand tRNA)
        if self._feature.strand is None:
            # This is the SeqFeature default (None), but the drawing code
            # only expects 0, +1 or -1.
            self.strand = 0
        else:
            self.strand = int(self._feature.strand)  # Feature strand
        if "color" in self._feature.qualifiers:  # Artemis color (if present)
            self.color = self._colortranslator.artemis_color(
                self._feature.qualifiers["color"][0]
            )
        self.name = self.type
        for qualifier in self.name_qualifiers:
            if qualifier in self._feature.qualifiers:
                self.name = self._feature.qualifiers[qualifier][0]
                break
        # Note will be 0 to N for origin wrapping feature on genome of length N
        self.start, self.end = min(bounds), max(bounds)

    def get_feature(self):
        """Return the unwrapped Bio.SeqFeature object."""
        return self._feature

    def set_colour(self, colour):
        """Backwards compatible variant of set_color(self, color) using UK spelling."""
        color = self._colortranslator.translate(colour)
        self.color = color

    def set_color(self, color):
        """Set the color in which the feature will be drawn.

        Arguments:
         - color    The color to draw the feature - either a colors.Color
           object, an RGB tuple of floats, or an integer corresponding a
           colors in colors.txt

        """
        # TODO - Make this into the set method for a color property?
        color = self._colortranslator.translate(color)
        self.color = color

    def __getattr__(self, name):
        """Get attribute by name.

        If the Feature class doesn't have the attribute called for,
        check in self._feature for it.
        """
        return getattr(self._feature, name)  # try to get the attribute from the feature


################################################################################
# RUN AS SCRIPT
################################################################################

if __name__ == "__main__":

    # Test code
    gdf = Feature()
