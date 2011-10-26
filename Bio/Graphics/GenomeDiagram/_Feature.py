# Copyright 2003-2008 by Leighton Pritchard.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
# Contact:       Leighton Pritchard, Scottish Crop Research Institute,
#                Invergowrie, Dundee, Scotland, DD2 5DA, UK
#                L.Pritchard@scri.ac.uk
################################################################################

""" Feature module

    Provides:

    o Feature - class to wrap Bio.SeqFeature objects with drawing information

    For drawing capabilities, this module uses reportlab to define colors:

    http://www.reportlab.com

    For dealing with biological information, the package uses BioPython:

    http://www.biopython.org
"""

# ReportLab imports
from reportlab.lib import colors

# GenomeDiagram imports
from _Colors import ColorTranslator

import string

class Feature(object):
    """ Class to wrap Bio.SeqFeature objects for GenomeDiagram

        Provides:

        Methods:

        o __init__(self, parent=None, feature_id=None, feature=None,
                 color=colors.lightgreen) Called when the feature is
                 instantiated

        o set_feature(self, feature) Wrap the passed feature

        o get_feature(self) Return the unwrapped Bio.SeqFeature object

        o set_color(self, color) Set the color in which the feature will
                be drawn (accepts multiple formats: reportlab color.Color()
                tuple and color.name, or integer representing Artemis color

        o get_color(self) Returns color.Color tuple of the feature's color

        o __getattr__(self, name) Catches attribute requests and passes them to
                the wrapped Bio.SeqFeature object

        Attributes:

        o parent    FeatureSet, container for the object

        o id        Unique id

        o color    color.Color, color to draw the feature

        o hide      Boolean for whether the feature will be drawn or not

        o sigil     String denoting the type of sigil to use for the feature.
                    Currently either "BOX" or "ARROW" are supported.

        o arrowhead_length  Float denoting length of the arrow head to be drawn,
                            relative to the bounding box height.  The arrow shaft
                            takes up the remainder of the bounding box's length.

        o arrowshaft_height  Float denoting length of the representative arrow
                             shaft to be drawn, relative to the bounding box height.
                             The arrow head takes the full height of the bound box.
         
        o name_qualifiers   List of Strings, describes the qualifiers that may
                    contain feature names in the wrapped Bio.SeqFeature object

        o label     Boolean, 1 if the label should be shown

        o label_font    String describing the font to use for the feature label

        o label_size    Int describing the feature label font size

        o label_color  color.Color describing the feature label color

        o label_angle   Float describing the angle through which to rotate the
                    feature label in degrees (default = 45, linear only)

        o label_position    String, 'start', 'end' or 'middle' denoting where
                    to place the feature label (linear only)

        o locations     List of tuples of (start, end) ints describing where the
                    feature and any subfeatures start and end

        o type      String denoting the feature type

        o name      String denoting the feature name

        o strand    Int describing the strand on which the feature is found

    """
    def __init__(self, parent=None, feature_id=None, feature=None,
                 color=colors.lightgreen, label=0, border=None, colour=None):
        """ __init__(self, parent=None, feature_id=None, feature=None,
                 color=colors.lightgreen, label=0)

            o parent    FeatureSet containing the feature

            o feature_id    Unique id for the feature

            o feature   Bio.SeqFeature object to be wrapped

            o color    color.Color Color to draw the feature (overridden
                       by backwards compatible argument with UK spelling,
                       colour).  Either argument is overridden if 'color'
                       is found in feature qualifiers

            o border   color.Color Color to draw the feature border, use
                       None for the same as the fill color, False for no border.

            o label     Boolean, 1 if the label should be shown
        """
        #Let the UK spelling (colour) override the USA spelling (color)
        if colour is not None:
            color = colour

        self._colortranslator = ColorTranslator()
        
        # Initialise attributes
        self.parent = parent
        self.id = feature_id        
        self.color = color            # default color to draw the feature
        self.border = border
        self._feature = None            # Bio.SeqFeature object to wrap
        self.hide = 0                   # show by default
        self.sigil = 'BOX'
        self.arrowhead_length = 0.5 # 50% of the box height
        self.arrowshaft_height = 0.4 # 40% of the box height
        self.name_qualifiers = ['gene', 'label', 'name', 'locus_tag', 'product']
        self.label = label
        self.label_font = 'Helvetica'
        self.label_size = 6
        self.label_color = colors.black
        self.label_angle = 45
        self.label_position = 'start'
        
        if feature is not None:
            self.set_feature(feature)

    def set_feature(self, feature):
        """ set_feature(self, feature)

            o feature   Bio.SeqFeature object to be wrapped

            Defines the Bio.SeqFeature object to be wrapped
        """
        self._feature = feature
        self.__process_feature()


    def __process_feature(self):
        """ __process_feature(self)

            Examine the feature to be wrapped, and set some of the Feature's
            properties accordingly
        """
        self.locations = []
        bounds = []
        if self._feature.sub_features == []:
            start = self._feature.location.nofuzzy_start
            end = self._feature.location.nofuzzy_end
            #if start > end and self.strand == -1:
            #    start, end = end, start
            self.locations.append((start, end))
            bounds += [start, end]
        else:
            for subfeature in self._feature.sub_features:
                start = subfeature.location.nofuzzy_start
                end = subfeature.location.nofuzzy_end
                #if start > end and self.strand == -1:
                #    start, end = end, start
                self.locations.append((start, end))                
                bounds += [start, end]
        self.type = str(self._feature.type)                     # Feature type
        #TODO - Strand can vary with subfeatures (e.g. mixed strand tRNA)
        if self._feature.strand is None:
            #This is the SeqFeature default (None), but the drawing code
            #only expects 0, +1 or -1.
            self.strand = 0
        else:
            self.strand = int(self._feature.strand)                 # Feature strand
        if 'color' in self._feature.qualifiers:                # Artemis color (if present)
            self.color = self._colortranslator.artemis_color( \
                                         self._feature.qualifiers['color'][0])
        self.name = self.type
        for qualifier in self.name_qualifiers:            
            if qualifier in self._feature.qualifiers:
                self.name = self._feature.qualifiers[qualifier][0]
                break
        #Note will be 0 to N for origin wrapping feature on genome of length N 
        self.start, self.end = min(bounds), max(bounds)


    def get_feature(self):
        """ get_feature(self) -> Bio.SeqFeature

            Returns the unwrapped Bio.SeqFeature object
        """
        return self._feature

    def set_colour(self, colour):
        """Backwards compatible variant of set_color(self, color) using UK spelling."""
        color = self._colortranslator.translate(colour)
        self.color = color

    def set_color(self, color):
        """ set_color(self, color)

            o color    The color to draw the feature - either a colors.Color
                       object, an RGB tuple of floats, or an integer
                       corresponding to colors in colors.txt
                           
            Set the color in which the feature will be drawn
        """
        #TODO - Make this into the set method for a color property?
        color = self._colortranslator.translate(color)
        self.color = color

    def __getattr__(self, name):
        """ __getattr__(self, name) -> various

            If the Feature class doesn't have the attribute called for,
            check in self._feature for it
        """
        return getattr(self._feature, name) # try to get the attribute from the feature


    
################################################################################
# RUN AS SCRIPT
################################################################################

if __name__ == '__main__':

    # Test code
    gdf = Feature()
