# Copyright 2003-2008 by Leighton Pritchard.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
# Contact:       Leighton Pritchard, Scottish Crop Research Institute,
#                Invergowrie, Dundee, Scotland, DD2 5DA, UK
#                L.Pritchard@scri.ac.uk
################################################################################

""" Track module

    Provides:

    o Track - Container for a single track on the diagram, containing
                FeatureSet and GraphSet objects

    For drawing capabilities, this module uses reportlab to draw and write
    the diagram:

    http://www.reportlab.com

    For dealing with biological information, the package expects BioPython
    objects:

    http://www.biopython.org

"""

# ReportLab imports
from reportlab.lib import colors

# GenomeDiagram imports
from _FeatureSet import FeatureSet
from _GraphSet import GraphSet

class Track(object):
    """ Track

        Provides:

        Methods:

        o __init__(self, name=None, ...) Called on instantiation

        o add_set(self, set)    Add a FeatureSet or GraphSet to the diagram

        o del_set(self, set_id) Delete a FeatureSet or GraphSet from the
                                diagram

        o get_sets(self)    Returns a list of the sets in the track

        o get_ids(self)     Returns a list of the ids for sets in the track

        o range(self)   Returns the base/position range covered by the data in
                        the track

        o to_string(self, verbose=0)    Returns a string describing the data in
                                        the track

        o __getitem__(self, key)    Returns the set with the passed id

        o __str__(self) Returns a formatted string describing the track

        Attributes:

        o height    Int describing the relative height to other trackscale_fontsizes in the
                    diagram

        o name      String describing the track

        o hide      Boolean, 0 if the track is not to be drawn
        
        o start, end    Integers (or None) specifying start/end to draw just
                        a partial track.

        o greytrack     Boolean, 1 if a grey background to the track is to be
                        drawn

        o greytrack_labels  Int describing how many track-identifying labels
                            should be placed on the track at regular intervals

        o greytrack_font    String describing the font to use for the greytrack
                            labels

        o greytrack_fontsize    Int describing the font size to display the
                                labels on the grey track

        o greytrack_font_rotation   Int describing the angle through which to
                                    rotate the grey track labels

        o greytrack_font_color     colors.Color describing the color to draw
                                    the grey track labels

        o scale     Boolean, 1 if a scale is to be drawn on the track

        o scale_format  String, defaults to None, when scale values are written
                        as numerals.  Setting this to 'SInt' invokes SI
                        unit-like multiples, such as Mbp, Kbp and so on.

        o scale_color  colors.Color to draw the elemnts of the scale

        o scale_font    String describing the font to use for the scale labels

        o scale_fontsize    Int describing the size of the scale label font

        o scale_fontangle   Int describing the angle at which to draw the scale
                            labels (linear only)

        o scale_ticks       Boolean, 1 if ticks should be drawn at all on the
                            scale

        o scale_largeticks  Float (0->1) describing the height of large
                            scale ticks relative to the track height.

        o scale_smallticks  Float (0->1) describing the height of large
                            scale ticks relative to the track height.

        o scale_largetick_interval  Int, describing the number of bases that
                                    should separate large ticks

        o scale_smalltick_interval  Int, describing the number of bases that
                                    should separate small ticks

        o scale_largetick_labels    Boolean describing whether position labels
                                    should be written over large ticks
                                    
        o scale_smalltick_labels    Boolean describing whether position labels
                                    should be written over small ticks

        o axis_labels       Boolean describing whether the value labels should
                            be placed on the Y axes
    """
    def __init__(self, name=None, height=1, hide=0, greytrack=0,
                 greytrack_labels=5, greytrack_fontsize=8,
                 greytrack_font='Helvetica', greytrack_font_rotation=0,
                 greytrack_font_color = colors.Color(0.6, 0.6, 0.6),
                 scale=1, scale_format=None, scale_color=colors.black,
                 scale_font='Helvetica', scale_fontsize=6,
                 scale_fontangle=45, scale_largeticks=0.5, scale_ticks=1,
                 scale_smallticks=0.3, scale_largetick_interval=1e6,
                 scale_smalltick_interval=1e4, scale_largetick_labels=1,
                 scale_smalltick_labels=0, axis_labels=1,
                 start=None, end=None,
                 greytrack_font_colour = None, scale_colour=None):
        """ __init__(self, name=None, ...)

            o height    Int describing the relative height to other tracks in the
                        diagram

            o name      String describing the track

            o hide      Boolean, 0 if the track is not to be drawn

            o greytrack     Boolean, 1 if a grey background to the track is to be
                            drawn

            o greytrack_labels  Int describing how many track-identifying labels
                                should be placed on the track at regular intervals

            o greytrack_font    String describing the font to use for the greytrack
                                labels

            o greytrack_fontsize    Int describing the font size to display the
                                    labels on the grey track

            o greytrack_font_rotation   Int describing the angle through which to
                                        rotate the grey track labels

            o greytrack_font_color     colors.Color describing the color to draw
                                       the grey track labels (overridden by
                                       backwards compatible argument with UK
                                       spelling, colour).

            o scale     Boolean, 1 if a scale is to be drawn on the track

            o scale_color  colors.Color to draw the elemnts of the scale
                           (overridden by backwards compatible argument with UK
                           spelling, colour).

            o scale_font    String describing the font to use for the scale labels

            o scale_fontsize    Int describing the size of the scale label font

            o scale_fontangle   Int describing the angle at which to draw the scale
                                labels (linear only)

            o scale_ticks       Boolean, 1 if ticks should be drawn at all on the
                                scale

            o scale_largeticks  Float (0->1) describing the height of large
                                scale ticks relative to the track height.

            o scale_smallticks  Float (0->1) describing the height of large
                                scale ticks relative to the track height.

            o scale_largetick_interval  Int, describing the number of bases that
                                        should separate large ticks

            o scale_smalltick_interval  Int, describing the number of bases that
                                        should separate small ticks

            o scale_largetick_labels    Boolean describing whether position labels
                                        should be written over large ticks
                                        
            o scale_smalltick_labels    Boolean describing whether position labels
                                        should be written over small ticks

            o name          String to help identify the track
            
            o height        Relative height to draw the track

            o axis_labels       Boolean describing whether the value labels should
                                be placed on the Y axes
        """
        #Let the UK spelling (colour) override the USA spelling (color)
        if greytrack_font_colour is not None:
            greytrack_font_color = greytrack_font_colour
        if scale_colour is not None:
            scale_color = scale_colour

        self._next_id = 0       # This will count sets as they are added to the track
        self._sets = {}         # Holds sets, keyed by unique ID

        # Assign attribute values from instantiation
        self.height = height
        if name is not None:
            self.name = str(name)
        else:
            self.name = "Track"
        self.hide = hide
        self.start = start
        self.end = end

        # Attributes for the grey track background and labels
        self.greytrack = greytrack
        self.greytrack_labels = greytrack_labels
        self.greytrack_fontsize = greytrack_fontsize
        self.greytrack_font = greytrack_font
        self.greytrack_font_rotation = greytrack_font_rotation
        self.greytrack_fontcolor = greytrack_font_color

        # Attributes for the track scale
        self.scale = scale
        self.scale_format = scale_format
        self.scale_color = scale_color
        self.scale_font = scale_font
        self.scale_fontsize = scale_fontsize
        self.scale_fontangle = scale_fontangle
        self.scale_ticks = scale_ticks
        self.scale_largeticks = scale_largeticks
        self.scale_smallticks = scale_smallticks
        self.scale_largetick_interval = scale_largetick_interval
        self.scale_smalltick_interval = scale_smalltick_interval
        self.scale_largetick_labels = scale_largetick_labels
        self.scale_smalltick_labels = scale_smalltick_labels
        self.axis_labels = axis_labels
        

    def add_set(self, set):
        """ add_set(self, set)

            o set       A FeatureSet or GraphSet object

            Add a preexisting FeatureSet or GraphSet object to the track
        """
        set.id = self._next_id          # Assign unique id to set
        set.parent = self               # Make set's parent this track
        self._sets[self._next_id] = set # Add set, keyed by unique id
        self._next_id += 1              # Increment unique set ids


    def new_set(self, type='feature', **args):
        """ new_set(self, type='feature') -> FeatureSet or GraphSet

            Create a new FeatureSet or GraphSet object, add it to the
            track, and return for user manipulation
        """
        type_dict = {'feature': FeatureSet,
                     'graph': GraphSet
                     }        
        set = type_dict[type]()
        for key in args:
            setattr(set, key, args[key])        
        set.id = self._next_id          # Assign unique id to set
        set.parent = self               # Make set's parent this track
        self._sets[self._next_id] = set # Add set, keyed by unique id
        self._next_id += 1              # Increment unique set ids
        return set


    def del_set(self, set_id):
        """ del_set(self, set_id)

            o set_id        The unique id for the set in this track

            Remove the set with the passed id from the track
        """
        del self._sets[set_id]


    def get_sets(self):
        """ get_sets(self) -> FeatureSet or GraphSet

            Return the sets contained in this track
        """
        return self._sets.values()


    def get_ids(self):
        """ get_ids(self) -> [int, int, ...]

            Return the ids of all sets contained in this track
        """
        return self._sets.keys()


    def range(self):
        """ range(self) -> (int, int)

            Returns the lowest and highest base (or mark) numbers as a tuple
        """
        lows, highs = [], []            # Holds set of low and high values from sets
        if self.start is not None:
            lows.append(self.start)
        if self.end is not None:
            highs.append(self.end)
        for set in self._sets.values():
            low, high = set.range()     # Get each set range
            lows.append(low)
            highs.append(high)
        if lows:
            low = min(lows)
        else:
            low = None
        if highs:
            high = max(highs)
        else:
            high = None
        return low, high  # Return lowest and highest values
    

    def to_string(self, verbose=0):
        """ to_string(self, verbose=0) -> ""

            o verbose       Boolean indicating whether a short or complete
                            account of the track is required

            Returns a formatted string with information about the track
        """
        if not verbose:             # Return the short description
            return "%s" % self      # Use __str__ method instead
        else:                       # Return the long desription
            outstr = ["\n<%s: %s>" % (self.__class__, self.name)]
            outstr.append("%d sets" % len(self._sets))
            for key in self._sets:
                outstr.append("set: %s" % self._sets[key])
            return "\n".join(outstr)      


    def __getitem__(self, key):
        """ __getitem__(self, key) -> int

            o key       The id of a set in the track

            Return the set with the passed id
        """
        return self._sets[key]


    def __str__(self):
        """ __str__(self) -> ""

            Returns a formatted string with information about the Track
        """
        outstr = ["\n<%s: %s>" % (self.__class__, self.name)]
        outstr.append("%d sets" % len(self._sets))
        return "\n".join(outstr)
    


################################################################################
# RUN AS SCRIPT
################################################################################

if __name__ == '__main__':

    # test code
    from Bio import SeqIO
    from Bio.SeqFeature import SeqFeature
    from _FeatureSet import FeatureSet
    from _GraphSet import GraphSet
    from random import normalvariate
    
    genbank_entry = SeqIO.read('/data/genomes/Bacteria/Nanoarchaeum_equitans/NC_005213.gbk', 'gb')

    gdfs1 = FeatureSet(0, 'Nanoarchaeum equitans CDS - CDS')
    gdfs2 = FeatureSet(1, 'Nanoarchaeum equitans CDS - gene')
    for feature in genbank_entry.features:
        if feature.type == 'CDS':
            gdfs1.add_feature(feature)
        if feature.type == 'gene':
            gdfs2.add_feature(feature)

    gdt = Track()
    gdt.add_set(gdfs1)
    gdt.add_set(gdfs2)

    graphdata = []
    for pos in xrange(1, len(genbank_entry.seq), 1000):
        graphdata.append((pos, normalvariate(0.5, 0.1)))
    gdgs = GraphSet(2, 'test data')
    gdgs.add_graph(graphdata, 'Test Data')
    gdt.add_set(gdgs)
    
    print gdt.get_ids()
    sets = gdt.get_sets()
    for set in sets:
        print set

    print gdt.get_element_limits()
