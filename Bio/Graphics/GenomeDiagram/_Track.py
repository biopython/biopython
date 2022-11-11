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

"""Track module.

Provides:
 - Track - Container for a single track on the diagram, containing
   FeatureSet and GraphSet objects

For drawing capabilities, this module uses reportlab to draw and write
the diagram: http://www.reportlab.com
"""


from reportlab.lib import colors

# GenomeDiagram imports
from ._FeatureSet import FeatureSet
from ._GraphSet import GraphSet

_grey = colors.Color(0.6, 0.6, 0.6)


class Track:
    """Track.

    Attributes:
     - height    Int describing the relative height to other trackscale_fontsizes
       in the diagram
     - name      String describing the track
     - hide      Boolean, 0 if the track is not to be drawn
     - start, end    Integers (or None) specifying start/end to draw just
       a partial track.
     - greytrack     Boolean, 1 if a grey background to the track is to be
       drawn
     - greytrack_labels  Int describing how many track-identifying labels
       should be placed on the track at regular intervals
     - greytrack_font    String describing the font to use for the greytrack
       labels
     - greytrack_fontsize    Int describing the font size to display the
       labels on the grey track
     - greytrack_font_rotation   Int describing the angle through which to
       rotate the grey track labels (Linear only)
     - greytrack_font_color     colors.Color describing the color to draw
       the grey track labels
     - scale     Boolean, 1 if a scale is to be drawn on the track
     - scale_format  String, defaults to None, when scale values are written
       as numerals.  Setting this to 'SInt' invokes SI
       unit-like multiples, such as Mbp, Kbp and so on.
     - scale_color  colors.Color to draw the elements of the scale
     - scale_font    String describing the font to use for the scale labels
     - scale_fontsize    Int describing the size of the scale label font
     - scale_fontangle   Int describing the angle at which to draw the scale
       labels (linear only)
     - scale_ticks       Boolean, 1 if ticks should be drawn at all on the
       scale
     - scale_largeticks  Float (0->1) describing the height of large
       scale ticks relative to the track height.
     - scale_smallticks  Float (0->1) describing the height of large
       scale ticks relative to the track height.
     - scale_largetick_interval  Int, describing the number of bases that
       should separate large ticks
     - scale_smalltick_interval  Int, describing the number of bases that
       should separate small ticks
     - scale_largetick_labels    Boolean describing whether position labels
       should be written over large ticks
     - scale_smalltick_labels    Boolean describing whether position labels
       should be written over small ticks
     - axis_labels       Boolean describing whether the value labels should
       be placed on the Y axes

    """

    def __init__(
        self,
        name=None,
        height=1,
        hide=0,
        greytrack=0,
        greytrack_labels=5,
        greytrack_fontsize=8,
        greytrack_font="Helvetica",
        greytrack_font_rotation=0,
        greytrack_font_color=_grey,
        scale=1,
        scale_format=None,
        scale_color=colors.black,
        scale_font="Helvetica",
        scale_fontsize=6,
        scale_fontangle=45,
        scale_largeticks=0.5,
        scale_ticks=1,
        scale_smallticks=0.3,
        scale_largetick_interval=1e6,
        scale_smalltick_interval=1e4,
        scale_largetick_labels=1,
        scale_smalltick_labels=0,
        axis_labels=1,
        start=None,
        end=None,
        greytrack_font_colour=None,
        scale_colour=None,
    ):
        """Initialize.

        Arguments:
         - height    Int describing the relative height to other tracks in the
           diagram
         - name      String describing the track
         - hide      Boolean, 0 if the track is not to be drawn
         - greytrack     Boolean, 1 if a grey background to the track is to be
           drawn
         - greytrack_labels  Int describing how many track-identifying labels
           should be placed on the track at regular intervals
         - greytrack_font    String describing the font to use for the greytrack
           labels
         - greytrack_fontsize    Int describing the font size to display the
           labels on the grey track
         - greytrack_font_rotation   Int describing the angle through which to
           rotate the grey track labels (Linear only)
         - greytrack_font_color     colors.Color describing the color to draw
           the grey track labels (overridden by backwards compatible argument
           with UK spelling, colour).
         - scale     Boolean, 1 if a scale is to be drawn on the track
         - scale_color  colors.Color to draw the elements of the scale
           (overridden by backwards compatible argument with UK
           spelling, colour).
         - scale_font    String describing the font to use for the scale labels
         - scale_fontsize    Int describing the size of the scale label font
         - scale_fontangle   Int describing the angle at which to draw the scale
           labels (linear only)
         - scale_ticks       Boolean, 1 if ticks should be drawn at all on the
           scale
         - scale_largeticks  Float (0->1) describing the height of large
           scale ticks relative to the track height.
         - scale_smallticks  Float (0->1) describing the height of large
           scale ticks relative to the track height.
         - scale_largetick_interval  Int, describing the number of bases that
           should separate large ticks
         - scale_smalltick_interval  Int, describing the number of bases that
           should separate small ticks
         - scale_largetick_labels    Boolean describing whether position labels
           should be written over large ticks
         - scale_smalltick_labels    Boolean describing whether position labels
           should be written over small ticks
         - name          String to help identify the track
         - height        Relative height to draw the track
         - axis_labels       Boolean describing whether the value labels should
           be placed on the Y axes

        """
        # Let the UK spelling (colour) override the USA spelling (color)
        if greytrack_font_colour is not None:
            greytrack_font_color = greytrack_font_colour
        if scale_colour is not None:
            scale_color = scale_colour

        self._next_id = 0  # This will count sets as they are added to the track
        self._sets = {}  # Holds sets, keyed by unique ID

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
        """Add a preexisting FeatureSet or GraphSet object to the track."""
        set.id = self._next_id  # Assign unique id to set
        set.parent = self  # Make set's parent this track
        self._sets[self._next_id] = set  # Add set, keyed by unique id
        self._next_id += 1  # Increment unique set ids

    def new_set(self, type="feature", **args):
        """Create a new FeatureSet or GraphSet object.

        Create a new FeatureSet or GraphSet object, add it to the
        track, and return for user manipulation
        """
        type_dict = {"feature": FeatureSet, "graph": GraphSet}
        set = type_dict[type]()
        for key in args:
            setattr(set, key, args[key])
        set.id = self._next_id  # Assign unique id to set
        set.parent = self  # Make set's parent this track
        self._sets[self._next_id] = set  # Add set, keyed by unique id
        self._next_id += 1  # Increment unique set ids
        return set

    def del_set(self, set_id):
        """Remove the set with the passed id from the track."""
        del self._sets[set_id]

    def get_sets(self):
        """Return the sets contained in this track."""
        return list(self._sets.values())

    def get_ids(self):
        """Return the ids of all sets contained in this track."""
        return list(self._sets.keys())

    def range(self):
        """Return the lowest and highest base (or mark) numbers as a tuple."""
        lows, highs = [], []  # Holds set of low and high values from sets
        if self.start is not None:
            lows.append(self.start)
        if self.end is not None:
            highs.append(self.end)
        for set in self._sets.values():
            low, high = set.range()  # Get each set range
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
        """Return a formatted string with information about the track.

        Arguments:
         - verbose - Boolean indicating whether a short or complete
           account of the track is required

        """
        if not verbose:  # Return the short description
            return f"{self}"  # Use __str__ method instead
        else:  # Return the long description
            outstr = [f"\n<{self.__class__}: {self.name}>"]
            outstr.append("%d sets" % len(self._sets))
            for key in self._sets:
                outstr.append(f"set: {self._sets[key]}")
            return "\n".join(outstr)

    def __getitem__(self, key):
        """Return the set with the passed id."""
        return self._sets[key]

    def __str__(self):
        """Return a formatted string with information about the Track."""
        outstr = [f"\n<{self.__class__}: {self.name}>"]
        outstr.append("%d sets" % len(self._sets))
        return "\n".join(outstr)
