# Copyright 2003-2008 by Leighton Pritchard.  All rights reserved.
# Revisions copyright 2008-2009 by Peter Cock.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
# Contact:       Leighton Pritchard, Scottish Crop Research Institute,
#                Invergowrie, Dundee, Scotland, DD2 5DA, UK
#                L.Pritchard@scri.ac.uk
################################################################################

""" LinearDrawer module

    Provides:

    o LinearDrawer -  Drawing object for linear diagrams

    For drawing capabilities, this module uses reportlab to draw and write
    the diagram:

    http://www.reportlab.com

    For dealing with biological information, the package expects BioPython
    objects:

    http://www.biopython.org
"""

# ReportLab imports
from reportlab.graphics.shapes import *
from reportlab.lib import colors

# GenomeDiagram imports
from _AbstractDrawer import AbstractDrawer, draw_box, draw_arrow
from _AbstractDrawer import intermediate_points, angle2trig
from _FeatureSet import FeatureSet
from _GraphSet import GraphSet

from math import ceil

class LinearDrawer(AbstractDrawer):
    """ LinearDrawer(AbstractDrawer)

        Inherits from:

        o AbstractDrawer

        Provides:

        Methods:

        o __init__(self, parent=None, pagesize='A3', orientation='landscape',
                 x=0.05, y=0.05, xl=None, xr=None, yt=None, yb=None,
                 start=None, end=None, tracklines=0, fragments=10,
                 fragment_size=0.9, track_size=0.75) Called on instantiation

        o set_page_size(self, pagesize, orientation)    Set the page size to the
                                                    passed size and orientation

        o set_margins(self, x, y, xl, xr, yt, yb)   Set the drawable area of the
                                                    page

        o set_bounds(self, start, end)  Set the bounds for the elements to be
                                        drawn

        o is_in_bounds(self, value)     Returns a boolean for whether the position
                                        is actually to be drawn

        o __len__(self)     Returns the length of sequence that will be drawn


        o draw(self)    Place the drawing elements on the diagram

        o init_fragments(self)  Calculate information
                                about sequence fragment locations on the drawing

        o set_track_heights(self)   Calculate information about the offset of
                                    each track from the fragment base
                                    
        o draw_test_tracks(self)    Add lines demarcating each track to the
                                    drawing

        o draw_track(self, track)   Return the contents of the passed track as
                                    drawing elements

        o draw_scale(self, track)   Return a scale for the passed track as
                                    drawing elements

        o draw_tick(self, tickpos, ctr, ticklen, track, draw_label) Return a
                                    tick line and possibly a label

        o draw_greytrack(self, track)   Return a grey background and superposed
                                        label for the passed track as drawing
                                        elements

        o draw_feature_set(self, set)   Return the features in the passed set as
                                        drawing elements

        o draw_feature(self, feature)   Return a single feature as drawing
                                        elements

        o get_feature_sigil(self, feature, x0, x1, fragment)    Return a single
                                        feature as its sigil in drawing elements

        o draw_graph_set(self, set)     Return the data in a set of graphs as
                                        drawing elements

        o draw_line_graph(self, graph)  Return the data in a graph as a line
                                        graph in drawing elements

        o draw_heat_graph(self, graph)  Return the data in a graph as a heat
                                        graph in drawing elements

        o draw_bar_graph(self, graph)   Return the data in a graph as a bar
                                        graph in drawing elements

        o canvas_location(self, base)   Return the fragment, and the offset from
                                        the left margin, of a passed position
                                        in the sequence, on the diagram.

        Attributes:

        o tracklines    Boolean for whether to draw lines dilineating tracks

        o pagesize      Tuple describing the size of the page in pixels

        o x0            Float X co-ord for leftmost point of drawable area

        o xlim          Float X co-ord for rightmost point of drawable area

        o y0            Float Y co-ord for lowest point of drawable area

        o ylim          Float Y co-ord for topmost point of drawable area

        o pagewidth     Float pixel width of drawable area

        o pageheight    Float pixel height of drawable area

        o xcenter       Float X co-ord of center of drawable area

        o ycenter       Float Y co-ord of center of drawable area

        o start         Int, base to start drawing from

        o end           Int, base to stop drawing at

        o length        Int, size of sequence to be drawn

        o fragments     Int, number of fragments into which to divide the
                        drawn sequence

        o fragment_size Float (0->1) the proportion of the fragment height to
                        draw in

        o track_size    Float (0->1) the proportion of the track height to
                        draw in

        o drawing       Drawing canvas

        o drawn_tracks  List of ints denoting which tracks are to be drawn

        o current_track_level   Int denoting which track is currently being
                                drawn

        o fragment_height   Float total fragment height in pixels

        o fragment_bases    Int total fragment length in bases

        o fragment_lines    Dictionary of top and bottom y-coords of fragment,
                            keyed by fragment number

        o fragment_limits   Dictionary of start and end bases of each fragment,
                            keyed by fragment number

        o track_offsets     Dictionary of number of pixels that each track top,
                            center and bottom is offset from the base of a
                            fragment, keyed by track

    """
    def __init__(self, parent=None, pagesize='A3', orientation='landscape',
                 x=0.05, y=0.05, xl=None, xr=None, yt=None, yb=None,
                 start=None, end=None, tracklines=0, fragments=10,
                 fragment_size=0.9, track_size=0.75):
        """ __init__(self, parent, pagesize='A3', orientation='landscape',
                     x=0.05, y=0.05, xl=None, xr=None, yt=None, yb=None,
                     start=None, end=None, tracklines=0, fragments=10,
                     fragment_size=0.9, track_size=0.75)

            o parent    Diagram object containing the data that the drawer
                        draws

            o pagesize  String describing the ISO size of the image, or a tuple
                        of pixels

            o orientation   String describing the required orientation of the
                            final drawing ('landscape' or 'portrait')

            o x         Float (0->1) describing the relative size of the X
                        margins to the page

            o y         Float (0->1) describing the relative size of the Y
                        margins to the page

            o xl        Float (0->1) describing the relative size of the left X
                        margin to the page (overrides x)

            o xl        Float (0->1) describing the relative size of the left X
                        margin to the page (overrides x)

            o xr        Float (0->1) describing the relative size of the right X
                        margin to the page (overrides x)

            o yt        Float (0->1) describing the relative size of the top Y
                        margin to the page (overrides y)

            o yb        Float (0->1) describing the relative size of the lower Y
                        margin to the page (overrides y)

            o start     Int, the position to begin drawing the diagram at

            o end       Int, the position to stop drawing the diagram at

            o tracklines    Boolean flag to show (or not) lines delineating tracks
                            on the diagram            

            o fragments Int, the number of equal fragments into which the
                        sequence should be divided for drawing

            o fragment_size Float(0->1) The proportion of the available height
                            for the fragment that should be taken up in drawing
                            
            o track_size    The proportion of the available track height that
                            should be taken up in drawing

        """
        # Use the superclass' instantiation method
        AbstractDrawer.__init__(self, parent, pagesize, orientation,
                                  x, y, xl, xr, yt, yb, start, end,
                                  tracklines)

        # Useful measurements on the page
        self.fragments = fragments
        self.fragment_size = fragment_size
        self.track_size = track_size


    def draw(self):
        """ draw(self)

            Draw a linear diagram of the data in the parent Diagram object
        """
        # Instantiate the drawing canvas
        self.drawing = Drawing(self.pagesize[0], self.pagesize[1])

        feature_elements = []           # holds feature elements
        feature_labels = []             # holds feature labels
        greytrack_bgs = []              # holds track background
        greytrack_labels = []           # holds track foreground labels
        scale_axes = []                 # holds scale axes
        scale_labels = []               # holds scale axis labels

        # Get the tracks to be drawn
        self.drawn_tracks = self._parent.get_drawn_levels()

        # Set fragment and track sizes
        self.init_fragments()
        self.set_track_heights()

        # Go through each track in the parent (if it is to be drawn) one by
        # one and collate the data as drawing elements
        for track_level in self.drawn_tracks: # only use tracks to be drawn
            self.current_track_level = track_level      # establish track level
            track = self._parent[track_level]           # get the track at that level
            gbgs, glabels = self.draw_greytrack(track)  # get greytrack elements
            greytrack_bgs.append(gbgs)
            greytrack_labels.append(glabels)
            features, flabels = self.draw_track(track)  # get feature and graph elements
            feature_elements.append(features)
            feature_labels.append(flabels)
            if track.scale:
                axes, slabels = self.draw_scale(track)      # get scale elements
                scale_axes.append(axes)
                scale_labels.append(slabels)

        # Groups listed in order of addition to page (from back to front)
        # Draw track backgrounds
        # Draw features and graphs
        # Draw scale axes
        # Draw scale labels
        # Draw feature labels
        # Draw track labels
        element_groups = [greytrack_bgs, feature_elements, scale_axes,
                         scale_labels, feature_labels, greytrack_labels]
        for element_group in element_groups:
            for element_list in element_group:
                [self.drawing.add(element) for element in element_list]
            
        if self.tracklines:             # Draw test tracks over top of diagram
            self.draw_test_tracks()


    def init_fragments(self):
        """ init_fragments(self)

            Initialises useful values for calculating the positioning of
            diagram elements
        """
        # Set basic heights, lengths etc
        self.fragment_height = 1.*self.pageheight/self.fragments     # total fragment height in pixels
        self.fragment_bases = ceil(1.*self.length/self.fragments)    # fragment length in bases

        # Key fragment base and top lines by fragment number
        self.fragment_lines = {}    # Holds bottom and top line locations of fragments, keyed by fragment number
        fragment_crop = (1-self.fragment_size)/2    # No of pixels to crop the fragment
        fragy = self.ylim           # Holder for current absolute fragment base
        for fragment in range(self.fragments):
            fragtop = fragy-fragment_crop * self.fragment_height     # top - crop
            fragbtm = fragy-(1-fragment_crop) * self.fragment_height # bottom + crop
            self.fragment_lines[fragment] = (fragbtm, fragtop)
            fragy -= self.fragment_height                   # next fragment base

        # Key base starts and ends for each fragment by fragment number
        self.fragment_limits = {}   # Holds first and last base positions in a fragment
        fragment_step = self.fragment_bases  # bases per fragment
        fragment_count = 0
        # Add start and end positions for each fragment to dictionary
        for marker in range(int(self.start), int(self.end), int(fragment_step)):
            self.fragment_limits[fragment_count] = (marker, marker+fragment_step)
            fragment_count += 1


    def set_track_heights(self):
        """ set_track_heights(self)

            Since tracks may not be of identical heights, the bottom and top
            offsets of each track relative to the fragment top and bottom is
            stored in a dictionary - self.track_offsets, keyed by track number
        """
        top_track = max(self.drawn_tracks)     # The 'highest' track number to draw

        trackunit_sum = 0           # Total number of 'units' for the tracks
        trackunits = {}             # The start and end units for each track, keyed by track number
        heightholder = 0            # placeholder variable
        for track in range(1, top_track+1):    # for all track numbers to 'draw'
            try:
                trackheight = self._parent[track].height    # Get track height
            except:
                trackheight = 1                             # ...or default to 1
            trackunit_sum += trackheight    # increment total track unit height
            trackunits[track] = (heightholder, heightholder+trackheight)
            heightholder += trackheight     # move to next height
        trackunit_height = 1.*self.fragment_height*self.fragment_size/trackunit_sum

        # Calculate top and bottom offsets for each track, relative to fragment
        # base
        track_offsets = {}      # The offsets from fragment base for each track
        track_crop = trackunit_height*(1-self.track_size)/2.    # 'step back' in pixels
        assert track_crop >= 0
        for track in trackunits:
            top = trackunits[track][1]*trackunit_height-track_crop  # top offset
            btm = trackunits[track][0]*trackunit_height+track_crop  # bottom offset
            ctr = btm+(top-btm)/2.                          # center offset
            track_offsets[track] = (btm, ctr, top)          
        self.track_offsets = track_offsets
        

    def draw_test_tracks(self):
        """ draw_test_tracks(self)

            Draw red lines indicating the top and bottom of each fragment,
            and blue ones indicating tracks to be drawn.
        """
        # Add lines for each fragment
        for fbtm, ftop in self.fragment_lines.values():
            self.drawing.add(Line(self.x0, ftop, self.xlim, ftop,
                                  strokeColor=colors.red))  # top line
            self.drawing.add(Line(self.x0, fbtm, self.xlim, fbtm,
                                  strokeColor=colors.red))  # bottom line

            # Add track lines for this fragment - but only for drawn tracks
            for track in self.drawn_tracks:
                trackbtm = fbtm + self.track_offsets[track][0]
                trackctr = fbtm + self.track_offsets[track][1]
                tracktop = fbtm + self.track_offsets[track][2]
                self.drawing.add(Line(self.x0, tracktop, self.xlim, tracktop,
                                      strokeColor=colors.blue))  # top line
                self.drawing.add(Line(self.x0, trackctr, self.xlim, trackctr,
                                      strokeColor=colors.green))  # center line
                self.drawing.add(Line(self.x0, trackbtm, self.xlim, trackbtm,
                                      strokeColor=colors.blue))  # bottom line
        

    def draw_track(self, track):
        """ draw_track(self, track) -> ([element, element,...], [element, element,...])

            o track     Track object

            Returns a tuple (list of elements in the track, list of labels in
            the track)
        """
        track_elements = []     # Holds elements from features and graphs
        track_labels = []       # Holds labels from features and graphs
        
        # Distribution dictionary for dealing with different set types
        set_methods = {FeatureSet: self.draw_feature_set,
                       GraphSet: self.draw_graph_set
                       }
        
        for set in track.get_sets():        # Draw the feature or graph sets
            elements, labels = set_methods[set.__class__](set)
            track_elements += elements
            track_labels += labels
        return track_elements, track_labels
        

    def draw_tick(self, tickpos, ctr, ticklen, track, draw_label):
        """ draw_tick(self, tickpos, ctr, ticklen) -> (element, element)

            o tickpos   Int, position of the tick on the sequence

            o ctr       Float, Y co-ord of the center of the track

            o ticklen   How long to draw the tick

            o track     Track, the track the tick is drawn on

            o draw_label    Boolean, write the tick label?

            Returns a drawing element that is the tick on the scale
        """
        assert self.start <= tickpos and tickpos <= self.end, \
               "Tick at %i, but showing %i to %i" \
               % (tickpos, self.start, self.end)
        fragment, tickx = self.canvas_location(tickpos) # Tick co-ordinates
        assert fragment >=0, \
               "Fragment %i, tickpos %i" % (fragment, tickpos)
        tctr = ctr + self.fragment_lines[fragment][0]   # Center line of the track
        tickx += self.x0                # Tick X co-ord
        ticktop = tctr + ticklen        # Y co-ord of tick top
        tick = Line(tickx, tctr, tickx, ticktop, strokeColor=track.scale_color)
        if draw_label: # Put tick position on as label
            if track.scale_format == 'SInt':
                if tickpos >= 1000000:
                    tickstring = str(tickpos//1000000) + " Mbp"
                elif tickpos >= 1000:
                    tickstring = str(tickpos//1000) + " Kbp"
                else:
                    tickstring = str(tickpos)
            else:
                tickstring = str(tickpos)
            label = String(0, 0, tickstring,  # Make label string
                   fontName=track.scale_font,
                   fontSize=track.scale_fontsize,
                   fillColor=track.scale_color)
            labelgroup = Group(label)
            rotation = angle2trig(track.scale_fontangle)
            labelgroup.transform = (rotation[0], rotation[1], rotation[2],
                                    rotation[3], tickx, ticktop)
        else:
            labelgroup = None
        return tick, labelgroup

    def draw_scale(self, track):
        """ draw_scale(self, track) -> ([element, element,...], [element, element,...])

            o track     Track object

            Returns a tuple of (list of elements in the scale, list of labels
            in the scale)
        """
        scale_elements = []     # Holds axes and ticks
        scale_labels = []       # Holds labels

        if not track.scale:     # No scale required, exit early
            return [], []

        # Get track location
        btm, ctr, top = self.track_offsets[self.current_track_level]
        trackheight = (top-ctr)

        # For each fragment, draw the scale for this track
        for fragment in range(self.fragments):
            tbtm = btm + self.fragment_lines[fragment][0]
            tctr = ctr + self.fragment_lines[fragment][0]
            ttop = top + self.fragment_lines[fragment][0]
            # X-axis
            if fragment == self.fragments - 1:
                frag, x = self.canvas_location(self.end)
                scale_elements.append(Line(self.x0, tctr, self.x0 + x, tctr,
                                       strokeColor=track.scale_color))
            else:
                scale_elements.append(Line(self.x0, tctr, self.xlim, tctr,
                                           strokeColor=track.scale_color))
            # Y-axis
            scale_elements.append(Line(self.x0, tbtm, self.x0, ttop,
                                       strokeColor=track.scale_color))
        if track.scale_ticks:   # Ticks are required on the scale
            # Draw large ticks
            #I want the ticks to be consistently positioned relative to
            #the start of the sequence (position 0), not relative to the
            #current viewpoint (self.start and self.end)

            ticklen = track.scale_largeticks * trackheight
            tickiterval = int(track.scale_largetick_interval)
            #Note that we could just start the list of ticks using
            #range(0,self.end,tickinterval) and the filter out the
            #ones before self.start - but this seems wasteful.
            #Using tickiterval * (self.start//tickiterval) is a shortcut.
            largeticks = [pos for pos \
                          in range(tickiterval * (self.start//tickiterval),
                                   int(self.end),
                                   tickiterval) \
                          if pos >= self.start]
            for tickpos in largeticks:
                tick, label = self.draw_tick(tickpos, ctr, ticklen,
                                             track,
                                             track.scale_largetick_labels)
                scale_elements.append(tick)
                if label is not None:   # If there's a label, add it
                    scale_labels.append(label)
            # Draw small ticks
            ticklen = track.scale_smallticks * trackheight
            tickiterval = int(track.scale_smalltick_interval)
            smallticks = [pos for pos \
                          in range(tickiterval * (self.start//tickiterval),
                                   int(self.end),
                                   tickiterval) \
                          if pos >= self.start]
            for tickpos in smallticks:
                tick, label = self.draw_tick(tickpos, ctr, ticklen,
                                             track,
                                             track.scale_smalltick_labels)
                scale_elements.append(tick)
                if label is not None:   # If there's a label, add it
                    scale_labels.append(label)

        # Check to see if the track contains a graph - if it does, get the
        # minimum and maximum values, and put them on the scale Y-axis
        if track.axis_labels:
            for set in track.get_sets():            # Check all sets...
                if set.__class__ is GraphSet:     # ...for a graph set
                    graph_label_min = []
                    graph_label_mid = []                    
                    graph_label_max = []
                    for graph in set.get_graphs():
                        quartiles = graph.quartiles()
                        minval, maxval = quartiles[0], quartiles[4]
                        if graph.center is None:
                            midval = (maxval + minval)/2.
                            graph_label_min.append("%.3f" % minval)
                            graph_label_max.append("%.3f" % maxval)
                        else:
                            diff = max((graph.center-minval),
                                       (maxval-graph.center))
                            minval = graph.center-diff
                            maxval = graph.center+diff
                            midval = graph.center
                            graph_label_mid.append("%.3f" % midval)
                            graph_label_min.append("%.3f" % minval)
                            graph_label_max.append("%.3f" % maxval)
                    for fragment in range(self.fragments):  # Add to all fragment axes
                        tbtm = btm + self.fragment_lines[fragment][0]
                        tctr = ctr + self.fragment_lines[fragment][0]
                        ttop = top + self.fragment_lines[fragment][0]
                        for val, pos in [(";".join(graph_label_min), tbtm),
                                         (";".join(graph_label_max), ttop),
                                         (";".join(graph_label_mid), tctr)]:
                            label = String(0, 0, val,
                                           fontName=track.scale_font,
                                           fontSize=track.scale_fontsize,
                                           fillColor=track.scale_color)
                            labelgroup = Group(label)
                            rotation = angle2trig(track.scale_fontangle)
                            labelgroup.transform = (rotation[0], rotation[1], rotation[2],
                                                    rotation[3], self.x0, pos)
                            scale_labels.append(labelgroup)

        return scale_elements, scale_labels
        
        

    def draw_greytrack(self, track):
        """ draw_greytrack(self) -> ([element, element,...], [element, element,...])

            o track     Track object

            Put in a grey background to the current track in all fragments,
            if track specifies that we should
        """        
        greytrack_bgs = []      # Holds grey track backgrounds
        greytrack_labels = []   # Holds grey foreground labels

        if not track.greytrack: # No greytrack required, return early
            return [], []

        # Get track location
        btm, ctr, top = self.track_offsets[self.current_track_level]

        # Add greytrack to all fragments for this track
        for fragment in range(self.fragments):
            tbtm = btm + self.fragment_lines[fragment][0]
            tctr = ctr + self.fragment_lines[fragment][0]
            ttop = top + self.fragment_lines[fragment][0]  
            box = draw_box((self.x0, tbtm), (self.xlim, ttop),  # Grey track bg
                           colors.Color(0.96,0.96, 0.96))       # is just a box
            greytrack_bgs.append(box)

            if track.greytrack_labels:  # If labels are required
                labelstep = (self.pagewidth)/track.greytrack_labels # how far apart should they be?
                label = String(0, 0, track.name,    # label contents
                               fontName=track.greytrack_font,
                               fontSize=track.greytrack_fontsize,
                               fillColor=track.greytrack_fontcolor)
                # Create a new labelgroup at each position the label is required
                for x in range(int(self.x0), int(self.xlim), int(labelstep)):
                    labelgroup = Group(label)
                    rotation = angle2trig(track.greytrack_font_rotation)
                    labelgroup.transform = (rotation[0], rotation[1], rotation[2],
                                            rotation[3], x, tbtm)
                    if not self.xlim-x <= labelstep:    # Don't overlap the end of the track
                        greytrack_labels.append(labelgroup)

        return greytrack_bgs, greytrack_labels


    def draw_feature_set(self, set):
        """ draw_feature_set(self, set) -> ([element, element,...], [element, element,...])

            o set       FeatureSet object

            Returns a tuple (list of elements describing features, list of
            labels for elements)
        """
        #print 'draw feature set'
        feature_elements = []   # Holds diagram elements belonging to the features
        label_elements = []     # Holds diagram elements belonging to feature labels 

        # Collect all the elements for the feature set
        for feature in set.get_features():
            if self.is_in_bounds(feature.start) or self.is_in_bounds(feature.end):
                features, labels = self.draw_feature(feature)   # get elements and labels
                feature_elements += features
                label_elements += labels

        return feature_elements, label_elements


    def draw_feature(self, feature):
        """ draw_feature(self, feature, parent_feature=None) -> ([element, element,...], [element, element,...])

            o feature           Feature containing location info

            Returns tuple of (list of elements describing single feature, list
            of labels for those elements)
        """
        if feature.hide:        # Feature hidden, don't draw it...
            return [], []
        
        feature_elements = []   # Holds diagram elements belonging to the feature
        label_elements = []     # Holds labels belonging to the feature

        # A single feature may be split into subfeatures, so loop over them
        for start, end in feature.locations:
            #print start, end, feature.name
            # Forward strand, start > end as it overlaps zero
            if start > end:
                locs = [(start, self.end), (self.start, end)]
            else:
                locs = [(start, end)]
            #print locs
            for locstart, locend in locs:
                #print locstart, locend, feature.name
                # Correct locations in case the feature overruns the drawn sequence
                if locstart < self.start:
                    locstart = self.start
                if locend > self.end:
                    locend = self.end

                feature_boxes = self.draw_feature_location(feature, locstart, locend)
                for box, label in feature_boxes:
                    feature_elements.append(box)
                    if label is not None:
                        label_elements.append(label)
            
        return feature_elements, label_elements
            
    def draw_feature_location(self, feature, locstart, locend):
        feature_boxes = []
        # Get start and end positions for feature/subfeatures
        start_fragment, start_offset = self.canvas_location(locstart)
        end_fragment, end_offset = self.canvas_location(locend)
            #print "start_fragment, start_offset", start_fragment, start_offset
            #print "end_fragment, end_offset", end_fragment, end_offset
            #print "start, end", locstart, locend

        # Note that there is a strange situation where a feature may be in
        # several parts, and one or more of those parts may end up being
        # drawn on a non-existent fragment.  So we check that the start and
        # end fragments do actually exist in terms of the drawing
        allowed_fragments = self.fragment_limits.keys()
        if start_fragment in allowed_fragments and end_fragment in allowed_fragments:
            
            #print feature.name, feature.start, feature.end, start_offset, end_offset
            if start_fragment == end_fragment:  # Feature is found on one fragment
                feature_box, label = self.get_feature_sigil(feature, start_offset,
                                                            end_offset, start_fragment)
                feature_boxes.append((feature_box, label))
                #feature_elements.append(feature_box)
                #if label is not None:   # There is a label for the feature
                #    label_elements.append(label)                        
            else: # Feature is split over two or more fragments
                fragment = start_fragment
                start = start_offset
                # The bit that runs up to the end of the first fragment,
                # and any bits that subsequently span whole fragments
                while self.fragment_limits[fragment][1] < locend:
                    #print fragment, self.fragment_limits[fragment][1], locend
                    feature_box, label = self.get_feature_sigil(feature, start,
                                                                self.pagewidth,
                                                                fragment)
                    
                    fragment += 1   # move to next fragment
                    start = 0       # start next sigil from start of fragment
                    feature_boxes.append((feature_box, label))
                    #feature_elements.append(feature_box)
                    #if label is not None:   # There's a label for the feature
                    #    label_elements.append(label)
                # The last bit of the feature
                #print locend, self.end, fragment
                #print self.fragment_bases, self.length
                feature_box, label = self.get_feature_sigil(feature, 0,
                                                            end_offset, fragment)
                feature_boxes.append((feature_box, label))
        #if locstart > locend:
        #    print locstart, locend, feature.strand, feature_boxes, feature.name
        return feature_boxes

        

    def get_feature_sigil(self, feature, x0, x1, fragment, **kwargs):
        """ get_feature_sigil(self, feature, x0, x1, fragment) -> (element, element, element)

            o feature       Feature object

            o x0            Start X co-ordinate on diagram

            o x1            End X co-ordinate on diagram

            o fragment      The fragment on which the feature appears

            Returns a drawable indicator of the feature, and any required label
            for it
        """
        # Establish co-ordinates for drawing
        x0, x1 = self.x0 + x0, self.x0 + x1
        btm, ctr, top = self.track_offsets[self.current_track_level]
        try:
            btm += self.fragment_lines[fragment][0]
            ctr += self.fragment_lines[fragment][0]
            top += self.fragment_lines[fragment][0]
        except:     # Only called if the method screws up big time
            print "We've got a screw-up"
            print self.start, self.end
            print self.fragment_bases
            print x0, x1
            for locstart, locend in feature.locations:
                print self.canvas_location(locstart)
                print self.canvas_location(locend)
            print 'FEATURE\n', feature
            1/0

        # Distribution dictionary for various ways of drawing the feature
        # Each method takes the corners of a containing box and a color
        # as argument
        draw_methods = {'BOX': draw_box,
                        'ARROW': draw_arrow,
                        }
        method = draw_methods[feature.sigil]
        kwargs['head_length_ratio'] = feature.arrowhead_length
        kwargs['shaft_height_ratio'] = feature.arrowshaft_height

        #Support for clickable links... needs ReportLab 2.4 or later
        #which added support for links in SVG output.
        if hasattr(feature, "url") :
            kwargs["hrefURL"] = feature.url
            kwargs["hrefTitle"] = feature.name

        strand = feature.strand
        
        # Get sigil for the feature, location dependent on the feature strand
        if strand == 1:
            sigil = method((x0, ctr), (x1, top), color=feature.color,
                           orientation='right', **kwargs)
        elif strand == -1:
            sigil = method((x1, btm), (x0, ctr), color=feature.color,
                           orientation='left', **kwargs)
        else:
            sigil = method((x0, btm), (x1, top), color=feature.color,
                           **kwargs)
        if feature.label:   # Feature requires a label
            label = String(0, 0, feature.name,
                           fontName=feature.label_font,
                           fontSize=feature.label_size,
                           fillColor=feature.label_color)
            labelgroup = Group(label)
            # Feature is on top, or covers both strands (location affects
            # the height and rotation of the label)
            if feature.strand in (0, 1):    
                rotation = angle2trig(feature.label_angle)
                if feature.label_position in ('start', "5'", 'left'):
                    pos = x0
                elif feature.label_position in ('middle', 'center', 'centre'):
                    pos = (x1 + x0)/2.
                else:
                    pos = x1
                labelgroup.transform = (rotation[0], rotation[1], rotation[2],
                                        rotation[3], pos, top)
            else:   # Feature on bottom strand
                rotation = angle2trig(feature.label_angle + 180)                        
                if feature.label_position in ('start', "5'", 'left'):
                    pos = x1
                elif feature.label_position in ('middle', 'center', 'centre'):
                    pos = (x1 + x0)/2.
                else:
                    pos = x0
                labelgroup.transform = (rotation[0], rotation[1], rotation[2],
                                        rotation[3], pos, btm)
        else:
            labelgroup = None
        return sigil, labelgroup

                
    def draw_graph_set(self, set):
        """ draw_graph_set(self, set) -> ([element, element,...], [element, element,...])

            o set       GraphSet object

            Returns tuple (list of graph elements, list of graph labels)
        """
        #print 'draw graph set'
        elements = []   # Holds graph elements

        # Distribution dictionary for how to draw the graph
        style_methods = {'line': self.draw_line_graph,
                         'heat': self.draw_heat_graph,
                         'bar': self.draw_bar_graph
                         }

        for graph in set.get_graphs():
            elements += style_methods[graph.style](graph)

        return elements, []


    def draw_line_graph(self, graph):
        """ draw_line_graph(self, graph) -> [element, element,...]

            o graph     Graph object

            Returns a line graph as a list of drawable elements 
        """
        #print '\tdraw_line_graph'
        line_elements = []                  # Holds drawable elements

        # Get graph data
        data_quartiles = graph.quartiles()
        minval, maxval = data_quartiles[0],data_quartiles[4]
        btm, ctr, top = self.track_offsets[self.current_track_level]
        trackheight = 0.5*(top-btm)
        datarange = maxval - minval
        if datarange == 0:
            datarange = trackheight
        data = graph[self.start:self.end]

        # midval is the value at which the x-axis is plotted, and is the
        # central ring in the track
        if graph.center is None:
            midval = (maxval + minval)/2.    
        else:
            midval = graph.center
        # Whichever is the greatest difference: max-midval or min-midval, is
        # taken to specify the number of pixel units resolved along the
        # y-axis
        resolution = max((midval-minval), (maxval-midval))

        # Start from first data point
        pos, val = data[0]
        lastfrag, lastx = self.canvas_location(pos)
        lastx += self.x0        # Start xy co-ords
        lasty = trackheight*(val-midval)/resolution + \
                self.fragment_lines[lastfrag][0] + ctr
        lastval = val
        # Add a series of lines linking consecutive data points
        for pos, val in data:   
            frag, x = self.canvas_location(pos)
            x += self.x0        # next xy co-ords
            y = trackheight*(val-midval)/resolution + \
                self.fragment_lines[frag][0] + ctr
            if frag == lastfrag:    # Points on the same fragment: draw the line
                line_elements.append(Line(lastx, lasty, x, y,
                                          strokeColor = graph.poscolor,
                                          strokeWidth = graph.linewidth))
            else:   # Points not on the same fragment, so interpolate
                tempval = 1.*(val-lastval)/(x-lastx)
                tempy = trackheight*(val-midval)/resolution + \
                        self.fragment_lines[lastfrag][0] + ctr
                line_elements.append(Line(lastx, lasty, self.xlim, tempy,
                                          strokeColor = graph.poscolor,
                                          strokeWidth = graph.linewidth))
                tempy = trackheight*(val-midval)/resolution + \
                        self.fragment_lines[frag][0] + ctr
                line_elements.append(Line(self.x0, tempy, x, y,
                                          strokeColor = graph.poscolor,
                                          strokeWidth = graph.linewidth))
            lastfrag, lastx, lasty, lastval = frag, x, y, val
            
        return line_elements


    def draw_heat_graph(self, graph):
        """ draw_heat_graph(self, graph) -> [element, element,...]

            o graph     Graph object

            Returns a list of drawable elements for the heat graph
        """
        #print '\tdraw_heat_graph'
        # At each point contained in the graph data, we draw a box that is the
        # full height of the track, extending from the midpoint between the
        # previous and current data points to the midpoint between the current
        # and next data points
        heat_elements = []  # Holds drawable elements for the graph

        # Get graph data and information
        data_quartiles = graph.quartiles()
        minval, maxval = data_quartiles[0],data_quartiles[4]
        midval = (maxval + minval)/2.    # mid is the value at the X-axis
        btm, ctr, top = self.track_offsets[self.current_track_level]
        trackheight = (top-btm)
        #print self.start, self.end
        newdata = intermediate_points(self.start, self.end,
                                      graph[self.start:self.end])
        #print newdata

        # Create elements on the graph, indicating a large positive value by
        # the graph's poscolor, and a large negative value by the graph's
        # negcolor attributes
        for pos0, pos1, val in newdata:
            fragment0, x0 = self.canvas_location(pos0)
            fragment1, x1 = self.canvas_location(pos1)
            x0, x1 = self.x0 + x0, self.x0 + x1     # account for margin
            #print 'x1 before:', x1
            
            # Calculate the heat color, based on the differential between
            # the value and the median value
            heat = colors.linearlyInterpolatedColor(graph.poscolor,
                                                    graph.negcolor,
                                                    maxval, minval, val)
            
            # Draw heat box
            if fragment0 == fragment1:  # Box is contiguous on one fragment
                if pos1 >= self.fragment_limits[fragment0][1]:
                    x1 = self.xlim
                ttop = top + self.fragment_lines[fragment0][0]
                tbtm = btm + self.fragment_lines[fragment0][0]
                #print 'equal', pos0, pos1, val
                #print pos0, pos1, fragment0, fragment1
                heat_elements.append(draw_box((x0, tbtm), (x1, ttop),
                                              color=heat, border=None))
            else:   # box is split over two or more fragments
                #if pos0 >= self.fragment_limits[fragment0][0]:
                #    fragment0 += 1
                fragment = fragment0
                start = x0
                while self.fragment_limits[fragment][1] <= pos1:
                    #print pos0, self.fragment_limits[fragment][1], pos1
                    ttop = top + self.fragment_lines[fragment][0]
                    tbtm = btm + self.fragment_lines[fragment][0]
                    heat_elements.append(draw_box((start, tbtm),
                                                  (self.xlim, ttop),
                                                  color=heat,
                                                  border=None))
                    fragment += 1
                    start = self.x0
                ttop = top + self.fragment_lines[fragment][0]
                tbtm = btm + self.fragment_lines[fragment][0]
                # Add the last part of the bar
                #print 'x1 after:', x1, '\n'
                heat_elements.append(draw_box((self.x0, tbtm), (x1, ttop),
                                              color=heat, border=None))                    

        return heat_elements


    def draw_bar_graph(self, graph):
        """ draw_bar_graph(self, graph) -> [element, element,...]

            o graph     Graph object

            Returns a list of drawable elements for a bar graph of the passed
            Graph object
        """
        #print '\tdraw_bar_graph'
        # At each point contained in the graph data, we draw a vertical bar
        # from the track center to the height of the datapoint value (positive
        # values go up in one color, negative go down in the alternative
        # color).
        bar_elements = []   # Holds drawable elements for the graph
        
        # Set the number of pixels per unit for the data
        data_quartiles = graph.quartiles()
        minval, maxval = data_quartiles[0],data_quartiles[4]
        btm, ctr, top = self.track_offsets[self.current_track_level]
        trackheight = 0.5*(top-btm)
        datarange = maxval - minval
        if datarange == 0:
            datarange = trackheight
        data = graph[self.start:self.end]
        # midval is the value at which the x-axis is plotted, and is the
        # central ring in the track
        if graph.center is None:
            midval = (maxval + minval)/2.    
        else:
            midval = graph.center

        # Convert data into 'binned' blocks, covering half the distance to the
        # next data point on either side, accounting for the ends of fragments
        # and tracks
        newdata = intermediate_points(self.start, self.end,
                                      graph[self.start:self.end])

        # Whichever is the greatest difference: max-midval or min-midval, is
        # taken to specify the number of pixel units resolved along the
        # y-axis
        resolution = max((midval-minval), (maxval-midval))
        if resolution == 0:
            resolution = trackheight
        
        # Create elements for the bar graph based on newdata
        for pos0, pos1, val in newdata:
            fragment0, x0 = self.canvas_location(pos0)
            fragment1, x1 = self.canvas_location(pos1)
            x0, x1 = self.x0 + x0, self.x0 + x1     # account for margin
            barval = trackheight*(val-midval)/resolution
            if barval >=0:  # Different colors for bars that extend above...
                barcolor = graph.poscolor
            else:           # ...or below the axis
                barcolor = graph.negcolor

            # Draw bar
            if fragment0 == fragment1:  # Box is contiguous
                if pos1 >= self.fragment_limits[fragment0][1]:
                    x1 = self.xlim
                tctr = ctr + self.fragment_lines[fragment0][0]
                barval += tctr                
                bar_elements.append(draw_box((x0, tctr), (x1, barval),
                                             color=barcolor))
            else:   # Box is split over two or more fragments                       
                fragment = fragment0
                #if pos0 >= self.fragment_limits[fragment0][0]:
                #    fragment += 1         
                start = x0
                while self.fragment_limits[fragment][1] < pos1:
                    tctr = ctr + self.fragment_lines[fragment][0]
                    thisbarval = barval + tctr
                    bar_elements.append(draw_box((start, tctr),
                                                 (self.xlim, thisbarval),
                                                 color=barcolor))
                    fragment += 1
                    start = self.x0
                tctr = ctr + self.fragment_lines[fragment1][0]
                barval += tctr
                # Add the last part of the bar
                bar_elements.append(draw_box((self.x0, tctr), (x1, barval),
                                             color=barcolor))

        return bar_elements

      
    def canvas_location(self, base):
        """ canvas_location(self, base) -> (int, float)

            o base      The base number on the genome sequence

            Returns the x-coordinate and fragment number of a base on the
            genome sequence, in the context of the current drawing setup
        """
        base = int(base - self.start)   # number of bases we are from the start
        fragment = int(base / self.fragment_bases)
        if fragment < 1:    # First fragment
            base_offset = base
            fragment = 0
        else:               # Calculate number of bases from start of fragment
            base_offset = base % self.fragment_bases
        # Calculate number of pixels from start of fragment
        x_offset = 1. * self.pagewidth * base_offset / self.fragment_bases
        return fragment, x_offset
        
