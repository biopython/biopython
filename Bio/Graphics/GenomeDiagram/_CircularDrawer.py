# Copyright 2003-2008 by Leighton Pritchard.  All rights reserved.
# Revisions copyright 2008-2017 by Peter Cock.
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

"""CircularDrawer module for GenomeDiagram."""

# ReportLab imports

from reportlab.graphics.shapes import Drawing, String, Group, Line, Circle, Polygon
from reportlab.lib import colors
from reportlab.graphics.shapes import ArcPath

# GenomeDiagram imports
from ._AbstractDrawer import AbstractDrawer, draw_polygon, intermediate_points
from ._AbstractDrawer import _stroke_and_fill_colors
from ._FeatureSet import FeatureSet
from ._GraphSet import GraphSet

from math import pi, cos, sin


class CircularDrawer(AbstractDrawer):
    """Object for drawing circular diagrams.

    Attributes:
     - tracklines    Boolean for whether to draw lines dilineating tracks
     - pagesize      Tuple describing the size of the page in pixels
     - x0            Float X co-ord for leftmost point of drawable area
     - xlim          Float X co-ord for rightmost point of drawable area
     - y0            Float Y co-ord for lowest point of drawable area
     - ylim          Float Y co-ord for topmost point of drawable area
     - pagewidth     Float pixel width of drawable area
     - pageheight    Float pixel height of drawable area
     - xcenter       Float X co-ord of center of drawable area
     - ycenter       Float Y co-ord of center of drawable area
     - start         Int, base to start drawing from
     - end           Int, base to stop drawing at
     - length        Size of sequence to be drawn
     - track_size    Float (0->1) the proportion of the track height to draw in
     - drawing       Drawing canvas
     - drawn_tracks  List of ints denoting which tracks are to be drawn
     - current_track_level   Int denoting which track is currently being drawn
     - track_offsets     Dictionary of number of pixels that each track top,
       center and bottom is offset from the base of a fragment, keyed by track
     - sweep     Float (0->1) the proportion of the circle circumference to
       use for the diagram
     - cross_track_links List of tuples each with four entries (track A,
       feature A, track B, feature B) to be linked.

    """

    def __init__(
        self,
        parent=None,
        pagesize="A3",
        orientation="landscape",
        x=0.05,
        y=0.05,
        xl=None,
        xr=None,
        yt=None,
        yb=None,
        start=None,
        end=None,
        tracklines=0,
        track_size=0.75,
        circular=1,
        circle_core=0.0,
        cross_track_links=None,
    ):
        """Create CircularDrawer object.

        Arguments:
         - parent    Diagram object containing the data that the drawer
           draws
         - pagesize  String describing the ISO size of the image, or a tuple
           of pixels
         - orientation   String describing the required orientation of the
           final drawing ('landscape' or 'portrait')
         - x         Float (0->1) describing the relative size of the X
           margins to the page
         - y         Float (0->1) describing the relative size of the Y
           margins to the page
         - xl        Float (0->1) describing the relative size of the left X
           margin to the page (overrides x)
         - xl        Float (0->1) describing the relative size of the left X
           margin to the page (overrides x)
         - xr        Float (0->1) describing the relative size of the right X
           margin to the page (overrides x)
         - yt        Float (0->1) describing the relative size of the top Y
           margin to the page (overrides y)
         - yb        Float (0->1) describing the relative size of the lower Y
           margin to the page (overrides y)
         - start     Int, the position to begin drawing the diagram at
         - end       Int, the position to stop drawing the diagram at
         - tracklines    Boolean flag to show (or not) lines delineating tracks
           on the diagram
         - track_size    The proportion of the available track height that
           should be taken up in drawing
         - circular      Boolean flaw to show whether the passed sequence is
           circular or not
         - circle_core   The proportion of the available radius to leave
           empty at the center of a circular diagram (0 to 1).
         - cross_track_links List of tuples each with four entries (track A,
           feature A, track B, feature B) to be linked.

        """
        # Use the superclass' instantiation method
        AbstractDrawer.__init__(
            self,
            parent,
            pagesize,
            orientation,
            x,
            y,
            xl,
            xr,
            yt,
            yb,
            start,
            end,
            tracklines,
            cross_track_links,
        )

        # Useful measurements on the page
        self.track_size = track_size
        self.circle_core = circle_core
        # Determine proportion of circumference around which information will be drawn
        if not circular:
            self.sweep = 0.9
        else:
            self.sweep = 1.0

    def set_track_heights(self):
        """Initialize track heights.

        Since tracks may not be of identical heights, the bottom and top
        radius for each track is stored in a dictionary - self.track_radii,
        keyed by track number
        """
        bot_track = min(min(self.drawn_tracks), 1)
        top_track = max(self.drawn_tracks)  # The 'highest' track to draw

        trackunit_sum = 0  # Total number of 'units' taken up by all tracks
        trackunits = {}  # Start and & units for each track keyed by track number
        heightholder = 0  # placeholder variable
        for track in range(bot_track, top_track + 1):  # track numbers to 'draw'
            try:
                trackheight = self._parent[track].height  # Get track height
            except Exception:  # TODO: ValueError? IndexError?
                trackheight = 1
            trackunit_sum += trackheight  # increment total track unit height
            trackunits[track] = (heightholder, heightholder + trackheight)
            heightholder += trackheight  # move to next height

        max_radius = 0.5 * min(self.pagewidth, self.pageheight)
        trackunit_height = max_radius * (1 - self.circle_core) / trackunit_sum
        track_core = max_radius * self.circle_core

        # Calculate top and bottom radii for each track
        self.track_radii = {}  # The inner, outer and center radii for each track
        track_crop = (
            trackunit_height * (1 - self.track_size) / 2.0
        )  # 'step back' in pixels
        for track in trackunits:
            top = trackunits[track][1] * trackunit_height - track_crop + track_core
            btm = trackunits[track][0] * trackunit_height + track_crop + track_core
            ctr = btm + (top - btm) / 2.0
            self.track_radii[track] = (btm, ctr, top)

    def draw(self):
        """Draw a circular diagram of the stored data."""
        # Instantiate the drawing canvas
        self.drawing = Drawing(self.pagesize[0], self.pagesize[1])

        feature_elements = []  # holds feature elements
        feature_labels = []  # holds feature labels
        greytrack_bgs = []  # holds track background
        greytrack_labels = []  # holds track foreground labels
        scale_axes = []  # holds scale axes
        scale_labels = []  # holds scale axis labels

        # Get tracks to be drawn and set track sizes
        self.drawn_tracks = self._parent.get_drawn_levels()
        self.set_track_heights()

        # Go through each track in the parent (if it is to be drawn) one by
        # one and collate the data as drawing elements
        for track_level in self._parent.get_drawn_levels():
            self.current_track_level = track_level
            track = self._parent[track_level]
            gbgs, glabels = self.draw_greytrack(track)  # Greytracks
            greytrack_bgs.append(gbgs)
            greytrack_labels.append(glabels)
            features, flabels = self.draw_track(track)  # Features and graphs
            feature_elements.append(features)
            feature_labels.append(flabels)
            if track.scale:
                axes, slabels = self.draw_scale(track)  # Scale axes
                scale_axes.append(axes)
                scale_labels.append(slabels)

        feature_cross_links = []
        for cross_link_obj in self.cross_track_links:
            cross_link_elements = self.draw_cross_link(cross_link_obj)
            if cross_link_elements:
                feature_cross_links.append(cross_link_elements)

        # Groups listed in order of addition to page (from back to front)
        # Draw track backgrounds
        # Draw feature cross track links
        # Draw features and graphs
        # Draw scale axes
        # Draw scale labels
        # Draw feature labels
        # Draw track labels
        element_groups = [
            greytrack_bgs,
            feature_cross_links,
            feature_elements,
            scale_axes,
            scale_labels,
            feature_labels,
            greytrack_labels,
        ]
        for element_group in element_groups:
            for element_list in element_group:
                [self.drawing.add(element) for element in element_list]

        if self.tracklines:
            # Draw test tracks over top of diagram
            self.draw_test_tracks()

    def draw_track(self, track):
        """Return list of track elements and list of track labels."""
        track_elements = []  # Holds elements for features and graphs
        track_labels = []  # Holds labels for features and graphs

        # Distribution dictionary for dealing with different set types
        set_methods = {FeatureSet: self.draw_feature_set, GraphSet: self.draw_graph_set}

        for set in track.get_sets():  # Draw the feature or graph sets
            elements, labels = set_methods[set.__class__](set)
            track_elements += elements
            track_labels += labels
        return track_elements, track_labels

    def draw_feature_set(self, set):
        """Return list of feature elements and list of labels for them."""
        # print 'draw feature set'
        feature_elements = []  # Holds diagram elements belonging to the features
        label_elements = []  # Holds diagram elements belonging to feature labels

        # Collect all the elements for the feature set
        for feature in set.get_features():
            if self.is_in_bounds(feature.start) or self.is_in_bounds(feature.end):
                features, labels = self.draw_feature(feature)
                feature_elements += features
                label_elements += labels

        return feature_elements, label_elements

    def draw_feature(self, feature):
        """Return list of feature elements and list of labels for them."""
        feature_elements = []  # Holds drawable elements for a single feature
        label_elements = []  # Holds labels for a single feature

        if feature.hide:  # Don't show feature: return early
            return feature_elements, label_elements

        start, end = self._current_track_start_end()
        # A single feature may be split into subfeatures, so loop over them
        for locstart, locend in feature.locations:
            if locend < start:
                continue
            locstart = max(locstart, start)
            if end < locstart:
                continue
            locend = min(locend, end)
            # Get sigil for the feature/ each subfeature
            feature_sigil, label = self.get_feature_sigil(feature, locstart, locend)
            feature_elements.append(feature_sigil)
            if label is not None:  # If there's a label
                label_elements.append(label)

        return feature_elements, label_elements

    def get_feature_sigil(self, feature, locstart, locend, **kwargs):
        """Return graphics for feature, and any required label for it.

        Arguments:
         - feature       Feature object
         - locstart      The start position of the feature
         - locend        The end position of the feature

        """
        # Establish the co-ordinates for the sigil
        btm, ctr, top = self.track_radii[self.current_track_level]

        startangle, startcos, startsin = self.canvas_angle(locstart)
        endangle, endcos, endsin = self.canvas_angle(locend)
        midangle, midcos, midsin = self.canvas_angle(float(locend + locstart) / 2)

        # Distribution dictionary for various ways of drawing the feature
        # Each method takes the inner and outer radii, the start and end angle
        # subtended at the diagram center, and the color as arguments
        draw_methods = {
            "BOX": self._draw_sigil_box,
            "OCTO": self._draw_sigil_cut_corner_box,
            "JAGGY": self._draw_sigil_jaggy,
            "ARROW": self._draw_sigil_arrow,
            "BIGARROW": self._draw_sigil_big_arrow,
        }

        # Get sigil for the feature, location dependent on the feature strand
        method = draw_methods[feature.sigil]
        kwargs["head_length_ratio"] = feature.arrowhead_length
        kwargs["shaft_height_ratio"] = feature.arrowshaft_height

        # Support for clickable links... needs ReportLab 2.4 or later
        # which added support for links in SVG output.
        if hasattr(feature, "url"):
            kwargs["hrefURL"] = feature.url
            kwargs["hrefTitle"] = feature.name

        sigil = method(
            btm,
            ctr,
            top,
            startangle,
            endangle,
            feature.strand,
            color=feature.color,
            border=feature.border,
            **kwargs
        )

        if feature.label:  # Feature needs a label
            # The spaces are a hack to force a little space between the label
            # and the edge of the feature
            label = String(
                0,
                0,
                " %s " % feature.name.strip(),
                fontName=feature.label_font,
                fontSize=feature.label_size,
                fillColor=feature.label_color,
            )
            labelgroup = Group(label)
            if feature.label_strand:
                strand = feature.label_strand
            else:
                strand = feature.strand
            if feature.label_position in ("start", "5'", "left"):
                # Position the label at the feature's start
                if strand != -1:
                    label_angle = startangle + 0.5 * pi  # Make text radial
                    sinval, cosval = startsin, startcos
                else:
                    label_angle = endangle + 0.5 * pi  # Make text radial
                    sinval, cosval = endsin, endcos
            elif feature.label_position in ("middle", "center", "centre"):
                # Position the label at the feature's midpoint
                label_angle = midangle + 0.5 * pi  # Make text radial
                sinval, cosval = midsin, midcos
            elif feature.label_position in ("end", "3'", "right"):
                # Position the label at the feature's end
                if strand != -1:
                    label_angle = endangle + 0.5 * pi  # Make text radial
                    sinval, cosval = endsin, endcos
                else:
                    label_angle = startangle + 0.5 * pi  # Make text radial
                    sinval, cosval = startsin, startcos
            elif startangle < pi:
                # Default to placing the label the bottom of the feature
                # as drawn on the page, meaning feature end on left half
                label_angle = endangle + 0.5 * pi  # Make text radial
                sinval, cosval = endsin, endcos
            else:
                # Default to placing the label on the bottom of the feature,
                # which means the feature end when on right hand half
                label_angle = startangle + 0.5 * pi  # Make text radial
                sinval, cosval = startsin, startcos
            if strand != -1:
                # Feature label on top
                radius = top
                if startangle < pi:  # Turn text round
                    label_angle -= pi
                else:
                    labelgroup.contents[0].textAnchor = "end"
            else:
                # Feature label on bottom
                radius = btm
                if startangle < pi:  # Turn text round and anchor end
                    label_angle -= pi
                    labelgroup.contents[0].textAnchor = "end"
            x_pos = self.xcenter + radius * sinval
            y_pos = self.ycenter + radius * cosval
            coslabel = cos(label_angle)
            sinlabel = sin(label_angle)
            labelgroup.transform = (
                coslabel,
                -sinlabel,
                sinlabel,
                coslabel,
                x_pos,
                y_pos,
            )
        else:
            # No label required
            labelgroup = None
        # if locstart > locend:
        #    print locstart, locend, feature.strand, sigil, feature.name
        # print locstart, locend, feature.name
        return sigil, labelgroup

    def draw_cross_link(self, cross_link):
        """Draw a cross-link between features."""
        startA = cross_link.startA
        startB = cross_link.startB
        endA = cross_link.endA
        endB = cross_link.endB

        if not self.is_in_bounds(startA) and not self.is_in_bounds(endA):
            return None
        if not self.is_in_bounds(startB) and not self.is_in_bounds(endB):
            return None

        if startA < self.start:
            startA = self.start
        if startB < self.start:
            startB = self.start
        if self.end < endA:
            endA = self.end
        if self.end < endB:
            endB = self.end

        trackobjA = cross_link._trackA(list(self._parent.tracks.values()))
        trackobjB = cross_link._trackB(list(self._parent.tracks.values()))
        assert trackobjA is not None
        assert trackobjB is not None
        if trackobjA == trackobjB:
            raise NotImplementedError()

        if trackobjA.start is not None:
            if endA < trackobjA.start:
                return
            startA = max(startA, trackobjA.start)
        if trackobjA.end is not None:
            if trackobjA.end < startA:
                return
            endA = min(endA, trackobjA.end)
        if trackobjB.start is not None:
            if endB < trackobjB.start:
                return
            startB = max(startB, trackobjB.start)
        if trackobjB.end is not None:
            if trackobjB.end < startB:
                return
            endB = min(endB, trackobjB.end)

        for track_level in self._parent.get_drawn_levels():
            track = self._parent[track_level]
            if track == trackobjA:
                trackA = track_level
            if track == trackobjB:
                trackB = track_level
        if trackA == trackB:
            raise NotImplementedError()

        startangleA, startcosA, startsinA = self.canvas_angle(startA)
        startangleB, startcosB, startsinB = self.canvas_angle(startB)
        endangleA, endcosA, endsinA = self.canvas_angle(endA)
        endangleB, endcosB, endsinB = self.canvas_angle(endB)

        btmA, ctrA, topA = self.track_radii[trackA]
        btmB, ctrB, topB = self.track_radii[trackB]

        if ctrA < ctrB:
            return [
                self._draw_arc_poly(
                    topA,
                    btmB,
                    startangleA,
                    endangleA,
                    startangleB,
                    endangleB,
                    cross_link.color,
                    cross_link.border,
                    cross_link.flip,
                )
            ]
        else:
            return [
                self._draw_arc_poly(
                    btmA,
                    topB,
                    startangleA,
                    endangleA,
                    startangleB,
                    endangleB,
                    cross_link.color,
                    cross_link.border,
                    cross_link.flip,
                )
            ]

    def draw_graph_set(self, set):
        """Return list of graph elements and list of their labels.

        Arguments:
         - set       GraphSet object

        """
        # print 'draw graph set'
        elements = []  # Holds graph elements

        # Distribution dictionary for how to draw the graph
        style_methods = {
            "line": self.draw_line_graph,
            "heat": self.draw_heat_graph,
            "bar": self.draw_bar_graph,
        }

        for graph in set.get_graphs():
            elements += style_methods[graph.style](graph)

        return elements, []

    def draw_line_graph(self, graph):
        """Return line graph as list of drawable elements.

        Arguments:
         - graph     GraphData object

        """
        line_elements = []  # holds drawable elements

        # Get graph data
        data_quartiles = graph.quartiles()
        minval, maxval = data_quartiles[0], data_quartiles[4]
        btm, ctr, top = self.track_radii[self.current_track_level]
        trackheight = 0.5 * (top - btm)
        datarange = maxval - minval
        if datarange == 0:
            datarange = trackheight

        start, end = self._current_track_start_end()
        data = graph[start:end]

        if not data:
            return []

        # midval is the value at which the x-axis is plotted, and is the
        # central ring in the track
        if graph.center is None:
            midval = (maxval + minval) / 2.0
        else:
            midval = graph.center
        # Whichever is the greatest difference: max-midval or min-midval, is
        # taken to specify the number of pixel units resolved along the
        # y-axis
        resolution = max((midval - minval), (maxval - midval))

        # Start from first data point
        pos, val = data[0]
        lastangle, lastcos, lastsin = self.canvas_angle(pos)
        # We calculate the track height
        posheight = trackheight * (val - midval) / resolution + ctr
        lastx = self.xcenter + posheight * lastsin  # start xy coords
        lasty = self.ycenter + posheight * lastcos
        for pos, val in data:
            posangle, poscos, possin = self.canvas_angle(pos)
            posheight = trackheight * (val - midval) / resolution + ctr
            x = self.xcenter + posheight * possin  # next xy coords
            y = self.ycenter + posheight * poscos
            line_elements.append(
                Line(
                    lastx,
                    lasty,
                    x,
                    y,
                    strokeColor=graph.poscolor,
                    strokeWidth=graph.linewidth,
                )
            )
            lastx, lasty, = x, y
        return line_elements

    def draw_bar_graph(self, graph):
        """Return list of drawable elements for a bar graph.

        Arguments:
         - graph     Graph object

        """
        # At each point contained in the graph data, we draw a vertical bar
        # from the track center to the height of the datapoint value (positive
        # values go up in one color, negative go down in the alternative
        # color).
        bar_elements = []

        # Set the number of pixels per unit for the data
        data_quartiles = graph.quartiles()
        minval, maxval = data_quartiles[0], data_quartiles[4]
        btm, ctr, top = self.track_radii[self.current_track_level]
        trackheight = 0.5 * (top - btm)
        datarange = maxval - minval
        if datarange == 0:
            datarange = trackheight
        data = graph[self.start : self.end]
        # midval is the value at which the x-axis is plotted, and is the
        # central ring in the track
        if graph.center is None:
            midval = (maxval + minval) / 2.0
        else:
            midval = graph.center

        # Convert data into 'binned' blocks, covering half the distance to the
        # next data point on either side, accounting for the ends of fragments
        # and tracks
        start, end = self._current_track_start_end()
        data = intermediate_points(start, end, graph[start:end])

        if not data:
            return []

        # Whichever is the greatest difference: max-midval or min-midval, is
        # taken to specify the number of pixel units resolved along the
        # y-axis
        resolution = max((midval - minval), (maxval - midval))
        if resolution == 0:
            resolution = trackheight

        # Create elements for the bar graph based on newdata
        for pos0, pos1, val in data:
            pos0angle, pos0cos, pos0sin = self.canvas_angle(pos0)
            pos1angle, pos1cos, pos1sin = self.canvas_angle(pos1)

            barval = trackheight * (val - midval) / resolution
            if barval >= 0:
                barcolor = graph.poscolor
            else:
                barcolor = graph.negcolor

            # Draw bar
            bar_elements.append(
                self._draw_arc(ctr, ctr + barval, pos0angle, pos1angle, barcolor)
            )
        return bar_elements

    def draw_heat_graph(self, graph):
        """Return list of drawable elements for the heat graph.

        Arguments:
         - graph     Graph object

        """
        # At each point contained in the graph data, we draw a box that is the
        # full height of the track, extending from the midpoint between the
        # previous and current data points to the midpoint between the current
        # and next data points
        heat_elements = []  # holds drawable elements

        # Get graph data
        data_quartiles = graph.quartiles()
        minval, maxval = data_quartiles[0], data_quartiles[4]
        midval = (maxval + minval) / 2.0  # mid is the value at the X-axis
        btm, ctr, top = self.track_radii[self.current_track_level]
        trackheight = top - btm

        start, end = self._current_track_start_end()
        data = intermediate_points(start, end, graph[start:end])

        # Create elements on the graph, indicating a large positive value by
        # the graph's poscolor, and a large negative value by the graph's
        # negcolor attributes
        for pos0, pos1, val in data:
            pos0angle, pos0cos, pos0sin = self.canvas_angle(pos0)
            pos1angle, pos1cos, pos1sin = self.canvas_angle(pos1)

            # Calculate the heat color, based on the differential between
            # the value and the median value
            heat = colors.linearlyInterpolatedColor(
                graph.poscolor, graph.negcolor, maxval, minval, val
            )

            # Draw heat box
            heat_elements.append(
                self._draw_arc(btm, top, pos0angle, pos1angle, heat, border=heat)
            )
        return heat_elements

    def draw_scale(self, track):
        """Return list of elements in the scale and list of their labels.

        Arguments:
         - track     Track object

        """
        scale_elements = []  # holds axes and ticks
        scale_labels = []  # holds labels

        if not track.scale:
            # no scale required, exit early
            return [], []

        # Get track locations
        btm, ctr, top = self.track_radii[self.current_track_level]
        trackheight = top - ctr

        # X-axis
        start, end = self._current_track_start_end()
        if track.start is not None or track.end is not None:
            # Draw an arc, leaving out the wedge
            p = ArcPath(strokeColor=track.scale_color, fillColor=None)
            startangle, startcos, startsin = self.canvas_angle(start)
            endangle, endcos, endsin = self.canvas_angle(end)
            p.addArc(
                self.xcenter,
                self.ycenter,
                ctr,
                90 - (endangle * 180 / pi),
                90 - (startangle * 180 / pi),
            )
            scale_elements.append(p)
            del p
            # Y-axis start marker
            x0, y0 = self.xcenter + btm * startsin, self.ycenter + btm * startcos
            x1, y1 = self.xcenter + top * startsin, self.ycenter + top * startcos
            scale_elements.append(Line(x0, y0, x1, y1, strokeColor=track.scale_color))
            # Y-axis end marker
            x0, y0 = self.xcenter + btm * endsin, self.ycenter + btm * endcos
            x1, y1 = self.xcenter + top * endsin, self.ycenter + top * endcos
            scale_elements.append(Line(x0, y0, x1, y1, strokeColor=track.scale_color))
        elif self.sweep < 1:
            # Draw an arc, leaving out the wedge
            p = ArcPath(strokeColor=track.scale_color, fillColor=None)
            # Note reportlab counts angles anti-clockwise from the horizontal
            # (as in mathematics, e.g. complex numbers and polar coordinates)
            # in degrees.
            p.addArc(
                self.xcenter,
                self.ycenter,
                ctr,
                startangledegrees=90 - 360 * self.sweep,
                endangledegrees=90,
            )
            scale_elements.append(p)
            del p
            # Y-axis start marker
            x0, y0 = self.xcenter, self.ycenter + btm
            x1, y1 = self.xcenter, self.ycenter + top
            scale_elements.append(Line(x0, y0, x1, y1, strokeColor=track.scale_color))
            # Y-axis end marker
            alpha = 2 * pi * self.sweep
            x0, y0 = self.xcenter + btm * sin(alpha), self.ycenter + btm * cos(alpha)
            x1, y1 = self.xcenter + top * sin(alpha), self.ycenter + top * cos(alpha)
            scale_elements.append(Line(x0, y0, x1, y1, strokeColor=track.scale_color))
        else:
            # Draw a full circle
            scale_elements.append(
                Circle(
                    self.xcenter,
                    self.ycenter,
                    ctr,
                    strokeColor=track.scale_color,
                    fillColor=None,
                )
            )

        start, end = self._current_track_start_end()
        if track.scale_ticks:  # Ticks are required on the scale
            # Draw large ticks
            # I want the ticks to be consistently positioned relative to
            # the start of the sequence (position 0), not relative to the
            # current viewpoint (self.start and self.end)

            ticklen = track.scale_largeticks * trackheight
            tickiterval = int(track.scale_largetick_interval)
            # Note that we could just start the list of ticks using
            # range(0,self.end,tickinterval) and the filter out the
            # ones before self.start - but this seems wasteful.
            # Using tickiterval * (self.start/tickiterval) is a shortcut.
            for tickpos in range(
                tickiterval * (self.start // tickiterval), int(self.end), tickiterval
            ):
                if tickpos <= start or end <= tickpos:
                    continue
                tick, label = self.draw_tick(
                    tickpos, ctr, ticklen, track, track.scale_largetick_labels
                )
                scale_elements.append(tick)
                if label is not None:  # If there's a label, add it
                    scale_labels.append(label)
            # Draw small ticks
            ticklen = track.scale_smallticks * trackheight
            tickiterval = int(track.scale_smalltick_interval)
            for tickpos in range(
                tickiterval * (self.start // tickiterval), int(self.end), tickiterval
            ):
                if tickpos <= start or end <= tickpos:
                    continue
                tick, label = self.draw_tick(
                    tickpos, ctr, ticklen, track, track.scale_smalltick_labels
                )
                scale_elements.append(tick)
                if label is not None:  # If there's a label, add it
                    scale_labels.append(label)

        # Check to see if the track contains a graph - if it does, get the
        # minimum and maximum values, and put them on the scale Y-axis
        # at 60 degree intervals, ordering the labels by graph_id
        startangle, startcos, startsin = self.canvas_angle(start)
        endangle, endcos, endsin = self.canvas_angle(end)
        if track.axis_labels:
            for set in track.get_sets():
                if set.__class__ is GraphSet:
                    # Y-axis
                    for n in range(7):
                        angle = n * 1.0471975511965976
                        if angle < startangle or endangle < angle:
                            continue
                        ticksin, tickcos = sin(angle), cos(angle)
                        x0, y0 = (
                            self.xcenter + btm * ticksin,
                            self.ycenter + btm * tickcos,
                        )
                        x1, y1 = (
                            self.xcenter + top * ticksin,
                            self.ycenter + top * tickcos,
                        )
                        scale_elements.append(
                            Line(x0, y0, x1, y1, strokeColor=track.scale_color)
                        )

                        graph_label_min = []
                        graph_label_max = []
                        graph_label_mid = []
                        for graph in set.get_graphs():
                            quartiles = graph.quartiles()
                            minval, maxval = quartiles[0], quartiles[4]
                            if graph.center is None:
                                midval = (maxval + minval) / 2.0
                                graph_label_min.append("%.3f" % minval)
                                graph_label_max.append("%.3f" % maxval)
                                graph_label_mid.append("%.3f" % midval)
                            else:
                                diff = max(
                                    (graph.center - minval), (maxval - graph.center)
                                )
                                minval = graph.center - diff
                                maxval = graph.center + diff
                                midval = graph.center
                                graph_label_mid.append("%.3f" % midval)
                                graph_label_min.append("%.3f" % minval)
                                graph_label_max.append("%.3f" % maxval)
                        xmid, ymid = (x0 + x1) / 2.0, (y0 + y1) / 2.0
                        for limit, x, y in [
                            (graph_label_min, x0, y0),
                            (graph_label_max, x1, y1),
                            (graph_label_mid, xmid, ymid),
                        ]:
                            label = String(
                                0,
                                0,
                                ";".join(limit),
                                fontName=track.scale_font,
                                fontSize=track.scale_fontsize,
                                fillColor=track.scale_color,
                            )
                            label.textAnchor = "middle"
                            labelgroup = Group(label)
                            labelgroup.transform = (
                                tickcos,
                                -ticksin,
                                ticksin,
                                tickcos,
                                x,
                                y,
                            )
                            scale_labels.append(labelgroup)

        return scale_elements, scale_labels

    def draw_tick(self, tickpos, ctr, ticklen, track, draw_label):
        """Return drawing element for a tick on the scale.

        Arguments:
         - tickpos   Int, position of the tick on the sequence
         - ctr       Float, Y co-ord of the center of the track
         - ticklen   How long to draw the tick
         - track     Track, the track the tick is drawn on
         - draw_label    Boolean, write the tick label?

        """
        # Calculate tick co-ordinates
        tickangle, tickcos, ticksin = self.canvas_angle(tickpos)
        x0, y0 = self.xcenter + ctr * ticksin, self.ycenter + ctr * tickcos
        x1, y1 = (
            self.xcenter + (ctr + ticklen) * ticksin,
            self.ycenter + (ctr + ticklen) * tickcos,
        )
        # Calculate height of text label so it can be offset on lower half
        # of diagram
        # LP: not used, as not all fonts have ascent_descent data in reportlab.pdfbase._fontdata
        # label_offset = _fontdata.ascent_descent[track.scale_font][0]*\
        #               track.scale_fontsize/1000.
        tick = Line(x0, y0, x1, y1, strokeColor=track.scale_color)
        if draw_label:
            # Put tick position on as label
            if track.scale_format == "SInt":
                if tickpos >= 1000000:
                    tickstring = str(tickpos // 1000000) + " Mbp"
                elif tickpos >= 1000:
                    tickstring = str(tickpos // 1000) + " Kbp"
                else:
                    tickstring = str(tickpos)
            else:
                tickstring = str(tickpos)
            label = String(
                0,
                0,
                tickstring,  # Make label string
                fontName=track.scale_font,
                fontSize=track.scale_fontsize,
                fillColor=track.scale_color,
            )
            if tickangle > pi:
                label.textAnchor = "end"
            # LP: This label_offset depends on ascent_descent data, which is not available for all
            # fonts, so has been deprecated.
            # if 0.5*pi < tickangle < 1.5*pi:
            #    y1 -= label_offset
            labelgroup = Group(label)
            labelgroup.transform = (1, 0, 0, 1, x1, y1)
        else:
            labelgroup = None
        return tick, labelgroup

    def draw_test_tracks(self):
        """Draw blue test tracks with grene line down their center."""
        # Add lines only for drawn tracks
        for track in self.drawn_tracks:
            btm, ctr, top = self.track_radii[track]
            self.drawing.add(
                Circle(
                    self.xcenter,
                    self.ycenter,
                    top,
                    strokeColor=colors.blue,
                    fillColor=None,
                )
            )  # top line
            self.drawing.add(
                Circle(
                    self.xcenter,
                    self.ycenter,
                    ctr,
                    strokeColor=colors.green,
                    fillColor=None,
                )
            )  # middle line
            self.drawing.add(
                Circle(
                    self.xcenter,
                    self.ycenter,
                    btm,
                    strokeColor=colors.blue,
                    fillColor=None,
                )
            )  # bottom line

    def draw_greytrack(self, track):
        """Drawing element for grey background to passed Track object."""
        greytrack_bgs = []  # Holds track backgrounds
        greytrack_labels = []  # Holds track foreground labels

        if not track.greytrack:  # No greytrack required, return early
            return [], []

        # Get track location
        btm, ctr, top = self.track_radii[self.current_track_level]

        start, end = self._current_track_start_end()
        startangle, startcos, startsin = self.canvas_angle(start)
        endangle, endcos, endsin = self.canvas_angle(end)

        # Make background
        if track.start is not None or track.end is not None:
            # Draw an arc, leaving out the wedge
            p = ArcPath(strokeColor=track.scale_color, fillColor=None)
            greytrack_bgs.append(
                self._draw_arc(
                    btm, top, startangle, endangle, colors.Color(0.96, 0.96, 0.96)
                )
            )
        elif self.sweep < 1:
            # Make a partial circle, a large arc box
            # This method assumes the correct center for us.
            greytrack_bgs.append(
                self._draw_arc(
                    btm, top, 0, 2 * pi * self.sweep, colors.Color(0.96, 0.96, 0.96)
                )
            )
        else:
            # Make a full circle (using a VERY thick linewidth)
            greytrack_bgs.append(
                Circle(
                    self.xcenter,
                    self.ycenter,
                    ctr,
                    strokeColor=colors.Color(0.96, 0.96, 0.96),
                    fillColor=None,
                    strokeWidth=top - btm,
                )
            )

        if track.greytrack_labels:
            # Labels are required for this track
            labelstep = self.length // track.greytrack_labels  # label interval
            for pos in range(self.start, self.end, labelstep):
                label = String(
                    0,
                    0,
                    track.name,  # Add a new label at
                    fontName=track.greytrack_font,  # each interval
                    fontSize=track.greytrack_fontsize,
                    fillColor=track.greytrack_fontcolor,
                )
                theta, costheta, sintheta = self.canvas_angle(pos)
                if theta < startangle or endangle < theta:
                    continue
                x, y = (
                    self.xcenter + btm * sintheta,
                    self.ycenter + btm * costheta,
                )  # start text halfway up marker
                labelgroup = Group(label)
                labelangle = (
                    self.sweep * 2 * pi * (pos - self.start) / self.length - pi / 2
                )
                if theta > pi:
                    label.textAnchor = "end"  # Anchor end of text to inner radius
                    labelangle += pi  # and reorient it
                cosA, sinA = cos(labelangle), sin(labelangle)
                labelgroup.transform = (cosA, -sinA, sinA, cosA, x, y)
                if not self.length - x <= labelstep:  # Don't overrun the circle
                    greytrack_labels.append(labelgroup)

        return greytrack_bgs, greytrack_labels

    def canvas_angle(self, base):
        """Given base-pair position, return (angle, cosine, sin) (PRIVATE)."""
        angle = self.sweep * 2 * pi * (base - self.start) / self.length
        return (angle, cos(angle), sin(angle))

    def _draw_sigil_box(
        self, bottom, center, top, startangle, endangle, strand, **kwargs
    ):
        """Draw BOX sigil (PRIVATE)."""
        if strand == 1:
            inner_radius = center
            outer_radius = top
        elif strand == -1:
            inner_radius = bottom
            outer_radius = center
        else:
            inner_radius = bottom
            outer_radius = top
        return self._draw_arc(
            inner_radius, outer_radius, startangle, endangle, **kwargs
        )

    def _draw_arc(
        self,
        inner_radius,
        outer_radius,
        startangle,
        endangle,
        color,
        border=None,
        colour=None,
        **kwargs
    ):
        """Return closed path describing an arc box (PRIVATE).

        Arguments:
         - inner_radius  Float distance of inside of arc from drawing center
         - outer_radius  Float distance of outside of arc from drawing center
         - startangle    Float angle subtended by start of arc at drawing center
           (in radians)
         - endangle      Float angle subtended by end of arc at drawing center
           (in radians)
         - color        colors.Color object for arc (overridden by backwards
           compatible argument with UK spelling, colour).

        Returns a closed path object describing an arced box corresponding to
        the passed values.  For very small angles, a simple four sided
        polygon is used.
        """
        # Let the UK spelling (colour) override the USA spelling (color)
        if colour is not None:
            color = colour

        strokecolor, color = _stroke_and_fill_colors(color, border)

        if abs(float(endangle - startangle)) > 0.01:
            # Wide arc, must use full curves
            p = ArcPath(strokeColor=strokecolor, fillColor=color, strokewidth=0)
            # Note reportlab counts angles anti-clockwise from the horizontal
            # (as in mathematics, e.g. complex numbers and polar coordinates)
            # but we use clockwise from the vertical.  Also reportlab uses
            # degrees, but we use radians.
            p.addArc(
                self.xcenter,
                self.ycenter,
                inner_radius,
                90 - (endangle * 180 / pi),
                90 - (startangle * 180 / pi),
                moveTo=True,
            )
            p.addArc(
                self.xcenter,
                self.ycenter,
                outer_radius,
                90 - (endangle * 180 / pi),
                90 - (startangle * 180 / pi),
                reverse=True,
            )
            p.closePath()
            return p
        else:
            # Cheat and just use a four sided polygon.
            # Calculate trig values for angle and coordinates
            startcos, startsin = cos(startangle), sin(startangle)
            endcos, endsin = cos(endangle), sin(endangle)
            x0, y0 = self.xcenter, self.ycenter  # origin of the circle
            x1, y1 = (x0 + inner_radius * startsin, y0 + inner_radius * startcos)
            x2, y2 = (x0 + inner_radius * endsin, y0 + inner_radius * endcos)
            x3, y3 = (x0 + outer_radius * endsin, y0 + outer_radius * endcos)
            x4, y4 = (x0 + outer_radius * startsin, y0 + outer_radius * startcos)
            return draw_polygon([(x1, y1), (x2, y2), (x3, y3), (x4, y4)], color, border)

    def _draw_arc_line(
        self, path, start_radius, end_radius, start_angle, end_angle, move=False
    ):
        """Add a list of points to a path object (PRIVATE).

        Assumes angles given are in degrees!

        Represents what would be a straight line on a linear diagram.
        """
        x0, y0 = self.xcenter, self.ycenter  # origin of the circle
        radius_diff = end_radius - start_radius
        angle_diff = end_angle - start_angle
        dx = 0.01  # heuristic
        a = start_angle * pi / 180
        if move:
            path.moveTo(x0 + start_radius * cos(a), y0 + start_radius * sin(a))
        else:
            path.lineTo(x0 + start_radius * cos(a), y0 + start_radius * sin(a))
        x = dx
        if 0.01 <= abs(dx):
            while x < 1:
                r = start_radius + x * radius_diff
                a = (
                    (start_angle + x * (angle_diff)) * pi / 180
                )  # to radians for sin/cos
                # print x0+r*cos(a), y0+r*sin(a)
                path.lineTo(x0 + r * cos(a), y0 + r * sin(a))
                x += dx
        a = end_angle * pi / 180
        path.lineTo(x0 + end_radius * cos(a), y0 + end_radius * sin(a))

    def _draw_arc_poly(
        self,
        inner_radius,
        outer_radius,
        inner_startangle,
        inner_endangle,
        outer_startangle,
        outer_endangle,
        color,
        border=None,
        flip=False,
        **kwargs
    ):
        """Return polygon path describing an arc."""
        strokecolor, color = _stroke_and_fill_colors(color, border)

        x0, y0 = self.xcenter, self.ycenter  # origin of the circle
        if (
            abs(inner_endangle - outer_startangle) > 0.01
            or abs(outer_endangle - inner_startangle) > 0.01
            or abs(inner_startangle - outer_startangle) > 0.01
            or abs(outer_startangle - outer_startangle) > 0.01
        ):
            # Wide arc, must use full curves
            p = ArcPath(
                strokeColor=strokecolor,
                fillColor=color,
                # default is mitre/miter which can stick out too much:
                strokeLineJoin=1,  # 1=round
                strokewidth=0,
            )
            # Note reportlab counts angles anti-clockwise from the horizontal
            # (as in mathematics, e.g. complex numbers and polar coordinates)
            # but we use clockwise from the vertical.  Also reportlab uses
            # degrees, but we use radians.
            i_start = 90 - (inner_startangle * 180 / pi)
            i_end = 90 - (inner_endangle * 180 / pi)
            o_start = 90 - (outer_startangle * 180 / pi)
            o_end = 90 - (outer_endangle * 180 / pi)
            p.addArc(x0, y0, inner_radius, i_end, i_start, moveTo=True, reverse=True)
            if flip:
                # Flipped, join end to start,
                self._draw_arc_line(p, inner_radius, outer_radius, i_end, o_start)
                p.addArc(x0, y0, outer_radius, o_end, o_start, reverse=True)
                self._draw_arc_line(p, outer_radius, inner_radius, o_end, i_start)
            else:
                # Not flipped, join start to start, end to end
                self._draw_arc_line(p, inner_radius, outer_radius, i_end, o_end)
                p.addArc(x0, y0, outer_radius, o_end, o_start, reverse=False)
                self._draw_arc_line(p, outer_radius, inner_radius, o_start, i_start)
            p.closePath()
            return p
        else:
            # Cheat and just use a four sided polygon.
            # Calculate trig values for angle and coordinates
            inner_startcos, inner_startsin = (
                cos(inner_startangle),
                sin(inner_startangle),
            )
            inner_endcos, inner_endsin = cos(inner_endangle), sin(inner_endangle)
            outer_startcos, outer_startsin = (
                cos(outer_startangle),
                sin(outer_startangle),
            )
            outer_endcos, outer_endsin = cos(outer_endangle), sin(outer_endangle)
            x1, y1 = (
                x0 + inner_radius * inner_startsin,
                y0 + inner_radius * inner_startcos,
            )
            x2, y2 = (
                x0 + inner_radius * inner_endsin,
                y0 + inner_radius * inner_endcos,
            )
            x3, y3 = (
                x0 + outer_radius * outer_endsin,
                y0 + outer_radius * outer_endcos,
            )
            x4, y4 = (
                x0 + outer_radius * outer_startsin,
                y0 + outer_radius * outer_startcos,
            )
            return draw_polygon(
                [(x1, y1), (x2, y2), (x3, y3), (x4, y4)],
                color,
                border,
                # default is mitre/miter which can stick out too much:
                strokeLineJoin=1,  # 1=round
            )

    def _draw_sigil_cut_corner_box(
        self,
        bottom,
        center,
        top,
        startangle,
        endangle,
        strand,
        color,
        border=None,
        corner=0.5,
        **kwargs
    ):
        """Draw OCTO sigil, box with corners cut off (PRIVATE)."""
        if strand == 1:
            inner_radius = center
            outer_radius = top
        elif strand == -1:
            inner_radius = bottom
            outer_radius = center
        else:
            inner_radius = bottom
            outer_radius = top

        strokecolor, color = _stroke_and_fill_colors(color, border)

        startangle, endangle = min(startangle, endangle), max(startangle, endangle)
        angle = float(endangle - startangle)

        middle_radius = 0.5 * (inner_radius + outer_radius)
        boxheight = outer_radius - inner_radius

        corner_len = min(0.5 * boxheight, 0.5 * boxheight * corner)
        shaft_inner_radius = inner_radius + corner_len
        shaft_outer_radius = outer_radius - corner_len

        cornerangle_delta = max(
            0.0, min(abs(boxheight) * 0.5 * corner / middle_radius, abs(angle * 0.5))
        )
        if angle < 0:
            cornerangle_delta *= -1  # reverse it

        # Calculate trig values for angle and coordinates
        startcos, startsin = cos(startangle), sin(startangle)
        endcos, endsin = cos(endangle), sin(endangle)
        x0, y0 = self.xcenter, self.ycenter  # origin of the circle
        p = ArcPath(
            strokeColor=strokecolor,
            fillColor=color,
            strokeLineJoin=1,  # 1=round
            strokewidth=0,
            **kwargs
        )
        # Inner curved edge
        p.addArc(
            self.xcenter,
            self.ycenter,
            inner_radius,
            90 - ((endangle - cornerangle_delta) * 180 / pi),
            90 - ((startangle + cornerangle_delta) * 180 / pi),
            moveTo=True,
        )
        # Corner edge - straight lines assumes small angle!
        # TODO - Use self._draw_arc_line(p, ...) here if we expose corner setting
        p.lineTo(x0 + shaft_inner_radius * startsin, y0 + shaft_inner_radius * startcos)
        p.lineTo(x0 + shaft_outer_radius * startsin, y0 + shaft_outer_radius * startcos)
        # Outer curved edge
        p.addArc(
            self.xcenter,
            self.ycenter,
            outer_radius,
            90 - ((endangle - cornerangle_delta) * 180 / pi),
            90 - ((startangle + cornerangle_delta) * 180 / pi),
            reverse=True,
        )
        # Corner edges
        p.lineTo(x0 + shaft_outer_radius * endsin, y0 + shaft_outer_radius * endcos)
        p.lineTo(x0 + shaft_inner_radius * endsin, y0 + shaft_inner_radius * endcos)
        p.closePath()
        return p

    def _draw_sigil_arrow(
        self, bottom, center, top, startangle, endangle, strand, **kwargs
    ):
        """Draw ARROW sigil (PRIVATE)."""
        if strand == 1:
            inner_radius = center
            outer_radius = top
            orientation = "right"
        elif strand == -1:
            inner_radius = bottom
            outer_radius = center
            orientation = "left"
        else:
            inner_radius = bottom
            outer_radius = top
            orientation = "right"  # backwards compatibility
        return self._draw_arc_arrow(
            inner_radius,
            outer_radius,
            startangle,
            endangle,
            orientation=orientation,
            **kwargs
        )

    def _draw_sigil_big_arrow(
        self, bottom, center, top, startangle, endangle, strand, **kwargs
    ):
        """Draw BIGARROW sigil, like ARROW but straddles the axis (PRIVATE)."""
        if strand == -1:
            orientation = "left"
        else:
            orientation = "right"
        return self._draw_arc_arrow(
            bottom, top, startangle, endangle, orientation=orientation, **kwargs
        )

    def _draw_arc_arrow(
        self,
        inner_radius,
        outer_radius,
        startangle,
        endangle,
        color,
        border=None,
        shaft_height_ratio=0.4,
        head_length_ratio=0.5,
        orientation="right",
        colour=None,
        **kwargs
    ):
        """Draw an arrow along an arc (PRIVATE)."""
        # Let the UK spelling (colour) override the USA spelling (color)
        if colour is not None:
            color = colour

        strokecolor, color = _stroke_and_fill_colors(color, border)

        # if orientation == 'right':
        #    startangle, endangle = min(startangle, endangle), max(startangle, endangle)
        # elif orientation == 'left':
        #    startangle, endangle = max(startangle, endangle), min(startangle, endangle)
        # else:
        startangle, endangle = min(startangle, endangle), max(startangle, endangle)
        if orientation != "left" and orientation != "right":
            raise ValueError(
                "Invalid orientation %r, should be 'left' or 'right'" % orientation
            )

        angle = float(endangle - startangle)  # angle subtended by arc
        middle_radius = 0.5 * (inner_radius + outer_radius)
        boxheight = outer_radius - inner_radius
        shaft_height = boxheight * shaft_height_ratio
        shaft_inner_radius = middle_radius - 0.5 * shaft_height
        shaft_outer_radius = middle_radius + 0.5 * shaft_height
        headangle_delta = max(
            0.0, min(abs(boxheight) * head_length_ratio / middle_radius, abs(angle))
        )
        if angle < 0:
            headangle_delta *= -1  # reverse it
        if orientation == "right":
            headangle = endangle - headangle_delta
        else:
            headangle = startangle + headangle_delta
        if startangle <= endangle:
            headangle = max(min(headangle, endangle), startangle)
        else:
            headangle = max(min(headangle, startangle), endangle)
        if not (
            startangle <= headangle <= endangle or endangle <= headangle <= startangle
        ):
            raise RuntimeError(
                "Problem drawing arrow, invalid positions. "
                "Start angle: %s, Head angle: %s, "
                "End angle: %s, Angle: %s" % (startangle, headangle, endangle, angle)
            )

        # Calculate trig values for angle and coordinates
        startcos, startsin = cos(startangle), sin(startangle)
        headcos, headsin = cos(headangle), sin(headangle)
        endcos, endsin = cos(endangle), sin(endangle)
        x0, y0 = self.xcenter, self.ycenter  # origin of the circle
        if 0.5 >= abs(angle) and abs(headangle_delta) >= abs(angle):
            # If the angle is small, and the arrow is all head,
            # cheat and just use a triangle.
            if orientation == "right":
                x1, y1 = (x0 + inner_radius * startsin, y0 + inner_radius * startcos)
                x2, y2 = (x0 + outer_radius * startsin, y0 + outer_radius * startcos)
                x3, y3 = (x0 + middle_radius * endsin, y0 + middle_radius * endcos)
            else:
                x1, y1 = (x0 + inner_radius * endsin, y0 + inner_radius * endcos)
                x2, y2 = (x0 + outer_radius * endsin, y0 + outer_radius * endcos)
                x3, y3 = (x0 + middle_radius * startsin, y0 + middle_radius * startcos)
            # return draw_polygon([(x1,y1),(x2,y2),(x3,y3)], color, border,
            #                    stroke_line_join=1)
            return Polygon(
                [x1, y1, x2, y2, x3, y3],
                strokeColor=border or color,
                fillColor=color,
                strokeLineJoin=1,  # 1=round, not mitre!
                strokewidth=0,
            )
        elif orientation == "right":
            p = ArcPath(
                strokeColor=strokecolor,
                fillColor=color,
                # default is mitre/miter which can stick out too much:
                strokeLineJoin=1,  # 1=round
                strokewidth=0,
                **kwargs
            )
            # Note reportlab counts angles anti-clockwise from the horizontal
            # (as in mathematics, e.g. complex numbers and polar coordinates)
            # but we use clockwise from the vertical.  Also reportlab uses
            # degrees, but we use radians.
            p.addArc(
                self.xcenter,
                self.ycenter,
                shaft_inner_radius,
                90 - (headangle * 180 / pi),
                90 - (startangle * 180 / pi),
                moveTo=True,
            )
            p.addArc(
                self.xcenter,
                self.ycenter,
                shaft_outer_radius,
                90 - (headangle * 180 / pi),
                90 - (startangle * 180 / pi),
                reverse=True,
            )
            if abs(angle) < 0.5:
                p.lineTo(x0 + outer_radius * headsin, y0 + outer_radius * headcos)
                p.lineTo(x0 + middle_radius * endsin, y0 + middle_radius * endcos)
                p.lineTo(x0 + inner_radius * headsin, y0 + inner_radius * headcos)
            else:
                self._draw_arc_line(
                    p,
                    outer_radius,
                    middle_radius,
                    90 - (headangle * 180 / pi),
                    90 - (endangle * 180 / pi),
                )
                self._draw_arc_line(
                    p,
                    middle_radius,
                    inner_radius,
                    90 - (endangle * 180 / pi),
                    90 - (headangle * 180 / pi),
                )
            p.closePath()
            return p
        else:
            p = ArcPath(
                strokeColor=strokecolor,
                fillColor=color,
                # default is mitre/miter which can stick out too much:
                strokeLineJoin=1,  # 1=round
                strokewidth=0,
                **kwargs
            )
            # Note reportlab counts angles anti-clockwise from the horizontal
            # (as in mathematics, e.g. complex numbers and polar coordinates)
            # but we use clockwise from the vertical.  Also reportlab uses
            # degrees, but we use radians.
            p.addArc(
                self.xcenter,
                self.ycenter,
                shaft_inner_radius,
                90 - (endangle * 180 / pi),
                90 - (headangle * 180 / pi),
                moveTo=True,
                reverse=True,
            )
            p.addArc(
                self.xcenter,
                self.ycenter,
                shaft_outer_radius,
                90 - (endangle * 180 / pi),
                90 - (headangle * 180 / pi),
                reverse=False,
            )
            # Note - two staight lines is only a good approximation for small
            # head angle, in general will need to curved lines here:
            if abs(angle) < 0.5:
                p.lineTo(x0 + outer_radius * headsin, y0 + outer_radius * headcos)
                p.lineTo(x0 + middle_radius * startsin, y0 + middle_radius * startcos)
                p.lineTo(x0 + inner_radius * headsin, y0 + inner_radius * headcos)
            else:
                self._draw_arc_line(
                    p,
                    outer_radius,
                    middle_radius,
                    90 - (headangle * 180 / pi),
                    90 - (startangle * 180 / pi),
                )
                self._draw_arc_line(
                    p,
                    middle_radius,
                    inner_radius,
                    90 - (startangle * 180 / pi),
                    90 - (headangle * 180 / pi),
                )
            p.closePath()
            return p

    def _draw_sigil_jaggy(
        self,
        bottom,
        center,
        top,
        startangle,
        endangle,
        strand,
        color,
        border=None,
        **kwargs
    ):
        """Draw JAGGY sigil (PRIVATE).

        Although we may in future expose the head/tail jaggy lengths, for now
        both the left and right edges are drawn jagged.
        """
        if strand == 1:
            inner_radius = center
            outer_radius = top
            teeth = 2
        elif strand == -1:
            inner_radius = bottom
            outer_radius = center
            teeth = 2
        else:
            inner_radius = bottom
            outer_radius = top
            teeth = 4

        # TODO, expose these settings?
        tail_length_ratio = 1.0
        head_length_ratio = 1.0

        strokecolor, color = _stroke_and_fill_colors(color, border)

        startangle, endangle = min(startangle, endangle), max(startangle, endangle)
        angle = float(endangle - startangle)  # angle subtended by arc
        height = outer_radius - inner_radius

        assert startangle <= endangle and angle >= 0
        if head_length_ratio and tail_length_ratio:
            headangle = max(
                endangle
                - min(height * head_length_ratio / (center * teeth), angle * 0.5),
                startangle,
            )
            tailangle = min(
                startangle
                + min(height * tail_length_ratio / (center * teeth), angle * 0.5),
                endangle,
            )
            # With very small features, can due to floating point calculations
            # violate the assertion below that start <= tail <= head <= end
            tailangle = min(tailangle, headangle)
        elif head_length_ratio:
            headangle = max(
                endangle - min(height * head_length_ratio / (center * teeth), angle),
                startangle,
            )
            tailangle = startangle
        else:
            headangle = endangle
            tailangle = min(
                startangle + min(height * tail_length_ratio / (center * teeth), angle),
                endangle,
            )

        if not startangle <= tailangle <= headangle <= endangle:
            raise RuntimeError(
                "Problem drawing jaggy sigil, invalid "
                "positions. Start angle: %s, "
                "Tail angle: %s, Head angle: %s, End angle %s, "
                "Angle: %s" % (startangle, tailangle, headangle, endangle, angle)
            )

        # Calculate trig values for angle and coordinates
        startcos, startsin = cos(startangle), sin(startangle)
        headcos, headsin = cos(headangle), sin(headangle)
        endcos, endsin = cos(endangle), sin(endangle)
        x0, y0 = self.xcenter, self.ycenter  # origin of the circle

        p = ArcPath(
            strokeColor=strokecolor,
            fillColor=color,
            # default is mitre/miter which can stick out too much:
            strokeLineJoin=1,  # 1=round
            strokewidth=0,
            **kwargs
        )
        # Note reportlab counts angles anti-clockwise from the horizontal
        # (as in mathematics, e.g. complex numbers and polar coordinates)
        # but we use clockwise from the vertical.  Also reportlab uses
        # degrees, but we use radians.
        p.addArc(
            self.xcenter,
            self.ycenter,
            inner_radius,
            90 - (headangle * 180 / pi),
            90 - (tailangle * 180 / pi),
            moveTo=True,
        )
        for i in range(0, teeth):
            p.addArc(
                self.xcenter,
                self.ycenter,
                inner_radius + i * height / teeth,
                90 - (tailangle * 180 / pi),
                90 - (startangle * 180 / pi),
            )
            # Curved line needed when drawing long jaggies
            self._draw_arc_line(
                p,
                inner_radius + i * height / teeth,
                inner_radius + (i + 1) * height / teeth,
                90 - (startangle * 180 / pi),
                90 - (tailangle * 180 / pi),
            )
        p.addArc(
            self.xcenter,
            self.ycenter,
            outer_radius,
            90 - (headangle * 180 / pi),
            90 - (tailangle * 180 / pi),
            reverse=True,
        )
        for i in range(0, teeth):
            p.addArc(
                self.xcenter,
                self.ycenter,
                outer_radius - i * height / teeth,
                90 - (endangle * 180 / pi),
                90 - (headangle * 180 / pi),
                reverse=True,
            )
            # Curved line needed when drawing long jaggies
            self._draw_arc_line(
                p,
                outer_radius - i * height / teeth,
                outer_radius - (i + 1) * height / teeth,
                90 - (endangle * 180 / pi),
                90 - (headangle * 180 / pi),
            )
        p.closePath()
        return p
