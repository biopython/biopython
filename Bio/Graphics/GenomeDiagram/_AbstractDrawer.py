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

"""AbstractDrawer module (considered to be a private module, the API may change!).

Provides:
 - AbstractDrawer - Superclass for methods common to the Drawer objects
 - page_sizes - Method that returns a ReportLab pagesize when passed
   a valid ISO size
 - draw_box - Method that returns a closed path object when passed
   the proper co-ordinates.  For HORIZONTAL boxes only.
 - angle2trig - Method that returns a tuple of values that are the
   vector for rotating a point through a passed angle,
   about an origin
 - intermediate_points - Method that returns a list of values intermediate
   between the points in a passed dataset

For drawing capabilities, this module uses reportlab to draw and write
the diagram: http://www.reportlab.com

For dealing with biological information, the package expects Biopython objects
like SeqFeatures.
"""

# ReportLab imports

from reportlab.lib import pagesizes
from reportlab.lib import colors
from reportlab.graphics.shapes import Polygon

from math import pi, sin, cos
from itertools import islice

################################################################################
# METHODS
################################################################################


# Utility method to translate strings to ISO page sizes
def page_sizes(size):
    """Convert size string into a Reportlab pagesize.

    Arguments:
     - size - A string representing a standard page size, eg 'A4' or 'LETTER'

    """
    sizes = {  # ReportLab pagesizes, keyed by ISO string
        "A0": pagesizes.A0,
        "A1": pagesizes.A1,
        "A2": pagesizes.A2,
        "A3": pagesizes.A3,
        "A4": pagesizes.A4,
        "A5": pagesizes.A5,
        "A6": pagesizes.A6,
        "B0": pagesizes.B0,
        "B1": pagesizes.B1,
        "B2": pagesizes.B2,
        "B3": pagesizes.B3,
        "B4": pagesizes.B4,
        "B5": pagesizes.B5,
        "B6": pagesizes.B6,
        "ELEVENSEVENTEEN": pagesizes.ELEVENSEVENTEEN,
        "LEGAL": pagesizes.LEGAL,
        "LETTER": pagesizes.LETTER,
    }
    try:
        return sizes[size]
    except KeyError:
        raise ValueError("%s not in list of page sizes" % size) from None


def _stroke_and_fill_colors(color, border):
    """Deal with  border and fill colors (PRIVATE)."""
    if not isinstance(color, colors.Color):
        raise ValueError("Invalid color %r" % color)

    if color == colors.white and border is None:
        # Force black border on white boxes with undefined border
        strokecolor = colors.black
    elif border is None:
        strokecolor = color  # use fill color
    elif border:
        if not isinstance(border, colors.Color):
            raise ValueError("Invalid border color %r" % border)
        strokecolor = border
    else:
        # e.g. False
        strokecolor = None

    return strokecolor, color


def draw_box(
    point1, point2, color=colors.lightgreen, border=None, colour=None, **kwargs
):
    """Draw a box.

    Arguments:
     - point1, point2 - coordinates for opposite corners of the box
       (x,y tuples)
     - color /colour - The color for the box (colour takes priority
       over color)
     - border - Border color for the box

    Returns a closed path object, beginning at (x1,y1) going round
    the four points in order, and filling with the passed color.
    """
    x1, y1 = point1
    x2, y2 = point2

    # Let the UK spelling (colour) override the USA spelling (color)
    if colour is not None:
        color = colour
        del colour

    strokecolor, color = _stroke_and_fill_colors(color, border)

    x1, y1, x2, y2 = min(x1, x2), min(y1, y2), max(x1, x2), max(y1, y2)
    return Polygon(
        [x1, y1, x2, y1, x2, y2, x1, y2],
        strokeColor=strokecolor,
        fillColor=color,
        strokewidth=0,
        **kwargs
    )


def draw_cut_corner_box(
    point1, point2, corner=0.5, color=colors.lightgreen, border=None, **kwargs
):
    """Draw a box with the corners cut off."""
    x1, y1 = point1
    x2, y2 = point2

    if not corner:
        return draw_box(point1, point2, color, border)
    elif corner < 0:
        raise ValueError("Arrow head length ratio should be positive")

    strokecolor, color = _stroke_and_fill_colors(color, border)

    boxheight = y2 - y1
    boxwidth = x2 - x1
    x_corner = min(boxheight * 0.5 * corner, boxwidth * 0.5)
    y_corner = min(boxheight * 0.5 * corner, boxheight * 0.5)

    points = [
        x1,
        y1 + y_corner,
        x1,
        y2 - y_corner,
        x1 + x_corner,
        y2,
        x2 - x_corner,
        y2,
        x2,
        y2 - y_corner,
        x2,
        y1 + y_corner,
        x2 - x_corner,
        y1,
        x1 + x_corner,
        y1,
    ]
    return Polygon(
        deduplicate(points),
        strokeColor=strokecolor,
        strokeWidth=1,
        strokeLineJoin=1,  # 1=round
        fillColor=color,
        **kwargs
    )


def draw_polygon(
    list_of_points, color=colors.lightgreen, border=None, colour=None, **kwargs
):
    """Draw polygon.

    Arguments:
     - list_of_point - list of (x,y) tuples for the corner coordinates
     - color / colour - The color for the box

    Returns a closed path object, beginning at (x1,y1) going round
    the four points in order, and filling with the passed colour.

    """
    # Let the UK spelling (colour) override the USA spelling (color)
    if colour is not None:
        color = colour
        del colour

    strokecolor, color = _stroke_and_fill_colors(color, border)

    xy_list = []
    for (x, y) in list_of_points:
        xy_list.append(x)
        xy_list.append(y)

    return Polygon(
        deduplicate(xy_list),
        strokeColor=strokecolor,
        fillColor=color,
        strokewidth=0,
        **kwargs
    )


def draw_arrow(
    point1,
    point2,
    color=colors.lightgreen,
    border=None,
    shaft_height_ratio=0.4,
    head_length_ratio=0.5,
    orientation="right",
    colour=None,
    **kwargs
):
    """Draw an arrow.

    Returns a closed path object representing an arrow enclosed by the
    box with corners at {point1=(x1,y1), point2=(x2,y2)}, a shaft height
    given by shaft_height_ratio (relative to box height), a head length
    given by head_length_ratio (also relative to box height), and
    an orientation that may be 'left' or 'right'.
    """
    x1, y1 = point1
    x2, y2 = point2

    if shaft_height_ratio < 0 or 1 < shaft_height_ratio:
        raise ValueError("Arrow shaft height ratio should be in range 0 to 1")
    if head_length_ratio < 0:
        raise ValueError("Arrow head length ratio should be positive")

    # Let the UK spelling (colour) override the USA spelling (color)
    if colour is not None:
        color = colour
        del colour

    strokecolor, color = _stroke_and_fill_colors(color, border)

    # Depending on the orientation, we define the bottom left (x1, y1) and
    # top right (x2, y2) coordinates differently, but still draw the box
    # using the same relative co-ordinates:
    xmin, ymin = min(x1, x2), min(y1, y2)
    xmax, ymax = max(x1, x2), max(y1, y2)
    if orientation == "right":
        x1, x2, y1, y2 = xmin, xmax, ymin, ymax
    elif orientation == "left":
        x1, x2, y1, y2 = xmax, xmin, ymin, ymax
    else:
        raise ValueError(
            "Invalid orientation %r, should be 'left' or 'right'" % orientation
        )

    # We define boxheight and boxwidth accordingly, and calculate the shaft
    # height from these.  We also ensure that the maximum head length is
    # the width of the box enclosure
    boxheight = y2 - y1
    boxwidth = x2 - x1
    shaftheight = boxheight * shaft_height_ratio
    headlength = min(abs(boxheight) * head_length_ratio, abs(boxwidth))
    if boxwidth < 0:
        headlength *= -1  # reverse it

    shafttop = 0.5 * (boxheight + shaftheight)
    shaftbase = boxheight - shafttop
    headbase = boxwidth - headlength
    midheight = 0.5 * boxheight

    points = [
        x1,
        y1 + shafttop,
        x1 + headbase,
        y1 + shafttop,
        x1 + headbase,
        y2,
        x2,
        y1 + midheight,
        x1 + headbase,
        y1,
        x1 + headbase,
        y1 + shaftbase,
        x1,
        y1 + shaftbase,
    ]

    return Polygon(
        deduplicate(points),
        strokeColor=strokecolor,
        # strokeWidth=max(1, int(boxheight/40.)),
        strokeWidth=1,
        # default is mitre/miter which can stick out too much:
        strokeLineJoin=1,  # 1=round
        fillColor=color,
        **kwargs
    )


def deduplicate(points):
    """Remove adjacent duplicate points.

    This is important for use with the Polygon class since reportlab has a
    bug with duplicate points.

    Arguments:
     - points - list of points [x1, y1, x2, y2,...]

    Returns a list in the same format with consecutive duplicates removed
    """
    assert len(points) % 2 == 0
    if len(points) < 2:
        return points
    newpoints = points[0:2]
    for x, y in zip(islice(points, 2, None, 2), islice(points, 3, None, 2)):
        if x != newpoints[-2] or y != newpoints[-1]:
            newpoints.append(x)
            newpoints.append(y)
    return newpoints


def angle2trig(theta):
    """Convert angle to a reportlab ready tuple.

    Arguments:
     - theta -  Angle in degrees, counter clockwise from horizontal

    Returns a representation of the passed angle in a format suitable
    for ReportLab rotations (i.e. cos(theta), sin(theta), -sin(theta),
    cos(theta) tuple)
    """
    c = cos(theta * pi / 180)
    s = sin(theta * pi / 180)
    return (c, s, -s, c)  # Vector for rotating point around an origin


def intermediate_points(start, end, graph_data):
    """Generate intermediate points describing provided graph data..

    Returns a list of (start, end, value) tuples describing the passed
    graph data as 'bins' between position midpoints.
    """
    newdata = []  # data in form (X0, X1, val)
    # add first block
    newdata.append(
        (
            start,
            graph_data[0][0] + (graph_data[1][0] - graph_data[0][0]) / 2.0,
            graph_data[0][1],
        )
    )
    # add middle set
    for index in range(1, len(graph_data) - 1):
        lastxval, lastyval = graph_data[index - 1]
        xval, yval = graph_data[index]
        nextxval, nextyval = graph_data[index + 1]
        newdata.append(
            (lastxval + (xval - lastxval) / 2.0, xval + (nextxval - xval) / 2.0, yval)
        )
    # add last block
    newdata.append((xval + (nextxval - xval) / 2.0, end, graph_data[-1][1]))
    return newdata


################################################################################
# CLASSES
################################################################################


class AbstractDrawer:
    """Abstract Drawer.

    Attributes:
     - tracklines    Boolean for whether to draw lines delineating tracks
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
     - cross_track_links List of tuples each with four entries (track A,
       feature A, track B, feature B) to be linked.

    """

    def __init__(
        self,
        parent,
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
        cross_track_links=None,
    ):
        """Create the object.

        Arguments:
         - parent    Diagram object containing the data that the drawer draws
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
         - cross_track_links List of tuples each with four entries (track A,
           feature A, track B, feature B) to be linked.

        """
        self._parent = parent  # The calling Diagram object

        # Perform 'administrative' tasks of setting up the page
        self.set_page_size(pagesize, orientation)  # Set drawing size
        self.set_margins(x, y, xl, xr, yt, yb)  # Set page margins
        self.set_bounds(start, end)  # Set limits on what will be drawn
        self.tracklines = tracklines  # Set flags
        if cross_track_links is None:
            cross_track_links = []
        else:
            self.cross_track_links = cross_track_links

    def set_page_size(self, pagesize, orientation):
        """Set page size of the drawing..

        Arguments:
         - pagesize      Size of the output image, a tuple of pixels (width,
           height, or a string in the reportlab.lib.pagesizes
           set of ISO sizes.
         - orientation   String: 'landscape' or 'portrait'

        """
        if isinstance(pagesize, str):  # A string, so translate
            pagesize = page_sizes(pagesize)
        elif isinstance(pagesize, tuple):  # A tuple, so don't translate
            pagesize = pagesize
        else:
            raise ValueError("Page size %s not recognised" % pagesize)
        shortside, longside = min(pagesize), max(pagesize)

        orientation = orientation.lower()
        if orientation not in ("landscape", "portrait"):
            raise ValueError("Orientation %s not recognised" % orientation)
        if orientation == "landscape":
            self.pagesize = (longside, shortside)
        else:
            self.pagesize = (shortside, longside)

    def set_margins(self, x, y, xl, xr, yt, yb):
        """Set page margins.

        Arguments:
         - x         Float(0->1), Absolute X margin as % of page
         - y         Float(0->1), Absolute Y margin as % of page
         - xl        Float(0->1), Left X margin as % of page
         - xr        Float(0->1), Right X margin as % of page
         - yt        Float(0->1), Top Y margin as % of page
         - yb        Float(0->1), Bottom Y margin as % of page

        Set the page margins as proportions of the page 0->1, and also
        set the page limits x0, y0 and xlim, ylim, and page center
        xorigin, yorigin, as well as overall page width and height
        """
        # Set left, right, top and bottom margins
        xmargin_l = xl or x
        xmargin_r = xr or x
        ymargin_top = yt or y
        ymargin_btm = yb or y

        # Set page limits, center and height/width
        self.x0, self.y0 = self.pagesize[0] * xmargin_l, self.pagesize[1] * ymargin_btm
        self.xlim, self.ylim = (
            self.pagesize[0] * (1 - xmargin_r),
            self.pagesize[1] * (1 - ymargin_top),
        )
        self.pagewidth = self.xlim - self.x0
        self.pageheight = self.ylim - self.y0
        self.xcenter, self.ycenter = (
            self.x0 + self.pagewidth / 2.0,
            self.y0 + self.pageheight / 2.0,
        )

    def set_bounds(self, start, end):
        """Set start and end points for the drawing as a whole.

        Arguments:
         - start - The first base (or feature mark) to draw from
         - end - The last base (or feature mark) to draw to

        """
        low, high = self._parent.range()  # Extent of tracks

        if start is not None and end is not None and start > end:
            start, end = end, start

        if start is None or start < 0:  # Check validity of passed args and
            start = 0  # default to 0
        if end is None or end < 0:
            end = high + 1  # default to track range top limit

        self.start, self.end = int(start), int(end)
        self.length = self.end - self.start + 1

    def is_in_bounds(self, value):
        """Check if given value is within the region selected for drawing.

        Arguments:
         - value - A base position

        """
        if value >= self.start and value <= self.end:
            return 1
        return 0

    def __len__(self):
        """Return the length of the region to be drawn."""
        return self.length

    def _current_track_start_end(self):
        track = self._parent[self.current_track_level]
        if track.start is None:
            start = self.start
        else:
            start = max(self.start, track.start)
        if track.end is None:
            end = self.end
        else:
            end = min(self.end, track.end)
        return start, end
