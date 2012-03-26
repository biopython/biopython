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

""" AbstractDrawer module (considered to be a private module, the API may change!)

    Provides:

    o AbstractDrawer -    Superclass for methods common to *Drawer objects

    o page_sizes -          Method that returns a ReportLab pagesize when passed
                            a valid ISO size

    o draw_box -            Method that returns a closed path object when passed
                            the proper co-ordinates.  For HORIZONTAL boxes only.

    o angle2trig -          Method that returns a tuple of values that are the
                            vector for rotating a point through a passed angle,
                            about an origin

    o intermediate_points - Method that returns a list of values intermediate
                            between the points in a passed dataset
    
    For drawing capabilities, this module uses reportlab to draw and write
    the diagram:

    http://www.reportlab.com

    For dealing with biological information, the package expects BioPython
    objects:

    http://www.biopython.org
"""

# ReportLab imports
from reportlab.lib import pagesizes
from reportlab.lib import colors
from reportlab.graphics.shapes import *

from math import pi

################################################################################
# METHODS
################################################################################
# Utility method to translate strings to ISO page sizes
def page_sizes(size):
    """ page_sizes(size)

        o size        A string representing a standard page size

        Returns a ReportLab pagesize when passed a valid size string
    """
    sizes = {'A0': pagesizes.A0,    # ReportLab pagesizes, keyed by ISO string
             'A1': pagesizes.A1,
             'A2': pagesizes.A2,
             'A3': pagesizes.A3,
             'A4': pagesizes.A4,
             'A5': pagesizes.A5,
             'A6': pagesizes.A6,
             'B0': pagesizes.B0,
             'B1': pagesizes.B1,
             'B2': pagesizes.B2,
             'B3': pagesizes.B3,
             'B4': pagesizes.B4,
             'B5': pagesizes.B5,
             'B6': pagesizes.B6,
             'ELEVENSEVENTEEN': pagesizes.ELEVENSEVENTEEN,
             'LEGAL': pagesizes.LEGAL,
             'LETTER': pagesizes.LETTER
             }
    try:
        return sizes[size]
    except:
        raise ValueError, "%s not in list of page sizes" % size


def draw_box(point1, point2,
             color=colors.lightgreen, border=None, colour=None,
             **kwargs):
    """ draw_box(self, (x1, y1), (x2, y2), (x3, y3), (x4, y4),
              color=colors.lightgreen)

        o point1, point2 Co-ordinates for opposite corners of the box
                         (x,y tuples)
        
        o color /colour       The color for the box
                              (colour takes priority over color)
                              
        o border              Border color for the box

        Returns a closed path object, beginning at (x1,y1) going round
        the four points in order, and filling with the passed color.            
    """
    x1, y1 = point1
    x2, y2 = point2

    #Let the UK spelling (colour) override the USA spelling (color)
    if colour is not None:
        color = colour
        del colour

    if not isinstance(color, colors.Color):
        raise ValueError("Invalid color %s" % repr(color))
    
    if color == colors.white and border is None:   # Force black border on 
        strokecolor = colors.black                 # white boxes with
    elif border is None:                           # undefined border, else
        strokecolor = color                        # use fill color
    elif border:
        if not isinstance(border, colors.Color):
            raise ValueError("Invalid border color %s" % repr(border))
        strokecolor = border
    else:
        #e.g. False
        strokecolor = None

    x1, y1, x2, y2 = min(x1, x2), min(y1, y2), max(x1, x2), max(y1, y2)
    return Polygon([x1, y1, x2, y1, x2, y2, x1, y2],
                   strokeColor=strokecolor,
                   fillColor=color,
                   strokewidth=0,
                   **kwargs)


def draw_polygon(list_of_points,
                 color=colors.lightgreen, border=None, colour=None,
                 **kwargs):
    """ draw_polygon(self, (x1, y1), (x2, y2), (x3, y3), (x4, y4)
              colour=colors.lightgreen)

        o list_of_point = list of (x,y) tuples for the corner coordinates
        
        o colour              The colour for the box

        Returns a closed path object, beginning at (x1,y1) going round
        the four points in order, and filling with the passed colour.          
    """
    #Let the UK spelling (colour) override the USA spelling (color)
    if colour is not None:
        color = colour
        del colour

    if color == colors.white and border is None:   # Force black border on 
        strokecolor = colors.black                 # white boxes with
    elif border is None:                           # undefined border, else
        strokecolor = color                        # use fill colour
    elif border:
        strokecolor = border
    else:
        #e.g. False
        strokecolor = None

    xy_list = []
    for (x,y) in list_of_points:
        xy_list.append(x)
        xy_list.append(y)

    return Polygon(xy_list,
                   strokeColor=strokecolor,
                   fillColor=color,
                   strokewidth=0,
                   **kwargs)


def draw_arrow(point1, point2, color=colors.lightgreen, border=None,
               shaft_height_ratio=0.4, head_length_ratio=0.5, orientation='right',
               colour=None, **kwargs):
    """ Returns a closed path object representing an arrow enclosed by the
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

    #Let the UK spelling (colour) override the USA spelling (color)
    if colour is not None:
        color = colour
        del colour

    if color == colors.white and border is None:   # Force black border on 
        strokecolor = colors.black                 # white boxes with
    elif border is None:                           # undefined border, else
        strokecolor = color                        # use fill colour
    elif border:
        strokecolor = border
    else:
        #e.g. False
        strokecolor = None

    # Depending on the orientation, we define the bottom left (x1, y1) and
    # top right (x2, y2) coordinates differently, but still draw the box
    # using the same relative co-ordinates:
    xmin, ymin = min(x1, x2), min(y1, y2)
    xmax, ymax = max(x1, x2), max(y1, y2)
    if orientation == 'right':
        x1, x2, y1, y2 = xmin, xmax, ymin, ymax
    elif orientation == 'left':
        x1, x2, y1, y2 = xmax, xmin, ymin, ymax
    else:
        raise ValueError("Invalid orientation %s, should be 'left' or 'right'" \
                         % repr(orientation))

    # We define boxheight and boxwidth accordingly, and calculate the shaft
    # height from these.  We also ensure that the maximum head length is
    # the width of the box enclosure
    boxheight = y2-y1
    boxwidth = x2-x1
    shaftheight = boxheight*shaft_height_ratio
    headlength = min(abs(boxheight)*head_length_ratio, abs(boxwidth))
    if boxwidth < 0:
        headlength *= -1 #reverse it


    shafttop = 0.5*(boxheight+shaftheight)
    shaftbase = boxheight-shafttop
    headbase = boxwidth-headlength
    midheight = 0.5*boxheight
    return Polygon([x1, y1+shafttop,
                    x1+headbase, y1+shafttop,
                    x1+headbase, y2,
                    x2, y1+midheight,
                    x1+headbase, y1,
                    x1+headbase, y1+shaftbase,
                    x1, y1+shaftbase],
                   strokeColor=strokecolor,
                   #strokeWidth=max(1, int(boxheight/40.)),
                   strokeWidth=1,
                   #default is mitre/miter which can stick out too much:
                   strokeLineJoin=1, #1=round
                   fillColor=color,
                   **kwargs)

def angle2trig(theta):
    """ angle2trig(angle)

        o theta     Angle in degrees, counter clockwise from horizontal

        Returns a representation of the passed angle in a format suitable
        for ReportLab rotations (i.e. cos(theta), sin(theta), -sin(theta),
        cos(theta) tuple)
    """
    c = cos(theta * pi / 180)
    s = sin(theta * pi / 180)
    return(c, s, -s, c)         # Vector for rotating point around an origin


def intermediate_points(start, end, graph_data):
    """ intermediate_points(start, end, graph_data)

        o graph_data

        o start

        o end

        Returns a list of (start, end, value) tuples describing the passed
        graph data as 'bins' between position midpoints.
    """
    #print start, end, len(graph_data)
    newdata = []    # data in form (X0, X1, val)
    # add first block
    newdata.append((start, graph_data[0][0]+(graph_data[1][0]-graph_data[0][0])/2.,
                    graph_data[0][1]))
    # add middle set
    for index in xrange(1, len(graph_data)-1):
        lastxval, lastyval = graph_data[index-1]
        xval, yval = graph_data[index]
        nextxval, nextyval = graph_data[index+1]
        newdata.append( (lastxval+(xval-lastxval)/2.,
                         xval+(nextxval-xval)/2., yval) )
    # add last block
    newdata.append( (xval+(nextxval-xval)/2.,
                         end, graph_data[-1][1]) )
    #print newdata[-1]
    #print newdata
    return newdata

################################################################################
# CLASSES
################################################################################

class AbstractDrawer(object):
    """ AbstractDrawer

        Provides:

        Methods:

        o __init__(self, parent, pagesize='A3', orientation='landscape',
                 x=0.05, y=0.05, xl=None, xr=None, yt=None, yb=None,
                 start=None, end=None, tracklines=0) Called on instantiation

        o set_page_size(self, pagesize, orientation)    Set the page size to the
                                                    passed size and orientation

        o set_margins(self, x, y, xl, xr, yt, yb)   Set the drawable area of the
                                                    page

        o set_bounds(self, start, end)  Set the bounds for the elements to be
                                        drawn

        o is_in_bounds(self, value)     Returns a boolean for whether the position
                                        is actually to be drawn

        o __len__(self)     Returns the length of sequence that will be drawn

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

        o length        Size of sequence to be drawn
        
        o cross_track_links List of tuples each with four entries (track A,
                            feature A, track B, feature B) to be linked.
    """
    def __init__(self, parent, pagesize='A3', orientation='landscape',
                 x=0.05, y=0.05, xl=None, xr=None, yt=None, yb=None,
                 start=None, end=None, tracklines=0, cross_track_links=None):
        """ __init__(self, parent, pagesize='A3', orientation='landscape',
                 x=0.05, y=0.05, xl=None, xr=None, yt=None, yb=None,
                 start=None, end=None, tracklines=0)

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

            o cross_track_links List of tuples each with four entries (track A,
                                feature A, track B, feature B) to be linked.
        """
        self._parent = parent   # The calling Diagram object

        # Perform 'administrative' tasks of setting up the page
        self.set_page_size(pagesize, orientation)   # Set drawing size
        self.set_margins(x, y, xl, xr, yt, yb)      # Set page margins
        self.set_bounds(start, end) # Set limits on what will be drawn
        self.tracklines = tracklines    # Set flags
        if cross_track_links is None:
            cross_track_links = []
        else:
            self.cross_track_links = cross_track_links
        
    def _set_xcentre(self, value):
        import warnings
        import Bio
        warnings.warn("The _set_xcentre method and .xcentre attribute are deprecated; please use the .xcenter attribute instead", Bio.BiopythonDeprecationWarning)
        self.xcenter = value
    xcentre = property(fget = lambda self : self.xcenter,
                       fset = _set_xcentre,
                       doc="Backwards compatible alias for xcenter (DEPRECATED)")

    def _set_ycentre(self, value):
        import warnings
        import Bio
        warnings.warn("The _set_ycentre method and .xcentre attribute are deprecated; please use the .ycenter attribute instead", Bio.BiopythonDeprecationWarning)
        self.ycenter = value
    ycentre = property(fget = lambda self : self.ycenter,
                       fset = _set_ycentre,
                       doc="Backwards compatible alias for ycenter (DEPRECATED)")

    def set_page_size(self, pagesize, orientation):
        """ set_page_size(self, pagesize, orientation)

            o pagesize      Size of the output image, a tuple of pixels (width,
                            height, or a string in the reportlab.lib.pagesizes
                            set of ISO sizes.

            o orientation   String: 'landscape' or 'portrait'

            Set the size of the drawing
        """
        if type(pagesize) == type('a'):     # A string, so translate
            pagesize = page_sizes(pagesize)
        elif type(pagesize) == type((1,2)): # A tuple, so don't translate
            pagesize = pagesize
        else:
            raise ValueError, "Page size %s not recognised" % pagesize        
        shortside, longside = min(pagesize), max(pagesize)

        orientation = orientation.lower()
        if orientation not in ('landscape', 'portrait'):
            raise ValueError, "Orientation %s not recognised" % orientation
        if orientation == 'landscape':
            self.pagesize = (longside, shortside)
        else:
            self.pagesize = (shortside, longside)


    def set_margins(self, x, y, xl, xr, yt, yb):
        """ set_margins(self, x, y, xl, xr, yt, yb)

            o x         Float(0->1), Absolute X margin as % of page

            o y         Float(0->1), Absolute Y margin as % of page

            o xl        Float(0->1), Left X margin as % of page

            o xr        Float(0->1), Right X margin as % of page

            o yt        Float(0->1), Top Y margin as % of page

            o yb        Float(0->1), Bottom Y margin as % of page

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
        self.x0, self.y0 = self.pagesize[0]*xmargin_l, self.pagesize[1]*ymargin_btm
        self.xlim, self.ylim = self.pagesize[0]*(1-xmargin_r), self.pagesize[1]*(1-ymargin_top)
        self.pagewidth = self.xlim-self.x0
        self.pageheight = self.ylim-self.y0
        self.xcenter, self.ycenter = self.x0+self.pagewidth/2., self.y0+self.pageheight/2.

            
    def set_bounds(self, start, end):
        """ set_bounds(self, start, end)

            o start     The first base (or feature mark) to draw from

            o end       The last base (or feature mark) to draw to

            Sets start and end points for the drawing as a whole
        """
        low, high = self._parent.range()  # Extent of tracks

        if start is not None and end is not None and start > end:
            start, end = end, start

        if start is None or start < 0:  # Check validity of passed args and 
            start = 0   # default to 0
        if end is None or end < 0:
            end = high + 1  # default to track range top limit
        
        self.start, self.end = int(start), int(end)
        self.length = self.end - self.start + 1


    def is_in_bounds(self, value):
        """ is_in_bounds(self, value)

            o value   A base position

            Returns 1 if the value is within the region selected for drawing
        """
        if value >= self.start and value <= self.end:
            return 1
        return 0


    def __len__(self):
        """ __len__(self)

            Returns the length of the region to be drawn
        """
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
