# Copyright 2003-2008 by Leighton Pritchard.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
# Contact:       Leighton Pritchard, Scottish Crop Research Institute,
#                Invergowrie, Dundee, Scotland, DD2 5DA, UK
#                L.Pritchard@scri.ac.uk
################################################################################

""" Diagram module

    Provides:

    o Diagram -   Container for information concerning the tracks to be
                    drawn in a diagram, and the interface for defining the
                    diagram (possibly split these functions in later version?)

    For drawing capabilities, this module uses reportlab to draw and write
    the diagram:

    http://www.reportlab.com

    For dealing with biological information, the package expects BioPython
    objects - namely SeqRecord ojbects containing SeqFeature objects.
"""

#------------------------------------------------------------------------------
# IMPORTS

# ReportLab
from reportlab.graphics import renderPS, renderPDF, renderSVG
try:
    from reportlab.graphics import renderPM
except ImportError:
    #This is an optional part of ReportLab, so may not be installed.
    renderPM=None
from reportlab.lib import pagesizes

# GenomeDiagram
from _LinearDrawer import LinearDrawer
from _CircularDrawer import CircularDrawer
from _Track import Track

# Builtins
import sys

from Bio.Graphics import _write

#------------------------------------------------------------------------------
# CLASSES

#------------------------------------------------------------
# Diagram

class Diagram(object):
    """ Diagram

        Provides:

        Attributes:

        o name         String, identifier for the diagram

        o tracks       List of Track objects comprising the diagram 

        o format       String, format of the diagram (circular/linear)

        o pagesize     String, the pagesize of output

        o orientation  String, the page orientation (landscape/portrait)

        o x            Float, the proportion of the page to take up with even 
                              X margins

        o y            Float, the proportion of the page to take up with even 
                              Y margins

        o xl           Float, the proportion of the page to take up with the 
                              left X margin

        o xr           Float, the proportion of the page to take up with the 
                              right X margin

        o yt           Float, the proportion of the page to take up with the 
                              top Y margin

        o yb           Float, the proportion of the page to take up with the 
                              bottom Y margin

        o circle_core  Float, the proportion of the available radius to leave
                       empty at the center of a circular diagram (0 to 1).

        o start        Int, the base/aa position to start the diagram at

        o end          Int, the base/aa position to end the diagram at

        o tracklines   Boolean, True if track guidelines are to be drawn

        o fragments    Int, for a linear diagram, the number of equal divisions
                                into which the sequence is divided

        o fragment_size Float, the proportion of the space available to each 
                                   fragment that should be used in drawing

        o track_size   Float, the proportion of the space available to each 
                                  track that should be used in drawing

        o circular     Boolean, True if the genome/sequence to be drawn is, in 
                                reality, circular.  

        Methods:

        o __init__(self, name=None) Called on instantiation

        o draw(self, format='circular', ...) Instructs the package to draw
            the diagram

        o write(self, filename='test1.ps', output='PS') Writes the drawn
            diagram to a specified file, in a specified format.

        o add_track(self, track, track_level) Adds a Track object to the
            diagram, with instructions to place it at a particular level on
            the diagram

        o del_track(self, track_level) Removes the track that is to be drawn
            at a particular level on the diagram

        o get_tracks(self) Returns the list of Track objects to be drawn
            contained in the diagram

        o renumber_tracks(self, low=1) Renumbers all tracks consecutively,
            optionally from a passed lowest number

        o get_levels(self) Returns a list of levels currently occupied by
            Track objects

        o get_drawn_levels(self) Returns a list of levels currently occupied
            by Track objects that will be shown in the drawn diagram (i.e.
            are not hidden)

        o range(self) Returns the lowest- and highest-numbered positions
            contained within features in all tracks on the diagram as a tuple.

        o __getitem__(self, key) Returns the track contained at the level of
            the passed key

        o __str__(self) Returns a formatted string describing the diagram

    """
    def __init__(self, name=None, format='circular', pagesize='A3', 
         orientation='landscape', x=0.05, y=0.05, xl=None, 
         xr=None, yt=None, yb=None, start=None, end=None, 
         tracklines=False, fragments=10, fragment_size=0.9, 
         track_size=0.75, circular=True, circle_core=0.0):
        """ __init__(self, name=None)

            o name  String describing the diagram

            o format    String: 'circular' or 'linear', depending on the sort of
                        diagram required

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

            o tracklines    Boolean flag to show (or not) lines delineating 
                        tracks on the diagram

            o fragments Int, for linear diagrams, the number of sections into
                        which to break the sequence being drawn

            o fragment_size     Float (0->1), for linear diagrams, describing 
                                the proportion of space in a fragment to take
                                up with tracks

            o track_size        Float (0->1) describing the proportion of space
                                in a track to take up with sigils

            o circular  Boolean flag to indicate whether the sequence being
                        drawn is circular
                        

        """
        self.tracks = {}   # Holds all Track objects, keyed by level
        self.name = name    # Description of the diagram
        # Diagram page setup attributes
        self.format = format
        self.pagesize = pagesize
        self.orientation = orientation
        self.x = x
        self.y = y
        self.xl = xl
        self.xr = xr
        self.yt = yt
        self.yb = yb
        self.start = start
        self.end = end
        self.tracklines = tracklines
        self.fragments = fragments
        self.fragment_size = fragment_size
        self.track_size = track_size
        self.circular = circular
        self.circle_core = circle_core
        self.cross_track_links = []

    def set_all_tracks(self, attr, value):
        """ set_all_tracks(self, attr, value)

            o attr      An attribute of the Track class

            o value     The value to set that attribute

            Set the passed attribute of all tracks in the set to the
            passed value
        """
        for track in self.tracks.values():
            if hasattr(track, attr):          # If the feature has the attribute
                if getattr(track, attr) != value:
                    setattr(track, attr, value)   # set it to the passed value     

    def draw(self, format=None, pagesize=None, orientation=None,
             x=None, y=None, xl=None, xr=None, yt=None, yb=None,
             start=None, end=None, tracklines=None, fragments=None,
             fragment_size=None, track_size=None, circular=None,
             circle_core=None, cross_track_links=None):
        """Draw the diagram, with passed parameters overriding existing attributes.
        """
        # Pass the parameters to the drawer objects that will build the 
        # diagrams.  At the moment, we detect overrides with an or in the 
        # Instantiation arguments, but I suspect there's a neater way to do 
        # this.
        if format == 'linear':
            drawer = LinearDrawer(self, pagesize or self.pagesize, 
                                  orientation or self.orientation, 
                                  x or self.x, y or self.y, xl or self.xl, 
                                  xr or self.xr, yt or self.yt, 
                                  yb or self.yb, start or self.start, 
                                  end or self.end, 
                                  tracklines or self.tracklines,
                                  fragments or self.fragments, 
                                  fragment_size or self.fragment_size, 
                                  track_size or self.track_size,
                                  cross_track_links or self.cross_track_links)
        else:
            drawer = CircularDrawer(self, pagesize or self.pagesize, 
                                    orientation or self.orientation, 
                                    x or self.x, y or self.y, xl or self.xl, 
                                    xr or self.xr, yt or self.yt, 
                                    yb or self.yb, start or self.start, 
                                    end or self.end, 
                                    tracklines or self.tracklines,
                                    track_size or self.track_size,
                                    circular or self.circular,
                                    circle_core or self.circle_core,
                                    cross_track_links or self.cross_track_links)
        drawer.draw()   # Tell the drawer to complete the drawing
        self.drawing = drawer.drawing  # Get the completed drawing
        
    def write(self, filename='test1.ps', output='PS', dpi=72):
        """ write(self, filename='test1.ps', output='PS', dpi=72)

            o filename      String indicating the name of the output file,
                            or a handle to write to.

            o output        String indicating output format, one of PS, PDF,
                            SVG, or provided the ReportLab renderPM module is
                            installed, one of the bitmap formats JPG, BMP,
                            GIF, PNG, TIFF or TIFF.  The format can be given
                            in upper or lower case.

            o dpi           Resolution (dots per inch) for bitmap formats.

            Write the completed drawing out to a file in a prescribed format

            No return value.
        """
        return _write(self.drawing, filename, output, dpi=dpi)
        
    def write_to_string(self, output='PS', dpi=72):
        """ write(self, output='PS')

            o output        String indicating output format, one of PS, PDF,
                            SVG, JPG, BMP, GIF, PNG, TIFF or TIFF (as
                            specified for the write method).

            o dpi           Resolution (dots per inch) for bitmap formats.

            Return the completed drawing as a string in a prescribed format
        """
        #The ReportLab drawToString method, which this function used to call,
        #just uses a cStringIO or StringIO handle with the drawToFile method.
        #In order to put all our complicated file format specific code in one
        #place we'll just use a StringIO handle here:
        from StringIO import StringIO
        handle = StringIO()
        self.write(handle, output, dpi)
        return handle.getvalue()

    def add_track(self, track, track_level):
        """ add_track(self, track, track_level)

            o track         Track object to draw

            o track_level   Int, the level at which the track will be drawn
                            (above an arbitrary baseline)

            Add a pre-existing Track to the diagram at a given level
        """
        if track is None:
            raise ValueError("Must specify track")
        if track_level not in self.tracks:     # No track at that level
            self.tracks[track_level] = track   # so just add it
        else:       # Already a track there, so shunt all higher tracks up one
            occupied_levels = self.get_levels() # Get list of occupied levels...
            occupied_levels.sort()              # ...sort it...
            occupied_levels.reverse()           # ...reverse it (highest first)
            for val in occupied_levels:
                # If track value >= that to be added
                if val >= track.track_level:
                    self.tracks[val+1] = self.tracks[val] # ...increment by 1
            self.tracks[track_level] = track   # And put the new track in
        self.tracks[track_level].track_level = track_level
                

    def new_track(self, track_level, **args):
        """ new_track(self, track_level) -> Track

            o track_level   Int, the level at which the track will be drawn
                            (above an arbitrary baseline)

            Add a new Track to the diagram at a given level and returns it for
            further user manipulation.
        """
        newtrack = Track()
        for key in args:
            setattr(newtrack, key, args[key])
        if track_level not in self.tracks:        # No track at that level
            self.tracks[track_level] = newtrack   # so just add it
        else:       # Already a track there, so shunt all higher tracks up one
            occupied_levels = self.get_levels() # Get list of occupied levels...
            occupied_levels.sort()              # ...sort it...
            occupied_levels.reverse()           # ...reverse (highest first)...
            for val in occupied_levels:     
                if val >= track_level:        # Track value >= that to be added
                    self.tracks[val+1] = self.tracks[val] # ..increment by 1
            self.tracks[track_level] = newtrack   # And put the new track in
        self.tracks[track_level].track_level = track_level
        return newtrack

            
    def del_track(self, track_level):
        """ del_track(self, track_level)

            o track_level   Int, the level of the track on the diagram to delete

            Remove the track at the passed level from the diagram
        """
        del self.tracks[track_level]


    def get_tracks(self):
        """ get_tracks(self) -> list

            Returns a list of the tracks contained in the diagram
        """
        return self.tracks.values()


    def move_track(self, from_level, to_level):
        """ move_track(self, from_level, to_level)

            o from_level    Int, the level at which the track to be moved is
                            found

            o to_level      Int, the level to move the track to

            Moves a track from one level on the diagram to another
        """
        aux = self.tracks[from_level]
        del self.tracks[from_level]
        self.add_track(aux, to_level)


    def renumber_tracks(self, low=1, step=1):
        """ renumber_tracks(self, low=1, step=1)

            o low       Int, the track number to start from

            o step      Int, the track interval for separation of tracks

            Reassigns all the tracks to run consecutively from the lowest
            value (low)
        """
        track = low                 # Start numbering from here
        levels = self.get_levels()  # 

        conversion = {}             # Holds new set of levels
        for level in levels:        # Starting at low...
            conversion[track] = self.tracks[level] # Add old tracks to new set
            conversion[track].track_level = track
            track += step                           # step interval
        self.tracks = conversion   # Replace old set of levels with new set

    def get_levels(self):
        """ get_levels(self) -> [int, int, ...]

            Return a sorted list of levels occupied by tracks in the diagram
        """
        levels = self.tracks.keys()
        levels.sort()
        return levels


    def get_drawn_levels(self):
        """ get_drawn_levels(self) -> [int, int, ...]

            Return a sorted list of levels occupied by tracks that are not
            explicitly hidden
        """
        drawn_levels = [key for key in self.tracks.keys() if \
                        not self.tracks[key].hide] # get list of shown levels
        drawn_levels.sort()
        return drawn_levels


    def range(self):
        """ range(self) -> (int, int)

            Returns the lowest and highest base (or mark) numbers containd in
            track features as a tuple
        """
        lows, highs = [], []
        for track in self.tracks.values(): # Get ranges for each track
            low, high = track.range()
            lows.append(low)
            highs.append(high)
        return (min(lows), max(highs))      # Return extremes from all tracks

    def __getitem__(self, key):
        """ __getitem__(self, key) -> Track

            o key       The id of a track in the diagram

            Return the Track object with the passed id
        """
        return self.tracks[key]

    def __str__(self):
        """ __str__(self) -> ""

            Returns a formatted string with information about the diagram
        """
        outstr = ["\n<%s: %s>" % (self.__class__, self.name)]
        outstr.append("%d tracks" % len(self.tracks))
        for level in self.get_levels():
            outstr.append("Track %d: %s\n" % (level, self.tracks[level]))
        outstr = '\n'.join(outstr)
        return outstr       

