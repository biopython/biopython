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
"""Provides a container for information concerning the tracks to be drawn in a diagram.

It also provides the interface for defining the diagram (possibly split these
functions in later version?).

For drawing capabilities, this module uses reportlab to draw and write the
diagram:

http://www.reportlab.com

For dealing with biological information, the package expects BioPython
objects - namely SeqRecord objects containing SeqFeature objects.
"""

try:
    from reportlab.graphics import renderPM
except ImportError:
    # This is an optional part of ReportLab, so may not be installed.
    renderPM = None

from ._LinearDrawer import LinearDrawer
from ._CircularDrawer import CircularDrawer
from ._Track import Track

from Bio.Graphics import _write


def _first_defined(*args):
    """Return the first non-null argument (PRIVATE)."""
    for arg in args:
        if arg is not None:
            return arg
    return None


class Diagram:
    """Diagram container.

    Arguments:
     - name           - a string, identifier for the diagram.
     - tracks         - a list of Track objects comprising the diagram.
     - format         - a string, format of the diagram 'circular' or
       'linear', depending on the sort of diagram required.
     - pagesize       - a string, the pagesize of output describing the ISO
       size of the image, or a tuple of pixels.
     - orientation    - a string describing the required orientation of the
       final drawing ('landscape' or 'portrait').
     - x              - a float (0->1), the proportion of the page to take
       up with even X margins t the page.
     - y              - a float (0->1), the proportion of the page to take
       up with even Y margins to the page.
     - xl             - a float (0->1), the proportion of the page to take
       up with the left X margin to the page (overrides x).
     - xr             - a float (0->1), the proportion of the page to take
       up with the right X margin to the page (overrides x).
     - yt             - a float (0->1), the proportion of the page to take
       up with the top Y margin to the page (overrides y).
     - yb             - a float (0->1), the proportion of the page to take
       up with the bottom Y margin to the page (overrides y).
     - circle_core    - a float, the proportion of the available radius to
       leave empty at the center of a circular diagram (0 to 1).
     - start          - an integer, the base/aa position to start the diagram at.
     - end            - an integer, the base/aa position to end the diagram at.
     - tracklines     - a boolean, True if track guidelines are to be drawn.
     - fragments      - and integer, for a linear diagram, the number of equal
       divisions into which the sequence is divided.
     - fragment_size  - a float (0->1), the proportion of the space
       available to each fragment that should be used in drawing.
     - track_size     - a float (0->1), the proportion of the space
       available to each track that should be used in drawing with sigils.
     - circular       - a boolean, True if the genome/sequence to be drawn
       is, in reality, circular.

    """

    def __init__(
        self,
        name=None,
        format="circular",
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
        tracklines=False,
        fragments=10,
        fragment_size=None,
        track_size=0.75,
        circular=True,
        circle_core=0.0,
    ):
        """Initialize.

        gdd = Diagram(name=None)
        """
        self.tracks = {}  # Holds all Track objects, keyed by level
        self.name = name  # Description of the diagram
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
        if fragment_size is not None:
            self.fragment_size = fragment_size
        else:
            if self.fragments == 1:
                # For single fragments, default to full height
                self.fragment_size = 1
            else:
                # Otherwise keep a 10% gap between fragments
                self.fragment_size = 0.9
        self.track_size = track_size
        self.circular = circular
        self.circle_core = circle_core
        self.cross_track_links = []
        self.drawing = None

    def set_all_tracks(self, attr, value):
        """Set the passed attribute of all tracks in the set to the passed value.

        Arguments:
         - attr    - An attribute of the Track class.
         - value   - The value to set that attribute.

        set_all_tracks(self, attr, value)
        """
        for track in self.tracks.values():
            if hasattr(track, attr):
                # If the feature has the attribute set it to the passed value
                setattr(track, attr, value)

    def draw(
        self,
        format=None,
        pagesize=None,
        orientation=None,
        x=None,
        y=None,
        xl=None,
        xr=None,
        yt=None,
        yb=None,
        start=None,
        end=None,
        tracklines=None,
        fragments=None,
        fragment_size=None,
        track_size=None,
        circular=None,
        circle_core=None,
        cross_track_links=None,
    ):
        """Draw the diagram, with passed parameters overriding existing attributes.

        gdd.draw(format='circular')
        """
        # Pass the parameters to the drawer objects that will build the
        # diagrams.  At the moment, we detect overrides with an or in the
        # Instantiation arguments, but I suspect there's a neater way to do
        # this.
        if format == "linear":
            drawer = LinearDrawer(
                self,
                _first_defined(pagesize, self.pagesize),
                _first_defined(orientation, self.orientation),
                _first_defined(x, self.x),
                _first_defined(y, self.y),
                _first_defined(xl, self.xl),
                _first_defined(xr, self.xr),
                _first_defined(yt, self.yt),
                _first_defined(yb, self.yb),
                _first_defined(start, self.start),
                _first_defined(end, self.end),
                _first_defined(tracklines, self.tracklines),
                _first_defined(fragments, self.fragments),
                _first_defined(fragment_size, self.fragment_size),
                _first_defined(track_size, self.track_size),
                _first_defined(cross_track_links, self.cross_track_links),
            )
        else:
            drawer = CircularDrawer(
                self,
                _first_defined(pagesize, self.pagesize),
                _first_defined(orientation, self.orientation),
                _first_defined(x, self.x),
                _first_defined(y, self.y),
                _first_defined(xl, self.xl),
                _first_defined(xr, self.xr),
                _first_defined(yt, self.yt),
                _first_defined(yb, self.yb),
                _first_defined(start, self.start),
                _first_defined(end, self.end),
                _first_defined(tracklines, self.tracklines),
                _first_defined(track_size, self.track_size),
                _first_defined(circular, self.circular),
                _first_defined(circle_core, self.circle_core),
                _first_defined(cross_track_links, self.cross_track_links),
            )
        drawer.draw()  # Tell the drawer to complete the drawing
        self.drawing = drawer.drawing  # Get the completed drawing

    def write(self, filename="test1.ps", output="PS", dpi=72):
        """Write the drawn diagram to a specified file, in a specified format.

        Arguments:
            - filename   - a string indicating the name of the output file,
              or a handle to write to.
            - output     - a string indicating output format, one of PS, PDF,
              SVG, or provided the ReportLab renderPM module is installed, one
              of the bitmap formats JPG, BMP, GIF, PNG, TIFF or TIFF.  The
              format can be given in upper or lower case.
            - dpi        - an integer. Resolution (dots per inch) for bitmap formats.

        Returns:
            No return value.

        write(self, filename='test1.ps', output='PS', dpi=72)

        """
        return _write(self.drawing, filename, output, dpi=dpi)

    def write_to_string(self, output="PS", dpi=72):
        """Return a byte string containing the diagram in the requested format.

        Arguments:
            - output    - a string indicating output format, one of PS, PDF,
              SVG, JPG, BMP, GIF, PNG, TIFF or TIFF (as specified for the write
              method).
            - dpi       - Resolution (dots per inch) for bitmap formats.

        Returns:
            Return the completed drawing as a bytes string in a prescribed
            format.

        """
        # The ReportLab drawToString method, which this function used to call,
        # originally just used a StringIO handle with the drawToFile method.
        #
        # TODO - Rename this method to include keyword bytes?
        from io import BytesIO

        handle = BytesIO()
        self.write(handle, output, dpi)
        return handle.getvalue()

    def add_track(self, track, track_level):
        """Add a Track object to the diagram.

        It also accepts instructions to place it at a particular level on the
        diagram.

        Arguments:
            - track          - Track object to draw.
            - track_level    - an integer. The level at which the track will be
              drawn (above an arbitrary baseline).

        add_track(self, track, track_level)
        """
        if track is None:
            raise ValueError("Must specify track")
        if track_level not in self.tracks:  # No track at that level
            self.tracks[track_level] = track  # so just add it
        else:  # Already a track there, so shunt all higher tracks up one
            occupied_levels = sorted(
                self.get_levels()
            )  # Get list of occupied levels...
            occupied_levels.reverse()  # ...reverse it (highest first)
            for val in occupied_levels:
                # If track value >= that to be added
                if val >= track.track_level:
                    self.tracks[val + 1] = self.tracks[val]  # ...increment by 1
            self.tracks[track_level] = track  # And put the new track in
        self.tracks[track_level].track_level = track_level

    def new_track(self, track_level, **args):
        """Add a new Track to the diagram at a given level.

        The track is returned for further user manipulation.

        Arguments:
            - track_level   - an integer. The level at which the track will be
              drawn (above an arbitrary baseline).

        new_track(self, track_level)
        """
        newtrack = Track()
        for key in args:
            setattr(newtrack, key, args[key])
        if track_level not in self.tracks:  # No track at that level
            self.tracks[track_level] = newtrack  # so just add it
        else:  # Already a track there, so shunt all higher tracks up one
            occupied_levels = sorted(
                self.get_levels()
            )  # Get list of occupied levels...
            occupied_levels.reverse()  # ...reverse (highest first)...
            for val in occupied_levels:
                if val >= track_level:
                    # Track value >= that to be added, increment by 1
                    self.tracks[val + 1] = self.tracks[val]
            self.tracks[track_level] = newtrack  # And put the new track in
        self.tracks[track_level].track_level = track_level
        return newtrack

    def del_track(self, track_level):
        """Remove the track to be drawn at a particular level on the diagram.

        Arguments:
            - track_level   - an integer. The level of the track on the diagram
              to delete.

        del_track(self, track_level)
        """
        del self.tracks[track_level]

    def get_tracks(self):
        """Return a list of the tracks contained in the diagram."""
        return list(self.tracks.values())

    def move_track(self, from_level, to_level):
        """Move a track from one level on the diagram to another.

        Arguments:
         - from_level   - an integer. The level at which the track to be
           moved is found.
         - to_level     - an integer. The level to move the track to.

        """
        aux = self.tracks[from_level]
        del self.tracks[from_level]
        self.add_track(aux, to_level)

    def renumber_tracks(self, low=1, step=1):
        """Renumber all tracks consecutively.

        Optionally from a passed lowest number.

        Arguments:
         - low     - an integer. The track number to start from.
         - step    - an integer. The track interval for separation of
           tracks.

        """
        track = low  # Start numbering from here
        levels = self.get_levels()

        conversion = {}  # Holds new set of levels
        for level in levels:  # Starting at low...
            conversion[track] = self.tracks[level]  # Add old tracks to new set
            conversion[track].track_level = track
            track += step  # step interval
        self.tracks = conversion  # Replace old set of levels with new set

    def get_levels(self):
        """Return a sorted list of levels occupied by tracks in the diagram."""
        return sorted(self.tracks)

    def get_drawn_levels(self):
        """Return a sorted list of levels occupied by tracks.

        These tracks are not explicitly hidden.
        """
        return sorted(key for key in self.tracks if not self.tracks[key].hide)

    def range(self):
        """Return lowest and highest base numbers from track features.

        Returned type is a tuple.
        """
        lows, highs = [], []
        for track in self.tracks.values():  # Get ranges for each track
            low, high = track.range()
            lows.append(low)
            highs.append(high)
        return min(lows), max(highs)  # Return extremes from all tracks

    def __getitem__(self, key):
        """Return the track contained at the level of the passed key."""
        return self.tracks[key]

    def __str__(self):
        """Return a formatted string describing the diagram."""
        outstr = ["\n<%s: %s>" % (self.__class__, self.name)]
        outstr.append("%d tracks" % len(self.tracks))
        for level in self.get_levels():
            outstr.append("Track %d: %s\n" % (level, self.tracks[level]))
        outstr = "\n".join(outstr)
        return outstr
