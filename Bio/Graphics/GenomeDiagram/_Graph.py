# Copyright 2003-2008 by Leighton Pritchard.  All rights reserved.
# Revisions copyright 2008-2009 by Peter Cock.
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

"""Graph module.

Provides:
 - GraphData - Contains data from which a graph will be drawn, and
   information about its presentation

For drawing capabilities, this module uses reportlab to draw and write
the diagram: http://www.reportlab.com
"""

# ReportLab imports

from reportlab.lib import colors

from math import sqrt


class GraphData:
    """Graph Data.

    Attributes:
     - id    Unique identifier for the data
     - data  Dictionary of describing the data, keyed by position
     - name  String describing the data
     - style String ('bar', 'heat', 'line') describing how to draw the data
     - poscolor     colors.Color for drawing high (some styles) or all
       values
     - negcolor     colors.Color for drawing low values (some styles)
     - linewidth     Int, thickness to draw the line in 'line' styles

    """

    def __init__(
        self,
        id=None,
        data=None,
        name=None,
        style="bar",
        color=colors.lightgreen,
        altcolor=colors.darkseagreen,
        center=None,
        colour=None,
        altcolour=None,
    ):
        """Initialize.

        Arguments:
         - id    Unique ID for the graph
         - data  List of (position, value) tuples
         - name  String describing the graph
         - style String describing the presentation style ('bar', 'line',
           'heat')
         - color   colors.Color describing the color to draw all or the
           'high' (some styles) values (overridden by backwards
           compatible argument with UK spelling, colour).
         - altcolor colors.Color describing the color to draw the 'low'
           values (some styles only) (overridden by backwards
           compatible argument with UK spelling, colour).
         - center Value at which x-axis crosses y-axis.

        """
        # Let the UK spelling (colour) override the USA spelling (color)
        if colour is not None:
            color = colour
        if altcolour is not None:
            altcolor = altcolour

        self.id = id  # Unique identifier for the graph
        self.data = {}  # holds values, keyed by sequence position
        if data is not None:
            self.set_data(data)
        self.name = name  # Descriptive string

        # Attributes describing how the graph will be drawn
        self.style = style  # One of 'bar', 'heat' or 'line'
        self.poscolor = color  # Color to draw all, or 'high' values
        self.negcolor = altcolor  # Color to draw 'low' values
        self.linewidth = 2  # linewidth to use in line graphs
        self.center = center  # value at which x-axis crosses y-axis

    def set_data(self, data):
        """Add data as a list of (position, value) tuples."""
        for pos, val in data:  # Fill data dictionary
            self.data[pos] = val

    def get_data(self):
        """Return data as a list of sorted (position, value) tuples."""
        data = []
        for xval in self.data:
            yval = self.data[xval]
            data.append((xval, yval))
        data.sort()
        return data

    def add_point(self, point):
        """Add a single point to the set of data as a (position, value) tuple."""
        pos, val = point
        self.data[pos] = val

    def quartiles(self):
        """Return (minimum, lowerQ, medianQ, upperQ, maximum) values as tuple."""
        data = sorted(self.data.values())
        datalen = len(data)
        return (
            data[0],
            data[datalen // 4],
            data[datalen // 2],
            data[3 * datalen // 4],
            data[-1],
        )

    def range(self):
        """Return range of data as (start, end) tuple.

        Returns the range of the data, i.e. its start and end points on
        the genome as a (start, end) tuple.
        """
        positions = sorted(self.data)  # i.e. dict keys
        # Return first and last positions in graph
        # print(len(self.data))
        return (positions[0], positions[-1])

    def mean(self):
        """Return the mean value for the data points (float)."""
        data = list(self.data.values())
        return sum(data) / len(data)

    def stdev(self):
        """Return the sample standard deviation for the data (float)."""
        data = list(self.data.values())
        m = self.mean()
        runtotal = 0.0
        for entry in data:
            runtotal += (entry - m) ** 2
        # This is sample standard deviation; population stdev would involve
        # division by len(data), rather than len(data)-1
        return sqrt(runtotal / (len(data) - 1))

    def __len__(self):
        """Return the number of points in the data set."""
        return len(self.data)

    def __getitem__(self, index):
        """Return data value(s) at the given position.

        Given an integer representing position on the sequence
        returns a float - the data value at the passed position.

        If a slice, returns graph data from the region as a list or
        (position, value) tuples. Slices with step are not supported.
        """
        if isinstance(index, int):
            return self.data[index]
        elif isinstance(index, slice):
            # TODO - Why does it treat the end points both as inclusive?
            # This doesn't match Python norms does it?
            low = index.start
            high = index.stop
            if index.step is not None and index.step != 1:
                raise ValueError
            outlist = []
            for pos in sorted(self.data):
                if pos >= low and pos <= high:
                    outlist.append((pos, self.data[pos]))
            return outlist
        else:
            raise TypeError("Need an integer or a slice")

    def __str__(self):
        """Return a string describing the graph data."""
        outstr = [f"\nGraphData: {self.name}, ID: {self.id}"]
        outstr.append("Number of points: %d" % len(self.data))
        outstr.append(f"Mean data value: {self.mean()}")
        outstr.append(f"Sample SD: {self.stdev():.3f}")
        outstr.append(
            "Minimum: %s\n1Q: %s\n2Q: %s\n3Q: %s\nMaximum: %s" % self.quartiles()
        )
        outstr.append("Sequence Range: %s..%s" % self.range())
        return "\n".join(outstr)
