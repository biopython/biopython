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

""" Graph module

    Provides:

    o GraphData - Contains data from which a graph will be drawn, and
                    information about its presentation

    For drawing capabilities, this module uses reportlab to draw and write
    the diagram:

    http://www.reportlab.com

    For dealing with biological information, the package expects BioPython
    objects:

    http://www.biopython.org
"""

# ReportLab imports
from __future__ import print_function

from reportlab.lib import colors

from math import sqrt


class GraphData(object):
    """ GraphData

        Provides:

        Methods:

        o __init__(self, id=None, data=None, name=None, style='bar',
                 color=colors.lightgreen, altcolor=colors.darkseagreen)
                 Called on instantiation

        o set_data(self, data)  Load the object with data to be plotted

        o get_data(self)    Returns the data to be plotted as a list of
                            (position, value) tuples

        o add_point(self, point)    Add a single point to the data set

        o quartiles(self)   Returns a tuple of the data quartiles

        o range(self)   Returns a tuple of the base range covered by the graph
                        data

        o mean(self)    Returns a float of the mean data point value

        o stdev(self)   Returns the sample standard deviation of the data values

        o __len__(self) Returns the length of sequence covered by the data

        o __getitem__(self, index)  Returns the value at the base specified,
                                    or graph data in the base range

        o __str__(self) Returns a formatted string describing the graph data

        Attributes:

        o id    Unique identifier for the data

        o data  Dictionary of describing the data, keyed by position

        o name  String describing the data

        o style String ('bar', 'heat', 'line') describing how to draw the data

        o poscolor     colors.Color for drawing high (some styles) or all
                        values

        o negcolor     colors.Color for drawing low values (some styles)

        o linewidth     Int, thickness to draw the line in 'line' styles
    """
    def __init__(self, id=None, data=None, name=None, style='bar',
                 color=colors.lightgreen, altcolor=colors.darkseagreen,
                 center=None, colour=None, altcolour=None):
        """__init__(self, id=None, data=None, name=None, style='bar',
                 color=colors.lightgreen, altcolor=colors.darkseagreen)

            o id    Unique ID for the graph

            o data  List of (position, value) tuples

            o name  String describing the graph

            o style String describing the presentation style ('bar', 'line',
                    'heat')

            o color   colors.Color describing the color to draw all or the
                      'high' (some styles) values (overridden by backwards
                      compatible argument with UK spelling, colour).

            o altcolor colors.Color describing the color to draw the 'low'
                       values (some styles only) (overridden by backwards
                       compatible argument with UK spelling, colour).

            o center Value at which x-axis crosses y-axis.

        """

        # Let the UK spelling (colour) override the USA spelling (color)
        if colour is not None:
            color = colour
        if altcolour is not None:
            altcolor = altcolour

        self.id = id            # Unique identifier for the graph
        self.data = {}          # holds values, keyed by sequence position
        if data is not None:
            self.set_data(data)
        self.name = name        # Descriptive string

        # Attributes describing how the graph will be drawn
        self.style = style          # One of 'bar', 'heat' or 'line'
        self.poscolor = color     # Color to draw all, or 'high' values
        self.negcolor = altcolor  # Color to draw 'low' values
        self.linewidth = 2          # linewidth to use in line graphs
        self.center = center        # value at which x-axis crosses y-axis

    def set_data(self, data):
        """ set_data(self, data)

            o data      List of (position, value) tuples

            Add data with a list of (position, value) tuples
        """
        for (pos, val) in data:     # Fill data dictionary
            self.data[pos] = val

    def get_data(self):
        """ get_data(self) -> [(int, float), (int, float), ...]

            Return data as a list of sorted (position, value) tuples
        """
        data = []
        for xval in self.data:
            yval = self.data[xval]
            data.append((xval, yval))
        data.sort()
        return data

    def add_point(self, point):
        """ add_point(self, point)

            o point     (position, value) tuple

            Add a single point to the set of data
        """
        pos, val = point
        self.data[pos] = val

    def quartiles(self):
        """ quartiles(self) -> (float, float, float, float, float)

            Returns the (minimum, lowerQ, medianQ, upperQ, maximum) values as
            a tuple
        """
        data = sorted(self.data.values())
        datalen = len(data)
        return(data[0], data[datalen//4], data[datalen//2],
               data[3*datalen//4], data[-1])

    def range(self):
        """ range(self) -> (int, int)

            Returns the range of the data, i.e. its start and end points on
            the genome as a (start, end) tuple
        """
        positions = sorted(self.data)  # i.e. dict keys
        # Return first and last positions in graph
        # print len(self.data)
        return (positions[0], positions[-1])

    def mean(self):
        """ mean(self) -> Float

            Returns the mean value for the data points
        """
        data = list(self.data.values())
        sum = 0.
        for item in data:
            sum += float(item)
        return sum/len(data)

    def stdev(self):
        """ stdev(self) -> Float

            Returns the sample standard deviation for the data
        """
        data = list(self.data.values())
        m = self.mean()
        runtotal = 0.
        for entry in data:
            runtotal += float((entry - m)**2)
        # This is sample standard deviation; population stdev would involve
        # division by len(data), rather than len(data)-1
        return sqrt(runtotal/(len(data)-1))

    def __len__(self):
        """ __len__(self) -> Int

            Returns the number of points in the data set
        """
        return len(self.data)

    def __getitem__(self, index):
        """ __getitem__(self, index) -> Float or list of tuples

            Given an integer representing position on the sequence
            returns a float - the data value at the passed position.

            If a slice, returns graph data from the region as a list or
            (position, value) tuples. Slices with step are not supported.

            Returns the data value at the passed position
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
                if pos >= low and pos <=high:
                    outlist.append((pos, self.data[pos]))
            return outlist
        else:
            raise TypeError("Need an integer or a slice")

    def __str__(self):
        """ __str__(self) -> ""

            Returns a string describing the graph data
        """
        outstr = ["\nGraphData: %s, ID: %s" % (self.name, self.id)]
        outstr.append("Number of points: %d" % len(self.data))
        outstr.append("Mean data value: %s" % self.mean())
        outstr.append("Sample SD: %.3f" % self.stdev())
        outstr.append("Minimum: %s\n1Q: %s\n2Q: %s\n3Q: %s\nMaximum: %s" % self.quartiles())
        outstr.append("Sequence Range: %s..%s" % self.range())
        return "\n".join(outstr)
