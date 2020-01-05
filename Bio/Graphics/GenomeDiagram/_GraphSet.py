# Copyright 2003-2008 by Leighton Pritchard.  All rights reserved.
# Revisions copyright 2008-2010 by Peter Cock.
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
#
# TODO: Make representation of Ymax and Ymin values at this level, so that
#       calculation of graph/axis drawing is simplified

"""GraphSet module.

Provides:
 - GraphSet - container for GraphData objects

For drawing capabilities, this module uses reportlab to draw and write
the diagram: http://www.reportlab.com
"""

# ReportLab imports

from reportlab.lib import colors

from ._Graph import GraphData


class GraphSet:
    """Graph Set.

    Attributes:
     - id        Unique identifier for the set
     - name      String describing the set

    """

    def __init__(self, name=None):
        """Initialize.

        Arguments:
         - name      String identifying the graph set sensibly

        """
        self.id = id  # Unique identifier for the set
        self._next_id = 0  # Holds unique ids for graphs
        self._graphs = {}  # Holds graphs, keyed by unique id
        self.name = name  # Holds description of graph

    def new_graph(
        self,
        data,
        name=None,
        style="bar",
        color=colors.lightgreen,
        altcolor=colors.darkseagreen,
        linewidth=1,
        center=None,
        colour=None,
        altcolour=None,
        centre=None,
    ):
        """Add a GraphData object to the diagram.

        Arguments:
         - data      List of (position, value) int tuples
         - name      String, description of the graph
         - style     String ('bar', 'heat', 'line') describing how the graph
           will be drawn
         - color    colors.Color describing the color to draw all or 'high'
           (some styles) data (overridden by backwards compatible
           argument with UK spelling, colour).
         - altcolor  colors.Color describing the color to draw 'low' (some
           styles) data (overridden by backwards compatible argument
           with UK spelling, colour).
         - linewidth     Float describing linewidth for graph
         - center        Float setting the value at which the x-axis
           crosses the y-axis (overridden by backwards
           compatible argument with UK spelling, centre)

        Add a GraphData object to the diagram (will be stored internally).
        """
        # Let the UK spelling (colour) override the USA spelling (color)
        if colour is not None:
            color = colour
        if altcolour is not None:
            altcolor = altcolour
        if centre is not None:
            center = centre

        id = self._next_id  # get id number
        graph = GraphData(id, data, name, style, color, altcolor, center)
        graph.linewidth = linewidth
        self._graphs[id] = graph  # add graph data
        self._next_id += 1  # increment next id
        return graph

    def del_graph(self, graph_id):
        """Remove a graph from the set, indicated by its id."""
        del self._graphs[graph_id]

    def get_graphs(self):
        """Return list of all graphs in the graph set, sorted by id.

        Sorting is to ensure reliable stacking.
        """
        return [self._graphs[id] for id in sorted(self._graphs)]

    def get_ids(self):
        """Return a list of all ids for the graph set."""
        return list(self._graphs.keys())

    def range(self):
        """Return the lowest and highest base (or mark) numbers as a tuple."""
        lows, highs = [], []
        for graph in self._graphs.values():
            low, high = graph.range()
            lows.append(low)
            highs.append(high)
        return (min(lows), max(highs))

    def data_quartiles(self):
        """Return (minimum, lowerQ, medianQ, upperQ, maximum) values as a tuple."""
        data = []
        for graph in self._graphs.values():
            data += list(graph.data.values())
        data.sort()
        datalen = len(data)
        return (
            data[0],
            data[datalen / 4],
            data[datalen / 2],
            data[3 * datalen / 4],
            data[-1],
        )

    def to_string(self, verbose=0):
        """Return a formatted string with information about the set.

        Arguments:
            - verbose - Flag indicating whether a short or complete account
              of the set is required

        """
        if not verbose:
            return "%s" % self
        else:
            outstr = ["\n<%s: %s>" % (self.__class__, self.name)]
            outstr.append("%d graphs" % len(self._graphs))
            for key in self._graphs:
                outstr.append("%s" % self._graphs[key])
            return "\n".join(outstr)

    def __len__(self):
        """Return the number of graphs in the set."""
        return len(self._graphs)

    def __getitem__(self, key):
        """Return a graph, keyed by id."""
        return self._graphs[key]

    def __str__(self):
        """Return a formatted string with information about the feature set."""
        outstr = ["\n<%s: %s>" % (self.__class__, self.name)]
        outstr.append("%d graphs" % len(self._graphs))
        outstr = "\n".join(outstr)
        return outstr
