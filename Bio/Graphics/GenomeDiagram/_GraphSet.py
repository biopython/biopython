# Copyright 2003-2008 by Leighton Pritchard.  All rights reserved.
# Revisions copyright 2008-2010 by Peter Cock.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
# Contact:       Leighton Pritchard, Scottish Crop Research Institute,
#                Invergowrie, Dundee, Scotland, DD2 5DA, UK
#                L.Pritchard@scri.ac.uk
################################################################################
#
# TODO: Make representation of Ymax and Ymin values at this level, so that
#       calculation of graph/axis drawing is simplified

""" GraphSet module

    Provides:

    o GraphSet - container for GraphData objects

    For drawing capabilities, this module uses reportlab to draw and write
    the diagram:

    http://www.reportlab.com

    For dealing with biological information, the package expects BioPython
    objects:

    http://www.biopython.org
"""

# ReportLab imports
from reportlab.lib import colors

from _Graph import GraphData

class GraphSet(object):
    """ GraphSet

        Provides:

        Methods:

        o __init__(self, set_id=None, name=None)    Called on instantiation

        o new_graph(self, data, name, style='bar', color=colors.lightgreen,
                  altcolor=colors.darkseagreen)    Create new graph in the set
                          from the passed data, with the passed parameters

        o del_graph(self, graph_id) Delete graph with the passed id

        o get_graphs(self)  Returns a list of all graphs

        o get_ids(self)     Returns a list of graph ids

        o range(self)       Returns the range covered by the graphs in the set

        o to_string(self, verbose=0)    Returns a string describing the set

        o __len__(self)     Returns the length of sequence covered by the set

        o __getitem__(self, key)    Returns the graph with the id of the passed key

        o __str__(self)     Returns a string describing the set

        Attributes:

        o id        Unique identifier for the set

        o name      String describing the set

    """
    def __init__(self, name=None):
        """ __init__(self, name=None)

            o name      String identifying the graph set sensibly
        """
        self.id = id            # Unique identifier for the set
        self._next_id = 0       # Holds unique ids for graphs
        self._graphs = {}       # Holds graphs, keyed by unique id
        self.name = name        # Holds description of graph


    def new_graph(self, data, name=None, style='bar', color=colors.lightgreen,
                  altcolor=colors.darkseagreen, linewidth=1, center=None,
                  colour=None, altcolour=None, centre=None):
        """ new_graph(self, data, name=None, style='bar', color=colors.lightgreen,
                  altcolor=colors.darkseagreen)

            o data      List of (position, value) int tuples

            o name      String, description of the graph

            o style     String ('bar', 'heat', 'line') describing how the graph
                        will be drawn

            o color    colors.Color describing the color to draw all or 'high'
                       (some styles) data (overridden by backwards compatible
                       argument with UK spelling, colour).

            o altcolor  colors.Color describing the color to draw 'low' (some
                        styles) data (overridden by backwards compatible argument
                        with UK spelling, colour).
            
            o linewidth     Float describing linewidth for graph

            o center        Float setting the value at which the x-axis
                            crosses the y-axis (overridden by backwards
                            compatible argument with UK spelling, centre)

            Add a GraphData object to the diagram (will be stored
            internally
        """
        #Let the UK spelling (colour) override the USA spelling (color)
        if colour is not None:
            color = colour
        if altcolour is not None:
            altcolor = altcolour
        if centre is not None:
            center = centre

        id = self._next_id                              # get id number
        graph = GraphData(id, data, name, style, color, altcolor, center)
        graph.linewidth = linewidth
        self._graphs[id] =  graph                       # add graph data
        self._next_id += 1                              # increment next id
        return graph


    def del_graph(self, graph_id):
        """ del_graph(self, graph_id)

            o graph_id        Identifying value of the graph

            Remove a graph from the set, indicated by its id
        """
        del self._graphs[graph_id]


    def get_graphs(self):
        """ get_graphs(self) -> [Graph, Graph, ...]

            Return a list of all graphs in the graph set, sorted by id (for
            reliable stacking...)
        """
        ids = self._graphs.keys()
        ids.sort()
        return [self._graphs[id] for id in ids]


    def get_ids(self):
        """ get_ids(self) -> [int, int, ...]

            Return a list of all ids for the graph set
        """
        return self._graphs.keys()


    def range(self):
        """ range(self) -> (int, int)

            Returns the lowest and highest base (or mark) numbers as a tuple
        """
        lows, highs = [], []
        for graph in self._graphs.values():
            low, high = graph.range()
            lows.append(low)
            highs.append(high)
        return (min(lows), max(highs))


    def data_quartiles(self):
        """ data_quartiles(self) -> (float, float, float, float, float)

            Returns the (minimum, lowerQ, medianQ, upperQ, maximum) values as
            a tuple
        """
        data = []
        for graph in self._graphs.values():
            data += graph.data.values()
        data.sort()
        datalen = len(data)
        return(data[0], data[datalen/4], data[datalen/2],
               data[3*datalen/4], data[-1])


    def to_string(self, verbose=0):
        """ to_string(self, verbose=0) -> ""

            o verbose       Flag indicating whether a short or complete account
                            of the set is required

            Returns a formatted string with information about the set
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
        """ __len__(self) -> int

            Return the number of graphs in the set
        """
        return len(self._graphs)


    def __getitem__(self, key):
        """ __getitem__(self, key) -> Graph

            Return a graph, keyed by id
        """
        return self._graphs[key]


    def __str__(self):
        """ __str__(self) -> ""

            Returns a formatted string with information about the feature set
        """
        outstr = ["\n<%s: %s>" % (self.__class__, self.name)]
        outstr.append("%d graphs" % len(self._graphs))
        outstr = "\n".join(outstr)
        return outstr


################################################################################
# RUN AS SCRIPT
################################################################################

if __name__ == '__main__':

    # Test code    
    gdgs = GraphSet(0, 'test data')

    testdata1 = [(1, 10), (5, 15), (10, 20), (20, 40)]
    testdata2 = [(250, .34), (251, .7), (252, .7), (253, .54), (254, .65)]

    gdgs.add_graph(testdata1, 'TestData 1')
    gdgs.add_graph(testdata2, 'TestData 2')

    print gdgs
