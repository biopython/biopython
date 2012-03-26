# Copyright 2002 by Tarjei Mikkelsen.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# get set abstraction for graph representation

class Graph(object):
    """A directed graph abstraction with labeled edges."""

    def __init__(self, nodes = []):
        """Initializes a new Graph object."""
        self._adjacency_list = {}    # maps parent -> set of child objects
        for n in nodes:
            self._adjacency_list[n] = set()
        self._label_map = {}         # maps label -> set of (parent, child) pairs
        self._edge_map = {}          # maps (parent, child) pair -> label

    def __eq__(self, g):
        """Returns true if g is equal to this graph."""
        return isinstance(g, Graph) and \
               (self._adjacency_list == g._adjacency_list) and \
               (self._label_map == g._label_map) and \
               (self._edge_map == g._edge_map)

    def __ne__(self, g):
        """Returns true if g is not equal to this graph."""
        return not self.__eq__(g)

    def __repr__(self):
        """Returns an unique string representation of this graph."""
        s = "<Graph: "
        keys = self._adjacency_list.keys()
        keys.sort()
        for key in keys:
            values = [(x,self._edge_map[(key,x)]) \
                      for x in self._adjacency_list[key].list()]
            values.sort()
            s = s + "(" + repr(key) + ": " + ",".join(map(repr, values)) + ")" 
        return s + ">"

    def __str__(self):
        """Returns a concise string description of this graph."""
        nodenum = len(self._adjacency_list.keys())
        edgenum = reduce(lambda x,y: x+y,
                         map(len, self._adjacency_list.values()))
        labelnum = len(self._label_map.keys())
        return "<Graph: " + \
               str(nodenum) + " node(s), " + \
               str(edgenum) + " edge(s), " + \
               str(labelnum) + " unique label(s)>"

    def add_node(self, node):
        """Adds a node to this graph."""
        if node not in self._adjacency_list:
            self._adjacency_list[node] = set()

    def add_edge(self, source, to, label = None):
        """Adds an edge to this graph."""
        if source not in self._adjacency_list:
            raise ValueError("Unknown <from> node: " + str(source))
        if to not in self._adjacency_list:
            raise ValueError("Unknown <to> node: " + str(to))
        if (source,to) in self._edge_map:
            raise ValueError(str(source) + " -> " + str(to) + " exists")
        self._adjacency_list[source].add(to)
        if label not in self._label_map:
            self._label_map[label] = set()
        self._label_map[label].add((source,to))
        self._edge_map[(source,to)] = label

    def child_edges(self, parent):
        """Returns a list of (child, label) pairs for parent."""
        if parent not in self._adjacency_list:
            raise ValueError("Unknown <parent> node: " + str(parent))
        return [(x, self._edge_map[(parent,x)]) \
                for x in sorted(self._adjacency_list[parent])]

    def children(self, parent):
        """Returns a list of unique children for parent."""
        return sorted(self._adjacency_list[parent])

    def edges(self, label):
        """Returns a list of all the edges with this label."""
        if label not in self._label_map:
            raise ValueError("Unknown label: " + str(label))
        return self._label_map[label].list()

    def labels(self):
        """Returns a list of all the edge labels in this graph."""
        return self._label_map.keys()

    def nodes(self):
        """Returns a list of the nodes in this graph."""
        return self._adjacency_list.keys()

    def parent_edges(self, child):
        """Returns a list of (parent, label) pairs for child."""
        if child not in self._adjacency_list:
            raise ValueError("Unknown <child> node: " + str(child))
        parents = []
        for parent, children in self._adjacency_list.iteritems():
            for x in children:
                if x is child:
                    parents.append((parent, self._edge_map[(parent, child)]))
        return sorted(parents)

    def parents(self, child):
        """Returns a list of unique parents for child."""
        return sorted(set([x[0] for x in self.parent_edges(child)]))

    def remove_node(self, node):
        """Removes node and all edges connected to it."""
        if node not in self._adjacency_list:
            raise ValueError("Unknown node: " + str(node))
        # remove node (and all out-edges) from adjacency list
        del self._adjacency_list[node]
        # remove all in-edges from adjacency list
        for n in self._adjacency_list.keys():
            self._adjacency_list[n] = set(x for x in self._adjacency_list[n] \
                                          if x is not node)
        # remove all refering pairs in label map
        for label in self._label_map.keys():
            lm = set(x for x in self._label_map[label] \
                     if (x[0] is not node) and (x[1] is not node))
            # remove the entry completely if the label is now unused
            if lm:
                self._label_map[label] = lm
            else:
                del self._label_map[label]
        # remove all refering entries in edge map
        for edge in self._edge_map.keys():
            if edge[0] is node or edge[1] is node:
                del self._edge_map[edge]
        
    def remove_edge(self, parent, child, label):
        """Removes edge. -- NOT IMPLEMENTED"""
        # hm , this is a multigraph - how should this be implemented?
        raise NotImplementedError("remove_edge is not yet implemented")



