# Copyright 2002 by Tarjei Mikkelsen.  All rights reserved.
# Revisions copyright 2018 by Maximilian Greil. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""get/set abstraction for graph representation."""

from functools import reduce


class Graph:
    """A directed graph abstraction with labeled edges."""

    def __init__(self, nodes=()):
        """Initialize a new Graph object."""
        self._adjacency_list = {}  # maps parent -> set of child objects
        for n in nodes:
            self._adjacency_list[n] = set()
        self._label_map = {}  # maps label -> set of (parent, child) pairs
        self._edge_map = {}  # maps (parent, child) pair -> label

    def __eq__(self, g):
        """Return true if g is equal to this graph."""
        return (
            isinstance(g, Graph)
            and self._adjacency_list == g._adjacency_list
            and self._label_map == g._label_map
            and self._edge_map == g._edge_map
        )

    def __repr__(self):
        """Return a unique string representation of this graph."""
        s = "<Graph: "
        for key in sorted(self._adjacency_list):
            values = sorted(
                (x, self._edge_map[(key, x)]) for x in list(self._adjacency_list[key])
            )
            s += f"({key!r}: {','.join(repr(v) for v in values)})"
        return s + ">"

    def __str__(self):
        """Return a concise string description of this graph."""
        nodenum = len(self._adjacency_list)
        edgenum = reduce(
            lambda x, y: x + y, [len(v) for v in self._adjacency_list.values()]
        )
        labelnum = len(self._label_map)
        return (
            "<Graph: "
            + str(nodenum)
            + " node(s), "
            + str(edgenum)
            + " edge(s), "
            + str(labelnum)
            + " unique label(s)>"
        )

    def add_node(self, node):
        """Add a node to this graph."""
        if node not in self._adjacency_list:
            self._adjacency_list[node] = set()

    def add_edge(self, source, to, label=None):
        """Add an edge to this graph."""
        if source not in self._adjacency_list:
            raise ValueError("Unknown <from> node: " + str(source))
        if to not in self._adjacency_list:
            raise ValueError("Unknown <to> node: " + str(to))
        if (source, to) in self._edge_map:
            raise ValueError(str(source) + " -> " + str(to) + " exists")
        self._adjacency_list[source].add(to)
        if label not in self._label_map:
            self._label_map[label] = set()
        self._label_map[label].add((source, to))
        self._edge_map[(source, to)] = label

    def child_edges(self, parent):
        """Return a list of (child, label) pairs for parent."""
        if parent not in self._adjacency_list:
            raise ValueError("Unknown <parent> node: " + str(parent))
        return [
            (x, self._edge_map[(parent, x)])
            for x in sorted(self._adjacency_list[parent])
        ]

    def children(self, parent):
        """Return a list of unique children for parent."""
        return sorted(self._adjacency_list[parent])

    def edges(self, label):
        """Return a list of all the edges with this label."""
        if label not in self._label_map:
            raise ValueError("Unknown label: " + str(label))
        return sorted(self._label_map[label])

    def labels(self):
        """Return a list of all the edge labels in this graph."""
        return sorted(self._label_map.keys())

    def nodes(self):
        """Return a list of the nodes in this graph."""
        return list(self._adjacency_list.keys())

    def parent_edges(self, child):
        """Return a list of (parent, label) pairs for child."""
        if child not in self._adjacency_list:
            raise ValueError("Unknown <child> node: " + str(child))
        parents = []
        for parent, children in self._adjacency_list.items():
            for x in children:
                if x == child:
                    parents.append((parent, self._edge_map[(parent, child)]))
        return sorted(parents)

    def parents(self, child):
        """Return a list of unique parents for child."""
        return sorted({x[0] for x in self.parent_edges(child)})

    def remove_node(self, node):
        """Remove node and all edges connected to it."""
        if node not in self._adjacency_list:
            raise ValueError("Unknown node: " + str(node))
        # remove node (and all out-edges) from adjacency list
        del self._adjacency_list[node]
        # remove all in-edges from adjacency list
        for n in self._adjacency_list.keys():
            self._adjacency_list[n] = {x for x in self._adjacency_list[n] if x != node}
        # remove all referring pairs in label map
        for label in list(self._label_map.keys()):  # we're editing this!
            lm = {
                x for x in self._label_map[label] if (x[0] != node) and (x[1] != node)
            }
            # remove the entry completely if the label is now unused
            if lm:
                self._label_map[label] = lm
            else:
                del self._label_map[label]
        # remove all referring entries in edge map
        for edge in list(self._edge_map.keys()):  # we're editing this!
            if edge[0] == node or edge[1] == node:
                del self._edge_map[edge]

    def remove_edge(self, parent, child, label):
        """Remove edge (NOT IMPLEMENTED)."""
        # hm , this is a multigraph - how should this be implemented?
        raise NotImplementedError("remove_edge is not yet implemented")
