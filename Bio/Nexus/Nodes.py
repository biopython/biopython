# Copyright 2005-2008 by Frank Kauff & Cymon J. Cox. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Linked list functionality for use in Bio.Nexus.

Provides functionality of a linked list.
Each node has one (or none) predecessor, and an arbitrary number of successors.
Nodes can store arbitrary data in a NodeData class.

Subclassed by Nexus.Trees to store phylogenetic trees.

Bug reports to Frank Kauff (fkauff@biologie.uni-kl.de)
"""

from typing import Dict, List, Optional


class ChainException(Exception):
    """Provision for the management of Chain exceptions."""


class NodeException(Exception):
    """Provision for the management of Node exceptions."""


class Chain:
    """Stores a list of nodes that are linked together."""

    def __init__(self) -> None:
        """Initialize a node chain."""
        self.chain: Dict[int, "Node"] = {}
        self.id = -1

    def _get_id(self) -> int:
        """Get a new id for a node in the chain (PRIVATE)."""
        self.id += 1
        return self.id

    def all_ids(self) -> List[int]:
        """Return a list of all node ids."""
        return list(self.chain.keys())

    def add(self, node: "Node", prev: Optional[int] = None) -> int:
        """Attach node to another."""
        if prev is not None and prev not in self.chain:
            raise ChainException("Unknown predecessor: " + str(prev))
        else:
            id = self._get_id()
            node.set_id(id)
            node.set_prev(prev)
            if prev is not None:
                self.chain[prev].add_succ(id)
            self.chain[id] = node
        return id

    def collapse(self, id):
        """Delete node from chain and relinks successors to predecessor."""
        if id not in self.chain:
            raise ChainException("Unknown ID: " + str(id))
        prev_id = self.chain[id].get_prev()
        self.chain[prev_id].remove_succ(id)
        succ_ids = self.chain[id].get_succ()
        for i in succ_ids:
            self.chain[i].set_prev(prev_id)
        self.chain[prev_id].add_succ(succ_ids)
        node = self.chain[id]
        self.kill(id)
        return node

    def kill(self, id):
        """Kill a node from chain without caring to what it is connected."""
        if id not in self.chain:
            raise ChainException("Unknown ID: " + str(id))
        else:
            del self.chain[id]

    def unlink(self, id):
        """Disconnect node from his predecessor."""
        if id not in self.chain:
            raise ChainException("Unknown ID: " + str(id))
        else:
            prev_id = self.chain[id].prev
            if prev_id is not None:
                self.chain[prev_id].succ.pop(self.chain[prev_id].succ.index(id))
            self.chain[id].prev = None
            return prev_id

    def link(self, parent, child):
        """Connect son to parent."""
        if child not in self.chain:
            raise ChainException("Unknown ID: " + str(child))
        elif parent not in self.chain:
            raise ChainException("Unknown ID: " + str(parent))
        else:
            self.unlink(child)
            self.chain[parent].succ.append(child)
            self.chain[child].set_prev(parent)

    def is_parent_of(self, parent, grandchild):
        """Check if grandchild is a subnode of parent."""
        if grandchild == parent or grandchild in self.chain[parent].get_succ():
            return True
        else:
            for sn in self.chain[parent].get_succ():
                if self.is_parent_of(sn, grandchild):
                    return True
            else:
                return False

    def trace(self, start, finish):
        """Return a list of all node_ids between two nodes (excluding start, including end)."""
        if start not in self.chain or finish not in self.chain:
            raise NodeException("Unknown node.")
        if not self.is_parent_of(start, finish) or start == finish:
            return []
        for sn in self.chain[start].get_succ():
            if self.is_parent_of(sn, finish):
                return [sn] + self.trace(sn, finish)


class Node:
    """A single node."""

    def __init__(self, data=None):
        """Represent a node with one predecessor and multiple successors."""
        self.id = None
        self.data = data
        self.prev = None
        self.succ = []

    def set_id(self, id):
        """Set the id of a node, if not set yet."""
        if self.id is not None:
            raise NodeException("Node id cannot be changed.")
        self.id = id

    def get_id(self):
        """Return the node's id."""
        return self.id

    def get_succ(self):
        """Return a list of the node's successors."""
        return self.succ

    def get_prev(self):
        """Return the id of the node's predecessor."""
        return self.prev

    def add_succ(self, id):
        """Add a node id to the node's successors."""
        if isinstance(id, type([])):
            self.succ.extend(id)
        else:
            self.succ.append(id)

    def remove_succ(self, id):
        """Remove a node id from the node's successors."""
        self.succ.remove(id)

    def set_succ(self, new_succ):
        """Set the node's successors."""
        if not isinstance(new_succ, type([])):
            raise NodeException("Node successor must be of list type.")
        self.succ = new_succ

    def set_prev(self, id):
        """Set the node's predecessor."""
        self.prev = id

    def get_data(self):
        """Return a node's data."""
        return self.data

    def set_data(self, data):
        """Set a node's data."""
        self.data = data
