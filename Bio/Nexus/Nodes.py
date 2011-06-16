# Copyright 2005-2008 by Frank Kauff & Cymon J. Cox. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.
#
# Nodes.py
# 
# Provides functionality of a linked list.
# Each node has one (or none) predecessor, and an arbitrary number of successors.
# Nodes can store arbitrary data in a NodeData class.
#
# Subclassed by Nexus.Trees to store phylogenetic trees.
#
# Bug reports to Frank Kauff (fkauff@biologie.uni-kl.de)
#

class ChainException(Exception):
    pass

class NodeException(Exception):
    pass

class Chain(object):
    """Stores a list of nodes that are linked together."""
    
    def __init__(self):
        """Initiates a node chain: (self)."""
        self.chain={}
        self.id=-1

    def _get_id(self):
        """Gets a new id for a node in the chain."""
        self.id+=1
        return self.id 
   
    def all_ids(self):
        """Return a list of all node ids."""
        return self.chain.keys()

    def add(self,node,prev=None):
        """Attaches node to another: (self, node, prev)."""
        if prev is not None and prev not in self.chain:
            raise ChainException('Unknown predecessor: '+str(prev))
        else:
            id=self._get_id()
            node.set_id(id)
            node.set_prev(prev)
            if prev is not None:
                self.chain[prev].add_succ(id)
            self.chain[id]=node
        return id

    def collapse(self,id):
        """Deletes node from chain and relinks successors to predecessor: collapse(self, id)."""
        if id not in self.chain:
            raise ChainException('Unknown ID: '+str(id))
        prev_id=self.chain[id].get_prev()
        self.chain[prev_id].remove_succ(id)
        succ_ids=self.chain[id].get_succ()
        for i in succ_ids:
            self.chain[i].set_prev(prev_id)
        self.chain[prev_id].add_succ(succ_ids)
        node=self.chain[id]
        self.kill(id)
        return node

    def kill(self,id):
        """Kills a node from chain without caring to what it is connected: kill(self,id)."""
        if id not in self.chain:
            raise ChainException('Unknown ID: '+str(id))
        else:
            del self.chain[id]

    def unlink(self,id):
        """Disconnects node from his predecessor: unlink(self,id)."""
        if id not in self.chain:
            raise ChainException('Unknown ID: '+str(id))
        else:
            prev_id=self.chain[id].prev
            if prev_id is not None:
                self.chain[prev_id].succ.pop(self.chain[prev_id].succ.index(id))
            self.chain[id].prev=None
            return prev_id

    def link(self, parent,child):
        """Connects son to parent: link(self,son,parent)."""
        if child not in self.chain:
            raise ChainException('Unknown ID: '+str(child))
        elif parent not in self.chain:
            raise ChainException('Unknown ID: '+str(parent))
        else:
            self.unlink(child)
            self.chain[parent].succ.append(child)
            self.chain[child].set_prev(parent)

    def is_parent_of(self,parent,grandchild):
        """Check if grandchild is a subnode of parent: is_parent_of(self,parent,grandchild)."""
        if grandchild==parent or grandchild in self.chain[parent].get_succ():
            return True
        else:
            for sn in self.chain[parent].get_succ():
                if self.is_parent_of(sn,grandchild):
                    return True
            else:
                return False

    def trace(self,start,finish):
        """Returns a list of all node_ids between two nodes (excluding start, including end): trace(start,end)."""
        if start not in self.chain or finish not in self.chain:
            raise NodeException('Unknown node.')
        if not self.is_parent_of(start,finish) or start==finish:
            return []
        for sn in self.chain[start].get_succ():
            if self.is_parent_of(sn,finish):
                return [sn]+self.trace(sn,finish)
                
class Node(object):
    """A single node."""

    def __init__(self,data=None):
        """Represents a node with one predecessor and multiple successors: (self, data=None)."""
        self.id=None
        self.data=data
        self.prev=None
        self.succ=[]

    def set_id(self,id):
        """Sets the id of a node, if not set yet: (self,id)."""
        if self.id is not None:
            raise NodeException('Node id cannot be changed.')
        self.id=id

    def get_id(self):
        """Returns the node's id: (self)."""
        return self.id

    def get_succ(self):
        """Returns a list of the node's successors: (self)."""
        return self.succ

    def get_prev(self):
        """Returns the id of the node's predecessor: (self)."""
        return self.prev

    def add_succ(self,id):
        """Adds a node id to the node's successors: (self,id)."""
        if isinstance(id,type([])):
            self.succ.extend(id)
        else:
            self.succ.append(id)

    def remove_succ(self,id):
        """Removes a node id from the node's successors: (self,id)."""
        self.succ.remove(id)

    def set_succ(self,new_succ):
        """Sets the node's successors: (self,new_succ)."""
        if not isinstance(new_succ,type([])):
            raise NodeException('Node successor must be of list type.')
        self.succ=new_succ

    def set_prev(self,id):
        """Sets the node's predecessor: (self,id)."""
        self.prev=id
    
    def get_data(self):
        """Returns a node's data: (self)."""
        return self.data

    def set_data(self,data):
        """Sets a node's data: (self,data)."""
        self.data=data
