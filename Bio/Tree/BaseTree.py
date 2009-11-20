# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Base classes for Bio.Tree objects.

All object representations for phylogenetic trees should derive from these base
classes in order to use the common methods defined on them.

The attributes on these classes map directly to the schema for PhyloDB, a BioSQL
extension for representing phylogenetic trees or networks. (Hopefully, this will
let us support database import/export through this module eventually.)

See: U{ http://biosql.org/wiki/Extensions }

"""
__docformat__ = "epytext en"

import re
from collections import deque
from itertools import izip


def trim_str(text, maxlen=60):
    assert isinstance(text, basestring), \
            "%s should be a string, not a %s" % (text, type(text))
    if len(text) > maxlen:
        return text[:maxlen-3] + '...'
    return text


# Factory functions to generalize searching for clades/nodes

def _identity_matcher(target):
    """Match each node (or its root) to the target object by identity."""
    def match(node):
        return (node is target or node is target.root)
    return match

def _class_matcher(target_cls):
    def match(node):
        return isinstance(node, target_cls)
    return match

def _attribute_matcher(**kwargs):
    """Match a node by specified attribute values.

    'terminal' is a special case: True restricts the search to external (leaf)
    nodes, False restricts to internal nodes, and None allows all tree elements
    to be searched, including phyloXML annotations.

    Otherwise, for a tree element to match the specification (i.e. for the
    function produced by _attribute_matcher to return True when given a tree
    element), it must have each of the attributes specified by the keys and
    match each of the corresponding values -- think 'and', not 'or', for
    multiple keys.
    """
    def match(node):
        for key, pattern in kwargs.iteritems():
            # Special case: restrict to internal/external/any nodes
            if key == 'terminal':
                if (pattern is None
                    or (hasattr(node, 'is_terminal')
                        and node.is_terminal() == pattern)):
                    continue
                return False
            # Nodes must match all other specified attributes
            if not hasattr(node, key):
                return False
            target = getattr(node, key)
            if isinstance(pattern, basestring):
                return (isinstance(target, basestring)
                        and re.match(pattern+'$', target))
            if isinstance(pattern, bool):
                return (pattern == bool(target))
            if isinstance(pattern, int):
                return (pattern == target)
            if pattern is None:
                return (target is None)
            raise TypeError('invalid query type: %s' % type(pattern))
        return True
    return match

def _function_matcher(matcher_func):
    """Safer attribute lookup -- returns False instead of raising"""
    def match(node):
        try:
            return matcher_func(node)
        except (LookupError, AttributeError, ValueError):
            return False
    return match

def _object_matcher(obj):
    """Retrieve a matcher function by passing an arbitrary object.

    i.e. passing a TreeElement such as a Node or Tree instance returns an
    identity matcher, passing a type such as the PhyloXML.Taxonomy class returns
    a class matcher, and passing a dictionary returns an attribute matcher.
    
    The resulting 'match' function returns True when given an object matching
    the specification (identity, type or attribute values), otherwise False.
    This is useful for writing functions that search the tree, and probably
    shouldn't be used directly by the end user.
    """
    if isinstance(obj, TreeElement):
        return _identity_matcher(obj)
    if isinstance(obj, type):
        return _class_matcher(obj)
    if isinstance(obj, dict):
        return _attribute_matcher(**obj)
    if callable(obj):
        return _function_matcher(obj)
    raise ValueError("%s (type %s) is not a valid type for comparison.")


# Class definitions

class TreeElement(object):
    """Base class for all Bio.Tree classes."""

    def __repr__(self):
        """Show this object's constructor with its primitive arguments."""
        s = '%s(%s)' % (self.__class__.__name__,
                        ', '.join("%s='%s'"
                                  % (key, trim_str(unicode(val)))
                            for key, val in self.__dict__.iteritems()
                            if val is not None
                            and type(val) in (str, int, float, bool, unicode)))
        return s.encode('utf-8')

    def __str__(self):
        return self.__class__.__name__


class Tree(TreeElement):
    """A phylogenetic tree.

    From PhyloDB:
    =============

        A tree basically is a namespace for nodes, and thereby implicitly for
        their relationships (edges). In this model, tree is also bit of misnomer
        because we try to support reticulating trees, i.e., networks, too, so
        arguably it should be called graph. Typically, this will be used for
        storing phylogenetic trees, sequence trees (a.k.a. gene trees) as much
        as species trees.

    @param clades: Sub-trees rooted directly under this tree's root.
    @type clades: list

    @param name:
        The name of the tree, in essence a label.
    @type name: str

    @param id:
        The identifier of the tree, if there is one.
    @type id: str

    @param rooted:
        Whether or not the tree is rooted. By default, a tree is assumed to be
        rooted.
    @type rooted: bool

    @param root:
        The starting node of the tree. If the tree is rooted, this will usually
        be the root node. Note that the root node(s) of a rooted tree must be
        stored in tree_root, too.
    @type root: list

    Not implemented:

    @param biodatabase_id:
        The namespace of the tree itself. Though trees are in a sense named
        containers themselves (namely for nodes), they also constitute (possibly
        identifiable!) data objects in their own right. Some data sources may
        only provide a single tree, so that assigning a namespace for the tree
        may seem excessive, but others, such as TreeBASE, contain many trees,
        just as sequence databanks contain many sequences. The choice of how to
        name a tree is up to the user; one may assign a default namespace (such
        as "biosql"), or create one named the same as the tree.

    """
    def __init__(self, root=None, clades=None, rooted=True, id=None, name=None):
        self._root = root or Node(tree=self)
        self.clades = clades or []  # each type Tree, "under" root
        self.rooted = rooted        # is_rooted=True
        self.id = id                #: identifier
        self.name = name or self.root.label

    # Properties may be overridden by subclasses

    def _get_nodes(self): return [c.root for c in self.clades]
    def _set_nodes(self, nodes): self.clades = [Tree(n) for n in nodes]
    def _del_nodes(self, x): self.clades = []
    nodes = property(_get_nodes, _set_nodes, _del_nodes)

    def _get_root(self): return self._root
    def _set_root(self, x): self._root = x
    def _del_root(self, x): self._root = Node(tree=self)
    root = property(_get_root, _set_root, _del_root)

    @classmethod
    def from_node(cls, node, **kwargs):
        """Create a new Tree object given a root node.

        Keyword arguments are the usual Tree constructor parameters.
        """
        return cls(node, **kwargs)

    # Plumbing

    def filter_search(self, filterfunc, breadth_first):
        """Perform a BFS or DFS through all elements in this tree.

        @return: generator of all elements for which 'filterfunc' is True.
        """
        def get_children(elem):
            singles = []
            lists = []
            # Sort attributes for consistent results
            for child in sorted(elem.__dict__.itervalues()):
                if child is None:
                    continue
                if isinstance(child, list):
                    lists.extend(child)
                else:
                    singles.append(child)
            return (x for x in singles + lists
                    if isinstance(x, TreeElement))

        Q = deque(get_children(self))
        if breadth_first:
            pop = deque.popleft
            extend = deque.extend
        else:
            pop = deque.pop
            extend = lambda q, s: deque.extend(q, reversed(tuple(s)))
        while Q:
            v = pop(Q)
            if filterfunc(v):
                yield v
            extend(Q, get_children(v))

    # TODO: write a unit test
    def get_path(self, target):
        """Find the direct path from the root to the given target.

        Returns an iterable of all clade origins along this path, ending with
        the given target.
        """
        # Only one path will work -- ignore weights and visits
        path = deque()
        match = _object_matcher(target)

        def check_in_path(v):
            if match(v):
                path.append(v)
                return True
            elif v.is_terminal():
                return False
            for child in v:
                if check_in_path(child):
                    path.append(v)
                    return True
                return False

        if not check_in_path(self):
            return None
        return reversed(path)

    def collapse(self, target):
        """Deletes target from chain and relinks successors to predecessor.

        Returns the predecessor clade.
        """
        path = list(self.get_path(target))
        if not path:
            raise ValueError("couldn't collapse %s in this tree" % target)
        if len(path) == 1:
            parent = self
        else:
            parent = path[-2]
        parent.clades.extend(
                parent.clades.pop(
                    parent.clades.index(target)
                    ).clades)
        return parent

    def common_ancestor(self, target1, target2):
        mrca = self
        for clade1, clade2 in izip(
                self.get_path(target1), 
                self.get_path(target2)): 
            if clade1 is clade2:
                mrca = clade1
            else:
                break
        return mrca
        # ENH: take arbitrary number of *targets
        # paths = [self.get_path(t) for t in targets]
        # mrca = self
        # for level in izip(paths):
        #     ref = level[0]
        #     for other in level[1:]:
        #         if ref is not other:
        #             break
        #     else:
        #         mrca = ref
        #     if ref is not mrca:
        #         break
        # return mrca

    # Porcelain

    def find(self, *args, **kwargs):
        """Return the first element found by find_all(), or None.

        This is also useful for checking whether any matching element exists in
        the tree.
        """
        hits = self.find_all(*args, **kwargs)
        try:
            return hits.next()
        except StopIteration:
            return None

    def find_all(self, cls=TreeElement, terminal=None, breadth_first=False,
            **kwargs):
        """Find all tree elements matching the given attributes.

        @param cls: 
            Specifies the class of the object to search for. Objects that
            inherit from this type will also match. (The default, TreeElement,
            matches any standard Bio.Tree type.)

        @param terminal:
            A boolean value to select for or against terminal nodes (a.k.a. leaf
            nodes). True searches for only terminal nodes, False excludes
            terminal nodes, and the default, None, searches both terminal and
            non-terminal nodes, as well as any tree elements lacking the
            'is_terminal' method.

        The arbitrary keyword arguments indicate the attribute name of the
        sub-element and the value to match: string, integer or boolean. Strings
        are evaluated as regular expression matches; integers are compared
        directly for equality, and booleans evaluate the attribute's truth value
        (True or False) before comparing. To handle nonzero floats, search with
        a boolean argument, then filter the result manually.

        If no keyword arguments are given, then just the class type is used for
        matching.

        The result is an iterable through all matching objects, by depth-first
        search. (Not necessarily the same order as the elements appear in the
        source file!)

        Example:

            >>> phx = TreeIO.read('phyloxml_examples.xml', 'phyloxml')
            >>> matches = phx.phylogenies[5].find_all(code='OCTVU')
            >>> matches.next()
            Taxonomy(code='OCTVU', scientific_name='Octopus vulgaris')

        """ 
        assert isinstance(cls, type), "cls argument must be a class or type"
        match_class = _class_matcher(cls)
        if terminal is not None:
            kwargs['terminal'] = terminal
        if kwargs:
            match_attr = _attribute_matcher(**kwargs)
            def is_matching_elem(elem):
                return (match_class(elem) and match_attr(elem))
        else:
            is_matching_elem = match_class

        return self.filter_search(is_matching_elem, breadth_first)

    def find_clades(self, cls=TreeElement, terminal=None, breadth_first=False,
            **kwargs):
        """Find each clade containing a matching element.

        That is, find each element as with find_all(), but return the
        corresponding clade object.
        """
        for clade in self.find_all(Tree,
                terminal=terminal, breadth_first=breadth_first):
            # Check whether any non-clade attributes/sub-elements match
            orig_clades = clade.__dict__.pop('clades')
            found = n.find(cls, **kwargs)
            clade.clades = orig_clades
            if found is not None:
                yield clade

    def get_terminals(self, breadth_first=False):
        """Iterate through all of this tree's terminal (leaf) nodes."""
        return self.find_all(Node, terminal=True, breadth_first=breadth_first)

    def is_terminal(self):
        """Returns True if the root of this tree is terminal."""
        return (not self.clades)

    def is_preterminal(self):
        """Returns True if all direct descendents are terminal."""
        if self.is_terminal():
            return False
        for clade in self.clades:
            if not clade.is_terminal():
                return False
        return True

    def count_terminals(self):
        """Counts the number of terminal (leaf) nodes within this tree."""
        counter = 0
        for i, leaf in enumerate(self.get_terminals()):
            counter = i
        return counter + 1

    def distance(self, target1, target2=None):
        """Calculate the sum of the branch lengths between two targets.

        If only one targets is specified, the other is the root of this tree.
        """
        if target2 is None:
            return sum(n.branch_length for n in self.get_path(target1)
                       if n.branch_length is not None)
        root = self.common_ancestor(target1, target2)
        return root.branch_length_to(target1) + root.branch_length_to(target2)

    def total_branch_length(self):
        """Calculate the sum of all the branch lengths in this tree."""
        return sum(node.branch_length
                   for node in self.find_all(branch_length=True))

    def trace(self, start, finish):
        """Returns a list of all tree elements between two targets.

        Excluding start, including end.
        """
        mrca = self.common_ancestor(start, finish)
        fromstart = list(mrca.get_path(start))[-2::-1]
        to = list(mrca.get_path(finish))
        return fromstart + [mrca] + to

    # Sequence-type behavior methods

    def __getitem__(self, index):
        """Get sub-trees by index (integer or slice)."""
        if isinstance(index, int) or isinstance(index, slice):
            return self.clades[index]
        ref = self
        for idx in index:
            ref = ref[idx]
        return ref

    def __iter__(self):
        """Iterate through this tree's direct sub-trees (clades)."""
        return iter(self.clades)

    def __len__(self):
        """Number of nodes/sub-trees directy under this tree's root."""
        return len(self.clades)


class Node(TreeElement):
    """A node in a tree.

    From PhyloDB:
    =============

        Typically, this will be a node in a phylogenetic tree, resembling either
        a nucleotide or protein sequence, or a taxon, or more generally an
        "operational taxonomic unit" (OTU).

    @param label:
         The label of a node. This may the latin binomial of the taxon, the
         accession number of a sequences, or any other construct that uniquely
         identifies the node within one tree.

    @param tree:
         The tree of which this node is a part of.
    """
    def __init__(self, tree=None, label=None, branch_length=None):
        self._tree = tree or Tree(root=self)     # tree_id, type Tree
        self._label = label
        self.branch_length = branch_length  # XXX or move this to Edge?

    # Properties may be overridden by subclasses

    def _get_label(self): return self._label
    def _set_label(self, x): self._label = x
    def _del_label(self, x): del self._label
    label = property(_get_label, _set_label, _del_label)

    def _get_tree(self): return self._tree
    def _set_tree(self, x): self._tree = x
    def _del_tree(self, x): self._tree = Tree(root=self)
    tree = property(_get_tree, _set_tree, _del_tree)

    def is_terminal(self):
        """Returns True if this is a terminal (leaf) node."""
        return (not self.tree.clades)


# Additional PhyloDB tables

class TreeQualifierValue(TreeElement):
    """Tree metadata as attribute/value pairs.

    From PhyloDB:
    =============

        Attribute names are from a controlled vocabulary (or ontology).

    @param tree_id:
	 The tree with which the metadata is being associated.

    @param term_id:
         The name of the metadate element as a term from a controlled vocabulary
         (or ontology).

    @param value:
         The value of the metadata element.

    @param rank:
         The index of the metadata value if there is more than one value for
         the same metadata element. If there is only one value, this may be left
         at the default of zero.
    """
    # value, rank
    # relations: tree, term
    pass


class TreeDbxref(TreeElement):
    """Secondary identifiers and other database cross-references for trees.

    From PhyloDB:
    =============

        There can only be one dbxref of a specific type for a tree.

    @param tree_id:
	The tree to which the database corss-reference is being assigned.

    @param dbxref_id:
	The database cross-reference being assigned to the tree.

    @param term_id:
        The type of the database cross-reference as a controlled vocabulary or
        ontology term. The type of a tree accession should be ''primary
        identifier''.
    """
    # relations: tree, term, dbxref
    pass


class NodeQualifierValue(TreeElement):
    """Tree (or network) node metadata as attribute/value pairs.

    From PhyloDB:
    =============

        Attribute names are from a controlled vocabulary (or ontology).

    @param node_id:
        The tree (or network) node to which the metadata is being associated.

    @param term_id:
        The name of the metadate element as a term from a controlled vocabulary
        (or ontology).

    @param value:
        The value of the attribute/value pair association of metadata (if
        applicable).

    @param rank:
        The index of the metadata value if there is more than one value for the
        same metadata element. If there is only one value, this may be left at
        the default of zero.

    """
    # value, rank
    # relations: node, term


class NodeDbxref(TreeElement):
    """Identifiers and other database cross-references for nodes. 

    From PhyloDB:
    =============

        There can only be one dbxref of a specific type for a node.

    @param node_id:
        The node to which the database cross-reference is being assigned.

    @param dbxref_id:
        The database cross-reference being assigned to the node.

    @param term_id:
        The type of the database cross-reference as a controlled vocabulary or
        ontology term. The type of a node identifier should be ''primary
        identifier''.
    """
    # relations: node, term, dbxref
    pass


class NodeBioentry(TreeElement):
    """Links tree nodes to sequences (or other bioentries).

    From PhyloDB:
    =============

        If the alignment is concatenated on molecular data, there will be more
        than one sequence, and rank can be used to order these appropriately.

    @param node_id:
        The node to which the bioentry is being linked.

    @param bioentry_id:
        The bioentry being linked to the node.

    @param rank:
        The index of this bioentry within the list of bioentries being linked to
        the node, if the order is significant. Typically, this will be used to
        represent the position of the respective sequence within the
        concatenated alignment, or the partition index.
    """
    # rank
    # relations: self, node, bioentry
    # may belong to a sequence
    pass


class NodeTaxon(TreeElement):
    """Links tree nodes to taxa.

    From PhyloDB:
    =============

        If the alignment is concatenated on molecular data, there may be more
        than one sequence, and these may not necessarily be from the same taxon
        (e.g., they might be from subspecies). Rank can be used to order these
        appropriately.';

    @param node_id:
	The node to which the taxon is being linked.

    @param taxon_id:
        The taxon being linked to the node.

    @param rank:
        The index of this taxon within the list of taxa being linked to the
        node, if the order is significant. Typically, this will be used to
        represent the position of the respective sequence within the
        concatenated alignment, or the partition index.
    """
    # rank
    # relations: self, node, taxon
    # may belong to a sequence
    pass


class TreeRoot(TreeElement):
    """Root node for a rooted tree.

    From PhyloDB:
    =============

        A phylogenetic analysis might suggest several alternative root nodes,
        with possible probabilities.

    @param tree_id:
        The tree for which the referenced node is a root node.

    @param node_id:
        The node that is a root for the referenced tree.

    @param is_alternate:
        True if the root node is the preferential (most likely) root node of the
        tree, and false otherwise.

    @param significance:
        The significance (such as likelihood, or posterior probability) with
        which the node is the root node. This only has meaning if the method
        used for reconstructing the tree calculates this value.
    """
    # is_alternate, significance
    # relations: self, tree, node
    # may belong to a sequence
    pass


class Edge(TreeElement):
    """An edge between two nodes in a tree (or graph).

    From PhyloDB:
    =============

    @param child_node_id:
        The endpoint node of the two nodes connected by a directed edge. In a
        phylogenetic tree, this is the descendant.

    @param parent_node_id:
        The startpoint node of the two nodes connected by a directed edge. In a
        phylogenetic tree, this is the ancestor.
    """
    # relations: self, child_node, parent_node
    # (distance/length, qualifier)
    # may belong to a sequence
    pass


class NodePath(TreeElement):
    """A path between two nodes in a tree (or graph).

    From PhyloDB:
    =============

        Two nodes A and B are connected by a (directed) path if B can be reached
        from A by following nodes that are connected by (directed) edges.

    @param child_node_id:
        The endpoint node of the two nodes connected by a (directed) path. In a
        phylogenetic tree, this is the descendant.

    @param parent_node_id:
        The startpoint node of the two nodes connected by a (directed) path. In
        a phylogenetic tree, this is the ancestor.

    @param path:
        The path from startpoint to endpoint as the series of nodes visited
        along the path. The nodes may be identified by label, or, typically more
        efficient, by their primary key, or left or right value. The latter or
        often smaller than the primary key, and hence consume less space. One
        may increase efficiency further by using a base-34 numeric
        representation (24 letters of the alphabet, plus 10 digits) instead of
        decimal (base-10) representation. The actual method used is not
        important, though it should be used consistently.

    @param distance:
        The distance (or length) of the path. The path between a node and itself
        has length zero, and length 1 between two nodes directly connected by an
        edge. If there is a path of length l between two nodes A and Z and an
        edge between Z and B, there is a path of length l+1 between nodes A and
        B.
    """
    # path, distance
    # relations: child_node, parent_node
    pass


class EdgeQualifierValue(TreeElement):
    """Edge metadata as attribute/value pairs.

    From PhyloDB:
    =============

        Attribute names are from a controlled vocabulary (or ontology).

    @param edge_id:
	The tree edge to which the metadata is being associated.

    @param term_id:
        The name of the metadate element as a term from a controlled vocabulary
        (or ontology).

    @param value:
        The value of the attribute/value pair association of metadata (if
        applicable).

    @param rank:
        The index of the metadata value if there is more than one value for the
        same metadata element. If there is only one value, this may be left at
        the default of zero.
    """
    # value, rank
    # relations: edge, term
    pass

