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


def trim_str(text, maxlen=60):
    assert isinstance(text, basestring), \
            "%s should be a string, not a %s" % (text, type(text))
    if len(text) > maxlen:
        return text[:maxlen-3] + '...'
    return text


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

    @param name:
        The name of the tree, in essence a label.

    @param identifier:
        The identifier of the tree, if there is one.

    @param is_rooted:
        Whether or not the tree is rooted. By default, a tree is assumed to be
        rooted.

    @param node_id:
        The starting node of the tree. If the tree is rooted, this will usually
        be the root node. Note that the root node(s) of a rooted tree must be
        stored in tree_root, too.

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
    def __init__(self, root=None, nodes=None, rooted=True, id=None, name=None):
        self.root = root or Node(tree=self)
        self.nodes = nodes or []    # each type Node, "under" root
        self.rooted = rooted        # is_rooted=True
        self.id = id                #: identifier
        self.name = name or self.root.label
        # relations: self, node (root_node), biodatabase
        # may belong to a sequence

    @classmethod
    def from_node(cls, node, **kwargs):
        """Create a new Tree object given a root node.

        Keyword arguments are the usual Tree constructor parameters.
        """
        return cls(node, **kwargs)

    # Plumbing
    def depth_first_search(self, node, filterfunc):
        """Perform a depth-first search through all nodes in this tree.
        
        @return: generator of all nodes for which 'filterfunc' is True.
        """
        singles = []
        lists = []
        # Sort attributes for consistent results
        for subnode in sorted(node.__dict__.itervalues()):
            if subnode is None:
                continue
            if isinstance(subnode, list):
                lists.extend(subnode)
            else:
                singles.append(subnode)
        for item in singles + lists:
            if isinstance(item, TreeElement):
                if filterfunc(item):
                    yield item
                for result in self.depth_first_search(item, filterfunc):
                    yield result

    def is_terminal(self):
        return (not self.nodes)

    def update_nested_set_index(self):
        """Calculate the left and right indexes for each node in the tree.

        See: Nested set representation.
        U{ http://www.oreillynet.com/pub/a/network/2002/11/27/bioconf.html }
        """
        # NB: resembles depth-first-search, but not quite
        # look at BioSQL's load_ncbi_taxonomy.pl
        pass

    # Porcelain
    def find(self, cls=None, terminal=None, **kwargs):
        """Find all sub-nodes matching the given attributes.

        The 'cls' argument specifies the class of the sub-node. Nodes that
        inherit from this type will also match. (The default, TreeElement,
        matches any standard Bio.Tree type.)

        The arbitrary keyword arguments indicate the attribute name of the
        sub-node and the value to match: string, integer or boolean. Strings are
        evaluated as regular expression matches; integers are compared directly
        for equality, and booleans evaluate the attribute's truth value (True or
        False) before comparing. To handle nonzero floats, search with a boolean
        argument, then filter the result manually.

        If no keyword arguments are given, then just the class type is used for
        matching.

        The result is an iterable through all matching objects, by depth-first
        search. (Not necessarily the same order as the elements appear in the
        source file!)

        Example:

            >>> phx = TreeIO.read('phyloxml_examples.xml', 'phyloxml')
            >>> matches = phx.phylogenies[5].find(code='OCTVU')
            >>> matches.next()
            Taxonomy(code='OCTVU', scientific_name='Octopus vulgaris')
        """ 
        if cls is None:
            cls = TreeElement

        def match_class(node):
            return isinstance(node, cls)

        def match_terminal(node):
            if hasattr(node, 'is_terminal') and not node.is_terminal():
                return (not terminal)
            return terminal

        def match_kwargs(node):
            for key, pattern in kwargs.iteritems():
                if not hasattr(node, key):
                    return False
                target = getattr(node, key)
                if (isinstance(pattern, basestring)
                        and isinstance(target, basestring)):
                    if not re.match(pattern, target):
                        return False
                elif isinstance(pattern, bool):
                    if pattern != bool(target):
                        return False
                elif isinstance(pattern, int):
                    if pattern != target:
                        return False
                else:
                    raise RuntimeError('invalid argument: ' + str(pattern))
            return True

        if terminal is None:
            if not kwargs :
                is_matching_node = match_class
            else:
                def is_matching_node(node):
                    return (match_class(node) and match_kwargs(node))
        elif not kwargs:
            def is_matching_node(node):
                return (match_class(node) and match_terminal(node))
        else:
            def is_matching_node(node):
                return (match_class(node)
                        and match_terminal(node)
                        and match_kwargs(node))

        return self.depth_first_search(self, is_matching_node)

    # From Bioperl's Bio::Tree::TreeI

    def total_branch_length(self):
        """Get the total length of this tree (sum of all branch lengths)."""
        raise NotImplementedError

    def get_leaf_nodes(self):
        """Request the taxa (leaves of the tree)."""
        raise NotImplementedError

    # From Bioperl's Bio::Tree::TreeFunctionsI

    # remove_node
    # get_lca (lowest common ancestor)
    # distance (between 2 nodes, specified however)
    # is_monophyletic
    # is_paraphyletic
    # reroot


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

    @param tree_id:
         The tree of which this node is a part of.

    @param left_idx:
         The left value of the nested set optimization structure for efficient
         hierarchical queries. Needs to be precomputed by a program, see J.
         Celko, SQL for Smarties.

    @param right_idx:
         The right value of the nested set optimization structure for efficient
         hierarchical queries. Needs to be precomputed by a program, see J.
         Celko, SQL for Smarties.
    """
    def __init__(self,
            tree=None, label=None,
            branch_length=None, left_idx=None, right_idx=None):
        self.tree = tree or Tree(root=self)     # tree_id, type Tree
        self.label = name
        self.branch_length = branch_length  # XXX or move this to Edge?
        self.left_idx = left_idx
        self.right_idx = right_idx
        # relations: self, tree
        # (id/identifier/label, parent/ancestor)
        # may belong to a sequence

    def is_terminal(self):
        return (not self.tree.nodes)

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

