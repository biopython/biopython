# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Base classes for Bio.Phylo objects.

All object representations for phylogenetic trees should derive from these base
classes in order to use the common methods defined on them.
"""
__docformat__ = "epytext en"

import itertools
import random
import re
from collections import deque


def trim_str(text, maxlen=60):
    assert isinstance(text, basestring), \
            "%s should be a string, not a %s" % (text, type(text))
    if len(text) > maxlen:
        return text[:maxlen-3] + '...'
    return text

# General tree-traversal algorithms

def _level_search(root, get_children):
    """Traverse a tree in breadth-first (level) order."""
    Q = deque([root])
    while Q:
        v = Q.popleft()
        yield v
        Q.extend(get_children(v))

def _preorder_search(root, get_children):
    """Traverse a tree in depth-first pre-order (parent before children)."""
    def dfs(elem):
        yield elem
        for v in get_children(elem):
            for u in dfs(v):
                yield u
    for elem in dfs(root):
        yield elem

def _postorder_search(root, get_children):
    """Traverse a tree in depth-first post-order (children before parent)."""
    def dfs(elem):
        for v in get_children(elem):
            for u in dfs(v):
                yield u
        yield elem
    for elem in dfs(root):
        yield elem

def _sorted_attrs(elem):
    """Get a flat list of elem's attributes, sorted for consistency."""
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

# Factory functions to generalize searching for clades/nodes

def _identity_matcher(target):
    """Match each node (or its root) to the target object by identity."""
    def match(node):
        return (node is target)
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
    """Base class for all Bio.Phylo classes."""

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
        if hasattr(self, 'name') and self.name:
            return self.name
        return self.__class__.__name__


class TreeMixin(object):
    """Methods for Tree- and Subtree-based classes.

    This lets Tree and Subtree support the same traversal and searching
    operations without requiring Subtree to inherit from Tree, so Subtree isn't
    required to have all of Tree's attributes -- just 'root' (a Subtree
    instance) and 'is_terminal()'.
    """
    # Traversal methods

    def _filter_search(self, filter_func, order, follow_attrs):
        """Perform a BFS or DFS traversal through all elements in this tree.

        @return: generator of all elements for which 'filter_func' is True.
        """
        order_opts = {'preorder': _preorder_search,
                      'postorder': _postorder_search,
                      'level': _level_search}
        try:
            order_func = order_opts[order]
        except KeyError:
            raise ValueError("Invalid order '%s'; must be one of: %s"
                             % (order, tuple(order_opts.keys())))
        if follow_attrs:
            get_children = _sorted_attrs
            root = self
        else:
            get_children = lambda elem: elem.clades
            root = self.root
        return itertools.ifilter(filter_func, order_func(root, get_children))

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

    def find_all(self, cls=TreeElement, terminal=None, order='preorder',
            **kwargs):
        """Find all tree elements matching the given attributes.

        @param cls: 
            Specifies the class of the object to search for. Objects that
            inherit from this type will also match. (The default, TreeElement,
            matches any standard Bio.Phylo type.)
        @type cls: type

        @param terminal:
            A boolean value to select for or against terminal nodes (a.k.a. leaf
            nodes). True searches for only terminal nodes, False excludes
            terminal nodes, and the default, None, searches both terminal and
            non-terminal nodes, as well as any tree elements lacking the
            'is_terminal' method.
        @type terminal: bool

        @param order:
            Tree traversal order: 'preorder' (default) is depth-first search,
            'postorder' is DFS with child nodes preceding parents, and 'level'
            is breadth-first search.
        @type order: string ('preorder'|'postorder'|'level')

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

            >>> from Bio.Phylo.IO import PhyloXMIO
            >>> phx = PhyloXMLIO.read('phyloxml_examples.xml')
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
        return self._filter_search(is_matching_elem, order, follow_attrs=True)

    def find_clades(self, cls=TreeElement, terminal=None, order='preorder',
            **kwargs):
        """Find each clade containing a matching element.

        That is, find each element as with find_all(), but return the
        corresponding clade object.
        """
        def match_attrs(elem):
            orig_clades = elem.__dict__.pop('clades')
            found = elem.find(cls=cls, **kwargs)
            elem.clades = orig_clades
            return (found is not None)
        if terminal is None:
            is_matching_elem = match_attrs
        else:
            def is_matching_elem(elem):
                return ((elem.is_terminal() == terminal)
                        and match_attrs(elem))
        return self._filter_search(is_matching_elem, order, follow_attrs=False)

    def get_path(self, target):
        """List the clades directly between the root and the given target.

        Returns a list of all clade objects along this path, ending with
        the given target, but excluding the root clade.
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
        if not check_in_path(self.root):
            return None
        return list(reversed(path))[1:]

    def get_terminals(self, order='preorder'):
        """Get a list of all of this tree's terminal (leaf) nodes."""
        return list(self.find_clades(terminal=True, order=order))

    def trace(self, start, finish):
        """List of all clade object between two targets in this tree.

        Excluding start, including finish.
        """
        mrca = self.common_ancestor(start, finish)
        fromstart = mrca.get_path(start)[-2::-1]
        to = mrca.get_path(finish)
        return fromstart + [mrca] + to

    # Information methods

    def common_ancestor(self, *targets):
        """Most recent common ancestor (clade) of all the given targets.

        Edge cases: 

            - If no target is given, returns self.root
            - If 1 target is given, returns the target
            - If any target is not found in this tree, raises an AssertionError
        """
        paths = [self.get_path(t) for t in targets]
        # Validation -- otherwise izip throws a spooky error below
        for p, t in zip(paths, targets):
            assert p is not None, "target %s is not in this tree" % repr(t)
        mrca = self.root
        for level in itertools.izip(*paths):
            ref = level[0]
            for other in level[1:]:
                if ref is not other:
                    break
            else:
                mrca = ref
            if ref is not mrca:
                break
        return mrca

    def count_terminals(self):
        """Counts the number of terminal (leaf) nodes within this tree."""
        counter = 0
        for i, leaf in enumerate(self.find_clades(terminal=True)):
            counter = i
        return counter + 1

    def is_bifurcating(self, node=None):
        """Return True if tree downstream of node is strictly bifurcating."""
        if node is not None:
            warnings.warn("use node.is_bifurcating() directly instead",
                          DeprecationWarning, stacklevel=2)
            return node.is_bifurcating()
        # Root can be trifurcating, because it has no ancestor
        if isinstance(self, BaseTree.Tree) and len(self.root) == 3:
            return (self.clade[0].is_bifurcating()
                    and self.clade[1].is_bifurcating()
                    and self.clade[2].is_bifurcating())
        if len(self.root) == 2:
            return (self.clade[0].is_bifurcating()
                    and self.clade[1].is_bifurcating())
        if len(self.root) == 0:
            return True
        return False

    def depths(self):
        """Create a mapping of tree clades to depths (by branch length).

        @return: dict of {clade: depth}
        """
        depths = {}
        def update_depths(node, curr_depth):
            depths[node] = curr_depth
            for child in node.clades:
                new_depth = curr_depth + (child.branch_length or 0)
                update_depths(child, new_depth)
        update_depths(self.root, 0)
        return depths

    def distance(self, target1, target2=None):
        """Calculate the sum of the branch lengths between two targets.

        If only one target is specified, the other is the root of this tree.
        """
        if target2 is None:
            return sum(n.branch_length for n in self.get_path(target1)
                       if n.branch_length is not None)
        mrca = self.common_ancestor(target1, target2)
        return mrca.distance(target1) + mrca.distance(target2)

    # TODO - unit test
    def is_bifurcating(self):
        """Return True if tree downstream of node is strictly bifurcating."""
        # Root can be trifurcating, because it has no ancestor
        if isinstance(self, BaseTree.Tree) and len(self.root) == 3:
            return (self.clade[0].is_bifurcating()
                    and self.clade[1].is_bifurcating()
                    and self.clade[2].is_bifurcating())
        if len(self.root) == 2:
            return (self.clade[0].is_bifurcating()
                    and self.clade[1].is_bifurcating())
        if len(self.root) == 0:
            return True
        return False

    # TODO - unit test
    def is_monophyletic(self, terminals):
        """Does taxon_list comprise a complete subclade of this clade?

        @return: common ancestor if subclades is monophyletic, otherwise False.
        """
        target_set = set(terminals)
        current = self.root
        while True:
            if set(current.get_terminals()) == target_set:
                return current
            # Try a narrower subclade
            for subclade in current.clades:
                if set(subclade.get_terminals()).issuperset(target_set):
                    current = subclade
                    break
                else:
                    return False

    def is_parent_of(self, target):
        """True if target is a descendent of this tree.

        Not required to be a direct descendent.
        """
        return (self.get_path(target) is not None)

    def is_preterminal(self):
        """True if all direct descendents are terminal."""
        if self.root.is_terminal():
            return False
        for clade in self.clades:
            if not clade.is_terminal():
                return False
        return True

    def total_branch_length(self):
        """Calculate the sum of all the branch lengths in this tree."""
        return sum(node.branch_length
                   for node in self.find_all(branch_length=True))

    # Tree manipulation methods

    def collapse(self, target):
        """Deletes target from the tree, relinking its children to its parent.

        @return: the parent clade.
        """
        path = list(self.get_path(target))
        if not path:
            raise ValueError("couldn't collapse %s in this tree" % target)
        if len(path) == 1:
            parent = self.root
        else:
            parent = path[-2]
        popped = parent.clades.pop(parent.clades.index(target))
        extra_length = popped.branch_length or 0
        for child in popped:
            child.branch_length += extra_length
        parent.clades.extend(popped.clades)
        return parent

    def ladderize(self, reverse=False):
        """Sort clades in-place according to the number of terminal nodes.

        Deepest clades are last by default. Use reverse=True to sort clades
        deepest-to-shallowest.
        """
        self.root.clades.sort(key=lambda c: c.count_terminals(),
                              reverse=reverse)
        for subclade in self.root.clades:
            subclade.ladderize(reverse=reverse)
        return

    # TODO - unit test
    def prune(self, target):
        """Prunes a terminal clade from the tree.

        If taxon is from a bifurcation, the connecting node will be collapsed
        and its branch length added to remaining terminal node. This might be no
        longer a meaningful value.

        @return parent clade of the pruned target
        """
        path = self.get_path(target)
        if not path:
            raise ValueError("can't prune the root clade")
        if not path[-1].is_terminal():
            raise ValueError("target must be terminal")
        if len(path) == 1:
            parent = self.root
        else:
            parent = path[-2]
        parent.clades.remove(path[-1])
        if len(parent) == 1:
            # We deleted one branch from a bifurcation
            if parent == self.root:
                # If we're at the root, move the root upwards
                newroot = parent.clades[0]
                newroot.branch_length = None
                parent = self.root = newroot
            else:
                # If we're not at the root, collapse this parent
                child = parent.clades[0]
                child.branch_length += (parent.branch_length or 0.0)
                if len(path) < 3:
                    grandparent = self.root
                else:
                    grandparent = path[-3]
                grandparent.clades.remove(parent)
                grandparent.clades.append(child)
                parent = grandparent
        return parent

    # TODO - unit test
    def split(self, n=2, branch_length=1.0):
        """Speciation: generate n (default 2) new descendants.

        New clades have the given branch_length and the same name as this
        clade's root plus an integer suffix (counting from 0).
        """
        subtree_cls = type(self.root)
        base_name = self.root.name or ''
        for i in range(n):
            clade = subtree_cls(name=base_name+str(i),
                                branch_length=branch_length)
            self.root.clades.append(clade)
        return self.root.clades[-n:]


class Tree(TreeElement, TreeMixin):
    """A phylogenetic tree, containing global info for the phylogeny.

    The structure and node-specific data is accessible through the 'root'
    subtree attached to the Tree instance.

    @param root:
        The starting node of the tree. If the tree is rooted, this will usually
        be the root node.
    @type root: Subtree

    @param rooted:
        Whether or not the tree is rooted. By default, a tree is assumed to be
        rooted.
    @type rooted: bool

    @param id: The identifier of the tree, if there is one.
    @type id: str

    @param name: The name of the tree, in essence a label.
    @type name: str
    """
    def __init__(self, root=None, rooted=True, id=None, name=None):
        self.root = root or Subtree()
        self.rooted = rooted
        self.id = id
        self.name = name

    @classmethod
    def from_subtree(cls, subtree, **kwargs):
        """Create a new Tree object given a subtree.

        Keyword arguments are the usual Tree constructor parameters.
        """
        return cls(subtree, **kwargs)

    @classmethod
    def randomized(cls, taxa, branch_length=1.0, branch_stdev=None):
        """Create a randomized bifurcating tree given a list of taxa.

        @param taxa: Either an integer specifying the number of taxa to create
            (automatically named taxon#), or an iterable of taxon names, as
            strings.

        @return: a tree of the same type as this class.
        """
        if isinstance(taxa, int):
            taxa = ['taxon%s' % (i+1) for i in xrange(taxa)]
        elif hasattr(taxa, '__iter__'):
            taxa = list(taxa)
        else:
            raise TypeError("taxa argument must be integer (# taxa) or "
                            "iterable of taxon names.")
        rtree = cls()
        terminals = [rtree.root]
        while len(terminals) < len(taxa):
            newsplit = random.choice(terminals)
            newterms = newsplit.split(branch_length=branch_length)
            if branch_stdev:
                # Add some noise to the branch lengths
                for nt in newterms:
                    nt.branch_length = max(0,
                            random.gauss(branch_length, branch_stdev))
            terminals.remove(newsplit)
            terminals.extend(newterms)
        # Distribute taxon labels randomly
        random.shuffle(taxa)
        for node, name in zip(terminals, taxa):
            node.name = name
        return rtree

    @property
    def clade(self):
        """The first subtree in this tree (not itself)."""
        return self.root

    def is_terminal(self):
        """Returns True if the root of this tree is terminal."""
        return (not self.root.clades)


class Subtree(TreeElement, TreeMixin):
    """A recursively defined subtree.

    @param branch_length:
        The length of the branch leading to the root node of this subtree.
    @type branch_length: str

    @param label: The label of a node.
    @type label: str

    @param clades: Sub-trees rooted directly under this tree's root.
    @type clades: list
    """
    def __init__(self, branch_length=None, name=None, clades=None):
        self.clades = clades or []
        self.name = name
        self.branch_length = branch_length

    def is_terminal(self):
        """Returns True if this is a terminal (leaf) node."""
        return (not self.clades)

    # Properties may be overridden by subclasses

    # XXX kind of superfluous
    # @property
    # def label(self):
    #     return str(self)

    @property
    def root(self):
        """Allow TreeMixin methods to traverse subtrees properly."""
        return self

    # Sequence-type behavior methods

    def __getitem__(self, index):
        """Get subtrees by index (integer or slice)."""
        if isinstance(index, int) or isinstance(index, slice):
            return self.clades[index]
        ref = self
        for idx in index:
            ref = ref[idx]
        return ref

    def __iter__(self):
        """Iterate through this tree's direct subtrees (clades)."""
        return iter(self.clades)

    def __len__(self):
        """Number of subtrees directy under the root."""
        return len(self.clades)

