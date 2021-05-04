# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Base classes for Bio.Phylo objects.

All object representations for phylogenetic trees should derive from these base
classes in order to use the common methods defined on them.
"""

import collections
import copy
import itertools
import random
import re
import warnings

from Bio import BiopythonDeprecationWarning


# General tree-traversal algorithms


def _level_traverse(root, get_children):
    """Traverse a tree in breadth-first (level) order (PRIVATE)."""
    Q = collections.deque([root])
    while Q:
        v = Q.popleft()
        yield v
        Q.extend(get_children(v))


def _preorder_traverse(root, get_children):
    """Traverse a tree in depth-first pre-order (parent before children) (PRIVATE)."""

    def dfs(elem):
        yield elem
        for v in get_children(elem):
            yield from dfs(v)

    yield from dfs(root)


def _postorder_traverse(root, get_children):
    """Traverse a tree in depth-first post-order (children before parent) (PRIVATE)."""

    def dfs(elem):
        for v in get_children(elem):
            yield from dfs(v)
        yield elem

    yield from dfs(root)


def _sorted_attrs(elem):
    """Get a flat list of elem's attributes, sorted for consistency (PRIVATE)."""
    singles = []
    lists = []
    # Sort attributes for consistent results
    for attrname, child in sorted(elem.__dict__.items(), key=lambda kv: kv[0]):
        if child is None:
            continue
        if isinstance(child, list):
            lists.extend(child)
        else:
            singles.append(child)
    return (x for x in singles + lists if isinstance(x, TreeElement))


# Factory functions to generalize searching for clades/nodes


def _identity_matcher(target):
    """Match a node to the target object by identity (PRIVATE)."""

    def match(node):
        return node is target

    return match


def _class_matcher(target_cls):
    """Match a node if it's an instance of the given class (PRIVATE)."""

    def match(node):
        return isinstance(node, target_cls)

    return match


def _string_matcher(target):
    def match(node):
        if isinstance(node, (Clade, Tree)):
            # Avoid triggering specialized or recursive magic methods
            return node.name == target
        return str(node) == target

    return match


def _attribute_matcher(kwargs):
    """Match a node by specified attribute values (PRIVATE).

    ``terminal`` is a special case: True restricts the search to external (leaf)
    nodes, False restricts to internal nodes, and None allows all tree elements
    to be searched, including phyloXML annotations.

    Otherwise, for a tree element to match the specification (i.e. for the
    function produced by ``_attribute_matcher`` to return True when given a tree
    element), it must have each of the attributes specified by the keys and
    match each of the corresponding values -- think 'and', not 'or', for
    multiple keys.
    """

    def match(node):
        if "terminal" in kwargs:
            # Special case: restrict to internal/external/any nodes
            kwa_copy = kwargs.copy()
            pattern = kwa_copy.pop("terminal")
            if pattern is not None and (
                not hasattr(node, "is_terminal") or node.is_terminal() != pattern
            ):
                return False
        else:
            kwa_copy = kwargs
        for key, pattern in kwa_copy.items():
            # Nodes must match all other specified attributes
            if not hasattr(node, key):
                return False
            target = getattr(node, key)
            if isinstance(pattern, str):
                return isinstance(target, str) and re.match(pattern + "$", target)
            if isinstance(pattern, bool):
                return pattern == bool(target)
            if isinstance(pattern, int):
                return pattern == target
            if pattern is None:
                return target is None
            raise TypeError("invalid query type: %s" % type(pattern))
        return True

    return match


def _function_matcher(matcher_func):
    """Safer attribute lookup -- returns False instead of raising an error (PRIVATE)."""

    def match(node):
        try:
            return matcher_func(node)
        except (LookupError, AttributeError, ValueError, TypeError):
            return False

    return match


def _object_matcher(obj):
    """Retrieve a matcher function by passing an arbitrary object (PRIVATE).

    Passing a ``TreeElement`` such as a ``Clade`` or ``Tree`` instance returns
    an identity matcher, passing a type such as the ``PhyloXML.Taxonomy`` class
    returns a class matcher, and passing a dictionary returns an attribute
    matcher.

    The resulting 'match' function returns True when given an object matching
    the specification (identity, type or attribute values), otherwise False.
    This is useful for writing functions that search the tree, and probably
    shouldn't be used directly by the end user.
    """
    if isinstance(obj, TreeElement):
        return _identity_matcher(obj)
    if isinstance(obj, type):
        return _class_matcher(obj)
    if isinstance(obj, str):
        return _string_matcher(obj)
    if isinstance(obj, dict):
        return _attribute_matcher(obj)
    if callable(obj):
        return _function_matcher(obj)
    raise ValueError(
        "%s (type %s) is not a valid type for comparison." % (obj, type(obj))
    )


def _combine_matchers(target, kwargs, require_spec):
    """Merge target specifications with keyword arguments (PRIVATE).

    Dispatch the components to the various matcher functions, then merge into a
    single boolean function.
    """
    if not target:
        if not kwargs:
            if require_spec:
                raise ValueError(
                    "you must specify a target object or keyword arguments."
                )
            return lambda x: True
        return _attribute_matcher(kwargs)
    match_obj = _object_matcher(target)
    if not kwargs:
        return match_obj
    match_kwargs = _attribute_matcher(kwargs)
    return lambda x: match_obj(x) and match_kwargs(x)


def _combine_args(first, *rest):
    """Convert ``[targets]`` or ``*targets`` arguments to a single iterable (PRIVATE).

    This helps other functions work like the built-in functions ``max`` and
    ``min``.
    """
    # Background: is_monophyletic takes a single list or iterable (like the
    # same method in Bio.Nexus.Trees); root_with_outgroup and common_ancestor
    # take separate arguments. This mismatch was in the initial release and I
    # didn't notice the inconsistency until after Biopython 1.55. I can think
    # of cases where either style is more convenient, so let's support both
    # (for backward compatibility and consistency between methods).
    if hasattr(first, "__iter__") and not (
        isinstance(first, TreeElement)
        or isinstance(first, type)
        or isinstance(first, str)
        or isinstance(first, dict)
    ):
        # terminals is an iterable of targets
        if rest:
            raise ValueError(
                "Arguments must be either a single list of "
                "targets, or separately specified targets "
                "(e.g. foo(t1, t2, t3)), but not both."
            )
        return first
    # terminals is a single target -- wrap in a container
    return itertools.chain([first], rest)


# Class definitions


class TreeElement:
    """Base class for all Bio.Phylo classes."""

    def __repr__(self):
        """Show this object's constructor with its primitive arguments."""

        def pair_as_kwarg_string(key, val):
            if isinstance(val, str):
                val = val[:57] + "..." if len(val) > 60 else val
                return "%s='%s'" % (key, val)
            return "%s=%s" % (key, val)

        return "%s(%s)" % (
            self.__class__.__name__,
            ", ".join(
                pair_as_kwarg_string(key, val)
                for key, val in sorted(self.__dict__.items())
                if val is not None and type(val) in (str, int, float, bool, str)
            ),
        )

    __str__ = __repr__


class TreeMixin:
    """Methods for Tree- and Clade-based classes.

    This lets ``Tree`` and ``Clade`` support the same traversal and searching
    operations without requiring Clade to inherit from Tree, so Clade isn't
    required to have all of Tree's attributes -- just ``root`` (a Clade
    instance) and ``is_terminal``.
    """

    # Traversal methods

    def _filter_search(self, filter_func, order, follow_attrs):
        """Perform a BFS or DFS traversal through all elements in this tree (PRIVATE).

        :returns: generator of all elements for which ``filter_func`` is True.

        """
        order_opts = {
            "preorder": _preorder_traverse,
            "postorder": _postorder_traverse,
            "level": _level_traverse,
        }
        try:
            order_func = order_opts[order]
        except KeyError:
            raise ValueError(
                "Invalid order '%s'; must be one of: %s" % (order, tuple(order_opts))
            ) from None

        if follow_attrs:
            get_children = _sorted_attrs
            root = self
        else:
            get_children = lambda elem: elem.clades  # noqa: E731
            root = self.root
        return filter(filter_func, order_func(root, get_children))

    def find_any(self, *args, **kwargs):
        """Return the first element found by find_elements(), or None.

        This is also useful for checking whether any matching element exists in
        the tree, and can be used in a conditional expression.
        """
        hits = self.find_elements(*args, **kwargs)
        try:
            return next(hits)
        except StopIteration:
            return None

    def find_elements(self, target=None, terminal=None, order="preorder", **kwargs):
        """Find all tree elements matching the given attributes.

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

        :Parameters:
            target : TreeElement instance, type, dict, or callable
                Specifies the characteristics to search for. (The default,
                TreeElement, matches any standard Bio.Phylo type.)
            terminal : bool
                A boolean value to select for or against terminal nodes (a.k.a.
                leaf nodes). True searches for only terminal nodes, False
                excludes terminal nodes, and the default, None, searches both
                terminal and non-terminal nodes, as well as any tree elements
                lacking the ``is_terminal`` method.
            order : {'preorder', 'postorder', 'level'}
                Tree traversal order: 'preorder' (default) is depth-first
                search, 'postorder' is DFS with child nodes preceding parents,
                and 'level' is breadth-first search.

        Examples
        --------
        >>> from Bio import Phylo
        >>> phx = Phylo.PhyloXMLIO.read('PhyloXML/phyloxml_examples.xml')
        >>> matches = phx.phylogenies[5].find_elements(code='OCTVU')
        >>> next(matches)
        Taxonomy(code='OCTVU', scientific_name='Octopus vulgaris')

        """
        if terminal is not None:
            kwargs["terminal"] = terminal
        is_matching_elem = _combine_matchers(target, kwargs, False)
        return self._filter_search(is_matching_elem, order, True)

    def find_clades(self, target=None, terminal=None, order="preorder", **kwargs):
        """Find each clade containing a matching element.

        That is, find each element as with find_elements(), but return the
        corresponding clade object. (This is usually what you want.)

        :returns: an iterable through all matching objects, searching
            depth-first (preorder) by default.

        """

        def match_attrs(elem):
            orig_clades = elem.__dict__.pop("clades")
            found = elem.find_any(target, **kwargs)
            elem.clades = orig_clades
            return found is not None

        if terminal is None:
            is_matching_elem = match_attrs
        else:

            def is_matching_elem(elem):
                return (elem.is_terminal() == terminal) and match_attrs(elem)

        return self._filter_search(is_matching_elem, order, False)

    def get_path(self, target=None, **kwargs):
        """List the clades directly between this root and the given target.

        :returns: list of all clade objects along this path, ending with the
            given target, but excluding the root clade.

        """
        # Only one path will work -- ignore weights and visits
        path = []
        match = _combine_matchers(target, kwargs, True)

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
        return path[-2::-1]

    def get_nonterminals(self, order="preorder"):
        """Get a list of all of this tree's nonterminal (internal) nodes."""
        return list(self.find_clades(terminal=False, order=order))

    def get_terminals(self, order="preorder"):
        """Get a list of all of this tree's terminal (leaf) nodes."""
        return list(self.find_clades(terminal=True, order=order))

    def trace(self, start, finish):
        """List of all clade object between two targets in this tree.

        Excluding ``start``, including ``finish``.
        """
        mrca = self.common_ancestor(start, finish)
        fromstart = mrca.get_path(start)[-2::-1]
        to = mrca.get_path(finish)
        return fromstart + [mrca] + to

    # Information methods

    def common_ancestor(self, targets, *more_targets):
        """Most recent common ancestor (clade) of all the given targets.

        Edge cases:
         - If no target is given, returns self.root
         - If 1 target is given, returns the target
         - If any target is not found in this tree, raises a ValueError

        """
        paths = [self.get_path(t) for t in _combine_args(targets, *more_targets)]
        # Validation -- otherwise izip throws a spooky error below
        for p, t in zip(paths, targets):
            if p is None:
                raise ValueError("target %r is not in this tree" % t)
        mrca = self.root
        for level in zip(*paths):
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
        """Count the number of terminal (leaf) nodes within this tree."""
        return sum(1 for clade in self.find_clades(terminal=True))

    def depths(self, unit_branch_lengths=False):  # noqa: D402
        """Create a mapping of tree clades to depths (by branch length).

        :Parameters:
            unit_branch_lengths : bool
                If True, count only the number of branches (levels in the tree).
                By default the distance is the cumulative branch length leading
                to the clade.

        :returns: dict of {clade: depth}, where keys are all of the Clade
            instances in the tree, and values are the distance from the root to
            each clade (including terminals).

        """  # noqa: D402
        if unit_branch_lengths:
            depth_of = lambda c: 1  # noqa: E731
        else:
            depth_of = lambda c: c.branch_length or 0  # noqa: E731
        depths = {}

        def update_depths(node, curr_depth):
            depths[node] = curr_depth
            for child in node.clades:
                new_depth = curr_depth + depth_of(child)
                update_depths(child, new_depth)

        update_depths(self.root, self.root.branch_length or 0)
        return depths

    def distance(self, target1, target2=None):
        """Calculate the sum of the branch lengths between two targets.

        If only one target is specified, the other is the root of this tree.
        """
        if target2 is None:
            return sum(
                n.branch_length
                for n in self.get_path(target1)
                if n.branch_length is not None
            )
        mrca = self.common_ancestor(target1, target2)
        return mrca.distance(target1) + mrca.distance(target2)

    def is_bifurcating(self):
        """Return True if tree downstream of node is strictly bifurcating.

        I.e., all nodes have either 2 or 0 children (internal or external,
        respectively). The root may have 3 descendents and still be considered
        part of a bifurcating tree, because it has no ancestor.
        """
        # Root can be trifurcating
        if isinstance(self, Tree) and len(self.root) == 3:
            return (
                self.root.clades[0].is_bifurcating()
                and self.root.clades[1].is_bifurcating()
                and self.root.clades[2].is_bifurcating()
            )
        if len(self.root) == 2:
            return (
                self.root.clades[0].is_bifurcating()
                and self.root.clades[1].is_bifurcating()
            )
        if len(self.root) == 0:
            return True
        return False

    def is_monophyletic(self, terminals, *more_terminals):
        """MRCA of terminals if they comprise a complete subclade, or False.

        I.e., there exists a clade such that its terminals are the same set as
        the given targets.

        The given targets must be terminals of the tree.

        To match both ``Bio.Nexus.Trees`` and the other multi-target methods in
        Bio.Phylo, arguments to this method can be specified either of two ways:
        (i) as a single list of targets, or (ii) separately specified targets,
        e.g. is_monophyletic(t1, t2, t3) -- but not both.

        For convenience, this method returns the common ancestor (MCRA) of the
        targets if they are monophyletic (instead of the value True), and False
        otherwise.

        :returns: common ancestor if terminals are monophyletic, otherwise False.

        """
        target_set = set(_combine_args(terminals, *more_terminals))
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

    def is_parent_of(self, target=None, **kwargs):
        """Check if target is a descendent of this tree.

        Not required to be a direct descendent.

        To check only direct descendents of a clade, simply use list membership
        testing: ``if subclade in clade: ...``
        """
        return self.get_path(target, **kwargs) is not None

    def is_preterminal(self):
        """Check if all direct descendents are terminal."""
        if self.root.is_terminal():
            return False
        for clade in self.root.clades:
            if not clade.is_terminal():
                return False
        return True

    def total_branch_length(self):
        """Calculate the sum of all the branch lengths in this tree."""
        return sum(node.branch_length for node in self.find_clades(branch_length=True))

    # Tree manipulation methods

    def collapse(self, target=None, **kwargs):
        """Delete target from the tree, relinking its children to its parent.

        :returns: the parent clade.

        """
        path = self.get_path(target, **kwargs)
        if not path:
            raise ValueError("couldn't collapse %s in this tree" % (target or kwargs))
        if len(path) == 1:
            parent = self.root
        else:
            parent = path[-2]
        popped = parent.clades.pop(parent.clades.index(path[-1]))
        extra_length = popped.branch_length or 0
        for child in popped:
            child.branch_length += extra_length
        parent.clades.extend(popped.clades)
        return parent

    def collapse_all(self, target=None, **kwargs):
        """Collapse all the descendents of this tree, leaving only terminals.

        Total branch lengths are preserved, i.e. the distance to each terminal
        stays the same.

        For example, this will safely collapse nodes with poor bootstrap
        support:

            >>> from Bio import Phylo
            >>> tree = Phylo.read('PhyloXML/apaf.xml', 'phyloxml')
            >>> print("Total branch length %0.2f" % tree.total_branch_length())
            Total branch length 20.44
            >>> tree.collapse_all(lambda c: c.confidence is not None and c.confidence < 70)
            >>> print("Total branch length %0.2f" % tree.total_branch_length())
            Total branch length 21.37

        This implementation avoids strange side-effects by using level-order
        traversal and testing all clade properties (versus the target
        specification) up front. In particular, if a clade meets the target
        specification in the original tree, it will be collapsed.  For example,
        if the condition is:

            >>> from Bio import Phylo
            >>> tree = Phylo.read('PhyloXML/apaf.xml', 'phyloxml')
            >>> print("Total branch length %0.2f" % tree.total_branch_length())
            Total branch length 20.44
            >>> tree.collapse_all(lambda c: c.branch_length < 0.1)
            >>> print("Total branch length %0.2f" % tree.total_branch_length())
            Total branch length 21.13

        Collapsing a clade's parent node adds the parent's branch length to the
        child, so during the execution of collapse_all, a clade's branch_length
        may increase. In this implementation, clades are collapsed according to
        their properties in the original tree, not the properties when tree
        traversal reaches the clade. (It's easier to debug.) If you want the
        other behavior (incremental testing), modifying the source code of this
        function is straightforward.
        """
        # Read the iterable into a list to protect against in-place changes
        matches = list(self.find_clades(target, False, "level", **kwargs))
        if not matches:
            # No matching nodes to collapse
            return
        # Skip the root node -- it can't be collapsed
        if matches[0] == self.root:
            matches.pop(0)
        for clade in matches:
            self.collapse(clade)

    def ladderize(self, reverse=False):
        """Sort clades in-place according to the number of terminal nodes.

        Deepest clades are last by default. Use ``reverse=True`` to sort clades
        deepest-to-shallowest.
        """
        self.root.clades.sort(key=lambda c: c.count_terminals(), reverse=reverse)
        for subclade in self.root.clades:
            subclade.ladderize(reverse=reverse)

    def prune(self, target=None, **kwargs):
        """Prunes a terminal clade from the tree.

        If taxon is from a bifurcation, the connecting node will be collapsed
        and its branch length added to remaining terminal node. This might be no
        longer be a meaningful value.

        :returns: parent clade of the pruned target

        """
        if "terminal" in kwargs and kwargs["terminal"]:
            raise ValueError("target must be terminal")
        path = self.get_path(target, terminal=True, **kwargs)
        if not path:
            raise ValueError("can't find a matching target below this root")
        if len(path) == 1:
            parent = self.root
        else:
            parent = path[-2]
        parent.clades.remove(path[-1])
        if len(parent) == 1:
            # We deleted a branch from a bifurcation
            if parent == self.root:
                # If we're at the root, move the root upwards
                # NB: This loses the length of the original branch
                newroot = parent.clades[0]
                newroot.branch_length = None
                parent = self.root = newroot
            else:
                # If we're not at the root, collapse this parent
                child = parent.clades[0]
                if child.branch_length is not None:
                    child.branch_length += parent.branch_length or 0.0
                if len(path) < 3:
                    grandparent = self.root
                else:
                    grandparent = path[-3]
                # Replace parent with child at the same place in grandparent
                index = grandparent.clades.index(parent)
                grandparent.clades.pop(index)
                grandparent.clades.insert(index, child)
                parent = grandparent
        return parent

    def split(self, n=2, branch_length=1.0):
        """Generate n (default 2) new descendants.

        In a species tree, this is a speciation event.

        New clades have the given branch_length and the same name as this
        clade's root plus an integer suffix (counting from 0). For example,
        splitting a clade named "A" produces sub-clades named "A0" and "A1".
        If the clade has no name, the prefix "n" is used for child nodes, e.g.
        "n0" and "n1".
        """
        clade_cls = type(self.root)
        base_name = self.root.name or "n"
        for i in range(n):
            clade = clade_cls(name=base_name + str(i), branch_length=branch_length)
            self.root.clades.append(clade)


class Tree(TreeElement, TreeMixin):
    """A phylogenetic tree, containing global info for the phylogeny.

    The structure and node-specific data is accessible through the 'root'
    clade attached to the Tree instance.

    :Parameters:
        root : Clade
            The starting node of the tree. If the tree is rooted, this will
            usually be the root node.
        rooted : bool
            Whether or not the tree is rooted. By default, a tree is assumed to
            be rooted.
        id : str
            The identifier of the tree, if there is one.
        name : str
            The name of the tree, in essence a label.

    """

    def __init__(self, root=None, rooted=True, id=None, name=None):
        """Initialize parameter for phylogenetic tree."""
        self.root = root or Clade()
        self.rooted = rooted
        self.id = id
        self.name = name

    @classmethod
    def from_clade(cls, clade, **kwargs):
        """Create a new Tree object given a clade.

        Keyword arguments are the usual ``Tree`` constructor parameters.
        """
        root = copy.deepcopy(clade)
        return cls(root, **kwargs)

    @classmethod
    def randomized(cls, taxa, branch_length=1.0, branch_stdev=None):
        """Create a randomized bifurcating tree given a list of taxa.

        :param taxa: Either an integer specifying the number of taxa to create
            (automatically named taxon#), or an iterable of taxon names, as
            strings.

        :returns: a tree of the same type as this class.

        """
        if isinstance(taxa, int):
            taxa = ["taxon%s" % (i + 1) for i in range(taxa)]
        elif hasattr(taxa, "__iter__"):
            taxa = list(taxa)
        else:
            raise TypeError(
                "taxa argument must be integer (# taxa) or iterable of taxon names."
            )
        rtree = cls()
        terminals = [rtree.root]
        while len(terminals) < len(taxa):
            newsplit = random.choice(terminals)
            newsplit.split(branch_length=branch_length)
            newterms = newsplit.clades
            if branch_stdev:
                # Add some noise to the branch lengths
                for nt in newterms:
                    nt.branch_length = max(0, random.gauss(branch_length, branch_stdev))
            terminals.remove(newsplit)
            terminals.extend(newterms)
        # Distribute taxon labels randomly
        random.shuffle(taxa)
        for node, name in zip(terminals, taxa):
            node.name = name
        return rtree

    @property
    def clade(self):
        """Return first clade in this tree (not itself)."""
        return self.root

    def as_phyloxml(self, **kwargs):
        """Convert this tree to a PhyloXML-compatible Phylogeny.

        This lets you use the additional annotation types PhyloXML defines, and
        save this information when you write this tree as 'phyloxml'.
        """
        from Bio.Phylo.PhyloXML import Phylogeny

        return Phylogeny.from_tree(self, **kwargs)

    def root_with_outgroup(
        self, outgroup_targets, *more_targets, outgroup_branch_length=None
    ):
        """Reroot this tree with the outgroup clade containing outgroup_targets.

        Operates in-place.

        Edge cases:
         - If ``outgroup == self.root``, no change
         - If outgroup is terminal, create new bifurcating root node with a
           0-length branch to the outgroup
         - If outgroup is internal, use the given outgroup node as the new
           trifurcating root, keeping branches the same
         - If the original root was bifurcating, drop it from the tree,
           preserving total branch lengths

        :param outgroup_branch_length: length of the branch leading to the
            outgroup after rerooting. If not specified (None), then:

            - If the outgroup is an internal node (not a single terminal taxon),
              then use that node as the new root.
            - Otherwise, create a new root node as the parent of the outgroup.

        """
        # This raises a ValueError if any target is not in this tree
        # Otherwise, common_ancestor guarantees outgroup is in this tree
        outgroup = self.common_ancestor(outgroup_targets, *more_targets)
        outgroup_path = self.get_path(outgroup)
        if len(outgroup_path) == 0:
            # Outgroup is the current root -- no change
            return

        prev_blen = outgroup.branch_length or 0.0

        if outgroup.is_terminal() or outgroup_branch_length is not None:
            # Create a new root with a 0-length branch to the outgroup
            outgroup.branch_length = outgroup_branch_length or 0.0
            new_root = self.root.__class__(
                branch_length=self.root.branch_length, clades=[outgroup]
            )
            # The first branch reversal (see the upcoming loop) is modified
            if len(outgroup_path) == 1:
                # No nodes between the original root and outgroup to rearrange.
                # Most of the code below will be skipped, but we still need
                # 'new_parent' pointing at the new root.
                new_parent = new_root
            else:
                parent = outgroup_path.pop(-2)
                # First iteration of reversing the path to the outgroup
                parent.clades.pop(parent.clades.index(outgroup))
                (prev_blen, parent.branch_length) = (
                    parent.branch_length,
                    prev_blen - outgroup.branch_length,
                )
                new_root.clades.insert(0, parent)
                new_parent = parent
        else:
            # Use the given outgroup node as the new (trifurcating) root
            new_root = outgroup
            new_root.branch_length = self.root.branch_length
            new_parent = new_root

        # Tracing the outgroup lineage backwards, reattach the subclades under a
        # new root clade. Reverse the branches directly above the outgroup in
        # the tree, but keep the descendants of those clades as they are.
        for parent in outgroup_path[-2::-1]:
            parent.clades.pop(parent.clades.index(new_parent))
            prev_blen, parent.branch_length = parent.branch_length, prev_blen
            new_parent.clades.insert(0, parent)
            new_parent = parent

        # Finally, handle the original root according to number of descendents
        old_root = self.root
        if outgroup in old_root.clades:
            assert len(outgroup_path) == 1
            old_root.clades.pop(old_root.clades.index(outgroup))
        else:
            old_root.clades.pop(old_root.clades.index(new_parent))
        if len(old_root) == 1:
            # Delete the old bifurcating root & add branch lengths
            ingroup = old_root.clades[0]
            if ingroup.branch_length:
                ingroup.branch_length += prev_blen
            else:
                ingroup.branch_length = prev_blen
            new_parent.clades.insert(0, ingroup)
            # ENH: If annotations are attached to old_root, do... something.
        else:
            # Keep the old trifurcating/polytomous root as an internal node
            old_root.branch_length = prev_blen
            new_parent.clades.insert(0, old_root)

        self.root = new_root
        self.rooted = True

    def root_at_midpoint(self):
        """Root the tree at the midpoint of the two most distant taxa.

        This operates in-place, leaving a bifurcating root. The topology of the
        tree is otherwise retained, though no guarantees are made about the
        stability of clade/node/taxon ordering.
        """
        # Identify the largest pairwise distance
        max_distance = 0.0
        tips = self.get_terminals()
        for tip in tips:
            self.root_with_outgroup(tip)
            new_max = max(self.depths().items(), key=lambda nd: nd[1])
            if new_max[1] > max_distance:
                tip1 = tip
                tip2 = new_max[0]
                max_distance = new_max[1]
        self.root_with_outgroup(tip1)
        # Depth to go from the ingroup tip toward the outgroup tip
        root_remainder = 0.5 * (max_distance - (self.root.branch_length or 0))
        assert root_remainder >= 0
        # Identify the midpoint and reroot there.
        # Trace the path to the outgroup tip until all of the root depth has
        # been traveled/accounted for.
        for node in self.get_path(tip2):
            root_remainder -= node.branch_length
            if root_remainder < 0:
                outgroup_node = node
                outgroup_branch_length = -root_remainder
                break
        else:
            raise ValueError("Somehow, failed to find the midpoint!")
        self.root_with_outgroup(
            outgroup_node, outgroup_branch_length=outgroup_branch_length
        )

    # Method assumed by TreeMixin

    def is_terminal(self):
        """Check if the root of this tree is terminal."""
        return not self.root.clades

    # Convention from SeqRecord and Alignment classes

    def __format__(self, format_spec):
        """Serialize the tree as a string in the specified file format.

        This method supports Python's ``format`` built-in function.

        :param format_spec: a lower-case string supported by ``Bio.Phylo.write``
            as an output file format.

        """
        if format_spec:
            from io import StringIO
            from Bio.Phylo import _io

            handle = StringIO()
            _io.write([self], handle, format_spec)
            return handle.getvalue()
        else:
            # Follow python convention and default to using __str__
            return str(self)

    def format(self, fmt=None, format=None):
        """Serialize the tree as a string in the specified file format.

        :param fmt: a lower-case string supported by ``Bio.Phylo.write``
            as an output file format.

        """
        if format is not None:
            if fmt is not None:
                raise ValueError("The ``format`` argument has been renamed to ``fmt``.")
            warnings.warn(
                "The ``format`` argument has been renamed to ``fmt``.",
                BiopythonDeprecationWarning,
            )
            fmt = format
        return self.__format__(fmt)

    # Pretty-printer for the entire tree hierarchy

    def __str__(self):
        """Return a string representation of the entire tree.

        Serialize each sub-clade recursively using ``repr`` to create a summary
        of the object structure.
        """
        TAB = "    "
        textlines = []

        def print_tree(obj, indent):
            """Recursively serialize sub-elements.

            This closes over textlines and modifies it in-place.
            """
            if isinstance(obj, (Tree, Clade)):
                # Avoid infinite recursion or special formatting from str()
                objstr = repr(obj)
            else:
                objstr = str(obj)
            textlines.append(TAB * indent + objstr)
            indent += 1
            for attr in obj.__dict__:
                child = getattr(obj, attr)
                if isinstance(child, TreeElement):
                    print_tree(child, indent)
                elif isinstance(child, list):
                    for elem in child:
                        if isinstance(elem, TreeElement):
                            print_tree(elem, indent)

        print_tree(self, 0)
        return "\n".join(textlines)


class Clade(TreeElement, TreeMixin):
    """A recursively defined sub-tree.

    :Parameters:
        branch_length : str
            The length of the branch leading to the root node of this clade.
        name : str
            The clade's name (a label).
        clades : list
            Sub-trees rooted directly under this tree's root.
        confidence : number
            Support.
        color : BranchColor
            The display color of the branch and descendents.
        width : number
            The display width of the branch and descendents.

    """

    def __init__(
        self,
        branch_length=None,
        name=None,
        clades=None,
        confidence=None,
        color=None,
        width=None,
    ):
        """Define parameters for the Clade tree."""
        self.branch_length = branch_length
        self.name = name
        self.clades = clades or []
        self.confidence = confidence
        self.color = color
        self.width = width

    @property
    def root(self):
        """Allow TreeMixin methods to traverse clades properly."""
        return self

    def is_terminal(self):
        """Check if this is a terminal (leaf) node."""
        return not self.clades

    # Sequence-type behavior methods

    def __getitem__(self, index):
        """Get clades by index (integer or slice)."""
        if isinstance(index, (int, slice)):
            return self.clades[index]
        ref = self
        for idx in index:
            ref = ref[idx]
        return ref

    def __iter__(self):
        """Iterate through this tree's direct descendent clades (sub-trees)."""
        return iter(self.clades)

    def __len__(self):
        """Return the number of clades directy under the root."""
        return len(self.clades)

    def __bool__(self):
        """Boolean value of an instance of this class (True).

        NB: If this method is not defined, but ``__len__``  is, then the object
        is considered true if the result of ``__len__()`` is nonzero. We want
        Clade instances to always be considered True.
        """
        return True

    def __str__(self):
        """Return name of the class instance."""
        if self.name:
            return self.name[:37] + "..." if len(self.name) > 40 else self.name
        return self.__class__.__name__

    # Syntax sugar for setting the branch color
    def _get_color(self):
        return self._color

    def _set_color(self, arg):
        if arg is None or isinstance(arg, BranchColor):
            self._color = arg
        elif isinstance(arg, str):
            if arg in BranchColor.color_names:
                # Known color name
                self._color = BranchColor.from_name(arg)
            elif arg.startswith("#") and len(arg) == 7:
                # HTML-style hex string
                self._color = BranchColor.from_hex(arg)
            else:
                raise ValueError("invalid color string %s" % arg)
        elif hasattr(arg, "__iter__") and len(arg) == 3:
            # RGB triplet
            self._color = BranchColor(*arg)
        else:
            raise ValueError("invalid color value %s" % arg)

    color = property(_get_color, _set_color, doc="Branch color.")


class BranchColor:
    """Indicates the color of a clade when rendered graphically.

    The color should be interpreted by client code (e.g. visualization
    programs) as applying to the whole clade, unless overwritten by the
    color(s) of sub-clades.

    Color values must be integers from 0 to 255.
    """

    color_names = {
        "red": (255, 0, 0),
        "r": (255, 0, 0),
        "yellow": (255, 255, 0),
        "y": (255, 255, 0),
        "green": (0, 128, 0),
        "g": (0, 128, 0),
        "cyan": (0, 255, 255),
        "c": (0, 255, 255),
        "blue": (0, 0, 255),
        "b": (0, 0, 255),
        "magenta": (255, 0, 255),
        "m": (255, 0, 255),
        "black": (0, 0, 0),
        "k": (0, 0, 0),
        "white": (255, 255, 255),
        "w": (255, 255, 255),
        # Names standardized in HTML/CSS spec
        # http://w3schools.com/html/html_colornames.asp
        "maroon": (128, 0, 0),
        "olive": (128, 128, 0),
        "lime": (0, 255, 0),
        "aqua": (0, 255, 255),
        "teal": (0, 128, 128),
        "navy": (0, 0, 128),
        "fuchsia": (255, 0, 255),
        "purple": (128, 0, 128),
        "silver": (192, 192, 192),
        "gray": (128, 128, 128),
        # More definitions from matplotlib/gcolor2
        "grey": (128, 128, 128),
        "pink": (255, 192, 203),
        "salmon": (250, 128, 114),
        "orange": (255, 165, 0),
        "gold": (255, 215, 0),
        "tan": (210, 180, 140),
        "brown": (165, 42, 42),
    }

    def __init__(self, red, green, blue):
        """Initialize BranchColor for a tree."""
        for color in (red, green, blue):
            assert (
                isinstance(color, int) and 0 <= color <= 255
            ), "Color values must be integers between 0 and 255."
        self.red = red
        self.green = green
        self.blue = blue

    @classmethod
    def from_hex(cls, hexstr):
        """Construct a BranchColor object from a hexadecimal string.

        The string format is the same style used in HTML and CSS, such as
        '#FF8000' for an RGB value of (255, 128, 0).
        """
        assert (
            isinstance(hexstr, str) and hexstr.startswith("#") and len(hexstr) == 7
        ), "need a 24-bit hexadecimal string, e.g. #000000"

        RGB = hexstr[1:3], hexstr[3:5], hexstr[5:]
        return cls(*[int("0x" + cc, base=16) for cc in RGB])

    @classmethod
    def from_name(cls, colorname):
        """Construct a BranchColor object by the color's name."""
        return cls(*cls.color_names[colorname])

    def to_hex(self):
        """Return a 24-bit hexadecimal RGB representation of this color.

        The returned string is suitable for use in HTML/CSS, as a color
        parameter in matplotlib, and perhaps other situations.

        Examples
        --------
        >>> bc = BranchColor(12, 200, 100)
        >>> bc.to_hex()
        '#0cc864'

        """
        return "#%02x%02x%02x" % (self.red, self.green, self.blue)

    def to_rgb(self):
        """Return a tuple of RGB values (0 to 255) representing this color.

        Examples
        --------
        >>> bc = BranchColor(255, 165, 0)
        >>> bc.to_rgb()
        (255, 165, 0)

        """
        return (self.red, self.green, self.blue)

    def __repr__(self):
        """Preserve the standard RGB order when representing this object."""
        return "%s(red=%d, green=%d, blue=%d)" % (
            self.__class__.__name__,
            self.red,
            self.green,
            self.blue,
        )

    def __str__(self):
        """Show the color's RGB values."""
        return "(%d, %d, %d)" % (self.red, self.green, self.blue)
