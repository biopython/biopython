# Copyright (C) 2013 by Yanbo Ye (yeyanbo289@gmail.com)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Classes and methods for finding consensus trees.

This module contains a ``_BitString`` class to assist the consensus tree
searching and some common consensus algorithms such as strict, majority rule and
adam consensus.
"""

import random
import itertools

from ast import literal_eval
from Bio.Phylo import BaseTree
from Bio.Align import MultipleSeqAlignment


class _BitString(str):
    """Helper class for binary string data (PRIVATE).

    Assistant class of binary string data used for storing and
    counting compatible clades in consensus tree searching. It includes
    some binary manipulation(&|^~) methods.

    _BitString is a sub-class of ``str`` object that only accepts two
    characters('0' and '1'), with additional functions for binary-like
    manipulation(&|^~). It is used to count and store the clades in
    multiple trees in consensus tree searching. During counting, the
    clades will be considered the same if their terminals(in terms of
    ``name`` attribute) are the same.

    For example, let's say two trees are provided as below to search
    their strict consensus tree::

        tree1: (((A, B), C),(D, E))
        tree2: ((A, (B, C)),(D, E))

    For both trees, a _BitString object '11111' will represent their
    root clade. Each '1' stands for the terminal clade in the list
    [A, B, C, D, E](the order might not be the same, it's determined
    by the ``get_terminal`` method of the first tree provided). For
    the clade ((A, B), C) in tree1 and (A, (B, C)) in tree2, they both
    can be represented by '11100'. Similarly, '11000' represents clade
    (A, B) in tree1, '01100' represents clade (B, C) in tree2, and '00011'
    represents clade (D, E) in both trees.

    So, with the ``_count_clades`` function in this module, finally we
    can get the clade counts and their _BitString representation as follows
    (the root and terminals are omitted)::

        clade   _BitString   count
        ABC     '11100'     2
        DE      '00011'     2
        AB      '11000'     1
        BC      '01100'     1

    To get the _BitString representation of a clade, we can use the following
    code snippet::

        # suppose we are provided with a tree list, the first thing to do is
        # to get all the terminal names in the first tree
        term_names = [term.name for term in trees[0].get_terminals()]
        # for a specific clade in any of the tree, also get its terminal names
        clade_term_names = [term.name for term in clade.get_terminals()]
        # then create a boolean list
        boolvals = [name in clade_term_names for name in term_names]
        # create the string version and pass it to _BitString
        bitstr = _BitString(''.join(map(str, map(int, boolvals))))
        # or, equivalently:
        bitstr = _BitString.from_bool(boolvals)

    To convert back::

        # get all the terminal clades of the first tree
        terms = [term for term in trees[0].get_terminals()]
        # get the index of terminal clades in bitstr
        index_list = bitstr.index_one()
        # get all terminal clades by index
        clade_terms = [terms[i] for i in index_list]
        # create a new calde and append all the terminal clades
        new_clade = BaseTree.Clade()
        new_clade.clades.extend(clade_terms)

    Examples
    --------
    >>> from Bio.Phylo.Consensus import _BitString
    >>> bitstr1 = _BitString('11111')
    >>> bitstr2 = _BitString('11100')
    >>> bitstr3 = _BitString('01101')
    >>> bitstr1
    _BitString('11111')
    >>> bitstr2 & bitstr3
    _BitString('01100')
    >>> bitstr2 | bitstr3
    _BitString('11101')
    >>> bitstr2 ^ bitstr3
    _BitString('10001')
    >>> bitstr2.index_one()
    [0, 1, 2]
    >>> bitstr3.index_one()
    [1, 2, 4]
    >>> bitstr3.index_zero()
    [0, 3]
    >>> bitstr1.contains(bitstr2)
    True
    >>> bitstr2.contains(bitstr3)
    False
    >>> bitstr2.independent(bitstr3)
    False
    >>> bitstr1.iscompatible(bitstr2)
    True
    >>> bitstr2.iscompatible(bitstr3)
    False

    """

    def __new__(cls, strdata):
        """Init from a binary string data."""
        if isinstance(strdata, str) and len(strdata) == strdata.count(
            "0"
        ) + strdata.count("1"):
            return str.__new__(cls, strdata)
        else:
            raise TypeError(
                "The input should be a binary string composed of '0' and '1'"
            )

    def __and__(self, other):
        selfint = literal_eval("0b" + self)
        otherint = literal_eval("0b" + other)
        resultint = selfint & otherint
        return _BitString(bin(resultint)[2:].zfill(len(self)))

    def __or__(self, other):
        selfint = literal_eval("0b" + self)
        otherint = literal_eval("0b" + other)
        resultint = selfint | otherint
        return _BitString(bin(resultint)[2:].zfill(len(self)))

    def __xor__(self, other):
        selfint = literal_eval("0b" + self)
        otherint = literal_eval("0b" + other)
        resultint = selfint ^ otherint
        return _BitString(bin(resultint)[2:].zfill(len(self)))

    def __rand__(self, other):
        selfint = literal_eval("0b" + self)
        otherint = literal_eval("0b" + other)
        resultint = otherint & selfint
        return _BitString(bin(resultint)[2:].zfill(len(self)))

    def __ror__(self, other):
        selfint = literal_eval("0b" + self)
        otherint = literal_eval("0b" + other)
        resultint = otherint | selfint
        return _BitString(bin(resultint)[2:].zfill(len(self)))

    def __rxor__(self, other):
        selfint = literal_eval("0b" + self)
        otherint = literal_eval("0b" + other)
        resultint = otherint ^ selfint
        return _BitString(bin(resultint)[2:].zfill(len(self)))

    def __repr__(self):
        return "_BitString(" + str.__repr__(self) + ")"

    def index_one(self):
        """Return a list of positions where the element is '1'."""
        return [i for i, n in enumerate(self) if n == "1"]

    def index_zero(self):
        """Return a list of positions where the element is '0'."""
        return [i for i, n in enumerate(self) if n == "0"]

    def contains(self, other):
        """Check if current bitstr1 contains another one bitstr2.

        That is to say, the bitstr2.index_one() is a subset of
        bitstr1.index_one().

        Examples:
            "011011" contains "011000", "011001", "000011"

        Be careful, "011011" also contains "000000". Actually, all _BitString
        objects contain all-zero _BitString of the same length.

        """
        xorbit = self ^ other
        return xorbit.count("1") == self.count("1") - other.count("1")

    def independent(self, other):
        """Check if current bitstr1 is independent of another one bitstr2.

        That is to say the bitstr1.index_one() and bitstr2.index_one() have
        no intersection.

        Be careful, all _BitString objects are independent of all-zero _BitString
        of the same length.
        """
        xorbit = self ^ other
        return xorbit.count("1") == self.count("1") + other.count("1")

    def iscompatible(self, other):
        """Check if current bitstr1 is compatible with another bitstr2.

        Two conditions are considered as compatible:
         1. bitstr1.contain(bitstr2) or vice versa;
         2. bitstr1.independent(bitstr2).

        """
        return self.contains(other) or other.contains(self) or self.independent(other)

    @classmethod
    def from_bool(cls, bools):
        return cls("".join(map(str, map(int, bools))))


def strict_consensus(trees):
    """Search strict consensus tree from multiple trees.

    :Parameters:
        trees : iterable
            iterable of trees to produce consensus tree.

    """
    trees_iter = iter(trees)
    first_tree = next(trees_iter)

    terms = first_tree.get_terminals()
    bitstr_counts, tree_count = _count_clades(itertools.chain([first_tree], trees_iter))

    # Store bitstrs for strict clades
    strict_bitstrs = [
        bitstr for bitstr, t in bitstr_counts.items() if t[0] == tree_count
    ]
    strict_bitstrs.sort(key=lambda bitstr: bitstr.count("1"), reverse=True)
    # Create root
    root = BaseTree.Clade()
    if strict_bitstrs[0].count("1") == len(terms):
        root.clades.extend(terms)
    else:
        raise ValueError("Taxons in provided trees should be consistent")
    # make a bitstr to clades dict and store root clade
    bitstr_clades = {strict_bitstrs[0]: root}
    # create inner clades
    for bitstr in strict_bitstrs[1:]:
        clade_terms = [terms[i] for i in bitstr.index_one()]
        clade = BaseTree.Clade()
        clade.clades.extend(clade_terms)
        for bs, c in bitstr_clades.items():
            # check if it should be the parent of current clade
            if bs.contains(bitstr):
                # remove old bitstring
                del bitstr_clades[bs]
                # update clade childs
                new_childs = [child for child in c.clades if child not in clade_terms]
                c.clades = new_childs
                # set current clade as child of c
                c.clades.append(clade)
                # update bitstring
                bs = bs ^ bitstr
                # update clade
                bitstr_clades[bs] = c
                break
        # put new clade
        bitstr_clades[bitstr] = clade
    return BaseTree.Tree(root=root)


def majority_consensus(trees, cutoff=0):
    """Search majority rule consensus tree from multiple trees.

    This is a extend majority rule method, which means the you can set any
    cutoff between 0 ~ 1 instead of 0.5. The default value of cutoff is 0 to
    create a relaxed binary consensus tree in any condition (as long as one of
    the provided trees is a binary tree). The branch length of each consensus
    clade in the result consensus tree is the average length of all counts for
    that clade.

    :Parameters:
        trees : iterable
            iterable of trees to produce consensus tree.

    """
    tree_iter = iter(trees)
    first_tree = next(tree_iter)

    terms = first_tree.get_terminals()
    bitstr_counts, tree_count = _count_clades(itertools.chain([first_tree], tree_iter))

    # Sort bitstrs by descending #occurrences, then #tips, then tip order
    bitstrs = sorted(
        bitstr_counts.keys(),
        key=lambda bitstr: (bitstr_counts[bitstr][0], bitstr.count("1"), str(bitstr)),
        reverse=True,
    )
    root = BaseTree.Clade()
    if bitstrs[0].count("1") == len(terms):
        root.clades.extend(terms)
    else:
        raise ValueError("Taxons in provided trees should be consistent")
    # Make a bitstr-to-clades dict and store root clade
    bitstr_clades = {bitstrs[0]: root}
    # create inner clades
    for bitstr in bitstrs[1:]:
        # apply majority rule
        count_in_trees, branch_length_sum = bitstr_counts[bitstr]
        confidence = 100.0 * count_in_trees / tree_count
        if confidence < cutoff * 100.0:
            break
        clade_terms = [terms[i] for i in bitstr.index_one()]
        clade = BaseTree.Clade()
        clade.clades.extend(clade_terms)
        clade.confidence = confidence
        clade.branch_length = branch_length_sum / count_in_trees
        bsckeys = sorted(bitstr_clades, key=lambda bs: bs.count("1"), reverse=True)

        # check if current clade is compatible with previous clades and
        # record its possible parent and child clades.
        compatible = True
        parent_bitstr = None
        child_bitstrs = []  # multiple independent childs
        for bs in bsckeys:
            if not bs.iscompatible(bitstr):
                compatible = False
                break
            # assign the closest ancestor as its parent
            # as bsckeys is sorted, it should be the last one
            if bs.contains(bitstr):
                parent_bitstr = bs
            # assign the closest descendant as its child
            # the largest and independent clades
            if (
                bitstr.contains(bs)
                and bs != bitstr
                and all(c.independent(bs) for c in child_bitstrs)
            ):
                child_bitstrs.append(bs)
        if not compatible:
            continue

        if parent_bitstr:
            # insert current clade; remove old bitstring
            parent_clade = bitstr_clades.pop(parent_bitstr)
            # update parent clade childs
            parent_clade.clades = [
                c for c in parent_clade.clades if c not in clade_terms
            ]
            # set current clade as child of parent_clade
            parent_clade.clades.append(clade)
            # update bitstring
            # parent = parent ^ bitstr
            # update clade
            bitstr_clades[parent_bitstr] = parent_clade

        if child_bitstrs:
            remove_list = []
            for c in child_bitstrs:
                remove_list.extend(c.index_one())
                child_clade = bitstr_clades[c]
                parent_clade.clades.remove(child_clade)
                clade.clades.append(child_clade)
            remove_terms = [terms[i] for i in remove_list]
            clade.clades = [c for c in clade.clades if c not in remove_terms]
        # put new clade
        bitstr_clades[bitstr] = clade
        if (len(bitstr_clades) == len(terms) - 1) or (
            len(bitstr_clades) == len(terms) - 2 and len(root.clades) == 3
        ):
            break
    return BaseTree.Tree(root=root)


def adam_consensus(trees):
    """Search Adam Consensus tree from multiple trees.

    :Parameters:
        trees : list
            list of trees to produce consensus tree.

    """
    clades = [tree.root for tree in trees]
    return BaseTree.Tree(root=_part(clades), rooted=True)


def _part(clades):
    """Recursive function for Adam Consensus algorithm (PRIVATE)."""
    new_clade = None
    terms = clades[0].get_terminals()
    term_names = [term.name for term in terms]
    if len(terms) == 1 or len(terms) == 2:
        new_clade = clades[0]
    else:
        bitstrs = {_BitString("1" * len(terms))}
        for clade in clades:
            for child in clade.clades:
                bitstr = _clade_to_bitstr(child, term_names)
                to_remove = set()
                to_add = set()
                for bs in bitstrs:
                    if bs == bitstr:
                        continue
                    elif bs.contains(bitstr):
                        to_add.add(bitstr)
                        to_add.add(bs ^ bitstr)
                        to_remove.add(bs)
                    elif bitstr.contains(bs):
                        to_add.add(bs ^ bitstr)
                    elif not bs.independent(bitstr):
                        to_add.add(bs & bitstr)
                        to_add.add(bs & bitstr ^ bitstr)
                        to_add.add(bs & bitstr ^ bs)
                        to_remove.add(bs)
                # bitstrs = bitstrs | to_add
                bitstrs ^= to_remove
                if to_add:
                    for ta in sorted(to_add, key=lambda bs: bs.count("1")):
                        independent = True
                        for bs in bitstrs:
                            if not ta.independent(bs):
                                independent = False
                                break
                        if independent:
                            bitstrs.add(ta)
        new_clade = BaseTree.Clade()
        for bitstr in sorted(bitstrs):
            indices = bitstr.index_one()
            if len(indices) == 1:
                new_clade.clades.append(terms[indices[0]])
            elif len(indices) == 2:
                bifur_clade = BaseTree.Clade()
                bifur_clade.clades.append(terms[indices[0]])
                bifur_clade.clades.append(terms[indices[1]])
                new_clade.clades.append(bifur_clade)
            elif len(indices) > 2:
                part_names = [term_names[i] for i in indices]
                next_clades = []
                for clade in clades:
                    next_clades.append(_sub_clade(clade, part_names))
                # next_clades = [clade.common_ancestor([clade.find_any(name=name) for name in part_names]) for clade in clades]
                new_clade.clades.append(_part(next_clades))
    return new_clade


def _sub_clade(clade, term_names):
    """Extract a compatible subclade that only contains the given terminal names (PRIVATE)."""
    term_clades = [clade.find_any(name) for name in term_names]
    sub_clade = clade.common_ancestor(term_clades)
    if len(term_names) != sub_clade.count_terminals():
        temp_clade = BaseTree.Clade()
        temp_clade.clades.extend(term_clades)
        for c in sub_clade.find_clades(terminal=False, order="preorder"):
            if c == sub_clade.root:
                continue
            childs = set(c.find_clades(terminal=True)) & set(term_clades)
            if childs:
                for tc in temp_clade.find_clades(terminal=False, order="preorder"):
                    tc_childs = set(tc.clades)
                    tc_new_clades = tc_childs - childs
                    if childs.issubset(tc_childs) and tc_new_clades:
                        tc.clades = list(tc_new_clades)
                        child_clade = BaseTree.Clade()
                        child_clade.clades.extend(list(childs))
                        tc.clades.append(child_clade)
        sub_clade = temp_clade
    return sub_clade


def _count_clades(trees):
    """Count distinct clades (different sets of terminal names) in the trees (PRIVATE).

    Return a tuple first a dict of bitstring (representing clade) and a tuple of its count of
    occurrences and sum of branch length for that clade, second the number of trees processed.

    :Parameters:
        trees : iterable
            An iterable that returns the trees to count

    """
    bitstrs = {}
    tree_count = 0
    for tree in trees:
        tree_count += 1
        clade_bitstrs = _tree_to_bitstrs(tree)
        for clade in tree.find_clades(terminal=False):
            bitstr = clade_bitstrs[clade]
            if bitstr in bitstrs:
                count, sum_bl = bitstrs[bitstr]
                count += 1
                sum_bl += clade.branch_length or 0
                bitstrs[bitstr] = (count, sum_bl)
            else:
                bitstrs[bitstr] = (1, clade.branch_length or 0)
    return bitstrs, tree_count


def get_support(target_tree, trees, len_trees=None):
    """Calculate branch support for a target tree given bootstrap replicate trees.

    :Parameters:
        target_tree : Tree
            tree to calculate branch support for.
        trees : iterable
            iterable of trees used to calculate branch support.
        len_trees : int
            optional count of replicates in trees. len_trees must be provided
            when len(trees) is not a valid operation.

    """
    term_names = sorted(term.name for term in target_tree.find_clades(terminal=True))
    bitstrs = {}

    size = len_trees
    if size is None:
        try:
            size = len(trees)
        except TypeError:
            raise TypeError(
                "Trees does not support len(trees), "
                "you must provide the number of replicates in trees "
                "as the optional parameter len_trees."
            ) from None

    for clade in target_tree.find_clades(terminal=False):
        bitstr = _clade_to_bitstr(clade, term_names)
        bitstrs[bitstr] = (clade, 0)
    for tree in trees:
        for clade in tree.find_clades(terminal=False):
            bitstr = _clade_to_bitstr(clade, term_names)
            if bitstr in bitstrs:
                c, t = bitstrs[bitstr]
                c.confidence = (t + 1) * 100.0 / size
                bitstrs[bitstr] = (c, t + 1)
    return target_tree


def bootstrap(msa, times):
    """Generate bootstrap replicates from a multiple sequence alignment (OBSOLETE).

    :Parameters:
        msa : MultipleSeqAlignment
            multiple sequence alignment to generate replicates.
        times : int
            number of bootstrap times.

    """
    length = len(msa[0])
    i = 0
    while i < times:
        i += 1
        item = None
        for j in range(length):
            col = random.randint(0, length - 1)
            if not item:
                item = msa[:, col : col + 1]
            else:
                item += msa[:, col : col + 1]
        yield item


def bootstrap_trees(alignment, times, tree_constructor):
    """Generate bootstrap replicate trees from a multiple sequence alignment.

    :Parameters:
        alignment : Alignment or MultipleSeqAlignment object
            multiple sequence alignment to generate replicates.
        times : int
            number of bootstrap times.
        tree_constructor : TreeConstructor
            tree constructor to be used to build trees.

    """
    if isinstance(alignment, MultipleSeqAlignment):
        length = len(alignment[0])
        for i in range(times):
            bootstrapped_alignment = None
            for j in range(length):
                col = random.randint(0, length - 1)
                if bootstrapped_alignment is None:
                    bootstrapped_alignment = alignment[:, col : col + 1]
                else:
                    bootstrapped_alignment += alignment[:, col : col + 1]
            tree = tree_constructor.build_tree(alignment)
            yield tree
    else:
        n, m = alignment.shape
        for i in range(times):
            cols = [random.randint(0, m - 1) for j in range(m)]
            tree = tree_constructor.build_tree(alignment[:, cols])
            yield tree


def bootstrap_consensus(alignment, times, tree_constructor, consensus):
    """Consensus tree of a series of bootstrap trees for a multiple sequence alignment.

    :Parameters:
        alignment : Alignment or MultipleSeqAlignment object
            Multiple sequence alignment to generate replicates.
        times : int
            Number of bootstrap times.
        tree_constructor : TreeConstructor
            Tree constructor to be used to build trees.
        consensus : function
            Consensus method in this module: ``strict_consensus``,
            ``majority_consensus``, ``adam_consensus``.

    """
    trees = bootstrap_trees(alignment, times, tree_constructor)
    tree = consensus(trees)
    return tree


def _clade_to_bitstr(clade, tree_term_names):
    """Create a BitString representing a clade, given ordered tree taxon names (PRIVATE)."""
    clade_term_names = {term.name for term in clade.find_clades(terminal=True)}
    return _BitString.from_bool((name in clade_term_names) for name in tree_term_names)


def _tree_to_bitstrs(tree):
    """Create a dict of a tree's clades to corresponding BitStrings (PRIVATE)."""
    clades_bitstrs = {}
    term_names = [term.name for term in tree.find_clades(terminal=True)]
    for clade in tree.find_clades(terminal=False):
        bitstr = _clade_to_bitstr(clade, term_names)
        clades_bitstrs[clade] = bitstr
    return clades_bitstrs


def _bitstring_topology(tree):
    """Generate a branch length dict for a tree, keyed by BitStrings (PRIVATE).

    Create a dict of all clades' BitStrings to the corresponding branch
    lengths (rounded to 5 decimal places).
    """
    bitstrs = {}
    for clade, bitstr in _tree_to_bitstrs(tree).items():
        bitstrs[bitstr] = round(clade.branch_length or 0.0, 5)
    return bitstrs


def _equal_topology(tree1, tree2):
    """Are two trees are equal in terms of topology and branch lengths (PRIVATE).

    (Branch lengths checked to 5 decimal places.)
    """
    term_names1 = {term.name for term in tree1.find_clades(terminal=True)}
    term_names2 = {term.name for term in tree2.find_clades(terminal=True)}
    return (term_names1 == term_names2) and (
        _bitstring_topology(tree1) == _bitstring_topology(tree2)
    )
