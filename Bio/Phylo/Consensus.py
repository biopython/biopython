# Copyright (C) 2013 by Yanbo Ye (yeyanbo289@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

""" Classes and methods for finding consensus trees.
    
    This module contains a ``BitString`` class to assist the consensus
    tree searching and some common consensus algorithms such as strict, 
    majority rule and adam consensus.
"""
from ast import literal_eval
from Bio.Phylo import BaseTree

class BitString(str):
    """Assistant class of binary string data used for storing and
     counting compatible clades in consensus tree searching. It includes
     some binary manipulation(&|^~) methods.

    BitString is a sub-class of ``str`` object that only accepts two
    characters('0' and '1'), with additional functions for binary-like
    manipulation(&|^~). It is used to count and store the clades in
    multiple trees in consensus tree searching. During counting, the
    clades will be considered the same if their terminals(in terms of
    ``name`` attribute) are the same.
    
    For example, let's say two trees are provided as below to search
    their strict consensus tree:
        
        tree1: (((A, B), C),(D, E))
        tree2: ((A, (B, C)),(D, E))

    For both trees, a BitString object '11111' will represent their
    root clade. Each '1' stands for the terminal clade in the list
    [A, B, C, D, E](the order might not be the same, it's determined
    by the ``get_terminal`` method of the first tree provided). For
    the clade ((A, B), C) in tree1 and (A, (B, C)) in tree2, they both
    can be represented by '11100'. Similarly, '11000' represents clade
    (A, B) in tree1, '01100' represents clade (B, C) in tree2, and '00011'
    represents clade (D, E) in both trees. 

    So, with the ``_count_clades`` function in this module, finnally we
    can get the clade counts and their BitString representation as 
    follows(the root and terminals are omitted):

        clade   BitString   count
        ABC     '11100'     2
        DE      '00011'     2
        AB      '11000'     1
        BC      '01100'     1

    To get the BitString representation of a clade, we can use the following
    code snippet:

        # suppose we are provided with a tree list, the first thing to do is 
        # to get all the terminal names in the first tree
        term_names = [term.name for term in trees[0].get_terminals()]
        # for a specific clade in any of the tree, also get its terminal names
        clade_term_names = [term.name for term in clade.get_terminals()]
        # then create a boolean list 
        boolvals = [name in clade_term_names for name in term_names]
        # create the string version and pass it to BitString  
        bitstr = BitString(''.join(map(str, map(int, boolvals))))

    To convert back:
        
        # get all the terminal clades of the first tree
        terms = [term for term in trees[0].get_terminals()]
        # get the index of terminal clades in bitstr
        index_list = bitstr.index_one()
        # get all terminal clades by index
        clade_terms = [terms[i] for i in index_list]
        # create a new calde and append all the terminal clades
        new_clade = BaseTree.Clade()
        new_clade.clades.extend(clade_terms)


    Example
    -------

    >>> from Bio.Phylo.Consensus import BitString
    >>> bitstr1 = BitString('11111')
    >>> bitstr2 = BitString('11100')
    >>> bitstr3 = BitString('01101')
    >>> bitstr1
    BitString('11111')
    >>> bitstr2 & bitstr3
    BitString('01100')
    >>> bitstr2 | bitstr3
    BitString('11101')
    >>> bitstr2 ^ bitstr3
    BitString('10001')
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
    >>> bitstr2.independent(bitstr4)
    True
    >>> bitstr1.iscompatible(bitstr2)
    True
    >>> bitstr2.iscompatible(bitstr3)
    False
    >>> bitstr2.iscompatible(bitstr4)
    True
     """

    def __new__(cls, strdata):
        """init from a binary string data"""
        if isinstance(strdata, str) and len(strdata) == strdata.count('0') + strdata.count('1'):
            return str.__new__(cls, strdata)
        else:
            raise TypeError("The input should be a binary string composed by '0' and '1'")

    def __and__(self, other):
        selfint = literal_eval('0b' + self)
        otherint = literal_eval('0b' + other)
        resultint = selfint & otherint
        return BitString(bin(resultint)[2:].zfill(len(self)))

    def __or__(self, other):
        selfint = literal_eval('0b' + self)
        otherint = literal_eval('0b' + other)
        resultint = selfint | otherint
        return BitString(bin(resultint)[2:].zfill(len(self)))

    def __xor__(self, other):
        selfint = literal_eval('0b' + self)
        otherint = literal_eval('0b' + other)
        resultint = selfint ^ otherint
        return BitString(bin(resultint)[2:].zfill(len(self)))


    def __rand__(self, other):
        selfint = literal_eval('0b' + self)
        otherint = literal_eval('0b' + other)
        resultint =  otherint & selfint
        return BitString(bin(resultint)[2:].zfill(len(self)))

    def __ror__(self, other):
        selfint = literal_eval('0b' + self)
        otherint = literal_eval('0b' + other)
        resultint = otherint | selfint
        return BitString(bin(resultint)[2:].zfill(len(self)))

    def __rxor__(self, other):
        selfint = literal_eval('0b' + self)
        otherint = literal_eval('0b' + other)
        resultint = otherint ^ selfint
        return BitString(bin(resultint)[2:].zfill(len(self)))

    def __repr__(self):
        return 'BitString(' + str.__repr__(self) + ')'

    def index_one(self):
        """Return a list of positions where the element is '1'"""
        return [i for i, n in enumerate(self) if n == '1']

    def index_zero(self):
        """Return a list of positions where the element is '0'"""
        return [i for i, n in enumerate(self) if n == '0']

    def contains(self, other):
        """Check if current bitstr1 contains another one bitstr2. That is
         to say the bitstr2.index_one() is a subset of bitstr1.index_one().
        Examples:
            "011011" contains "011000", "011001", "000011"
        
        Be careful, "011011" also contains "000000". Actually, all BitString
        objects contain all-zero BitString of the same length. 
        """
        xorbit = self ^ other
        if xorbit.count('1') == self.count('1') - other.count('1'):
            return True
        else:
            return False

    def independent(self, other):
        """Check if current bitstr1 is independent of another one bitstr2.
        That is to say the bitstr1.index_one() and bitstr2.index_one() have
        no intersection.

        Be careful, all BitString objects are independent of all-zero BitString
        of the same length. 
        """
        xorbit = self ^ other
        if xorbit.count('1') == self.count('1') + other.count('1'):
            return True
        else:
            return False

    def iscompatible(self, other):
        """Check if current bitstr1 is compatible with another bitstr2. Two
        conditions are considered as compatible: 
        1. bitstr1.contain(bitstr2) or vise versa;
        2. bitstr1.independent(bitstr2).
        """
        return self.contains(other) or other.contains(self) or self.independent(other)


def strict_consensus(trees):
    """Search stric consensus tree from multiple trees"""
    terms = trees[0].get_terminals()
    bitstr_counts = _count_clades(trees)
    # store bitstrs for strict clades
    strict_bitstrs = []
    for bitstr, t in bitstr_counts.items():
        if t[0] == len(trees):
            strict_bitstrs.append(bitstr)
    # sort bitstrs and create root
    strict_bitstrs.sort(key=lambda bitstr: bitstr.count('1'), reverse=True)
    root = BaseTree.Clade()
    if strict_bitstrs[0].count('1') == len(terms):
        root.clades.extend(terms)
    else:
        raise ValueError('Taxons in provided trees should be consistent')
    # make a bitstr to clades dict and store root clade
    bitstr_clades = {strict_bitstrs[0]: root}
    # create inner clades
    for bitstr in strict_bitstrs[1:]:
        index_list = bitstr.index_one()
        clade_terms = [terms[i] for i in index_list]
        clade = BaseTree.Clade()
        clade.clades.extend(clade_terms)
        for bs, c in bitstr_clades.items():
            # check if it should be the parent of current clade
            if bs.contains(bitstr):
                # remove old bitstring
                del bitstr_clades[bs]
                # update clade childs
                childs = c.clades
                new_childs = [child for child in childs if child not in clade_terms]
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
    create a relaxed binary consensus tree in any condition(as long as one 
    of the provided trees is a binary tree). The branch length of each consensus
    clade in the result consensus tree is the average length of all counts
    for that clade.
    """
    terms = trees[0].get_terminals()
    bitstr_counts = _count_clades(trees)
    bitstrs = bitstr_counts.keys()
    # sort bitstrs from higher count to lower
    bitstrs.sort(key=lambda bitstr: (bitstr_counts[bitstr][0], bitstr.count('1')), reverse=True)
    root = BaseTree.Clade()
    if bitstrs[0].count('1') == len(terms):
        root.clades.extend(terms)
    else:
        raise ValueError('Taxons in provided trees should be consistent')
    # make a bitstr to clades dict and store root clade
    bitstr_clades = {bitstrs[0]: root}
    # create inner clades
    for bitstr in bitstrs[1:]:
        # apply majority rule
        confidence = bitstr_counts[bitstr][0] * 100.0 / len(trees)
        if confidence < cutoff * 100.0:
            break
        index_list = bitstr.index_one()
        clade_terms = [terms[i] for i in index_list]
        clade = BaseTree.Clade()
        clade.clades.extend(clade_terms)
        clade.confidence = confidence
        clade.branch_length = bitstr_counts[bitstr][1] * 1.0 / bitstr_counts[bitstr][0]
        bsckeys = bitstr_clades.keys()
        bsckeys.sort(key=lambda bs: bs.count('1'), reverse=True)

        # check if current clade is compatible with previous clades and
        # record it's possible parent and child clades.
        compatible = True
        parent_bitstr = None
        child_bitstrs = [] # multiple independent childs
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
            if bitstr.contains(bs) and bs != bitstr and all([c.independent(bs) for c in child_bitstrs]):
                child_bitstrs.append(bs)
        if not compatible:
            continue

        if parent_bitstr:
            # insert current clade
            parent_clade = bitstr_clades[parent_bitstr]
            # remove old bitstring
            del bitstr_clades[parent_bitstr]
            # update parent clade childs
            parent_childs = parent_clade.clades
            new_parent_childs = [c for c in parent_childs if c not in clade_terms]
            parent_clade.clades = new_parent_childs
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
            new_clade_childs = [c for c in clade.clades if c not in remove_terms]
            clade.clades = new_clade_childs
        # put new clade
        bitstr_clades[bitstr] = clade
        if len(bitstr_clades) == len(terms) - 1:
            break
        elif len(bitstr_clades) == len(terms) - 2 and len(root.clades) == 3:
            break
    return BaseTree.Tree(root=root)


def adam_consensus(trees):
    """Search Adam Consensus tree from multiple trees"""
    clades = [tree.root for tree in trees]
    return BaseTree.Tree(root=_part(clades), rooted=True)


def _part(clades):
    new_clade = None
    terms = [term for term in clades[0].get_terminals()]
    term_names = [term.name for term in clades[0].get_terminals()]
    if len(terms) == 1 or len(terms) == 2:
        new_clade = clades[0]
    else:
        bitstrs = set([BitString('1' * len(terms))])
        for clade in clades:
            for childs in clade.clades:
                clade_term_names = [term.name for term in childs.get_terminals()]
                boolvals = [name in clade_term_names for name in term_names]
                bitstr = BitString(''.join(map(str, map(int, boolvals))))
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
                bitstrs = bitstrs | to_add
                bitstrs = bitstrs ^ to_remove
        new_clade = BaseTree.Clade()
        for bitstr in bitstrs:
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
                next_clades = [clade.common_ancestor([clade.find_any(name=name) for name in part_names]) for clade in clades]
                new_clade.clades.append(_part(next_clades))
    return new_clade


def _count_clades(trees):
    """Count all clades(considered as the same if terminals are the same)
    in the trees and return a dict of bitstring(representing clade) and 
    a tuple of its time of presence and sum of branch length for that clade."""
    term_names = [term.name for term in trees[0].get_terminals()]
    bitstrs = {}
    for tree in trees:
        for clade in tree.get_nonterminals():
            clade_term_names = [term.name for term in clade.get_terminals()]
            boolvals = [name in clade_term_names for name in term_names]
            bitstr = BitString(''.join(map(str, map(int, boolvals))))
            if bitstr in bitstrs.keys():
                time, sum_bl = bitstrs[bitstr]
                time += 1
                sum_bl += clade.branch_length
                bitstrs[bitstr] = (time, sum_bl)
            else:
                bitstrs[bitstr] = (1, clade.branch_length)
    return bitstrs


def get_support(target_tree, trees):
    """Calculate branch support given a target tree and a list of bootstrap
    replicate trees"""
    term_names = [term.name for term in target_tree.get_terminals()]
    bitstrs = {}
    size = len(trees)
    for clade in target_tree.get_nonterminals():
        clade_term_names = [term.name for term in clade.get_terminals()]
        boolvals = [name in clade_term_names for name in term_names]
        bitstr = BitString(''.join(map(str, map(int, boolvals))))
        bitstrs[bitstr] = (clade, 0)
    for tree in trees:
        for clade in tree.get_nonterminals():
            clade_term_names = [term.name for term in clade.get_terminals()]
            boolvals = [name in clade_term_names for name in term_names]
            bitstr = BitString(''.join(map(str, map(int, boolvals))))
            if bitstr in bitstrs.keys():
                c,t = bitstrs[bitstr]
                c.confidence = (t + 1) * 100.0 / size
                bitstrs[bitstr] = (c, t + 1)
    return target_tree