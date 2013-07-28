# Copyright (C) 2013 by Yanbo Ye (yeyanbo289@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Classes and methods for finding consensus trees"""
from ast import literal_eval
from Bio.Phylo import BaseTree

class BitString(str):
    """Assistant class of binary string data used for storing and
     counting compatible clades in consensus tree searching. It includes
     some binary manipulation(&|^~) methods."""

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
        """
        xorbit = self ^ other
        if xorbit.count('1') == self.count('1') - other.count('1'):
            return True
        else:
            return False


def strict_consensus(trees):
    """Search stric consensus tree from a multiple trees"""
    terms = trees[0].get_terminals()
    bitstr_counts = _count_clades(trees)
    # store bitstrs for strict clades
    strict_bitstrs = []
    for bitstr, count in bitstr_counts.items():
        if count == len(trees):
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


def majority_consensus(trees):
    pass


def adam_consensus(trees):
    pass


def _count_clades(trees):
    """Count all clades(considered as the same if terminals are the same)
    in the trees and return a dict of bitstring(representing clade) and its
    time of presence"""
    term_names = [term.name for term in trees[0].get_terminals()]
    bitstrs = {}
    for tree in trees:
        for clade in tree.get_nonterminals():
            clade_term_names = [term.name for term in clade.get_terminals()]
            boolvals = [name in clade_term_names for name in term_names]
            bitstr = BitString(''.join(map(str, map(int, boolvals))))
            if bitstr in bitstrs.keys():
                bitstrs[bitstr] = bitstrs[bitstr] + 1
            else:
                bitstrs[bitstr] = 1
    return bitstrs
