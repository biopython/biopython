# Copyright 2018 Aleksandra Jarmolinska. Based on code by Thomas Sicheritz-Ponten,
# Cecilia Alsmark, Markus Piotrowski and Peter Cock. All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Functions for dealing with GC/AT content."""


from __future__ import print_function
from Bio.Data.IUPACData import ambiguous_dna_values


def _calc_gc_content_values():
    gc_u = {}
    gc_a = {}
    unamb = "GCTASWN"
    for b, opts in ambiguous_dna_values.items():
        d = gc_u if b in unamb else gc_a
        d[b] = float(opts.count("C") + opts.count("G")) / len(opts)
        d[b.lower()] = float(opts.count("C") + opts.count("G")) / len(opts)
    return gc_u, gc_a


_GC_CONTENT_VALUE, _GC_CONTENT_VALUE_AMBG = _calc_gc_content_values()


def _calc_base_values(base):
    d = {}
    for b, opts in ambiguous_dna_values.items():
        d[b] = float(opts.count(base)) / len(opts)
        d[b.lower()] = float(opts.count(base)) / len(opts)
    return d


def calc_gc_content(seq, interpret_all=False):
    """Calculate G+C content, returns a fraction (float between 0 and 1).

    Copes with mixed case sequences, and with the ambiguous nucleotides
    S (can be either G or C) and N (can be any base: A or C or G or T )
    when counting the G and C content.

    When interpret_all is True, includes probability of other ambiguous
    to be C or G:
    - for M,K,R,Y it is 0.5
    - for B,V it is 0.66
    - for D,H it is 0.33.

    The percentage is calculated against the full length, e.g.:

    >>> from Bio.SeqUtils.NucUtils import calc_gc_content
    >>> calc_gc_content("ACTGS")
    0.6

    >>> from Bio.SeqUtils.NucUtils import calc_gc_content
    >>> calc_gc_content("ACTGN")
    0.5

    >>> from Bio.SeqUtils.NucUtils import calc_gc_content
    >>> calc_gc_content("ACTGBHNN")
    0.375

    >>> from Bio.SeqUtils.NucUtils import calc_gc_content
    >>> calc_gc_content("ACTGBHNN", True)
    0.5

    Note that this will return zero for an empty sequence.
    """
    len_seq = len(seq)
    if not len_seq:
        return 0.

    gc = sum(_GC_CONTENT_VALUE.get(x,
             _GC_CONTENT_VALUE_AMBG.get(x, 0.) if interpret_all else 0.)
             for x in seq)
    if not gc:
        return None

    return gc / len_seq


def calc_at_content(seq, interpret_all=False):
    """Calculate A+T content, return a fraction (float between 0 and 1).

    Copes with mixed case sequences, and with the ambiguous nucleotides
    W (can be either A or T) and N (can be any base: A or C or G or T )
    when counting the A and T content.

    When interpret_all is True, includes probability of other ambiguous
    to be A or T:
    - for M,K,R,Y it is 0.5
    - for B,V it is 0.33
    - for D,H it is 0.66.

    The percentage is calculated against the full length, e.g.:

    >>> from Bio.SeqUtils.NucUtils import calc_at_content
    >>> calc_at_content("ACTGS")
    0.4

    >>> from Bio.SeqUtils.NucUtils import calc_at_content
    >>> calc_at_content("ACTGN")
    0.5

    >>> from Bio.SeqUtils.NucUtils import calc_at_content
    >>> calc_at_content("ACTGBHNN")
    0.625

    >>> from Bio.SeqUtils.NucUtils import calc_at_content
    >>> calc_at_content("ACTGBHNN", True)
    0.5

    Note that this will return zero for an empty sequence.
    """
    len_seq = len(seq)
    if not len_seq:
        return 0.

    at = sum(1 - _GC_CONTENT_VALUE.get(x,
             _GC_CONTENT_VALUE_AMBG.get(x, 0.) if interpret_all else 0.)
             for x in seq)
    if not at:
        return None

    return at / len_seq


def calc_gc_skew(seq, interpret_all=False):
    """Calculate GC skew (G-C)/(G+C) for the whole sequence.

    Returns the ratio (float).

    When interpret_all is True includes also ambiguous nucleotides,
    as their probability of being G or C:
    - G is 0.5 in S, K and R, 0.33 in B, D, and V and 0.25 in N
    - C is 0.5 in S, M and Y, 0.33 in B, H, and V and 0.25 in N.

    Returns zero for an empty sequence.
    Returns None if there is no G+C content.
    """
    g_dict = _calc_base_values("G")

    c_dict = _calc_base_values("C")

    len_seq = len(seq)
    if not len_seq:
        return 0

    g = sum(g_dict.get(x, 0.) if interpret_all else x in "Gg" for x in seq)
    c = sum(c_dict.get(x, 0.) if interpret_all else x in "Cc" for x in seq)

    if not g + c:
        return None

    skew = (g - c) / float(g + c)

    return skew


def calc_at_skew(seq, interpret_all=False):
    """Calculate AT skew (A-T)/(A+T) for the whole sequence.

    Returns the ratio (float).

    When interpret_all is True includes also ambiguous nucleotides,
    as their probability of being G or C:
    - A is 0.5 in W, M and R, 0.33 in D, H, and V and 0.25 in N
    - C is 0.5 in W, K and Y, 0.33 in B, D, and H and 0.25 in N.

    Returns zero for an empty sequence.
    Returns None if there is no G+C content.
    """
    a_dict = _calc_base_values("A")

    t_dict = _calc_base_values("T")

    len_seq = len(seq)
    if not len_seq:
        return 0

    a = sum(a_dict.get(x, 0.) if interpret_all else x in "Aa" for x in seq)
    t = sum(t_dict.get(x, 0.) if interpret_all else x in "Tt" for x in seq)

    if not a + t:
        return None

    skew = (a - t) / float(a + t)

    return skew
