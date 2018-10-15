# Copyright 2018 Aleksandra Jarmolinska. Based on code by Thomas Sicheritz-Ponten,
# Cecilia Alsmark, Markus Piotrowski and Peter Cock. All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Functions for dealing with GC/AT content."""


from __future__ import print_function

_GC_CONTENT_VALUE = {
    "G": 1.,
    "C": 1.,
    "A": 0.,
    "T": 0.,
    "S": 1.,
    "W": 0.,
    "N": 0.5
}
_GC_CONTENT_VALUE_AMBG = {
    "M": 0.5,
    "K": 0.5,
    "R": 0.5,
    "Y": 0.5,

    "B": 0.66,
    "D": 0.33,
    "H": 0.33,
    "V": 0.66,
}


def calc_gc_content(seq, interpret_all=False):
    """Calculate G+C content, returns a fraction (float between 0 and 1).

    Copes with mixed case sequences, and with the ambiguous nucleotides S 
    (can be either G or C) and N (can be any base: A or C or G or T ) 
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
    >>> calc_gc_content("ACTGB")
    0.4

    >>> from Bio.SeqUtils.NucUtils import calc_gc_content
    >>> calc_gc_content("ACTGB", True)
    0.532

    Note that this will return zero for an empty sequence.
    """
    len_seq = len(seq)
    if not len_seq:
        return 0.

    gc = sum(_GC_CONTENT_VALUE.get(x.upper(),
             _GC_CONTENT_VALUE_AMBG.get(x.upper(), 0.) if interpret_all else 0.)
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
    >>> calc_at_content("ACTGBNNN")
    0.5625

    >>> from Bio.SeqUtils.NucUtils import calc_at_content
    >>> calc_at_content("ACTGBNNN", True)
    0.48

    Note that this will return zero for an empty sequence.
    """
    len_seq = len(seq)
    if not len_seq:
        return 0.

    at = sum(1 - _GC_CONTENT_VALUE.get(x.upper(),
             _GC_CONTENT_VALUE_AMBG.get(x.upper(), 0.) if interpret_all else 0.)
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
    g_dict = {"G": 1., "S": 0.5, "K": 0.5, "R": 0.5,
              "B": 0.33, "D": 0.33,
              "V": 0.33, "N": 0.25}

    c_dict = {"C": 1., "S": 0.5, "M": 0.5, "Y": 0.5,
              "B": 0.33, "H": 0.33,
              "V": 0.33, "N": 0.25}

    len_seq = len(seq)
    if not len_seq:
        return 0

    g = sum(g_dict.get(x.upper(), 0.) if interpret_all else x.upper() == "G" for x in seq)
    c = sum(c_dict.get(x.upper(), 0.) if interpret_all else x.upper() == "C" for x in seq)

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
    a_dict = {"A": 1., "W": 0.5, "M": 0.5, "R": 0.5,
              "H": 0.33, "D": 0.33,
              "V": 0.33, "N": 0.25}
    t_dict = {"T": 1., "S": 0.5, "K": 0.5, "Y": 0.5,
              "B": 0.33, "H": 0.33,
              "D": 0.33, "N": 0.25}

    len_seq = len(seq)
    if not len_seq:
        return 0

    a = sum(a_dict.get(x.upper(), 0.) if interpret_all else x.upper() == "A" for x in seq)
    t = sum(t_dict.get(x.upper(), 0.) if interpret_all else x.upper() == "T" for x in seq)

    if not a + t:
        return None

    skew = (a - t) / float(a + t)

    return skew
