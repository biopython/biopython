# Copyright 2001 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""This internal module provides code to support dynamic programming
code.  Do not use this module unless you know what you are doing!

Functions:
rint     Round a floating point up a few digits.

calc_affine_penalty  Calculate an affine gap penalty.
find_global_best     Find the best score from a global alignment.
find_local_best      Find the best score from a local alignment.

pad_until_equal      Pad two strings with chars until they're equal length.
lpad_until_equal     Left pad two strings until they're equal length.
pad                  Pad a string with a characters.
lpad                 Left pad a string with characters.

Classes:
identity_match    Calculate match scores based on identity.
dictionary_match  Calculate match scores from a dictionary lookup.
affine_penalty    Calculate a linear gap penalty.
no_penalty        No gap penalty.

"""
from Bio.Tools import listfns

_PRECISION = 1000
def rint(x, precision=_PRECISION):
    return int(x * precision + 0.5)

class identity_match:
    """identity_match([match][, mismatch]) -> match_fn

    Create a match function for use in an alignment.  match and
    mismatch are the scores to give when two residues are equal or
    unequal.  By default, match is 1 and mismatch is 0.

    """
    def __init__(self, match=1, mismatch=0):
        self.match = match
        self.mismatch = mismatch
    def __call__(self, charA, charB):
        if charA == charB:
            return self.match
        return self.mismatch

class dictionary_match:
    """dictionary_match(score_dict[, symmetric]) -> match_fn

    Create a match function for use in an alignment.  score_dict is a
    dictionary where the keys are tuples (residue 1, residue 2) and
    the values are the match scores between those residues.  symmetric
    is a flag that indicates whether the scores are symmetric.  If
    true, then if (res 1, res 2) doesn't exist, I will use the score
    at (res 2, res 1).

    """
    def __init__(self, score_dict, symmetric=0):
        self.score_dict = score_dict
        self.symmetric = symmetric
    def __call__(self, charA, charB):
        if symmetric and not score_dict.has_key((charA, charB)):
            # If the score dictionary is symmetric, then look up the
            # score both ways.
            charB, charA = charA, charB
        return score_dict[(charA, charB)]

class affine_penalty:
    """affine_penalty(open, extend[, count_first]) -> gap_fn

    Create a gap function for use in an alignment.

    """
    def __init__(self, open, extend, count_first=0):
        if open > 0 or extend > 0:
            raise ValueError, "Gap penalties should be negative."
        self.open, self.extend = open, extend
        self.count_first = count_first
    def __call__(self, length, sequence, index):
        return calc_affine_penalty(
            length, self.open, self.extend, self.count_first)

class no_penalty:
    """no_penalty() -> gap_fn

    Create a gap function for use in an alignment.

    """
    def __call__(self):
        return 0

def pad_until_equal(s1, s2, char):
    ls1, ls2 = len(s1), len(s2)
    if ls1 < ls2:
        s1 = pad(s1, char, ls2-ls1)
    elif ls2 < ls1:
        s2 = pad(s2, char, ls1-ls2)
    return s1, s2

def lpad_until_equal(s1, s2, char):
    ls1, ls2 = len(s1), len(s2)
    if ls1 < ls2:
        s1 = lpad(s1, char, ls2-ls1)
    elif ls2 < ls1:
        s2 = lpad(s2, char, ls1-ls2)
    return s1, s2

def pad(s, char, n):
    return s + char*n

def lpad(s, char, n):
    return char*n + s

def calc_affine_penalty(length, open, extend, count_first):
    if length <= 0:
        return 0
    penalty = open + extend * length
    if not count_first:
        penalty -= extend
    return penalty

def find_global_best(sequenceA, sequenceB, 
                     score_matrix, penalize_end_gaps, gap_A_fn, gap_B_fn):
    # gap_A_fn and gap_B_fn only necessary if penalize_end_gaps
    nrows, ncols = len(score_matrix), len(score_matrix[0])
    best_score = best_score_rint = None
    best_indexes = []
    
    # Search all rows in the last column.
    for row in range(nrows):
        # Find the score, penalizing end gaps if necessary.
        score = score_matrix[row][ncols-1]
        if penalize_end_gaps:
            score += gap_B_fn(nrows-row-1, sequenceB, ncols)
        # Check to see whether this score exceeds the previous best.
        score_rint = rint(score)
        if best_score_rint is None or score_rint > best_score_rint:
            best_score, best_score_rint = score, score_rint
            best_indexes = [(row, ncols-1)]
        elif best_score_rint == score_rint:
            best_indexes.append((row, ncols-1))

    # Search all columns in the last row.
    for col in range(ncols-1):
        score = score_matrix[nrows-1][col]
        if penalize_end_gaps:
            score += gap_A_fn(ncols-col-1, sequenceA, nrows)
        score_rint = rint(score)
        if best_score_rint is None or score_rint > best_score_rint:
            best_score, best_score_rint = score, score_rint
            best_indexes = [(nrows-1, col)]
        elif best_score_rint == score_rint:
            best_indexes.append((nrows-1, col))

    return best_score, best_indexes

def find_local_best(score_matrix):
    # local alignment looks everywhere
    nrows, ncols = len(score_matrix), len(score_matrix[0])
    best_score = best_score_rint = None
    best_indexes = []
    for row in range(nrows):
        for col in range(ncols):
            score = score_matrix[row][col]
            score_rint = rint(score)
            if best_score_rint is None or score_rint > best_score_rint:
                best_score, best_score_rint = score, score_rint
                best_indexes = [(row, col)]
            elif score_rint == best_score_rint:
                best_indexes.append((row, col))
    return best_score, best_indexes

def clean_alignments(alignments):
    alignments.sort()
    i = 0
    while i < len(alignments):
        # Get rid of duplicates.
        if i < len(alignments)-1 and alignments[i] == alignments[i+1]:
            del alignments[i]
            continue
        seqA, seqB, score, begin, end = alignments[i]
        # Make sure end is set reasonably.
        if end == None:
            end = len(seqA)
        elif end < 0:
            end = end + len(seqA)
        # If there's no alignment here, get rid of it.
        if begin >= end:
            del alignments[i]
            continue
        alignments[i] = seqA, seqB, score, begin, end
        i += 1
    return alignments
