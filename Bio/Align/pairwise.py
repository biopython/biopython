# Copyright 2001 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""This module implements some dynamic programming code to align
sequences.  The algorithm is general and allows things like non-linear
gap penalties and position-dependent gap penalties.  For a faster
but more restrictive algorithm, see the gotoh module.

Functions:
align_global   Do a global alignment between two sequences.
align_local    Do a local alignment between two sequences.

identity_match    Calculate match scores based on identity.
dictionary_match  Calculate match scores from a dictionary lookup.

affine_penalty    Use linear gap penalties.
no_penalty        Use no gap penalties.

"""
from support import *

# Because of differences in how this and the gotoh algorithms
# generate score matrices, they may stop off at different points in
# the tracebacks during a local alignment.  For example, if you align:
# AAA to ATAA with gap penalty (-1, -1) and match function (1, -1),
# gotoh will get:        general will get:
# -AAA                   A-AA
# --||                   ||||
# ATAA                   ATAA
#
# The score matrix in the general algorithm never gets a 0, because
# the gap is implicit in the traceback.

def align_global(sequenceA, sequenceB, gap_A_fn, gap_B_fn,
                 match_fn=None, **keywds):
    """align_global(sequenceA, sequenceB, gap_A_fn, gap_B_fn[, match_fn]) ->
    alignments

    Do a global alignment on sequenceA and sequenceB.  gap_A_fn and
    gap_B_fn are callback functions that take the gap length, the
    sequence, and the index of the start of the gap and returns a gap
    penalty, which should be negative.  match_fn is a callback
    function that takes two characters and returns the score.  By
    default, will give a score of 1 for matches and 0 for mismatches.

    alignments is a list of tuples (seqA, seqB, score, begin, end).
    seqA and seqB are strings showing the alignment between the
    sequences.  score is the score of the alignment.  begin and end
    are indexes into seqA and seqB that indicate the where the
    alignment occurs.  Here, it should be the whole sequence.

    """
    if match_fn is None:
        match_fn = identity_match(1, 0)
    params = [
        ('global_alignment', 1),
        ('penalize_end_gaps', 1),
        ('gap_char', '-'),
        ]
    for name, default in params:
        keywds[name] = keywds.get(name, default)
    alignments = _align(sequenceA, sequenceB, match_fn, gap_A_fn, gap_B_fn,
                        **keywds)
    return alignments

def align_local(sequenceA, sequenceB, gap_A_fn, gap_B_fn,
                 match_fn=None, **keywds):
    """align_local(sequenceA, sequenceB, gap_A_fn, gap_B_fn[, match_fn]) ->
    alignments

    Do a local alignment on sequenceA and sequenceB.  gap_A_fn and
    gap_B_fn are callback functions that take the gap length, the
    sequence, and the index of the start of the gap and returns a gap
    penalty, which should be negative.  match_fn is a callback
    function that takes two characters and returns the score.  By
    default, will give a score of 1 for matches and 0 for mismatches.

    alignments is a list of tuples (seqA, seqB, score, begin, end).
    seqA and seqB are strings showing the alignment between the
    sequences.  score is the score of the alignment.  begin and end
    are indexes into seqA and seqB that indicate the local region
    where the alignment occurs.

    """
    if match_fn is None:
        match_fn = identity_match(1, 0)
    params = [
        ('global_alignment', 0),
        ('penalize_end_gaps', 0),
        ('gap_char', '-'),
        ]
    for name, default in params:
        keywds[name] = keywds.get(name, default)
    return _align(sequenceA, sequenceB, match_fn, gap_A_fn, gap_B_fn, **keywds)

def _align(sequenceA, sequenceB, match_fn, gap_A_fn, gap_B_fn,
           global_alignment, penalize_end_gaps, gap_char):
    score_matrix, traceback_matrix = _make_score_matrix(
        sequenceA, sequenceB, match_fn, gap_A_fn, gap_B_fn,
        global_alignment, penalize_end_gaps)

    score, indexes = _find_best_score(
        score_matrix, sequenceA, sequenceB, gap_A_fn, gap_B_fn,
        global_alignment, penalize_end_gaps)

##     for i in range(len(score_matrix)):
##         for j in range(len(score_matrix[i])):
##             print score_matrix[i][j],
##         print
##     print score, indexes

    alignments = _recover_alignments(sequenceA, sequenceB, score, indexes,
                                     match_fn,
                                     score_matrix, traceback_matrix,
                                     global_alignment, penalize_end_gaps,
                                     gap_char)
    return alignments

def _make_score_matrix(sequenceA, sequenceB, match_fn, gap_A_fn, gap_B_fn,
                       global_alignment, penalize_end_gaps):
    lenA, lenB = len(sequenceA), len(sequenceB)
    # Create the score and traceback matrices.  These should be in the
    # shape:
    # sequenceA (down) x sequenceB (across)
    score_matrix, traceback_matrix = [], []
    for i in range(lenA):
        score_matrix.append([None] * lenB)
        traceback_matrix.append([[None]] * lenB)

    if not sequenceA or not sequenceB:
        return score_matrix, traceback_matrix

    # The borders are special cases.  Handle them first.
    for i in range(lenA):
        # Align the first residue in sequenceB to the ith residue in
        # sequence A.  This is like opening up i gaps at the beginning
        # of sequence B.
        score = match_fn(sequenceA[i], sequenceB[0])
        if penalize_end_gaps:
            score += gap_B_fn(i, sequenceB, 0)
        score_matrix[i][0] = score
    for i in range(1, lenB):
        score = match_fn(sequenceA[0], sequenceB[i])
        if penalize_end_gaps:
            score += gap_A_fn(i, sequenceA, 0)
        score_matrix[0][i] = score

    # fill in the score_matrix
    for row in range(1, lenA):
        for col in range(1, lenB):
            # The score by extending the alignment, no gaps.
            best_score = score_matrix[row-1][col-1]
            best_score_rint = rint(best_score)
            best_indexes = [(row-1, col-1)]

            # Look for alignment in previous row.  Open gap in sequenceA.
            for i in range(0, col-1):
                score = score_matrix[row-1][i] + gap_A_fn(
                    col-1-i, sequenceA, i)
                score_rint = rint(score)
                if score_rint == best_score_rint:
                    best_indexes.append((row-1, i))
                elif score_rint > best_score_rint:
                    best_score, best_score_rint = score, score_rint
                    best_indexes = [(row-1, i)]
            
            # Look for alignment in previous col.  Open gap in sequenceB.
            for i in range(0, row-1):
                score = score_matrix[i][col-1] + gap_B_fn(
                    row-1-i, sequenceB, i)
                score_rint = rint(score)
                if score_rint == best_score_rint:
                    best_indexes.append((i, col-1))
                elif score_rint > best_score_rint:
                    best_score, best_score_rint = score, score_rint
                    best_indexes = [(i, col-1)]

            score_matrix[row][col] = best_score + \
                                     match_fn(sequenceA[row], sequenceB[col])
            if not global_alignment and score_matrix[row][col] < 0:
                score_matrix[row][col] = 0
            traceback_matrix[row][col] = best_indexes
    return score_matrix, traceback_matrix

def _recover_alignments(sequenceA, sequenceB, score, indexes,
                        match_fn,
                        score_matrix, traceback_matrix, global_alignment,
                        penalize_end_gaps, gap_char):
    # Recover the alignments.  This is a recursive procedure, but it's
    # implemented here iteratively with a stack.
    lenA, lenB = len(sequenceA), len(sequenceB)
    tracebacks = [] # list of (seq1, seq2, score, begin, end)
    in_process = [] # list of ([same as tracebacks], prev_pos, next_pos)
    
    # Initialize the in_process stack
    for row, col in indexes:
        if global_alignment:
            begin, end = None, None
        else:
            begin, end = None, -max(lenA-row, lenB-col)+1
            if not end:
                end = None
        in_process.append(
            ('', '', score, begin, end, (lenA, lenB), (row, col)))
    while in_process:
        seqA, seqB, score, begin, end, prev_pos, next_pos = in_process.pop()
        prevA, prevB = prev_pos
        if next_pos is None:
            prevlen = len(seqA)
            # add the rest of the sequences
            seqA = sequenceA[:prevA] + seqA
            seqB = sequenceB[:prevB] + seqB
            lseqA, lseqB = len(seqA), len(seqB)
            # add the rest of the gaps
            if lseqA < lseqB:
                seqA = lpad(seqA, gap_char, lseqB-lseqA)
            elif len(seqB) < len(seqA):
                seqB = lpad(seqB, gap_char, lseqA-lseqB)
            # Now make sure begin is set.
            if begin is None:
                if global_alignment:
                    begin = 0
                else:
                    begin = len(seqA) - prevlen
            tracebacks.append((seqA, seqB, score, begin, end))
        else:
            nextA, nextB = next_pos
            nseqA, nseqB = prevA-nextA, prevB-nextB
            maxseq = max(nseqA, nseqB)
            ngapA, ngapB = maxseq-nseqA, maxseq-nseqB
            seqA = sequenceA[nextA:nextA+nseqA] + gap_char*ngapA + seqA
            seqB = sequenceB[nextB:nextB+nseqB] + gap_char*ngapB + seqB
            prev_pos = next_pos
            # local alignment stops early if score falls < 0
            if not global_alignment and score_matrix[nextA][nextB] <= 0:
                begin = max(prevA, prevB)
                in_process.append(
                    (seqA, seqB, score, begin, end, prev_pos, None))
            else:
                for next_pos in traceback_matrix[nextA][nextB]:
                    in_process.append(
                        (seqA, seqB, score, begin, end, prev_pos, next_pos))
    return clean_alignments(tracebacks)

def _find_best_score(score_matrix, sequenceA, sequenceB, gap_A_fn, gap_B_fn,
                     global_alignment, penalize_end_gaps):
    if global_alignment:
        score, indexes = find_global_best(
            sequenceA, sequenceB,
            score_matrix, penalize_end_gaps, gap_A_fn, gap_B_fn)
    else:
        score, indexes = find_local_best(score_matrix)
    return score, indexes
