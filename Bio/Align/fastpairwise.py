# Copyright 2001 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""This module implements Gotoh's O(MN) algorithm for aligning
sequences.

Gotoh O.  An improved algorithm for matching biological sequences.  J
Mol Biol. 1982 Dec 15;162(3):705-8.

Please refer to the original paper for technical details on the
working of the algorithm.

Functions:
align_global   Do a global alignment between two sequences.
align_local    Do a local alignment between two sequences.

identity_match    Calculate match scores based on identity.
dictionary_match  Calculate match scores from a dictionary lookup.

"""
from support import *

# These constants indicate the path of an optimal alignment through a
# traceback matrix.
GO_D = 1   # go diagonally
GO_P = 2   # go up
GO_Q = 4   # go left

def align_global(sequenceA, sequenceB, open, extend, match_fn=None, **keywds):
    """align_global(sequenceA, sequenceB, open, extend[, match_fn]) ->
    alignments

    Do a global alignment on sequenceA and sequenceB.  open and extend
    are the gap penalties to apply to the sequences.  They should be
    negative.  match_fn is a callback function that takes two
    characters and returns the score.  By default, will give a score
    of 1 for matches and 0 for mismatches.

    alignments is a list of tuples (seqA, seqB, score, begin, end).
    seqA and seqB are strings showing the alignment between the
    sequences.  score is the score of the alignment.  begin and end
    are indexes into seqA and seqB that indicate the where the
    alignment occurs.  Here, it should be the whole sequence.

    """
    if match_fn is not None:
        keywds['match_fn'] = match_fn
    params = [
        ('match_fn', identity_match(1, 0)),
        ('gap_penalty_A', (open, extend)),
        ('gap_penalty_B', (open, extend)),
        ('count_first', 0),
        ('global_alignment', 1),
        ('penalize_end_gaps', 1),
        ('gap_char', '-'),
        ]
    for name, default in params:
        keywds[name] = keywds.get(name, default)
    alignments = _align(sequenceA, sequenceB, **keywds)
    return alignments
    

def align_local(sequenceA, sequenceB, open, extend, match_fn=None, **keywds):
    """align_local(sequenceA, sequenceB, open, extend[, match_fn]) ->
    alignments

    Do a local alignment on sequenceA and sequenceB.  open and extend
    are the gap penalties to apply to the sequences.  They should be
    negative.  match_fn is a callback function that takes two
    characters and returns the score.  By default, will give a score
    of 1 for matches and 0 for mismatches.

    alignments is a list of tuples (seqA, seqB, score, begin, end).
    seqA and seqB are strings showing the alignment between the
    sequences.  score is the score of the alignment.  begin and end
    are indexes into seqA and seqB that indicate the local region
    where the alignment occurs.

    """
    if match_fn is not None:
        keywds['match_fn'] = match_fn
    params = [
        ('match_fn', identity_match(1, 0)),
        ('gap_penalty_A', (open, extend)),
        ('gap_penalty_B', (open, extend)),
        ('count_first', 0),
        ('global_alignment', 0),
        ('penalize_end_gaps', 0),
        ('gap_char', '-'),
        ]
    for name, default in params:
        keywds[name] = keywds.get(name, default)
    return _align(sequenceA, sequenceB, **keywds)

def _align(sequenceA, sequenceB, match_fn, gap_penalty_A, gap_penalty_B,
           count_first, global_alignment, penalize_end_gaps, gap_char):
    # - count_first is whether to count an extension on the first gap
    # open.  If false, a gap of 1 is penalized by "open".  Otherwise,
    # it's penalized "open + extend".
    # - gap_penalty_A/B should be tuples of (open, extend) penalties.
    # - global_alignment is whether the alignment should be global.
    # - penalize_end_gaps is whether the gaps at the beginning and
    # end of the alignments should be counted.
    # - gap_char is the character to use for gaps.

    # First, check to make sure the gap penalties are reasonable.
    open_A, extend_A = gap_penalty_A
    open_B, extend_B = gap_penalty_B
    if open_A > 0 or extend_A > 0 or open_B > 0 or extend_B > 0:
        raise ValueError, "Gap penalties must be negative"

    score_matrix, direction_matrix = _make_score_matrix(
        sequenceA, sequenceB, open_A, extend_A, open_B, extend_B,
        count_first, penalize_end_gaps, match_fn)

    # Look for the proper starting point
    score, indexes = _find_best_score(
        score_matrix, sequenceA, sequenceB, 
        open_A, extend_A, open_B, extend_B,
        global_alignment, penalize_end_gaps, count_first)

##     for i in range(len(score_matrix)):
##         for j in range(len(score_matrix[i])):
##             print score_matrix[i][j],
##         print
##     print

##     for i in range(len(direction_matrix)):
##         for j in range(len(direction_matrix[i])):
##             print direction_matrix[i][j],
##         print
##     print
##     print score, indexes

    # Recover the alignments and return them.
    alignments = _recover_alignments(
        sequenceA, sequenceB, score, indexes,
        score_matrix, direction_matrix,
        global_alignment, gap_char)
    return alignments

def _make_score_matrix(sequenceA, sequenceB,
                       open_A, extend_A, open_B, extend_B,
                       count_first, penalize_end_gaps, match_fn):
    # Make a matrix sequenceA (down) x sequenceB (across)
    M, N = len(sequenceA), len(sequenceB)
    Dmatrix, Pmatrix, Qmatrix, direction_matrix = [], [], [], []
    for i in range(M):
        Dmatrix.append([0] * N)             # going diagonal, best score
        Pmatrix.append([0] * N)             # going up
        Qmatrix.append([0] * N)             # going left
        direction_matrix.append([0] * N)    # for the traceback
    
    # Cache some commonly used gap penalties.
    first_B_gap = calc_affine_penalty(1, open_B, extend_B, count_first)
    first_A_gap = calc_affine_penalty(1, open_A, extend_A, count_first)

    # Now fill in the score matrix.
    for i in range(M):
        for j in range(N):
            if i:
                Pscore = max(
                    Pmatrix[i-1][j] + extend_B,
                    Dmatrix[i-1][j] + first_B_gap
                    )
            else:
##                 if penalize_end_gaps:
##                     Pscore = calc_affine_penalty(j, open_B, extend_B,
##                                                  count_first)
##                 else:
##                     Pscore = 0
                Pscore = -10000
            if j:
                Qscore = max(
                    Qmatrix[i][j-1] + extend_A,
                    Dmatrix[i][j-1] + first_A_gap
                    )
            else:
##                 if penalize_end_gaps:
##                     Qscore = calc_affine_penalty(i, open_A, extend_A,
##                                                  count_first)
##                 else:
##                     Qscore = 0
                Qscore = -10000
            Dscore = match_fn(sequenceA[i], sequenceB[j])
            if i and j:
                Dscore += Dmatrix[i-1][j-1]
            elif penalize_end_gaps:
                if not i:
                    Dscore += calc_affine_penalty(
                        j, open_B, extend_B, count_first)
                else:
                    Dscore += calc_affine_penalty(
                        i, open_A, extend_A, count_first)
                
            Pscore_rint, Qscore_rint, Dscore_rint = map(
                rint, (Pscore, Qscore, Dscore))
            maxscore = max(Pscore, Qscore, Dscore)
            maxscore_rint = max(Pscore_rint, Qscore_rint, Dscore_rint)

            direction = direction_matrix[i][j]
            if Pscore_rint == maxscore_rint and i:
                direction |= GO_P
            if Qscore_rint == maxscore_rint and j:
                direction |= GO_Q
            if Dscore_rint == maxscore_rint:
                direction |= GO_D
            direction_matrix[i][j] = direction

            Pmatrix[i][j] = Pscore
            Qmatrix[i][j] = Qscore
            Dmatrix[i][j] = maxscore
    return Dmatrix, direction_matrix

def _find_best_score(score_matrix, sequenceA, sequenceB,
                     open_A, extend_A, open_B, extend_B,
                     global_alignment, penalize_end_gaps, count_first):
    if global_alignment:
        if penalize_end_gaps:
            gap_A_fn = affine_penalty(open_A, extend_A, count_first)
            gap_B_fn = affine_penalty(open_B, extend_B, count_first)
            score, indexes = find_global_best(sequenceA, sequenceB,
                score_matrix, 1, gap_A_fn, gap_B_fn)
        else:
            score, indexes = find_global_best(
                score_matrix, 0, None, None)
    else:
        score, indexes = find_local_best(score_matrix)
    return score, indexes

def _recover_alignments(sequenceA, sequenceB, score, indexes, 
                        score_matrix, direction_matrix,
                        global_alignment, gap_char):
    # Recover the alignments.  This is a recursive procedure, but it's
    # implemented here iteratively with a stack.
    tracebacks = []   # list of (seqA, seqB, score, begin, end)
    in_process = []   # list of ([same as tracebacks], row, col)

    # Initialize the in_process stack.
    for row, col in indexes:
        seqA, seqB = sequenceA[row+1:], sequenceB[col+1:]
        seqA, seqB = pad_until_equal(seqA, seqB, gap_char)
        if global_alignment:
            begin, end = 0, None
        else:
            begin, end = 0, -len(seqA)
            if not end:
                end = None
        in_process.append((seqA, seqB, score, begin, end, row, col))
        
    while in_process:
        seqA, seqB, score, begin, end, row, col = in_process.pop()
        if row < 0 and col < 0:
            # This one is done.  Put it into the tracebacks list and
            # continue.
            tracebacks.append((seqA, seqB, score, begin, end))
            continue
        if row < 0:
            nseqA = gap_char + seqA
            nseqB = sequenceB[col] + seqB
            if not global_alignment and begin < col+1:
                begin = col + 1
            in_process.append((nseqA, nseqB, score, begin, end, row, col-1))
        elif col < 0:
            nseqA = sequenceA[row] + seqA
            nseqB = gap_char + seqB
            if not global_alignment and begin < row+1:
                begin = row + 1
            in_process.append((nseqA, nseqB, score, begin, end, row-1, col))
        elif not global_alignment and score_matrix[row][col] <= 0:
            begin = max(row, col) + 1
            nseqA = sequenceA[:row+1] + seqA
            nseqB = sequenceB[:col+1] + seqB
            nseqA, nseqB = lpad_until_equal(nseqA, nseqB, gap_char)
            in_process.append((nseqA, nseqB, score, begin, end, -1, -1))
        else:
            # when I see GO_D, add the current chars to the alignment
            if direction_matrix[row][col] & GO_D:
                nseqA = sequenceA[row] + seqA
                nseqB = sequenceB[col] + seqB
                in_process.append(
                    (nseqA, nseqB, score, begin, end, row-1, col-1))
            # when I see GO_P, add the current seqA row, and a gap
            if direction_matrix[row][col] & GO_P:
                nseqA = sequenceA[row] + seqA
                nseqB = gap_char + seqB
                in_process.append(
                    (nseqA, nseqB, score, begin, end, row-1, col))
            # when I see GO_Q, add the current seqB col, and a gap
            if direction_matrix[row][col] & GO_Q:
                nseqA = gap_char + seqA
                nseqB = sequenceB[col] + seqB
                in_process.append(
                    (nseqA, nseqB, score, begin, end, row, col-1))
    return clean_alignments(tracebacks)
