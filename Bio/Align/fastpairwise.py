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
from types import *

from support import *



MAX_ALIGNMENTS = 1000


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
        ('score_cutoff_fn', None),
        ('forgiveness', 0),
        ]
    for name, default in params:
        keywds[name] = keywds.get(name, default)
    return _align(sequenceA, sequenceB, **keywds)

def _default_score_cutoff_fn(score, pos, all_scores):
    return rint(score) == rint(all_scores[0])

def _align(sequenceA, sequenceB, match_fn, gap_penalty_A, gap_penalty_B,
           count_first, global_alignment, penalize_end_gaps, gap_char,
           score_cutoff_fn, forgiveness):
    # - count_first is whether to count an extension on the first gap
    # open.  If false, a gap of 1 is penalized by "open".  Otherwise,
    # it's penalized "open + extend".
    # - gap_penalty_A/B should be tuples of (open, extend) penalties.
    # - global_alignment is whether the alignment should be global.
    # - penalize_end_gaps is whether the gaps at the beginning and
    # end of the alignments should be counted.
    # - gap_char is the character to use for gaps.
    # - score_cutoff_fn is a callback that determines.  It takes the
    # current score, the (row, col), and a list of all scores (sorted
    # from best to worst) and returns a boolean for whether to recover
    # this alignment.
    # - forgiveness is the some number of points the alignments are
    # allowed to deviate from the best score
    
    if not sequenceA or not sequenceB:
        return []
    # First, check to make sure the gap penalties are reasonable.
    open_A, extend_A = gap_penalty_A
    open_B, extend_B = gap_penalty_B
    if open_A > 0 or extend_A > 0 or open_B > 0 or extend_B > 0:
        raise ValueError, "Gap penalties must be negative"

    score_matrix, direction_matrix = _make_score_matrix_optimized(
        sequenceA, sequenceB, open_A, extend_A, open_B, extend_B,
        count_first, global_alignment, penalize_end_gaps, match_fn,
        forgiveness)

    # Look for the proper starting point
    starts = _find_start(score_matrix, sequenceA, sequenceB, 
                         open_A, extend_A, open_B, extend_B,
                         global_alignment, penalize_end_gaps, count_first)
    if score_cutoff_fn is None:
        score_cutoff_fn = _default_score_cutoff_fn
    scores = [x[0] for x in starts]
    scores.sort()
    scores.reverse()
    i = 0
    while i < len(starts):
        score, pos = starts[i]
        if not score_cutoff_fn(score, pos, scores):
            del starts[i]
        else:
            i += 1

    # Recover the alignments and return them.
    alignments = _recover_alignments(sequenceA, sequenceB, starts,
                                     score_matrix, direction_matrix,
                                     global_alignment, gap_char)
    return alignments

def _make_score_matrix_optimized(*args):
    sequenceA, sequenceB, match_fn = args[0], args[1], args[-2]
    
    # See if we can make assumptions about the inputs
    if type(sequenceA) is StringType and type(sequenceB) is StringType and \
       type(match_fn) is InstanceType and \
       match_fn.__class__ is identity_match:
        fn = _make_score_matrix_faster
    else:
        fn = _make_score_matrix
    return fn(*args)

def _make_score_matrix_faster(*args):
    return _make_score_matrix(*args)

def _make_score_matrix(sequenceA, sequenceB,
                       open_A, extend_A, open_B, extend_B,
                       count_first, global_alignment, penalize_end_gaps,
                       match_fn, forgiveness):
    forgiveness_rint = rint(forgiveness)
    # Make a matrix sequenceA (down) x sequenceB (across)
    M, N = len(sequenceA), len(sequenceB)
    Dmatrix, Pmatrix, Qmatrix, direction_matrix = [], [], [], []
    for i in range(M):
        Dmatrix.append([0] * N)             # going diagonal, best score
        Pmatrix.append([0] * N)             # going up
        Qmatrix.append([0] * N)             # going left
        direction_matrix.append([0] * N)    # for the traceback
    
    # Cache some commonly used gap penalties.
    first_A_gap = calc_affine_penalty(1, open_A, extend_A, count_first)
    first_B_gap = calc_affine_penalty(1, open_B, extend_B, count_first)

    # Now fill in the score matrix.
    for i in range(M):
        for j in range(N):
            if i > 0 and i < M-1:
                Pscore = max(
                    Pmatrix[i-1][j] + extend_B,
                    Dmatrix[i-1][j] + first_B_gap
                    )
            else:
                # This is not an airtight solution.  This is meant to
                # disallow the P move.  However, if the user chooses
                # perverse parameters, this score may actually be
                # good.  Ideally, this should be replaced with
                # something out of bound, like None.  However, it
                # causes the code to get much more messy.
                Pscore = -10000
            if j > 0 and j < N-1:
                Qscore = max(
                    Qmatrix[i][j-1] + extend_A,
                    Dmatrix[i][j-1] + first_A_gap
                    )
            else:
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

            if not global_alignment:
                Pscore, Qscore, Dscore = \
                        max(Pscore, 0), max(Qscore, 0), max(Dscore, 0)

            Pscore_rint, Qscore_rint, Dscore_rint = map(
                rint, (Pscore, Qscore, Dscore))
            maxscore = max(Pscore, Qscore, Dscore)
            maxscore_rint = max(Pscore_rint, Qscore_rint, Dscore_rint)

            Pmatrix[i][j] = Pscore
            Qmatrix[i][j] = Qscore
            Dmatrix[i][j] = maxscore

            direction = 0
            if maxscore_rint-Pscore_rint <= forgiveness_rint and i:
                direction |= GO_P
            if maxscore_rint-Qscore_rint <= forgiveness_rint and j:
                direction |= GO_Q
            if maxscore_rint-Dscore_rint <= forgiveness_rint:
                direction |= GO_D
            direction_matrix[i][j] = direction

    return Dmatrix, direction_matrix

def _find_start(score_matrix, sequenceA, sequenceB,
                open_A, extend_A, open_B, extend_B,
                global_alignment, penalize_end_gaps, count_first):
    if global_alignment:
        if penalize_end_gaps:
            gap_A_fn = affine_penalty(open_A, extend_A, count_first)
            gap_B_fn = affine_penalty(open_B, extend_B, count_first)
            starts = find_global_start(sequenceA, sequenceB,
                                       score_matrix, 1, gap_A_fn, gap_B_fn)
        else:
            starts = find_global_start(sequenceA, sequenceB,
                                       score_matrix, 0, None, None)
    else:
        starts = find_local_start(score_matrix)
    return starts

##def _find_best_score(score_matrix, sequenceA, sequenceB,
##                     open_A, extend_A, open_B, extend_B,
##                     global_alignment, penalize_end_gaps, count_first):
##    if global_alignment:
##        if penalize_end_gaps:
##            gap_A_fn = affine_penalty(open_A, extend_A, count_first)
##            gap_B_fn = affine_penalty(open_B, extend_B, count_first)
##            score, indexes = find_global_best(sequenceA, sequenceB,
##                score_matrix, 1, gap_A_fn, gap_B_fn)
##        else:
##            score, indexes = find_global_best(sequenceA, sequenceB,
##                score_matrix, 0, None, None)
##    else:
##        score, indexes = find_local_best(score_matrix)
##    return score, indexes

def _print_matrix(matrix):
    nums = []
    for m in matrix:
        nums.extend(m)
    numstrs = ["%g" % x for x in nums]
    lens = map(len, numstrs)
    ndigits = max(lens)
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            print "%*g " % (ndigits, matrix[i][j]),
        print


def _recover_alignments(sequenceA, sequenceB, starts,
                        score_matrix, direction_matrix,
                        global_alignment, gap_char):
    # Recover the alignments.  This is a recursive procedure, but it's
    # implemented here iteratively with a stack.
    tracebacks = []   # list of (seqA, seqB, score, begin, end)
    in_process = []   # list of ([same as tracebacks], row, col)

    #print "SCORE"; _print_matrix(score_matrix)
    #print "DIRECTION"; _print_matrix(direction_matrix)
    #print "INDEXES", indexes

    # sequenceA and sequenceB may be sequences, including strings,
    # lists, or list-like objects.  In order to preserve the type of
    # the object, we need to use slices on the sequences instead of
    # indexes.  For example, sequenceA[row] may return a type that's
    # not compatible with sequenceA, e.g. if sequenceA is a list and
    # sequenceA[row] is a string.  Thus, avoid using indexes and use
    # slices, e.g. sequenceA[row:row+1].  Assume that client-defined
    # sequence classes preserve these semantics.

    # Initialize the in_process stack.
    for score, (row, col) in starts:
        seqA, seqB = sequenceA[row+1:], sequenceB[col+1:]
        seqA, seqB = pad_until_equal(seqA, seqB, gap_char)
        if global_alignment:
            # Set the begin to 0.  I don't know how long the alignment
            # is going to be right now, so set end to None.  Figure
            # out the length of the alignment later.
            begin, end = 0, None
        else:
            # Since the length of the alignment is unknown, set an
            # offset from the end of the alignment.
            begin, end = 0, -len(seqA)
            if not end:
                end = None
        x = seqA, seqB, score, begin, end, row, col, []
        in_process.append(x)

    # This is a potentially exponential process.  This is particularly
    # prone to happen when there are no gap penalties, because there
    # are many ways things can be unaligned.  To solve this, we can
    # cache the paths that lead to different alignments.
    path_cache = {}  # (row, col) -> list of (seqA, seqB, begin, taillen)

    while in_process and len(tracebacks) < MAX_ALIGNMENTS:
        seqA, seqB, score, begin, end, row, col, path = in_process.pop()
        path = path + [((row, col), len(seqA))]
        if path_cache.has_key((row, col)):
            path.pop()     # This is already cached.  Don't repeat.
            # Need to set the begin properly, if necessary.
            if not global_alignment:
                if row < 0:
                    begin = col + 1
                elif col < 0:
                    begin = row + 1
            for cseqA, cseqB, cbegin, taillen in path_cache[(row, col)]:
                nbegin = begin
                nseqA = cseqA[:-taillen] + seqA
                nseqB = cseqB[:-taillen] + seqB
                if nbegin < cbegin:
                    nbegin = cbegin
                x = nseqA, nseqB, score, nbegin, end, -1, -1, path
                in_process.append(x)
        elif row < 0 and col < 0:
            # This one is done.  Put it into the tracebacks list and
            # continue.
            tracebacks.append((seqA, seqB, score, begin, end))
            path.pop()   # Don't cache the last step.
            # Update the cache.
            for p, taillen in path:
                if not path_cache.has_key(p):
                    path_cache[p] = []
                path_cache[p].append((seqA, seqB, begin, taillen))
            continue
        elif row < 0:
            nseqA = gap_char + seqA
            nseqB = sequenceB[col:col+1] + seqB
            if not global_alignment:
                begin = col + 1
            x = nseqA, nseqB, score, begin, end, row, col-1, path
            in_process.append(x)
        elif col < 0:
            nseqA = sequenceA[row:row+1] + seqA
            nseqB = gap_char + seqB
            if not global_alignment:
                begin = row + 1
            x = nseqA, nseqB, score, begin, end, row-1, col, path
            in_process.append(x)
        elif not global_alignment and score_matrix[row][col] <= 0:
            # Stop the alignment early.
            begin = max(row, col) + 1
            nseqA = sequenceA[:row+1] + seqA
            nseqB = sequenceB[:col+1] + seqB
            nseqA, nseqB = lpad_until_equal(nseqA, nseqB, gap_char)
            x = nseqA, nseqB, score, begin, end, -1, -1, path
            in_process.append(x)
        else:
            # when I see GO_D, add the current chars to the alignment
            if direction_matrix[row][col] & GO_D:
                nseqA = sequenceA[row:row+1] + seqA
                nseqB = sequenceB[col:col+1] + seqB
                x = nseqA, nseqB, score, begin, end, row-1, col-1, path
                in_process.append(x)
            # when I see GO_P, add the current seqA row, and a gap
            if direction_matrix[row][col] & GO_P:
                nseqA = sequenceA[row:row+1] + seqA
                nseqB = gap_char + seqB
                x = nseqA, nseqB, score, begin, end, row-1, col, path
                in_process.append(x)
            # when I see GO_Q, add the current seqB col, and a gap
            if direction_matrix[row][col] & GO_Q:
                nseqA = gap_char + seqA
                nseqB = sequenceB[col:col+1] + seqB
                x = nseqA, nseqB, score, begin, end, row, col-1, path
                in_process.append(x)
    #print "TRACEBACKS", tracebacks
    return clean_alignments(tracebacks)


# Try and load C implementations of functions.  If I can't,
# then just ignore and use the pure python implementations.
try:
    #raise ImportError
    import cfastpairwise
except ImportError:
    pass
else:
    _make_score_matrix = cfastpairwise._make_score_matrix
    _make_score_matrix_faster = cfastpairwise._make_score_matrix_faster
    
