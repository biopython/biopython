# Copyright 2000, 2004 by Brad Chapman.
# Revisions copyright 2010-2013, 2015-2018 by Peter Cock.
# All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Code for dealing with sequence alignments.

One of the most important things in this module is the MultipleSeqAlignment
class, used in the Bio.AlignIO module.

"""

import warnings

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Please install numpy if you want to use Bio.Align. "
        "See http://www.numpy.org/"
    ) from None

from Bio import BiopythonDeprecationWarning
from Bio.Align import _aligners
from Bio.Align._base import Alignment, Alignments, LazyAlignments
from Bio.Seq import Seq, MutableSeq, reverse_complement

# Import errors may occur here if a compiled aligners.c file
# (_aligners.pyd or _aligners.so) is missing or if the user is
# importing from within the Biopython source tree, see PR #2007:
# https://github.com/biopython/biopython/pull/2007


class PairwiseAlignments(LazyAlignments):
    """Implements an iterator over pairwise alignments returned by the aligner.

    This class also supports indexing, which is fast for increasing indices,
    but may be slow for random access of a large number of alignments.

    Note that pairwise aligners can return an astronomical number of alignments,
    even for relatively short sequences, if they align poorly to each other. We
    therefore recommend to first check the number of alignments, accessible as
    len(alignments), which can be calculated quickly even if the number of
    alignments is very large.
    """

    def __init__(self, seqA, seqB, score, paths):
        """Initialize a new PairwiseAlignments object.

        Arguments:
         - seqA  - The first sequence, as a plain string, without gaps.
         - seqB  - The second sequence, as a plain string, without gaps.
         - score - The alignment score.
         - paths - An iterator over the paths in the traceback matrix;
                   each path defines one alignment.

        You would normally obtain a PairwiseAlignments object by calling
        aligner.align(seqA, seqB), where aligner is a PairwiseAligner object.
        """
        self.sequences = [seqA, seqB]
        self.score = score
        self.paths = paths
        self._index = -1

    def __len__(self):
        """Return the number of alignments."""
        return len(self.paths)

    def __getitem__(self, index):
        if isinstance(index, slice):
            self._load()
            alignments = Alignments()
            items = self.__getitem__(index)
            alignments.extend(items)
            return alignments
        if index == self._index:
            return self.alignment
        if index < self._index:
            self.paths.reset()
            self._index = -1
        while self._index < index:
            try:
                alignment = next(self)
            except StopIteration:
                raise IndexError("index out of range") from None
        return alignment

    def __iter__(self):
        self.paths.reset()
        self._index = -1
        return self

    def __next__(self):
        path = next(self.paths)
        self._index += 1
        coordinates = numpy.array(path)
        alignment = Alignment(self.sequences, coordinates)
        alignment.score = self.score
        self.alignment = alignment
        return alignment

    def clear(self):
        # black complains if we inherit the parent's docstring  # fmt: skip
        del self.paths
        self._index = 0
        self.__class__ = Alignments


class PairwiseAligner(_aligners.PairwiseAligner):
    """Performs pairwise sequence alignment using dynamic programming.

    This provides functions to get global and local alignments between two
    sequences.  A global alignment finds the best concordance between all
    characters in two sequences.  A local alignment finds just the
    subsequences that align the best.

    To perform a pairwise sequence alignment, first create a PairwiseAligner
    object.  This object stores the match and mismatch scores, as well as the
    gap scores.  Typically, match scores are positive, while mismatch scores
    and gap scores are negative or zero.  By default, the match score is 1,
    and the mismatch and gap scores are zero.  Based on the values of the gap
    scores, a PairwiseAligner object automatically chooses the appropriate
    alignment algorithm (the Needleman-Wunsch, Smith-Waterman, Gotoh, or
    Waterman-Smith-Beyer global or local alignment algorithm).

    Calling the "score" method on the aligner with two sequences as arguments
    will calculate the alignment score between the two sequences.
    Calling the "align" method on the aligner with two sequences as arguments
    will return a generator yielding the alignments between the two
    sequences.

    Some examples:

    >>> from Bio import Align
    >>> aligner = Align.PairwiseAligner()
    >>> alignments = aligner.align("TACCG", "ACG")
    >>> for alignment in sorted(alignments):
    ...     print("Score = %.1f:" % alignment.score)
    ...     print(alignment)
    ...
    Score = 3.0:
    TACCG
    -|-||
    -A-CG
    <BLANKLINE>
    Score = 3.0:
    TACCG
    -||-|
    -AC-G
    <BLANKLINE>

    Specify the aligner mode as local to generate local alignments:

    >>> aligner.mode = 'local'
    >>> alignments = aligner.align("TACCG", "ACG")
    >>> for alignment in sorted(alignments):
    ...     print("Score = %.1f:" % alignment.score)
    ...     print(alignment)
    ...
    Score = 3.0:
    TACCG
     |-||
     A-CG
    <BLANKLINE>
    Score = 3.0:
    TACCG
     ||-|
     AC-G
    <BLANKLINE>

    Do a global alignment.  Identical characters are given 2 points,
    1 point is deducted for each non-identical character.

    >>> aligner.mode = 'global'
    >>> aligner.match_score = 2
    >>> aligner.mismatch_score = -1
    >>> for alignment in aligner.align("TACCG", "ACG"):
    ...     print("Score = %.1f:" % alignment.score)
    ...     print(alignment)
    ...
    Score = 6.0:
    TACCG
    -||-|
    -AC-G
    <BLANKLINE>
    Score = 6.0:
    TACCG
    -|-||
    -A-CG
    <BLANKLINE>

    Same as above, except now 0.5 points are deducted when opening a
    gap, and 0.1 points are deducted when extending it.

    >>> aligner.open_gap_score = -0.5
    >>> aligner.extend_gap_score = -0.1
    >>> aligner.target_end_gap_score = 0.0
    >>> aligner.query_end_gap_score = 0.0
    >>> for alignment in aligner.align("TACCG", "ACG"):
    ...     print("Score = %.1f:" % alignment.score)
    ...     print(alignment)
    ...
    Score = 5.5:
    TACCG
    -|-||
    -A-CG
    <BLANKLINE>
    Score = 5.5:
    TACCG
    -||-|
    -AC-G
    <BLANKLINE>

    The alignment function can also use known matrices already included in
    Biopython:

    >>> from Bio.Align import substitution_matrices
    >>> aligner = Align.PairwiseAligner()
    >>> aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    >>> alignments = aligner.align("KEVLA", "EVL")
    >>> alignments = list(alignments)
    >>> print("Number of alignments: %d" % len(alignments))
    Number of alignments: 1
    >>> alignment = alignments[0]
    >>> print("Score = %.1f" % alignment.score)
    Score = 13.0
    >>> print(alignment)
    KEVLA
    -|||-
    -EVL-
    <BLANKLINE>

    You can also set the value of attributes directly during construction
    of the PairwiseAligner object by providing them as keyword arguments:

    >>> aligner = Align.PairwiseAligner(mode='global', match_score=2, mismatch_score=-1)
    >>> for alignment in aligner.align("TACCG", "ACG"):
    ...     print("Score = %.1f:" % alignment.score)
    ...     print(alignment)
    ...
    Score = 6.0:
    TACCG
    -||-|
    -AC-G
    <BLANKLINE>
    Score = 6.0:
    TACCG
    -|-||
    -A-CG
    <BLANKLINE>

    """

    def __init__(self, **kwargs):
        """Initialize a new PairwiseAligner with the keyword arguments as attributes.

        Loops over the keyword arguments and sets them as attributes on the object.
        """
        super().__init__()
        for name, value in kwargs.items():
            setattr(self, name, value)

    def __setattr__(self, key, value):
        if key not in dir(_aligners.PairwiseAligner):
            # To prevent confusion, don't allow users to create new attributes.
            # On CPython, __slots__ can be used for this, but currently
            # __slots__ does not behave the same way on PyPy at least.
            raise AttributeError("'PairwiseAligner' object has no attribute '%s'" % key)
        _aligners.PairwiseAligner.__setattr__(self, key, value)

    def align(self, seqA, seqB, strand="+"):
        """Return the alignments of two sequences using PairwiseAligner."""
        if isinstance(seqA, (Seq, MutableSeq)):
            sA = bytes(seqA)
        else:
            sA = seqA
        if strand == "+":
            sB = seqB
        else:  # strand == "-":
            sB = reverse_complement(seqB, inplace=False)
        if isinstance(sB, (Seq, MutableSeq)):
            sB = bytes(sB)
        score, paths = _aligners.PairwiseAligner.align(self, sA, sB, strand)
        alignments = PairwiseAlignments(seqA, seqB, score, paths)
        return alignments

    def score(self, seqA, seqB, strand="+"):
        """Return the alignments score of two sequences using PairwiseAligner."""
        if isinstance(seqA, (Seq, MutableSeq)):
            seqA = bytes(seqA)
        if strand == "-":
            seqB = reverse_complement(seqB, inplace=False)
        if isinstance(seqB, (Seq, MutableSeq)):
            seqB = bytes(seqB)
        return _aligners.PairwiseAligner.score(self, seqA, seqB, strand)

    def __getstate__(self):
        state = {
            "wildcard": self.wildcard,
            "target_internal_open_gap_score": self.target_internal_open_gap_score,
            "target_internal_extend_gap_score": self.target_internal_extend_gap_score,
            "target_left_open_gap_score": self.target_left_open_gap_score,
            "target_left_extend_gap_score": self.target_left_extend_gap_score,
            "target_right_open_gap_score": self.target_right_open_gap_score,
            "target_right_extend_gap_score": self.target_right_extend_gap_score,
            "query_internal_open_gap_score": self.query_internal_open_gap_score,
            "query_internal_extend_gap_score": self.query_internal_extend_gap_score,
            "query_left_open_gap_score": self.query_left_open_gap_score,
            "query_left_extend_gap_score": self.query_left_extend_gap_score,
            "query_right_open_gap_score": self.query_right_open_gap_score,
            "query_right_extend_gap_score": self.query_right_extend_gap_score,
            "mode": self.mode,
        }
        if self.substitution_matrix is None:
            state["match_score"] = self.match_score
            state["mismatch_score"] = self.mismatch_score
        else:
            state["substitution_matrix"] = self.substitution_matrix
        return state

    def __setstate__(self, state):
        self.wildcard = state["wildcard"]
        self.target_internal_open_gap_score = state["target_internal_open_gap_score"]
        self.target_internal_extend_gap_score = state[
            "target_internal_extend_gap_score"
        ]
        self.target_left_open_gap_score = state["target_left_open_gap_score"]
        self.target_left_extend_gap_score = state["target_left_extend_gap_score"]
        self.target_right_open_gap_score = state["target_right_open_gap_score"]
        self.target_right_extend_gap_score = state["target_right_extend_gap_score"]
        self.query_internal_open_gap_score = state["query_internal_open_gap_score"]
        self.query_internal_extend_gap_score = state["query_internal_extend_gap_score"]
        self.query_left_open_gap_score = state["query_left_open_gap_score"]
        self.query_left_extend_gap_score = state["query_left_extend_gap_score"]
        self.query_right_open_gap_score = state["query_right_open_gap_score"]
        self.query_right_extend_gap_score = state["query_right_extend_gap_score"]
        self.mode = state["mode"]
        substitution_matrix = state.get("substitution_matrix")
        if substitution_matrix is None:
            self.match_score = state["match_score"]
            self.mismatch_score = state["mismatch_score"]
        else:
            self.substitution_matrix = substitution_matrix


class PairwiseAlignment(Alignment):
    """Represents a pairwise sequence alignment.

    Internally, the pairwise alignment is stored as the path through
    the traceback matrix, i.e. a tuple of pairs of indices corresponding
    to the vertices of the path in the traceback matrix.
    """

    def __init__(self, target, query, path, score):
        """Initialize a new PairwiseAlignment object.

        Arguments:
         - target  - The first sequence, as a plain string, without gaps.
         - query   - The second sequence, as a plain string, without gaps.
         - path    - The path through the traceback matrix, defining an
                     alignment.
         - score   - The alignment score.

        You would normally obtain a PairwiseAlignment object by iterating
        over a PairwiseAlignments object.
        """
        warnings.warn(
            "The PairwiseAlignment class is deprecated; please use the "
            "Alignment class instead.  Note that the coordinates attribute of "
            "an Alignment object is a numpy array and the transpose of the "
            "path attribute of a PairwiseAlignment object.",
            BiopythonDeprecationWarning,
        )
        sequences = [target, query]
        coordinates = numpy.array(path).transpose()
        super().__init__(sequences, coordinates)
        self.score = score


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
