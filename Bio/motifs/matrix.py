# Copyright 2013 by Michiel de Hoon.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Support for various forms of sequence motif matrices.

Implementation of frequency (count) matrices, position-weight matrices,
and position-specific scoring matrices.
"""

import math

try:
    import numpy as np
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.motifs.matrix."
    )

from Bio.Seq import Seq

from . import _pwm


def _calculate(score_dict, sequence, m):
    """Calculate scores using C code (PRIVATE)."""
    n = len(sequence)
    # Create the numpy arrays here; the C module then does not rely on numpy
    # Use a float32 for the scores array to save space
    scores = np.empty(n - m + 1, np.float32)
    logodds = np.array(
        [[score_dict[letter][i] for letter in "ACGT"] for i in range(m)], float
    )
    _pwm.calculate(sequence, logodds, scores)
    return scores


class GenericPositionMatrix(dict):
    """Base class for the support of position matrix operations."""

    def __init__(self, alphabet, values):
        """Initialize the class."""
        self.length = None
        for letter in alphabet:
            if self.length is None:
                self.length = len(values[letter])
            elif self.length != len(values[letter]):
                raise Exception("data has inconsistent lengths")
            self[letter] = list(values[letter])
        self.alphabet = alphabet

    def __str__(self):
        """Return a string containing nucleotides and counts of the alphabet in the Matrix."""
        words = ["%6d" % i for i in range(self.length)]
        line = "   " + " ".join(words)
        lines = [line]
        for letter in self.alphabet:
            words = ["%6.2f" % value for value in self[letter]]
            line = "%c: " % letter + " ".join(words)
            lines.append(line)
        text = "\n".join(lines) + "\n"
        return text

    def __getitem__(self, key):
        """Return the position matrix of index key."""
        if isinstance(key, tuple):
            if len(key) == 2:
                key1, key2 = key
                if isinstance(key1, slice):
                    start1, stop1, stride1 = key1.indices(len(self.alphabet))
                    indices1 = range(start1, stop1, stride1)
                    letters1 = [self.alphabet[i] for i in indices1]
                    dim1 = 2
                elif isinstance(key1, int):
                    letter1 = self.alphabet[key1]
                    dim1 = 1
                elif isinstance(key1, tuple):
                    letters1 = [self.alphabet[i] for i in key1]
                    dim1 = 2
                elif isinstance(key1, str):
                    if len(key1) == 1:
                        letter1 = key1
                        dim1 = 1
                    else:
                        raise KeyError(key1)
                else:
                    raise KeyError("Cannot understand key %s", str(key1))
                if isinstance(key2, slice):
                    start2, stop2, stride2 = key2.indices(self.length)
                    indices2 = range(start2, stop2, stride2)
                    dim2 = 2
                elif isinstance(key2, int):
                    index2 = key2
                    dim2 = 1
                else:
                    raise KeyError("Cannot understand key %s", str(key2))
                if dim1 == 1 and dim2 == 1:
                    return dict.__getitem__(self, letter1)[index2]
                elif dim1 == 1 and dim2 == 2:
                    values = dict.__getitem__(self, letter1)
                    return tuple(values[index2] for index2 in indices2)
                elif dim1 == 2 and dim2 == 1:
                    d = {}
                    for letter1 in letters1:
                        d[letter1] = dict.__getitem__(self, letter1)[index2]
                    return d
                else:
                    d = {}
                    for letter1 in letters1:
                        values = dict.__getitem__(self, letter1)
                        d[letter1] = [values[_] for _ in indices2]
                    if sorted(letters1) == self.alphabet:
                        return self.__class__(self.alphabet, d)
                    else:
                        return d
            elif len(key) == 1:
                key = key[0]
            else:
                raise KeyError("keys should be 1- or 2-dimensional")
        if isinstance(key, slice):
            start, stop, stride = key.indices(len(self.alphabet))
            indices = range(start, stop, stride)
            letters = [self.alphabet[i] for i in indices]
            dim = 2
        elif isinstance(key, int):
            letter = self.alphabet[key]
            dim = 1
        elif isinstance(key, tuple):
            letters = [self.alphabet[i] for i in key]
            dim = 2
        elif isinstance(key, str):
            if len(key) == 1:
                letter = key
                dim = 1
            else:
                raise KeyError(key)
        else:
            raise KeyError("Cannot understand key %s", str(key))
        if dim == 1:
            return dict.__getitem__(self, letter)
        elif dim == 2:
            d = {}
            for letter in letters:
                d[letter] = dict.__getitem__(self, letter)
            return d
        else:
            raise RuntimeError("Should not get here")

    @property
    def consensus(self):
        """Return the consensus sequence."""
        sequence = ""
        for i in range(self.length):
            maximum = -math.inf
            for letter in self.alphabet:
                count = self[letter][i]
                if count > maximum:
                    maximum = count
                    sequence_letter = letter
            sequence += sequence_letter
        return Seq(sequence)

    @property
    def anticonsensus(self):
        """Return the anticonsensus sequence."""
        sequence = ""
        for i in range(self.length):
            minimum = math.inf
            for letter in self.alphabet:
                count = self[letter][i]
                if count < minimum:
                    minimum = count
                    sequence_letter = letter
            sequence += sequence_letter
        return Seq(sequence)

    @property
    def degenerate_consensus(self):
        """Return the degenerate consensus sequence."""
        # Following the rules adapted from
        # D. R. Cavener: "Comparison of the consensus sequence flanking
        # translational start sites in Drosophila and vertebrates."
        # Nucleic Acids Research 15(4): 1353-1361. (1987).
        # The same rules are used by TRANSFAC.
        degenerate_nucleotide = {
            "A": "A",
            "C": "C",
            "G": "G",
            "T": "T",
            "AC": "M",
            "AG": "R",
            "AT": "W",
            "CG": "S",
            "CT": "Y",
            "GT": "K",
            "ACG": "V",
            "ACT": "H",
            "AGT": "D",
            "CGT": "B",
            "ACGT": "N",
        }
        sequence = ""
        for i in range(self.length):

            def get(nucleotide):
                return self[nucleotide][i]

            nucleotides = sorted(self, key=get, reverse=True)
            counts = [self[c][i] for c in nucleotides]
            # Follow the Cavener rules:
            if counts[0] > sum(counts[1:]) and counts[0] > 2 * counts[1]:
                key = nucleotides[0]
            elif 4 * sum(counts[:2]) > 3 * sum(counts):
                key = "".join(sorted(nucleotides[:2]))
            elif counts[3] == 0:
                key = "".join(sorted(nucleotides[:3]))
            else:
                key = "ACGT"
            nucleotide = degenerate_nucleotide.get(key, key)
            sequence += nucleotide
        return Seq(sequence)

    @property
    def gc_content(self):
        """Compute the fraction GC content."""
        alphabet = self.alphabet
        gc_total = 0.0
        total = 0.0
        for i in range(self.length):
            for letter in alphabet:
                if letter in "CG":
                    gc_total += self[letter][i]
                total += self[letter][i]
        return gc_total / total

    def reverse_complement(self):
        """Compute reverse complement."""
        values = {}
        if self.alphabet == "ACGU":
            values["A"] = self["U"][::-1]
            values["U"] = self["A"][::-1]
        else:
            values["A"] = self["T"][::-1]
            values["T"] = self["A"][::-1]
        values["G"] = self["C"][::-1]
        values["C"] = self["G"][::-1]
        alphabet = self.alphabet
        return self.__class__(alphabet, values)


class FrequencyPositionMatrix(GenericPositionMatrix):
    """Class for the support of frequency calculations on the Position Matrix."""

    def normalize(self, pseudocounts=None):
        """Create and return a position-weight matrix by normalizing the counts matrix.

        If pseudocounts is None (default), no pseudocounts are added
        to the counts.

        If pseudocounts is a number, it is added to the counts before
        calculating the position-weight matrix.

        Alternatively, the pseudocounts can be a dictionary with a key
        for each letter in the alphabet associated with the motif.
        """
        counts = {}
        if pseudocounts is None:
            for letter in self.alphabet:
                counts[letter] = [0.0] * self.length
        elif isinstance(pseudocounts, dict):
            for letter in self.alphabet:
                counts[letter] = [float(pseudocounts[letter])] * self.length
        else:
            for letter in self.alphabet:
                counts[letter] = [float(pseudocounts)] * self.length
        for i in range(self.length):
            for letter in self.alphabet:
                counts[letter][i] += self[letter][i]
        # Actual normalization is done in the PositionWeightMatrix initializer
        return PositionWeightMatrix(self.alphabet, counts)


class PositionWeightMatrix(GenericPositionMatrix):
    """Class for the support of weight calculations on the Position Matrix."""

    def __init__(self, alphabet, counts):
        """Initialize the class."""
        GenericPositionMatrix.__init__(self, alphabet, counts)
        for i in range(self.length):
            total = sum(float(self[letter][i]) for letter in alphabet)
            for letter in alphabet:
                self[letter][i] /= total
        for letter in alphabet:
            self[letter] = tuple(self[letter])

    def log_odds(self, background=None):
        """Return the Position-Specific Scoring Matrix.

        The Position-Specific Scoring Matrix (PSSM) contains the log-odds
        scores computed from the probability matrix and the background
        probabilities. If the background is None, a uniform background
        distribution is assumed.
        """
        values = {}
        alphabet = self.alphabet
        if background is None:
            background = dict.fromkeys(self.alphabet, 1.0)
        else:
            background = dict(background)
        total = sum(background.values())
        for letter in alphabet:
            background[letter] /= total
            values[letter] = []
        for i in range(self.length):
            for letter in alphabet:
                b = background[letter]
                if b > 0:
                    p = self[letter][i]
                    if p > 0:
                        logodds = math.log(p / b, 2)
                    else:
                        logodds = -math.inf
                else:
                    p = self[letter][i]
                    if p > 0:
                        logodds = math.inf
                    else:
                        logodds = math.nan
                values[letter].append(logodds)
        pssm = PositionSpecificScoringMatrix(alphabet, values)
        return pssm


class PositionSpecificScoringMatrix(GenericPositionMatrix):
    """Class for the support of Position Specific Scoring Matrix calculations."""

    def calculate(self, sequence):
        """Return the PWM score for a given sequence for all positions.

        Notes:
         - the sequence can only be a DNA sequence
         - the search is performed only on one strand
         - if the sequence and the motif have the same length, a single
           number is returned
         - otherwise, the result is a one-dimensional list or numpy array

        """
        # TODO - Code itself tolerates ambiguous bases (as NaN).
        if sorted(self.alphabet) != ["A", "C", "G", "T"]:
            raise ValueError(
                "PSSM has wrong alphabet: %s - Use only with DNA motifs" % self.alphabet
            )

        # NOTE: The C code handles mixed case input as this could be large
        # (e.g. contig or chromosome), so requiring it be all upper or lower
        # case would impose an overhead to allocate the extra memory.
        sequence = str(sequence)
        m = self.length
        scores = _calculate(self, sequence, m)

        if len(scores) == 1:
            return scores[0]
        else:
            return scores

    def search(self, sequence, threshold=0.0, both=True, chunksize=10 ** 6):
        """Find hits with PWM score above given threshold.

        A generator function, returning found hits in the given sequence
        with the pwm score higher than the threshold.
        """
        sequence = sequence.upper()
        seq_len = len(sequence)
        motif_l = self.length
        chunk_starts = np.arange(0, seq_len, chunksize)
        if both:
            rc = self.reverse_complement()
        for chunk_start in chunk_starts:
            subseq = sequence[chunk_start : chunk_start + chunksize + motif_l - 1]
            pos_scores = self.calculate(subseq)
            pos_ind = pos_scores >= threshold
            pos_positions = np.where(pos_ind)[0] + chunk_start
            pos_scores = pos_scores[pos_ind]
            if both:
                neg_scores = rc.calculate(subseq)
                neg_ind = neg_scores >= threshold
                neg_positions = np.where(neg_ind)[0] + chunk_start
                neg_scores = neg_scores[neg_ind]
            else:
                neg_positions = np.empty((0), dtype=int)
                neg_scores = np.empty((0), dtype=int)
            chunk_positions = np.append(pos_positions, neg_positions - seq_len)
            chunk_scores = np.append(pos_scores, neg_scores)
            order = np.argsort(np.append(pos_positions, neg_positions))
            chunk_positions = chunk_positions[order]
            chunk_scores = chunk_scores[order]
            yield from zip(chunk_positions, chunk_scores)

    @property
    def max(self):
        """Maximal possible score for this motif.

        returns the score computed for the consensus sequence.
        """
        score = 0.0
        letters = self.alphabet
        for position in range(0, self.length):
            score += max(self[letter][position] for letter in letters)
        return score

    @property
    def min(self):
        """Minimal possible score for this motif.

        returns the score computed for the anticonsensus sequence.
        """
        score = 0.0
        letters = self.alphabet
        for position in range(0, self.length):
            score += min(self[letter][position] for letter in letters)
        return score

    @property
    def gc_content(self):
        """Compute the GC-ratio."""
        raise Exception("Cannot compute the %GC composition of a PSSM")

    def mean(self, background=None):
        """Return expected value of the score of a motif."""
        if background is None:
            background = dict.fromkeys(self.alphabet, 1.0)
        else:
            background = dict(background)
        total = sum(background.values())
        for letter in self.alphabet:
            background[letter] /= total
        sx = 0.0
        for i in range(self.length):
            for letter in self.alphabet:
                logodds = self[letter, i]
                if math.isnan(logodds):
                    continue
                if math.isinf(logodds) and logodds < 0:
                    continue
                b = background[letter]
                p = b * math.pow(2, logodds)
                sx += p * logodds
        return sx

    def std(self, background=None):
        """Return standard deviation of the score of a motif."""
        if background is None:
            background = dict.fromkeys(self.alphabet, 1.0)
        else:
            background = dict(background)
        total = sum(background.values())
        for letter in self.alphabet:
            background[letter] /= total
        variance = 0.0
        for i in range(self.length):
            sx = 0.0
            sxx = 0.0
            for letter in self.alphabet:
                logodds = self[letter, i]
                if math.isnan(logodds):
                    continue
                if math.isinf(logodds) and logodds < 0:
                    continue
                b = background[letter]
                p = b * math.pow(2, logodds)
                sx += p * logodds
                sxx += p * logodds * logodds
            sxx -= sx * sx
            variance += sxx
        variance = max(variance, 0)  # to avoid roundoff problems
        return math.sqrt(variance)

    def dist_pearson(self, other):
        """Return the similarity score based on pearson correlation for the given motif against self.

        We use the Pearson's correlation of the respective probabilities.
        """
        if self.alphabet != other.alphabet:
            raise ValueError("Cannot compare motifs with different alphabets")

        max_p = -2
        for offset in range(-self.length + 1, other.length):
            if offset < 0:
                p = self.dist_pearson_at(other, -offset)
            else:  # offset>=0
                p = other.dist_pearson_at(self, offset)
            if max_p < p:
                max_p = p
                max_o = -offset
        return 1 - max_p, max_o

    def dist_pearson_at(self, other, offset):
        """Return the similarity score based on pearson correlation at the given offset."""
        letters = self.alphabet
        sx = 0.0  # \sum x
        sy = 0.0  # \sum y
        sxx = 0.0  # \sum x^2
        sxy = 0.0  # \sum x \cdot y
        syy = 0.0  # \sum y^2
        norm = max(self.length, offset + other.length) * len(letters)
        for pos in range(min(self.length - offset, other.length)):
            xi = [self[letter, pos + offset] for letter in letters]
            yi = [other[letter, pos] for letter in letters]
            sx += sum(xi)
            sy += sum(yi)
            sxx += sum(x * x for x in xi)
            sxy += sum(x * y for x, y in zip(xi, yi))
            syy += sum(y * y for y in yi)
        sx /= norm
        sy /= norm
        sxx /= norm
        sxy /= norm
        syy /= norm
        numerator = sxy - sx * sy
        denominator = math.sqrt((sxx - sx * sx) * (syy - sy * sy))
        return numerator / denominator

    def distribution(self, background=None, precision=10 ** 3):
        """Calculate the distribution of the scores at the given precision."""
        from .thresholds import ScoreDistribution

        if background is None:
            background = dict.fromkeys(self.alphabet, 1.0)
        else:
            background = dict(background)
        total = sum(background.values())
        for letter in self.alphabet:
            background[letter] /= total
        return ScoreDistribution(precision=precision, pssm=self, background=background)
