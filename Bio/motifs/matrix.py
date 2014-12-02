# Copyright 2013 by Michiel de Hoon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Implementation of frequency (count) matrices, position-weight matrices,
and position-specific scoring matrices.
"""

import math

from Bio._py3k import range

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

# Hack for Python 2.5, isnan and isinf were new in Python 2.6
try:
    from math import isnan as _isnan
except ImportError:
    def _isnan(value):
        # This is tricky due to cross platform float differences
        if str(value).lower() == "nan":
            return True
        return value != value
try:
    from math import isinf as _isinf
except ImportError:
    def _isinf(value):
        # This is tricky due to cross platform float differences
        if str(value).lower().endswith("inf"):
            return True
        return False
# Hack for Python 2.5 on Windows:
try:
    _nan = float("nan")
except ValueError:
    _nan = 1e1000 / 1e1000


class GenericPositionMatrix(dict):

    def __init__(self, alphabet, values):
        self.length = None
        for letter in alphabet.letters:
            if self.length is None:
                self.length = len(values[letter])
            elif self.length!=len(values[letter]):
                raise Exception("data has inconsistent lengths")
            self[letter] = list(values[letter])
        self.alphabet = alphabet
        self._letters = sorted(self.alphabet.letters)

    def __str__(self):
        words = ["%6d" % i for i in range(self.length)]
        line = "   " + " ".join(words)
        lines = [line]
        for letter in self._letters:
            words = ["%6.2f" % value for value in self[letter]]
            line = "%c: " % letter + " ".join(words)
            lines.append(line)
        text = "\n".join(lines) + "\n"
        return text

    def __getitem__(self, key):
        if isinstance(key, tuple):
            if len(key)==2:
                key1, key2 = key
                if isinstance(key1, slice):
                    start1, stop1, stride1 = key1.indices(len(self._letters))
                    indices1 = range(start1, stop1, stride1)
                    letters1 = [self._letters[i] for i in indices1]
                    dim1 = 2
                elif isinstance(key1, int):
                    letter1 = self._letters[key1]
                    dim1 = 1
                elif isinstance(key1, tuple):
                    letters1 = [self._letters[i] for i in key1]
                    dim1 = 2
                elif isinstance(key1, str):
                    if len(key1)==1:
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
                if dim1==1 and dim2==1:
                    return dict.__getitem__(self, letter1)[index2]
                elif dim1==1 and dim2==2:
                    values = dict.__getitem__(self, letter1)
                    return tuple(values[index2] for index2 in indices2)
                elif dim1==2 and dim2==1:
                    d = {}
                    for letter1 in letters1:
                        d[letter1] = dict.__getitem__(self, letter1)[index2]
                    return d
                else:
                    d = {}
                    for letter1 in letters1:
                        values = dict.__getitem__(self, letter1)
                        d[letter1] = [values[index2] for index2 in indices2]
                    if sorted(letters1)==self._letters:
                        return self.__class__(self.alphabet, d)
                    else:
                        return d
            elif len(key)==1:
                key = key[0]
            else:
                raise KeyError("keys should be 1- or 2-dimensional")
        if isinstance(key, slice):
            start, stop, stride = key.indices(len(self._letters))
            indices = range(start, stop, stride)
            letters = [self._letters[i] for i in indices]
            dim = 2
        elif isinstance(key, int):
            letter = self._letters[key]
            dim = 1
        elif isinstance(key, tuple):
            letters = [self._letters[i] for i in key]
            dim = 2
        elif isinstance(key, str):
            if len(key)==1:
                letter = key
                dim = 1
            else:
                raise KeyError(key)
        else:
            raise KeyError("Cannot understand key %s", str(key))
        if dim==1:
            return dict.__getitem__(self, letter)
        elif dim==2:
            d = {}
            for letter in letters:
                d[letter] = dict.__getitem__(self, letter)
            return d
        else:
            raise RuntimeError("Should not get here")

    @property
    def consensus(self):
        """Returns the consensus sequence."""
        sequence = ""
        for i in range(self.length):
            try:
                maximum = float("-inf")
            except ValueError:
                # On Python 2.5 or older that was handled in C code,
                # and failed on Windows XP 32bit
                maximum = - 1E400
            for letter in self.alphabet.letters:
                count = self[letter][i]
                if count > maximum:
                    maximum = count
                    sequence_letter = letter
            sequence += sequence_letter
        return Seq(sequence, self.alphabet)

    @property
    def anticonsensus(self):
        sequence = ""
        for i in range(self.length):
            try:
                minimum = float("inf")
            except ValueError:
                # On Python 2.5 or older that was handled in C code,
                # and failed on Windows XP 32bit
                minimum = 1E400
            for letter in self.alphabet.letters:
                count = self[letter][i]
                if count < minimum:
                    minimum = count
                    sequence_letter = letter
            sequence += sequence_letter
        return Seq(sequence, self.alphabet)

    @property
    def degenerate_consensus(self):
        # Following the rules adapted from
        # D. R. Cavener: "Comparison of the consensus sequence flanking
        # translational start sites in Drosophila and vertebrates."
        # Nucleic Acids Research 15(4): 1353-1361. (1987).
        # The same rules are used by TRANSFAC.
        degenerate_nucleotide = {
            'A': 'A',
            'C': 'C',
            'G': 'G',
            'T': 'T',
            'AC': 'M',
            'AG': 'R',
            'AT': 'W',
            'CG': 'S',
            'CT': 'Y',
            'GT': 'K',
            'ACG': 'V',
            'ACT': 'H',
            'AGT': 'D',
            'CGT': 'B',
            'ACGT': 'N',
        }
        sequence = ""
        for i in range(self.length):
            def get(nucleotide):
                return self[nucleotide][i]
            nucleotides = sorted(self, key=get, reverse=True)
            counts = [self[c][i] for c in nucleotides]
            # Follow the Cavener rules:
            if counts[0] >= sum(counts[1:]) and counts[0] >= 2*counts[1]:
                key = nucleotides[0]
            elif 4*sum(counts[:2]) > 3*sum(counts):
                key = "".join(sorted(nucleotides[:2]))
            elif counts[3]==0:
                key = "".join(sorted(nucleotides[:3]))
            else:
                key = "ACGT"
            nucleotide = degenerate_nucleotide[key]
            sequence += nucleotide
        return Seq(sequence, alphabet=IUPAC.ambiguous_dna)

    @property
    def gc_content(self):
        """Compute the fraction GC content."""
        alphabet = self.alphabet
        gc_total = 0.0
        total = 0.0
        for i in range(self.length):
            for letter in alphabet.letters:
                if letter in 'CG':
                    gc_total += self[letter][i]
                total += self[letter][i]
        return gc_total / total

    def reverse_complement(self):
        values = {}
        values["A"] = self["T"][::-1]
        values["T"] = self["A"][::-1]
        values["G"] = self["C"][::-1]
        values["C"] = self["G"][::-1]
        alphabet = self.alphabet
        return self.__class__(alphabet, values)


class FrequencyPositionMatrix(GenericPositionMatrix):

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
            for letter in self.alphabet.letters:
                counts[letter] = [0.0] * self.length
        elif isinstance(pseudocounts, dict):
            for letter in self.alphabet.letters:
                counts[letter] = [float(pseudocounts[letter])] * self.length
        else:
            for letter in self.alphabet.letters:
                counts[letter] = [float(pseudocounts)] * self.length
        for i in range(self.length):
            for letter in self.alphabet.letters:
                counts[letter][i] += self[letter][i]
        # Actual normalization is done in the PositionWeightMatrix initializer
        return PositionWeightMatrix(self.alphabet, counts)


class PositionWeightMatrix(GenericPositionMatrix):

    def __init__(self, alphabet, counts):
        GenericPositionMatrix.__init__(self, alphabet, counts)
        for i in range(self.length):
            total = sum(float(self[letter][i]) for letter in alphabet.letters)
            for letter in alphabet.letters:
                self[letter][i] /= total
        for letter in alphabet.letters:
            self[letter] = tuple(self[letter])

    def log_odds(self, background=None):
        """Returns the Position-Specific Scoring Matrix.

        The Position-Specific Scoring Matrix (PSSM) contains the log-odds
        scores computed from the probability matrix and the background
        probabilities. If the background is None, a uniform background
        distribution is assumed.
        """
        values = {}
        alphabet = self.alphabet
        if background is None:
            background = dict.fromkeys(self._letters, 1.0)
        else:
            background = dict(background)
        total = sum(background.values())
        for letter in alphabet.letters:
            background[letter] /= total
            values[letter] = []
        for i in range(self.length):
            for letter in alphabet.letters:
                b = background[letter]
                if b > 0:
                    p = self[letter][i]
                    if p > 0:
                        logodds = math.log(p/b, 2)
                    else:
                        # TODO - Ensure this has unittest coverage!
                        try:
                            logodds = float("-inf")
                        except ValueError:
                            # On Python 2.5 or older that was handled in C code,
                            # and failed on Windows XP 32bit
                            logodds = - 1E400
                else:
                    p = self[letter][i]
                    if p > 0:
                        logodds = float("inf")
                    else:
                        logodds = _nan
                values[letter].append(logodds)
        pssm = PositionSpecificScoringMatrix(alphabet, values)
        return pssm


class PositionSpecificScoringMatrix(GenericPositionMatrix):

    def calculate(self, sequence):
        """Returns the PWM score for a given sequence for all positions.

        Notes:

         - the sequence can only be a DNA sequence
         - the search is performed only on one strand
         - if the sequence and the motif have the same length, a single
           number is returned
         - otherwise, the result is a one-dimensional list or numpy array
        """
        # TODO - Code itself tolerates ambiguous bases (as NaN).
        if not isinstance(self.alphabet, IUPAC.IUPACUnambiguousDNA):
            raise ValueError("PSSM has wrong alphabet: %s - Use only with DNA motifs"
                                 % self.alphabet)
        if not isinstance(sequence.alphabet, IUPAC.IUPACUnambiguousDNA):
            raise ValueError("Sequence has wrong alphabet: %r - Use only with DNA sequences"
                                 % sequence.alphabet)

        # TODO - Force uppercase here and optimise switch statement in C
        # by assuming upper case?
        sequence = str(sequence)
        m = self.length
        n = len(sequence)

        scores = []
        # check if the fast C code can be used
        try:
            import _pwm
        except ImportError:
            # use the slower Python code otherwise
            # The C code handles mixed case so Python version must too:
            sequence = sequence.upper()
            for i in range(n-m+1):
                score = 0.0
                for position in range(m):
                    letter = sequence[i+position]
                    try:
                        score += self[letter][position]
                    except KeyError:
                        score = _nan
                        break
                scores.append(score)
        else:
            # get the log-odds matrix into a proper shape
            # (each row contains sorted (ACGT) log-odds values)
            logodds = [[self[letter][i] for letter in "ACGT"] for i in range(m)]
            scores = _pwm.calculate(sequence, logodds)
        if len(scores)==1:
            return scores[0]
        else:
            return scores

    def search(self, sequence, threshold=0.0, both=True):
        """Find hits with PWM score above given threshold.

        A generator function, returning found hits in the given sequence
        with the pwm score higher than the threshold.
        """
        sequence = sequence.upper()
        n = len(sequence)
        m = self.length
        if both:
            rc = self.reverse_complement()
        for position in range(0, n-m+1):
            s = sequence[position:position+m]
            score = self.calculate(s)
            if score > threshold:
                yield (position, score)
            if both:
                score = rc.calculate(s)
                if score > threshold:
                    yield (position-n, score)

    @property
    def max(self):
        """Maximal possible score for this motif.

        returns the score computed for the consensus sequence.
        """
        score = 0.0
        letters = self._letters
        for position in range(0, self.length):
            score += max(self[letter][position] for letter in letters)
        return score

    @property
    def min(self):
        """Minimal possible score for this motif.

        returns the score computed for the anticonsensus sequence.
        """
        score = 0.0
        letters = self._letters
        for position in range(0, self.length):
            score += min(self[letter][position] for letter in letters)
        return score

    @property
    def gc_content(self):
        raise Exception("Cannot compute the %GC composition of a PSSM")

    def mean(self, background=None):
        """Expected value of the score of a motif."""
        if background is None:
            background = dict.fromkeys(self._letters, 1.0)
        else:
            background = dict(background)
        total = sum(background.values())
        for letter in self._letters:
            background[letter] /= total
        sx = 0.0
        for i in range(self.length):
            for letter in self._letters:
                logodds = self[letter, i]
                if _isnan(logodds):
                    continue
                if _isinf(logodds) and logodds < 0:
                    continue
                b = background[letter]
                p = b * math.pow(2, logodds)
                sx += p * logodds
        return sx

    def std(self, background=None):
        """Standard deviation of the score of a motif."""
        if background is None:
            background = dict.fromkeys(self._letters, 1.0)
        else:
            background = dict(background)
        total = sum(background.values())
        for letter in self._letters:
            background[letter] /= total
        variance = 0.0
        for i in range(self.length):
            sx = 0.0
            sxx = 0.0
            for letter in self._letters:
                logodds = self[letter, i]
                if _isnan(logodds):
                    continue
                if _isinf(logodds) and logodds < 0:
                    continue
                b = background[letter]
                p = b * math.pow(2, logodds)
                sx += p*logodds
                sxx += p*logodds*logodds
            sxx -= sx*sx
            variance += sxx
        variance = max(variance, 0) # to avoid roundoff problems
        return math.sqrt(variance)

    def dist_pearson(self, other):
        """Return the similarity score based on pearson correlation for the given motif against self.

        We use the Pearson's correlation of the respective probabilities.
        """
        if self.alphabet != other.alphabet:
            raise ValueError("Cannot compare motifs with different alphabets")

        max_p=-2
        for offset in range(-self.length+1, other.length):
            if offset<0:
                p = self.dist_pearson_at(other, -offset)
            else:  # offset>=0
                p = other.dist_pearson_at(self, offset)
            if max_p<p:
                max_p=p
                max_o=-offset
        return 1-max_p, max_o

    def dist_pearson_at(self, other, offset):
        letters = self._letters
        sx = 0.0   # \sum x
        sy = 0.0   # \sum y
        sxx = 0.0  # \sum x^2
        sxy = 0.0  # \sum x \cdot y
        syy = 0.0  # \sum y^2
        norm=max(self.length, offset+other.length)*len(letters)
        for pos in range(min(self.length-offset, other.length)):
            xi = [self[letter, pos+offset] for letter in letters]
            yi = [other[letter, pos] for letter in letters]
            sx += sum(xi)
            sy += sum(yi)
            sxx += sum(x*x for x in xi)
            sxy += sum(x*y for x, y in zip(xi, yi))
            syy += sum(y*y for y in yi)
        sx /= norm
        sy /= norm
        sxx /= norm
        sxy /= norm
        syy /= norm
        numerator = sxy - sx*sy
        denominator = math.sqrt((sxx-sx*sx)*(syy-sy*sy))
        return numerator/denominator

    def distribution(self, background=None, precision=10**3):
        """calculate the distribution of the scores at the given precision."""
        from .thresholds import ScoreDistribution
        if background is None:
            background = dict.fromkeys(self._letters, 1.0)
        else:
            background = dict(background)
        total = sum(background.values())
        for letter in self._letters:
            background[letter] /= total
        return ScoreDistribution(precision=precision, pssm=self, background=background)
