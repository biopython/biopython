# Copyright 2003-2009 by Bartek Wilczynski.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Implementation of sequence motifs (PRIVATE).
"""
from Bio.Seq import Seq
from Bio.SubsMat import FreqTable
from Bio.Alphabet import IUPAC
import math

import warnings
from Bio import BiopythonExperimentalWarning


class GenericPositionMatrix(dict):

    def __init__(self, alphabet, values):
        self.length = None
        for letter in alphabet.letters:
            if self.length is None:
                self.length = len(values[letter])
            elif self.length!=len(values[letter]):
                raise Exception("Inconsistent lengths found in dictionary")
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
        """Returns the consensus sequence.
        """
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
        return Seq(sequence, alphabet = IUPAC.ambiguous_dna)

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
        """
        create and return a position-weight matrix by normalizing the counts matrix.

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
        for i in xrange(self.length):
            for letter in self.alphabet.letters:
                counts[letter][i] += self[letter][i]
        # Actual normalization is done in the PositionWeightMatrix initializer
        return PositionWeightMatrix(self.alphabet, counts)


class PositionWeightMatrix(GenericPositionMatrix):

    def __init__(self, alphabet, counts):
        GenericPositionMatrix.__init__(self, alphabet, counts)
        for i in xrange(self.length):
            total = sum([float(self[letter][i]) for letter in alphabet.letters])
            for letter in alphabet.letters:
                self[letter][i] /= total
        for letter in alphabet.letters:
            self[letter] = tuple(self[letter])

    def log_odds(self, background=None):
        """
        returns the Position-Specific Scoring Matrix.

        The Position-Specific Scoring Matrix (PSSM) contains the log-odds
        scores computed from the probability matrix and the background
        probabilities. If the background is None, a uniform background
        distribution is assumed.
        """
        values = {}
        alphabet = self.alphabet
        if background is None:
            background = {}
            for letter in alphabet.letters:
                background[letter] = 1.0
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
                        #TODO - Ensure this has unittest coverage!
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
                        logodds = float("nan")
                values[letter].append(logodds)
        pssm = PositionSpecificScoringMatrix(alphabet, values)
        pssm._background = background
        return pssm


class PositionSpecificScoringMatrix(GenericPositionMatrix):

    def calculate(self, sequence):
        """
        returns the PWM score for a given sequence for all positions.

        - the sequence can only be a DNA sequence
        - the search is performed only on one strand
        - if the sequence and the motif have the same length, a single
          number is returned
        - otherwise, the result is a one-dimensional list or numpy array
        """
        if self.alphabet!=IUPAC.unambiguous_dna:
            raise ValueError("Wrong alphabet! Use only with DNA motifs")
        if sequence.alphabet!=IUPAC.unambiguous_dna:
            raise ValueError("Wrong alphabet! Use only with DNA sequences")

        sequence = str(sequence)
        m = self.length
        n = len(sequence)

        scores = []
        # check if the fast C code can be used
        try:
            import _pwm
        except ImportError:
            # use the slower Python code otherwise
            for i in xrange(n-m+1):
                score = 0.0
                for position in xrange(m):
                    letter = sequence[i+position]
                    score += self[letter][position]
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
        """
        a generator function, returning found hits in a given sequence with the pwm score higher than the threshold
        """
        sequence = sequence.upper()
        n = len(sequence)
        m = self.length
        if both:
            rc = self.reverse_complement()
        for position in xrange(0,n-m+1):
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
        for position in xrange(0,self.length):
            score += max([self[letter][position] for letter in letters])
        return score

    @property
    def min(self):
        """Minimal possible score for this motif.

        returns the score computed for the anticonsensus sequence.
        """
        score = 0.0
        letters = self._letters
        for position in xrange(0,self.length):
            score += min([self[letter][position] for letter in letters])
        return score

    @property
    def mean(self):
        """Expected value of the score of a motif.

        returns None if the expected value is undefined"""
        if self._background is None:
            return None
        background = self._background
        sx = 0.0
        for i in range(self.length):
            for letter in self._letters:
                logodds = self[letter,i]
                b = background[letter]
                p = b * math.pow(2,logodds)
                sx += p * logodds
        return sx
        

    @property
    def std(self):
        """Standard deviation of the score of a motif.

        returns None if the standard deviation is undefined"""
        if self._background is None:
            return None
        background = self._background
        variance = 0.0
        for i in range(self.length):
            sx = 0.0
            sxx = 0.0
            for letter in self._letters:
                logodds = self[letter,i]
                b = background[letter]
                p = b * math.pow(2,logodds)
                sx += p*logodds
                sxx += p*logodds*logodds
            sxx -= sx*sx
            variance += sxx
        variance = max(variance, 0) # to avoid roundoff problems
        return math.sqrt(variance)

    def dist_pearson(self, other):
        """
        return the similarity score based on pearson correlation for the given motif against self.

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
        return 1-max_p,max_o

    def dist_pearson_at(self, other, offset):
        letters = self._letters
        sx  = 0.0   # \sum x
        sy  = 0.0   # \sum y
        sxx = 0.0  # \sum x^2
        sxy = 0.0  # \sum x \cdot y
        syy = 0.0  # \sum y^2
        norm=max(self.length,offset+other.length)*len(letters)
        for pos in range(min(self.length-offset, other.length)):
            xi = [self[letter,pos+offset] for letter in letters]
            yi = [other[letter,pos] for letter in letters]
            sx += sum(xi)
            sy += sum(yi)
            sxx += sum([x*x for x in xi])
            sxy += sum([x*y for x,y in zip(xi,yi)])
            syy += sum([y*y for y in yi])
        sx /= norm
        sy /= norm
        sxx /= norm
        sxy /= norm
        syy /= norm
        numerator = sxy - sx*sy
        denominator = math.sqrt((sxx-sx*sx)*(syy-sy*sy))
        return numerator/denominator


class Motif(object):
    """
    A class representing sequence motifs.
    """
    def __init__(self, alphabet=None, instances=None, counts=None):
        self._pwm_is_current = False
        self._pwm = []
        self._log_odds_is_current = False
        self._log_odds = []
        self.length=None
        self.info=None
        self.name=""
        if counts is not None and instances is not None:
            raise Exception(ValueError,
                "Specify either instances or counts, don't specify both")
        elif counts is not None:
            warnings.warn("This is experimental code, and may change in future versions", BiopythonExperimentalWarning)
            if alphabet is None:
                alphabet = IUPAC.unambiguous_dna
            for letter in counts:
                length = len(counts[letter])
                if self.length is None:
                    self.length = length
                elif self.length!=length:
                    raise Exception("counts matrix has inconsistent lengths")
            self.instances = None
            self.counts = FrequencyPositionMatrix(alphabet, counts)
        elif instances is not None:
            warnings.warn("This is experimental code, and may change in future versions", BiopythonExperimentalWarning)
            self.instances = []
            for instance in instances:
                if alphabet is None:
                    alphabet=instance.alphabet
                elif alphabet != instance.alphabet:
                    raise ValueError("Alphabets are inconsistent")
                if self.length is None:
                    self.length = len(instance)
                elif self.length != len(instance):
                    message = "All instances should have the same length (%d found, %d expected)" % (len(instance), self.length)
                    raise ValueError(message)
                self.instances.append(instance)
            if alphabet.letters is None:
                # If we didn't get a meaningful alphabet from the instances,
                # assume it is DNA.
                alphabet = IUPAC.unambiguous_dna
            counts = {}
            for letter in alphabet.letters:
                counts[letter] = [0] * self.length
            for instance in self.instances:
                for position, letter in enumerate(instance):
                    counts[letter][position] += 1
            self.counts = FrequencyPositionMatrix(alphabet, counts)
        else:
            self.counts = None
            self.instances = None
            if alphabet is None:
                alphabet = IUPAC.unambiguous_dna
        self.alphabet = alphabet
        if self.length is None:
            self.__mask = ()
        else:
            self.__mask = (1,) * self.length
        self.background=dict((n, 1.0/len(self.alphabet.letters))
                             for n in self.alphabet.letters)
        self.beta=1.0

    @property
    def has_instances(self):
        """Legacy property, check if m.instances is None instead (DEPRECATED)."""
        from Bio import BiopythonDeprecationWarning
        warnings.warn("Instead of 'm.has_instances' use 'm.instances is not None'",
                      BiopythonDeprecationWarning)
        return self.instances is not None

    @property
    def has_counts(self):
        """Legacy property, check if m.counts is None instead (DEPRECATED)."""
        from Bio import BiopythonDeprecationWarning
        warnings.warn("Instead of 'm.has_counts' use 'm.counts is not None'",
                      BiopythonDeprecationWarning)
        return self.counts is not None

    def _check_length(self, len):
        warnings.warn("""\
This function is now obsolete, and will be deprecated and removed
in a future release of Biopython.
""", PendingDeprecationWarning)
        if self.length is None:
            self.length = len
        elif self.length != len:
            print "len",self.length,self.instances, len
            raise ValueError("You can't change the length of the motif")

    def _check_alphabet(self,alphabet):
        warnings.warn("""\
This function is now obsolete, and will be deprecated and removed
in a future release of Biopython.
""", PendingDeprecationWarning)
        if self.alphabet is None:
            self.alphabet=alphabet
        elif self.alphabet != alphabet:
                raise ValueError("Wrong Alphabet")

    def add_instance(self,instance):
        """
        adds new instance to the motif
        """
        warnings.warn("""\
This function is now obsolete, and will be deprecated and removed
in a future release of Biopython. Instead of adding instances to an existing
Motif object, please create a new Motif object.
""", PendingDeprecationWarning)
        self._check_alphabet(instance.alphabet)
        self._check_length(len(instance))
        if self.counts is not None:
            for i in range(self.length):
                let=instance[i]
                self.counts[let][i]+=1

        if self.instances is not None or self.counts is None:
            if self.instances is None:
                self.instances = []
            self.instances.append(instance)

        self._pwm_is_current = False
        self._log_odds_is_current = False

    def set_mask(self,mask):
        """
        sets the mask for the motif

        The mask should be a string containing asterisks in the position of significant columns and spaces in other columns
        """
        warnings.warn("""\
This function is now obsolete, and will be deprecated and removed
in a future release of Biopython. As a replacement, instead of
>>> motif.set_mask(mask)
please use
>>> motif.mask = mask
""", PendingDeprecationWarning)
        self._check_length(len(mask))
        self.__mask=[]
        for char in mask:
            if char=="*":
                self.__mask.append(1)
            elif char==" ":
                self.__mask.append(0)
            else:
                raise ValueError("Mask should contain only '*' or ' ' and not a '%s'"%char)
        self.__mask = tuple(self.__mask)

    def __get_mask(self):
        return self.__mask

    def __set_mask(self, mask):
        if mask is None:
            self.__mask = (1,) * self.length
        elif len(mask)!=self.length:
            raise ValueError("The length (%d) of the mask is inconsistent with the length (%d) of the motif", (len(mask), self.length))
        elif isinstance(mask, str):
            self.__mask=[]
            for char in mask:
                if char=="*":
                    self.__mask.append(1)
                elif char==" ":
                    self.__mask.append(0)
                else:
                    raise ValueError("Mask should contain only '*' or ' ' and not a '%s'"%char)
            self.__mask = tuple(self.__mask)
        else:
            self.__mask = tuple(int(bool(c)) for c in mask)

    mask = property(__get_mask, __set_mask)
    del __get_mask
    del __set_mask

    def pwm(self,laplace=True):
        """
        returns the PWM computed for the set of instances

        if laplace=True (default), pseudocounts equal to self.background multiplied by self.beta are added to all positions.
        """

        warnings.warn("""\
This function is now obsolete, and will be deprecated and removed
in a future release of Biopython. As a replacement, instead of
>>> motif.pwm()
use
>>> pwm = motif.counts.normalize()
See the documentation of motif.counts.normalize and pwm.log_odds for
details on treatment of pseudocounts and background probabilities.
""", PendingDeprecationWarning)
        if self._pwm_is_current:
            return self._pwm
        #we need to compute new pwm
        self._pwm = []
        for i in xrange(self.length):
            dict = {}
            #filling the dict with 0's
            for letter in self.alphabet.letters:
                if laplace:
                    dict[letter]=self.beta*self.background[letter]
                else:
                    dict[letter]=0.0
            if self.counts is not None:
                #taking the raw counts
                for letter in self.alphabet.letters:
                    dict[letter]+=self.counts[letter][i]
            elif self.instances is not None:
                #counting the occurences of letters in instances
                for seq in self.instances:
                    #dict[seq[i]]=dict[seq[i]]+1
                    try:
                        dict[seq[i]]+=1
                    except KeyError:  # we need to ignore non-alphabet letters
                        pass
            self._pwm.append(FreqTable.FreqTable(dict,FreqTable.COUNT,self.alphabet))
        self._pwm_is_current=1
        return self._pwm

    def log_odds(self,laplace=True):
        """
        returns the log odds matrix computed for the set of instances
        """
        warnings.warn("""\
This function is now obsolete, and will be deprecated and removed
in a future release of Biopython. As a replacement, instead of
>>> motif.log_odds()
use
>>> pwm = motif.counts.normalize()
>>> pssm = pwm.log_odds()
See the documentation of motif.counts.normalize and pwm.log_odds for
details on treatment of pseudocounts and background probabilities.
""", PendingDeprecationWarning)
        if self._log_odds_is_current:
            return self._log_odds
        #we need to compute new pwm
        self._log_odds = []
        pwm=self.pwm(laplace)
        for i in xrange(self.length):
            d = {}
            for a in self.alphabet.letters:
                    d[a]=math.log(pwm[i][a]/self.background[a],2)
            self._log_odds.append(d)
        self._log_odds_is_current=1
        return self._log_odds

    def ic(self):
        """Method returning the information content of a motif.
        """
        warnings.warn("""\
This function is now obsolete, and will be deprecated and removed
in a future release of Biopython. As a replacement, instead of
>>> motif.ic()
please use
>>> pwm = motif.counts.normalize()
>>> pwm.ic()
Please be aware though that by default, motif.counts.normalize()
does not use pseudocounts, while motif.ic() does. See the documentation
of motif.counts.normalize for more details.
""", PendingDeprecationWarning)
        res=0
        pwm=self.pwm()
        for i in range(self.length):
            res+=2
            for a in self.alphabet.letters:
                if pwm[i][a]!=0:
                    res+=pwm[i][a]*math.log(pwm[i][a],2)
        return res

    def exp_score(self,st_dev=False):
        """
        Computes expected score of motif's instance and its standard deviation
        """
        warnings.warn("""\
This function is now obsolete, and will be deprecated and removed
in a future release of Biopython. As a replacement, instead of
>>> motif.exp_score()
please use
>>> pwm = motif.counts.normalize()
>>> pwm.exp_score()
See the documentation of motif.counts.normalize for details on treatment of
pseudocounts.
""", PendingDeprecationWarning)
        exs=0.0
        var=0.0
        pwm=self.pwm()
        for i in range(self.length):
            ex1=0.0
            ex2=0.0
            for a in self.alphabet.letters:
                if pwm[i][a]!=0:
                    ex1+=pwm[i][a]*(math.log(pwm[i][a],2)-math.log(self.background[a],2))
                    ex2+=pwm[i][a]*(math.log(pwm[i][a],2)-math.log(self.background[a],2))**2
            exs+=ex1
            var+=ex2-ex1**2
        if st_dev:
            return exs,math.sqrt(var)
        else:
            return exs

    def search_instances(self,sequence):
        """
        a generator function, returning found positions of instances of the motif in a given sequence
        """
        if self.instances is None:
            raise ValueError("This motif has no instances")
        for pos in xrange(0,len(sequence)-self.length+1):
            for instance in self.instances:
                if str(instance) == str(sequence[pos:pos+self.length]):
                    yield(pos,instance)
                    break # no other instance will fit (we don't want to return multiple hits)

    def score_hit(self,sequence,position,normalized=0,masked=0):
        """
        give the pwm score for a given position
        """
        warnings.warn("""\
This function is now obsolete, and will be deprecated and removed
in a future release of Biopython. As a replacement, instead of
>>> motif.score_hit(sequence, position)
please use
>>> pwm = motif.counts.normalize()
>>> pssm = pwm.log_odds()
>>> s = sequence[position:positon+len(pssm)]
>>> pssm.calculate(s)
See the documentation of motif.counts.normalize() and pwm.log_odds
for details on the treatment of pseudocounts and background probabilities.
""", PendingDeprecationWarning)
        lo=self.log_odds()
        score = 0.0
        for pos in xrange(self.length):
            a = sequence[position+pos]
            if not masked or self.__mask[pos]:
                try:
                    score += lo[pos][a]
                except:
                    pass
        if normalized:
            if not masked:
                score/=self.length
            else:
                score/=len([x for x in self.__mask if x])
        return score

    def search_pwm(self,sequence,normalized=0,masked=0,threshold=0.0,both=True):
        """
        a generator function, returning found hits in a given sequence with the pwm score higher than the threshold
        """
        warnings.warn("""\
This function is now obsolete, and will be deprecated and removed
in a future release of Biopython. As a replacement, instead of
>>> motif.score_hit(sequence, position)
please use
>>> pwm = motif.counts.normalize()
>>> pssm = pwm.log_odds()
>>> pssm.search(sequence)
See the documentation of motif.counts.normalize() and pwm.log_odds
for details on treatment of pseudocounts and background probabilities.
""", PendingDeprecationWarning)
        raise Exception
        if both:
            rc = self.reverse_complement()

        sequence=sequence.upper()
        n = len(sequence)
        for pos in xrange(0,n-self.length+1):
            score = self.score_hit(sequence,pos,normalized,masked)
            if score > threshold:
                yield (pos,score)
            if both:
                rev_score = rc.score_hit(sequence,pos,normalized,masked)
                if rev_score > threshold:
                    yield (pos-n,rev_score)

    def dist_pearson(self, motif, masked = 0):
        """
        return the similarity score based on pearson correlation for the given motif against self.

        We use the Pearson's correlation of the respective probabilities.
        """
        warnings.warn("""\
This function is now obsolete, and will be deprecated and removed
in a future release of Biopython. As a replacement, instead of
>>> motif1.dist_pearson(motif2)
please use
>>> pwm1 = motif1.counts.normalize()
>>> pwm2 = motif2.counts.normalize()
>>> pwm1.dist_pearson(pwm2)
Please see the documentation of motif.counts.normalize and
pwm.dist_pearson for more details.
""", PendingDeprecationWarning)

        if self.alphabet != motif.alphabet:
            raise ValueError("Cannot compare motifs with different alphabets")

        max_p=-2
        for offset in range(-self.length+1,motif.length):
            if offset<0:
                p = self.dist_pearson_at(motif,-offset)
            else:  # offset>=0
                p = motif.dist_pearson_at(self,offset)

            if max_p<p:
                max_p=p
                max_o=-offset
        return 1-max_p,max_o

    def dist_pearson_at(self,motif,offset):
        warnings.warn("""\
This function is now obsolete, and will be deprecated and removed
in a future release of Biopython.""", PendingDeprecationWarning)
        sxx = 0  # \sum x^2
        sxy = 0  # \sum x \cdot y
        sx = 0   # \sum x
        sy = 0   # \sum y
        syy = 0  # \sum y^2
        norm=max(self.length,offset+motif.length)

        for pos in range(max(self.length,offset+motif.length)):
            for l in self.alphabet.letters:
                xi = self[pos][l]
                yi = motif[pos-offset][l]
                sx = sx + xi
                sy = sy + yi
                sxx = sxx + xi * xi
                syy = syy + yi * yi
                sxy = sxy + xi * yi

        norm *= len(self.alphabet.letters)
        s1 = (sxy - sx*sy*1.0/norm)
        s2 = (norm*sxx - sx*sx*1.0)*(norm*syy- sy*sy*1.0)
        return s1/math.sqrt(s2)

    def dist_product(self,other):
        """
        A similarity measure taking into account a product probability of generating overlapping instances of two motifs
        """
        warnings.warn("""\
This function is now obsolete, and will be deprecated and removed
in a future release of Biopython.""", PendingDeprecationWarning)
        max_p=0.0
        for offset in range(-self.length+1,other.length):
            if offset<0:
                p = self.dist_product_at(other,-offset)
            else:  # offset>=0
                p = other.dist_product_at(self,offset)
            if max_p<p:
                max_p=p
                max_o=-offset
        return 1-max_p/self.dist_product_at(self,0),max_o

    def dist_product_at(self,other,offset):
        warnings.warn("""\
This function is now obsolete, and will be deprecated and removed
in a future release of Biopython.""", PendingDeprecationWarning)
        s=0
        for i in range(max(self.length,offset+other.length)):
            f1=self[i]
            f2=other[i-offset]
            for n,b in self.background.iteritems():
                s+=b*f1[n]*f2[n]
        return s/i

    def dist_dpq(self,other):
        warnings.warn("""\
This function is now obsolete, and will be deprecated and removed
in a future release of Biopython.""", PendingDeprecationWarning)
        r"""Calculates the DPQ distance measure between motifs.

        It is calculated as a maximal value of DPQ formula (shown using LaTeX
        markup, familiar to mathematicians):

        \sqrt{\sum_{i=1}^{alignment.len()} \sum_{k=1}^alphabet.len() \
        \{ m1[i].freq(alphabet[k])*log_2(m1[i].freq(alphabet[k])/m2[i].freq(alphabet[k])) +
           m2[i].freq(alphabet[k])*log_2(m2[i].freq(alphabet[k])/m1[i].freq(alphabet[k]))
        }

        over possible non-spaced alignments of two motifs.  See this reference:

        D. M Endres and J. E Schindelin, "A new metric for probability
        distributions", IEEE transactions on Information Theory 49, no. 7
        (July 2003): 1858-1860.
        """

        min_d=float("inf")
        min_o=-1
        d_s=[]
        for offset in range(-self.length+1,other.length):
            #print "%2.3d"%offset,
            if offset<0:
                d = self.dist_dpq_at(other,-offset)
                overlap = self.length+offset
            else:  # offset>=0
                d = other.dist_dpq_at(self,offset)
                overlap = other.length-offset
            overlap = min(self.length,other.length,overlap)
            out = self.length+other.length-2*overlap
            #print d,1.0*(overlap+out)/overlap,d*(overlap+out)/overlap
            #d = d/(2*overlap)
            d = (d/(out+overlap))*(2*overlap+out)/(2*overlap)
            #print d
            d_s.append((offset,d))
            if min_d> d:
                min_d=d
                min_o=-offset
        return min_d,min_o  # ,d_s

    def dist_dpq_at(self,other,offset):
        warnings.warn("""\
This function is now obsolete, and will be deprecated and removed
in a future release of Biopython.""", PendingDeprecationWarning)
        """
        calculates the dist_dpq measure with a given offset.

        offset should satisfy 0<=offset<=len(self)
        """
        def dpq(f1,f2,alpha):
            s=0
            for n in alpha.letters:
                avg=(f1[n]+f2[n])/2
                s+=f1[n]*math.log(f1[n]/avg,2)+f2[n]*math.log(f2[n]/avg,2)
            return math.sqrt(s)

        s=0
        for i in range(max(self.length,offset+other.length)):
            f1=self[i]
            f2=other[i-offset]
            s+=dpq(f1,f2,self.alphabet)
        return s

    def _read(self,stream):
        """Reads the motif from the stream (in AlignAce format).

        the self.alphabet variable must be set beforehand.
        If the last line contains asterisks it is used for setting mask
        """
        warnings.warn("This function is now obsolete, and will be deprecated and removed in a future release of Biopython. As a replacement, please use Bio.Motif.parse instead.", PendingDeprecationWarning)

        while 1:
            ln = stream.readline()
            if "*" in ln:
                self.set_mask(ln.strip("\n\c"))
                break
            self.add_instance(Seq(ln.strip(),self.alphabet))

    def __str__(self,masked=False):
        """ string representation of a motif.
        """
        string = ""
        if self.instances is not None:
            for inst in self.instances:
                string += str(inst) + "\n"

        if masked:
            for i in xrange(self.length):
                if self.__mask[i]:
                    string += "*"
                else:
                    string += " "
            string += "\n"
        return string

    def __len__(self):
        """return the length of a motif

        Please use this method (i.e. invoke len(m)) instead of referring to m.length directly.
        """
        if self.length is None:
            return 0
        else:
            return self.length

    def _write(self,stream):
        """
        writes the motif to the stream
        """
        warnings.warn("""\
This function is now obsolete, and will be deprecated and removed
in a future release of Biopython.""")

        stream.write(self.__str__())

    def _to_fasta(self):
        """
        FASTA representation of motif
        """
        if self.instances is None:
            warnings.warn("""\
Creating simulated instances of a motif is now obsolete. This functionality
will be deprecated and removed in a future release of Biopython.""", PendingDeprecationWarning)
            alpha="".join(self.alphabet.letters)
            #col[i] is a column taken from aligned motif instances
            col=[]
            instances=[]
            s = sum(self.counts[nuc][0] for nuc in self.alphabet.letters)
            for i in range(self.length):
                col.append("")
                for n in self.alphabet.letters:
                    col[i] += n * self.counts[n][i]
                if len(col[i])<s:
                    warnings.warn(UserWarning, "WARNING, column too short (%d; expected %d)" % (len(col[i]), s))
                    col[i]+=(alpha*s)[:(s-len(col[i]))]
            #iterate over instances
            for i in range(s):
                instance = ""  # start with empty seq
                for j in range(self.length):  # iterate over positions
                    instance+=col[j][i]
                instance = Seq(instance, self.alphabet)
                instances.append(instance)
        else:
            instances = self.instances
        text = ""
        for i, instance in enumerate(instances):
            text += ">instance%d\n%s\n" % (i, instance)
        return text

    def reverse_complement(self):
        """
        Gives the reverse complement of the motif
        """
        if self.instances is not None:
            instances = []
            for instance in self.instances:
                instance = instance.reverse_complement()
                instances.append(instance)
            res = Motif(instances=instances)
        else:  # has counts
            alphabet = self.alphabet
            res = Motif(alphabet)
            res.counts={}
            res.counts["A"]=self.counts["T"][::-1]
            res.counts["T"]=self.counts["A"][::-1]
            res.counts["G"]=self.counts["C"][::-1]
            res.counts["C"]=self.counts["G"][::-1]
            res.length=self.length
        res.__mask = self.__mask[::-1]
        return res

    def _from_jaspar_pfm(self,stream,make_instances=False):
        """
        reads the motif from Jaspar .pfm file

        The instances are fake, but the pwm is accurate.
        """
        warnings.warn("This function is now obsolete, and will be deprecated and removed in a future release of Biopython. Please use the 'pfm' format instead of the 'jaspar-pfm' format", PendingDeprecationWarning)
        return self._from_horiz_matrix(stream,letters="ACGT",make_instances=make_instances)

    def _from_vert_matrix(self,stream,letters=None,make_instances=False):
        """reads a vertical count matrix from stream and fill in the counts.
        """
        warnings.warn("This function is now obsolete, and will be deprecated and removed in a future release of Biopython.", PendingDeprecationWarning)
        self.counts = {}
        if letters is None:
            letters=self.alphabet.letters
        self.length=0
        for i in letters:
            self.counts[i]=[]
        for ln in stream.readlines():
            rec=map(float,ln.strip().split())
            for k,v in zip(letters,rec):
                self.counts[k].append(v)
            self.length+=1
        self.__mask = (1,) * self.length
        if make_instances is True:
            self.make_instances_from_counts()
        return self

    def _from_horiz_matrix(self,stream,letters=None,make_instances=False):
        """reads a horizontal count matrix from stream and fill in the counts.
        """
        warnings.warn("This function is now obsolete, and will be deprecated and removed in a future release of Biopython.", PendingDeprecationWarning)
        if letters is None:
            letters=self.alphabet.letters
        self.counts = {}

        for i in letters:
            ln = stream.readline().strip().split()
            #if there is a letter in the beginning, ignore it
            if ln[0]==i:
                ln=ln[1:]
            #print ln
            try:
                self.counts[i]=map(int,ln)
            except ValueError:  # not integers
                self.counts[i]=map(float,ln)  # map(lambda s: int(100*float(s)),ln)
            #print counts[i]

        s = sum(self.counts[nuc][0] for nuc in letters)
        l = len(self.counts[letters[0]])
        self.length=l
        self.__mask = (1,) * l
        if make_instances is True:
            self.make_instances_from_counts()
        return self

    def make_instances_from_counts(self):
        """Creates "fake" instances for a motif created from a count matrix.

        In case the sums of counts are different for different columns, the
        shorter columns are padded with background.
        """
        warnings.warn("This function is now obsolete, and will be deprecated and removed in a future release of Biopython.", PendingDeprecationWarning)
        alpha="".join(self.alphabet.letters)
        #col[i] is a column taken from aligned motif instances
        col=[]
        self.instances=[]
        s = sum(map(lambda nuc: self.counts[nuc][0],self.alphabet.letters))
        for i in range(self.length):
            col.append("")
            for n in self.alphabet.letters:
                col[i] = col[i]+ (n*(self.counts[n][i]))
            if len(col[i])<s:
                print "WARNING, column too short",len(col[i]),s
                col[i]+=(alpha*s)[:(s-len(col[i]))]
            #print i,col[i]
        #iterate over instances
        for i in range(s):
            inst=""  # start with empty seq
            for j in range(self.length):  # iterate over positions
                inst+=col[j][i]
            #print i,inst
            inst=Seq(inst,self.alphabet)
            self.add_instance(inst)
        return self.instances

    def make_counts_from_instances(self):
        """Creates the count matrix for a motif with instances.

        """
        warnings.warn("This function is now obsolete, and will be deprecated and removed in a future release of Biopython.", PendingDeprecationWarning)
        raise Exception(remove)
        #make strings for "columns" of motifs
        #col[i] is a column taken from aligned motif instances
        counts={}
        for a in self.alphabet.letters:
            counts[a]=[]
        s = len(self.instances)
        for i in range(self.length):
            ci = dict((a,0) for a in self.alphabet.letters)
            for inst in self.instances:
                ci[inst[i]]+=1
            for a in self.alphabet.letters:
                counts[a].append(ci[a])
        self.counts=counts
        return counts

    def _from_jaspar_sites(self,stream):
        """
        reads the motif from Jaspar .sites file

        The instances and pwm are OK.
        """
        warnings.warn("This function is now obsolete, and will be deprecated and removed in a future release of Biopython. Please use the 'sites' format instead of the 'jaspar-sites' format", PendingDeprecationWarning)
        # Probably this should be in a separate submodule of Bio.Motif
        self.instances = []
        for line in stream:
            if not line.startswith(">"):
                break
            # line contains the header ">...."
            # now read the actual sequence
            line = stream.next()
            instance = ""
            for c in line.strip():
                if c==c.upper():
                    instance += c
            instance = Seq(instance, self.alphabet)
            length = len(instance)
            if self.length is None:
                self.length = length
            elif length!=self.length:
                raise ValueError("Inconsistent motif lengths found")
            self.instances.append(instance)
        self._pwm_is_current = False
        self._log_odds_is_current = False
        self.__mask = (1,) * length
        return self

    def __getitem__(self,index):
        """Returns the probability distribution over symbols at a given position, padding with background.

        If the requested index is out of bounds, the returned distribution comes from background.
        """
        warnings.warn("""\
This function is now obsolete, and will be deprecated and removed
in a future release of Biopython. Instead of
>>> motif[i]
please use
>>> pwm = motif.counts.normalize()
>>> pwm[:,i]
""", PendingDeprecationWarning)
        if index in range(self.length):
            return self.pwm()[index]
        else:
            return self.background

    class ConsensusSeq(Seq):
        # This is an ugly hack that allows us raise a warning if a user attempts
        # to call the method motif.consensus() instead of accessing the property
        # motif.consensus. It can be removed after consensus() as a method has
        # been removed from Biopython. Same thing for motif.anticonsensus.
        def __call__(self):
            warnings.warn("""\
Motif.consensus and Motif.anticonsensus are now properties instead of methods.
Please use yourmotif.consensus instead of yourmotif.consensus(), and
yourmotif.anticonsensus instead of yourmotif.anticonsensus()""",
                PendingDeprecationWarning)
            return self
    ConsensusSeq.__name__ = Seq.__name__

    @property
    def consensus(self):
        """Returns the consensus sequence.
        """
        sequence = self.counts.consensus
        return Motif.ConsensusSeq(str(sequence), sequence.alphabet)

    @property
    def anticonsensus(self):
        """returns the least probable pattern to be generated from this motif.
        """
        sequence = self.counts.anticonsensus
        return Motif.ConsensusSeq(str(sequence), sequence.alphabet)

    @property
    def degenerate_consensus(self):
        """Following the rules adapted from
D. R. Cavener: "Comparison of the consensus sequence flanking
translational start sites in Drosophila and vertebrates."
Nucleic Acids Research 15(4): 1353-1361. (1987).
The same rules are used by TRANSFAC."""
        return self.counts.degenerate_consensus

    def max_score(self):
        """Maximal possible score for this motif.

        returns the score computed for the consensus sequence.
        """
        warnings.warn("""\
This function is now deprecated. Instead of
>>> motif.max_score()
please use
>>> pwm = motif.counts.normalize()
>>> pssm = pwm.log_odds()
>>> pssm.max_score()
""",
                PendingDeprecationWarning)
        return self.score_hit(self.consensus,0)

    def min_score(self):
        """Minimal possible score for this motif.

        returns the score computed for the anticonsensus sequence.
        """
        warnings.warn("""\
This function is now deprecated. Instead of
>>> motif.min_score()
please use
>>> pwm = motif.counts.normalize()
>>> pssm = pwm.log_odds()
>>> pssm.min_score()
""",
                PendingDeprecationWarning)
        return self.score_hit(self.anticonsensus,0)

    def weblogo(self,fname,format="PNG",version="2.8.2", **kwds):
        """
        uses the Berkeley weblogo service to download and save a weblogo of
        itself

        requires an internet connection.
        The parameters from **kwds are passed directly to the weblogo server.

        Currently, this method uses WebLogo 2.8.2 by default. This is likely
        to change in a future release of Biopython in favor of WebLogo 3.3.
        To use WebLogo 3.3 now, use the 'version="3"' argument when calling
        this method.

        These are the arguments and their default values passed to
        WebLogo 2.8.2; see their website at http://weblogo.berkeley.edu
        for more information:

            'logowidth' : '18',
            'logoheight' : '5',
            'logounits' : 'cm',
            'kind' : 'AUTO',
            'firstnum' : "1",
            'command' : 'Create Logo',
            'smallsamplecorrection' : "on",
            'symbolsperline' : 32,
            'res' : '96',
            'res_units' : 'ppi',
            'antialias' : 'on',
            'title' : '',
            'barbits' : '',
            'xaxis': 'on',
            'xaxis_label' : '',
            'yaxis': 'on',
            'yaxis_label' : '',
            'showends' : 'on',
            'shrink' : '0.5',
            'fineprint' : 'on',
            'ticbits' : '1',
            'colorscheme' : 'DEFAULT',
            'color1' : 'black',
            'color2' : 'blue',
            'color3' : 'red',
            'color4' : 'black',
            'color5' : 'purple',
            'color6' : 'orange',

        These are the arguments and their default values passed to
        WebLogo 3.3; see their website at http://weblogo.threeplusone.com
        for more information:

            'stack_width' : 'medium',
            'stack_per_line' : '40',
            'alphabet' : 'alphabet_dna',
            'ignore_lower_case' : True,
            'unit_name' : "bits",
            'first_index' : '1',
            'logo_start' : '1',
            'logo_end': str(self.length),
            'composition' : "comp_auto",
            'percentCG' : '',
            'scale_width' : True,
            'show_errorbars' : True,
            'logo_title' : '',
            'logo_label' : '',
            'show_xaxis': True,
            'xaxis_label': '',
            'show_yaxis': True,
            'yaxis_label': '',
            'yaxis_scale': 'auto',
            'yaxis_tic_interval' : '1.0',
            'show_ends' : True,
            'show_fineprint' : True,
            'symbols0': '',
            'symbols1': '',
            'symbols2': '',
            'symbols3': '',
            'symbols4': '',
            'color0': '',
            'color1': '',
            'color2': '',
            'color3': '',
            'color4': '',

        """
        import urllib
        import urllib2
        version = str(version)
        if version.startswith("2"):
            warnings.warn("""\
Please use WebLogo 3 instead of the older WebLogo 2.8.2.
To use WebLogo 3, use the 'version="3"' argument when calling this method.
By default, this method uses WebLogo 2.8.2; this default will be changed
to WebLogo 3 in a future release of Biopython.
Note that the arguments to WebLogo 3 can be different from those to
WebLogo 2.8.2. See the WebLogo website (http://weblogo.threeplusone.com)
for more information. You can also use help() on this method to see the
arguments to WebLogo 2.8.2 and WebLogo 3.""",
                PendingDeprecationWarning)
            al= self._to_fasta()
            url = 'http://weblogo.berkeley.edu/logo.cgi'
            values = {'sequence' : al,
                      'format' : format,
                      'logowidth' : '18',
                      'logoheight' : '5',
                      'logounits' : 'cm',
                      'kind' : 'AUTO',
                      'firstnum' : "1",
                      'command' : 'Create Logo',
                      'smallsamplecorrection' : "on",
                      'symbolsperline' : 32,
                      'res' : '96',
                      'res_units' : 'ppi',
                      'antialias' : 'on',
                      'title' : '',
                      'barbits' : '',
                      'xaxis': 'on',
                      'xaxis_label' : '',
                      'yaxis': 'on',
                      'yaxis_label' : '',
                      'showends' : 'on',
                      'shrink' : '0.5',
                      'fineprint' : 'on',
                      'ticbits' : '1',
                      'colorscheme' : 'DEFAULT',
                      'color1' : 'black',
                      'color2' : 'blue',
                      'color3' : 'red',
                      'color4' : 'black',
                      'color5' : 'purple',
                      'color6' : 'orange',
                      }
            for k,v in kwds.iteritems():
                values[k]=str(v)
        else:
            frequencies= self._to_transfac()
            url = 'http://weblogo.threeplusone.com/create.cgi'
            values = {'sequences' : frequencies,
                      'format' : format.lower(),
                      'stack_width' : 'medium',
                      'stack_per_line' : '40',
                      'alphabet' : 'alphabet_dna',
                      'ignore_lower_case' : True,
                      'unit_name' : "bits",
                      'first_index' : '1',
                      'logo_start' : '1',
                      'logo_end': str(self.length),
                      'composition' : "comp_auto",
                      'percentCG' : '',
                      'scale_width' : True,
                      'show_errorbars' : True,
                      'logo_title' : '',
                      'logo_label' : '',
                      'show_xaxis': True,
                      'xaxis_label': '',
                      'show_yaxis': True,
                      'yaxis_label': '',
                      'yaxis_scale': 'auto',
                      'yaxis_tic_interval' : '1.0',
                      'show_ends' : True,
                      'show_fineprint' : True,
                      'symbols0': '',
                      'symbols1': '',
                      'symbols2': '',
                      'symbols3': '',
                      'symbols4': '',
                      'color0': '',
                      'color1': '',
                      'color2': '',
                      'color3': '',
                      'color4': '',
                      }
            for k,v in kwds.iteritems():
                if type(values[k])==bool:
                    if not v:
                        v = ""
                values[k]=str(v)

        data = urllib.urlencode(values)
        req = urllib2.Request(url, data)
        response = urllib2.urlopen(req)
        f=open(fname,"w")
        im=response.read()

        f.write(im)
        f.close()

    def _to_transfac(self):
        """Write the representation of a motif in TRANSFAC format
        """
        from Bio.Motif import TRANSFAC
        multiple_value_keys = TRANSFAC.Motif.multiple_value_keys
        sections = (('AC', 'AS',), # Accession
                    ('ID',),       # ID
                    ('DT', 'CO'),  # Date, copyright
                    ('NA',),       # Name
                    ('DE',),       # Short factor description
                    ('TY',),       # Type
                    ('OS', 'OC'),  # Organism
                    ('HP', 'HC'),  # Superfamilies, subfamilies
                    ('BF',),       # Binding factors
                    ('P0',),       # Frequency matrix
                    ('BA',),       # Statistical basis
                    ('BS',),       # Factor binding sites
                    ('CC',),       # Comments
                    ('DR',),       # External databases
                    ('OV', 'PV',), # Versions
                   )
        lines = []
        for section in sections:
            blank = False
            for key in section:
                if key=='P0':
                    # Frequency matrix
                    length = self.length
                    if length==0:
                        continue
                    sequence = self.degenerate_consensus
                    line = "P0      A      C      G      T"
                    lines.append(line)
                    for i in range(length):
                        line = "%02.d %6.20g %6.20g %6.20g %6.20g      %s" % (
                                             i+1,
                                             self.counts['A'][i],
                                             self.counts['C'][i],
                                             self.counts['G'][i],
                                             self.counts['T'][i],
                                             sequence[i],
                                            )
                        lines.append(line)
                    blank = True
                else:
                    try:
                        value = self.get(key)
                    except AttributeError:
                        value = None
                    if value is not None:
                        if key in multiple_value_keys:
                            for v in value:
                                line = "%s  %s" % (key, v)
                                lines.append(line)
                        else:
                            line = "%s  %s" % (key, value)
                            lines.append(line)
                        blank = True
                if key=='PV':
                    # References
                    try:
                        references = self.references
                    except AttributeError:
                        pass
                    else:
                        keys = ("RN", "RX", "RA", "RT", "RL")
                        for reference in references:
                            for key in keys:
                                value = reference.get(key)
                                if value is None:
                                    continue
                                line = "%s  %s" % (key, value)
                                lines.append(line)
                                blank = True
            if blank:
                line = 'XX'
                lines.append(line)
        # Finished; glue the lines together
        line = "//"
        lines.append(line)
        text = "\n".join(lines) + "\n"
        return text

    def _to_vertical_matrix(self,letters=None):
        """Return string representation of the motif as  a matrix.
        """
        warnings.warn("""\
This function is now obsolete, and will be deprecated and removed
in a future release of Biopython.
""", PendingDeprecationWarning)
        if letters is None:
            letters=self.alphabet.letters
        self._pwm_is_current=False
        pwm=self.pwm(laplace=False)
        res=""
        for i in range(self.length):
            res+="\t".join([str(pwm[i][a]) for a in letters])
            res+="\n"
        return res

    def _to_horizontal_matrix(self,letters=None,normalized=True):
        """Return string representation of the motif as  a matrix.
        """
        warnings.warn("""\
This function is now obsolete, and will be deprecated and removed
in a future release of Biopython.
""", PendingDeprecationWarning)
        if letters is None:
            letters=self.alphabet.letters
        res=""
        if normalized:  # output PWM
            self._pwm_is_current=False
            mat=self.pwm(laplace=False)
            for a in letters:
                res+="\t".join([str(mat[i][a]) for i in range(self.length)])
                res+="\n"
        else:  # output counts
            if self.counts is None:
                self.make_counts_from_instances()
            mat=self.counts
            for a in letters:
                res+="\t".join([str(mat[a][i]) for i in range(self.length)])
                res+="\n"
        return res


    def _to_jaspar_pfm(self):
        """Returns the pfm representation of the motif
        """
        letters = "ACGT"
        counts = self.counts
        length = self.length
        lines = []
        for letter in letters:
            terms = [str(counts[letter][i]) for i in range(length)]
            line = "\t".join(terms) + "\n"
            lines.append(line)
        # Finished; glue the lines together
        text = "".join(lines)
        return text


    def format(self,format):
        """Returns a string representation of the Motif in a given format

        Currently supported fromats:
         - pfm : JASPAR Position Frequency Matrix
         - transfac : TRANSFAC like files
         - fasta : FASTA file with instances
        """

        formatters={
            "jaspar-pfm":   self._to_jaspar_pfm,
            "pfm":   self._to_jaspar_pfm,
            "transfac":     self._to_transfac,
            "fasta" :       self._to_fasta,
            }

        if format=="jaspar-pfm":
            warnings.warn("Please uses 'pfm' instead of 'jaspar-pfm'.",
PendingDeprecationWarning)
        try:
            return formatters[format]()
        except KeyError:
            raise ValueError("Wrong format type")


    def scanPWM(self,seq):
        """Matrix of log-odds scores for a nucleotide sequence.

        scans a nucleotide sequence and returns the matrix of log-odds
        scores for all positions.

        - the result is a one-dimensional list or numpy array
        - the sequence can only be a DNA sequence
        - the search is performed only on one strand
        """
        warnings.warn("""\
This function is now obsolete, and will be deprecated and removed
in a future release of Biopython. As a replacement, instead of
>>> motif.scanPWM(sequence)
use
>>> pwm = motif.counts.normalize()
>>> pssm = pwm.log_odds()
>>> pssm.calculate(sequence)
See the documentation of motif.counts.normalize, pwm.log_odds, and
pssm.calculate for details.
""", PendingDeprecationWarning)
        if self.alphabet!=IUPAC.unambiguous_dna:
            raise ValueError("Wrong alphabet! Use only with DNA motifs")
        if seq.alphabet!=IUPAC.unambiguous_dna:
            raise ValueError("Wrong alphabet! Use only with DNA sequences")

        seq = str(seq)

        # check if the fast C code can be used
        try:
            import _pwm
        except ImportError:
            # use the slower Python code otherwise
            return self._pwm_calculate(seq)

        # get the log-odds matrix into a proper shape
        # (each row contains sorted (ACGT) log-odds values)
        logodds=[[y[1] for y in sorted(x.items())] for x in self.log_odds()]
        return _pwm.calculate(seq, logodds)


    def _pwm_calculate(self, sequence):
        warnings.warn("""\
This function is now obsolete, and will be deprecated and removed
in a future release of Biopython. As a replacement, instead of
>>> motif._pwm_calculate(sequence)
use
>>> pwm = motif.counts.normalize()
>>> pssm = pwm.log_odds()
>>> pssm.calculate(sequence)
See the documentation of motif.counts.normalize, pwm.log_odds, and
pssm.calculate for details.
""", PendingDeprecationWarning)
        logodds = self.log_odds()
        m = len(logodds)
        s = len(sequence)
        n = s - m + 1
        result = [None] * n
        for i in xrange(n):
            score = 0.0
            for j in xrange(m):
                c = sequence[i+j]
                temp = logodds[j].get(c)
                if temp is None:
                    break
                score += temp
            else:
                result[i] = score
        return result
