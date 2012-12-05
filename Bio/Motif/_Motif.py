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

class Motif(object):
    """
    A class representing sequence motifs.
    """
    def __init__(self, alphabet=None, instances=None, counts=None):
        self.mask = []
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
            self.counts = counts
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
            self.counts = {}
            for letter in alphabet.letters:
                self.counts[letter] = [0] * self.length
            for instance in self.instances:
                for position, letter in enumerate(instance):
                    self.counts[letter][position] += 1
        else:
            self.counts = None
            self.instances = None
            if alphabet is None:
                alphabet = IUPAC.unambiguous_dna
        self.alphabet = alphabet
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
        if self.length is None:
            self.length = len
        elif self.length != len:
            print "len",self.length,self.instances, len
            raise ValueError("You can't change the length of the motif")

    def _check_alphabet(self,alphabet):
        if self.alphabet is None:
            self.alphabet=alphabet
        elif self.alphabet != alphabet:
                raise ValueError("Wrong Alphabet")

    def add_instance(self,instance):
        """
        adds new instance to the motif
        """
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
        self._check_length(len(mask))
        self.mask=[]
        for char in mask:
            if char=="*":
                self.mask.append(1)
            elif char==" ":
                self.mask.append(0)
            else:
                raise ValueError("Mask should contain only '*' or ' ' and not a '%s'"%char)

    def pwm(self,laplace=True):
        """
        returns the PWM computed for the set of instances

        if laplace=True (default), pseudocounts equal to self.background multiplied by self.beta are added to all positions.
        """

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
                    except KeyError: #we need to ignore non-alphabet letters
                        pass
            self._pwm.append(FreqTable.FreqTable(dict,FreqTable.COUNT,self.alphabet))
        self._pwm_is_current=1
        return self._pwm

    def log_odds(self,laplace=True):
        """
        returns the logg odds matrix computed for the set of instances
        """
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

    def ic(self, background=None):
        """\
Returns the information content of a motif.

By default, a uniform background is used. To specify a non-uniform
background, use the 'background' argument to pass a dictionary containing
the probability of each letter in the alphabet associated with the motif
under the background distribution.
        """
        result=0
        if background==None:
            background = {}
            for a in self.alphabet.letters:
                background[a] = 1.0
        total = sum(background.values())
        for a in self.alphabet.letters:
            background[a] /= total
        for a in self.alphabet.letters:
            if background[a]!=0:
                result-=background[a]*math.log(background[a],2)
        result *= self.length
        pwm=self.pwm()
        for i in range(self.length):
            for a in self.alphabet.letters:
                if pwm[i][a]!=0:
                    result+=pwm[i][a]*math.log(pwm[i][a],2)
        return result

    def exp_score(self,st_dev=False):
        """
        Computes expected score of motif's instance and its standard deviation
        """
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
        lo=self.log_odds()
        score = 0.0
        for pos in xrange(self.length):
            a = sequence[position+pos]
            if not masked or self.mask[pos]:
                try:
                    score += lo[pos][a]
                except:
                    pass
        if normalized:
            if not masked:
                score/=self.length
            else:
                score/=len([x for x in self.mask if x])
        return score

    def search_pwm(self,sequence,normalized=0,masked=0,threshold=0.0,both=True):
        """
        a generator function, returning found hits in a given sequence with the pwm score higher than the threshold
        """
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

        if self.alphabet != motif.alphabet:
            raise ValueError("Cannot compare motifs with different alphabets")

        max_p=-2
        for offset in range(-self.length+1,motif.length):
            if offset<0:
                p = self.dist_pearson_at(motif,-offset)
            else: #offset>=0
                p = motif.dist_pearson_at(self,offset)

            if max_p<p:
                max_p=p
                max_o=-offset
        return 1-max_p,max_o

    def dist_pearson_at(self,motif,offset):
        sxx = 0 # \sum x^2
        sxy = 0 # \sum x \cdot y
        sx = 0  # \sum x
        sy = 0  # \sum y
        syy = 0 # \sum x^2
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
        A similarity measure taking into account a product probability of generating overlaping instances of two motifs
        """
        max_p=0.0
        for offset in range(-self.length+1,other.length):
            if offset<0:
                p = self.dist_product_at(other,-offset)
            else: #offset>=0
                p = other.dist_product_at(self,offset)
            if max_p<p:
                max_p=p
                max_o=-offset
        return 1-max_p/self.dist_product_at(self,0),max_o

    def dist_product_at(self,other,offset):
        s=0
        for i in range(max(self.length,offset+other.length)):
            f1=self[i]
            f2=other[i-offset]
            for n,b in self.background.iteritems():
                s+=b*f1[n]*f2[n]
        return s/i

    def dist_dpq(self,other):
        r"""Calculates the DPQ distance measure between motifs.

        It is calculated as a maximal value of DPQ formula (shown using LaTeX
        markup, familiar to mathematicians):

        \sqrt{\sum_{i=1}^{alignment.len()} \sum_{k=1}^alphabet.len() \
        \{ m1[i].freq(alphabet[k])*log_2(m1[i].freq(alphabet[k])/m2[i].freq(alphabet[k])) +
           m2[i].freq(alphabet[k])*log_2(m2[i].freq(alphabet[k])/m1[i].freq(alphabet[k]))
        }

        over possible non-spaced alignemts of two motifs.  See this reference:

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
            else: #offset>=0
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
        return min_d,min_o#,d_s

    def dist_dpq_at(self,other,offset):
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
                if self.mask[i]:
                    string += "*"
                else:
                    string += " "
            string += "\n"
        return string

    def __len__(self):
        """return the length of a motif

        Please use this method (i.e. invoke len(m)) instead of refering to the m.length directly.
        """
        if self.length is None:
            return 0
        else:
            return self.length

    def _write(self,stream):
        """
        writes the motif to the stream
        """

        stream.write(self.__str__())

    def _to_fasta(self):
        """
        FASTA representation of motif
        """
        warnings.warn("""\
This function is now obsolete, and will be deprecated and removed in
a future release of Biopython.""", PendingDeprecationWarning)
        if self.instances is None:
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
                instance="" #start with empty seq
                for j in range(self.length): #iterate over positions
                    instance+=col[j][i]
                instance = Seq(instance, self.alphabet)
                instances.append(instance)
        else:
            instances = self.instances
        string = ""
        for i, instance in enumerate(instances):
            string += ">instance%d\n%s\n "% (i, instance)
        return string

    def reverse_complement(self):
        """
        Gives the reverse complement of the motif
        """
        alphabet = self.alphabet
        if self.instances is not None:
            instances = []
            for instance in self.instances:
                instance = instance.reverse_complement()
                instances.append(instance)
            res = Motif(alphabet, instances)
        else: # has counts
            res = Motif(alphabet)
            res.counts={}
            res.counts["A"]=self.counts["T"][:]
            res.counts["T"]=self.counts["A"][:]
            res.counts["G"]=self.counts["C"][:]
            res.counts["C"]=self.counts["G"][:]
            res.counts["A"].reverse()
            res.counts["C"].reverse()
            res.counts["G"].reverse()
            res.counts["T"].reverse()
            res.length=self.length
        res.mask = self.mask
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
        self.set_mask("*"*self.length)
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
            except ValueError: #not integers
                self.counts[i]=map(float,ln) #map(lambda s: int(100*float(s)),ln)
            #print counts[i]

        s = sum(self.counts[nuc][0] for nuc in letters)
        l = len(self.counts[letters[0]])
        self.length=l
        self.set_mask("*"*l)
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
            inst="" #start with empty seq
            for j in range(self.length): #iterate over positions
                inst+=col[j][i]
            #print i,inst
            inst=Seq(inst,self.alphabet)
            self.add_instance(inst)
        return self.instances

    def make_counts_from_instances(self):
        """Creates the count matrix for a motif with instances.

        """
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
        self.set_mask("*"*len(instance))
        return self

    def __getitem__(self,index):
        """Returns the probability distribution over symbols at a given position, padding with background.

        If the requested index is out of bounds, the returned distribution comes from background.
        """
        if index in range(self.length):
            return self.pwm()[index]
        else:
            return self.background

    @property
    def consensus(self):
        """Returns the consensus sequence of a motif.
        """
        res=""
        for i in range(self.length):
            max_f=0
            max_n="X"
            for n in sorted(self[i]):
                if self[i][n]>max_f:
                    max_f=self[i][n]
                    max_n=n
            res+=max_n
        sequence = ConsensusSeq(res,self.alphabet)
        # This is an ugly hack that allows us raise a warning if a user
        # attempts to call the method motif.consensus() instead of accessing
        # the property motif.consensus.
        return sequence

    @property
    def anticonsensus(self):
        """returns the least probable pattern to be generated from this motif.
        """
        res=""
        for i in range(self.length):
            min_f=10.0
            min_n="X"
            for n in sorted(self[i]):
                if self[i][n]<min_f:
                    min_f=self[i][n]
                    min_n=n
            res+=min_n
        sequence = ConsensusSeq(res,self.alphabet)
        # This is an ugly hack that allows us raise a warning if a user
        # attempts to call the method motif.anticonsensus() instead of
        # accessing the property motif.anticonsensus.
        return sequence

    @property
    def degenerate_consensus(self):
        """Following the rules adapted from
D. R. Cavener: "Comparison of the consensus sequence flanking
translational start sites in Drosophila and vertebrates."
Nucleic Acids Research 15(4): 1353-1361. (1987).
The same rules are used by TRANSFAC."""
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
        res = ""
        for i in range(self.length):
            def get(nucleotide):
                return self.counts[nucleotide][i]
            nucleotides = sorted(self.counts, key=get, reverse=True)
            counts = [self.counts[c][i] for c in nucleotides]
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
            res += nucleotide
        sequence = Seq(res, alphabet = IUPAC.ambiguous_dna)
        return sequence

    def max_score(self):
        """Maximal possible score for this motif.

        returns the score computed for the consensus sequence.
        """
        return self.score_hit(self.consensus,0)

    def min_score(self):
        """Minimal possible score for this motif.

        returns the score computed for the anticonsensus sequence.
        """
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
        res="XX\nTY Motif\n" #header
        try:
            res+="ID %s\n"%self.name
        except:
            pass
        res+="BF undef\nP0"
        for a in self.alphabet.letters:
            res+=" %s"%a
        res+="\n"
        if self.counts is None:
            self.make_counts_from_instances()
        for i in range(self.length):
            if i<9:
                res+="0%d"%(i+1)
            else:
                res+="%d"%(i+1)
            for a in self.alphabet.letters:
                res+=" %d"%self.counts[a][i]
            res+="\n"
        res+="XX\n"
        return res

    def _to_vertical_matrix(self,letters=None):
        """Return string representation of the motif as  a matrix.
        """
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
        if letters is None:
            letters=self.alphabet.letters
        res=""
        if normalized: #output PWM
            self._pwm_is_current=False
            mat=self.pwm(laplace=False)
            for a in letters:
                res+="\t".join([str(mat[i][a]) for i in range(self.length)])
                res+="\n"
        else: #output counts
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
        return self._to_horizontal_matrix(normalized=False,letters="ACGT")

    def format(self,format):
        """Returns a string representation of the Motif in a given format

        Currently supported fromats:
         - jaspar-pfm : JASPAR Position Frequency Matrix
         - transfac : TRANSFAC like files
         - fasta : FASTA file with instances
        """

        formatters={
            "jaspar-pfm":   self._to_jaspar_pfm,
            "transfac":     self._to_transfac,
            "fasta" :       self._to_fasta,
            }

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
