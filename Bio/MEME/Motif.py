# Copyright 2004 by Jason A. Hackney.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio import Seq
from Bio.Alphabet import IUPAC
from math import sqrt

class Motif:
    """A generic motif class.
    
    A Motif has information about the alphabet used in the motif,
    the length of the motif, the motif name, the number of occurrences,
    and the individual instances of the motif. Each of the instances
    is an object of the Instance class.
    
    Methods:
    add_instance(instance): adds a new instance of the motif.
    get_instance_by_name(name): returns an instance with the specified name.
    """
    def __init__(self):
        self.instances = []
        self.alphabet = None
        self.length = 0
        self.num_occurrences = 0
        self.name = ""
        self.consensus = ""
    
    def add_instance (self, instance):
        if isinstance(instance, Instance):
            if self.length:
                instance._length(self.length)
            if self.name:
                instance._motifname(self.name)
            self.instances.append(instance)
    
    def get_instance_by_name (self,name):
        for i in set.instances:
            if i.seqname == name:
                return i
        return None
    
    def _alphabet (self, alphabet):
        if alphabet == IUPAC.unambiguous_dna or alphabet == IUPAC.protein or alphabet == IUPAC.ambiguous_dna:
            self.alphabet = alphabet
        else:
            return -1
    
    def _length (self, length):
        if type(length) == int:
            self.length = length
        else:
            length = int(length)
            self.length = length
    
    def _name (self, name):
        self.name = name    
    
    def _consensus (self, consensus):
        if self.alphabet:
            self.consensus = Seq.Seq(consensus, self.alphabet)
        else:
            self.consensus = consensus
    
    def _numoccurrences (self, number):
        if type(number) == int:
            self.num_occurrences = number
        else:
            number = int(number)
            self.num_occurrences = number
    


class MEMEMotif (Motif):
    """A subclass of Motif used in parsing MEME (and MAST) output.
    
    This sublcass defines functions and data specific to MEME motifs. 
    This includes the evalue for a motif and the PSSM of the motif.
    
    Methods:
    add_instance_from_values (name = 'default', pvalue = 1, sequence = 'ATA', start = 0, strand = +): create a new instance of the motif with the specified values.
    add_to_pssm (position): add a new position to the pssm. The position should be a list of nucleotide/amino acid frequencies
    compare_motifs (other_motif): returns the maximum correlation between this motif and other_motif
    """
    def __init__ (self):
        Motif.__init__(self)
        self.evalue = 0.0
        self.pssm = []
    
    def add_instance_from_values (self, name = 'default', pvalue = 1, sequence = 'ATA', start = 0, strand = '+'):
        inst = Instance()
        inst._sequence(sequence)
        inst._pvalue(pvalue)
        inst._seqname(name)
        inst._start(start)
        inst._strand(strand)
        if self.length:
            inst._length(self.length)
        if self.name:
            inst._motifname(self.name)
        self.instances.append(inst)
    
    def _evalue (self, evalue):
        if type(evalue) == float:
            self.evalue = evalue
        else:
            evalue = float(evalue)
            self.evalue = evalue
    
    def add_to_pssm (self, thisposition):
        self.pssm.append(thisposition)
    
    def compare_motifs (self, motif):
        if isinstance(motif, MEMEMotif):
            mylen = len(self.pssm)
            yourlen = len(motif.pssm)
            myr = None
            if mylen == yourlen:
                myr = _corr(self.pssm, motif.pssm)
            elif mylen < yourlen:
                diff = yourlen - mylen
                for i in range(0,diff):
                    yourpssm = motif.pssm[i:mylen+i]
                    r = _corr(self.pssm, yourpssm)
                    if r > myr:
                        myr = r
            else:
                diff = mylen - yourlen
                for i in range(0, diff):
                    mypssm = self.pssm[i:yourlen+i]
                    r = _corr(mypssm, motif.pssm)
                    if r > myr:
                        myr = r
            return myr
        else:
            sys.stderr.write(str(m2))
            return None
    


class Instance:
    """A class describing the instances of a motif, and the data thereof. 
    """
    def __init__ (self):
        self.sequence = None
        self.sequence_name = ""
        self.start = 0
        self.pvalue = 1.0
        self.strand = 0
        self.length = 0
        self.motif_name = ""
        
    def _sequence (self, seq):
        self.sequence = seq
    
    def _seqname (self, name):
        self.sequence_name = name
        
    def _motifname (self, name):
        self.motif_name = name
    
    def _start (self,start):
        start = int(start)
        self.start = start
    
    def _pvalue (self,pval):
        pval = float(pval)
        self.pvalue = pval
    
    def _score (self, score):
        score = float(score)
        self.score = score
    
    def _strand (self, strand):
        self.strand = strand
    
    def _length (self, length):
        self.length = length
    


def _corr (x,y):
    sx = 0
    sy = 0
    sxx = 0
    syy = 0
    sxy = 0
    #make the two-dimensional lists one dimensional
    if type(x[0] == list):
        x = reduce(lambda a,b: a+b, x)
    if type(y[0] == list):
        y = reduce(lambda a,b: a+b, y)
    sq = lambda t: t*t
    length = len(x)
    for a,b in zip(x,y):
        sx += a
        sy += b
        sxx += sq(a)
        syy += sq(b)
        sxy += a*b
    s1 = sxy - sx * sy * 1.0/length
    s2 = (sxx - sx * sy * 1.0/length) * (syy - sx * sy * 1.0/length)
    r = s1 / sqrt (s2)
    return r
    
