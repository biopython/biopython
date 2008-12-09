#Copyright 2008 by Bartek Wilczynski, adapted from Bio.MEME.Motif by Jason A. Hackney
# Copyright 2004 by Jason A. Hackney.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio import Seq
from Bio.Alphabet import IUPAC
from math import sqrt
import sys
from Motif import Motif

class MEMEMotif (Motif):
    """A subclass of Motif used in parsing MEME (and MAST) output.
    
    This sublcass defines functions and data specific to MEME motifs. 
    This includes the evalue for a motif and the PSSM of the motif.
    
    Methods:
    add_instance_from_values (name = 'default', pvalue = 1, sequence = 'ATA', start = 0, strand = +): create a new instance of the motif with the specified values.
    add_to_pssm (position): add a new position to the pssm. The position should be a list of nucleotide/amino acid frequencies
    add_to_logodds (position): add a new position to the log odds matrix. The position should be a tuple of log odds values for the nucleotide/amino acid at that position.
    compare_motifs (other_motif): returns the maximum correlation between this motif and other_motif
    """
    def __init__ (self):
        Motif.__init__(self)
        self.evalue = 0.0
    
    def _numoccurrences (self, number):
        if type(number) == int:
            self.num_occurrences = number
        else:
            number = int(number)
            self.num_occurrences = number

    def get_instance_by_name (self,name):
        for i in self.instances:
            if i.sequence_name == name:
                return i
        return None

    def add_instance_from_values (self, name = 'default', pvalue = 1, sequence = 'ATA', start = 0, strand = '+'):
        inst = MEMEInstance(sequence,self.alphabet)
        inst._pvalue(pvalue)
        inst._seqname(name)
        inst._start(start)
        inst._strand(strand)
        if self.length:
            inst._length(self.length)
        if self.name:
            inst._motifname(self.name)
        self.add_instance(inst)
    
    def _evalue (self, evalue):
        if type(evalue) == float:
            self.evalue = evalue
        else:
            evalue = float(evalue)
            self.evalue = evalue
    

class MEMEInstance(Seq.Seq):
    """A class describing the instances of a MEME motif, and the data thereof. 
    """
    def __init__ (self,*args,**kwds):
        Seq.Seq.__init__(self,*args,**kwds)
        self.sequence_name = ""
        self.start = 0
        self.pvalue = 1.0
        self.strand = 0
        self.length = 0
        self.motif_name = ""
        
    
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
    
