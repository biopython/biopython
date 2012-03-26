# Copyright 2008 by Bartek Wilczynski
# Adapted from  Bio.MEME.Parser by Jason A. Hackney.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio.Alphabet import IUPAC
from Bio import Seq
import re
from math import sqrt
import sys
from Bio.Motif import Motif



def read(handle):
    """Parses the text output of the MEME program into MEME.Record object.
    
    Example:
    
    >>> f = open("meme.output.txt")
    >>> from Bio.Motif.Parsers import MEME
    >>> record = MEME.read(f)
    >>> for motif in record.motifs:
    ...     for instance in motif.instances:
    ...         print instance.motif_name, instance.sequence_name, instance.strand, instance.pvalue
    
    """
    record = MEMERecord()
    __read_version(record, handle)
    __read_datafile(record, handle)
    __read_alphabet(record, handle)
    __read_sequence_names(record, handle)
    __read_command(record, handle)
    for line in handle:
        if line.startswith('MOTIF  1'):
            break
    else:
        raise ValueError('Unexpected end of stream')
    while True:
        motif = __create_motif(line)
        motif.alphabet = record.alphabet
        record.motifs.append(motif)
        __read_motif_name(motif, handle)
        __read_motif_sequences(motif, handle, 'revcomp' in record.command)
        __skip_unused_lines(handle)
        try:
            line = handle.next()
        except StopIteration:
            raise ValueError('Unexpected end of stream: Expected to find new motif, or the summary of motifs')
        if line.startswith("SUMMARY OF MOTIFS"):
            break
        if not line.startswith('MOTIF'):
            raise ValueError("Line does not start with 'MOTIF':\n%s" % line)
    return record


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
        else:
            inst._length(len(sequence))
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
    

class MEMERecord(object):
    """A class for holding the results of a MEME run.
    
    A MEMERecord is an object that holds the results from running
    MEME. It implements no methods of its own.
        
    """
    def __init__ (self):
        """__init__ (self)"""
        self.motifs = []
        self.version = ""
        self.datafile = ""
        self.command = ""
        self.alphabet = None
        self.sequence_names = []
        
    def get_motif_by_name (self, name):
        for m in self.motifs:
            if m.name == name:
                return m


# Everything below is private


def __read_version(record, handle):
    for line in handle:
        if line.startswith('MEME version'):
            break
    else:
        raise ValueError("Improper input file. File should contain a line starting MEME version.")
    line = line.strip()
    ls = line.split()
    record.version = ls[2]


def __read_datafile(record, handle):
    for line in handle:
        if line.startswith('TRAINING SET'):
            break
    else:
        raise ValueError("Unexpected end of stream: 'TRAINING SET' not found.")
    try:
        line = handle.next()
    except StopIteration:
        raise ValueError("Unexpected end of stream: Expected to find line starting with '****'")
    if not line.startswith('****'):
        raise ValueError("Line does not start with '****':\n%s" % line)
    try:
        line = handle.next()
    except StopIteration:
        raise ValueError("Unexpected end of stream: Expected to find line starting with 'DATAFILE'")
    if not line.startswith('DATAFILE'):
        raise ValueError("Line does not start with 'DATAFILE':\n%s" % line)
    line = line.strip()
    line = line.replace('DATAFILE= ','')
    record.datafile = line


def __read_alphabet(record, handle):
    try:
        line = handle.next()
    except StopIteration:
        raise ValueError("Unexpected end of stream: Expected to find line starting with 'ALPHABET'")
    if not line.startswith('ALPHABET'):
        raise ValueError("Line does not start with 'ALPHABET':\n%s" % line)
    line = line.strip()
    line = line.replace('ALPHABET= ','')
    if line == 'ACGT':
        al = IUPAC.unambiguous_dna
    else:
        al = IUPAC.protein
    record.alphabet = al


def __read_sequence_names(record, handle):
    try:
        line = handle.next()
    except StopIteration:
        raise ValueError("Unexpected end of stream: Expected to find line starting with 'Sequence name'")
    if not line.startswith('Sequence name'):
        raise ValueError("Line does not start with 'Sequence name':\n%s" % line)
    try:
        line = handle.next()
    except StopIteration:
        raise ValueError("Unexpected end of stream: Expected to find line starting with '----'")
    if not line.startswith('----'):
        raise ValueError("Line does not start with '----':\n%s" % line)
    for line in handle:
        if line.startswith('***'):
            break
        line = line.strip()
        ls = line.split()
        record.sequence_names.append(ls[0])
        if len(ls) == 6:
            record.sequence_names.append(ls[3])
    else:
        raise ValueError("Unexpected end of stream: Expected to find line starting with '***'")


def __read_command(record, handle):
    for line in handle:
        if line.startswith('command:'):
            break
    else:
        raise ValueError("Unexpected end of stream: Expected to find line starting with 'command'")
    line = line.strip()
    line = line.replace('command: ','')
    record.command = line


def __create_motif(line):
    line = line[5:].strip()
    ls = line.split()
    motif = MEMEMotif()
    motif.length = int(ls[3])
    motif._numoccurrences(ls[6])
    motif._evalue(ls[12])
    return motif


def __read_motif_name(motif, handle):
    for line in handle:
        if 'sorted by position p-value' in line:
            break
    else:
        raise ValueError('Unexpected end of stream: Failed to find motif name')
    line = line.strip()
    ls = line.split()
    name = " ".join(ls[0:2])
    motif.name=name


def __read_motif_sequences(motif, handle, rv):
    try:
        line = handle.next()
    except StopIteration:
        raise ValueError('Unexpected end of stream: Failed to find motif sequences')
    if not line.startswith('---'):
        raise ValueError("Line does not start with '---':\n%s" % line)
    try:
        line = handle.next()
    except StopIteration:
        raise ValueError("Unexpected end of stream: Expected to find line starting with 'Sequence name'")
    if not line.startswith('Sequence name'):
        raise ValueError("Line does not start with 'Sequence name':\n%s" % line)
    try:
        line = handle.next()
    except StopIteration:
        raise ValueError('Unexpected end of stream: Failed to find motif sequences')
    if not line.startswith('---'):
        raise ValueError("Line does not start with '---':\n%s" % line)
    for line in handle:
        if line.startswith('---'):
            break
        line = line.strip()
        ls = line.split()
        if rv:
            #seq = Seq.Seq(ls[5], record.alphabet)
            motif.add_instance_from_values(name = ls[0], sequence = ls[5], start = ls[2], pvalue = ls[3], strand = ls[1])
        else:
            #seq = Seq.Seq(ls[4], record.alphabet)
            motif.add_instance_from_values(name = ls[0], sequence = ls[4], start = ls[1], pvalue = ls[2])
    else:
        raise ValueError('Unexpected end of stream')


def __skip_unused_lines(handle):
    for line in handle:
        if line.startswith('log-odds matrix'):
            break
    else:
        raise ValueError("Unexpected end of stream: Expected to find line starting with 'log-odds matrix'")
    for line in handle:
        if line.startswith('---'):
            break
    else:
        raise ValueError("Unexpected end of stream: Expected to find line starting with '---'")
    for line in handle:
        if line.startswith('letter-probability matrix'):
            break
    else:
        raise ValueError("Unexpected end of stream: Expected to find line starting with 'letter-probability matrix'")
    for line in handle:
        if line.startswith('---'):
            break
    else:
        raise ValueError("Unexpected end of stream: Expected to find line starting with '---'")
    for line in handle:
        if line.startswith('Time'):
            break
    else:
        raise ValueError("Unexpected end of stream: Expected to find line starting with 'Time'")
    try:
        line = handle.next()
    except StopIteration:
        raise ValueError('Unexpected end of stream: Expected to find blank line')
    if line.strip():
        raise ValueError("Expected blank line, but got:\n%s" % line)
    try:
        line = handle.next()
    except StopIteration:
        raise ValueError("Unexpected end of stream: Expected to find line starting with '***'")
    if not line.startswith('***'):
        raise ValueError("Line does not start with '***':\n%s" % line)
    for line in handle:
        if line.strip():
            break
    else:
        raise ValueError("Unexpected end of stream: Expected to find line starting with '***'")
    if not line.startswith('***'):
        raise ValueError("Line does not start with '***':\n%s" % line)

