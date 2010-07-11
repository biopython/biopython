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
    

class MEMERecord:
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
    line = line.strip()
    ls = line.split()
    motif = MEMEMotif()
    motif.length = int(ls[4])
    motif._numoccurrences(ls[7])
    motif._evalue(ls[13])
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


# Everything below is obsolete.


from Bio import File
from Bio.ParserSupport import *


class MEMEParser (AbstractParser):
    """A parser for the text output of the MEME program (OBSOLETE).
    Parses the output into an object of the MEMERecord class.
    
    Methods:
    parse (handle): parses the contents of the file handle passed to it.
    
    Example:
    
    >>>f = open("meme.output.txt")
    >>>parser = MEMEParser()
    >>>meme_record = parser.parse(f)
    >>>for motif in meme_record.motifs:
    ...    for instance in motif.instances:
    ...        print instance.motif_name, instance.sequence_name, instance.strand, instance.pvalue
    
    This class is OBSOLETE; please use the read() function in this module
    instead.
    """
    def __init__ (self):
        """__init__ (self)"""
        self._scanner = _MEMEScanner()
        self._consumer = _MEMEConsumer()
    
    def parse (self, handle):
        """parse (self, handle)"""
        self._scanner.feed(handle, self._consumer)
        return self._consumer.data
    


class _MEMEScanner:
    """Scanner for MEME output (OBSOLETE).
    
    Methods:
    feed
    
    This class is OBSOLETE; please use the read() function in this module
    instead.
    """
    
    def feed (self, handle, consumer):
        """
        Feeds in MEME output for scanning. handle should
        implement the readline method. consumer is 
        a Consumer object that can receive the salient events.
        """
        if isinstance(handle, File.UndoHandle):
            uhandle = handle
        else:
            uhandle = File.UndoHandle(handle)
        
        self._scan_header(uhandle, consumer)
        self._scan_motifs    (uhandle, consumer)
    
    def _scan_header(self, uhandle, consumer):
        try:
            read_and_call_until(uhandle, consumer.noevent, contains = 'MEME version')
        except ValueError:
            raise ValueError("Improper input file. File should contain a line starting MEME version.")
        read_and_call(uhandle, consumer._version, start = 'MEME version')
        read_and_call_until(uhandle, consumer.noevent, start = 'TRAINING SET')
        read_and_call(uhandle, consumer.noevent, start = 'TRAINING SET')
        read_and_call(uhandle, consumer.noevent, start = '****')
        read_and_call(uhandle, consumer._datafile, start = 'DATAFILE')
        read_and_call(uhandle, consumer._alphabet, start = 'ALPHABET')
        read_and_call(uhandle, consumer.noevent, start = 'Sequence name')
        read_and_call(uhandle, consumer.noevent, start = '----')
        read_and_call_until(uhandle, consumer._sequence_name, start = '***')
        read_and_call_until(uhandle, consumer.noevent, start = 'command:')
        read_and_call(uhandle, consumer._commandline, start = 'command:')
        read_and_call_until(uhandle, consumer.noevent, start = 'MOTIF  1')
    
    def _scan_motifs(self, uhandle, consumer):
        while 1:
            read_and_call(uhandle, consumer._add_motif_with_info, start = 'MOTIF')
            read_and_call_until(uhandle, consumer.noevent, contains = 'sorted by position p-value')
            read_and_call(uhandle, consumer.motif_name, contains = 'sorted by position p-value')
            read_and_call(uhandle, consumer.noevent, start = '---')
            read_and_call(uhandle, consumer.noevent, start = 'Sequence name')
            read_and_call(uhandle, consumer.noevent, start = '---')
            read_and_call_until(uhandle, consumer.add_instance, start = '---')
            read_and_call_until(uhandle, consumer.noevent, start = 'log-odds matrix')
            read_and_call(uhandle, consumer.noevent)
            read_and_call_until(uhandle, consumer.add_to_logodds, start = '---')
            read_and_call_until(uhandle, consumer.noevent, start = 'letter-probability matrix')
            read_and_call(uhandle, consumer.noevent, start = 'letter-probability matrix')    
            read_and_call_until(uhandle, consumer.add_to_pssm, start = '---')
            read_and_call_until(uhandle, consumer.noevent, start = 'Time')
            read_and_call(uhandle, consumer.noevent, start = 'Time')
            read_and_call(uhandle, consumer.noevent, blank = 1)
            read_and_call(uhandle, consumer.noevent, start = '***')
            read_and_call_while(uhandle, consumer.noevent, blank = 1)
            read_and_call(uhandle, consumer.noevent, start = '***')
            line = safe_peekline(uhandle)
            if line.startswith("SUMMARY OF MOTIFS"):
                break
    


class _MEMEConsumer:
    """
    Consumer that can receive events from MEME Scanner (OBSOLETE).
    
    This is the Consumer object that should be passed to the 
    MEME Scanner.
    
    This class is OBSOLETE; please use the read() function in this module
    instead.
    """
    
    def __init__ (self):
        self.current_motif = None
        self.sequence_names = []
        self.data = MEMERecord()
    
    def _version (self, line):
        line = line.strip()
        ls = line.split()
        self.data.version = ls[2]
    
    def _datafile (self, line):
        line = line.strip()
        line = line.replace('DATAFILE= ','')
        self.data.datafile = line
    
    def _alphabet (self, line):
        line = line.strip()
        line = line.replace('ALPHABET= ','')
        if line == 'ACGT':
            al = IUPAC.unambiguous_dna
        else:
            al = IUPAC.protein
        self.data.alphabet = al
    
    def _sequence_name (self, line):
        line = line.strip()
        ls = line.split()
        self.data.sequence_names.append(ls[0])
        if len(ls) == 6:
            self.data.sequence_names.append(ls[3])
    
    def _commandline (self, line):
        line = line.strip()
        line = line.replace('command: ','')
        self.data.command = line
    
    def _add_motif_with_info (self, line):
        line = line.strip()
        ls = line.split()
        motif = MEMEMotif()
        #motif.length=ls[4]
        motif._numoccurrences(ls[7])
        motif._evalue(ls[13])
        motif.alphabet=self.data.alphabet
        self.data.motifs.append(motif)
        self.current_motif = motif
    
    def motif_name (self, line):
        line = line.strip()
        ls = line.split()
        name = " ".join(ls[0:2])
        self.current_motif.name=name
    
    def add_instance (self, line):
        line = line.strip()
        ls = line.split()
        if self.data.command.find('revcomp') != -1:
            #seq = Seq.Seq(ls[5], self.data.alphabet)
            self.current_motif.add_instance_from_values(name = ls[0], sequence = ls[5], start = ls[2], pvalue = ls[3], strand = ls[1])
        else:
            #seq = Seq.Seq(ls[4], self.data.alphabet)
            self.current_motif.add_instance_from_values(name = ls[0], sequence = ls[4], start = ls[1], pvalue = ls[2])
    
    def add_to_pssm (self, line):
        pass
    
    def add_to_logodds (self, line):
        pass
    
    def noevent (self,line):
        pass
    


class _MASTConsumer:
    """
    Consumer that can receive events from _MASTScanner (OBSOLETE).
    
    A _MASTConsumer parses lines from a mast text output file.
    The motif match diagrams are parsed using line buffering. 
    Each of the buffering functions have a dummy variable that is
    required for testing using the Bio.ParserSupport.TaggingConsumer.
    If this variable isn't there, the TaggingConsumer barfs. In
    the _MASTScanner, None is passed in the place of this variable.
    
    This class is OBSOLETE; please use the read() function in the module
    Bio.Motif.Parsers.MAST instead.
    """
    def __init__ (self):
        self.data = MASTRecord()
        self._current_seq = ""
        self._line_buffer = []
        self._buffer_size = 0
        self._buffered_seq_start = 0
    
    def _version (self, line):
        line = line.strip()
        ls = line.split()
        self.data._version(ls[2])
    
    def _database (self, line):
        line = line.strip()
        ls = line.split()
        self.data._database(ls[1])
        al = ""
        if ls[2] == '(nucleotide)':
            al = IUPAC.unambiguous_dna
            self.data._alphabet(al)        
        else:
            al = IUPAC.protein
            self.data._alphabet(al)
        
    def _add_motif (self, line):
        line = line.strip()
        ls = line.split()
        m = MEMEMotif()
        m.alphabet=self.data.alphabet
        m.length=ls[1]
        name = ls[0]
        m.name=name
        m.add_instance(ls[2])
        self.data._add_motif(m)
    
    def _add_match_diagram (self, line):
        line = line.strip()
        ls = line.split()
        self.data._add_diagram_for_sequence(ls[1], self._current_seq)
        ds = ls[1].split('_')
        i = 0
        start = 0
        for i in range(0,len(ds)):
            if ds[i].find('[') != -1 or ds[i].find('<') != -1:
                inst = MEMEInstance()
                inst._seqname (self._current_seq)
                inst._start (start)
                r = re.compile('\d+')
                mn = r.findall(ds[i])[0]
                if ds[i].find('-') != -1:
                    inst.strand = '-'
                else:
                    inst.strand = '+'
                motif = self.data.get_motif_by_name(mn)
                motif.add_instance(inst)
                start += motif.length
            else:
                start += int(ds[i])            
    
    def _add_sequence_match_with_diagram (self, line):
        line = line.strip()
        ls = line.split()
        self.data._add_sequence(ls[0])
        self.data._add_diagram_for_sequence(ls[2],ls[0])
        ds = ls[2].split('_')
        i = 0
        start = 0
        for i in range(0,len(ds)):
            if ds[i].find('+') != -1 or ds[i].find('-') != -1:
                inst = MEMEInstance()
                inst._seqname (ls[0])
                inst._start (start)
                r = re.compile('\d+')
                mn = r.findall(ds[i])[0]
                if ds[i].find('-') != -1:
                    inst.strand = '-'
                else:
                    inst.strand = '+'
                motif = self.data.get_motif_by_name(mn)
                motif.add_instance(inst)
                start += motif.length
            else:
                start += int(ds[i])            
    
    def _add_diagram_from_buffer (self, dummy):
        line = ""
        for l in self._line_buffer:
            line += l.strip()
        ls = line.split()
        self.data._add_diagram_for_sequence(ls[1], self._current_seq)
        ds = ls[1].split('_')
        i = 0
        start = 0
        for i in range(0,len(ds)):
            if ds[i].find('[') != -1 or ds[i].find('<') != -1:
                inst = MEMEInstance()
                inst._seqname (self._current_seq)
                inst._start (start)
                r = re.compile('\d+')
                mn = r.findall(ds[i])[0]
                if ds[i].find('-') != -1:
                    inst.strand = '-'
                else:
                    inst.strand = '+'
                motif = self.data.get_motif_by_name(mn)
                motif.add_instance(inst)
                start += motif.length
            else:
                start += int(ds[i])            
    
    def _set_current_seq (self, line):
        line = line.strip()
        self._current_seq = line
        if not self.data.sequences.count(line):
            self.data.sequences.append(line)
    
    def _add_line_to_buffer (self, line):
        line = line.strip()
        if not line.startswith('*****'):
            self._line_buffer.append(line)
        else:
            return -1
    
    def _parse_buffer (self, dummy):
        """Parses the line buffer to get e-values for each instance of a motif.
        This buffer parser is the most likely point of failure for the 
        MASTParser.
        """
        insts = self.data.get_motif_matches_for_sequence(self._current_seq)    
        if len(insts) > 0:
            
            fullSeq = self._line_buffer[self._buffer_size-1]
            pvals = self._line_buffer[1].split()
            p = 0
            lpval = len(pvals)
            while p < lpval:
                if pvals[p].count('e') > 1:
                #Break blocks up by e and parse into valid floats. This only 
                #works if there are no e-values greater than 1e-5.
                    pvs = []
                    spe = pvals[p].split('e')
                    spe.reverse()
                    dotind = spe[1].find('.')
                    if dotind == -1:
                        thispval = spe[1][-1] + 'e' + spe[0]
                    else:
                        thispval = spe[1][dotind-1:] + 'e' + spe[0]
                    pvs.append(thispval)
                    for spi in range(2,len(spe)):
                        dotind = spe[spi].find('.')
                        prevdotind = spe[spi-1].find('.')
                        if dotind != -1:
                            if prevdotind == -1:
                                thispval = spe[spi][dotind-1:] + 'e' + spe[spi-1][:-1]
                            else:
                                thispval = spe[spi][dotind-1:] + 'e' + spe[spi-1][0:prevdotind-1]
                        else:
                            if prevdotind == -1:
                                thispval = spe[spi][-1] + 'e' + spe[spi-1][:-1]
                            else:
                                thispval = spe[spi][-1] + 'e' + spe[spi-1][0:prevdotind-1]
                        pvs.append(thispval)
                    pvs.reverse()
                    if p > 0:
                        pvals = pvals[0:p] + pvs + pvals[p+1:]
                    else:
                        pvals = pvs + pvals[p+1:]
                    lpval = len(pvals)
                p += 1
            i = 0
            if len(pvals) != len(insts):
                sys.stderr.write("Failure to parse p-values for " + self._current_seq +  ":  " + self._line_buffer[1] + " to: " + str(pvals) + "\n")
                pvals = []
#            else:
#                sys.stderr.write('These are just fine' + self._current_seq + ': ' + self._line_buffer[1] + " to: " + str(pvals) + "\n")
            for i in range(0,len(insts)):
                inst = insts[i]
                start = inst.start - self._buffered_seq_start + 1
                thisSeq = fullSeq[start:start+inst.length]
                thisSeq = Seq.Seq(thisSeq, self.data.alphabet)
                inst._sequence(thisSeq)
                if pvals:
                    inst._pvalue(float(pvals[i]))

    def _blank_buffer (self, dummy):
        self._line_buffer = []
        self._buffer_size = 0
    
    def _collapse_buffer(self, dummy):
        if self._buffer_size == 0:
            if len(self._line_buffer) > 0:
                self._buffer_size = len(self._line_buffer)
                ll = self._line_buffer[self._buffer_size-1].split()
                self._line_buffer[self._buffer_size-1] = ll[1]
                self._buffered_seq_start = int(ll[0])
        else:
            i = 0
            for i in range(self._buffer_size, len(self._line_buffer)-1):
                    self._line_buffer[i-self._buffer_size] = self._line_buffer[i-self._buffer_size] + self._line_buffer[i].strip()
            ll = self._line_buffer[len(self._line_buffer)-1].split()
            if int(ll[0]) == self._buffered_seq_start + len(self._line_buffer[self._buffer_size-1]):
                self._line_buffer[self._buffer_size-1] += ll[1]
            else:
                differ = int(ll[0]) - (self._buffered_seq_start + len(self._line_buffer[self._buffer_size-1]))
                self._line_buffer[self._buffer_size-1] += "N"*differ
                self._line_buffer[self._buffer_size-1] += ll[1]
            self._line_buffer = self._line_buffer[0:self._buffer_size]
    
    def _add_motif_match (self, line):
        line = line.strip()
        if line.find('[') != -1 or line.find('<') != -1:
            pass
        elif line.find('e') != -1:
            pass
        elif line.find('+') != -1:
            pass
    
    def noevent (self, line):
        pass
    

class MASTParser(AbstractParser):
    """
    Parser for MAST text output (OBSOLETE).
    HTML output cannot be parsed, yet. Returns a MASTRecord
    
    A MASTParser takes a file handle for a MAST text output file and 
    returns a MASTRecord, containing the hits between motifs and 
    sequences. The parser does some unusual line buffering to parse out 
    match diagrams. Really complex diagrams often lead to an error message 
    and p-values not being parsed for a given line.
    
    Methods:
    parse (handle): parses the data from the file handle passed to it.
    
    Example:
    
    >>>f = open("mast_file.txt")
    >>>parser = MASTParser()
    >>>mast_record = parser.parse(f)
    >>>for motif in mast_record.motifs:
    ...    for instance in motif.instances:
    ...        print instance.motif_name, instance.sequence_name, instance.strand, instance.pvalue
    
    This class is OBSOLETE; please use the read() function in the module
    Bio.Motif.Parsers.MAST instead.
    """
    def __init__ (self):
        self._consumer = _MASTConsumer()
        self._scanner = _MASTScanner()
    
    def parse (self, handle):
        self._scanner.feed(handle, self._consumer)
        return self._consumer.data
    


class _MASTScanner:
    """
    Scanner for MAST text output (OBSOLETE).
    
    This class is OBSOLETE; please use the read() function in the module
    Bio.Motif.Parsers.MAST instead.
    """
    def feed (self, handle, consumer):
        if isinstance(handle, File.UndoHandle):
            uhandle = handle
        else:
            uhandle = File.UndoHandle(handle)
            
        self._scan_header(uhandle, consumer)
        self._scan_matches(uhandle, consumer)
        self._scan_annotated_matches(uhandle, consumer)
    
    def _scan_header (self, uhandle, consumer):
        try:
            read_and_call_until(uhandle, consumer.noevent, contains = "MAST version")
        except ValueError:
            raise ValueError("Improper input file. Does not begin with a line with 'MAST version'")
        read_and_call(uhandle, consumer._version, contains = 'MAST version')
        read_and_call_until(uhandle, consumer.noevent, start = 'DATABASE AND MOTIFS')
        read_and_call(uhandle, consumer.noevent, start = 'DATABASE')
        read_and_call(uhandle, consumer.noevent, start = '****')
        read_and_call(uhandle, consumer._database, contains = 'DATABASE')
        read_and_call_until(uhandle, consumer.noevent, contains = 'MOTIF WIDTH')
        read_and_call(uhandle, consumer.noevent, contains = 'MOTIF')
        read_and_call(uhandle, consumer.noevent, contains = '----')
        read_and_call_until(uhandle, consumer._add_motif, blank = 1)
        read_and_call_until(uhandle, consumer.noevent, start = 'SECTION II:')
    
    def _scan_matches (self, uhandle, consumer):
        read_and_call_until(uhandle, consumer.noevent, start = 'SEQUENCE NAME')
        read_and_call(uhandle, consumer.noevent, start = 'SEQUENCE NAME')
        read_and_call(uhandle, consumer.noevent, start = '---')
#        read_and_call_until(uhandle, consumer._add_sequence_match_with_diagram, blank = 1)
        read_and_call_until(uhandle, consumer.noevent, blank = 1)
        read_and_call(uhandle, consumer.noevent, blank = 1)
    
    def _scan_annotated_matches (self, uhandle, consumer):
        read_and_call_until(uhandle, consumer.noevent, start = 'SECTION III:')
        read_and_call(uhandle, consumer.noevent, start = 'SECTION III:')
        read_and_call_until(uhandle, consumer.noevent, start = '****')
        read_and_call(uhandle, consumer.noevent, start = '****')
        read_and_call_until(uhandle, consumer.noevent, start = '*****')
        read_and_call(uhandle, consumer.noevent)
        read_and_call_while(uhandle, consumer.noevent, blank = 1)
        readMatches = 1
        while readMatches == 1:
            if consumer._current_seq:
                if consumer._buffer_size != 0:
                    consumer._parse_buffer(None)
                consumer._blank_buffer(None)
            read_and_call(uhandle, consumer._set_current_seq)
            read_and_call_until(uhandle, consumer.noevent, start = '  DIAGRAM')
            read_and_call_until(uhandle, consumer._add_line_to_buffer, blank = 1)
            consumer._add_diagram_from_buffer(None)
            consumer._blank_buffer(None)
            read_and_call(uhandle, consumer.noevent, blank = 1)
            while 1:
                line = safe_peekline(uhandle)
                if line.startswith('****'):
                    consumer._parse_buffer(None)
                    readMatches = 0
                    break
                read_and_call_until(uhandle, consumer._add_line_to_buffer, blank = 1)
                read_and_call(uhandle, consumer.noevent, blank = 1)
                consumer._collapse_buffer(None)
                if attempt_read_and_call(uhandle, consumer.noevent, blank = 1):
                    break
                elif attempt_read_and_call(uhandle, consumer.noevent, start = '*****'):
                    consumer._parse_buffer(None)
                    consumer._blank_buffer(None)
                    readMatches = 0
                    break
    


class MASTRecord:
    """The class for holding the results from a MAST run (OBSOLETE).
    
    A MASTRecord holds data about matches between motifs and sequences.
    The motifs held by the MASTRecord are objects of the class MEMEMotif.
    
    Methods:
    get_motif_matches_for_sequence(sequence_name): returns all of the
    motif matches within a given sequence. The matches are objects of
    the class MEMEInstance
    get_motif_matches (motif_name): returns all of the matches for a motif
    in the sequences searched. The matches returned are of class 
    MEMEInstance
    get_motif_by_name (motif_name): returns a MEMEMotif with the given
    name.
    
    This class is OBSOLETE; please use the read() function in the module
    Bio.Motif.Parsers.MAST instead.
    """
    def __init__ (self):
        self.sequences = []
        self.version = ""
        self.matches = []
        self.database = ""
        self.diagrams = {}
        self.alphabet = None
        self.motifs = []
    
    def _version (self, version):
        self.version = version
    
    def _alphabet (self, alphabet):
        if alphabet == IUPAC.protein or alphabet == IUPAC.ambiguous_dna or alphabet == IUPAC.unambiguous_dna:
            self.alphabet = alphabet
        else:
            return -1
    
    def _database(self, database):
        self.database = database
    
    def get_motif_matches_for_sequence (self, seq):
        insts = []
        for m in self.motifs:
            for i in m.instances:
                if i.sequence_name == seq:
                    insts.append(i)
        insts.sort(lambda x,y: cmp(x.start, y.start))
        return insts
    
    def get_motif_matches (self, motif):
        m = self.get_motif_by_name (motif.name)
        return m.instances
    
    def _add_diagram_for_sequence (self, diagram, seq):
        self.diagrams[seq] = diagram
    
    def _add_match (self, match):
        self.matches.append(match)
    
    def _add_sequence (self, sequence):
        self.sequences.append(sequence)
    
    def _add_motif (self, motif):
        self.motifs.append(motif)
    
    def get_motif_by_name (self, name):
        for m in self.motifs:
            if m.name == name:
                return m
    
