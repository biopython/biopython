# Copyright 2008 by Bartek Wilczynski
# Adapted from  Bio.MEME.Parser by Jason A. Hackney.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio.Alphabet import IUPAC
from Bio import File
from Bio.ParserSupport import *
from Bio import Seq
import re
import warnings
from Bio.Motif.Parsers.MEME import MEMEMotif, MEMEInstance


def read(handle):
    """read(handle)"""
    parser = MASTParser()
    record = parser.parse(handle)
    return record

class _MASTConsumer:
    """
    Consumer that can receive events from _MASTScanner.
    
    A _MASTConsumer parses lines from a mast text output file.
    The motif match diagrams are parsed using line buffering. 
    Each of the buffering functions have a dummy variable that is
    required for testing using the Bio.ParserSupport.TaggingConsumer.
    If this variable isn't there, the TaggingConsumer barfs. In
    the _MASTScanner, None is passed in the place of this variable.
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
        m.length=int(ls[1])
        name = ls[0]
        m.name=name
        # m.add_instance(ls[2])
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
        return
        # code below doesn't work yet
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
                warnings.warn("Failure to parse p-values for %s: %s to: %s" % (self._current_seq, self._line_buffer[1], str(pvals)))
                pvals = []
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
    Parser for MAST text output. HTML output cannot be parsed, yet. Returns a MASTRecord
    
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
    >>>    for instance in motif.instances:
    >>>        print instance.motif_name, instance.sequence_name, instance.strand, instance.pvalue
    """
    def __init__ (self):
        self._consumer = _MASTConsumer()
        self._scanner = _MASTScanner()
    
    def parse (self, handle):
        self._scanner.feed(handle, self._consumer)
        return self._consumer.data
    


class _MASTScanner:
    """
    Scanner for MAST text output. 
        
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
    """The class for holding the results from a MAST run.
    
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
