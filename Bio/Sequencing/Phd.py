"""
Parser for PHD files output by PHRED and used by PHRAP and
CONSED.

Works fine with PHRED 0.020425.c

Version 1.1, 03/09/2004
written by Cymon J. Cox (cymon@duke.edu) and Frank Kauff (fkauff@duke.edu)
Comments, bugs, problems, suggestions to one uf us are welcome!

Uses the Biopython Parser interface for parsing: ParserSupport.py

"""

import os
from types import *

from Bio import File
from Bio import Index
from Bio import Seq
from Bio import SeqRecord
from Bio.ParserSupport import *
from Bio.Alphabet import IUPAC

CKEYWORDS=['CHROMAT_FILE','ABI_THUMBPRINT','PHRED_VERSION','CALL_METHOD',\
        'QUALITY_LEVELS','TIME','TRACE_ARRAY_MIN_INDEX','TRACE_ARRAY_MAX_INDEX',\
        'TRIM','TRACE_PEAK_AREA_RATIO','CHEM','DYE']

class Record:
    """Hold information from a PHD file

    """
    def __init__(self):
        self.file_name = ''
        self.comments={}
        for kw in CKEYWORDS:
            self.comments[kw.lower()]=None
        self.sites = []
        self.seq = ''
        self.seq_trimmed = ''


class Iterator:
    """Iterates over a file of multiple PHD records
    
    Methods: 
    next    Return the next record from the stream, or None.
    """

    def __init__(self, handle, parser=None):
        """__init__(self, handle, parser=None)
        
        Create a new iterator.  handle is a file-like object.  parser
        is an optional Parser object to change the results into another form.
        If set to None, then the raw contents of the file will be returned.
        """

        if type(handle) is not FileType and type(handle) is not InstanceType:
            raise ValueError, "I expected a file handle or file-like object"
        self._uhandle = File.UndoHandle(handle)
        self._parser = parser

    def next(self):
        """next(self) -> object

        Return the next PHD record from the file. If no more records
        return None.
        """

        lines = []
        while 1:
            line = self._uhandle.readline()
            if not line:
                break
            # If a new record, then put the line back and stop.
            if lines and line[:14] == 'BEGIN_SEQUENCE':
                self._uhandle.saveline(line)
                break
            lines.append(line)

        if not lines:
            return None

        data = ''.join(lines)
        if self._parser is not None:
            return self._parser.parse(File.StringHandle(data))
        return data
    
    def __iter__(self):
        return iter(self.next, None)

class RecordParser(AbstractParser):
    """Parses PHD file data into a Record object

    """
    def __init__(self):
        self._scanner = _Scanner()
        self._consumer = _RecordConsumer()

    def parse(self, handle):
        if isinstance(handle, File.UndoHandle):
            uhandle = handle
        else:
            uhandle = File.UndoHandle(handle)
        self._scanner.feed(uhandle, self._consumer)
        return self._consumer.data


class _Scanner:
    """Scans a PHD-formatted file
    
    Methods:
    feed - Feed one PHD record.
    """
    def feed(self, handle, consumer):
        """feed(self, handle, consumer)

        Feed in PHD data for scanning.  handle is a file-like object
        containing PHD data.  consumer is a Consumer object that will
        receive events as the PHD data is scanned.
        """
        assert isinstance(handle, File.UndoHandle), \
            "handle must be an UndoHandle"
        if handle.peekline():
            self._scan_record(handle, consumer)

    def _scan_record(self, uhandle, consumer):
        self._scan_begin_sequence(uhandle, consumer)
        self._scan_comments(uhandle, consumer)
        self._scan_dna(uhandle, consumer)
        consumer.end_sequence()

    def _scan_begin_sequence(self, uhandle, consumer):
        read_and_call(uhandle, consumer.begin_sequence, start = 'BEGIN_SEQUENCE')
    
    def _scan_comments(self, uhandle, consumer):
        
        read_and_call_while(uhandle, consumer.noevent, blank=1)
        read_and_call(uhandle, consumer.noevent,  start = 'BEGIN_COMMENT')
        read_and_call_while(uhandle, consumer.noevent, blank=1)
       
        while 1:
            for kw in CKEYWORDS:
                if attempt_read_and_call(uhandle,getattr(consumer,kw.lower()),start=kw+':'):
                    break   # recognized keyword: end for loop and do another while
            else:
                break       # no keywords found: end while loop
            
        read_and_call_while(uhandle, consumer.noevent, blank=1)
        read_and_call(uhandle, consumer.noevent, start = 'END_COMMENT')

    def _scan_dna(self, uhandle, consumer):
        while 1:
            line = uhandle.readline()
            if is_blank_line(line) or line == 'BEGIN_DNA\n':
                continue
            elif line == 'END_DNA\n':
                break
            consumer.read_dna(line)
        
        
class _RecordConsumer(AbstractConsumer):
    """Consumer that converts a PHD record to a Record object

    """
    def __init__(self):
        self.data = None

    def begin_sequence(self, line):
        self.data = Record()
        self.data.file_name = line[15:].rstrip() 

    def end_sequence(self):
        self.data.seq = Seq.Seq(''.join([n[0] for n in self.data.sites]), IUPAC.IUPACAmbiguousDNA()) 
        first = self.data.comments['trim'][0]
        last = self.data.comments['trim'][1]
        self.data.seq_trimmed = Seq.Seq(self.data.seq.tostring()[first : last], IUPAC.IUPACAmbiguousDNA())

    def chromat_file(self, line):
        self.data.comments['chromat_file'] = line[13:-1].strip()

    def abi_thumbprint(self, line):
        self.data.comments['abi_thumbprint'] = int(line[15:-1].strip())

    def phred_version(self, line):
        self.data.comments['phred_version'] = line[14:-1].strip()

    def call_method(self, line):
        self.data.comments['call_method'] = line[12:-1].strip()

    def quality_levels(self, line):
        self.data.comments['quality_levels'] = int(line[15:-1].strip())

    def time(self, line):
        self.data.comments['time'] = line[5:-1].strip()
        
    def trace_array_min_index(self, line):
        self.data.comments['trace_array_min_index'] = int(line[22:-1].strip())
        
    def trace_array_max_index(self, line):
        self.data.comments['trace_array_max_index'] = int(line[22:-1].strip())
    
    def trim(self, line):
        first, last, prob = line[5:-1].split()
        self.data.comments['trim'] = (int(first), int(last), float(prob))
    
    def trace_peak_area_ratio(self, line):
        self.data.comments['trace_peak_area_ratio'] = float(line[22:-1].strip())
    
    def chem(self, line):
        self.data.comments['chem'] = line[5:-1].strip()

    def dye(self, line):
        self.data.comments['dye'] = line[4:-1].strip()
        
    def read_dna(self, line):
        base, quality, location = line.split()
        self.data.sites.append((base, quality, location))

