# Copyright 2004 by Cymon J. Cox and Frank Kauff.  All rights reserved.
# Copyright 2008 by Michiel de Hoon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
Parser for PHD files output by PHRED and used by PHRAP and CONSED.

This module can be used used directly which will return Record objects
which should contain all the original data in the file.

Alternatively, using Bio.SeqIO with the "phd" format will call this module
internally.  This will give SeqRecord objects for each contig sequence.
"""

from Bio import Seq
from Bio.Alphabet import generic_dna

CKEYWORDS=['CHROMAT_FILE','ABI_THUMBPRINT','PHRED_VERSION','CALL_METHOD',\
        'QUALITY_LEVELS','TIME','TRACE_ARRAY_MIN_INDEX','TRACE_ARRAY_MAX_INDEX',\
        'TRIM','TRACE_PEAK_AREA_RATIO','CHEM','DYE']

class Record:
    """Hold information from a PHD file."""
    def __init__(self):
        self.file_name = ''
        self.comments={}
        for kw in CKEYWORDS:
            self.comments[kw.lower()]=None
        self.sites = []
        self.seq = ''
        self.seq_trimmed = ''


def read(handle):
    """Reads the next PHD record from the file, returning it as a Record object.

    This function reads PHD file data line by line from the handle,
    and returns a single Record object.
    """
    for line in handle:
        if line.startswith("BEGIN_SEQUENCE"):
            record = Record()
            record.file_name = line[15:].rstrip() 
            break
    else:
        return # No record found

    for line in handle:
        if line.startswith("BEGIN_COMMENT"):
            break
    else:
        raise ValueError("Failed to find BEGIN_COMMENT line")
       
    for line in handle:
        line = line.strip()
        if not line:
            continue
        if line=="END_COMMENT":
            break
        keyword, value = line.split(":", 1)
        keyword = keyword.lower()
        value = value.strip()
        if keyword in ('chromat_file',
                       'phred_version',
                       'call_method',
                       'chem',
                       'dye',
                       'time',
                       'basecaller_version',
                       'trace_processor_version'):
            record.comments[keyword] = value
        elif keyword in ('abi_thumbprint',
                         'quality_levels',
                         'trace_array_min_index',
                         'trace_array_max_index'):
            record.comments[keyword] = int(value)
        elif keyword=='trace_peak_area_ratio':
            record.comments[keyword] = float(value)
        elif keyword=='trim':
            first, last, prob = value.split()
            record.comments[keyword] = (int(first), int(last), float(prob))
    else:
        raise ValueError("Failed to find END_COMMENT line")

    for line in handle:
        if line.startswith('BEGIN_DNA'):
            break
    else:
        raise ValueError("Failed to find BEGIN_DNA line")

    for line in handle:
        if line.startswith('END_DNA'):
            break
        else:
            base, quality, location = line.split()
            record.sites.append((base, quality, location))

    for line in handle:
        if line.startswith("END_SEQUENCE"):
            break
    else:
        raise ValueError("Failed to find END_SEQUENCE line")

    record.seq = Seq.Seq(''.join([n[0] for n in record.sites]), generic_dna)
    if record.comments['trim'] is not None:
        first, last = record.comments['trim'][:2]
        record.seq_trimmed = record.seq[first:last]

    return record

def parse(handle):
    """Iterates over a file returning multiple PHD records.

    The data is read line by line from the handle. The handle can be a list
    of lines, an open file, or similar; the only requirement is that we can
    iterate over the handle to retrieve lines from it.

    Typical usage:

    records = parse(handle)
    for record in records:
        # do something with the record object
    """
    while True:
        record = read(handle)
        if not record:
            return
        yield record


# ---------- Everything below is deprecated

from Bio import File
from Bio.ParserSupport import *


class Iterator:
    """Iterates over a file of multiple PHD records (DEPRECATED).
    
    Methods: 
    next    Return the next record from the stream, or None.
    """

    def __init__(self, handle, parser=None):
        """Create a new iterator.
        
        Create a new iterator.  handle is a file-like object.  parser
        is an optional Parser object to change the results into another form.
        If set to None, then the raw contents of the file will be returned.
        """
        import warnings
        warnings.warn("Bio.Sequencing.Phd.Iterator is deprecated. Please use Bio.Sequencing.Phd.parse(handle) instead of Bio.Sequencing.Phd.Iterator(handle, RecordParser())", DeprecationWarning)
        self._uhandle = File.UndoHandle(handle)
        self._parser = parser

    def next(self):
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
        """Iterate over the PHY file record by record."""
        return iter(self.next, None)

class RecordParser(AbstractParser):
    """Parses PHD file data into a Record object (DEPRECATED)."""
    def __init__(self):
        import warnings
        warnings.warn("Bio.Sequencing.Phd.RecordParser is deprecated. Please use Bio.Sequencing.Phd.read(handle) instead of Bio.Sequencing.Phd.RecordParser().parse(handle)", DeprecationWarning)
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
    """Scans a PHD-formatted file (DEPRECATED).
    
    Methods:
    feed - Feed one PHD record.
    """
    def feed(self, handle, consumer):
        """Reads in PDH data from the handle for scanning.

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
    """Consumer that converts a PHD record to a Record object (DEPRECATED)."""
    def __init__(self):
        self.data = None

    def begin_sequence(self, line):
        self.data = Record()
        self.data.file_name = line[15:].rstrip() 

    def end_sequence(self):
        self.data.seq = Seq.Seq(''.join([n[0] for n in self.data.sites]), generic_dna)
        if self.data.comments['trim'] is not None :
            first = self.data.comments['trim'][0]
            last = self.data.comments['trim'][1]
            self.data.seq_trimmed = self.data.seq[first:last]

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

if __name__ == "__main__" :
    print "Quick self test"
    #Test the iterator,
    handle = open("../../Tests/Phd/phd1")
    records = parse(handle)
    for record in records:
        print record.file_name, len(record.seq)
    handle.close()
    print "Done"
