# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Fasta.py

This module provides code to work with FASTA-formatted sequences.


Classes:
Record     Holds information from a FASTA record.
Scanner    Scans a FASTA-format stream.
Consumer   Standard consumer that converts a FASTA record into a Record object.
Iterator   An iterator over records in a FASTA file.

Functions:
parse      Parse a handle with FASTA-formatted data into a list of Records.

"""

# To do:
# Random access for FASTA entries
# index file?

from Bio import File
from Bio.ParserSupport import *

class Record:
    """Holds information from a FASTA record.

    Members:
    title       Title line ('>' character not included).
    sequence    The sequence.
    
    """
    def __init__(self, colwidth=80):
        """__init__(self, colwidth=80)

        Create a new Record.  colwidth specifies the number of residues
        to put on each line when generating FASTA format.

        """
        self.title = ''
        self.sequence = ''
        self._colwidth = colwidth
        
    def __str__(self):
        s = []
        s.append('>%s' % self.title)
        i = 0
        while i < len(self.sequence):
            s.append(self.sequence[i:i+self._colwidth])
            i = i + self._colwidth
        return string.join(s, '\n')

class Scanner:
    """Scans a FASTA-formatted file.

    Methods:
    feed   Feed in one FASTA record.

    """

    def feed(self, handle, consumer):
        """feed(self, handle, consumer)

        Feed in FASTA data for scanning.  handle is a file-like object
        containing FASTA data.  consumer is a Consumer object that will
        recieve events as the FASTA data is scanned.

        """
        if isinstance(handle, File.UndoHandle):
            uhandle = handle
        else:
            uhandle = File.UndoHandle(handle)
        
        if not is_blank_line(uhandle.peekline()):
            self._scan_record(uhandle, consumer)

    def _scan_record(self, uhandle, consumer):
        consumer.start_sequence()
        self._scan_title(uhandle, consumer)
        self._scan_sequence(uhandle, consumer)
        consumer.end_sequence()

    def _scan_title(self, uhandle, consumer):
        read_and_call(uhandle, consumer.title, start='>')

    def _scan_sequence(self, uhandle, consumer):
        while 1:
            line = uhandle.readline()
            if is_blank_line(line):
                break
            elif line[0] == '>':
                uhandle.saveline(line)
                break
            consumer.sequence(line)

class Consumer(AbstractConsumer):
    """Standard consumer that converts a FASTA record to a Record object.

    Members:
    record    Record with FASTA data.

    """
    def __init__(self):
        self.record = None
    
    def start_sequence(self):
        self.record = Record()

    def end_sequence(self):
        pass
    
    def title(self, line):
        assert line[0] == '>'
        self.record.title = string.rstrip(line[1:])

    def sequence(self, line):
        # This can be optimized
        seq = string.rstrip(line)
        self.record.sequence = self.record.sequence + seq

class Iterator:
    """An iterator that returns one record at a time from a FASTA file.

    Methods:
    next   Return the next record from the stream, or None.

    """
    def __init__(self, handle):
        # XXX quick hack to check my code
        assert not isinstance(handle, File.UndoHandle), "oops, legacy stuff"
            
        self._uhandle = File.UndoHandle(handle)
        self._scanner = Scanner()

    def next(self):
        """next(self) -> string

        Return the next FASTA record from the file.  If no more records,
        return the empty string.

        """
        lines = []
        while 1:
            line = self._uhandle.readline()
            if not line:
                break
            if line[0] == '>' and lines:
                self._uhandle.saveline(line)
                break
            lines.append(line)
        return string.join(lines, '')

def parse(handle):
    """parse(handle) -> list of Records

    Parse FASTA-formatted data from handle to a list of Records.
    Warning: this will convert all the data in handle into an in-memory
    data structure.  Do not call this if the amount of data will
    exceed the amount of memory in your machine!

    """
    scanner = Scanner()
    consumer = Consumer()

    iter = Iterator(handle)
    records = []
    while 1:
        r = iter.next()
        if not r:
            break
        scanner.feed(File.StringHandle(r), consumer)
        records.append(consumer.record)
    return records
