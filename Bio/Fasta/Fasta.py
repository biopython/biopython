# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Fasta.py

This module provides code to work with FASTA-formatted sequences.


Classes:
Record     Holds information from a FASTA record.
Scanner    Scans a FASTA-format file.
Consumer   Standard consumer that converts a FASTA stream to a list of Records.

Functions:
parse      Parse a handle with FASTA-formatted data into a list of Records.

"""

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
    feed     Feed in FASTA data for scanning.

    """

    def feed(self, handle, consumer):
        """feed(self, handle, consumer)

        Feed in FASTA data for scanning.  handle is a file-like object
        that contains the keyword information.  consumer is a Consumer
        object that will receive events as the report is scanned.

        """
        ohandle = OopsHandle(handle)
        while not is_blank_line(ohandle.peekline()):   # Am I done yet?
            self._scan_record(ohandle, consumer)

    def _scan_record(self, ohandle, consumer):
        consumer.start_sequence()
        self._scan_title(ohandle, consumer)
        self._scan_sequence(ohandle, consumer)
        consumer.end_sequence()

    def _scan_title(self, ohandle, consumer):
        read_and_call(ohandle, consumer.title, start='>')

    def _scan_sequence(self, ohandle, consumer):
        while 1:
            line = ohandle.readline()
            if is_blank_line(line):
                break
            elif line[0] == '>':
                ohandle.saveline(line)
                break
            consumer.sequence(line)

class Consumer(AbstractConsumer):
    """Standard consumer that converts a FASTA stream to a list of Records.

    Members:
    records     List of Records.

    """
    def __init__(self):
        self.records = []
    
    def start_sequence(self):
        self.records.append(Record())

    def end_sequence(self):
        pass
    
    def title(self, line):
        assert line[0] == '>'
        self.records[-1].title = string.rstrip(line[1:])

    def sequence(self, line):
        # This can be optimized
        seq = string.rstrip(line)
        self.records[-1].sequence = self.records[-1].sequence + seq
    

def parse(handle):
    """parse(handle) -> list of Records

    Parse FASTA-formatted data from handle to a list of Records.
    Warning: this will convert all the data in handle into an in-memory
    data structure.  Do not call this if the amount of data will
    exceed the amount of memory in your machine!

    """
    c = Consumer()
    Scanner().feed(handle, c)
    return c.records
    
