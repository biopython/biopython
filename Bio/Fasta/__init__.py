# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# 2/2003: Iddo: replaced '\n' with os.linesep for compatibility along
# platforms using different EOL characters
"""Fasta package

This module provides code to work with FASTA-formatted sequences.


Classes:
Record             Holds FASTA sequence data.
Iterator           Iterates over sequence data in a FASTA file.
Dictionary         Accesses a FASTA file using a dictionary interface.
RecordParser       Parses FASTA sequence data into a Record object.
SequenceParser     Parses FASTA sequence data into a Sequence object.

_Scanner           Scans a FASTA-format stream.
_RecordConsumer    Consumes FASTA data to a Record object.
_SequenceConsumer  Consumes FASTA data to a Sequence object.


Functions:
index_file         Index a FASTA file for a Dictionary.

"""
from types import *
import string
import os

__all__ = [
    'FastaAlign',
    ]


from Bio import File
from Bio import Index
from Bio import Seq
from Bio import SeqRecord
from Bio.ParserSupport import *
from Bio import Alphabet

class Record:
    """Holds information from a FASTA record.

    Members:
    title       Title line ('>' character not included).
    sequence    The sequence.
    
    """
    def __init__(self, colwidth=60):
        """__init__(self, colwidth=60)

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
        return string.join(s, os.linesep)

class Iterator:
    """Returns one record at a time from a FASTA file.

    Methods:
    next   Return the next record from the stream, or None.

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

        Return the next FASTA record from the file.  If no more records,
        return None.

        """
        lines = []
        while 1:
            line = self._uhandle.readline()
            if not line:
                break
            if not line.strip():   # ignore blank lines
                continue
            if line[0] == '>' and lines:
                self._uhandle.saveline(line)
                break
            lines.append(line)
            
        if not lines:
            return None
            
        data = string.join(lines, '')
        if self._parser is not None:
            return self._parser.parse(File.StringHandle(data))
        return data

    def __iter__(self):
        return iter(self.next, None)

class Dictionary:
    """Accesses a FASTA file using a dictionary interface.

    """
    __filename_key = '__filename'
    
    def __init__(self, indexname, parser=None, filename = None):
        """__init__(self, indexname, parser=None)

        Open a Fasta Dictionary.  indexname is the name of the
        index for the dictionary.  The index should have been created
        using the index_file function.  parser is an optional Parser
        object to change the results into another form.  If set to None,
        then the raw contents of the file will be returned.
        
        filename specifies the name of the file that this index references. 
        This is useful in cases where the file has been moved since indexing.
        If no filename is supplied (the default) the filename stored in the
        index will be used.
        """
        self._index = Index.Index(indexname)
        if filename:
            self._handle = open(filename)
        else:
            self._handle = open(self._index[Dictionary.__filename_key])
        self._parser = parser

    def __len__(self):
        return len(self._index)

    def __getitem__(self, key):
        start, len = self._index[key]
        self._handle.seek(start)
        data = self._handle.read(len)
        if self._parser is not None:
            return self._parser.parse(File.StringHandle(data))
        return data

    def __getattr__(self, name):
        return getattr(self._index, name)

    def keys(self):
        """Override keys to return only valid keys.

        The keys which is called normally on this dictionary will return
        internal values such as '__filename'. This just strips those out
        before returning.
        """
        all_keys = self._index.keys()
        to_remove = []
        for key in all_keys:
            # if the key is an internal element with a '__', strip it
            if key[:2] == '__':
                to_remove.append(key)

        for key in to_remove:
            all_keys.remove(key)

        return all_keys

class RecordParser(AbstractParser):
    """Parses FASTA sequence data into a Record object.

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

class SequenceParser:
    """Parses FASTA sequence data into a Sequence object.

    """
    def __init__(self, alphabet = Alphabet.generic_alphabet, title2ids = None):
        """Initialize a Scanner and Sequence Consumer.

        Arguments:
        o alphabet - The alphabet of the sequences to be parsed. If not
        passed, this will be set as generic_alphabet.
        o title2ids - A function that, when given the title of the FASTA
        file (without the beginning >), will return the name, id and
        description (in that order) for the record. If this is not given,
        then the entire title line will be used as the description.
        """
        self._scanner = _Scanner()
        self._consumer = _SequenceConsumer(alphabet, title2ids)

    def parse(self, handle):
        if isinstance(handle, File.UndoHandle):
            uhandle = handle
        else:
            uhandle = File.UndoHandle(handle)
        self._scanner.feed(uhandle, self._consumer)
        return self._consumer.data


class _Scanner:
    """Scans a FASTA-formatted file.

    Methods:
    feed   Feed in one FASTA record.

    """
    def feed(self, handle, consumer):
        """feed(self, handle, consumer)

        Feed in FASTA data for scanning.  handle is a file-like object
        containing FASTA data.  consumer is a Consumer object that will
        receive events as the FASTA data is scanned.

        """
        assert isinstance(handle, File.UndoHandle), \
               "handle must be an UndoHandle"
        if handle.peekline():
            self._scan_record(handle, consumer)

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

class _RecordConsumer(AbstractConsumer):
    """Consumer that converts a FASTA record to a Record object.

    Members:
    data    Record with FASTA data.
    
    """
    def __init__(self):
        self.data = None
    
    def start_sequence(self):
        self.data = Record()

    def end_sequence(self):
        pass
    
    def title(self, line):
        assert line[0] == '>'
        self.data.title = string.rstrip(line[1:])

    def sequence(self, line):
        # This can be optimized
        seq = string.rstrip(line)
        self.data.sequence = self.data.sequence + seq

class _SequenceConsumer(AbstractConsumer):
    """Consumer that converts a FASTA record to a Sequence object.

    Members:
    data    Sequence with FASTA data.

    """
    def __init__(self, alphabet = Alphabet.generic_alphabet, title2ids = None):
        """Initialize the Consumer.

        Arguments:
        o alphabet - The alphabet of the sequences we will be creating.
        o title2ids - A function that will convert the title of a FASTA
        record into the id, name and description information.
        """
        self.data = None
        self.alphabet = alphabet
        self.title2ids = title2ids
    
    def start_sequence(self):
        seq = Seq.Seq('', self.alphabet)
        self.data = SeqRecord.SeqRecord(seq)

    def end_sequence(self):
        pass
    
    def title(self, line):
        assert line[0] == '>'
        if self.title2ids:
            id, name, descr = self.title2ids(string.rstrip(line[1:]))
            self.data.id = id
            self.data.name = name
            self.data.description = descr
        else:
            self.data.description = string.rstrip(line[1:])

    def sequence(self, line):
        # This can be optimized
        seq = Seq.Seq(string.rstrip(line), self.alphabet)
        self.data.seq = self.data.seq + seq

def index_file(filename, indexname, rec2key=None):
    """index_file(filename, indexname, rec2key=None)

    Index a FASTA file.  filename is the name of the file.
    indexname is the name of the dictionary.  rec2key is an
    optional callback that takes a Record and generates a unique key
    (e.g. the accession number) for the record.  If not specified,
    the sequence title will be used.

    """
    if not os.path.exists(filename):
        raise ValueError, "%s does not exist" % filename

    index = Index.Index(indexname, truncate=1)
    index[Dictionary._Dictionary__filename_key] = filename
    
    iter = Iterator(open(filename), parser=RecordParser())
    while 1:
        start = iter._uhandle.tell()
        rec = iter.next()
        length = iter._uhandle.tell() - start
        
        if rec is None:
            break
        if rec2key is not None:
            key = rec2key(rec)
        else:
            key = rec.title
            
        if not key:
            raise KeyError, "empty sequence key was produced"
        elif index.has_key(key):
            raise KeyError, "duplicate key %s found" % key

        index[key] = start, length
