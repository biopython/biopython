# Copyright 2000 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
This module is OBSOLETE.
Most of the functionality in this module has moved to Bio.ExPASy.Prodoc;
please see

Bio.ExPASy.Prodoc.read          To read a Prodoc file containing one entry.
Bio.ExPASy.Prodoc.parse         Iterates over entries in a Prodoc file.
Bio.ExPASy.Prodoc.Record        Holds Prodoc data.
Bio.ExPASy.Prodoc.Reference     Holds data from a Prodoc reference.

The other functions and classes in Bio.Prosite.Prodoc (including
Bio.Prosite.Prodoc.index_file and Bio.Prosite.Prodoc.Dictionary) are
considered deprecated, and were not moved to Bio.ExPASy.Prodoc. If you use
this functionality, please contact the Biopython developers at
biopython-dev@biopython.org to avoid permanent removal of this module from
Biopython.




This module provides code to work with the prosite.doc file from
Prosite, available at http://www.expasy.ch/prosite/.

Tested with:
Release 15.0, July 1998
Release 16.0, July 1999
Release 20.22, 13 November 2007


Functions:
parse              Iterates over entries in a Prodoc file.
index_file         Index a Prodoc file for a Dictionary.
_extract_record    Extract Prodoc data from a web page.


Classes:
Record             Holds Prodoc data.
Reference          Holds data from a Prodoc reference.
Dictionary         Accesses a Prodoc file using a dictionary interface.
RecordParser       Parses a Prodoc record into a Record object.

_Scanner           Scans Prodoc-formatted data.
_RecordConsumer    Consumes Prodoc data to a Record object.
"""

import warnings
warnings.warn("This module is OBSOLETE. Most of the functionality in this module has moved to Bio.ExPASy.Prodoc.", PendingDeprecationWarning)

from types import *
import os
import sgmllib
from Bio import File
from Bio import Index
from Bio.ParserSupport import *

def parse(handle):
    import cStringIO
    parser = RecordParser()
    text = ""
    for line in handle:
        text += line
        if line[:5] == '{END}':
            handle = cStringIO.StringIO(text)
            record = parser.parse(handle)
            text = ""
            yield record

def read(handle):
    parser = RecordParser()
    record = parser.parse(handle)
    # We should have reached the end of the record by now
    remainder = handle.read()
    if remainder:
        raise ValueError("More than one Prodoc record found")
    return record


# It may be a good idea to rewrite read(), parse() at some point to avoid
# using the old-style "parser = RecordParser(); parser.parse(handle)" approach.

class Record:
    """Holds information from a Prodoc record.

    Members:
    accession      Accession number of the record.
    prosite_refs   List of tuples (prosite accession, prosite name).
    text           Free format text.
    references     List of reference objects.

    """
    def __init__(self):
        self.accession = ''
        self.prosite_refs = []
        self.text = ''
        self.references = []

class Reference:
    """Holds information from a Prodoc citation.

    Members:
    number     Number of the reference. (string)
    authors    Names of the authors.
    citation   Describes the citation.

    """
    def __init__(self):
        self.number = ''
        self.authors = ''
        self.citation = ''

class Dictionary:
    """Accesses a Prodoc file using a dictionary interface.

    """
    __filename_key = '__filename'
    
    def __init__(self, indexname, parser=None):
        """__init__(self, indexname, parser=None)

        Open a Prodoc Dictionary.  indexname is the name of the
        index for the dictionary.  The index should have been created
        using the index_file function.  parser is an optional Parser
        object to change the results into another form.  If set to None,
        then the raw contents of the file will be returned.

        """
        self._index = Index.Index(indexname)
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

class RecordParser(AbstractParser):
    """Parses Prodoc data into a Record object.

    """
    def __init__(self):
        self._scanner = _Scanner()
        self._consumer = _RecordConsumer()

    def parse(self, handle):
        self._scanner.feed(handle, self._consumer)
        return self._consumer.data

class _Scanner:
    """Scans Prodoc-formatted data.

    Tested with:
    Release 15.0, July 1998
    
    """
    def feed(self, handle, consumer):
        """feed(self, handle, consumer)

        Feed in Prodoc data for scanning.  handle is a file-like
        object that contains prosite data.  consumer is a
        Consumer object that will receive events as the report is scanned.

        """
        if isinstance(handle, File.UndoHandle):
            uhandle = handle
        else:
            uhandle = File.UndoHandle(handle)

        while 1:
            line = uhandle.peekline()
            if not line:
                break
            elif is_blank_line(line):
                # Skip blank lines between records
                uhandle.readline()
                continue
            else:
                self._scan_record(uhandle, consumer)
            
    def _scan_record(self, uhandle, consumer):
        consumer.start_record()

        self._scan_accession(uhandle, consumer)
        self._scan_prosite_refs(uhandle, consumer)
        read_and_call(uhandle, consumer.noevent, start='{BEGIN}')
        self._scan_text(uhandle, consumer)
        self._scan_refs(uhandle, consumer)
        self._scan_copyright(uhandle, consumer)
        read_and_call(uhandle, consumer.noevent, start='{END}')

        consumer.end_record()

    def _scan_accession(self, uhandle, consumer):
        read_and_call(uhandle, consumer.accession, start='{PDOC')

    def _scan_prosite_refs(self, uhandle, consumer):
        while attempt_read_and_call(uhandle, consumer.prosite_reference,
                                    start='{PS'):
            pass

    def _scan_text(self, uhandle, consumer):
        while 1:
            line = safe_readline(uhandle)
            if (line[0] == '[' and line[3] == ']' and line[4] == ' ') or \
               line[:5] == '{END}':
                uhandle.saveline(line)
                break
            consumer.text(line)

    def _scan_refs(self, uhandle, consumer):
        while 1:
            line = safe_readline(uhandle)
            if line[:5] == '{END}' or is_blank_line(line):
                uhandle.saveline(line)
                break
            consumer.reference(line)

    def _scan_copyright(self, uhandle, consumer):
        # Cayte Lindner found some PRODOC records with the copyrights
        # appended at the end.  We'll try and recognize these.
        read_and_call_while(uhandle, consumer.noevent, blank=1)
        if attempt_read_and_call(uhandle, consumer.noevent, start='+----'):
            read_and_call_until(uhandle, consumer.noevent, start='+----')
            read_and_call(uhandle, consumer.noevent, start='+----')
        read_and_call_while(uhandle, consumer.noevent, blank=1)

class _RecordConsumer(AbstractConsumer):
    """Consumer that converts a Prodoc record to a Record object.

    Members:
    data    Record with Prodoc data.

    """
    def __init__(self):
        self.data = None
        
    def start_record(self):
        self.data = Record()
        
    def end_record(self):
        self._clean_data()

    def accession(self, line):
        line = line.rstrip()
        if line[0] != '{' or line[-1] != '}':
            raise ValueError("I don't understand accession line\n%s" % line)
        acc = line[1:-1]
        if acc[:4] != 'PDOC':
            raise ValueError("Invalid accession in line\n%s" % line)
        self.data.accession = acc

    def prosite_reference(self, line):
        line = line.rstrip()
        if line[0] != '{' or line[-1] != '}':
            raise ValueError("I don't understand accession line\n%s" % line)
        acc, name = line[1:-1].split('; ')
        self.data.prosite_refs.append((acc, name))
    
    def text(self, line):
        self.data.text = self.data.text + line
    
    def reference(self, line):
        if line[0] == '[' and line[3] == ']':  # new reference
            self._ref = Reference()
            self._ref.number = line[1:3].strip()
            if line[1] == 'E':
                # If it's an electronic reference, then the URL is on the
                # line, instead of the author.
                self._ref.citation = line[4:].strip()
            else:
                self._ref.authors = line[4:].strip()
            self.data.references.append(self._ref)
        elif line[:4] == '    ':
            if not self._ref:
                raise ValueError("Unnumbered reference lines\n%s" % line)
            self._ref.citation = self._ref.citation + line[5:]
        else:
            raise Exception("I don't understand the reference line\n%s" % line)

    def _clean_data(self):
        # get rid of trailing newlines
        for ref in self.data.references:
            ref.citation = ref.citation.rstrip()
            ref.authors = ref.authors.rstrip()
    
def index_file(filename, indexname, rec2key=None):
    """index_file(filename, indexname, rec2key=None)

    Index a Prodoc file.  filename is the name of the file.
    indexname is the name of the dictionary.  rec2key is an
    optional callback that takes a Record and generates a unique key
    (e.g. the accession number) for the record.  If not specified,
    the id name will be used.

    """
    import os
    if not os.path.exists(filename):
        raise ValueError("%s does not exist" % filename)

    index = Index.Index(indexname, truncate=1)
    index[Dictionary._Dictionary__filename_key] = filename

    handle = open(filename)
    records = parse(handle)
    end = 0L
    for record in records:
        start = end
        end = handle.tell()
        length = end - start

        if rec2key is not None:
            key = rec2key(record)
        else:
            key = record.accession
            
        if not key:
            raise KeyError("empty key was produced")
        elif key in index:
            raise KeyError("duplicate key %s found" % key)

        index[key] = start, length
