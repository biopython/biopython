# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""KeyWList.py

This module provides code to work with the keywlist.txt file from
SwissProt.
http://www.expasy.ch/sprot/sprot-top.html


Classes:
ListParser        Parses a keywlist.txt file into a list of keywords.

_Scanner          Scans the keywlist.txt file.
_ListConsumer     Consumes keywlist data to a list.


Functions:
extract_keywords  Return the keywords from a keywlist.txt file.

"""

from types import *

from Bio import File
from Bio.ParserSupport import *

class ListParser(AbstractParser):
    """Parses keywlist.txt data into a list of keywords.

    """
    def __init__(self):
        self._scanner = _Scanner()
        self._consumer = _ListConsumer()

    def parse(self, handle):
        self._scanner.feed(handle, self._consumer)
        return self._consumer.keywords


class _Scanner:
    """Scan the keywlist.txt file included with the SwissProt distribution.

    Tested with:
    Release 37
    Release 38
    """

    def feed(self, handle, consumer):
        """feed(self, handle, consumer)

        Feed in the keywlist.txt file for scanning.  handle is a file-like
        object that contains keyword information.  consumer is a
        Consumer object that will receive events as the report is scanned.

        """
        if isinstance(handle, File.UndoHandle):
            uhandle = handle
        else:
            uhandle = File.UndoHandle(handle)
        
        self._scan_header(uhandle, consumer)
        self._scan_keywords(uhandle, consumer)
        self._scan_footer(uhandle, consumer)

    def _scan_header(self, uhandle, consumer):
        consumer.start_header()
        
        read_and_call(uhandle, consumer.noevent, start='----')
        read_and_call(uhandle, consumer.noevent, blank=1)
        read_and_call(uhandle, consumer.noevent, contains="SWISS-PROT")
        read_and_call(uhandle, consumer.noevent, contains="Release")
        read_and_call(uhandle, consumer.noevent, blank=1)
        read_and_call(uhandle, consumer.noevent, start='----')

        read_and_call(uhandle, consumer.noevent, blank=1)
        read_and_call(uhandle, consumer.noevent, start='List of keywords')
        read_and_call(uhandle, consumer.noevent, blank=1)
        read_and_call(uhandle, consumer.noevent, start='----')

        while 1:
            if attempt_read_and_call(uhandle, consumer.noevent, start='----'):
                break
            read_and_call(uhandle, consumer.noevent, blank=0)

        read_and_call(uhandle, consumer.noevent, start='Document name')
        read_and_call(uhandle, consumer.noevent, start='----')
        read_and_call(uhandle, consumer.noevent, blank=1)
        
        consumer.end_header()

    def _scan_keywords(self, uhandle, consumer):
        consumer.start_keywords()

        # SwissProt38 starts with lines:
        # Keyword
        # ______________________________________
        #
        # Check and see if it's release 38, and parse it.
        if attempt_read_and_call(uhandle, consumer.noevent, start='Keyword'):
            read_and_call(uhandle, consumer.noevent, start='____')

        while 1:
            if not attempt_read_and_call(uhandle, consumer.keyword, blank=0):
                break
        read_and_call(uhandle, consumer.noevent, blank=1)
        
        consumer.end_keywords()

    def _scan_footer(self, uhandle, consumer):
        consumer.start_footer()

        read_and_call(uhandle, consumer.noevent, start='----')
        while 1:
            if attempt_read_and_call(uhandle, consumer.noevent, start='----'):
                break
            read_and_call(uhandle, consumer.copyright, blank=0)

        consumer.end_footer()

class _ListConsumer(AbstractConsumer):
    """Consumer that converts a keywlist.txt file into a list of keywords.

    Members:
    keywords    List of keywords.

    """
    def __init__(self):
        self.keywords = None

    def start_keywords(self):
        self.keywords = []

    def keyword(self, line):
        self.keywords.append(string.rstrip(line))

def extract_keywords(keywlist_handle):
    """extract_keywords(keywlist_handle) -> list of keywords

    Return the keywords from a keywlist.txt file.

    """
    if type(keywlist_handle) is not FileType and \
       type(keywlist_handle) is not InstanceType:
        raise ValueError, "I expected a file handle or file-like object"
    return ListParser().parse(keywlist_handle)
