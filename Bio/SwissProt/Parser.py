# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""SwissProt Parser

This module provides code to parse file from SwissProt.
http://www.expasy.ch/sprot/sprot-top.html


Classes:
KeyWListScanner   Scans the keywlist.txt file.

"""

from Bio.ParserSupport import *

class KeyWListScanner:
    """Scan the keywlist.txt file included with the SwissProt distribution.

    This scanner produces the following events:
    header
    keywords
        keyword
    footer
        copyright

    Tested with:
    Release 37
    Release 38
    """

    def feed(self, handle, consumer):
        """feed(self, handle, consumer)

        Feed in the keywlist.txt file for scanning.  handle is a file-like
        object that contains the keyword information.  consumer is a
        Consumer object that will receive events as the report is scanned.

        """
        ohandle = OopsHandle(handle)

        self._scan_header(ohandle, consumer)
        self._scan_keywords(ohandle, consumer)
        self._scan_footer(ohandle, consumer)

    def _scan_header(self, ohandle, consumer):
        consumer.start_header()
        
        read_and_call(ohandle, consumer.noevent, start='----')
        read_and_call(ohandle, consumer.noevent, blank=1)
        read_and_call(ohandle, consumer.noevent, contains="SWISS-PROT")
        read_and_call(ohandle, consumer.noevent, contains="Release")
        read_and_call(ohandle, consumer.noevent, blank=1)
        read_and_call(ohandle, consumer.noevent, start='----')

        read_and_call(ohandle, consumer.noevent, blank=1)
        read_and_call(ohandle, consumer.noevent, start='List of keywords')
        read_and_call(ohandle, consumer.noevent, blank=1)
        read_and_call(ohandle, consumer.noevent, start='----')

        while 1:
            if attempt_read_and_call(ohandle, consumer.noevent, start='----'):
                break
            read_and_call(ohandle, consumer.noevent, blank=0)

        read_and_call(ohandle, consumer.noevent, start='Document name')
        read_and_call(ohandle, consumer.noevent, start='----')
        read_and_call(ohandle, consumer.noevent, blank=1)
        
        consumer.end_header()

    def _scan_keywords(self, ohandle, consumer):
        consumer.start_keywords()

        # SwissProt38 starts with lines:
        # Keyword
        # ______________________________________
        #
        # Check and see if it's release 38, and parse it.
        if attempt_read_and_call(ohandle, consumer.noevent, start='Keyword'):
            read_and_call(ohandle, consumer.noevent, start='____')

        while 1:
            if not attempt_read_and_call(ohandle, consumer.keyword, blank=0):
                break
        read_and_call(ohandle, consumer.noevent, blank=1)
        
        consumer.end_keywords()

    def _scan_footer(self, ohandle, consumer):
        consumer.start_footer()

        read_and_call(ohandle, consumer.noevent, start='----')
        while 1:
            if attempt_read_and_call(ohandle, consumer.noevent, start='----'):
                break
            read_and_call(ohandle, consumer.copyright, blank=0)

        consumer.end_footer()
