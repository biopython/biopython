# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Fasta.py

This module provides code to work with FASTA-formatted sequences.


Classes:
Scanner     Scans a FASTA-format file.

"""

from Bio.ParserSupport import *

class Scanner:
    """Scans a FASTA-formatted file.

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
            if line[0] == '>':
                ohandle.saveline(line)
                break
            elif is_blank_line(line):
                break
            consumer.sequence(line)
