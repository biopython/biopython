# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Enzyme.py

This module provides code to work with the enzyme.dat file from
Enzyme.
http://www.expasy.ch/enzyme/


Classes:
_Scanner     Scans Enzyme data.

"""

from Bio import File
from Bio.ParserSupport import *

class _Scanner:
    """Scans Enzyme data.

    Tested with:
    Release 24
    """

    def feed(self, handle, consumer):
        """feed(self, handle, consumer)

        Feed in Enzyme data for scanning.  handle is a file-like object
        that contains keyword information.  consumer is a Consumer
        object that will receive events as the report is scanned.

        """
        if isinstance(handle, File.UndoHandle):
            uhandle = handle
        else:
            uhandle = File.UndoHandle(handle)

        while not is_blank_line(uhandle.peekline()):   # Am I done yet?
            self._scan_record(uhandle, consumer)

    def _scan_record(self, uhandle, consumer):
        # The first record is just copyright information embedded in
        # comments.  Check to see if I'm at the first record.  If so,
        # then just scan the comments and the terminator.
        consumer.start_record()
        line = uhandle.peekline()
        if line[:2] == 'CC':
            self._scan_cc(uhandle, consumer)
            self._scan_terminator(uhandle, consumer)
        else:
            for fn in self._scan_fns:
                fn(self, uhandle, consumer)
        consumer.end_record()

    def _scan_line(self, line_type, uhandle, event_fn,
                   exactly_one=None, one_or_more=None, any_number=None,
                   up_to_one=None):
        # Callers must set exactly one of exactly_one, one_or_more, or
        # any_number to a true value.  I do not explicitly check to
        # make sure this function is called correctly.
        
        # This does not guarantee any parameter safety, but I
        # like the readability.  The other strategy I tried was have
        # parameters min_lines, max_lines.
        
        if exactly_one or one_or_more:
            read_and_call(uhandle, event_fn, start=line_type)
        if one_or_more or any_number:
            while 1:
                if not attempt_read_and_call(uhandle, event_fn,
                                             start=line_type):
                    break
        if up_to_one:
            attempt_read_and_call(uhandle, event_fn, start=line_type)

    def _scan_id(self, uhandle, consumer):
        self._scan_line('ID', uhandle, consumer.identification, exactly_one=1)

    def _scan_de(self, uhandle, consumer):
        self._scan_line('DE', uhandle, consumer.description, one_or_more=1)
    
    def _scan_an(self, uhandle, consumer):
        self._scan_line('AN', uhandle, consumer.alternate_name, any_number=1)
    
    def _scan_ca(self, uhandle, consumer):
        self._scan_line('CA', uhandle, consumer.catalytic_activity,
                        any_number=1)
    
    def _scan_cf(self, uhandle, consumer):
        self._scan_line('CF', uhandle, consumer.cofactor, any_number=1)

    def _scan_cc(self, uhandle, consumer):
        self._scan_line('CC', uhandle, consumer.comment, any_number=1)
    
    def _scan_di(self, uhandle, consumer):
        self._scan_line('DI', uhandle, consumer.disease, any_number=1)
    
    def _scan_pr(self, uhandle, consumer):
        self._scan_line('PR', uhandle, consumer.prosite_reference,
                        any_number=1)
    
    def _scan_dr(self, uhandle, consumer):
        self._scan_line('DR', uhandle, consumer.databank_reference,
                        any_number=1)

    def _scan_terminator(self, uhandle, consumer):
        self._scan_line('//', uhandle, consumer.terminator, exactly_one=1)
    
    _scan_fns = [
        _scan_id,
        _scan_de,
        _scan_an,
        _scan_ca,
        _scan_cf,
        _scan_cc,
        _scan_di,
        _scan_pr,
        _scan_dr,
        _scan_terminator
        ]
