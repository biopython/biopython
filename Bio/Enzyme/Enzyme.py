# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Enzyme.py

This module provides code to work with the enzyme.dat file from
Enzyme.
http://www.expasy.ch/enzyme/


Classes:
Scanner     Scans Enzyme data.

"""

from Bio.ParserSupport import *

class Scanner:
    """Scans Enzyme data.

    Tested with:
    XXX ??
    """

    def feed(self, handle, consumer):
        """feed(self, handle, consumer)

        Feed in Enzyme data for scanning.  handle is a file-like object
        that contains the keyword information.  consumer is a Consumer
        object that will receive events as the report is scanned.

        """
        ohandle = OopsHandle(handle)
        while not is_blank_line(ohandle.peekline()):   # Am I done yet?
            self._scan_record(ohandle, consumer)

    def _scan_record(self, ohandle, consumer):
        # The first record is just copyright information embedded in
        # comments.  Check to see if I'm at the first record.  If so,
        # then just scan the comments and the terminator.
        consumer.start_record()
        line = ohandle.peekline()
        if line[:2] == 'CC':
            self._scan_cc(ohandle, consumer)
            self._scan_terminator(ohandle, consumer)
        else:
            for fn in self._scan_fns:
                fn(self, ohandle, consumer)
        consumer.end_record()

    def _scan_line(self, line_type, ohandle, event_fn,
                   exactly_one=None, one_or_more=None, any_number=None,
                   up_to_one=None):
        # Callers must set exactly one of exactly_one, one_or_more, or
        # any_number to a true value.  I do not explicitly check to
        # make sure this function is called correctly.
        
        # This does not guarantee any parameter safety, but I
        # like the readability.  The other strategy I tried was have
        # parameters min_lines, max_lines.
        
        if exactly_one or one_or_more:
            read_and_call(ohandle, event_fn, start=line_type)
        if one_or_more or any_number:
            while 1:
                if not attempt_read_and_call(ohandle, event_fn,
                                             start=line_type):
                    break
        if up_to_one:
            attempt_read_and_call(ohandle, event_fn, start=line_type)

    def _scan_id(self, ohandle, consumer):
        self._scan_line('ID', ohandle, consumer.identification, exactly_one=1)

    def _scan_de(self, ohandle, consumer):
        self._scan_line('DE', ohandle, consumer.description, one_or_more=1)
    
    def _scan_an(self, ohandle, consumer):
        self._scan_line('AN', ohandle, consumer.alternate_name, any_number=1)
    
    def _scan_ca(self, ohandle, consumer):
        self._scan_line('CA', ohandle, consumer.catalytic_activity,
                        any_number=1)
    
    def _scan_cf(self, ohandle, consumer):
        self._scan_line('CF', ohandle, consumer.cofactor, any_number=1)

    def _scan_cc(self, ohandle, consumer):
        self._scan_line('CC', ohandle, consumer.comment, any_number=1)
    
    def _scan_di(self, ohandle, consumer):
        self._scan_line('DI', ohandle, consumer.disease, any_number=1)
    
    def _scan_pr(self, ohandle, consumer):
        self._scan_line('PR', ohandle, consumer.prosite_reference,
                        any_number=1)
    
    def _scan_dr(self, ohandle, consumer):
        self._scan_line('DR', ohandle, consumer.databank_reference,
                        any_number=1)

    def _scan_terminator(self, ohandle, consumer):
        self._scan_line('//', ohandle, consumer.terminator, exactly_one=1)
    
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
