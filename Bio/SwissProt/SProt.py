# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""SProt.py

This module provides code to work with the sprotXX.dat file from
SwissProt.
http://www.expasy.ch/sprot/sprot-top.html


Classes:
Scanner     Scans SwissProt-formatted data.

"""

from Bio.ParserSupport import *

class Scanner:
    """Scans SwissProt-formatted data.

    Tested with:
    Release 37
    Release 38
    """

    def feed(self, handle, consumer):
        """feed(self, handle, consumer)

        Feed in SwissProt data for scanning.  handle is a file-like
        object that contains the keyword information.  consumer is a
        Consumer object that will receive events as the report is scanned.

        """
        ohandle = OopsHandle(handle)

        while not is_blank_line(ohandle.peekline()):
            self._scan_record(ohandle, consumer)

    def _scan_record(self, ohandle, consumer):
        consumer.start_record()
        for fn in self._scan_fns:
            fn(self, ohandle, consumer)

            # In Release 38, ID N33_HUMAN has a DR buried within comments.
            # Check for this and do more comments, if necessary.
            # XXX handle this better
            if fn is self._scan_dr.im_func:
                self._scan_cc(ohandle, consumer)
                self._scan_dr(ohandle, consumer)
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

    def _scan_ac(self, ohandle, consumer):
        self._scan_line('AC', ohandle, consumer.accession, exactly_one=1)
    
    def _scan_dt(self, ohandle, consumer):
        self._scan_line('DT', ohandle, consumer.date, exactly_one=1)
        self._scan_line('DT', ohandle, consumer.date, exactly_one=1)
        self._scan_line('DT', ohandle, consumer.date, exactly_one=1)

    def _scan_de(self, ohandle, consumer):
        self._scan_line('DE', ohandle, consumer.description, one_or_more=1)
    
    def _scan_gn(self, ohandle, consumer):
        self._scan_line('GN', ohandle, consumer.gene_name, any_number=1)
    
    def _scan_os(self, ohandle, consumer):
        self._scan_line('OS', ohandle, consumer.organism_species,
                        one_or_more=1)
    
    def _scan_og(self, ohandle, consumer):
        self._scan_line('OG', ohandle, consumer.organelle, any_number=1)
    
    def _scan_oc(self, ohandle, consumer):
        self._scan_line('OC', ohandle, consumer.organism_classification,
                        one_or_more=1)

    def _scan_reference(self, ohandle, consumer):
        while 1:
            if safe_peekline(ohandle)[:2] != 'RN':
                break
            self._scan_rn(ohandle, consumer)
            self._scan_rp(ohandle, consumer)
            self._scan_rc(ohandle, consumer)
            self._scan_rx(ohandle, consumer)
            self._scan_ra(ohandle, consumer)
            self._scan_rt(ohandle, consumer)
            self._scan_rl(ohandle, consumer)
    
    def _scan_rn(self, ohandle, consumer):
        self._scan_line('RN', ohandle, consumer.reference_number,
                        exactly_one=1)
    
    def _scan_rp(self, ohandle, consumer):
        self._scan_line('RP', ohandle, consumer.reference_position,
                        exactly_one=1)
    
    def _scan_rc(self, ohandle, consumer):
        self._scan_line('RC', ohandle, consumer.reference_comment,
                        any_number=1)
    
    def _scan_rx(self, ohandle, consumer):
        self._scan_line('RX', ohandle, consumer.reference_cross_reference,
                        up_to_one=1)
    
    def _scan_ra(self, ohandle, consumer):
        self._scan_line('RA', ohandle, consumer.reference_author,
                        one_or_more=1)
    
    def _scan_rt(self, ohandle, consumer):
        self._scan_line('RT', ohandle, consumer.reference_title,
                        any_number=1)
    
    def _scan_rl(self, ohandle, consumer):
        self._scan_line('RL', ohandle, consumer.reference_location,
                        one_or_more=1)
    
    def _scan_cc(self, ohandle, consumer):
        self._scan_line('CC', ohandle, consumer.comment, any_number=1)
    
    def _scan_dr(self, ohandle, consumer):
        self._scan_line('DR', ohandle, consumer.database_cross_reference,
                        any_number=1)
    
    def _scan_kw(self, ohandle, consumer):
        self._scan_line('KW', ohandle, consumer.keyword, any_number=1)
    
    def _scan_ft(self, ohandle, consumer):
        self._scan_line('FT', ohandle, consumer.feature_table, any_number=1)
    
    def _scan_sq(self, ohandle, consumer):
        self._scan_line('SQ', ohandle, consumer.sequence_header, exactly_one=1)
    
    def _scan_sequence_data(self, ohandle, consumer):
        self._scan_line('  ', ohandle, consumer.sequence_data, one_or_more=1)
    
    def _scan_terminator(self, ohandle, consumer):
        self._scan_line('//', ohandle, consumer.terminator, exactly_one=1)
    
    _scan_fns = [
        _scan_id,
        _scan_ac,
        _scan_dt,
        _scan_de,
        _scan_gn,
        _scan_os,
        _scan_og,
        _scan_oc,
        _scan_reference,
        _scan_cc,
        _scan_dr,
        _scan_kw,
        _scan_ft,
        _scan_sq,
        _scan_sequence_data,
        _scan_terminator
        ]
