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

from Bio import File
from Bio.ParserSupport import *

class Scanner:
    """Scans SwissProt-formatted data.

    Tested with:
    Release 37
    Release 38
    """

    def feed(self, uhandle, consumer):
        """feed(self, uhandle, consumer)

        Feed in SwissProt data for scanning.  uhandle must be an
        UndoHandle that contains the keyword information.  consumer is a
        Consumer object that will receive events as the report is scanned.

        """
        assert isinstance(uhandle, File.UndoHandle), \
               "uhandle must be an instance of Bio.File.UndoHandle"

        while not is_blank_line(uhandle.peekline()):
            self._scan_record(uhandle, consumer)

    def _scan_record(self, uhandle, consumer):
        consumer.start_record()
        for fn in self._scan_fns:
            fn(self, uhandle, consumer)

            # In Release 38, ID N33_HUMAN has a DR buried within comments.
            # Check for this and do more comments, if necessary.
            # XXX handle this better
            if fn is self._scan_dr.im_func:
                self._scan_cc(uhandle, consumer)
                self._scan_dr(uhandle, consumer)
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

    def _scan_ac(self, uhandle, consumer):
        self._scan_line('AC', uhandle, consumer.accession, exactly_one=1)
    
    def _scan_dt(self, uhandle, consumer):
        self._scan_line('DT', uhandle, consumer.date, exactly_one=1)
        self._scan_line('DT', uhandle, consumer.date, exactly_one=1)
        self._scan_line('DT', uhandle, consumer.date, exactly_one=1)

    def _scan_de(self, uhandle, consumer):
        self._scan_line('DE', uhandle, consumer.description, one_or_more=1)
    
    def _scan_gn(self, uhandle, consumer):
        self._scan_line('GN', uhandle, consumer.gene_name, any_number=1)
    
    def _scan_os(self, uhandle, consumer):
        self._scan_line('OS', uhandle, consumer.organism_species,
                        one_or_more=1)
    
    def _scan_og(self, uhandle, consumer):
        self._scan_line('OG', uhandle, consumer.organelle, any_number=1)
    
    def _scan_oc(self, uhandle, consumer):
        self._scan_line('OC', uhandle, consumer.organism_classification,
                        one_or_more=1)

    def _scan_reference(self, uhandle, consumer):
        while 1:
            if safe_peekline(uhandle)[:2] != 'RN':
                break
            self._scan_rn(uhandle, consumer)
            self._scan_rp(uhandle, consumer)
            self._scan_rc(uhandle, consumer)
            self._scan_rx(uhandle, consumer)
            self._scan_ra(uhandle, consumer)
            self._scan_rt(uhandle, consumer)
            self._scan_rl(uhandle, consumer)
    
    def _scan_rn(self, uhandle, consumer):
        self._scan_line('RN', uhandle, consumer.reference_number,
                        exactly_one=1)
    
    def _scan_rp(self, uhandle, consumer):
        self._scan_line('RP', uhandle, consumer.reference_position,
                        exactly_one=1)
    
    def _scan_rc(self, uhandle, consumer):
        self._scan_line('RC', uhandle, consumer.reference_comment,
                        any_number=1)
    
    def _scan_rx(self, uhandle, consumer):
        self._scan_line('RX', uhandle, consumer.reference_cross_reference,
                        up_to_one=1)
    
    def _scan_ra(self, uhandle, consumer):
        self._scan_line('RA', uhandle, consumer.reference_author,
                        one_or_more=1)
    
    def _scan_rt(self, uhandle, consumer):
        self._scan_line('RT', uhandle, consumer.reference_title,
                        any_number=1)
    
    def _scan_rl(self, uhandle, consumer):
        self._scan_line('RL', uhandle, consumer.reference_location,
                        one_or_more=1)
    
    def _scan_cc(self, uhandle, consumer):
        self._scan_line('CC', uhandle, consumer.comment, any_number=1)
    
    def _scan_dr(self, uhandle, consumer):
        self._scan_line('DR', uhandle, consumer.database_cross_reference,
                        any_number=1)
    
    def _scan_kw(self, uhandle, consumer):
        self._scan_line('KW', uhandle, consumer.keyword, any_number=1)
    
    def _scan_ft(self, uhandle, consumer):
        self._scan_line('FT', uhandle, consumer.feature_table, any_number=1)
    
    def _scan_sq(self, uhandle, consumer):
        self._scan_line('SQ', uhandle, consumer.sequence_header, exactly_one=1)
    
    def _scan_sequence_data(self, uhandle, consumer):
        self._scan_line('  ', uhandle, consumer.sequence_data, one_or_more=1)
    
    def _scan_terminator(self, uhandle, consumer):
        self._scan_line('//', uhandle, consumer.terminator, exactly_one=1)
    
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
