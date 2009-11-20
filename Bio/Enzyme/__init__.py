# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Module to work with enzyme.dat file (DEPRECATED).

This module provides code to work with the enzyme.dat file from
Enzyme (OBSOLETE as of Biopython version 1.50).
http://www.expasy.ch/enzyme/

The functionality of Bio.Enzyme has moved to Bio.ExPASy.ExPASy;
please use that module instead of Bio.Enzyme. Bio.Enzyme is now
deprecated and will be removed in a future release of Biopython.
"""

import warnings
warnings.warn("Bio.Enzyme is deprecated, and will be removed in a"\
              " future release of Biopython. Most of the functionality "
              " is now provided by Bio.ExPASy.Enzyme.  If you want to "
              " continue to use Bio.Enzyme, please get in contact "
              " via the mailing lists to avoid its permanent removal from"\
              " Biopython.", DeprecationWarning)

from Bio import File
from Bio.ParserSupport import *

class _Scanner:
    """Scans Enzyme data (PRIVATE).

    Tested with:
    Release 33
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
class DataRecord:
    def __init__(self,tr_code='',sw_code=''):
        self.tr_code = tr_code
        self.sw_code = sw_code
    
    def __str__(self):
        return self.tr_code + ", " + self.sw_code
        
class EnzymeRecord:
    def __init__(self):
        self.ID = ''
        self.DE = []
        self.AN = []
        self.CA = ''
        self.CF = []
        self.CC = []   # one comment per line
        self.DI = []
        self.PR = []
        self.DR = []
    
    def __repr__(self):
        if self.ID:
            if self.DE:
                return "%s (%s, %s)" % (self.__class__.__name__, 
                                        self.ID, self.DE[0])
            else:
                return "%s (%s)" % (self.__class__.__name__, 
                                       self.ID)
        else:
            return "%s ( )" % (self.__class__.__name__)
            
    def __str__(self):
        output = "ID: " + self.ID
        output += " DE: " + repr(self.DE)
        output += " AN: " + repr(self.AN)
        output += " CA: '" + self.CA + "'"
        output += " CF: " + repr(self.CF)
        output += " CC: " + repr(self.CC)
        output += " DI: " + repr(self.DI)
        output += " PR: " + repr(self.PR)
        output += " DR: %d Records" % len(self.DR)
        
        return output
        
class RecordParser(AbstractParser):
	def __init__(self):
		self._scanner = _Scanner()
		self._consumer = _RecordConsumer()

	def parse(self, handle):
		if isinstance(handle, File.UndoHandle):
			uhandle = handle
		else:
			uhandle = File.UndoHandle(handle)
			self._scanner.feed(uhandle, self._consumer)
		return self._consumer.enzyme_record

class Iterator:
	def __init__(self, handle, parser=None):
		self._uhandle = File.UndoHandle(handle)

	def next(self):
		self._parser = RecordParser()
		lines = []
		while 1:
			line = self._uhandle.readline()
			if not line: break
			if line[:2] == '//':
				break
			lines.append(line)
		if not lines:
			return None
		lines.append('//')
		data = string.join(lines,'')
		if self._parser is not None:
			return self._parser.parse(File.StringHandle(data))
		return data

        def __iter__(self):
                return iter(self.next, None)

class _RecordConsumer(AbstractConsumer):
    def __init__(self):
        self.enzyme_record = EnzymeRecord()
    def identification(self, id_info):
        self.enzyme_record.ID = id_info.split()[1]
    def description(self,de_info):
        self.enzyme_record.DE.append(de_info[2:].strip())
    def alternate_name(self,an_info):
        self.enzyme_record.AN.append(an_info[2:].strip())
    def catalytic_activity(self, ca_info):
        self.enzyme_record.CA = string.join([self.enzyme_record.CA,ca_info[2:].strip()],'')
    def cofactor(self, cf_info):
        self.enzyme_record.CF.append(cf_info[2:].strip())
    def comment(self, cc_info):
        cc = cc_info[2:].strip()
        if cc.startswith("-!-"):
            self.enzyme_record.CC.append(cc[len("-!-"):].strip())
        else:
            # The header is all CC, but doesn't start with -!-
            if self.enzyme_record.CC:
                pre_cc = self.enzyme_record.CC.pop()
            else:
                pre_cc = ""
            new_cc = pre_cc + " " + cc
            self.enzyme_record.CC.append(new_cc)
    def disease(self, di_info):
        self.enzyme_record.DI.append(di_info[2:].strip())
        
    def prosite_reference(self,pr_info):
        self.enzyme_record.PR.append(pr_info.split(';')[1].strip())

    def databank_reference(self,dr_info):
        good_data = dr_info[2:].strip()
        pair_data = good_data.split(';')
        for pair in pair_data:
            if not pair: continue
            data_record = DataRecord()
            t1, t2 = pair.split(',')
            data_record.tr_code, data_record.sw_code = \
                t1.strip(), t2.strip()
            self.enzyme_record.DR.append(data_record)

	def terminator(self,schwarzenegger):
		pass # Hasta la Vista, baby!
