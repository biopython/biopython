# Copyright 2003 Iddo Friedberg. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import string
"""A parser for the NCBI blastpgp version 2.2.5 output format. Currently only supports
the '-m 9' option, (table w/ annotations).
Returns a BlastTableRec instance
"""

class BlastTableEntry:
   def __init__(self,in_rec):
      bt_fields = in_rec.split()
      self.qid = bt_fields[0].split('|')
      self.sid = bt_fields[1].split('|')
      self.pid = string.atof(bt_fields[2])
      self.ali_len = string.atoi(bt_fields[3])
      self.mis = string.atoi(bt_fields[4])
      self.gaps = string.atoi(bt_fields[5])
      self.q_bounds = (string.atoi(bt_fields[6]), string.atoi(bt_fields[7]))
      self.s_bounds = (string.atoi(bt_fields[8]), string.atoi(bt_fields[9]))
      self.e_value = string.atof(bt_fields[10])
      self.bit_score = string.atof(bt_fields[11])
      
class BlastTableRec:
   def __init__(self):
      self.program = None
      self.version = None
      self.date = None
      self.iteration = None
      self.query = None
      self.database = None
      self.entries = []
   def add_entry(self, entry):
      self.entries.append(entry)

reader_keywords = {'BLASTP': 'version',
                   'Iteration': 'iteration',
                   'Query': 'query',
                   'Database': 'database',
                   'Fields': 'fields'}
class BlastTableReader:
   def __init__(self, handle):
      self.handle = handle
      inline = self.handle.readline()
      # zip forward to start of record
      while inline and inline.find('BLASTP') == -1:
         inline = self.handle.readline()
      self._lookahead = inline
      self._n = 0
      self._in_header = 1
   def next(self):
      self.table_record = BlastTableRec()
      self._n += 1
      inline = self._lookahead
      if not inline:
         return None
      while inline:
         if inline[0] == '#':
            if self._in_header:
               self._in_header = self._consume_header(inline)
            else:
               break
         else:
            self._consume_entry(inline)
            self._in_header = 0

         inline = self.handle.readline()
      self._lookahead = inline
      self._in_header = 1
      return self.table_record
         
   def _consume_entry(self, inline):
      current_entry = BlastTableEntry(inline)
      self.table_record.add_entry(current_entry)
   def _consume_header(self, inline):
      for keyword in reader_keywords.keys():
         if inline.find(keyword) > -1:
            in_header = self._Parse('_parse_%s' % reader_keywords[keyword],inline)
            break
      return in_header
   def _parse_version(self, inline):
      program, version, date = inline.split()[1:]
      self.table_record.program = program
      self.table_record.version = version
      self.table_record.date = date
      return 1
   def _parse_iteration(self, inline):
      self.table_record.iteration = string.atoi(inline.split()[2])
      return 1
   def _parse_query(self, inline):
      self.table_record.query = inline.split()[2:]
      return 1
   def _parse_database(self, inline):
      self.table_record.database = inline.split()[2]
      return 1
   def _parse_fields(self, inline):
      return 0
   def _Parse(self, method_name, inline):
      return getattr(self,method_name)(inline)
      


