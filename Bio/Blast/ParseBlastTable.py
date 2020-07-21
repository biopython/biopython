# Copyright 2003 Iddo Friedberg. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""A parser for the NCBI blastpgp version 2.2.5 output format.

Currently only supports the '-m 9' option, (table w/ annotations).
Returns a BlastTableRec instance
"""


class BlastTableEntry:
    """Container for Blast Table Entry, the field values from the table."""

    def __init__(self, in_rec):
        """Initialize the class."""
        bt_fields = in_rec.split()
        self.qid = bt_fields[0].split("|")
        self.sid = bt_fields[1].split("|")
        self.pid = float(bt_fields[2])
        self.ali_len = int(bt_fields[3])
        self.mis = int(bt_fields[4])
        self.gaps = int(bt_fields[5])
        self.q_bounds = (int(bt_fields[6]), int(bt_fields[7]))
        self.s_bounds = (int(bt_fields[8]), int(bt_fields[9]))
        self.e_value = float(bt_fields[10])
        self.bit_score = float(bt_fields[11])


class BlastTableRec:
    """Container for Blast Table record, list of Blast Table Entries."""

    def __init__(self):
        """Initialize the class."""
        self.program = None
        self.version = None
        self.date = None
        self.iteration = None
        self.query = None
        self.database = None
        self.entries = []

    def add_entry(self, entry):
        """Add entry to Blast Table."""
        self.entries.append(entry)


class BlastTableReader:
    """Reader for the output of blastpgp."""

    reader_keywords = {
        "BLASTP": "version",
        "Iteration": "iteration",
        "Query": "query",
        "Database": "database",
        "Fields": "fields",
    }

    def __init__(self, handle):
        """Initialize the class."""
        self.handle = handle
        inline = self.handle.readline()
        # zip forward to start of record
        while inline and "BLASTP" not in inline:
            inline = self.handle.readline()
        self._lookahead = inline
        self._n = 0
        self._in_header = 1

    def __next__(self):
        """Return the next record when iterating over the file."""
        self.table_record = BlastTableRec()
        self._n += 1
        inline = self._lookahead
        if not inline:
            return None
        while inline:
            if inline[0] == "#":
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
        for keyword in self.reader_keywords:
            if keyword in inline:
                return self._Parse("_parse_%s" % self.reader_keywords[keyword], inline)

    def _parse_version(self, inline):
        program, version, date = inline.split()[1:]
        self.table_record.program = program
        self.table_record.version = version
        self.table_record.date = date
        return 1

    def _parse_iteration(self, inline):
        self.table_record.iteration = int(inline.split()[2])
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
        return getattr(self, method_name)(inline)
