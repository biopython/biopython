# Copyright 2001 by Katharine Lindner.  All rights reserved.
# Copyright 2006 by PeterC.  All rights reserved.
# Copyright 2007 by Michiel de Hoon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Parser for files from NCBI's Gene Expression Omnibus (GEO).

http://www.ncbi.nlm.nih.gov/geo/
"""

import Record


def _read_key_value(line):
    words = line[1:].split("=", 1)
    try:
        key, value = words
        value = value.strip()
    except ValueError:
        key = words[0]
        value = ""
    key = key.strip()
    return key, value


def parse(handle):
    record = None
    for line in handle:
        line = line.strip('\n').strip('\r')
        if not line: continue # Ignore empty lines
        c = line[0]
        if c=='^':
            if record: yield record
            record = Record.Record()
            record.entity_type, record.entity_id = _read_key_value(line)
        elif c=='!':
            if line in ('!Sample_table_begin',
                        '!Sample_table_end',
                        '!Platform_table_begin',
                        '!Platform_table_end'):
                continue
            key, value = _read_key_value(line)
            if key in record.entity_attributes:
                if type(record.entity_attributes[key])==list:
                    record.entity_attributes[key].append(value)
                else:
                    existing = record.entity_attributes[key]
                    record.entity_attributes[key] = [existing, value]
            else:
                record.entity_attributes[key] = value
        elif c=='#':
            key, value = _read_key_value(line)
            assert key not in record.col_defs
            record.col_defs[key] = value
        else:
            row = line.split("\t")
            record.table_rows.append(row)
    yield record


class Iterator:
    """Iterator interface to move over a file of Geo entries one at a time.

    Uses the fact that each GEO record begins with a line starting ^ (caret).
    """
    def __init__(self, handle, parser = None):
        """Initialize the iterator.

        Arguments:
        o handle - A handle with GEO entries to iterate through.
        returning them. If None, then the raw entry will be returned.
        """
        import warnings
        warnings.warn("Bio.Geo.Iterator(handle, parser) is deprecated. Please use Bio.Geo.parse(handle) instead. It also returns an iterator.""",
              DeprecationWarning)
        self.records = parse(handle)

    def next(self):
        return self.records.next()

    def __iter__(self):
        return iter(self.next, None)
