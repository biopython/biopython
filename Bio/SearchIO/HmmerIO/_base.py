# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO base classes for HMMER-related code."""

import re

from Bio._py3k import _as_bytes, _bytes_to_string
from Bio.SearchIO._index import SearchIndexer
from Bio.SearchIO._utils import read_forward


class _BaseHmmerTextIndexer(SearchIndexer):

    """Base indexer class for HMMER plain text output."""

    def __iter__(self):
        handle = self._handle
        handle.seek(0)
        start_offset = handle.tell()

        while True:
            line = read_forward(handle)
            end_offset = handle.tell()

            if line.startswith(self.qresult_start):
                regx = re.search(self.regex_id, line)
                qresult_key = regx.group(1).strip()
                # qresult start offset is the offset of this line
                # (starts with the start mark)
                start_offset = end_offset - len(line)
            elif line.startswith(self.qresult_end):
                yield _bytes_to_string(qresult_key), start_offset, 0
                start_offset = end_offset
            elif not line:
                break

    def get_raw(self, offset):
        handle = self._handle
        qresult_raw = _as_bytes('')

        # read header first
        handle.seek(0)
        while True:
            line = handle.readline()
            if line.startswith(self.qresult_start):
                break
            qresult_raw += line

        # and read the qresult raw string
        handle.seek(offset)
        while True:
            # preserve whitespace, don't use read_forward
            line = handle.readline()
            qresult_raw += line

            # break when we've reached qresult end
            if line.startswith(self.qresult_end):
                break

        return qresult_raw
