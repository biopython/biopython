# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO base classes for HMMER-related code."""

from Bio._py3k import _as_bytes
from Bio.SearchIO._index import SearchIndexer


__docformat__ = "restructuredtext en"


class _BaseHmmerTextIndexer(SearchIndexer):

    """Base indexer class for HMMER plain text output."""

    def __init__(self, *args, **kwargs):
        super(_BaseHmmerTextIndexer, self).__init__(*args, **kwargs)
        self._preamble = _as_bytes('')

    def get_raw(self, offset):
        handle = self._handle
        qresult_raw = _as_bytes('')

        # read header first
        if not self._preamble:
            handle.seek(0)
            while True:
                line = handle.readline()
                if line.startswith(self.qresult_start):
                    break
                qresult_raw += line
        else:
            qresult_raw += self._preamble

        # and read the qresult raw string
        handle.seek(offset)
        while True:
            # preserve whitespace, don't use read_forward
            line = handle.readline()
            qresult_raw += line

            # break when we've reached qresult end
            if line.startswith(self.qresult_end) or not line:
                break

        return qresult_raw
