# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.SearchIO base classes for HMMER-related code."""

from Bio.SearchIO._index import SearchIndexer


class _BaseHmmerTextIndexer(SearchIndexer):
    """Base indexer class for HMMER plain text output."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._preamble = b""

    def get_raw(self, offset):
        """Return the raw record from the file as a bytes string."""
        handle = self._handle
        qresult_raw = b""

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
