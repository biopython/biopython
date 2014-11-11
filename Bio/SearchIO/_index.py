# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Custom indexing for Bio.SearchIO objects."""

from Bio._py3k import StringIO
from Bio._py3k import _bytes_to_string
from Bio import bgzf
from Bio.File import _IndexedSeqFileProxy, _open_for_random_access


__docformat__ = "restructuredtext en"


class SearchIndexer(_IndexedSeqFileProxy):
    """Base class for file format specific random access.

    Subclasses for each file format should define '_parser' and optionally
    'get_raw' methods.
    """

    def __init__(self, filename, **kwargs):
        self._handle = _open_for_random_access(filename)
        self._kwargs = kwargs

    def _parse(self, handle):
        return next(iter(self._parser(handle, **self._kwargs)))

    def get(self, offset):
        return self._parse(StringIO(_bytes_to_string(self.get_raw(offset))))
