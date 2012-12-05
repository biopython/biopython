# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

#TODO: factor out this module with SeqIO's _index to stay DRY

"""Custom indexing for Bio.SearchIO objects."""

from StringIO import StringIO
from Bio._py3k import _bytes_to_string
from Bio import bgzf
from Bio.File import _IndexedSeqFileProxy


class SearchIndexer(_IndexedSeqFileProxy):
    """Base class for file format specific random access.

    Subclasses for each file format should define '_parser' and optionally
    'get_raw' methods.
    """

    def __init__(self, filename, **kwargs):
        h = open(filename, 'rb')
        try:
            self._handle = bgzf.BgzfReader(mode="rb", fileobj=h)
        except ValueError, e:
            assert "BGZF" in str(e)
            #Not a BGZF file
            h.seek(0)
            self._handle = h
        self._kwargs = kwargs

    def _parse(self, handle):
        return iter(self._parser(handle, **self._kwargs)).next()

    def get(self, offset):
        return self._parse(StringIO(_bytes_to_string(self.get_raw(offset))))
