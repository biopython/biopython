# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# Revisions copyright 2012-2016 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Custom indexing for Bio.SearchIO objects."""

from io import StringIO

from Bio.File import _IndexedSeqFileProxy, _open_for_random_access


class SearchIndexer(_IndexedSeqFileProxy):
    """Base class for file format specific random access.

    Subclasses for each file format should define '_parser' and optionally
    'get_raw' methods.
    """

    def __init__(self, filename, **kwargs):
        """Initialize the class."""
        self._handle = _open_for_random_access(filename)
        self._kwargs = kwargs

    def _parse(self, handle):
        """Pass handle and arguments to the next iterable (PRIVATE)."""
        return next(iter(self._parser(handle, **self._kwargs)))

    def get(self, offset):
        """Get offset and convert it from bytes to string."""
        return self._parse(StringIO(self.get_raw(offset).decode()))
