# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO parser for HMMER table output format."""

from Bio._py3k import _as_bytes, _bytes_to_string
from Bio.SearchIO._objects import QueryResult, Hit, HSP
from Bio.SearchIO._index import SearchIndexer


def hmmer_tab_iterator(handle):
    """Generator function to parse HMMER table output as QueryResult objects.

    handle -- Handle to the file.

    """
    for qresult in HmmerTabIterator(handle):
        yield qresult


def read_forward(handle, strip=True):
    """Reads through whitespaces, returns the first non-whitespace line."""
    while True:
        line = handle.readline()
        # return the line if it has characters
        if line and line.strip():
            if strip:
                return line.strip()
            else:
                return line
        # or if has no characters (EOF)
        elif not line:
            return line


class HmmerTabIterator(object):

    """Parser for the HMMER table format."""

    def __init__(self, handle):
        self.handle = handle
        self.line = read_forward(self.handle)

    def __iter__(self):
        pass

    def parse_result_row(self):
        """Returns a dictionary of parsed row values."""

    def parse_qresult(self):
        """Generator function that returns QueryResult objects."""


class HmmerTabIndexer(SearchIndexer):

    """Indexer class for HMMER table output."""

    def __init__(self, *args, **kwargs):
        SearchIndexer.__init__(self, *args, **kwargs)
        # set parser for on-the-fly parsing
        self._parser = hmmer_tab_iterator
        self._handle.seek(0)

    def __iter__(self):
        """Iterates over the file handle; yields key, start offset, and length."""
        handle = self._handle
        handle.seek(0)
        start_offset = handle.tell()

    def get_raw(self, offset):
        """Returns the raw string of a QueryResult object from the given offset."""
        handle = self._handle
        handle.seek(offset)
        qresult_raw = ''

        return qresult_raw


def _test():
    """Run the Bio.SearchIO.HmmerIO.hmmertab module's doctests.

    This will try and locate the unit tests directory, and run the doctests
    from there in order that the relative paths used in the examples work.
    """
    import doctest
    import os

    test_dir = 'Tests'

    if os.path.isdir(os.path.join('..', '..', test_dir)):
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join('..', '..', test_dir))
        doctest.testmod()
        os.chdir(cur_dir)
        print "Done"


# if not used as a module, run the doctest
if __name__ == "__main__":
    _test()
