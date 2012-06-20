# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO parser for HMMER domain table output format."""

from itertools import chain

from Bio._py3k import _as_bytes, _bytes_to_string
from Bio.SearchIO._objects import QueryResult, Hit, HSP
from Bio.SearchIO._index import SearchIndexer


def hmmer_domtab_hmmhit_iterator(handle):
    """Generator function to parse HMMER domain table output as QueryResults.

    handle -- Handle to the file.

    This iterator assumes that the HMM coordinates are search hit coordinates.

    """
    for qresult in HmmerDomtabIterator(handle, hmm_as_hit=True):
        yield qresult


def hmmer_domtab_hmmquery_iterator(handle):
    """Generator function to parse HMMER domain table output as QueryResults.

    handle -- Handle to the file.

    This iterator assumes that the HMM coordinates are search query
    coordinates.

    """
    for qresult in HmmerDomtabIterator(handle, hmm_as_query=False):
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


class HmmerDomtabHmmhitIterator(object):

    """Parser for the HMMER domain table format that assumes HMM profile
    coordinaes are hit coordinates."""

    def __init__(self, handle):
        self.handle = handle
        self.line = read_forward(self.handle)

    def __iter__(self):
        # stop iterating if it's an empty file
        if not self.line:
            raise StopIteration
        # if line starts with '#', it's a header line
        # and we want to read through that
        else:
            if self.line.startswith('#'):
                while True:
                    self.line = read_forward(self.handle)
                    # break out of loop when it's not the header lines anymore
                    if not self.line.startswith('#'):
                        break
            # stop iterating if we only have headers
            if not self.line:
                raise StopIteration
            for qresult in self.parse_qresult():
                yield qresult

    def parse_result_row(self):
        """Returns a dictionary of parsed row values."""

        return {'qresult': qresult, 'hit': hit, 'hsp': hsp}

    def parse_qresult(self):
        """Generator function that returns QueryResult objects."""
        qid_cache = ''


class HmmerDomtabHmmhitIndexer(SearchIndexer):

    """Indexer class for HMMER domain table output that assumes HMM profile
    coordinates are hit coordinates."""

    def __init__(self, *args, **kwargs):
        SearchIndexer.__init__(self, *args, **kwargs)
        # set parser for on-the-fly parsing
        self._parser = hmmer_domtab_iterator
        self._handle.seek(0)

    def __iter__(self):
        """Iterates over the file handle; yields key, start offset, and length."""
        handle = self._handle
        handle.seek(0)
        split_char = _as_bytes(' ')
        qresult_key = None

        # read through header
        while True:
            start_offset = handle.tell()
            line = read_forward(handle, strip=False)
            if not line.startswith('#'):
                break

        # and index the qresults
        #while True:
        #    end_offset = handle.tell()

        #    if not line:
        #        break
        #    if qresult_key is None:

    def get_raw(self, offset):
        """Returns the raw string of a QueryResult object from the given offset."""
        handle = self._handle
        handle.seek(offset)
        split_char = _as_bytes(' ')
        qresult_key = None
        qresult_raw = ''

        return qresult_raw


class HmmerDomtabHmmhitWriter(object):

    """Writer for hmmer-domtab output format which writes hit coordinates
    as HMM profile coordinates."""

    def __init__(self, handle):
        self.handle = handle

    def write_file(self, qresults):
        """Writes to the handle.

        Returns a tuple of how many QueryResult, Hit, and HSP objects were written.

        """
        handle = self.handle
        qresult_counter, hit_counter, hsp_counter = 0, 0, 0

        try:
            first_qresult = qresults.next()
        except StopIteration:
            handle.write(self.build_header())
        else:
            # write header
            handle.write(self.build_header(first_qresult))
            # and then the qresults
            for qresult in chain([first_qresult], qresults):
                if qresult:
                    handle.write(self.build_row(qresult))
                    qresult_counter += 1
                    hit_counter += len(qresult)
                    hsp_counter += sum([len(hit) for hit in qresult])

        return qresult_counter, hit_counter, hsp_counter

    def build_header(self, first_qresult=None):
        """Returns the header string of a HMMER table output."""

        return header

    def build_row(self, qresult):
        """Returns a string or one row or more of the QueryResult object."""

        return rows


def _test():
    """Run the Bio.SearchIO.HmmerIO.hmmerdomtab module's doctests.

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
