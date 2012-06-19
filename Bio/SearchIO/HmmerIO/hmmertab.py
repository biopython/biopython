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
        # stop iterating if it's an empty file
        if not self.line:
            raise StopIteration
        # if line starts with '#', it's a header line
        # and we want to read through that
        elif self.line.startswith('#'):
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
        assert self.line
        cols = filter(None, self.line.strip().split(' '))
        # if len(cols) > 18, we have extra description columns
        # combine them all into one string in the 19th column
        if len(cols) > 18:
            cols[18] = ' '.join(cols[18:])
        else:
            cols[18] = ''

        # assign parsed column data into qresult, hit, and hsp dicts
        qresult = {}
        qresult['id'] = cols[2]                  # query name
        qresult['acc'] = cols[3]                 # query accession
        hit = {}
        hit['id'] = cols[0]                      # target name
        hit['acc'] = cols[1]                     # target accession
        hit['evalue'] = cols[4]                  # evalue (full sequence)
        hit['bitscore'] = cols[5]                # score (full sequence)
        hit['bias'] = cols[6]                    # bias (full sequence)
        hit['domain_exp_num'] = cols[10]         # exp
        hit['reg_num'] = cols[11]                # reg
        hit['multidomain_num'] = cols[12]        # clu
        hit['env_overlap_num'] = cols[13]        # ov
        hit['env_num'] = cols[14]                # env
        hit['domain_obs_num'] = cols[15]         # dom
        hit['domain_report_num'] = cols[16]      # rep
        hit['domain_inclusion_num'] = cols[17]   # inc
        hit['desc'] = cols[18]                   # description of target
        hsp = {}
        hsp['evalue'] = cols[7]                  # evalue (best 1 domain)
        hsp['bitscore'] = cols[8]                # score (best 1 domain)
        hsp['bias'] = cols[9]                    # bias (best 1 domain)

        return {'qresult': qresult, 'hit': hit, 'hsp': hsp}

    def parse_qresult(self):
        """Generator function that returns QueryResult objects."""
        qid_cache = ''
        while True:
            # only parse the result row if it's not EOF
            if self.line:
                parsed = self.parse_result_row()
                qresult_id = parsed['qresult']['id']

            # a new qresult is created whenever qid_cache != qresult_id
            if qid_cache != qresult_id:
                # yield qresult if qid_cache is filled
                if qid_cache:
                    yield qresult
                qid_cache = qresult_id
                qresult = QueryResult(qresult_id)
                for attr, value in parsed['qresult'].items():
                    setattr(qresult, attr, value)
            # when we've reached EOF, try yield any remaining qresult and break
            elif not self.line:
                yield qresult
                break

            # create Hit and set its attributes
            hit_id = parsed['hit']['id']
            hit = Hit(hit_id, qresult_id)
            for attr, value in parsed['hit'].items():
                setattr(hit, attr, value)

            # create HSP and set its attributes
            hsp = HSP(hit_id, qresult_id)
            for attr, value in parsed['hsp'].items():
                setattr(hsp, attr, value)

            # since domain tab formats only have 1 HSP per line
            # we don't have to worry about appending other HSPs to the Hit
            hit.append(hsp)
            qresult.append(hit)

            self.line = read_forward(self.handle)


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
