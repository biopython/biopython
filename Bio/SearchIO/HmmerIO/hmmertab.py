# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO parser for HMMER table output format."""

from itertools import chain

from Bio._py3k import _as_bytes, _bytes_to_string
from Bio.SearchIO._objects import QueryResult, Hit, HSP
from Bio.SearchIO._index import SearchIndexer


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
        # read through the header if it exists
        while self.line.startswith('#'):
            self.line = read_forward(self.handle)
        # if we have result rows, parse it
        if self.line:
            for qresult in self.parse_qresult():
                yield qresult

    def parse_result_row(self):
        """Returns a dictionary of parsed row values."""
        cols = filter(None, self.line.strip().split(' '))
        # if len(cols) > 19, we have extra description columns
        # combine them all into one string in the 19th column
        if len(cols) > 19:
            cols[18] = ' '.join(cols[18:])
        # if it's < 19, we have no description columns, so use an empty string
        # instead
        elif len(cols) < 19:
            cols.append('')
            assert len(cols) == 19

        # assign parsed column data into qresult, hit, and hsp dicts
        qresult = {}
        qresult['id'] = cols[2]                     # query name
        qresult['acc'] = cols[3]                    # query accession
        hit = {}
        hit['id'] = cols[0]                         # target name
        hit['acc'] = cols[1]                        # target accession
        hit['evalue'] = float(cols[4])              # evalue (full sequence)
        hit['bitscore'] = float(cols[5])            # score (full sequence)
        hit['bias'] = float(cols[6])                # bias (full sequence)
        hit['domain_exp_num'] = float(cols[10])     # exp
        hit['region_num'] = int(cols[11])           # reg
        hit['cluster_num'] = int(cols[12])          # clu
        hit['overlap_num'] = int(cols[13])          # ov
        hit['env_num'] = int(cols[14])              # env
        hit['domain_obs_num'] = int(cols[15])       # dom
        hit['domain_reported_num'] = int(cols[16])  # rep
        hit['domain_included_num'] = int(cols[17])  # inc
        hit['desc'] = cols[18]                      # description of target
        hsp = {}
        hsp['evalue'] = float(cols[7])              # evalue (best 1 domain)
        hsp['bitscore'] = float(cols[8])            # score (best 1 domain)
        hsp['bias'] = float(cols[9])                # bias (best 1 domain)
        # strand is always 0, since HMMER now only handles protein
        hsp['hit_strand'] = hsp['query_strand'] = 0

        return qresult, hit, hsp

    def parse_qresult(self):
        """Generator function that returns QueryResult objects."""
        qid_cache = ''
        while True:
            # only parse the result row if it's not EOF
            if self.line:
                qres_attrs, hit_attrs, hsp_attrs = self.parse_result_row()
                qresult_id = qres_attrs['id']

            # a new qresult is created whenever qid_cache != qresult_id
            if qid_cache != qresult_id:
                # yield qresult if qid_cache is filled
                if qid_cache:
                    yield qresult
                qid_cache = qresult_id
                qresult = QueryResult(qresult_id)
                for attr, value in qres_attrs.items():
                    setattr(qresult, attr, value)
            # when we've reached EOF, try yield any remaining qresult and break
            elif not self.line:
                yield qresult
                break

            # create Hit and set its attributes
            hit_id = hit_attrs['id']
            hit = Hit(hit_id, qresult_id)
            for attr, value in hit_attrs.items():
                setattr(hit, attr, value)

            # create HSP and set its attributes
            hsp = HSP(hit_id, qresult_id)
            for attr, value in hsp_attrs.items():
                setattr(hsp, attr, value)

            # since domain tab formats only have 1 HSP per line
            # we don't have to worry about appending other HSPs to the Hit
            hit.append(hsp)
            qresult.append(hit)

            self.line = read_forward(self.handle)


class HmmerTabIndexer(SearchIndexer):

    """Indexer class for HMMER table output."""

    _parser = HmmerTabIterator
    # denotes column location for query identifier
    _query_id_idx = 2

    def __iter__(self):
        """Iterates over the file handle; yields key, start offset, and length."""
        handle = self._handle
        handle.seek(0)
        split_char = _as_bytes(' ')
        query_id_idx = self._query_id_idx
        qresult_key = None
        # set line with initial mock value, to emulate header
        line = '#'

        # read through header
        while line.startswith('#'):
            start_offset = handle.tell()
            line = read_forward(handle, strip=False)

        # and index the qresults
        while True:
            end_offset = handle.tell()

            if not line:
                break
            if qresult_key is None:
                qresult_key = filter(None, \
                        line.strip().split(split_char))[query_id_idx]
            else:
                curr_key = filter(None, \
                        line.strip().split(split_char))[query_id_idx]

                if curr_key != qresult_key:
                    yield _bytes_to_string(qresult_key), start_offset, \
                            end_offset - start_offset
                    qresult_key = curr_key
                    start_offset = end_offset - len(line)

            line = read_forward(handle, strip=False)
            if not line:
                yield _bytes_to_string(qresult_key), start_offset, \
                        end_offset - start_offset
                break

    def get_raw(self, offset):
        """Returns the raw string of a QueryResult object from the given offset."""
        handle = self._handle
        handle.seek(offset)
        split_char = _as_bytes(' ')
        query_id_idx = self._query_id_idx
        qresult_key = None
        qresult_raw = ''

        while True:
            line = read_forward(handle, strip=False)
            if not line:
                break
            if qresult_key is None:
                qresult_key = filter(None, \
                        line.strip().split(split_char))[query_id_idx]
            else:
                curr_key = filter(None, \
                        line.strip().split(split_char))[query_id_idx]
                if curr_key != qresult_key:
                    break
            qresult_raw += line

        return qresult_raw


class HmmerTabWriter(object):

    """Writer for hmmer-tab output format."""

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

        # calculate whitespace required
        # adapted from HMMER's source: src/p7_tophits.c#L1083
        if first_qresult is not None:
            #qnamew = max(20, len(first_qresult.id))
            qnamew = 20 # why doesn't the above work?
            tnamew = max(20, len(first_qresult[0].id))
            qaccw = max(10, len(first_qresult.acc))
            taccw = max(10, len(first_qresult[0].acc))
        else:
            qnamew, tnamew, qaccw, taccw = 20, 20, 10, 10

        header = ''
        header += "#%*s %22s %22s %33s\n" % \
                (tnamew + qnamew + taccw + qaccw + 2, "", \
                "--- full sequence ----", "--- best 1 domain ----", \
                "--- domain number estimation ----")
        header += "#%-*s %-*s %-*s %-*s %9s %6s %5s %9s %6s %5s %5s %3s " \
                "%3s %3s %3s %3s %3s %3s %s\n" % (tnamew-1, " target name", \
                        taccw, "accession",  qnamew, "query name", qaccw, \
                        "accession",  "  E-value", " score", " bias", \
                        "  E-value", " score", " bias", "exp", \
                        "reg", "clu", " ov", "env", "dom", "rep", \
                        "inc", "description of target")
        header += "#%*s %*s %*s %*s %9s %6s %5s %9s %6s %5s %5s %3s %3s " \
                "%3s %3s %3s %3s %3s %s\n" % (tnamew-1, "-------------------", \
                taccw, "----------", qnamew, "--------------------", qaccw, \
                "----------", "---------", "------", "-----", "---------", \
                "------", "-----", "---", "---", "---", "---", "---", "---", \
                "---", "---", "---------------------")

        return header

    def build_row(self, qresult):
        """Returns a string or one row or more of the QueryResult object."""
        rows = ''

        # calculate whitespace required
        # adapted from HMMER's source: src/p7_tophits.c#L1083
        qnamew = max(20, len(qresult.id))
        tnamew = max(20, len(qresult[0].id))
        qaccw = max(10, len(qresult.acc))
        taccw = max(10, len(qresult[0].acc))

        for hit in qresult:
            rows += "%-*s %-*s %-*s %-*s %9.2g %6.1f %5.1f %9.2g %6.1f %5.1f " \
            "%5.1f %3d %3d %3d %3d %3d %3d %3d %s\n" % (tnamew, hit.id, taccw, \
            hit.acc, qnamew, qresult.id, qaccw, qresult.acc, hit.evalue, \
            hit.bitscore, hit.bias, hit[0].evalue, hit[0].bitscore, \
            hit[0].bias, hit.domain_exp_num, hit.region_num, hit.cluster_num, \
            hit.overlap_num, hit.env_num, hit.domain_obs_num, \
            hit.domain_reported_num, hit.domain_included_num, hit.desc)

        return rows


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
