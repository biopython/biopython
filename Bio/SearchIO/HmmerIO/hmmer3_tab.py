# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO parser for HMMER table output format."""

from itertools import chain

from Bio._py3k import _as_bytes, _bytes_to_string
from Bio.Alphabet import generic_protein
from Bio.SearchIO._index import SearchIndexer
from Bio.SearchIO._model import QueryResult, Hit, HSP, HSPFragment


__all__ = ['Hmmer3TabParser', 'Hmmer3TabIndexer', 'Hmmer3TabWriter']

__docformat__ = "restructuredtext en"


class Hmmer3TabParser(object):

    """Parser for the HMMER table format."""

    def __init__(self, handle):
        self.handle = handle
        self.line = self.handle.readline()

    def __iter__(self):
        header_mark = '#'
        # read through the header if it exists
        while self.line.startswith(header_mark):
            self.line = self.handle.readline()
        # if we have result rows, parse it
        if self.line:
            for qresult in self._parse_qresult():
                yield qresult

    def _parse_row(self):
        """Returns a dictionary of parsed row values."""
        cols = [x for x in self.line.strip().split(' ') if x]
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
        qresult['accession'] = cols[3]              # query accession
        hit = {}
        hit['id'] = cols[0]                         # target name
        hit['accession'] = cols[1]                  # target accession
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
        hit['description'] = cols[18]               # description of target
        hsp = {}
        hsp['evalue'] = float(cols[7])              # evalue (best 1 domain)
        hsp['bitscore'] = float(cols[8])            # score (best 1 domain)
        hsp['bias'] = float(cols[9])                # bias (best 1 domain)
        # strand is always 0, since HMMER now only handles protein
        frag = {}
        frag['hit_strand'] = frag['query_strand'] = 0
        frag['alphabet'] = generic_protein

        return {'qresult': qresult, 'hit': hit, 'hsp': hsp, 'frag': frag}

    def _parse_qresult(self):
        """Generator function that returns QueryResult objects."""
        # state values, determines what to do for each line
        state_EOF = 0
        state_QRES_NEW = 1
        state_QRES_SAME = 3
        # initial value dummies
        qres_state = None
        file_state = None
        prev_qid = None
        cur, prev = None, None
        # container for Hit objects, used to create QueryResult
        hit_list = []

        while True:
            # store previous line's parsed values for all lines after the first
            if cur is not None:
                prev = cur
                prev_qid = cur_qid
            # only parse the result row if it's not EOF
            # NOTE: we are not parsing the extra '#' lines appended to the end
            # of hmmer31b1 tabular results since storing them in qresult
            # objects means we can not do a single-pass parsing
            if self.line and not self.line.startswith('#'):
                cur = self._parse_row()
                cur_qid = cur['qresult']['id']
            else:
                file_state = state_EOF
                # mock value for cur_qid, since we have nothing to parse
                cur_qid = None

            if prev_qid != cur_qid:
                qres_state = state_QRES_NEW
            else:
                qres_state = state_QRES_SAME

            if prev is not None:
                # since domain tab formats only have 1 Hit per line
                # we always create HSPFragment, HSP, and Hit per line
                prev_hid = prev['hit']['id']

                # create fragment and HSP and set their attributes
                frag = HSPFragment(prev_hid, prev_qid)
                for attr, value in prev['frag'].items():
                    setattr(frag, attr, value)
                hsp = HSP([frag])
                for attr, value in prev['hsp'].items():
                    setattr(hsp, attr, value)

                # create Hit and set its attributes
                hit = Hit([hsp])
                for attr, value in prev['hit'].items():
                    setattr(hit, attr, value)
                hit_list.append(hit)

                # create qresult and yield if we're at a new qresult or at EOF
                if qres_state == state_QRES_NEW or file_state == state_EOF:
                    qresult = QueryResult(hit_list, prev_qid)
                    for attr, value in prev['qresult'].items():
                        setattr(qresult, attr, value)
                    yield qresult
                    # if we're at EOF, break
                    if file_state == state_EOF:
                        break
                    hit_list = []

            self.line = self.handle.readline()


class Hmmer3TabIndexer(SearchIndexer):

    """Indexer class for HMMER table output."""

    _parser = Hmmer3TabParser
    # denotes column location for query identifier
    _query_id_idx = 2

    def __iter__(self):
        """Iterates over the file handle; yields key, start offset, and length."""
        handle = self._handle
        handle.seek(0)
        query_id_idx = self._query_id_idx
        qresult_key = None
        header_mark = _as_bytes('#')
        split_mark = _as_bytes(' ')
        # set line with initial mock value, to emulate header
        line = header_mark

        # read through header
        while line.startswith(header_mark):
            start_offset = handle.tell()
            line = handle.readline()

        # and index the qresults
        while True:
            end_offset = handle.tell()

            if not line:
                break

            cols = [x for x in line.strip().split(split_mark) if x]
            if qresult_key is None:
                qresult_key = cols[query_id_idx]
            else:
                curr_key = cols[query_id_idx]

                if curr_key != qresult_key:
                    adj_end = end_offset - len(line)
                    yield _bytes_to_string(qresult_key), start_offset, \
                            adj_end - start_offset
                    qresult_key = curr_key
                    start_offset = adj_end

            line = handle.readline()
            if not line:
                yield _bytes_to_string(qresult_key), start_offset, \
                        end_offset - start_offset
                break

    def get_raw(self, offset):
        """Returns the raw string of a QueryResult object from the given offset."""
        handle = self._handle
        handle.seek(offset)
        query_id_idx = self._query_id_idx
        qresult_key = None
        qresult_raw = _as_bytes('')
        split_mark = _as_bytes(' ')

        while True:
            line = handle.readline()
            if not line:
                break
            cols = [x for x in line.strip().split(split_mark) if x]
            if qresult_key is None:
                qresult_key = cols[query_id_idx]
            else:
                curr_key = cols[query_id_idx]
                if curr_key != qresult_key:
                    break
            qresult_raw += line

        return qresult_raw


class Hmmer3TabWriter(object):

    """Writer for hmmer3-tab output format."""

    def __init__(self, handle):
        self.handle = handle

    def write_file(self, qresults):
        """Writes to the handle.

        Returns a tuple of how many QueryResult, Hit, and HSP objects were written.

        """
        handle = self.handle
        qresult_counter, hit_counter, hsp_counter, frag_counter = 0, 0, 0, 0

        try:
            first_qresult = next(qresults)
        except StopIteration:
            handle.write(self._build_header())
        else:
            # write header
            handle.write(self._build_header(first_qresult))
            # and then the qresults
            for qresult in chain([first_qresult], qresults):
                if qresult:
                    handle.write(self._build_row(qresult))
                    qresult_counter += 1
                    hit_counter += len(qresult)
                    hsp_counter += sum(len(hit) for hit in qresult)
                    frag_counter += sum(len(hit.fragments) for hit in qresult)

        return qresult_counter, hit_counter, hsp_counter, frag_counter

    def _build_header(self, first_qresult=None):
        """Returns the header string of a HMMER table output."""

        # calculate whitespace required
        # adapted from HMMER's source: src/p7_tophits.c#L1083
        if first_qresult is not None:
            # qnamew = max(20, len(first_qresult.id))
            qnamew = 20 # why doesn't the above work?
            tnamew = max(20, len(first_qresult[0].id))
            qaccw = max(10, len(first_qresult.accession))
            taccw = max(10, len(first_qresult[0].accession))
        else:
            qnamew, tnamew, qaccw, taccw = 20, 20, 10, 10

        header = "#%*s %22s %22s %33s\n" % \
                (tnamew + qnamew + taccw + qaccw + 2, "",
                "--- full sequence ----", "--- best 1 domain ----",
                "--- domain number estimation ----")
        header += "#%-*s %-*s %-*s %-*s %9s %6s %5s %9s %6s %5s %5s %3s " \
                "%3s %3s %3s %3s %3s %3s %s\n" % (tnamew-1, " target name",
                        taccw, "accession", qnamew, "query name", qaccw,
                        "accession", "  E-value", " score", " bias",
                        "  E-value", " score", " bias", "exp",
                        "reg", "clu", " ov", "env", "dom", "rep",
                        "inc", "description of target")
        header += "#%*s %*s %*s %*s %9s %6s %5s %9s %6s %5s %5s %3s %3s " \
                "%3s %3s %3s %3s %3s %s\n" % (tnamew-1, "-------------------",
                taccw, "----------", qnamew, "--------------------", qaccw,
                "----------", "---------", "------", "-----", "---------",
                "------", "-----", "---", "---", "---", "---", "---", "---",
                "---", "---", "---------------------")

        return header

    def _build_row(self, qresult):
        """Returns a string or one row or more of the QueryResult object."""
        rows = ''

        # calculate whitespace required
        # adapted from HMMER's source: src/p7_tophits.c#L1083
        qnamew = max(20, len(qresult.id))
        tnamew = max(20, len(qresult[0].id))
        qaccw = max(10, len(qresult.accession))
        taccw = max(10, len(qresult[0].accession))

        for hit in qresult:
            rows += "%-*s %-*s %-*s %-*s %9.2g %6.1f %5.1f %9.2g %6.1f %5.1f " \
            "%5.1f %3d %3d %3d %3d %3d %3d %3d %s\n" % (tnamew, hit.id, taccw,
            hit.accession, qnamew, qresult.id, qaccw, qresult.accession, hit.evalue,
            hit.bitscore, hit.bias, hit.hsps[0].evalue, hit.hsps[0].bitscore,
            hit.hsps[0].bias, hit.domain_exp_num, hit.region_num, hit.cluster_num,
            hit.overlap_num, hit.env_num, hit.domain_obs_num,
            hit.domain_reported_num, hit.domain_included_num, hit.description)

        return rows


# if not used as a module, run the doctest
if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()
