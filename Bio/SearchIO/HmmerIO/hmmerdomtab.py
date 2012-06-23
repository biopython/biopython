# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO parser for HMMER domain table output format."""

from itertools import chain

from Bio.SearchIO._objects import QueryResult, Hit, HSP
from hmmertab import HmmerTabIterator, HmmerTabIndexer


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
    for qresult in HmmerDomtabIterator(handle, hmm_as_hit=False):
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


class HmmerDomtabIterator(HmmerTabIterator):

    """Parser for the HMMER domain table format that assumes HMM profile
    coordinates are hit coordinates."""

    def __init__(self, handle, hmm_as_hit):
        HmmerTabIterator.__init__(self, handle)
        self.hmm_as_hit = hmm_as_hit

    def parse_result_row(self):
        """Returns a dictionary of parsed row values."""
        assert self.line
        cols = filter(None, self.line.strip().split(' '))
        # if len(cols) > 22, we have extra description columns
        # combine them all into one string in the 19th column
        if len(cols) > 22:
            cols[22] = ' '.join(cols[22:])
        else:
            cols[22] = ''

        # assign parsed column data into qresult, hit, and hsp dicts
        qresult = {}
        qresult['id'] = cols[3]                 # query name
        qresult['acc'] = cols[4]                # query accession
        qresult['seq_len'] = cols[5]            # qlen
        hit = {}
        hit['id'] = cols[0]                     # target name
        hit['acc'] = cols[1]                    # target accession
        hit['seq_len'] = cols[2]                # tlen
        hit['evalue'] = cols[6]                 # evalue
        hit['bitscore'] = cols[7]               # score
        hit['bias'] = cols[8]                   # bias
        hit['desc'] = cols[22]                  # description of target
        hsp = {}
        hsp['domain_index'] = cols[9]           # # (domain number)
        # not parsing cols[10] since it's basically len(hit)
        hsp['evalue_cond'] = cols[11]           # c-evalue
        hsp['evalue'] = cols[12]                # i-evalue
        hsp['bitscore'] = cols[13]              # score
        hsp['bias'] = cols[14]                  # bias
        hsp['hit_from'] = int(cols[15]) - 1     # hmm from
        hsp['hit_to'] = int(cols[16]) - 1       # hmm to
        hsp['query_from'] = int(cols[17]) - 1   # ali from
        hsp['query_to'] = int(cols[18]) - 1     # ali to
        hsp['env_from'] = int(cols[19]) - 1     # env from
        hsp['env_to'] = int(cols[20]) - 1       # env to
        hsp['acc_avg'] = cols[21]               # acc

        # switch hmm<-->ali coordinates if hmm is not hit
        if not self.hmm_as_hit:
            hsp['hit_to'], hsp['query_to'] = \
                    hsp['query_to'], hsp['hit_to']
            hsp['hit_from'], hsp['query_from'] = \
                    hsp['query_from'], hsp['hit_from']

        return qresult, hit, hsp

    def parse_qresult(self):
        """Generator function that returns QueryResult objects."""
        qid_cache = None
        hid_cache = None
        # flag for denoting whether the query is the same as the previous line
        # or not
        same_query = False
        while True:
            # only parse the result row if it's not EOF
            if self.line:
                qres_parsed, hit_parsed, hsp_parsed = self.parse_result_row()
                qresult_id = qres_parsed['id']

            # a new qresult is created whenever qid_cache != qresult_id
            if qid_cache != qresult_id:
                # append the last hit and yield qresult if qid_cache is filled
                if qid_cache is not None:
                    qresult.append(hit)
                    yield qresult
                    same_query = False
                qid_cache = qresult_id
                qresult = QueryResult(qresult_id)
                for attr, value in qres_parsed.items():
                    setattr(qresult, attr, value)
            # when we've reached EOF, try yield any remaining qresult and break
            elif not self.line:
                qresult.append(hit)
                yield qresult
                break
            # otherwise, we must still be in the same query, so set the flag
            elif not same_query:
                same_query = True

            hit_id = hit_parsed['id']
            # a new hit is created whenever hid_cache != hit_id
            if hid_cache != hit_id:
                # if we're in the same query, append the previous line's hit
                if same_query:
                    qresult.append(hit)
                hid_cache = hit_id
                hit = Hit(hit_id, qresult_id)
                for attr, value in hit_parsed.items():
                    setattr(hit, attr, value)

            # each line is basically a different HSP, so we always add it to
            # any hit object we have
            hsp = HSP(hit_id, qresult_id)
            for attr, value in hsp_parsed.items():
                setattr(hsp, attr, value)
            hit.append(hsp)

            self.line = read_forward(self.handle)


class HmmerDomtabHmmhitIndexer(HmmerTabIndexer):

    """Indexer class for HMMER domain table output that assumes HMM profile
    coordinates are hit coordinates."""

    def __init__(self, *args, **kwargs):
        HmmerTabIndexer.__init__(self, *args, **kwargs)
        # set parser for on-the-fly parsing
        self._parser = hmmer_domtab_hmmhit_iterator
        self._query_id_idx = 3


class HmmerDomtabHmmqueryIndexer(HmmerTabIndexer):

    """Indexer class for HMMER domain table output that assumes HMM profile
    coordinates are query coordinates."""

    def __init__(self, *args, **kwargs):
        HmmerTabIndexer.__init__(self, *args, **kwargs)
        # set parser for on-the-fly parsing
        self._parser = hmmer_domtab_hmmquery_iterator
        self._query_id_idx = 3


class HmmerDomtabHmmhitWriter(object):

    """Writer for hmmer-domtab output format which writes hit coordinates
    as HMM profile coordinates."""

    def __init__(self, handle):
        self.handle = handle
        self.hmm_as_hit = True

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
        """Returns the header string of a domain HMMER table output."""

        # calculate whitespace required
        # adapted from HMMER's source: src/p7_tophits.c#L1157
        if first_qresult:
            #qnamew = max(20, len(first_qresult.id))
            qnamew = 20
            tnamew = max(20, len(first_qresult[0].id))
            try:
                qaccw = max(10, len(first_qresult.acc))
                taccw = max(10, len(first_qresult[0].acc))
            except AttributeError:
                qaccw, taccw = 10, 10
        else:
            qnamew, tnamew, qaccw, taccw = 20, 20, 10, 10

        header = ''
        header += "#%*s %22s %40s %11s %11s %11s\n" % \
                (tnamew+qnamew-1+15+taccw+qaccw, "", "--- full sequence ---", \
                "-------------- this domain -------------", "hmm coord", \
                "ali coord", "env coord")
        header += "#%-*s %-*s %5s %-*s %-*s %5s %9s %6s %5s %3s %3s %9s " \
                "%9s %6s %5s %5s %5s %5s %5s %5s %5s %4s %s\n" % (tnamew-1, \
                " target name", taccw, "accession", "tlen", qnamew, \
                "query name", qaccw, "accession", "qlen", "E-value", "score", \
                "bias", "#", "of", "c-Evalue", "i-Evalue", "score", "bias", \
                "from", "to", "from", "to", "from", "to", "acc", \
                "description of target")
        header += "#%*s %*s %5s %*s %*s %5s %9s %6s %5s %3s %3s %9s %9s " \
                "%6s %5s %5s %5s %5s %5s %5s %5s %4s %s\n" % (tnamew-1, \
                "-------------------", taccw, "----------", "-----", \
                qnamew, "--------------------", qaccw, "----------", \
                "-----", "---------", "------", "-----", "---", "---", \
                "---------", "---------", "------", "-----", "-----", "-----", \
                "-----", "-----", "-----", "-----", "----", \
                "---------------------")

        return header

    def build_row(self, qresult):
        """Returns a string or one row or more of the QueryResult object."""
        rows = ''

        # calculate whitespace required
        # adapted from HMMER's source: src/p7_tophits.c#L1083
        qnamew = max(20, len(qresult.id))
        tnamew = max(20, len(qresult[0].id))
        try:
            qaccw = max(10, len(qresult.acc))
            taccw = max(10, len(qresult[0].acc))
            qresult_acc = qresult.acc
        except AttributeError:
            qaccw, taccw = 10, 10
            qresult_acc = '-'

        for hit in qresult:
            
            # try to get hit accession
            try:
                hit_acc = hit.acc
            except AttributeError:
                hit_acc = '-'

            for hsp in hit:
                if self.hmm_as_hit:
                    hmm_to = hsp.hit_to + 1
                    hmm_from = hsp.hit_from + 1
                    ali_to = hsp.query_to + 1
                    ali_from = hsp.query_from + 1
                else:
                    hmm_to = hsp.query_to + 1
                    hmm_from = hsp.query_from + 1
                    ali_to = hsp.hit_to + 1
                    ali_from = hsp.hit_from + 1

                rows += "%-*s %-*s %5d %-*s %-*s %5d %9.2g %6.1f %5.1f %3d %3d" \
                " %9.2g %9.2g %6.1f %5.1f %5d %5d %5ld %5ld %5d %5d %4.2f %s\n" % \
                (tnamew, hit.id, taccw, hit_acc, hit.seq_len, qnamew, qresult.id, \
                qaccw, qresult_acc, qresult.seq_len, hit.evalue, hit.bitscore, \
                hit.bias, hsp.domain_index, len(hit), hsp.evalue_cond, hsp.evalue, \
                hsp.bitscore, hsp.bias, hmm_from, hmm_to, ali_from, ali_to, \
                hsp.env_from + 1, hsp.env_to + 1, hsp.acc_avg, hit.desc)

        return rows


class HmmerDomtabHmmqueryWriter(HmmerDomtabHmmhitWriter):

    """Writer for hmmer-domtab output format which writes query coordinates
    as HMM profile coordinates."""

    def __init__(self, handle):
        self.handle = handle
        self.hmm_as_hit = False


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
