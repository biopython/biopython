# Copyright 2014 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO parser for HMMER pfam table output format."""

from itertools import chain

from Bio.Alphabet import generic_protein
from Bio.SearchIO._model import QueryResult, Hit, HSP, HSPFragment
from Bio._py3k import OrderedDict


class Hmmer3PfamtabParser(object):

    """Base hmmer3-pfamtab iterator."""

    def __init__(self, handle):
        self.handle = handle
        self.line = self.handle.readline()

    def __iter__(self):
        while not self.line.startswith('# Sequence scores'):
            self.line = self.handle.readline()
            if not self.line:
                raise StopIteration

        for qresult in self._parse_qresult():
            yield qresult

    def _parse_seq_row(self):
        cols = [x for x in self.line.strip().split() if x]
        # adjust for description with multiple or no words
        if len(cols) > 7:
            cols[6] = ' '.join(cols[6:])
            del cols[7:]
        elif len(cols) < 7:
            cols.append('')
        assert len(cols) == 7
        # assign values per column
        hit = {}
        hit['id'] = cols[0]
        hit['bitscore'] = cols[1]
        hit['evalue'] = cols[2]
        hit['__n'] = cols[3]
        hit['domain_exp_num'] = cols[4]
        hit['bias'] = cols[5]
        hit['description'] = cols[6]

        return hit

    def _parse_domain_row(self):
        cols = [x for x in self.line.strip().split() if x]
        if len(cols) > 12:
            cols[11] = ' '.join(cols[11:])
            del cols[12:]
        elif len(cols) < 12:
            cols.append('')
        assert len(cols) == 12
        hsp, frag = {}, {}
        hsp['hit_id'] = cols[0]
        hsp['bitscore'] = cols[1]
        hsp['evalue'] = cols[2]
        hsp['domain_index'] = cols[3]
        hsp['bias'] = cols[4]
        hsp['env_start'] = int(cols[5]) - 1
        hsp['env_end'] = int(cols[6])
        hsp['hit_description'] = cols[11]
        frag['hit_strand'] = frag['query_strand'] = 0
        frag['query_start'] = int(cols[7]) - 1
        frag['query_end'] = int(cols[8])
        frag['hit_start'] = int(cols[9]) - 1
        frag['hit_end'] = int(cols[10])
        frag['alphabet'] = generic_protein

        # switch hmm<-->ali coordinates if hmm is not hit
        if not self.hmm_as_hit:
            frag['hit_end'], frag['query_end'] = \
                    frag['query_end'], frag['hit_end']
            frag['hit_start'], frag['query_start'] = \
                    frag['query_start'], frag['hit_start']

        return {'hsp': hsp, 'frag': frag}

    def _parse_qresult(self, base_id='query_'):
        # counter for query IDs
        # NOTE: this is because pfamtab does not contain query IDs
        idx = 0
        while True:
            assert self.line.startswith('# Sequence scores')
            hits = OrderedDict()
            hsps = {}
            # read through sequence scores header (5 lines)
            for _ in range(5):
                self.line = self.handle.readline()
            # short circuit when there are no hits
            if self.line and not self.line.strip():
                # finish through the domain scores header
                for _ in range(6):
                    self.line = self.handle.readline()
                yield QueryResult(id=base_id + str(idx))
                idx += 1
            elif not self.line:
                break
            else:
                # parse sequence scores table
                while not self.line.startswith('#'):
                    raw_hit = self._parse_seq_row()
                    hits[raw_hit['id']] = raw_hit
                    self.line = self.handle.readline()
                    # also break when line is empty
                    if self.line and not self.line.strip():
                        # advance 1 line first
                        self.line = self.handle.readline()
                        break
                # parse domain scores table
                assert self.line.startswith('# Domain scores')
                for _ in range(5):
                    self.line = self.handle.readline()
                while not self.line.startswith('#'):
                    cur = self._parse_domain_row()
                    frag = HSPFragment(cur['hsp']['hit_id'], base_id + str(idx))
                    for attr, value in cur['frag'].items():
                        setattr(frag, attr, value)
                    hsp = HSP([frag])
                    for attr, value in cur['hsp'].items():
                        setattr(hsp, attr, value)
                    if hsp.hit_id in hsps:
                        hsps[hsp.hit_id].append(hsp)
                    else:
                        hsps[hsp.hit_id] = [hsp]
                    # move on to next line
                    self.line = self.handle.readline()
                # create hits and yield query
                hit_list = []
                for hit_id, raw_hit in hits.items():
                    hit = Hit(hsps[hit_id])
                    for attr, value in raw_hit.items():
                        setattr(hit, attr, value)
                    hit_list.append(hit)
                yield QueryResult(hit_list)
                idx += 1

                # '#\n' (or equivalent) means we have reached the end
                # of all results
                # we can discard all the other metrics
                if self.line.strip() == '#':
                    for _ in self.handle:
                        pass
                    raise StopIteration


class Hmmer3PfamtabHmmhitParser(Hmmer3PfamtabParser):

    """Parser for the HMMER pfam table format that assumes HMM profile
    coordinates are hit coordinates."""

    hmm_as_hit = True


class Hmmer3PfamtabHmmqueryParser(Hmmer3PfamtabParser):

    """Parser for the HMMER pfam table format that assumes HMM profile
    coordinates are query coordinates."""

    hmm_as_hit = False


class Hmmer3PfamtabHmmhitIndexer(object):

    """Indexer class for HMMER pfam table output that assumes HMM profile
    coordinates are hit coordinates."""

    _parser = Hmmer3PfamtabHmmhitParser
    _query_id_idx = 3


class Hmmer3PfamtabHmmqueryIndexer(object):

    """Indexer class for HMMER pfam table output that assumes HMM profile
    coordinates are query coordinates."""

    _parser = Hmmer3PfamtabHmmqueryParser
    _query_id_idx = 3


class Hmmer3PfamtabHmmhitWriter(object):

    """Writer for hmmer3-pfamtab output format which writes hit coordinates
    as HMM profile coordinates."""

    hmm_as_hit = True

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

        header = "#%*s %22s %40s %11s %11s %11s\n" % \
                (tnamew+qnamew-1+15+taccw+qaccw, "", "--- full sequence ---",
                "-------------- this domain -------------", "hmm coord",
                "ali coord", "env coord")
        header += "#%-*s %-*s %5s %-*s %-*s %5s %9s %6s %5s %3s %3s %9s " \
                "%9s %6s %5s %5s %5s %5s %5s %5s %5s %4s %s\n" % (tnamew-1,
                " target name", taccw, "accession", "tlen", qnamew,
                "query name", qaccw, "accession", "qlen", "E-value", "score",
                "bias", "#", "of", "c-Evalue", "i-Evalue", "score", "bias",
                "from", "to", "from", "to", "from", "to", "acc",
                "description of target")
        header += "#%*s %*s %5s %*s %*s %5s %9s %6s %5s %3s %3s %9s %9s " \
                "%6s %5s %5s %5s %5s %5s %5s %5s %4s %s\n" % (tnamew-1,
                "-------------------", taccw, "----------", "-----",
                qnamew, "--------------------", qaccw, "----------",
                "-----", "---------", "------", "-----", "---", "---",
                "---------", "---------", "------", "-----", "-----", "-----",
                "-----", "-----", "-----", "-----", "----",
                "---------------------")

        return header

    def _build_row(self, qresult):
        """Returns a string or one row or more of the QueryResult object."""
        rows = ''

        # calculate whitespace required
        # adapted from HMMER's source: src/p7_tophits.c#L1083
        qnamew = max(20, len(qresult.id))
        tnamew = max(20, len(qresult[0].id))
        try:
            qaccw = max(10, len(qresult.accession))
            taccw = max(10, len(qresult[0].accession))
            qresult_acc = qresult.accession
        except AttributeError:
            qaccw, taccw = 10, 10
            qresult_acc = '-'

        for hit in qresult:

            # try to get hit accession
            try:
                hit_acc = hit.accession
            except AttributeError:
                hit_acc = '-'

            for hsp in hit.hsps:
                if self.hmm_as_hit:
                    hmm_to = hsp.hit_end
                    hmm_from = hsp.hit_start + 1
                    ali_to = hsp.query_end
                    ali_from = hsp.query_start + 1
                else:
                    hmm_to = hsp.query_end
                    hmm_from = hsp.query_start + 1
                    ali_to = hsp.hit_end
                    ali_from = hsp.hit_start + 1

                rows += "%-*s %-*s %5d %-*s %-*s %5d %9.2g %6.1f %5.1f %3d %3d" \
                " %9.2g %9.2g %6.1f %5.1f %5d %5d %5ld %5ld %5d %5d %4.2f %s\n" % \
                (tnamew, hit.id, taccw, hit_acc, hit.seq_len, qnamew, qresult.id,
                qaccw, qresult_acc, qresult.seq_len, hit.evalue, hit.bitscore,
                hit.bias, hsp.domain_index, len(hit.hsps), hsp.evalue_cond, hsp.evalue,
                hsp.bitscore, hsp.bias, hmm_from, hmm_to, ali_from, ali_to,
                hsp.env_start + 1, hsp.env_end, hsp.acc_avg, hit.description)

        return rows


class Hmmer3PfamtabHmmqueryWriter(Hmmer3PfamtabHmmhitWriter):

    """Writer for hmmer3-pfamtab output format which writes query coordinates
    as HMM profile coordinates."""

    hmm_as_hit = False


# if not used as a module, run the doctest
if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()
