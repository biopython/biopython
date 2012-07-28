# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO abstract base parser for Exonerate standard output format."""

from Bio._py3k import _bytes_to_string
from Bio.SearchIO._objects import QueryResult, Hit, HSP, BatchHSP
from Bio.SearchIO._index import SearchIndexer


# strand char-value mapping
_STRAND_MAP = {'+': 1, '-': -1, '.': 0}


def _absorb_hit(qresult, hit):
    """Absorbs the given hit by merging it with an existing hit."""
    try:
        qresult.append(hit)
    except ValueError:
        assert hit.id in qresult
        for hsp_item in hit.batch_hsps:
            qresult[hit.id].append(hsp_item)

    return qresult


def _create_batch_hsp(hid, qid, hspd):
    """Returns a list of HSP objects from the given parsed HSP values."""
    hsps = []

    # we are iterating over query_ranges, but hit_ranges works just as well
    for idx, qcoords in enumerate(hspd['query_ranges']):
        # get sequences, create object
        hseqlist = hspd.get('hit')
        hseq = '' if hseqlist is None else hseqlist[idx]
        qseqlist = hspd.get('query')
        qseq = '' if qseqlist is None else qseqlist[idx]
        hsp = HSP(hid, qid, hseq, qseq)
        # coordinates
        hsp.query_start = qcoords[0]
        hsp.query_end = qcoords[1]
        hsp.hit_start = hspd['hit_ranges'][idx][0]
        hsp.hit_end = hspd['hit_ranges'][idx][1]
        # alignment annotation
        try:
            aln_annot = hspd.get('alignment_annotation', {})
            for key, value in aln_annot.items():
                hsp.alignment_annotation[key] = value[idx]
        except IndexError:
            pass
        # strands
        hsp.query_strand = hspd['query_strand']
        hsp.hit_strand = hspd['hit_strand']
        # and append the hsp object to the list
        hsps.append(hsp)

    ghsp = BatchHSP(hsps)
    # set batchhsp-specific attributes
    for attr in ('score', 'hit_scodon_ranges', 'query_scodon_ranges', \
            'hit_intron_ranges', 'query_intron_ranges', 'hit_ner_ranges', \
            'query_ner_ranges', 'query_strand', 'hit_strand', 'model', \
            'vulgar_comp', 'cigar_comp'):
        if attr in hspd:
            setattr(ghsp, attr, hspd[attr])

    return ghsp


def _parse_hit_or_query_line(line):
    """Parse the 'Query:' line of exonerate alignment outputs."""
    try:
        mark, id, desc = line.split(' ', 2)
    except ValueError: # no desc
        mark, id = line.split(' ', 1)
        desc = ''

    return id, desc


class BaseExonerateIterator(object):

    """Abstract iterator for exonerate format."""

    _ALN_MARK = None

    def __init__(self, handle):
        self.handle = handle
        self.has_c4_alignment = False

    def __iter__(self):
        # read line until the first alignment block or cigar/vulgar lines
        while True:
            self.line = self.handle.readline()
            # flag for human-readable alignment block
            if self.line.startswith('C4 Alignment:') and not \
                    self.has_c4_alignment:
                self.has_c4_alignment = True
            if self.line.startswith('C4 Alignment:') or \
                    self.line.startswith('vulgar:') or \
                    self.line.startswith('cigar:'):
                break
            elif not self.line or self.line.startswith('-- completed '):
                raise StopIteration

        for qresult in self.parse_qresult():
            qresult.program = 'exonerate'
            # HACK: so that all descriptions are set
            qresult.description = qresult.description
            for hit in qresult:
                hit.description = hit.description
            yield qresult

    def read_until(self, bool_func):
        """Reads the file handle until the given bool function returns True."""
        while True:
            if not self.line or bool_func(self.line):
                return
            else:
                self.line = self.handle.readline()

    def parse_alignment_block(self, qres_dict, hit_dict, hsp_dict):
        raise NotImplementedError("Subclass must implement this")

    def parse_alignment_header(self, aln_header, qresult, hit, hsp):
        for line in aln_header:
            # query line
            if line.startswith('Query:'):
                qresult['id'], qresult['description']  = \
                        _parse_hit_or_query_line(line)
            # target line
            elif line.startswith('Target:'):
                hit['id'], hit['description'] = \
                        _parse_hit_or_query_line(line)
            # model line
            elif line.startswith('Model:'):
                qresult['model'] = line.split(' ', 1)[1]
            # score line
            elif line.startswith('Raw score:'):
                hsp['score'] = line.split(' ', 2)[2]
            # query range line
            elif line.startswith('Query range:'):
                # line is always 'Query range: \d+ -> \d+', so we can pluck
                # the numbers directly
                hsp['query_start'], hsp['query_end'] = line.split(' ', 4)[2:5:2]
            # hit range line
            elif line.startswith('Target range:'):
                # same logic with query range
                hsp['hit_start'], hsp['hit_end'] = line.split(' ', 4)[2:5:2]

        # determine strand
        if qresult['description'].endswith(':[revcomp]'):
            hsp['query_strand'] = '-'
            qresult['description'] = qresult['description'].replace(':[revcomp]', '')
        elif qresult['model'].startswith('protein2'):
            hsp['query_strand'] = '.'
        else:
            hsp['query_strand'] = '+'
        if hit['description'].endswith(':[revcomp]'):
            hsp['hit_strand'] = '-'
            hit['description'] = hit['description'].replace(':[revcomp]', '')
        else:
            hsp['hit_strand'] = '+'

        # NOTE: we haven't processed the coordinates types
        # and the strands are not yet Biopython's standard (1 / -1 / 0)
        # since it's easier if we do the conversion later

        return qresult, hit, hsp

    def parse_qresult(self):
        qid_cache = None
        hid_cache = None
        same_query = False # flag for tracking active query
        # if the file has c4 alignments, use that as the alignment mark
        if self.has_c4_alignment:
            self._ALN_MARK = 'C4 Alignment:'

        while True:
            self.read_until(lambda line: line.startswith(self._ALN_MARK))
            # only parse the result row if it's not EOF
            if self.line:
                assert self.line.startswith(self._ALN_MARK), self.line
                # create temp dicts for storing parsed values
                qres_dict, hit_dict, hsp_dict = {}, {}, {}
                # if the file has c4 alignments, try to parse the header
                if self.has_c4_alignment:
                    self.read_until(lambda line: line.strip().startswith('Query:'))
                    # collect header lines and parse them
                    aln_header = []
                    while not self.line == '\n':
                        aln_header.append(self.line.strip())
                        self.line = self.handle.readline()
                    # parse header values into the temp dicts
                    qres_dict, hit_dict, hsp_dict = \
                            self.parse_alignment_header(aln_header,qres_dict, \
                            hit_dict, hsp_dict)

                qres_parsed, hit_parsed, hsp_parsed = \
                        self.parse_alignment_block(qres_dict, hit_dict, hsp_dict)
                qresult_id = qres_parsed['id']

            # a new qresult is created whenever qid_cache != qresult_id
            if qid_cache != qresult_id:
                # append the last hit and yield qresult if qid_cache is filled
                if qid_cache is not None:
                    yield _absorb_hit(qresult, hit)
                    same_query = False
                qid_cache = qresult_id
                hid_cache = None
                qresult = QueryResult(qresult_id)
                for attr, value in qres_parsed.items():
                    setattr(qresult, attr, value)
            # when we've reached EOF, try yield any remaining qresult and break
            elif not self.line or self.line.startswith('-- completed '):
                yield _absorb_hit(qresult, hit)
                break
            # otherwise, we must still be in the same query, so set the flag
            elif not same_query:
                same_query = True

            hit_id = hit_parsed['id']
            # a new hit is created whenever hid_cache != hit_id
            if hid_cache != hit_id:
                # if hit is already in qresult, merge them
                # for some reason, exonerate doesn't always group hits together (?)
                if same_query:
                    qresult = _absorb_hit(qresult, hit)
                hid_cache = hit_id
                hit = Hit(hit_id, qresult_id)
                for attr, value in hit_parsed.items():
                    setattr(hit, attr, value)

            # create the HSP objects from a single parsed HSP results,
            # group them in one BatchHSP object, and append to Hit
            batch_hsp = _create_batch_hsp(hit_id, qresult_id, hsp_parsed)
            hit.append(batch_hsp)

            if not self.has_c4_alignment:
                self.line = self.handle.readline()


class BaseExonerateIndexer(SearchIndexer):

    """Indexer class for Exonerate plain text."""

    _parser = None # should be defined by subclass
    _query_mark = None # this one too

    def get_qresult_id(self, pos):
        raise NotImplementedError("Should be defined by subclass")

    def __iter__(self):
        """Iterates over the file handle; yields key, start offset, and length."""
        handle = self._handle
        handle.seek(0)
        qresult_key = None

        while True:
            start_offset = handle.tell()
            line = handle.readline()
            if line.startswith(self._query_mark):
                if qresult_key is None:
                    qresult_key = self.get_qresult_id(start_offset)
                    qresult_offset = start_offset
                else:
                    curr_key = self.get_qresult_id(start_offset)
                    if curr_key != qresult_key:
                        yield _bytes_to_string(qresult_key), qresult_offset, \
                                start_offset - qresult_offset
                        qresult_key = curr_key
                        qresult_offset = start_offset
                        handle.seek(qresult_offset)
            elif not line:
                yield _bytes_to_string(qresult_key), qresult_offset, \
                        start_offset - qresult_offset
                break


def _test():
    """Run the Bio.SearchIO.ExonerateIO._base module's doctests.

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
