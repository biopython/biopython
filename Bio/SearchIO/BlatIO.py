# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO parser for BLAT output formats.

This module adds support for parsing BLAT outputs. BLAT (BLAST-Like Alignment
Tool) is a sequence homology search program initially built for annotating
the human genome.

Specifically, this module supports the following BLAT output formats:

  - PSL - 'blat-psl'
  - PSLX - 'blat-pslx'

More information are available through these links:
  
  - Publication: http://genome.cshlp.org/content/12/4/656
  - User guide: http://genome.ucsc.edu/goldenPath/help/blatSpec.html
  - Source download: http://www.soe.ucsc.edu/~kent/src
  - Executable download: http://hgdownload.cse.ucsc.edu/admin/exe/
  - Blat score calculation: http://genome.ucsc.edu/FAQ/FAQblat.html#blat4

"""

import re
from math import log

from Bio._py3k import _as_bytes, _bytes_to_string
from Bio.SearchIO._objects import QueryResult, Hit, BatchHSP, HSP
from Bio.SearchIO._index import SearchIndexer


# precompile regex patterns
_RE_ROW_CHECK = re.compile(r'^\d+\s+\d+\s+\d+\s+\d+')


def _append_hit(qresult, hit):
    """Appends Hit to the given QueryResult objects; combines HSP if Hit exists."""
    if hit not in qresult:
        qresult.append(hit)
    else:
        for hsp in hit.batch_hsps:
            qresult[hit.id].append(hsp)


def _list_from_csv(csv_string, caster=None):
    """Transforms the given comma-separated string into a list.

    Arguments:
    csv_string -- Comma-separated string to transform.
    caster -- Cast function to use on each list item.

    """
    filtered = [x for x in filter(None, csv_string.split(','))]
    if caster is None:
        return filtered
    else:
        return [caster(x) for x in filtered]


def _reorient_starts(starts, blksizes, seqlen, strand):
    """Reorients block starts into the opposite strand's coordinates.

    Arguments:
    starts -- List of integers, start coordinates.
    start -- Integer, 'Q start' or 'T start' column
    blksizes -- List of integers, block sizes.
    seqlen -- Integer of total sequence length.
    strand -- Integer denoting sequence strand.

    """
    assert len(starts) == len(blksizes), \
            "Unequal start coordinates and block sizes list (%r vs %r)" \
            % (len(starts), len(blksizes))
    # see: http://genome.ucsc.edu/goldenPath/help/blatSpec.html
    # no need to reorient if it's already the positive strand
    if strand >= 0:
        return starts
    else:
        # the plus-oriented coordinate is calculated by this:
        # plus_coord = length - minus_coord - block_size
        return [seqlen - start - blksize for \
                start, blksize in zip(starts, blksizes)]


def _is_protein(psl):
    # check if query is protein or not
    # adapted from http://genome.ucsc.edu/FAQ/FAQblat.html#blat4
    if len(psl['strand']) == 2:
        if psl['strand'][1] == '+':
            return psl['tend'] == psl['tstarts'][-1] + \
                    3 * psl['blocksizes'][-1]
        elif psl['strand'][1] == '-':
            return psl['tstart'] == psl['tsize'] - \
                    (psl['tstarts'][-1] + 3 * psl['blocksizes'][-1])

    return False


def _calc_millibad(psl, is_protein):
    # calculates millibad
    # adapted from http://genome.ucsc.edu/FAQ/FAQblat.html#blat4
    size_mul = 3 if is_protein else 1
    millibad = 0

    qali_size = size_mul * (psl['qend'] - psl['qstart'])
    tali_size = psl['tend'] - psl['tstart']
    ali_size = min(qali_size, tali_size)
    if ali_size <= 0:
        return 0

    size_dif = qali_size - tali_size
    size_dif = 0 if size_dif < 0 else size_dif

    total = size_mul * (psl['matches'] + psl['repmatches'] + psl['mismatches'])
    if total != 0:
        millibad = (1000 * (psl['mismatches'] * size_mul + psl['qnuminsert'] + \
                round(3 * log(1 + size_dif)))) / total

    return millibad


def _calc_score(psl, is_protein):
    # calculates score
    # adapted from http://genome.ucsc.edu/FAQ/FAQblat.html#blat4
    size_mul = 3 if is_protein else 1
    return size_mul * (psl['matches'] + (psl['repmatches'] >> 1)) - \
            size_mul * psl['mismatches'] - psl['qnuminsert'] - psl['tnuminsert']


def _create_batch_hsp(hid, qid, psl):
    # protein flag
    is_protein = _is_protein(psl)
    # strand
    #if query is protein, strand is 0
    if is_protein:
        qstrand = 0
    else:
        qstrand = 1 if psl['strand'][0] == '+' else -1
    # try to get hit strand, if it exists
    try:
        hstrand = 1 if psl['strand'][1] == '+' else -1
    except IndexError:
        hstrand = 1  # hit strand defaults to plus

    # query block starts
    qstarts = _reorient_starts(psl['qstarts'], \
            psl['blocksizes'], psl['qsize'], qstrand)
    # hit block starts
    if len(psl['strand']) == 2:
        hstarts = _reorient_starts(psl['tstarts'], \
                psl['blocksizes'], psl['tsize'], hstrand)
    else:
        hstarts = psl['tstarts']
    # set query and hit coords
    # this assumes each block has no gaps (which seems to be the case)
    assert len(qstarts) == len(hstarts) == len(psl['blocksizes'])
    query_ranges = zip(qstarts, [x + y for x, y in \
            zip(qstarts, psl['blocksizes'])])
    hit_ranges = zip(hstarts, [x + y for x, y in \
            zip(hstarts, psl['blocksizes'])])
    # check length of sequences and coordinates, all must match
    if psl['tseqs'] and psl['qseqs']:
        assert len(psl['tseqs']) == len(psl['qseqs']) == \
                len(query_ranges) == len(hit_ranges)
    else:
        assert len(query_ranges) == len(hit_ranges)

    hsps = []
    # iterating over query_ranges, but hit_ranges works just as well
    for idx, qcoords in enumerate(query_ranges):
        hseqlist = psl.get('tseqs')
        hseq = '' if not hseqlist else hseqlist[idx]
        qseqlist = psl.get('qseqs')
        qseq = '' if not qseqlist else qseqlist[idx]
        hsp = HSP(hid, qid, hseq, qseq)
        # set coordinates
        hsp.query_start = qcoords[0]
        hsp.query_end = qcoords[1]
        hsp.hit_start = hit_ranges[idx][0]
        hsp.hit_end = hit_ranges[idx][1]
        # and strands
        hsp.query_strand = qstrand
        hsp.hit_strand = hstrand
        hsps.append(hsp)

    # create batch hsp object
    ghsp = BatchHSP(hid, qid, hsps)
    # check if start and end are set correctly
    assert ghsp.query_start == psl['qstart']
    assert ghsp.query_end == psl['qend']
    assert ghsp.hit_start == psl['tstart']
    assert ghsp.hit_end == psl['tend']
    # and check block spans as well
    assert ghsp.query_spans == ghsp.hit_spans == psl['blocksizes']
    # set its attributes
    ghsp.match_num = psl['matches']
    ghsp.mismatch_num = psl['mismatches']
    ghsp.match_rep_num = psl['repmatches']
    ghsp.n_num = psl['ncount']
    ghsp.query_gapopen_num = psl['qnuminsert']
    ghsp.query_gap_num = psl['qbaseinsert']
    ghsp.hit_gapopen_num = psl['tnuminsert']
    ghsp.hit_gap_num = psl['tbaseinsert']

    ghsp.query_strand = qstrand
    ghsp.hit_strand = hstrand
    ghsp.ident_num = psl['matches'] + psl['repmatches']
    ghsp.gapopen_num = psl['qnuminsert'] + psl['tnuminsert']
    ghsp.gap_num = psl['qbaseinsert'] + psl['tbaseinsert']
    ghsp.query_is_protein = is_protein
    ghsp.ident_pct = 100.0 - _calc_millibad(psl, is_protein) * 0.1
    ghsp.score = _calc_score(psl, is_protein)
    # helper flag, for writing
    ghsp._has_hit_strand = len(psl['strand']) == 2

    return ghsp


class BlatPslIterator(object):

    """Parser for the BLAT PSL format."""

    def __init__(self, handle):
        self.handle = handle
        self.line = self.handle.readline()

    def __iter__(self):
        # break out if it's an empty file
        if not self.line:
            raise StopIteration

        # read through header
        # this assumes that the result row match the regex
        while not re.search(_RE_ROW_CHECK, self.line.strip()):
            self.line = self.handle.readline()
            if not self.line:
                raise StopIteration

        # parse into query results
        for qresult in self._parse_qresult():
            qresult.program = 'blat'
            yield qresult

    def _parse_qresult(self):
        """Generator function that returns QueryResult objects."""
        qid_cache = None
        hid_cache = None
        # flag for denoting whether the query is the same as the previous line
        # or not
        same_query = False
        while True:
            # only parse the result row if it's not EOF
            if self.line:
                cols = filter(None, self.line.strip().split('\t'))
                self._validate_cols(cols)
                psl =  self._parse_cols(cols)
                qresult_id = psl['qname']

            # a new qresult is created whenever qid_cache != qresult_id
            if qid_cache != qresult_id:
                # append the last hit and yield qresult if qid_cache is filled
                if qid_cache is not None:
                    qresult.append(hit)
                    yield qresult
                    same_query = False
                qid_cache = qresult_id
                hid_cache = None
                qresult = QueryResult(qresult_id)
                qresult.seq_len = psl['qsize']
            # when we've reached EOF, try yield any remaining qresult and break
            elif not self.line:
                _append_hit(qresult, hit)
                yield qresult
                break
            # otherwise, we must still be in the same query, so set the flag
            elif not same_query:
                same_query = True

            hit_id = psl['tname']
            # a new hit is created whenever hid_cache != hit_id
            if hid_cache != hit_id:
                # if we're in the same query, append the previous line's hit
                if same_query:
                    _append_hit(qresult, hit)
                hid_cache = hit_id
                hit = Hit(hit_id, qresult_id)
                hit.seq_len = psl['tsize']

            # create the HSP objects from a single parsed HSP results,
            # group them in one BatchHSP object, and append to Hit
            batch_hsp = _create_batch_hsp(hit_id, qresult_id, psl)
            hit.append(batch_hsp)

            self.line = self.handle.readline()

    def _parse_cols(self, cols):
        """Returns a dictionary of parsed column values."""
        psl = {}

        psl['qname'] = cols[9]                            # qName
        psl['qsize'] = int(cols[10])                      # qSize
        psl['tname'] = cols[13]                           # tName
        psl['tsize'] = int(cols[14])                      # tSize
        psl['matches'] = int(cols[0])                     # matches
        psl['mismatches'] = int(cols[1])                  # misMatches
        psl['repmatches'] = int(cols[2])                  # repMatches
        psl['ncount'] = int(cols[3])                      # nCount
        psl['qnuminsert'] = int(cols[4])                  # qNumInsert
        psl['qbaseinsert'] = int(cols[5])                 # qBaseInsert
        psl['tnuminsert'] = int(cols[6])                  # tNumInsert
        psl['tbaseinsert'] = int(cols[7])                 # tBaseInsert
        psl['strand'] = cols[8]                           # strand
        psl['qstart'] = int(cols[11])                     # qStart
        psl['qend'] = int(cols[12])                       # qEnd
        psl['tstart'] = int(cols[15])                     # tStart
        psl['tend'] = int(cols[16])                       # tEnd
        psl['blockcount'] = int(cols[17])                 # blockCount
        psl['blocksizes'] = _list_from_csv(cols[18], int) # blockSizes
        psl['qstarts'] = _list_from_csv(cols[19], int)    # qStarts
        psl['tstarts'] = _list_from_csv(cols[20], int)    # tStarts
        # PSL doesn't have any sequences; these are needed for instantiating
        # BatchHSP objects
        psl['qseqs'] = []
        psl['tseqs'] = []

        return psl

    def _validate_cols(self, cols):
        assert len(cols) == 21, "Invalid PSL line: %r. " \
        "Expected 21 tab-separated columns, found %i" % (self.line, len(cols))


class BlatPslxIterator(BlatPslIterator):

    """Parser for the BLAT PSLX format."""

    def _validate_cols(self, cols):
        assert len(cols) == 23, "Invalid PSLX line: %r. " \
        "Expected 23 tab-separated columns, found %i" % (self.line, len(cols))

    def _parse_cols(self, cols):
        psl = BlatPslIterator._parse_cols(self, cols)
        # append seqs to hsp dict
        psl['qseqs'] = _list_from_csv(cols[21])    # query sequence
        psl['tseqs'] = _list_from_csv(cols[22])    # hit sequence

        return psl


class BlatPslIndexer(SearchIndexer):

    """Indexer class for BLAT PSL output."""

    _parser = BlatPslIterator

    def __iter__(self):
        """Iterates over the file handle; yields key, start offset, and length."""
        handle = self._handle
        handle.seek(0)
        split_char = _as_bytes('\t')
        # denotes column location for query identifier
        query_id_idx = 9
        qresult_key = None

        start_offset = handle.tell()
        line = handle.readline()
        # read through header
        # this assumes that the result row match the regex
        while not re.search(_RE_ROW_CHECK, line.strip()):
            start_offset = handle.tell()
            line = handle.readline()
            if not line:
                raise StopIteration

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

            line = handle.readline()
            if not line:
                yield _bytes_to_string(qresult_key), start_offset, \
                        end_offset - start_offset
                break

    def get_raw(self, offset):
        """Returns the raw string of a QueryResult object from the given offset."""
        handle = self._handle
        handle.seek(offset)
        split_char = _as_bytes('\t')
        query_id_idx = 9
        qresult_key = None
        qresult_raw = ''

        while True:
            line = handle.readline()
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


class BlatPslxIndexer(BlatPslIndexer):

    """Indexer class for BLAT PSL output."""

    _parser = BlatPslxIterator


class BlatPslWriter(object):

    """Writer for the blat-psl format."""

    fmt = 'psl'

    def __init__(self, handle, header=False):
        self.handle = handle
        # flag for writing header or not
        self.header = header

    def write_file(self, qresults):
        handle = self.handle
        qresult_counter, hit_counter, hsp_counter = 0, 0, 0

        if self.header:
            handle.write(self.build_header())

        for qresult in qresults:
            if qresult:
                handle.write(self.build_row(qresult))
                qresult_counter += 1
                hit_counter += len(qresult)
                hsp_counter += sum([len(hit) for hit in qresult])

        return qresult_counter, hit_counter, hsp_counter

    def build_header(self):
        # for now, always use the psLayout version 3
        header = 'psLayout version 3\n'

        # adapted from BLAT's source: lib/psl.c#L496
        header += "\nmatch\tmis- \trep. \tN's\tQ gap\tQ gap\tT gap\tT " \
        "gap\tstrand\tQ        \tQ   \tQ    \tQ  \tT        \tT   \tT    " \
        "\tT  \tblock\tblockSizes \tqStarts\t tStarts\n     " \
        "\tmatch\tmatch\t   \tcount\tbases\tcount\tbases\t      \tname     " \
        "\tsize\tstart\tend\tname     \tsize\tstart\tend\tcount" \
        "\n%s\n" % ('-' * 159)

        return header

    def build_row(self, qresult):
        """Returns a string or one row or more of the QueryResult object."""
        # For now, our writer writes the row according to the order in
        # the QueryResult and Hit objects.
        # This is different from BLAT's native output, where the rows are
        # grouped by strand.
        # Should we tweak the behavior to better mimic the native output?
        qresult_lines = []

        for hit in qresult:
            for hsp in hit.batch_hsps:

                line = []
                line.append(hsp.match_num)
                line.append(hsp.mismatch_num)
                line.append(hsp.match_rep_num)
                line.append(hsp.n_num)
                line.append(hsp.query_gapopen_num)
                line.append(hsp.query_gap_num)
                line.append(hsp.hit_gapopen_num)
                line.append(hsp.hit_gap_num)

                # check spans
                assert hsp.query_spans == hsp.hit_spans
                block_sizes = hsp.query_spans

                # set strand and starts
                if hsp.query_strand >= 0: # since it may be a protein seq
                    strand = '+'
                else:
                    strand = '-'
                qstarts = _reorient_starts([x[0] for x in hsp.query_ranges], \
                        hsp.query_spans, qresult.seq_len, hsp.query_strand)

                if hsp.hit_strand == 1:
                    hstrand = 1
                    # only write hit strand if it was present in the source file
                    if hsp._has_hit_strand:
                        strand += '+'
                else:
                    hstrand = -1
                    strand += '-'
                hstarts = _reorient_starts([x[0] for x in hsp.hit_ranges], \
                        hsp.hit_spans, hit.seq_len, hstrand)

                line.append(strand)
                line.append(qresult.id)
                line.append(qresult.seq_len)
                line.append(hsp.query_start)
                line.append(hsp.query_end)
                line.append(hit.id)
                line.append(hit.seq_len)
                line.append(hsp.hit_start)
                line.append(hsp.hit_end)
                line.append(len(hsp))
                line.append(','.join((str(x) for x in block_sizes)) + ',')
                line.append(','.join((str(x) for x in qstarts)) + ',')
                line.append(','.join((str(x) for x in hstarts)) + ',')

                if self.fmt == 'pslx':
                    line.append(','.join((str(x.seq) for x in hsp.queries)) + ',')
                    line.append(','.join((str(x.seq) for x in hsp.hits)) + ',')

                qresult_lines.append('\t'.join((str(x) for x in line)))

        return '\n'.join(qresult_lines) + '\n'


class BlatPslxWriter(BlatPslWriter):

    """Writer for the blat-pslx format."""

    fmt = 'pslx'


def _test():
    """Run the Bio.SearchIO.BlatIO module's doctests.

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
