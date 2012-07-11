# Copyright 2008-2011 by Peter Cock.
# Copyright 2012 by Wibowo Arindrarto.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO parser for Bill Pearson's FASTA tools.

This module adds support for parsing FASTA outputs. FASTA is a suite of
programs that finds regions of local or global similarity between protein
or nucleotide sequences, either by searching databases or identifying
local duplications.

The following FASTA output format is supported: format 10 (triggered by the
-m 10 flag).

More information are available through these links:
  - Website: http://fasta.bioch.virginia.edu/fasta_www2/fasta_list2.shtml
  - User guide: http://fasta.bioch.virginia.edu/fasta_www2/fasta_guide.pdf

"""

import re

from Bio.Alphabet import generic_dna, generic_protein
from Bio.File import UndoHandle
from Bio.SearchIO._objects import QueryResult, Hit, HSP
from Bio.SearchIO._index import SearchIndexer
from Bio._py3k import _bytes_to_string

# precompile regex patterns
# regex for program name
_RE_FLAVS = re.compile(r't?fast[afmsxy]|pr[sf][sx]|lalign|[gs]?[glso]search')
# regex for sequence ID and length
_RE_ID_DESC_SEQLEN = re.compile(r'>>>(.+?)\s+(.*?) *- (\d+) (?:aa|nt)$')
# regex for qresult, hit, or hsp attribute value
_RE_ATTR = re.compile(r'^; [a-z]+(_[ \w-]+):\s+(.*)$')

# attribute name mappings
_HSP_ATTR_MAP = {
    '_initn': ('initn_score', int),
    '_init1': ('init1_score', int),
    '_opt': ('opt_score', int),
    '_s-w opt': ('opt_score', int),
    '_z-score': ('z_score', float),
    '_bits': ('bitscore', float),
    '_expect': ('evalue', float),
    '_score': ('sw_score', int),
    '_ident': ('ident_pct', float),
    '_sim': ('pos_pct', float),
}

# state flags
STATE_NONE = 0
STATE_QUERY_BLOCK = 1
STATE_HIT_BLOCK = 2
STATE_CONS_BLOCK = 3


def _set_qresult_hits(qresult, hit_rows):
    """Helper function for appending Hits without alignments into QueryResults."""
    for hit_row in hit_rows:
        hit_id, remainder = hit_row.split(' ', 1)
        # TODO: parse hit and hsp properties properly; by dealing with:
        #   - any character in the description (brackets, spaces, etc.)
        #   - possible [f] or [r] presence (for frame info)
        #   - possible presence of E2() column
        #   - possible incomplete hit_id due to column length limit
        # The current method only looks at the Hit ID, none of the things above
        if hit_id not in qresult:
            hit = Hit(hit_id, qresult.id)
            hsp = HSP(hit_id, qresult.id)
            hit.append(hsp)
            qresult.append(hit)

    return qresult


def _set_hsp_seqs(hsp, hseq, qseq, annot, program, strand):
    """Helper function for the main parsing code (PRIVATE).

    Arguments:
    hsp -- HSP object whose properties are to be set.
    hseq -- String of raw Hit sequence.
    qseq -- String of raw Query sequence.
    annot -- Dictionary containing HSP annotation.
    program -- String of program name.
    strand -- String of frame, 'f' or 'r', from Fasta frame.

    """
    if 'tfast' not in program:
        qseq = _extract_alignment_region(qseq, annot['query'])
        hseq = _extract_alignment_region(hseq, annot['hit'])
    assert len(hseq) == len(qseq), annot

    for seq, seq_type in zip([hseq, qseq], ['hit', 'query']):
        # get and set start and end coordinates
        start = int(annot[seq_type]['_start'])
        end = int(annot[seq_type]['_stop'])

        setattr(hsp, seq_type + '_start', min(start, end) - 1)
        setattr(hsp, seq_type + '_end', max(start, end))
        # set seq and alphabet
        setattr(hsp, seq_type, seq)

        if seq_type == 'query':
            setattr(hsp.query.seq, 'alphabet', _get_alphabet(seq, \
                    annot[seq_type]))

    # set hsp strand properties
    if hsp.query.seq.alphabet is not generic_protein:
        if strand == 'f':
            hsp.query_strand = 1
        else:
            hsp.query_strand = -1
    else:
        hsp.query_strand = 0
        hsp.hit_strand = 0

    # set hsp alignment length
    hsp.ali_len = len(hsp.query)

    return hsp


def _get_alphabet(seq, annot):
    """Helper function to determine sequence alphabet."""
    if annot['_type'] == 'D':
        alphabet = generic_dna
    elif annot['_type'] == 'p':
        alphabet =  generic_protein
    #TODO: should we implement this?
    #if '-' in seq:
    #    alphabet = Gapped(alphabet)

    return alphabet


def _extract_alignment_region(seq, annot):
    """Helper function for the main parsing code (PRIVATE).

    To get the actual pairwise alignment sequences, we must first
    translate the un-gapped sequence based coordinates into positions
    in the gapped sequence (which may have a flanking region shown
    using leading - characters).  To date, I have never seen any
    trailing flanking region shown in the m10 file, but the
    following code should also cope with that.

    Note that this code seems to work fine even when the "sq_offset"
    entries are prsent as a result of using the -X command line option.
    """
    seq_stripped = seq.strip('-')
    disp_start = int(annot['_display_start'])
    start = int(annot['_start'])
    stop = int(annot['_stop'])

    if start <= stop:
        start = start - disp_start
        stop = stop - disp_start + 1
    else:
        start = disp_start - start
        stop = disp_start - stop + 1
    stop += seq_stripped.count('-')
    assert 0 <= start and start < stop and stop <= len(seq_stripped), \
           "Problem with sequence start/stop,\n%s[%i:%i]\n%s" \
           % (seq, start, stop, annot)
    return seq_stripped[start:stop]


class FastaM10Iterator(object):

    """Iterator for the Fasta -m10 output."""

    def __init__(self, handle, parse_hit_table=False):
        self.handle = UndoHandle(handle)
        self._preamble = self.parse_preamble()

    def __iter__(self):
        for qresult in self.parse_qresult():
            # re-set desc, for hsp query description
            qresult.desc = qresult.desc
            yield qresult

    def parse_preamble(self):
        """Parses the Fasta preamble for Fasta flavor and version."""
        preamble = {}

        while True:
            self.line = self.handle.readline()
            # this should be the line just before the first qresult
            if self.line.startswith('Query'):
                break
            # try to match for version line
            elif self.line.startswith(' version'):
                preamble['version'] = self.line.split(' ')[2]
            else:
                # try to match for flavor line
                flav_match = re.match(_RE_FLAVS, self.line.lower())
                if flav_match:
                    preamble['program'] = flav_match.group(0)

        return preamble

    def parse_hit_table(self):
        """Parses hit table rows."""
        hit_rows = []
        while self.line and not self.line.strip():
            hit_rows.append(self.line.strip())
            self.line = self.handle.readline()
        return hit_rows

    def parse_qresult(self):
        # initial qresult value
        qresult = None
        hit_rows = []
        while True:
            # parse hit table if flag is set
            if self.parse_hit_table and \
                    self.line.startswith('The best scores are:'):
                # move to the first row
                self.line = self.handle.readline()
                hit_rows = self.parse_hit_table()

            # this marks the end of a query or the file altogether
            elif self.line.strip() in [">>>///"] or not self.line:
                if self.parse_hit_table:
                    yield _set_qresult_hits(qresult, hit_rows)
                else:
                    yield qresult
                break

            # this marks the beginning of a new query
            elif not self.line.startswith('>>>') and '>>>' in self.line:
                # if qresult is filled, yield it first before creating a new one
                if qresult is not None:
                    if self.parse_hit_table:
                        yield _set_qresult_hits(qresult, hit_rows)
                    else:
                        yield qresult
                regx = re.search(_RE_ID_DESC_SEQLEN, self.line)
                query_id = regx.group(1)
                seq_len = regx.group(3)
                desc = regx.group(2)
                qresult = QueryResult(query_id)
                qresult.seq_len = int(seq_len)
                # get target from the next line
                qresult.target = filter(None, \
                        self.handle.peekline().split(' '))[1].strip()
                if desc is not None:
                    qresult.desc = desc
                # set values from preamble
                for key, value in self._preamble.items():
                    setattr(qresult, key, value)

            # this marks the beginning of the query info and its hits + hsps
            elif self.line.startswith('>>>') and not \
                    self.line.strip() == '>>><<<':
                assert self.line[3:].startswith(qresult.id), self.line
                for hit, strand in self.parse_hit(query_id):
                    # re-set desc, for hsp hit description
                    hit.desc = hit.desc
                    # if hit is not in qresult, append it
                    try:
                        qresult.append(hit)
                    # otherwise, it might be the same hit with a different strand
                    except ValueError:
                        # make sure strand is different and then append hsp to
                        # existing hit
                        for hsp in hit:
                            assert strand != hsp.query_strand
                            qresult[hit.id].append(hsp)

            self.line = self.handle.readline()

    def parse_hit(self, query_id):
        while True:
            self.line = self.handle.readline()
            if self.line.startswith('>>'):
                break

        hit = None
        strand = None
        # parse the hits
        while True:
            peekline = self.handle.peekline()
            # yield hit if we've reached the start of a new query or
            # the end of the search
            if peekline.strip() in [">>><<<", ">>>///"] or \
                    (not peekline.startswith('>>>') and '>>>' in peekline):
                # append last hseq line
                if state == STATE_HIT_BLOCK:
                    hseq += self.line.strip()
                elif state == STATE_CONS_BLOCK:
                    hit[-1].alignment_annotation['homology'] += self.line.strip('\n')
                # process HSP alignment and coordinates
                hit[-1] = _set_hsp_seqs(hit[-1], hseq, qseq, hsp_annot, \
                        self._preamble['program'], strand)
                yield hit, strand
                break
            elif self.line.startswith('>'):
                # yield hit and create a new one if we're still in the same query
                if self.line.startswith('>>'):
                    # try yielding,  if hit is not None
                    try:
                        hit[-1] = _set_hsp_seqs(hit[-1], hseq, qseq, hsp_annot, \
                                self._preamble['program'], strand)
                        yield hit, strand
                    except TypeError:
                        assert hit is None
                    # try to get the hit id and desc, and handle cases without descs
                    try:
                        hit_id, hit_desc = self.line[2:].strip().split(' ', 1)
                    except ValueError:
                        hit_id = self.line[2:].strip().split(' ', 1)[0]
                        hit_desc = ''
                    hit = Hit(hit_id, query_id)
                    hit.desc = hit_desc
                    # create the HSP object for Hit
                    hsp = HSP(hit_id, query_id)
                    hit.append(hsp)
                    # set or reset the state to none
                    state = STATE_NONE
                    hsp_annot = {'query':{}, 'hit': {}}
                # create and append a new HSP if line starts with '>--'
                elif self.line.startswith('>--'):
                    # set seq attributes of previous hsp
                    hit[-1] = _set_hsp_seqs(hit[-1], hseq, qseq, hsp_annot, \
                            self._preamble['program'], strand)
                    # and create a new one
                    hsp = HSP(hit_id, query_id)
                    hit.append(hsp)
                    # set the state ~ none yet
                    state = STATE_NONE
                    hsp_annot = {'query':{}, 'hit': {}}
                # this is either query or hit data in the HSP, depending on the state
                elif self.line.startswith('>'):
                    if state == STATE_NONE:
                        # make sure it's the correct query
                        assert query_id.startswith(self.line[1:].split(' ')[0]), \
                                "%r vs %r" % (query_id, self.line)
                        state = STATE_QUERY_BLOCK
                        qseq = ''
                    elif state == STATE_QUERY_BLOCK:
                        # make sure it's the correct hit
                        assert hit_id.startswith(self.line[1:].split(' ')[0])
                        state = STATE_HIT_BLOCK
                        hseq = ''
                else:
                    raise ValueError("Unexpected line: %r" % self.line)
            # set the Hit and HSP attribute values if line starts with ';'
            elif self.line.startswith(';'):
                # check for conservation block
                if self.line.startswith('; al_cons'):
                    state = STATE_CONS_BLOCK
                    hsp.alignment_annotation = {}
                    hsp.alignment_annotation['homology'] = ''
                else:
                    # Fasta outputs do not make a clear distinction between Hit
                    # and HSPs, so we check the attribute names to determine
                    # whether it belongs to a Hit or HSP
                    regx = re.search(_RE_ATTR, self.line.strip())
                    name = regx.group(1)
                    value = regx.group(2)

                    # for values before the '>...' query block
                    if state == STATE_NONE:
                        if name in _HSP_ATTR_MAP:
                            attr_name, caster = _HSP_ATTR_MAP[name]
                            if caster is not str:
                                value = caster(value)
                            if name in ['_ident', '_sim']:
                                value *= 100
                            setattr(hit[-1], attr_name, value)
                        # store strand
                        elif name == '_frame':
                            strand = value
                    # otherwise, pool the values for processing later
                    elif state == STATE_QUERY_BLOCK:
                        hsp_annot['query'][name] = value
                    elif state == STATE_HIT_BLOCK:
                        if name == '_len':
                            hit.seq_len = int(value)
                        else:
                            hsp_annot['hit'][name] = value
                    # for values in the hit block
                    else:
                        raise ValueError("Unexpected line: %r" % self.line)
            # otherwise, it must be lines containing the sequences
            else:
                assert '>' not in self.line
                # if we're in hit, parse into hsp.hit
                if state == STATE_HIT_BLOCK:
                    hseq += self.line.strip()
                elif state == STATE_QUERY_BLOCK:
                    qseq += self.line.strip()
                elif state == STATE_CONS_BLOCK:
                    hit[-1].alignment_annotation['homology'] += self.line.strip('\n')
                # we should not get here!
                else:
                    raise ValueError("Unexpected line: %r" % self.line)

            self.line = self.handle.readline()


class FastaM10Indexer(SearchIndexer):

    """Indexer class for FASTA m10 output."""

    _parser = FastaM10Iterator

    def __init__(self, *args, **kwargs):
        SearchIndexer.__init__(self, *args, **kwargs)
        self._handle = UndoHandle(self._handle)

    def __iter__(self):
        handle = self._handle
        handle.seek(0)
        start_offset = handle.tell()
        qresult_key = None

        while True:
            line = handle.readline()
            peekline = handle.peekline()
            end_offset = handle.tell()

            if not line.startswith('>>>') and '>>>' in line:
                regx = re.search(_RE_ID_DESC_SEQLEN, line)
                qresult_key = regx.group(1)
                start_offset = end_offset - len(line)
            # yield whenever we encounter a new query
            elif not peekline.startswith('>>>') and '>>>' in peekline and qresult_key is not None:
                yield _bytes_to_string(qresult_key), start_offset, \
                        end_offset - start_offset
                start_offset = end_offset
            # or we arrive at the end of the search
            elif not line:
                yield _bytes_to_string(qresult_key), start_offset, \
                        end_offset - start_offset
                break

    def get_raw(self, offset):
        handle = self._handle
        qresult_raw = ''

        # read header first
        handle.seek(0)
        while True:
            line = handle.readline()
            peekline = handle.peekline()
            qresult_raw += line
            if not peekline.startswith('>>>') and '>>>' in peekline:
                break

        # and read the qresult raw string
        handle.seek(offset)
        while True:
            # preserve whitespace, don't use read_forward
            line = handle.readline()
            peekline = handle.peekline()
            qresult_raw += line

            # break when we've reached qresult end
            if (not peekline.startswith('>>>') and '>>>' in peekline) or \
                    not line:
                break

        # append mock end marker to qresult_raw, since it's not always present
        return qresult_raw + '>>><<<\n'


def _test():
    """Run the Bio.SearchIO.FastaIO module's doctests.

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
