# Adapted from Bio.AlignIO.FastaIO copyright 2008-2011 by Peter Cock.
# Copyright 2012 by Wibowo Arindrarto.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO support for Bill Pearson's FASTA tools.

This module adds support for parsing FASTA outputs. FASTA is a suite of
programs that finds regions of local or global similarity between protein
or nucleotide sequences, either by searching databases or identifying
local duplications.

Bio.SearchIO.FastaIO was tested on the following FASTA flavors and versions:

    - flavors: fasta, ssearch, tfastx
    - versions: 35, 36

Other flavors and/or versions may introduce some bugs. Please file a bug report
if you see such problems to Biopython's bug tracker.

More information on FASTA are available through these links:
  - Website: http://fasta.bioch.virginia.edu/fasta_www2/fasta_list2.shtml
  - User guide: http://fasta.bioch.virginia.edu/fasta_www2/fasta_guide.pdf


Supported Formats
=================

Bio.SearchIO.FastaIO supports parsing and indexing FASTA outputs triggered by
the -m 10 flag. Other formats that mimic other programs (e.g. the BLAST tabular
format using the -m 8 flag) may be parseable but using SearchIO's other parsers
(in this case, using the 'blast-tab' parser).


fasta-m10
=========

Note that in FASTA -m 10 outputs, HSPs from different strands are considered to
be from different hits. They are listed as two separate entries in the hit
table. FastaIO recognizes this and will group HSPs with the same hit ID into a
single Hit object, regardless of strand.

FASTA also sometimes output extra sequences adjacent to the HSP match. These
extra sequences are discarded by FastaIO. Only regions containing the actual
sequence match are extracted.

The following object attributes are provided:

+-----------------+-------------------------+----------------------------------+
| Object          | Attribute               | Value                            |
+=================+=========================+==================================+
| QueryResult     | description             | query sequence description       |
|                 +-------------------------+----------------------------------+
|                 | id                      | query sequence ID                |
|                 +-------------------------+----------------------------------+
|                 | program                 | FASTA flavor                     |
|                 +-------------------------+----------------------------------+
|                 | seq_len                 | full length of query sequence    |
|                 +-------------------------+----------------------------------+
|                 | target                  | target search database           |
|                 +-------------------------+----------------------------------+
|                 | version                 | FASTA version                    |
+-----------------+-------------------------+----------------------------------+
| Hit             | seq_len                 | full length of the hit sequence  |
+-----------------+-------------------------+----------------------------------+
| HSP             | bitscore                | *_bits line                      |
|                 +-------------------------+----------------------------------+
|                 | evalue                  | *_expect line                    |
|                 +-------------------------+----------------------------------+
|                 | ident_pct               | *_ident line                     |
|                 +-------------------------+----------------------------------+
|                 | init1_score             | *_init1 line                     |
|                 +-------------------------+----------------------------------+
|                 | initn_score             | *_initn line                     |
|                 +-------------------------+----------------------------------+
|                 | opt_score               | *_opt line, *_s-w opt line       |
|                 +-------------------------+----------------------------------+
|                 | pos_pct                 | *_sim line                       |
|                 +-------------------------+----------------------------------+
|                 | sw_score                | *_score line                     |
|                 +-------------------------+----------------------------------+
|                 | z_score                 | *_z-score line                   |
+-----------------+-------------------------+----------------------------------+
| HSPFragment     | aln_annotation          | al_cons block, if present        |
| (also via HSP)  +-------------------------+----------------------------------+
|                 | hit                     | hit sequence                     |
|                 +-------------------------+----------------------------------+
|                 | hit_end                 | hit sequence end coordinate      |
|                 +-------------------------+----------------------------------+
|                 | hit_start               | hit sequence start coordinate    |
|                 +-------------------------+----------------------------------+
|                 | hit_strand              | hit sequence strand              |
|                 +-------------------------+----------------------------------+
|                 | query                   | query sequence                   |
|                 +-------------------------+----------------------------------+
|                 | query_end               | query sequence end coordinate    |
|                 +-------------------------+----------------------------------+
|                 | query_start             | query sequence start coordinate  |
|                 +-------------------------+----------------------------------+
|                 | query_strand            | query sequence strand            |
+-----------------+-------------------------+----------------------------------+

"""

import re

from Bio._py3k import _as_bytes, _bytes_to_string
from Bio.Alphabet import generic_dna, generic_protein
from Bio.File import UndoHandle
from Bio.SearchIO._index import SearchIndexer
from Bio.SearchIO._model import QueryResult, Hit, HSP, HSPFragment


__all__ = ['FastaM10Parser', 'FastaM10Indexer']


# precompile regex patterns
# regex for program name
_RE_FLAVS = re.compile(r't?fast[afmsxy]|pr[sf][sx]|lalign|[gs]?[glso]search')
# regex for sequence ID and length
_PTR_ID_DESC_SEQLEN = r'>>>(.+?)\s+(.*?) *- (\d+) (?:aa|nt)$'
_RE_ID_DESC_SEQLEN = re.compile(_PTR_ID_DESC_SEQLEN)
_RE_ID_DESC_SEQLEN_IDX = re.compile(_as_bytes(_PTR_ID_DESC_SEQLEN))
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
_STATE_NONE = 0
_STATE_QUERY_BLOCK = 1
_STATE_HIT_BLOCK = 2
_STATE_CONS_BLOCK = 3


def _set_qresult_hits(qresult, hit_rows=[]):
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
            frag = HSPFragment(hit_id, qresult.id)
            hsp = HSP([frag])
            hit = Hit([hsp])
            qresult.append(hit)

    return qresult


def _set_hsp_seqs(hsp, parsed, program):
    """Helper function for the main parsing code.

    Arguments:
    hsp -- HSP object whose properties are to be set.
    parsed -- Dictionary containing parsed values for HSP attributes.
    program -- String of program name.

    """
    # get aligned sequences and check if they have equal lengths
    for seq_type in ('hit', 'query'):
        if 'tfast' not in program:
            parsed[seq_type]['seq'] = _extract_alignment(parsed[seq_type])
    assert len(parsed['query']['seq']) == len(parsed['hit']['seq']), parsed

    # query and hit sequence types must be the same
    assert parsed['query']['_type'] == parsed['hit']['_type']
    type_val = parsed['query']['_type'] # hit works fine too
    alphabet = generic_dna if type_val == 'D' else generic_protein
    setattr(hsp.fragment, 'alphabet', alphabet)

    for seq_type in ('hit', 'query'):
        # get and set start and end coordinates
        start = int(parsed[seq_type]['_start'])
        end = int(parsed[seq_type]['_stop'])

        setattr(hsp.fragment, seq_type + '_start', min(start, end) - 1)
        setattr(hsp.fragment, seq_type + '_end', max(start, end))
        # set seq and alphabet
        setattr(hsp.fragment, seq_type, parsed[seq_type]['seq'])

        if alphabet is not generic_protein:
            # get strand from coordinate; start <= end is plus
            # start > end is minus
            if start <= end:
                setattr(hsp.fragment, seq_type + '_strand', 1)
            else:
                setattr(hsp.fragment, seq_type + '_strand', -1)
        else:
            setattr(hsp.fragment, seq_type + '_strand', 0)


def _extract_alignment(parsed_hsp):
    """Helper function for the main parsing code.

    To get the actual pairwise alignment sequences, we must first
    translate the un-gapped sequence based coordinates into positions
    in the gapped sequence (which may have a flanking region shown
    using leading - characters).  To date, I have never seen any
    trailing flanking region shown in the m10 file, but the
    following code should also cope with that.

    Note that this code seems to work fine even when the "sq_offset"
    entries are prsent as a result of using the -X command line option.
    """
    seq = parsed_hsp['seq']
    seq_stripped = seq.strip('-')
    disp_start = int(parsed_hsp['_display_start'])
    start = int(parsed_hsp['_start'])
    stop = int(parsed_hsp['_stop'])

    if start <= stop:
        start = start - disp_start
        stop = stop - disp_start + 1
    else:
        start = disp_start - start
        stop = disp_start - stop + 1
    stop += seq_stripped.count('-')
    assert 0 <= start and start < stop and stop <= len(seq_stripped), \
           "Problem with sequence start/stop,\n%s[%i:%i]\n%s" \
           % (seq, start, stop, parsed_hsp)
    return seq_stripped[start:stop]


class FastaM10Parser(object):
    """Parser for Bill Pearson's FASTA suite's -m 10 output."""

    def __init__(self, handle, __parse_hit_table=False):
        self.handle = UndoHandle(handle)
        self._preamble = self._parse_preamble()

    def __iter__(self):
        for qresult in self._parse_qresult():
            # re-set desc, for hsp query description
            qresult.description = qresult.description
            yield qresult

    def _parse_preamble(self):
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

    def __parse_hit_table(self):
        """Parses hit table rows."""
        # move to the first row
        self.line = self.handle.readline()
        # parse hit table until we see an empty line
        hit_rows = []
        while self.line and not self.line.strip():
            hit_rows.append(self.line.strip())
            self.line = self.handle.readline()
        return hit_rows

    def _parse_qresult(self):
        # initial qresult value
        qresult = None
        hit_rows = []
        # state values
        state_QRES_NEW = 1
        state_QRES_HITTAB = 3
        state_QRES_CONTENT = 5
        state_QRES_END = 7

        while True:

            # one line before the hit table
            if self.line.startswith('The best scores are:'):
                qres_state = state_QRES_HITTAB
            # the end of a query or the file altogether
            elif self.line.strip() == '>>>///' or not self.line:
                qres_state = state_QRES_END
            # the beginning of a new query
            elif not self.line.startswith('>>>') and '>>>' in self.line:
                qres_state = state_QRES_NEW
            # the beginning of the query info and its hits + hsps
            elif self.line.startswith('>>>') and not \
                    self.line.strip() == '>>><<<':
                qres_state = state_QRES_CONTENT
            # default qres mark
            else:
                qres_state = None

            if qres_state is not None:
                if qres_state == state_QRES_HITTAB:
                    # parse hit table if flag is set
                    hit_rows = self.__parse_hit_table()

                elif qres_state == state_QRES_END:
                    yield _set_qresult_hits(qresult, hit_rows)
                    break

                elif qres_state == state_QRES_NEW:
                    # if qresult is filled, yield it first
                    if qresult is not None:
                        yield _set_qresult_hits(qresult, hit_rows)
                    regx = re.search(_RE_ID_DESC_SEQLEN, self.line)
                    query_id = regx.group(1)
                    seq_len = regx.group(3)
                    desc = regx.group(2)
                    qresult = QueryResult(query_id)
                    qresult.seq_len = int(seq_len)
                    # get target from the next line
                    self.line = self.handle.readline()
                    qresult.target = list(filter(None,
                            self.line.split(' ')))[1].strip()
                    if desc is not None:
                        qresult.description = desc
                    # set values from preamble
                    for key, value in self._preamble.items():
                        setattr(qresult, key, value)

                elif qres_state == state_QRES_CONTENT:
                    assert self.line[3:].startswith(qresult.id), self.line
                    for hit, strand in self._parse_hit(query_id):
                        # re-set desc, for hsp hit description
                        hit.description = hit.description
                        # if hit is not in qresult, append it
                        try:
                            qresult.append(hit)
                        # otherwise, it might be the same hit with a different strand
                        except ValueError:
                            # make sure strand is different and then append hsp to
                            # existing hit
                            for hsp in hit.hsps:
                                assert strand != hsp.query_strand
                                qresult[hit.id].append(hsp)

            self.line = self.handle.readline()

    def _parse_hit(self, query_id):
        while True:
            self.line = self.handle.readline()
            if self.line.startswith('>>'):
                break

        strand = None
        hsp_list = []
        while True:
            peekline = self.handle.peekline()
            # yield hit if we've reached the start of a new query or
            # the end of the search
            if peekline.strip() in [">>><<<", ">>>///"] or \
                    (not peekline.startswith('>>>') and '>>>' in peekline):
                # append last parsed_hsp['hit']['seq'] line
                if state == _STATE_HIT_BLOCK:
                    parsed_hsp['hit']['seq'] += self.line.strip()
                elif state == _STATE_CONS_BLOCK:
                    hsp.aln_annotation['homology'] += \
                            self.line.strip('\n')
                # process HSP alignment and coordinates
                _set_hsp_seqs(hsp, parsed_hsp, self._preamble['program'])
                hit = Hit(hsp_list)
                hit.description = hit_desc
                hit.seq_len = seq_len
                yield hit, strand
                hsp_list = []
                break
            # yield hit and create a new one if we're still in the same query
            elif self.line.startswith('>>'):
                # try yielding,  if we have hsps
                if hsp_list:
                    _set_hsp_seqs(hsp, parsed_hsp, self._preamble['program'])
                    hit = Hit(hsp_list)
                    hit.description = hit_desc
                    hit.seq_len = seq_len
                    yield hit, strand
                    hsp_list = []
                # try to get the hit id and desc, and handle cases without descs
                try:
                    hit_id, hit_desc = self.line[2:].strip().split(' ', 1)
                except ValueError:
                    hit_id = self.line[2:].strip().split(' ', 1)[0]
                    hit_desc = ''
                # create the HSP object for Hit
                frag = HSPFragment(hit_id, query_id)
                hsp = HSP([frag])
                hsp_list.append(hsp)
                # set or reset the state to none
                state = _STATE_NONE
                parsed_hsp = {'query':{}, 'hit': {}}
            # create and append a new HSP if line starts with '>--'
            elif self.line.startswith('>--'):
                # set seq attributes of previous hsp
                _set_hsp_seqs(hsp, parsed_hsp, self._preamble['program'])
                # and create a new one
                frag = HSPFragment(hit_id, query_id)
                hsp = HSP([frag])
                hsp_list.append(hsp)
                # set the state ~ none yet
                state = _STATE_NONE
                parsed_hsp = {'query':{}, 'hit': {}}
            # this is either query or hit data in the HSP, depending on the state
            elif self.line.startswith('>'):
                if state == _STATE_NONE:
                    # make sure it's the correct query
                    assert query_id.startswith(self.line[1:].split(' ')[0]), \
                            "%r vs %r" % (query_id, self.line)
                    state = _STATE_QUERY_BLOCK
                    parsed_hsp['query']['seq'] = ''
                elif state == _STATE_QUERY_BLOCK:
                    # make sure it's the correct hit
                    assert hit_id.startswith(self.line[1:].split(' ')[0])
                    state = _STATE_HIT_BLOCK
                    parsed_hsp['hit']['seq'] = ''
            # check for conservation block
            elif self.line.startswith('; al_cons'):
                state = _STATE_CONS_BLOCK
                hsp.fragment.aln_annotation['homology'] = ''
            elif self.line.startswith(';'):
                # Fasta outputs do not make a clear distinction between Hit
                # and HSPs, so we check the attribute names to determine
                # whether it belongs to a Hit or HSP
                regx = re.search(_RE_ATTR, self.line.strip())
                name = regx.group(1)
                value = regx.group(2)

                # for values before the '>...' query block
                if state == _STATE_NONE:
                    if name in _HSP_ATTR_MAP:
                        attr_name, caster = _HSP_ATTR_MAP[name]
                        if caster is not str:
                            value = caster(value)
                        if name in ['_ident', '_sim']:
                            value *= 100
                        setattr(hsp, attr_name, value)
                # otherwise, pool the values for processing later
                elif state == _STATE_QUERY_BLOCK:
                    parsed_hsp['query'][name] = value
                elif state == _STATE_HIT_BLOCK:
                    if name == '_len':
                        seq_len = int(value)
                    else:
                        parsed_hsp['hit'][name] = value
                # for values in the hit block
                else:
                    raise ValueError("Unexpected line: %r" % self.line)
            # otherwise, it must be lines containing the sequences
            else:
                assert '>' not in self.line
                # if we're in hit, parse into hsp.hit
                if state == _STATE_HIT_BLOCK:
                    parsed_hsp['hit']['seq'] += self.line.strip()
                elif state == _STATE_QUERY_BLOCK:
                    parsed_hsp['query']['seq'] += self.line.strip()
                elif state == _STATE_CONS_BLOCK:
                    hsp.fragment.aln_annotation['homology'] += \
                            self.line.strip('\n')
                # we should not get here!
                else:
                    raise ValueError("Unexpected line: %r" % self.line)

            self.line = self.handle.readline()


class FastaM10Indexer(SearchIndexer):
    """Indexer class for Bill Pearson's FASTA suite's -m 10 output."""

    _parser = FastaM10Parser

    def __init__(self, filename):
        SearchIndexer.__init__(self, filename)
        self._handle = UndoHandle(self._handle)

    def __iter__(self):
        handle = self._handle
        handle.seek(0)
        start_offset = handle.tell()
        qresult_key = None
        query_mark = _as_bytes('>>>')

        while True:
            line = handle.readline()
            peekline = handle.peekline()
            end_offset = handle.tell()

            if not line.startswith(query_mark) and query_mark in line:
                regx = re.search(_RE_ID_DESC_SEQLEN_IDX, line)
                qresult_key = _bytes_to_string(regx.group(1))
                start_offset = end_offset - len(line)
            # yield whenever we encounter a new query or at the end of the file
            if qresult_key is not None:
                if (not peekline.startswith(query_mark)
                        and query_mark in peekline) or not line:
                    yield qresult_key, start_offset, end_offset - start_offset
                    if not line:
                        break
                    start_offset = end_offset

    def get_raw(self, offset):
        handle = self._handle
        qresult_raw = _as_bytes('')
        query_mark = _as_bytes('>>>')

        # read header first
        handle.seek(0)
        while True:
            line = handle.readline()
            peekline = handle.peekline()
            qresult_raw += line
            if not peekline.startswith(query_mark) and query_mark in peekline:
                break

        # and read the qresult raw string
        handle.seek(offset)
        while True:
            # preserve whitespace, don't use read_forward
            line = handle.readline()
            peekline = handle.peekline()
            qresult_raw += line

            # break when we've reached qresult end
            if (not peekline.startswith(query_mark) and query_mark in peekline) or \
                    not line:
                break

        # append mock end marker to qresult_raw, since it's not always present
        return qresult_raw + _as_bytes('>>><<<\n')


# if not used as a module, run the doctest
if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()
