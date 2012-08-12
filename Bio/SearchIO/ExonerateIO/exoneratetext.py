# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO parser for Exonerate plain text output format."""

import re
from itertools import chain

from Bio._py3k import _as_bytes, _bytes_to_string

from _base import _BaseExonerateParser, _BaseExonerateIndexer, _STRAND_MAP, \
        _parse_hit_or_query_line
from exoneratevulgar import parse_vulgar_comp, _RE_VULGAR


__all__ = ['ExonerateTextParser', 'ExonerateTextIndexer']


_RE_ALN_ROW = re.compile(r'\s*\d+\s+: (.*) :\s+\d+')
_RE_EXON = re.compile(r'[atgc ]{2,}?[<>]+ \w+ Intron \d+ [<>]+[atgc ]{2,}?')
_RE_EXON_LEN = re.compile(r'(?:(\d+) bp // (\d+) bp)|(?:(\d+) bp)')
_RE_NER = re.compile(r'--<\s+\d+\s+>--')
_RE_NER_LEN = re.compile(r'--<\s+(\d+)\s+>--')
_RE_SCODON_START = re.compile(r'\{(\w{1,2})\}$')
_RE_SCODON_END = re.compile(r'^\{(\w{1,2})\}')


def _flip_codons(codon_seq, target_seq):
    """Flips the codon characters from one seq to another."""
    a, b = '', ''
    for char1, char2 in zip(codon_seq, target_seq):
        # no need to do anything if the codon seq line has nothing
        if char1 == ' ':
            a += char1
            b += char2
        else:
            a += char2
            b += char1

    return a, b


def _get_block_coords(parsed_seq, has_ner=False):
    """Returns a list of start, end coordinates for each given block in the sequence."""
    start = 0
    coords = []
    if not has_ner:
        splitter = _RE_EXON
    else:
        splitter = _RE_NER
    for block in re.split(splitter, parsed_seq):
        start += parsed_seq[start:].find(block)
        end = start + len(block)
        coords.append((start, end))

    return coords


def _get_inter_coords(coords, strand=1):
    """From the given pairs of coordinates, returns a list of pairs
    covering the intervening ranges."""
    # adapted from Python's itertools guide
    # if strand is -1, adjust coords to the ends and starts are chained
    if strand == -1:
        sorted_coords = [(max(a, b), min(a, b)) for a, b in coords]
        inter_coords = list(chain(*sorted_coords))[1:-1]
        return zip(inter_coords[1::2], inter_coords[::2])
    else:
        inter_coords = list(chain(*coords))[1:-1]
        return zip(inter_coords[::2], inter_coords[1::2])


def _stitch_rows(raw_rows):
    """Stitches together the parsed alignment rows and returns them in a list."""
    # deal with possible codon surprise!
    # (i.e. alignments with codons using cdna2genome model)
    # by creating additional rows to contain the codons
    try:
        max_len = max([len(x) for x in raw_rows])
        for row in raw_rows:
            assert len(row) == max_len
    except AssertionError:
        for idx, row in enumerate(raw_rows):
            if len(row) != max_len:
                # codons must be present in the query and hit (so +2)
                assert len(row) + 2 == max_len
                # add additional empty lines to contain codons
                raw_rows[idx] = [' ' * len(row[0])] + row + [' ' * len(row[0])]

    cmbn_rows = []
    for idx, row in enumerate(raw_rows[0]):
        cmbn_row = ''.join(aln_row[idx] for aln_row in raw_rows)
        cmbn_rows.append(cmbn_row)

    return cmbn_rows


def _get_row_idx(row_len):
    """Returns a dictionary of row indices for parsing alignment blocks."""
    idx = {}
    # 3 lines, usually in dna vs dna models
    if row_len == 3:
        idx['query'] = 0
        idx['midline'] = 1
        idx['hit'] = 2
        idx['qannot'], idx['hannot'] = None, None
    # 4 lines, in protein vs dna models
    elif row_len == 4:
        idx['query'] = 0
        idx['midline'] = 1
        idx['hit'] = 2
        idx['hannot'] = 3
        idx['qannot'] = None
    # 5 lines, translated dna vs translated dna
    elif row_len == 5:
        # set sequence indexes
        idx['qannot'] = 0
        idx['query'] = 1
        idx['midline'] = 2
        idx['hit'] = 3
        idx['hannot'] = 4
    else:
        raise ValueError("Unexpected row count in alignment block: "
                "%i" % row_len)
    return idx


def _get_blocks(rows, coords, idx):
    """Returns a list of dictionaries of sequences split by the coordinates."""
    for idx_name in ('query', 'hit', 'midline', 'qannot', 'hannot'):
        assert idx_name in idx
    blocks = []
    for start, end in coords:
        block = {}
        # get seqs according to index
        block['query'] = rows[idx['query']][start:end]
        block['hit'] = rows[idx['hit']][start:end]
        block['homology'] = rows[idx['midline']][start:end]
        if idx['qannot'] is not None:
            block['query_annotation'] = rows[idx['qannot']][start:end]
        if idx['hannot'] is not None:
            block['hit_annotation'] = rows[idx['hannot']][start:end]
        blocks.append(block)

    return blocks


def _fill_coords(hsp, seq_type, inter_lens):
    """Fill the block coordinates of the given hsp dictionary."""

    # manually fill the first coord
    seq_step = 1 if hsp[seq_type + 'strand'] >= 0 else -1
    fstart = hsp[seq_type + 'start']
    fend = fstart + len(\
            hsp[seq_type[:-1]][0].replace('-','').replace('>', \
            '').replace('<', '')) * seq_step
    coords = [(fstart, fend)]
    # and start from the second block, after the first inter seq
    for idx, block in enumerate(hsp[seq_type[:-1]][1:]):
        bstart = coords[-1][1] + inter_lens[idx] * seq_step
        bend = bstart + seq_step * \
                len(block.replace('-', ''))
        coords.append((bstart, bend))

    # adjust the coords so the smallest is [0], if strand is -1
    # couldn't do this above since we need the initial ordering
    if seq_step != 1:
        for i in range(len(coords)):
            coords[i] = coords[i][1], coords[i][0]
    hsp[seq_type + 'ranges'] = coords

    return hsp


class ExonerateTextParser(_BaseExonerateParser):

    """Parser for Exonerate plain text output."""

    _ALN_MARK = 'C4 Alignment:'

    def parse_alignment_block(self, header):
        qresult = header['qresult']
        hit = header['hit']
        hsp = header['hsp']
        # check for values that must have been set by previous methods
        for val_name in ('query_start', 'query_end' ,'hit_start', 'hit_end', \
                'query_strand', 'hit_strand'):
            assert val_name in hsp, hsp

        # get the alignment rows
        raw_aln_blocks, vulgar_comp = self._read_alignment()
        # and stitch them so we have the full sequences in single strings
        # cmbn_rows still has split codon markers
        cmbn_rows = _stitch_rows(raw_aln_blocks)
        row_idx = _get_row_idx(len(cmbn_rows))

        if len(cmbn_rows) == 5:
            # the real aligned sequence is always the 'outer' one, so we want
            # to flip them with their 'inner' pairs
            # flip query sequence
            cmbn_rows[0], cmbn_rows[1] = \
                    _flip_codons(cmbn_rows[0], cmbn_rows[1])
            # flip hit sequence
            cmbn_rows[4], cmbn_rows[3] = \
                    _flip_codons(cmbn_rows[4], cmbn_rows[3])

        # get the sequence blocks
        # we use the query row as reference for block coords
        block_ref = cmbn_rows[row_idx['query']]
        has_ner = 'NER' in qresult['model'].upper()
        seq_coords = _get_block_coords(block_ref, has_ner)
        tmp_seq_blocks = _get_blocks(cmbn_rows, seq_coords, row_idx)

        # get split codon temp coords for later use
        # this result in pairs of base movement for both ends of each row
        scodon_coords = {'query': [], 'hit': []}
        for seq_type in scodon_coords:
            scoords = []
            for block in tmp_seq_blocks:
                m_start = re.search(_RE_SCODON_START, block[seq_type])
                m_end = re.search(_RE_SCODON_END, block[seq_type])
                if m_start:
                    m_start = len(m_start.group(1))
                    scoords.append((m_start, 0))
                else:
                    scoords.append((0, 0))
                if m_end:
                    m_end = len(m_end.group(1))
                    scoords.append((0, m_end))
                else:
                    scoords.append((0, 0))
            # discard first element, assuming that we never start alns with codons
            scodon_coords[seq_type] = scoords

        # remove the split codon markers
        seq_blocks = []
        for seq_block in tmp_seq_blocks:
            for line_name in seq_block:
                seq_block[line_name] = \
                        seq_block[line_name].replace('{', '').replace('}', '')
            seq_blocks.append(seq_block)

        # adjust strands
        hsp['query_strand'] = _STRAND_MAP[hsp['query_strand']]
        hsp['hit_strand'] = _STRAND_MAP[hsp['hit_strand']]
        # cast coords into ints
        hsp['query_start'] = int(hsp['query_start'])
        hsp['query_end'] = int(hsp['query_end'])
        hsp['hit_start'] = int(hsp['hit_start'])
        hsp['hit_end'] = int(hsp['hit_end'])
        # cast score into ints
        hsp['score'] = int(hsp['score'])
        # set sequences
        hsp['query'] = [x['query'] for x in seq_blocks]
        hsp['hit'] = [x['hit'] for x in seq_blocks]
        hsp['alignment_annotation'] = {}
        # get the annotations if they exist
        for annot_type in ('homology', 'query_annotation', 'hit_annotation'):
            try:
                hsp['alignment_annotation'][annot_type] = \
                        [x[annot_type] for x in seq_blocks]
            except KeyError:
                pass

        # use vulgar coordinates if vulgar line is present and return
        if vulgar_comp is not None:
            hsp = parse_vulgar_comp(hsp, vulgar_comp)
            return {'qresult': qresult, 'hit': hit, 'hsp': hsp}

        # otherwise we need to get the coordinates from the alignment
        # get the intervening blocks first, so we can use them
        # to adjust the coordinates
        inter_coords = _get_inter_coords(seq_coords)
        inter_blocks = _get_blocks(cmbn_rows, inter_coords, row_idx)
        # if model is not ner, scan for possible intron lengths, for later use
        if not has_ner:
            # returns a three-component tuple of intron lengths
            # first two component filled == intron in hit and query
            # last component filled == intron in hit or query
            parsed_lens = re.findall(_RE_EXON_LEN, \
                    cmbn_rows[row_idx['midline']])
        # compute start and end coords for each block
        for seq_type in ('query_', 'hit_'):

            # ner blocks and intron blocks require different adjustments
            if not has_ner:
                # set opposite type, for setting introns
                opp_type = 'hit_' if seq_type == 'query_' else 'query_'
                # list of flags to denote if an intron follows a block
                # it reads e.g. this line:
                # "ATGTT{TT}  >>>> Target Intron 1 >>>>  {G}TGTGTGTACATT"
                # and sets the opposing sequence type's intron (since this
                # line is present on the opposite sequence type line)
                has_intron_after = ['Intron' in x[seq_type[:-1]] for x in \
                        inter_blocks]
                assert len(has_intron_after) == len(parsed_lens)
                # create list containing coord adjustments incorporating
                # intron lengths
                inter_lens = []
                for flag, parsed_len in zip(has_intron_after, parsed_lens):
                    if flag:
                        # joint introns
                        if all(parsed_len[:2]):
                            # intron len is [0] if opp_type is query, otherwise
                            # it's [1]
                            intron_len = int(parsed_len[0]) if opp_type == \
                                    'query_' else int(parsed_len[1])
                        # single hit/query introns
                        elif parsed_len[2]:
                            intron_len = int(parsed_len[2])
                        else:
                            raise ValueError("Unexpected intron parsing " \
                                    "result: %r" % parsed_len)
                    else:
                        intron_len = 0

                    inter_lens.append(intron_len)
                # check that inter_lens's length is len opp_type block - 1
                assert len(inter_lens) == len(hsp[opp_type[:-1]])-1, \
                        "%r vs %r" % (len(inter_lens), len(hsp[opp_type[:-1]])-1)
            else:
                inter_lens = [int(x) for x in \
                        re.findall(_RE_NER_LEN, cmbn_rows[row_idx[seq_type[:-1]]])]
                opp_type = seq_type

            # fill the hsp query and hit coordinates
            hsp = _fill_coords(hsp, opp_type, inter_lens)
            strand = 1 if hsp[opp_type + 'strand'] >= 0 else -1
            # and fill the intervening ranges' values
            if not has_ner:
                # set split codon coordinates
                scodons = []
                for idx in range(len(scodon_coords[seq_type[:-1]])):
                    pair = scodon_coords[opp_type[:-1]][idx]
                    if not any(pair):
                        continue
                    else:
                        assert not all(pair)
                    a, b = pair
                    anchor_pair = hsp[opp_type + 'ranges'][idx // 2]
                    if a:
                        func = max if strand == 1 else min
                        anchor = func(anchor_pair)
                        start_c, end_c = anchor + a * strand * -1, anchor
                    elif b:
                        func = min if strand == 1 else max
                        anchor = func(anchor_pair)
                        start_c, end_c = anchor + b * strand, anchor
                    scodons.append((min(start_c, end_c), max(start_c, end_c)))
                hsp[opp_type + 'split_codons'] = scodons

        # now that we've finished parsing coords, we can set the hit and start
        # coord according to Biopython's convention (start <= end)
        for seq_type in ('query_', 'hit_'):
            if hsp[seq_type + 'strand'] == -1:
                n_start = seq_type + 'start'
                n_end = seq_type + 'end'
                hsp[n_start], hsp[n_end] = hsp[n_end], hsp[n_start]

        return {'qresult': qresult, 'hit': hit, 'hsp': hsp}

    def _read_alignment(self):
        """Reads the raw alignment block strings, returns them in a list."""
        raw_aln_blocks = []
        # flag to check whether we're in an aligment row
        in_aln_row = False
        # flag for vulgar line, if present, we can parse coordinates from it
        vulgar_comp = None
        while True:

            match = re.search(_RE_ALN_ROW, self.line.strip())
            # if we have a match, set flags and values
            if match and not in_aln_row:
                start_idx = self.line.index(match.group(1))
                row_len = len(match.group(1))
                in_aln_row = True
                raw_aln_block = []
            # if we're in an alignment row, grab the sequence
            if in_aln_row:
                raw_aln_block.append(self.line[start_idx:start_idx+row_len])
            # reset flags and values if the line matches, we're in an alignment
            # row, and there are more than 1 line in rows
            if match and in_aln_row and len(raw_aln_block) > 1:
                raw_aln_blocks.append(raw_aln_block)
                start_idx = None
                row_len = None
                in_aln_row = False

            self.line = self.handle.readline()
            # try to parse vulgar line if present
            if self.line.startswith('vulgar'):
                vulgar = re.search(_RE_VULGAR, self.line)
                vulgar_comp = vulgar.group(10)
            if not self.line or self.line.startswith(self._ALN_MARK):
                # HACK: this is so that the parse_qresult method does not
                # yield the objects before appending the last HSP. We are doing
                # this to keep the parser compatible with outputs without
                # human-readable alignment outputs. This also relies on the
                # fact that repeated readline() always returns '' on EOF.
                if not self.line:
                    self.line = 'mock'
                break

        return raw_aln_blocks, vulgar_comp


class ExonerateTextIndexer(_BaseExonerateIndexer):

    """Indexer class for Exonerate plain text."""

    _parser = ExonerateTextParser
    _query_mark = _as_bytes('C4 Alignment')

    def get_qresult_id(self, pos):
        """Returns the query ID from the nearest "Query:" line."""
        handle = self._handle
        handle.seek(pos)
        sentinel = _as_bytes('Query:')

        while True:
            line = handle.readline().strip()
            if line.startswith(sentinel):
                break
            if not line:
                raise StopIteration
        qid, desc = _parse_hit_or_query_line(_bytes_to_string(line))

        return qid

    def get_raw(self, offset):
        """Returns the raw string of a QueryResult object from the given offset."""
        handle = self._handle
        handle.seek(offset)
        qresult_key = None
        qresult_raw = _as_bytes('')

        while True:
            line = handle.readline()
            if not line:
                break
            elif line.startswith(self._query_mark):
                cur_pos = handle.tell()
                if qresult_key is None:
                    qresult_key = self.get_qresult_id(cur_pos)
                else:
                    curr_key = self.get_qresult_id(cur_pos)
                    if curr_key != qresult_key:
                        break
                handle.seek(cur_pos)
            qresult_raw += line

        return qresult_raw


def _test():
    """Run the Bio.SearchIO.ExonerateIO.exoneratetext module's doctests.

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
