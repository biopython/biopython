# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO parser for Exonerate standard output format."""

# Known issues & gotchas:
# - The cigar parser does not use the extended cigar string; only supports MID
# - Cigar and vulgar parsing results will most likely be different, due to the
#   different type of data stored by both formats

import re

from Bio.SearchIO._objects import QueryResult, Hit, SegmentedHSP
from Bio.SearchIO._index import SearchIndexer


# precompile regex for parsing vulgar and cigar lines
_RE_VULGAR = re.compile(r"""^vulgar:\s+
        (\S+)\s+(\d+)\s+(\d+)\s+([\+-\.])\s+  # query: ID, start, end, strand
        (\S+)\s+(\d+)\s+(\d+)\s+([\+-\.])\s+  # hit: ID, start, end, strand
        (\d+)(\s+.*)$                         # score, vulgar components
        """, re.VERBOSE)

_RE_VCOMP = re.compile(r"""
        \s+(\S+) # vulgar label (C/M: codon/match, G: gap, N: ner, 5/3: splice
                 #               site, I: intron, S: split codon, F: frameshift)
        \s+(\d+) # how many residues to advance in query sequence
        \s+(\d+) # how many residues to advance in hit sequence
        """, re.VERBOSE)

_RE_CIGAR = re.compile(r"""^cigar:\s+
        (\S+)\s+(\d+)\s+(\d+)\s+([\+-\.])\s+  # query: ID, start, end, strand
        (\S+)\s+(\d+)\s+(\d+)\s+([\+-\.])\s+  # hit: ID, start, end, strand
        (\d+)(\s+.*)$                         # score, vulgar components
        """, re.VERBOSE)


# strand char-value mapping
_STRAND_MAP = {'+': 1, '-': -1, '.': 0}


def _absorb_hit(qresult, hit):
    """Absorbs the given hit by merging it with an existing hit."""
    try:
        qresult.append(hit)
    except ValueError:
        assert hit.id in qresult
        for hsp_item in hit:
            qresult[hit.id].append(hsp_item)

    return qresult


class BaseExonerateIterator(object):

    """Abstract iterator for exonerate format."""

    _ALN_MARK = None

    def __init__(self, handle):
        self.handle = handle
        self.line = handle.readline()
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
            yield qresult

    def read_until(self, bool_func):
        """Reads the file handle until the given bool function returns True."""
        while True:
            if not self.line or bool_func(self.line):
                return
            else:
                self.line = self.handle.readline()

    def _parse_hit_or_query_line(line):
        # get id and desc
        try:
            mark, id, desc = line.split(' ', 2)
        except ValueError: # no desc
            mark, id = line.split(' ', 1)
            desc = ''

        return id, desc

    _parse_hit_or_query_line = staticmethod(_parse_hit_or_query_line)

    def parse_alignment_header(self, aln_header, qresult, hit, hsp):
        for line in aln_header:
            # query line
            if line.startswith('Query:'):
                qresult['id'], qresult['desc']  = \
                        BaseExonerateIterator._parse_hit_or_query_line(line)
            # target line
            elif line.startswith('Target:'):
                hit['id'], hit['desc'] = \
                        BaseExonerateIterator._parse_hit_or_query_line(line)
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
        if qresult['desc'].endswith(':[revcomp]'):
            hsp['query_strand'] = '-'
            qresult['desc'] = qresult['desc'].replace(':[revcomp]', '')
        elif qresult['model'].startswith('protein2'):
            hsp['query_strand'] = '.'
        else:
            hsp['query_strand'] = '+'
        if hit['desc'].endswith(':[revcomp]'):
            hsp['hit_strand'] = '-'
            hit['desc'] = hit['desc'].replace(':[revcomp]', '')
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

        while True:
            # if the file has c4 alignments, use that as the alignment mark
            if self.has_c4_alignment:
                self._ALN_MARK = 'C4 Alignment:'

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

            # each line is basically a different HSP, so we always add it to
            # any hit object we have
            hsp = SegmentedHSP(hit_id, qresult_id)
            for attr, value in hsp_parsed.items():
                setattr(hsp, attr, value)
            hit.append(hsp)

            self.line = self.handle.readline()


class ExonerateVulgarIterator(BaseExonerateIterator):

    """Iterator for exonerate vulgar strings."""

    _ALN_MARK = 'vulgar'

    def parse_alignment_block(self, qresult, hit, hsp):
        self.read_until(lambda line: line.startswith('vulgar'))
        vulgars = re.search(_RE_VULGAR, self.line)
        # if the file has c4 alignments
        # check if vulgar values match our previously parsed header values
        if self.has_c4_alignment:
            assert qresult['id'] == vulgars.group(1)
            assert hsp['query_start'] == vulgars.group(2)
            assert hsp['query_end'] == vulgars.group(3)
            assert hsp['query_strand'] == vulgars.group(4)
            assert hit['id'] == vulgars.group(5)
            assert hsp['hit_start'] == vulgars.group(6)
            assert hsp['hit_end'] == vulgars.group(7)
            assert hsp['hit_strand'] == vulgars.group(8)
            assert hsp['score'] == vulgars.group(9)
        else:
            qresult['id'] = vulgars.group(1)
            hsp['query_start'] = vulgars.group(2)
            hsp['query_end'] = vulgars.group(3)
            hsp['query_strand'] = vulgars.group(4)
            hit['id'] = vulgars.group(5)
            hsp['hit_start'] = vulgars.group(6)
            hsp['hit_end'] = vulgars.group(7)
            hsp['hit_strand'] = vulgars.group(8)
            hsp['score'] = vulgars.group(9)

        # cast coords into ints
        hsp['query_start'] = int(hsp['query_start'])
        hsp['query_end'] = int(hsp['query_end'])
        hsp['hit_start'] = int(hsp['hit_start'])
        hsp['hit_end'] = int(hsp['hit_end'])
        # adjust strands
        hsp['hit_strand'] = _STRAND_MAP[hsp['hit_strand']]
        hsp['query_strand'] = _STRAND_MAP[hsp['query_strand']]
        # cast score into int
        hsp['score'] = int(hsp['score'])
        # store vulgar line and parse it
        hsp['vulgar'] = vulgars.group(10)
        hsp = ExonerateVulgarIterator.parse_vulgar(hsp)

        return qresult, hit, hsp

    def parse_vulgar(hsp):
        """Parses the vulgar line present in the hsp dictionary."""
        # containers for block coordinates
        hsp['query_starts'], hsp['query_ends'], \
                hsp['hit_starts'], hsp['hit_ends'] = \
                [hsp['query_start']], [], [hsp['hit_start']], []
        # containers for split codons
        hsp['query_split_codons'], hsp['hit_split_codons'] = [], []
        # containers for introns
        hsp['query_introns'], hsp['hit_introns'] = [], []
        # containers for ner blocks
        hsp['query_ners'], hsp['hit_ners'] = [], []
        # sentinels for tracking query and hit positions
        qpos, hpos = hsp['query_start'], hsp['hit_start']
        # multiplier for determining sentinel movement
        qmove = 1 if hsp['query_strand'] >= 0 else -1
        hmove = 1 if hsp['hit_strand'] >= 0 else -1

        vcomps = re.findall(_RE_VCOMP, hsp['vulgar'])
        for idx, match in enumerate(vcomps):
            label, qstep, hstep = match[0], int(match[1]), int(match[2])
            # check for label, must be recognized
            assert label in 'MCGF53INS', "Unexpected vulgar label: %r" % label
            # match, codon, gaps, or frameshift
            # for now, consider frameshift as gaps
            if label in 'MCGFS':
                # if the previous comp is not an MCGFS block, it's the
                # start of a new block
                if vcomps[idx-1][0] not in 'MCGFS':
                    hsp['query_starts'].append(qpos)
                    hsp['hit_starts'].append(hpos)
            # other labels
            # store the values in the hsp dict as a tuple of (start, stop)
            if label in '53INS':
                # get start and stop from parsed values
                qstart, hstart = qpos, hpos
                qend = qstart + qstep * qmove
                hend = hstart + hstep * hmove
                # adjust the start-stop coords
                sqstart, sqend = min(qstart, qend), max(qstart, qend)
                shstart, shend = min(hstart, hend), max(hstart, hend)
                # then decide which list to store these values into
                # splice sites (5' and 3') are grouped into introns
                # note that here, the 5', 3', and intron coordinates live
                # in separate tuples even though they're part of the same
                # intron. we'll merge them later on after sorting
                if label in '53I':
                    qlist = hsp['query_introns']
                    hlist = hsp['hit_introns']
                # ner blocks
                elif label == 'N':
                    qlist = hsp['query_ners']
                    hlist = hsp['hit_ners']
                # split codons
                # XXX: is it possible to have a frameshift that introduces
                # a codon split? If so, this may need a different treatment..
                elif label == 'S':
                    qlist = hsp['query_split_codons']
                    hlist = hsp['hit_split_codons']
                # and store the values
                qlist.append((sqstart, sqend))
                hlist.append((shstart, shend))

            # move sentinels accordingly
            qpos += qstep * qmove
            hpos += hstep * hmove

            # append to ends if the next comp is not an MCGFS block or
            # if it's the last comp
            if idx == len(vcomps)-1 or \
                    (label in 'MCGFS' and vcomps[idx+1][0] not in 'MCGFS'):
                    hsp['query_ends'].append(qpos)
                    hsp['hit_ends'].append(hpos)

        # adjust coordinates
        for seq_type in ('query_', 'hit_'):
            strand = hsp[seq_type + 'strand']
            # switch coordinates if strand is < 0
            if strand < 0:
                # switch the starts and ends
                hsp[seq_type + 'start'], hsp[seq_type + 'end'] = \
                        hsp[seq_type + 'end'], hsp[seq_type + 'start']
                hsp[seq_type + 'starts'], hsp[seq_type + 'ends'] = \
                        hsp[seq_type + 'ends'], hsp[seq_type + 'starts']

            # merge adjacent 5', 3', and introns into single intron blocks
            introns = []
            for start, end in hsp[seq_type + 'introns']:
                if strand >= 0:
                    if not introns or introns[-1][1] != start:
                        introns.append((start, end))
                    # merge if the end coord of the previous intron is the same
                    # as the current intron
                    elif introns[-1][1] == start:
                        introns[-1] = (introns[-1][0], end)
                else:
                    # merging is slightly different if the strand is -
                    if not introns or introns[-1][0] != end:
                        introns.append((start, end))
                    elif introns[-1][0] == end:
                        introns[-1] = (start, introns[-1][1])
            # set the merged coords back to hsp dict
            hsp[seq_type + 'introns'] = introns

        return hsp

    parse_vulgar = staticmethod(parse_vulgar)


class ExonerateVulgarIndexer(SearchIndexer):

    """Indexer class for exonerate vulgar lines."""

    def __init__(self, *args, **kwargs):
        pass


class ExonerateCigarIterator(BaseExonerateIterator):

    """Iterator for exonerate cigar strings."""

    _ALN_MARK = 'cigar'

    def parse_alignment_block(self, qresult, hit, hsp):
        self.read_until(lambda line: line.startswith('cigar'))
        cigars = re.search(_RE_CIGAR, self.line)
        # if the file has c4 alignments
        # check if cigar values match our previously parsed header values
        if self.has_c4_alignment:
            assert qresult['id'] == cigars.group(1)
            assert hsp['query_start'] == cigars.group(2)
            assert hsp['query_end'] == cigars.group(3)
            assert hsp['query_strand'] == cigars.group(4)
            assert hit['id'] == cigars.group(5)
            assert hsp['hit_start'] == cigars.group(6)
            assert hsp['hit_end'] == cigars.group(7)
            assert hsp['hit_strand'] == cigars.group(8)
            assert hsp['score'] == cigars.group(9)
        else:
            qresult['id'] = cigars.group(1)
            hsp['query_start'] = cigars.group(2)
            hsp['query_end'] = cigars.group(3)
            hsp['query_strand'] = cigars.group(4)
            hit['id'] = cigars.group(5)
            hsp['hit_start'] = cigars.group(6)
            hsp['hit_end'] = cigars.group(7)
            hsp['hit_strand'] = cigars.group(8)
            hsp['score'] = cigars.group(9)

        # cast coords into ints
        hsp['query_start'] = int(hsp['query_start'])
        hsp['query_end'] = int(hsp['query_end'])
        hsp['hit_start'] = int(hsp['hit_start'])
        hsp['hit_end'] = int(hsp['hit_end'])
        # adjust strands
        hsp['query_strand'] = _STRAND_MAP[hsp['query_strand']]
        hsp['hit_strand'] = _STRAND_MAP[hsp['hit_strand']]
        if hsp['query_strand'] < 0:
            hsp['query_start'], hsp['query_end'] = hsp['query_end'], \
                    hsp['query_start']
        if hsp['hit_strand'] < 0:
            hsp['hit_start'], hsp['hit_end'] = hsp['hit_end'], \
                    hsp['hit_start']
        # cast score into int
        hsp['score'] = int(hsp['score'])
        # store cigar line and parse it
        hsp['cigar'] = cigars.group(10)
        # container for block coordinates
        # should be present for SegmentedHSPs
        hsp['query_starts'] = [hsp['query_start']]
        hsp['query_ends'] = [hsp['query_end']]
        hsp['hit_starts'] = [hsp['hit_start']]
        hsp['hit_ends'] = [hsp['hit_end']]

        return qresult, hit, hsp


def _test():
    """Run the Bio.SearchIO.ExonIO module's doctests.

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
