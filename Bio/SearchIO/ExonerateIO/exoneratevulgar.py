# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO parser for Exonerate vulgar output format."""

import re

from Bio.SearchIO._index import SearchIndexer

from _base import BaseExonerateIterator, _STRAND_MAP


# precompile regex
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


class ExonerateVulgarIterator(BaseExonerateIterator):

    """Iterator for Exonerate vulgar strings."""

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
        hsp['query_scodon_coords'], hsp['hit_scodon_coords'] = [], []
        # containers for introns
        hsp['query_intron_coords'], hsp['hit_intron_coords'] = [], []
        # containers for ner blocks
        hsp['query_ner_coords'], hsp['hit_ner_coords'] = [], []
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
                    qlist = hsp['query_intron_coords']
                    hlist = hsp['hit_intron_coords']
                # ner blocks
                elif label == 'N':
                    qlist = hsp['query_ner_coords']
                    hlist = hsp['hit_ner_coords']
                # split codons
                # XXX: is it possible to have a frameshift that introduces
                # a codon split? If so, this may need a different treatment..
                elif label == 'S':
                    qlist = hsp['query_scodon_coords']
                    hlist = hsp['hit_scodon_coords']
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
            for start, end in hsp[seq_type + 'intron_coords']:
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


def _test():
    """Run the Bio.SearchIO.ExonerateIO.exoneratevulgar module's doctests.

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
