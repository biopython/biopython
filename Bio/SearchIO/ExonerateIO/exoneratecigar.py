# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO parser for Exonerate cigar output format."""

import re

from Bio.SearchIO._index import SearchIndexer

from _base import BaseExonerateIterator, _STRAND_MAP


# precompile regex
_RE_CIGAR = re.compile(r"""^cigar:\s+
        (\S+)\s+(\d+)\s+(\d+)\s+([\+-\.])\s+  # query: ID, start, end, strand
        (\S+)\s+(\d+)\s+(\d+)\s+([\+-\.])\s+  # hit: ID, start, end, strand
        (\d+)(\s+.*)$                         # score, vulgar components
        """, re.VERBOSE)


class ExonerateCigarIterator(BaseExonerateIterator):

    """Iterator for Exonerate cigar strings."""

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


class ExonerateCigarIndexer(SearchIndexer):

    """Indexer class for exonerate cigar lines."""

    def __init__(self, *args, **kwargs):
        pass


def _test():
    """Run the Bio.SearchIO.ExonerateIO.exoneratecigar module's doctests.

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
