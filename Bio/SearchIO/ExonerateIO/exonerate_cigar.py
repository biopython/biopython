# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO parser for Exonerate cigar output format."""

import re

from Bio._py3k import _as_bytes, _bytes_to_string

from ._base import _BaseExonerateParser, _STRAND_MAP
from .exonerate_vulgar import ExonerateVulgarIndexer


__all__ = ['ExonerateCigarParser', 'ExonerateCigarIndexer']


# precompile regex
_RE_CIGAR = re.compile(r"""^cigar:\s+
        (\S+)\s+(\d+)\s+(\d+)\s+([\+-\.])\s+  # query: ID, start, end, strand
        (\S+)\s+(\d+)\s+(\d+)\s+([\+-\.])\s+  # hit: ID, start, end, strand
        (\d+)(\s+.*)$                         # score, vulgar components
        """, re.VERBOSE)


class ExonerateCigarParser(_BaseExonerateParser):

    """Parser for Exonerate cigar strings."""

    _ALN_MARK = 'cigar'

    def parse_alignment_block(self, header):
        qresult = header['qresult']
        hit = header['hit']
        hsp = header['hsp']
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

        # adjust strands
        hsp['query_strand'] = _STRAND_MAP[hsp['query_strand']]
        hsp['hit_strand'] = _STRAND_MAP[hsp['hit_strand']]
        # cast coords into ints
        qstart = int(hsp['query_start'])
        qend = int(hsp['query_end'])
        hstart = int(hsp['hit_start'])
        hend = int(hsp['hit_end'])
        # set coords (start <= end)
        hsp['query_start'] = min(qstart, qend)
        hsp['query_end'] = max(qstart, qend)
        hsp['hit_start'] = min(hstart, hend)
        hsp['hit_end'] = max(hstart, hend)
        # cast score into int
        hsp['score'] = int(hsp['score'])
        # store cigar components
        hsp['cigar_comp'] = cigars.group(10)
        # HACK: since we can't really figure out exactly when a
        # HSP starts or ends, we set the entire alignment as one HSP
        hsp['query_ranges'] = [(hsp['query_start'], hsp['query_end'])]
        hsp['hit_ranges'] = [(hsp['hit_start'], hsp['hit_end'])]

        return {'qresult': qresult, 'hit': hit, 'hsp': hsp}


class ExonerateCigarIndexer(ExonerateVulgarIndexer):

    """Indexer class for exonerate cigar lines."""

    _parser = ExonerateCigarParser
    _query_mark = _as_bytes('cigar')

    def get_qresult_id(self, pos):
        """Returns the query ID of the nearest cigar line."""
        handle = self._handle
        handle.seek(pos)
        # get line, check if it's a vulgar line, and get query ID
        line = handle.readline()
        assert line.startswith(self._query_mark), line
        id = re.search(_RE_CIGAR, _bytes_to_string(line))
        return id.group(1)


# if not used as a module, run the doctest
if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()
