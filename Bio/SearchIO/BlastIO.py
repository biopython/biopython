# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO parser for BLAST+ output formats.

This module adds support for parsing BLAST+ outputs. BLAST+ is a rewrite
of NCBI's legacy BLAST (Basic Local Alignment Search Tool), based on the
NCBI C++ toolkit. The BLAST+ suite is available as command line programs
or on NCBI's web page.

Specifically, this module supports the following BLAST+ output formats:

  - XML        - 'blast-xml'
  - Tabular    - 'blast-tab'
  - Plain text - 'blast-text'

And the following BLAST+ programs: blastn, blastp, blastx, tblastn, tblastx

More information are available through these links:
  - Publication: http://www.biomedcentral.com/1471-2105/10/421
  - Web interface: http://blast.ncbi.nlm.nih.gov/
  - User guide: http://www.ncbi.nlm.nih.gov/books/NBK1762/

For legacy BLAST outputs, see the Bio.Blast module.

"""

from Bio.SearchIO._objects import Result, Hit, HSP, SearchIndexer


def blast_xml_iterator(handle):
    """Generator function to parse BLAST+ XML output as Result objects.

    handle -- Handle to the file, or the filename as string.

    """
    # TODO: improve parser performance
    from Bio.Blast.NCBIXML import parse
    for query in parse(handle):
        # HACK: query is the first word before any space
        query_id = query.query.split(' ')[0]
        program = query.application
        target = query.database

        # meta information about search, stored prior to any hits in the file
        meta = {
            'program_version': query.version,
            'program_reference': query.reference,
            'parameter_score_match': query.sc_mismatch,
            'parameter_score_mismatch': query.sc_mismatch,
            # only defined in blastp, blastx, tblastx
            'parameter_matrix': query.matrix,
            'parameter_evalue': float(query.expect),
            'parameter_gap_open': query.gap_penalties[0],
            'parameter_gap_extend': query.gap_penalties[1],
        }

        result = Result(query_id, program=program, target=target, meta=meta)
        for hit in query.alignments:
            hit_id = hit.hit_id

            hsps = []
            for hsp in hit.hsps:
                query_seq = hsp.query
                hit_seq = hsp.sbjct
                hsp_obj = HSP(hit_id, query_id, hit_seq, query_seq)

                # set other parsed hsp attributes
                hsp_attrs = {
                    'bitscore': hsp.bits,
                    'bitscore_raw': hsp.score,
                    'evalue': hsp.expect,
                    'frame': hsp.frame,
                    'gap_num': hsp.gaps,
                    'hit_start_idx': hsp.sbjct_start,
                    'hit_stop_idx': hsp.sbjct_end,
                    'homology': hsp.match,
                    'identity_num': hsp.identities,
                    'positive_num': hsp.positives,
                    'query_start_idx': hsp.query_start,
                    'query_stop_idx': hsp.query_end,
                }
                for attr in hsp_attrs:
                    setattr(hsp_obj, attr, hsp_attrs[attr])

                # append hsp to temporary list
                hsps.append(hsp_obj)

            # create hit object with the hsps
            hit_obj = Hit(hit_id, query_id, hsps)

            # set other parsed hit attributes
            hit_attrs = {
                'full_length': hit.length
            }
            for attr in hit_attrs:
                setattr(hit_obj, attr, hit_attrs[attr])

            result.append(hit_obj)

        yield result


def blast_tabular_iterator(handle):
    """Generator function to parse BLAST+ tabular output as Result objects.

    handle -- Handle to the file, or the filename as string.

    This method accepts the tabular output variants with or without headers.
    If the handle points to the tabular variant file with headers, it can
    parse arbitrary tabs. However, is the tabular file does not have any
    headers, then it will raise an Exception if the tab columns are not
    the default ones.

    """



def blast_text_iterator(handle):
    """Generator function to parse BLAST+ plain text output as Result objects.

    handle -- Handle to the file, or the filename as string.

    """



class BlastXmlIndexer(SearchIndexer):

    """Indexer class for BLAST+ XML output."""

    def __init__(self, handle):
        pass



class BlastTabularIndexer(SearchIndexer):

    """Indexer class for BLAST+ tabular output."""

    def __init__(self, handle):
        pass



class BlastTextIndexer(SearchIndexer):

    """Indexer class for BLAST+ plain text output."""

    def __init__(self, handle):
        pass



def _test():
    """Run the Bio.SearchIO.BlastIO module's doctests.

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
