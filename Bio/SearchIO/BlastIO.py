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

from Bio.SearchIO._objects import QueryResult, Hit, HSP, SearchIndexer


def _safe_float(value):
    """Float caster that handles None."""
    if value is not None:
        return float(value)
    return value


def blast_xml_iterator(handle):
    """Generator function to parse BLAST+ XML output as QueryResult objects.

    handle -- Handle to the file, or the filename as string.

    """
    # TODO: improve parser performance
    from Bio.Blast.NCBIXML import parse
    for record in parse(handle):
        # HACK: query is the first word before any space
        parsed_id = record.query.split(' ', 1)
        try:
            query_id, query_description = parsed_id
        except ValueError:
            query_id, query_description = parsed_id[0], None
        program = record.application
        target = record.database

        # meta information about search, stored prior to any hits in the file
        meta = {
            # present in all blast flavors
            'program_version': record.version,
            'program_reference': record.reference,
            'param_score_match': _safe_float(record.sc_match),
            'param_score_mismatch': _safe_float(record.sc_mismatch),
            # only defined in blastp, blastx, tblastx
            'param_matrix': record.matrix,
            'param_evalue_threshold': _safe_float(record.expect),
            'param_filter': record.filter,
            'param_gap_open': _safe_float(record.gap_penalties[0]),
            'param_gap_extend': _safe_float(record.gap_penalties[1]),
        }

        qresult = QueryResult(query_id, program=program, target=target, meta=meta)

        # set attributes of the QueryResult object
        qresult_attrs = {
            'description': query_description,
            'query_length': record.query_length,
            'stat_db_sequences': record.database_sequences,
            'stat_db_length': record.database_length,
            'stat_eff_search_space': record.effective_search_space,
            'stat_kappa': record.ka_params[1],
            'stat_lambda': record.ka_params[0],
            'stat_entropy': record.ka_params[2],
        }
        for attr in qresult_attrs:
            setattr(qresult, attr, qresult_attrs[attr])

        # fill the QueryResult object with Hits
        for hit in record.alignments:
            # HACK: hit id the first word in hit_def before any space
            #parsed_id = hit.hit_def.split(' ', 1)
            #try:
            #    hit_id, hit_description = parsed_id
            #except ValueError:
            #    hit_id, hit_description = parsed_id[0], None
            hit_id = hit.hit_id
            hit_description = hit.hit_def

            hsps = []
            for hsp in hit.hsps:
                query_seq = hsp.query
                hit_seq = hsp.sbjct
                hsp_obj = HSP(hit_id, query_id, hit_seq, query_seq)

                # set attributes of the HSP object
                hsp_attrs = {
                    'bitscore': hsp.bits,
                    'bitscore_raw': hsp.score,
                    'evalue': hsp.expect,
                    'gap_num': hsp.gaps,
                    'hit_frame': hsp.frame[1],
                    'hit_from': hsp.sbjct_start,
                    'hit_to': hsp.sbjct_end,
                    'homology': hsp.match,
                    'identity_num': hsp.identities,
                    'positive_num': hsp.positives,
                    'query_frame': hsp.frame[0],
                    'query_from': hsp.query_start,
                    'query_to': hsp.query_end,
                }
                for attr in hsp_attrs:
                    setattr(hsp_obj, attr, hsp_attrs[attr])

                # append hsp to temporary list
                hsps.append(hsp_obj)

            # create hit object with the hsps
            hit_obj = Hit(hit_id, query_id, hsps)

            # set attributes of the Hit object
            hit_attrs = {
                'description': hit_description,
                'seq_length': hit.length,
            }
            for attr in hit_attrs:
                setattr(hit_obj, attr, hit_attrs[attr])

            qresult.append(hit_obj)

        yield qresult


def blast_tabular_iterator(handle):
    """Generator function to parse BLAST+ tabular output as QueryResult objects.

    handle -- Handle to the file, or the filename as string.

    This method accepts the tabular output variants with or without headers.
    If the handle points to the tabular variant file with headers, it can
    parse arbitrary tabs. However, is the tabular file does not have any
    headers, then it will raise an Exception if the tab columns are not
    the default ones.

    """



def blast_text_iterator(handle):
    """Generator function to parse BLAST+ plain text output as QueryResult objects.

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
