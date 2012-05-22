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


def blast_xml_iterator(handle):
    """Generator function to parse BLAST+ XML output as QueryResult objects.

    handle -- Handle to the file, or the filename as string.

    """
    # for more info: http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.mod.dtd
    try:
        from xml.etree import cElementTree as ET
    except ImportError:
        from xml.etree import ElementTree as ET

    def _parse_hit(root_hit_elem):
        """Iterator that transforms Iteration_hits XML elements into Hit objects.

        Arguments:
        root_hit_elem -- Element object of the Iteration_hits tag.

        """
        # Hit level processing
        # Hits are stored in the Iteration_hits tag, with the following
        # DTD
        # <!ELEMENT Hit (
        #        Hit_num,
        #        Hit_id,
        #        Hit_def,
        #        Hit_accession,
        #        Hit_len,
        #        Hit_hsps?)>

        if root_hit_elem is None:
            root_hit_elem = []

        for hit_elem in root_hit_elem:

            # create empty hit object
            hit_id = hit_elem.find('Hit_id').text
            hit = Hit(hit_id, query_id)

            hit.description = hit_elem.find('Hit_def').text
            hit.accession = hit_elem.find('Hit_accession').text
            hit.seq_length = int(hit_elem.find('Hit_len').text)

            for hsp in _parse_hsp(hit_elem.find('Hit_hsps'), hit_id):
                hit.append(hsp)

            # delete element after we finish parsing it
            hit_elem.clear()
            yield hit

    def _parse_hsp(root_hsp_elem, hit_id):
        """Iterator that transforms Hit_hsps XML elements into HSP objects.

        Arguments:
        root_hsp_elem -- Element object of the Hit_hsps tag.
        hit_id -- String of hit ID to initialize all HSP objects.

        """
        # Hit_hsps DTD:
        # <!ELEMENT Hsp (
        #        Hsp_num,
        #        Hsp_bit-score,
        #        Hsp_score,
        #        Hsp_evalue,
        #        Hsp_query-from,
        #        Hsp_query-to,
        #        Hsp_hit-from,
        #        Hsp_hit-to,
        #        Hsp_pattern-from?,
        #        Hsp_pattern-to?,
        #        Hsp_query-frame?,
        #        Hsp_hit-frame?,
        #        Hsp_identity?,
        #        Hsp_positive?,
        #        Hsp_gaps?,
        #        Hsp_align-len?,
        #        Hsp_density?,
        #        Hsp_qseq,
        #        Hsp_hseq,
        #        Hsp_midline?)>

        # if value is None, feed the loop below and empty list
        if root_hsp_elem is None:
            root_hsp_elem = []

        for hsp_elem in root_hsp_elem:
            # get the hit, and query tags first as they are required
            # to initialize the HSP object
            hit_seq = hsp_elem.find('Hsp_hseq').text
            query_seq = hsp_elem.find('Hsp_qseq').text
            hsp = HSP(hit_id, query_id, hit_seq, query_seq)

            hsp.bitscore = float(hsp_elem.find('Hsp_bit-score').text)
            hsp.bitscore_raw = float(hsp_elem.find('Hsp_score').text)
            hsp.evalue = float(hsp_elem.find('Hsp_evalue').text)
            hsp.query_from = int(hsp_elem.find('Hsp_query-from').text)
            hsp.query_to = int(hsp_elem.find('Hsp_query-to').text)
            hsp.hit_from = int(hsp_elem.find('Hsp_hit-from').text)
            hsp.hit_to = int(hsp_elem.find('Hsp_hit-to').text)
            hsp.query_frame = int(hsp_elem.find('Hsp_query-frame').text)
            hsp.hit_frame = int(hsp_elem.find('Hsp_hit-frame').text)
            hsp.identity_num = int(hsp_elem.find('Hsp_identity').text)
            hsp.positive_num = int(hsp_elem.find('Hsp_positive').text)
            hsp.gap_num = int(hsp_elem.find('Hsp_gaps').text)
            hsp.homology = hsp_elem.find('Hsp_midline').text

            # optional attributes
            try:
                hsp.pattern_from = int(hsp_elem.find('Hsp_pattern-from').text)
            except AttributeError:
                hsp.pattern_from = None

            try:
                hsp.pattern_to = int(hsp_elem.find('Hsp_pattern-to').text)
            except AttributeError:
                hsp.pattern_to = None

            try:
                hsp.density = float(hsp_elem.find('Hsp_density').text)
            except AttributeError:
                hsp.density = None

            # delete element after we finish parsing it
            hsp_elem.clear()
            yield hsp

    # dictionary for containing all information prior to the first query
    meta = {}

    # dictionary for mapping tag name and meta key name
    _tag_meta_map = {
        'BlastOutput_program': 'program',
        'BlastOutput_db': 'target',
        'BlastOutput_version': 'program_version',
        'BlastOutput_reference': 'program_reference',
        'BlastOutput_query-ID': 'query_id',
        'BlastOutput_query-def': 'query_description',
        'BlastOutput_query-len': 'query_length',
        'Parameters_matrix': 'param_matrix',
        'Parameters_expect': 'param_evalue_threshold',
        'Parameters_gap-open': 'param_gap_open',
        'Parameters_gap-extend': 'param_gap_extend',
        'Parameters_filter': 'param_filter',
        'Parameters_sc-match': 'param_score_match',
        'Parameters_sc-mismatch': 'param_score_mismatch',
    }

    xml_iter = ET.iterparse(handle)

    # parse the preamble part (anything prior to the first result)
    for event, elem in xml_iter:
        # get the tag values and store them into meta
        if event == 'end' and elem.tag in _tag_meta_map:
            meta_key = _tag_meta_map[elem.tag]

            if meta_key  == 'param_evalue_threshold':
                meta[meta_key] = float(elem.text)
            elif meta_key in ['param_gap_open', 'param_gap_extend', \
                    'param_score_match', 'param_score_mismatch']:
                meta[meta_key] = int(elem.text)
            else:
                meta[meta_key] = elem.text
            # delete element after we finish parsing it
            elem.clear()
            continue
        break

    # we only want the version number, sans the program name
    meta['program_version'] = meta['program_version'].split(' ', 1)[1]

    # parse the queries
    for event, qresult_elem in xml_iter:
        # </Iteration> marks the end of a single query
        # which means we can process it
        if event == 'end' and qresult_elem.tag == 'Iteration':

            # get program and target
            program = meta['program']
            target = meta['target']

            # we'll use the order outlined in the BLAST XML schema
            # <!ELEMENT Iteration (
            #        Iteration_iter-num,
            #        Iteration_query-ID?,
            #        Iteration_query-def?,
            #        Iteration_query-len?,
            #        Iteration_hits?,
            #        Iteration_stat?,
            #        Iteration_message?)>

            # we need the query ID to instantiate the QueryResult object
            query_id = qresult_elem.find('Iteration_query-ID').text

            qresult = QueryResult(query_id, program=program, target=target, meta=meta)
            qresult.description = qresult_elem.find('Iteration_query-def').text
            qresult.query_length = int(qresult_elem.find('Iteration_query-len').text)

            # statistics are stored Iteration_stat's 'grandchildren' with the
            # following DTD
            # <!ELEMENT Statistics (
            #        Statistics_db-num,
            #        Statistics_db-len,
            #        Statistics_hsp-len,
            #        Statistics_eff-space,
            #        Statistics_kappa,
            #        Statistics_lambda,
            #        Statistics_entropy)>

            stat_elem = qresult_elem.find('Iteration_stat').find('Statistics')

            qresult.stat_db_num = int(stat_elem.find('Statistics_db-num').text)
            qresult.stat_db_len = float(stat_elem.find('Statistics_db-len').text)
            qresult.stat_eff_space = float(stat_elem.find('Statistics_eff-space').text)
            qresult.stat_kappa = float(stat_elem.find('Statistics_kappa').text)
            qresult.stat_lambda = float(stat_elem.find('Statistics_lambda').text)
            qresult.stat_entropy = float(stat_elem.find('Statistics_entropy').text)

            for hit in _parse_hit(qresult_elem.find('Iteration_hits')):
                # only append the Hit object if we have HSPs
                if hit:
                    qresult.append(hit)

            # delete element after we finish parsing it
            qresult_elem.clear()
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
