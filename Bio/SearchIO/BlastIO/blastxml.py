# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO parser for BLAST+ XML output formats."""

import re

from Bio.SearchIO._objects import QueryResult, Hit, HSP
from Bio.SearchIO._index import SearchIndexer


def blast_xml_iterator(handle):
    """Generator function to parse BLAST+ XML output as QueryResult objects.

    handle -- Handle to the file, or the filename as string.

    """
    # for more info: http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.mod.dtd
    try:
        from xml.etree import cElementTree as ET
    except ImportError:
        from xml.etree import ElementTree as ET

    # element - optional qresult attribute name mapping
    _elem_qresult_opt = {
        'Statistics_db-num': 'stat_db_num',
        'Statistics_db-len': 'stat_db_len',
        'Statistics_eff-space': 'stat_eff_space',
        'Statistics_kappa': 'stat_kappa',
        'Statistics_lambda': 'stat_lambda',
        'Statistics_entropy': 'stat_entropy',
    }

    # element - hit attribute name mapping
    _elem_hit = {
        'Hit_def': 'desc',
        'Hit_accession': 'acc',
        'Hit_len': 'seq_len',
    }

    # element - hsp attribute name mapping
    _elem_hsp = {
        'Hsp_bit-score': 'bitscore',
        'Hsp_score': 'bitscore_raw',
        'Hsp_evalue': 'evalue',
        'Hsp_query-from': 'query_from',
        'Hsp_query-to': 'query_to',
        'Hsp_hit-from': 'hit_from',
        'Hsp_hit-to': 'hit_to',
        'Hsp_query-frame': 'query_frame',
        'Hsp_hit-frame': 'hit_frame',
        'Hsp_identity': 'ident_num',
        'Hsp_positive': 'pos_num',
        'Hsp_gaps': 'gap_num',
        'Hsp_align-len': 'init_len',
        #'Hsp_midline': 'homology',
        'Hsp_pattern-from': 'pattern_from',
        'Hsp_pattern-to': 'pattern_to',
        'Hsp_density': 'density',
    }

    def _get_elem_data(elem, name):
        """XML element text retriever that casts and handles None."""
        try:
            return elem.find(name).text
        except AttributeError:
            return None

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

        # feed the loop below an empty list so iteration still works
        if root_hit_elem is None:
            root_hit_elem = []

        for hit_elem in root_hit_elem:

            # create empty hit object
            hit_id = hit_elem.find('Hit_id').text
            hit = Hit(hit_id, query_id)

            for hit_tag in _elem_hit:
                setattr(hit, _elem_hit[hit_tag], \
                        _get_elem_data(hit_elem, hit_tag))

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

            for hsp_tag in _elem_hsp:
                setattr(hsp, _elem_hsp[hsp_tag], \
                        _get_elem_data(hsp_elem, hsp_tag))

            # delete element after we finish parsing it
            hsp_elem.clear()
            yield hsp

    # dictionary for containing all information prior to the first query
    meta = {}
    # dictionary for mapping tag name and meta key name
    _elem_meta = {
        'BlastOutput_db': 'target',
        'BlastOutput_program': 'program',
        'BlastOutput_version': 'version',
        'BlastOutput_reference': 'reference',
        'Parameters_expect': 'param_evalue_threshold',
        'Parameters_entrez-query': 'param_entrez_query',
        'Parameters_filter': 'param_filter',
        'Parameters_gap-extend': 'param_gap_extend',
        'Parameters_gap-open': 'param_gap_open',
        'Parameters_include': 'param_include',
        'Parameters_matrix': 'param_matrix',
        'Parameters_pattern': 'param_pattern',
        'Parameters_sc-match': 'param_score_match',
        'Parameters_sc-mismatch': 'param_score_mismatch',
    }

    # these are fallback tags that store information on the first query
    # outside the <Iteration> tag
    # only used if query_{ID,def,len} is not found in <Iteration>
    # (seen in legacy Blast <2.2.14)
    _fallback = {}
    _elem_fallback_map = {
        'BlastOutput_query-ID': 'query_id',
        'BlastOutput_query-def': 'query_desc',
        'BlastOutput_query-len': 'query_len',
    }

    # compile RE patterns
    _re_version = re.compile(r'\d+\.\d+\.\d+\+?')

    xml_iter = ET.iterparse(handle, events=('start', 'end'))

    # parse the preamble part (anything prior to the first result)
    for event, elem in xml_iter:
        # get the tag values, cast appropriately, store into meta
        if event == 'end' and elem.tag in _elem_meta:
            meta_key = _elem_meta[elem.tag]

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
        # capture fallback values
        # these are used only if the first <Iteration> does not have any
        # ID, ref, or len.
        elif event == 'end' and elem.tag in _elem_fallback_map:
            fallback_key = _elem_fallback_map[elem.tag]
            _fallback[fallback_key] = elem.text
            elem.clear()
            continue

        if event == 'start' and elem.tag == 'BlastOutput_iterations':
            break

    # we only want the version number, sans the program name or date
    meta['version'] = re.search(_re_version, \
            meta['version']).group(0)

    # parse the queries
    for event, qresult_elem in xml_iter:
        # </Iteration> marks the end of a single query
        # which means we can process it
        if event == 'end' and qresult_elem.tag == 'Iteration':

            # we'll use the following schema
            # <!ELEMENT Iteration (
            #        Iteration_iter-num,
            #        Iteration_query-ID?,
            #        Iteration_query-def?,
            #        Iteration_query-len?,
            #        Iteration_hits?,
            #        Iteration_stat?,
            #        Iteration_message?)>

            # assign query attributes with fallbacks
            try:
                query_id = qresult_elem.find('Iteration_query-ID').text
            except AttributeError:
                query_id = _fallback['query_id']
            try:
                description = qresult_elem.find('Iteration_query-def').text
            except AttributeError:
                description = _fallback['query_desc']
            try:
                query_len = qresult_elem.find('Iteration_query-len').text
            except AttributeError:
                query_len = _fallback['query_len']

            # create qresult and assign its attributes
            qresult = QueryResult(query_id)
            qresult.desc = description
            qresult.seq_len = query_len
            for meta_attr in meta:
                setattr(qresult, meta_attr, meta[meta_attr])

            # statistics are stored in Iteration_stat's 'grandchildren' with the
            # following DTD
            # <!ELEMENT Statistics (
            #        Statistics_db-num,
            #        Statistics_db-len,
            #        Statistics_hsp-len,
            #        Statistics_eff-space,
            #        Statistics_kappa,
            #        Statistics_lambda,
            #        Statistics_entropy)>

            stat_iter_elem = qresult_elem.find('Iteration_stat')
            if stat_iter_elem is not None:
                stat_elem = stat_iter_elem.find('Statistics')

                for stat_tag in _elem_qresult_opt:
                    setattr(qresult, _elem_qresult_opt[stat_tag], \
                            _get_elem_data(stat_elem, stat_tag))

            for hit in _parse_hit(qresult_elem.find('Iteration_hits')):
                # only append the Hit object if we have HSPs
                if hit:
                    qresult.append(hit)

            # delete element after we finish parsing it
            qresult_elem.clear()
            yield qresult


class BlastXmlIndexer(SearchIndexer):

    """Indexer class for BLAST+ XML output."""

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
