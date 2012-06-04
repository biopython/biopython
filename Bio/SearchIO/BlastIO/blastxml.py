# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO parser for BLAST+ XML output formats."""
# for more info: http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.mod.dtd

import re
try:
    from xml.etree import cElementTree as ET
except ImportError:
    from xml.etree import ElementTree as ET

from Bio.SearchIO._objects import QueryResult, Hit, HSP
from Bio.SearchIO._index import SearchIndexer


# compile RE patterns
_RE_VERSION = re.compile(r'\d+\.\d+\.\d+\+?')


# element - optional qresult attribute name mapping
_ELEM_QRESULT_OPT = {
    'Statistics_db-num': 'stat_db_num',
    'Statistics_db-len': 'stat_db_len',
    'Statistics_eff-space': 'stat_eff_space',
    'Statistics_kappa': 'stat_kappa',
    'Statistics_lambda': 'stat_lambda',
    'Statistics_entropy': 'stat_entropy',
}

# element - hit attribute name mapping
_ELEM_HIT = {
    'Hit_def': 'desc',
    'Hit_accession': 'acc',
    'Hit_len': 'seq_len',
}

# element - hsp attribute name mapping
_ELEM_HSP = {
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

# dictionary for mapping tag name and meta key name
_ELEM_META = {
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
_ELEM_QRESULTfallback = {
    'BlastOutput_query-ID': 'id',
    'BlastOutput_query-def': 'desc',
    'BlastOutput_query-len': 'len',
}


def blast_xml_iterator(handle):
    """Generator function to parse BLAST+ XML output as QueryResult objects.

    handle -- Handle to the file.

    """
    for qresult in BlastXmlIterator(handle):
        yield qresult


class BlastXmlIterator(object):

    def __init__(self, handle):
        self.xmlhandle = ET.iterparse(handle, events=('start', 'end'))
        self._meta, self._fallback = self.parse_preamble()

    def __iter__(self):
        for qresult in self.parse_qresult():
            yield qresult

    def parse_hit(self, root_hit_elem, query_id):
        """Generator that transforms Iteration_hits XML elements into Hit objects.

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
            hit_id = hit_elem.findtext('Hit_id')
            hit = Hit(hit_id, query_id)

            for hit_tag in _ELEM_HIT:
                setattr(hit, _ELEM_HIT[hit_tag], hit_elem.findtext(hit_tag))

            for hsp in self.parse_hsp(hit_elem.find('Hit_hsps'), query_id, \
                    hit_id):
                hit.append(hsp)

            # delete element after we finish parsing it
            hit_elem.clear()
            yield hit

    def parse_hsp(self, root_hsp_elem, query_id, hit_id):
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
            hit_seq = hsp_elem.findtext('Hsp_hseq')
            query_seq = hsp_elem.findtext('Hsp_qseq')
            hsp = HSP(hit_id, query_id, hit_seq, query_seq)

            for hsp_tag in _ELEM_HSP:
                setattr(hsp, _ELEM_HSP[hsp_tag], hsp_elem.findtext(hsp_tag))

            # delete element after we finish parsing it
            hsp_elem.clear()
            yield hsp

    def parse_preamble(self):
        """Parses all tag data prior to the first query result."""
        # dictionary for containing all information prior to the first query
        meta = {}
        # dictionary for fallback information
        fallback = {}

        # parse the preamble part (anything prior to the first result)
        for event, elem in self.xmlhandle:
            # get the tag values, cast appropriately, store into meta
            if event == 'end' and elem.tag in _ELEM_META:
                meta_key = _ELEM_META[elem.tag]

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
            elif event == 'end' and elem.tag in _ELEM_QRESULTfallback:
                fallback_key = _ELEM_QRESULTfallback[elem.tag]
                fallback[fallback_key] = elem.text
                elem.clear()
                continue

            if event == 'start' and elem.tag == 'BlastOutput_iterations':
                break

        # we only want the version number, sans the program name or date
        meta['version'] = re.search(_RE_VERSION, \
                meta['version']).group(0)

        return meta, fallback

    def parse_qresult(self):
        """Parses query results."""
        # parse the queries
        for event, qresult_elem in self.xmlhandle:
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
                    query_id = self._fallback['id']
                try:
                    description = qresult_elem.find('Iteration_query-def').text
                except AttributeError:
                    description = self._fallback['desc']
                try:
                    query_len = qresult_elem.find('Iteration_query-len').text
                except AttributeError:
                    query_len = self._fallback['len']

                # create qresult and assign its attributes
                qresult = QueryResult(query_id)
                qresult.desc = description
                qresult.seq_len = query_len
                for meta_attr in self._meta:
                    setattr(qresult, meta_attr, self._meta[meta_attr])

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

                    for stat_tag in _ELEM_QRESULT_OPT:
                        setattr(qresult, _ELEM_QRESULT_OPT[stat_tag], \
                                stat_elem.findtext(stat_tag))

                for hit in self.parse_hit(qresult_elem.find('Iteration_hits'), \
                        query_id):
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
