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

from Bio._py3k import _as_bytes, _bytes_to_string
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
_ELEM_QRESULT_FALLBACK = {
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
                value = hsp_elem.findtext(hsp_tag)
                # adjust 'from' and 'to' coordinates to 0-based ones
                if value is not None and ('-from' in hsp_tag or '-to' \
                        in hsp_tag):
                    value = int(value) - 1
                setattr(hsp, _ELEM_HSP[hsp_tag], value)

            # set the homology characters into alignment_annotation dict
            hm_chars = hsp_elem.findtext('Hsp_midline')
            hsp.alignment_annotation = {}
            hsp.alignment_annotation['homology'] = hm_chars

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
            elif event == 'end' and elem.tag in _ELEM_QRESULT_FALLBACK:
                fallback_key = _ELEM_QRESULT_FALLBACK[elem.tag]
                fallback[fallback_key] = elem.text
                elem.clear()
                continue

            if event == 'start' and elem.tag == 'Iteration':
                break

        # we only want the version number, sans the program name or date
        if meta.get('version') is not None:
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
                query_id = qresult_elem.findtext('Iteration_query-ID')
                if query_id is None:
                    query_id = self._fallback['id']

                query_desc = qresult_elem.findtext('Iteration_query-def')
                if query_desc is None:
                    query_desc = self._fallback['desc']

                query_len = qresult_elem.findtext('Iteration_query-len')
                if query_len is None:
                    query_len = self._fallback['len']

                # handle blast searches against databases with Blast's IDs
                # TODO: handle Blast IDs of legacy blast suite?
                if query_id.startswith('Query_'):
                    id_desc = query_desc.split(' ', 1)
                    query_id = id_desc[0]
                    try:
                        query_desc = id_desc[1]
                    except IndexError:
                        query_desc = ''

                # create qresult and assign its attributes
                qresult = QueryResult(query_id)
                qresult.desc = query_desc
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

                        # handle blast searches against databases with Blast's IDs
                        if hit.id.startswith('gnl|BL_ORD_ID|'):
                            real_id = hit.desc.split(' ')[0]
                            # only change the ID if it's not yet present in qresult
                            if real_id not in qresult:
                                id_desc = hit.desc.split(' ', 1)
                                hit.id = id_desc[0]
                                hit.desc = id_desc[1]

                        qresult.append(hit)

                # delete element after we finish parsing it
                qresult_elem.clear()
                yield qresult


class BlastXmlIndexer(SearchIndexer):

    """Indexer class for BLAST+ XML output."""

    def __init__(self, *args, **kwargs):
        SearchIndexer.__init__(self, *args, **kwargs)
        self._parser = blast_xml_iterator
        self.qstart_mark = _as_bytes('<Iteration>')
        self.qend_mark = _as_bytes('</Iteration>')
        self.block_size = 16384
        # TODO: better way to do this?
        iter_obj = BlastXmlIterator(self._handle)
        self._meta, self._fallback = iter_obj._meta, iter_obj._fallback

    def get_offsets(self, string, sub, offset=0):
        idx = string.find(sub, offset)
        while idx >= 0:
            yield idx
            idx = string.find(sub, idx + 1)

    def __iter__(self):
        qstart_mark = self.qstart_mark
        qend_mark = self.qend_mark
        block_size = self.block_size
        handle = self._handle
        handle.seek(0)
        re_desc = re.compile(r'<Iteration_query-ID>(.*?)'
                '</Iteration_query-ID>\s+?<Iteration_query-def>'
                '(.*?)</Iteration_query-def>')
        re_desc_end = re.compile(r'</Iteration_query-def>')
        counter = 0

        while True:
            # read file content (every block size) and store into block
            block = handle.read(block_size)
            # for each iteration start mark
            for qstart_idx in self.get_offsets(block, qstart_mark):
                # get the id of the query (re.search gets the 1st result)
                regx = re.search(re_desc, block[qstart_idx:])
                try:
                    qstart_desc = regx.group(2)
                    qstart_id = regx.group(1)
                # handle cases where the iteration query ID tag lies in a block
                # split
                except AttributeError:
                    # if we still haven't found the query ID, use the
                    # fallback value
                    if re.search(re_desc_end, block[qstart_idx:]):
                        qstart_desc = self._fallback['desc']
                        qstart_id = self._fallback['id']
                    else:
                        # extend the cached read
                        block_ext = block + handle.read(block_size)
                        regx = re.search(re_desc, block_ext[qstart_idx:])
                        qstart_desc = regx.group(2)
                        qstart_id = regx.group(1)
                        # set file pointer to the position before block_ext
                        handle.seek((counter + 1) * block_size)

                # now for getting the length
                # try finding it in the block
                qlen = block[qstart_idx:].find(qend_mark)

                # if not there, loop until query end is found
                block_ext = block
                while qlen < 0:
                    ext = handle.read(block_size)
                    if not ext:
                        raise ValueError("Query end for %r not found" % qstart_id)
                    block_ext += ext
                    qlen = block_ext[qstart_idx:].find(qend_mark)

                # adjust for read blocks
                qstart_idx = counter * block_size + qstart_idx
                qlen = len(qend_mark) + qlen
                # move pointer to original position
                handle.seek((counter + 1) * block_size)
                if qstart_id.startswith('Query_'):
                    qstart_id = qstart_desc.split(' ', 1)[0]
                # yield key, offset, length
                yield qstart_id, qstart_idx, qlen

            counter += 1

            if not block:
                break

    def get_raw(self, offset):
        qend_mark = self.qend_mark
        block_size = self.block_size
        handle = self._handle
        handle.seek(offset)
        counter = 0
        qresult_raw = ''

        while True:
            block = handle.read(block_size)

            # if we reach EOF without encountering any query end mark
            if not block:
                raise ValueError("Query end not found")

            qresult_raw += block
            qend_idx = qresult_raw.find(qend_mark)

            # if a match is found, return the raw qresult string
            if qend_idx > 0:
                return qresult_raw[:qend_idx + len(qend_mark)]
            # otherwise, increment the counter and go on to the next iteration
            counter += 1

    def get(self, offset):
        qresult = SearchIndexer.get(self, offset)
        # set qresult object attributes with values from preamble
        for key, value in self._meta.items():
            setattr(qresult, key, value)
        return qresult


class BlastXmlWriter(object):

    """Writer for blast-xml output format."""

    def __init__(self, handle):
        self.handle = handle

    def write_file(self, qresults):
        pass

    def build_preamble(self):
        pass

    def build_qresults(self):
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
