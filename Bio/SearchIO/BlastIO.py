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

import re

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

    def _get_elem_data(elem, name, caster=None):
        """XML element text retriever that casts and handles None."""
        try:
            if caster is not None:
                return caster(elem.find(name).text)
            else:
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

            hit.desc = hit_elem.find('Hit_def').text
            hit.acc = hit_elem.find('Hit_accession').text
            hit.seq_len = int(hit_elem.find('Hit_len').text)

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

            # optional attributes
            hsp.query_frame = _get_elem_data(hsp_elem, 'Hsp_query-frame', int)
            hsp.hit_frame = _get_elem_data(hsp_elem, 'Hsp_hit-frame', int)
            hsp.ident_num = _get_elem_data(hsp_elem, 'Hsp_identity', int)
            hsp.pos_num = _get_elem_data(hsp_elem, 'Hsp_positive', int)
            hsp.gap_num = _get_elem_data(hsp_elem, 'Hsp_gaps', int)
            hsp.homology = _get_elem_data(hsp_elem, 'Hsp_midline')
            hsp.pattern_from = _get_elem_data(hsp_elem, 'Hsp_pattern-from', int)
            hsp.pattern_to = _get_elem_data(hsp_elem, 'Hsp_pattern-to', int)
            hsp.density = _get_elem_data(hsp_elem, 'Hsp_density', float)

            # delete element after we finish parsing it
            hsp_elem.clear()
            yield hsp

    # dictionary for containing all information prior to the first query
    meta = {}
    # dictionary for mapping tag name and meta key name
    _tag_meta_map = {
        'BlastOutput_db': 'target',
        'BlastOutput_program': 'program',
        'BlastOutput_reference': 'program_reference',
        'BlastOutput_version': 'program_version',
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
    _tag_fallback_map = {
        'BlastOutput_query-ID': 'query_id',
        'BlastOutput_query-def': 'query_desc',
        'BlastOutput_query-len': 'query_len',
    }

    # compile RE patterns
    _re_version = re.compile(r'\d+\.\d+\.\d+\+?')

    xml_iter = ET.iterparse(handle)

    # parse the preamble part (anything prior to the first result)
    for event, elem in xml_iter:
        # get the tag values, cast appropriately, store into meta
        if elem.tag in _tag_meta_map:
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
        # capture fallback values
        # these are used only if the first <Iteration> does not have any
        # ID, ref, or len.
        elif elem.tag in _tag_fallback_map:
            fallback_key = _tag_fallback_map[elem.tag]
            _fallback[fallback_key] = elem.text
            elem.clear()
            continue

        if elem.tag == 'BlastOutput_param':
            break

    # we only want the version number, sans the program name or date
    meta['program_version'] = re.search(_re_version, \
            meta['program_version']).group(0)

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

            # get the values required for QueryResult instantiation
            try:
                query_id = qresult_elem.find('Iteration_query-ID').text
            except AttributeError:
                query_id = _fallback['query_id']
            program = meta['program']
            target = meta['target']
            qresult = QueryResult(query_id, program=program, target=target, \
                    meta=meta)

            # assign attributes, with fallbacks
            try:
                description = qresult_elem.find('Iteration_query-def').text
            except AttributeError:
                description = _fallback['query_desc']
            try:
                query_len = int(qresult_elem.find('Iteration_query-len').text)
            except AttributeError:
                query_len = int(_fallback['query_len'])

            qresult.desc = description
            qresult.seq_len = query_len

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

                qresult.stat_db_num = int(stat_elem.find('Statistics_db-num').text)
                qresult.stat_db_len = int(stat_elem.find('Statistics_db-len').text)
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
    # column to class attribute map
    _column_qresult = {
        'query id': 'id',               # qseqid
        'query length': 'seq_len'       # qlen
    }
    _column_hit = {
        'subject id': 'id',             # sseqid
        'subject acc.': 'acc',          # sacc
        'subject length': 'seq_len',    # slen
    }
    _column_hsp = {
        'bit score': 'bitscore',        # bitscore
        'score': 'bitscore_raw',        # score
        'evalue': 'evalue',             # evalue
        'alignment length': 'len',      # length
        'identical': 'ident_num',       # nident
        '% identity': 'ident_pct',      # pident
        'positives': 'pos_num',         # positive
        '% positives': 'pos_pct',       # ppos
        'mismatches': 'mismatch_num',   # mismatch
        'gaps': 'gap_num',              # gaps
        'q. start': 'query_from',       # qstart
        'q. end': 'query_to',           # qend
        's. start': 'hit_from',         # sstart
        's. end': 'hit_to',             # send
        'query frame': 'query_frame',   # qframe
        'subject frame': 'hit_frame',   # sframe
        'query seq': 'query',           # qseq
        'subject seq': 'hit',           # sseq
        'gap opens': 'gap_opens',       # gap opens
    }
    # ignored columns (for now) are:
    # query gi -- qgi
    # query acc -- qacc
    # query acc.ver -- qaccver
    # subject ids --  sallseqid
    # subject gi -- sgi
    # subject gis -- sallgi
    # subject acc.ver -- saccver
    # alignment length -- length
    # query/sbjct frames -- frames
    # BTOP -- btop

    # column order in the non-commented tabular output variant
    # values must be keys inside the column-attribute maps above
    _default_order = ['query id', 'subject id', '% identity', \
            'alignment length', 'mismatches', 'gap opens', 'q. start', \
            'q. end', 's. start', 's. end', 'evalue', 'bit score']

    def _parse_row(line, column_order):
        # returns a dict of assigned var names to level names
        columns = line.strip().split('\t')
        assert len(columns) == len(column_order)
        qresult, hit, hsp = {}, {}, {}

        for idx, value in enumerate(columns):
            attr_name = column_order[idx]

            if attr_name in _column_qresult:
                qresult[_column_qresult[attr_name]] = value

            elif attr_name in _column_hit:
                hit[_column_hit[attr_name]] = value

            elif attr_name in _column_hsp:
                hsp[_column_hsp[attr_name]] = value

            else:
                raise ValueError("Column '%s' not supported in SearchIO." % \
                        attr_name)

        return {'qresult': qresult, 'hit': hit, 'hsp': hsp}

    while True:
        line = handle.readline()
        if line.startswith('#'):
            #column_order = ...
            pass
        else:
            column_order = _default_order
            break

    parsed = _parse_row(line, column_order)
    while True:

        # create qresult object, setattr with parsed values
        qid_cache = parsed['qresult']['id']
        qresult = QueryResult(qid_cache)
        for qresult_attr in parsed['qresult']:
            setattr(qresult, qresult_attr, parsed['qresult'][qresult_attr])

        # append hit to qresult while qresult.id is the same as the
        # previous row
        while True:

            # create hit object, setattr with parsed values
            hid_cache = parsed['hit']['id']
            hit = Hit(hid_cache, qid_cache)
            for hit_attr in parsed['hit']:
                setattr(hit, hit_attr, parsed['hit'][hit_attr])

            # append hsp to hit while hit.id and qresult.id are the same
            # as the previous row
            while True:

                # create hsp object, setattr with parsed values, append to hit
                hsp = HSP(hid_cache, qid_cache)
                for hsp_attr in parsed['hsp']:
                    setattr(hsp, hsp_attr, parsed['hsp'][hsp_attr])
                hit.append(hsp)

                # read next line and parse it if it exists
                line = handle.readline()
                if line:
                    parsed = _parse_row(line, column_order)
                # if line doesn't exist (file end), break out of loop
                else:
                    break
                # if hit.id or qresult.id is different, break out of loop
                if hid_cache != parsed['hit']['id'] or \
                        qid_cache != parsed['qresult']['id']:
                    break

            # append hsp-filled hit into qresult
            qresult.append(hit)
            # if qresult.id is different compared to the previous line
            # break out of loop
            if not line or qid_cache != parsed['qresult']['id']:
                break

        yield qresult

        if not line:
            break


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
