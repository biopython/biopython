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
import warnings

from Bio.SearchIO._objects import QueryResult, Hit, HSP, SearchIndexer


# longname-shortname map
# maps the column names shown in a commented output to its short name
# (the one used in the command line)
_LONG_SHORT_MAP = {
    'query id': 'qseqid',
    'query acc.': 'qacc',
    'query acc.ver': 'qaccver',
    'query length': 'qlen',
    'subject id': 'sseqid',
    'subject acc.': 'sacc',
    'subject acc.ver': 'saccver',
    'subject length': 'slen',
    'alignment length': 'length',
    'bit score': 'bitscore',
    'score': 'score',
    'evalue': 'evalue',
    'identical': 'nident',
    '% identity': 'pident',
    'positives': 'positive',
    '% positives': 'ppos',
    'mismatches': 'mismatch',
    'gaps': 'gaps',
    'q. start': 'qstart',
    'q. end': 'qend',
    's. start': 'sstart',
    's. end': 'send',
    'query frame': 'qframe',
    'sbjct frame': 'sframe',
    'query/sbjct frames': 'frames',
    'query seq': 'qseq',
    'subject seq': 'sseq',
    'gap opens': 'gapopen',
    # unsupported columns
    'query gi': 'qgi',
    'subject ids': 'sallseqid',
    'subject gi': 'sgi',
    'subject gis': 'sallgi',
    'BTOP': 'btop',
}

# column to class attribute map
_COLUMN_QRESULT = {
    'qseqid': 'id',
    'qacc': 'acc',
    'qaccver': 'acc_ver',
    'qlen': 'seq_len',
}
_COLUMN_HIT = {
    'sseqid': 'id',
    'sacc': 'acc',
    'saccver': 'acc_ver',
    'slen': 'seq_len',
}
_COLUMN_HSP = {
    'length': 'init_len',
    'bitscore': 'bitscore',
    'score': 'bitscore_raw',
    'evalue': 'evalue',
    'nident': 'ident_num',
    'pident': 'ident_pct',
    'positive': 'pos_num',
    'ppos': 'pos_pct',
    'mismatch': 'mismatch_num',
    'gaps': 'gap_num',
    'qstart': 'query_from',
    'qend': 'query_to',
    'sstart': 'hit_from',
    'send': 'hit_to',
    'qframe': 'query_frame',
    'sframe': 'hit_frame',
    'frames': 'frames',
    'qseq': 'query',
    'sseq': 'hit',
    'gapopen': 'gapopen_num',
}
_SUPPORTED_FIELDS = _COLUMN_QRESULT.keys() + _COLUMN_HIT.keys() + \
        _COLUMN_HSP.keys()
# ignored columns (for now) are:
# query gi -- qgi
# subject ids --  sallseqid
# subject gi -- sgi
# subject gis -- sallgi
# BTOP -- btop

# column order in the non-commented tabular output variant
# values must be keys inside the column-attribute maps above
_DEFAULT_FIELDS = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', \
        'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
# one field from each of the following sets must exist in order for the
# parser to work
_MIN_QUERY_FIELDS = set(['qseqid', 'qacc', 'qaccver'])
_MIN_HIT_FIELDS = set(['sseqid', 'sacc', 'saccver'])


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


def blast_tab_iterator(handle):
    """Generator function to parse BLAST+ tabular output as QueryResult objects.

    handle -- Handle to the file, or the filename as string.

    This method accepts the tabular output variants with or without headers.
    If the handle points to the tabular variant file with headers, it can
    parse arbitrary tabs. However, is the tabular file does not have any
    headers, then it will raise an Exception if the tab columns are not
    the default ones.

    """
    for qresult in BlastTabIterator(handle):
        yield qresult


class BlastTabIterator(object):

    """Parser for the Blast tabular format."""

    def __init__(self, handle, fields=_DEFAULT_FIELDS, line=None):
        self.handle = handle
        self.fields = fields
        self.cache = {}
        if line is None:
            self.line = self.read_forward()
        else:
            self.line = line

    def __iter__(self):
        # stop iteration if file has no lines
        if self.line is None:
            raise StopIteration
        # if line starts with '#', it's a commented file
        # so we'll parse the comments
        elif self.line.startswith('# '):
            while True:
                # parse comments
                comments = self.parse_comments()
                if not comments:
                    break
                # no fields means the query has no result, so we'll yield
                # and empty query object with the parsed comments
                elif 'fields' not in comments:
                    qresult = QueryResult('')
                    for key, value in comments.items():
                        setattr(qresult, key, value)
                    yield qresult
                # otherwise, we'll use the plain qresult parser to parse
                # the result lines
                else:
                    # set fields according to the parsed comments
                    self.fields = comments['fields']
                    for qresult in self.parse_qresult():
                        for key, value in comments.items():
                            setattr(qresult, key, value)
                        yield qresult
                # if there's no more lines, break
                if not self.line:
                    break
        # otherwise it's a noncommented file
        # and we can parse the result lines directly
        else:
            for qresult in self.parse_qresult():
                yield qresult

    def read_forward(self):
        """Reads through whitespaces, returns the first non-whitespace line."""
        while True:
            line = self.handle.readline()
            # if line has characters and stripping does not remove them,
            # return the line
            if line and line.strip():
                return line.strip()
            # if line ends, return None
            elif not line:
                return

    def parse_comments(self):
        """Returns a dictionary containing tab file comments."""

        comments = {}
        while True:
            # parse program and version
            # example: # BLASTX 2.2.26+
            if 'BLAST' in self.line and 'processed' not in self.line:
                program_line = self.line[len(' #'):].split(' ')
                comments['program'] = program_line[0].lower()
                comments['version'] = program_line[1]
                self.line = self.read_forward()
            # parse query id and description (if available)
            # example: # Query: gi|356995852 Mus musculus POU domain
            elif 'Query' in self.line:
                query_line = self.line[len('# Query: '):].split(' ')
                comments['id'] = query_line[0]
                if len(query_line) > 1:
                    comments['desc'] = ' '.join(query_line[1:])
                self.line = self.read_forward()
            # parse target database
            # example: # Database: db/minirefseq_protein
            elif 'Database' in self.line:
                comments['target'] = self.line[len('# Database: '):]
                self.line = self.read_forward()
            # parse RID (from remote searches)
            elif 'RID' in self.line:
                comments['rid'] = self.line[len('# RID: '):]
                self.line = self.read_forward()
            # parse column order, required for parsing the result lines
            # example: # Fields: query id, query gi, query acc., query length
            elif 'Fields' in self.line:
                raw_field_str = self.line[len('# Fields: '):]
                long_fields = raw_field_str.split(', ')
                fields = [_LONG_SHORT_MAP[long_name] for long_name in long_fields]
                # warn if there are unsupported columns
                for field in fields:
                    if field not in _SUPPORTED_FIELDS:
                        message = "Warning: field '%s' is not yet " \
                                "supported by SearchIO. The data in the " \
                                "corresponding column will be ignored." % field
                        warnings.warn(message)
                # if set(fields) has a null intersection with minimum required
                # fields for hit and query, raise an exception
                if set(fields).isdisjoint(_MIN_QUERY_FIELDS) or \
                        set(fields).isdisjoint(_MIN_HIT_FIELDS):
                    raise ValueError("Required field is not found.")
                comments['fields'] = fields
                self.line = self.read_forward()
            # if the line has these strings, it's either the end of a comment
            # or the end of a file, so we return all the comments we've parsed
            elif ' hits found' in self.line or 'processed' in self.line:
                self.line = self.read_forward()
                return comments
            # otherwise, keep on reading the lines
            else:
                self.line = self.read_forward()

            if not self.line:
                return

    def parse_result_row(self):
        """Returns a dictionary of parsed row values."""
        # returns a dict of assigned var names to level names
        fields = self.fields
        columns = self.line.strip().split('\t')
        assert len(fields) == len(columns), "Expected %i columns, found: " \
            "%i" % (len(fields), len(columns))
        qresult, hit, hsp = {}, {}, {}

        for idx, value in enumerate(columns):
            attr_name = fields[idx]

            if attr_name in _COLUMN_QRESULT:
                qresult[_COLUMN_QRESULT[attr_name]] = value

            elif attr_name in _COLUMN_HIT:
                hit[_COLUMN_HIT[attr_name]] = value

            elif attr_name in _COLUMN_HSP:
                hsp[_COLUMN_HSP[attr_name]] = value
            # make sure that any unhandled field is not supported
            else:
                assert attr_name not in _SUPPORTED_FIELDS

        return {'qresult': qresult, 'hit': hit, 'hsp': hsp}

    def get_id(self, parsed):
        """Returns the value used for a QueryResult or Hit ID from a parsed row."""
        # use 'id', with 'acc' and 'acc_ver' fallbacks
        # one of these must have a value since we've checked whether
        # they exist or not when parsing the comments
        id_cache = parsed.get('id')
        if id_cache is None:
            id_cache = parsed.get('acc')
        if id_cache is None:
            id_cache = parsed.get('acc_ver')

        return id_cache

    def parse_qresult(self):
        """Generator function that returns QueryResult objects.

        Argument:
        fields -- List of field short names or the result lines.

        """
        while True:
            parsed = self.parse_result_row()
            # create qresult object, setattr with parsed values
            qid_cache = self.get_id(parsed['qresult'])

            qresult = QueryResult(qid_cache)
            for qresult_attr in parsed['qresult']:
                setattr(qresult, qresult_attr, parsed['qresult'][qresult_attr])

            # read over the next lines and add the parsed values into qresult
            # as long as the id is the same
            for hit in self.parse_hit_hsp(parsed, qid_cache):
                qresult.append(hit)

            yield qresult

            if not self.line or self.line.startswith('# '):
                break

    def parse_hit_hsp(self, parsed, qid_cache):
        """Generator function that returns Hit objects.

        Argument:
        parsed -- Dictionary of parsed result line.
        qid_cache -- String of QueryResult ID for the iteration.
        fields -- List of field short names or the result lines.

        """
        # return hits while qresult.id is the same as the previous row
        while True:

            # create hit object, setattr with parsed values
            hid_cache = self.get_id(parsed['hit'])

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
                # try to set hit_frame and/or query_frame if frames
                # attribute is set
                if not hasattr(hsp, 'query_frame') and hasattr(hsp, 'frames'):
                    setattr(hsp, 'query_frame', hsp.frames.split('/')[0])
                if not hasattr(hsp, 'hit_frame') and hasattr(hsp, 'frames'):
                    setattr(hsp, 'hit_frame', hsp.frames.split('/')[1])
                hit.append(hsp)

                # read next line and parse it if it exists
                self.line = self.read_forward()
                # if line doesn't exist (file end), break out of loop
                if not self.line or self.line.startswith('#'):
                    break
                else:
                    parsed = self.parse_result_row()
                # if hit.id or qresult.id is different, break out of loop
                if hid_cache != self.get_id(parsed['hit']) or \
                        qid_cache != self.get_id(parsed['qresult']):
                    break

            # append hsp-filled hit into qresult
            yield hit

            # if qresult.id is different compared to the previous line
            # break out of loop
            if not self.line or qid_cache != self.get_id(parsed['qresult']) or \
                    self.line.startswith('#'):
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
