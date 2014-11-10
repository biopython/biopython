# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO parser for BLAST+ tab output format, with or without comments."""

import re

from Bio._py3k import _as_bytes, _bytes_to_string
from Bio._py3k import basestring

from Bio.SearchIO._index import SearchIndexer
from Bio.SearchIO._model import QueryResult, Hit, HSP, HSPFragment


__all__ = ['BlastTabIndexer', 'BlastTabParser', 'BlastTabWriter']

__docformat__ = "restructuredtext en"


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
    'query gi': 'qgi',
    'subject ids': 'sallseqid',
    'subject gi': 'sgi',
    'subject gis': 'sallgi',
    'BTOP': 'btop',
    'subject accs.': 'sallacc',
    'subject tax ids': 'staxids',
    'subject sci names': 'sscinames',
    'subject com names': 'scomnames',
    'subject blast names': 'sblastnames',
    'subject super kingdoms': 'sskingdoms',
    'subject title': 'stitle',
    'subject titles': 'salltitles',
    'subject strand': 'sstrand',
    '% subject coverage': 'qcovs',
    '% hsp coverage': 'qcovhsp',
}

# function to create a list from semicolon-delimited string
# used in BlastTabParser._parse_result_row
_list_semicol = lambda x: x.split(';')
_list_diamond = lambda x: x.split('<>')
# column to class attribute map
_COLUMN_QRESULT = {
    'qseqid': ('id', str),
    'qacc': ('accession', str),
    'qaccver': ('accession_version', str),
    'qlen': ('seq_len', int),
    'qgi': ('gi', str),
}
_COLUMN_HIT = {
    'sseqid': ('id', str),
    'sallseqid': ('id_all', _list_semicol),
    'sacc': ('accession', str),
    'saccver': ('accession_version', str),
    'sallacc': ('accession_all', _list_semicol),
    'sgi': ('gi', str),
    'sallgi': ('gi_all', str),
    'slen': ('seq_len', int),
    'staxids': ('tax_ids', _list_semicol),
    'sscinames': ('sci_names', _list_semicol),
    'scomnames': ('com_names', _list_semicol),
    'sblastnames': ('blast_names', _list_semicol),
    'sskingdoms': ('super_kingdoms', _list_semicol),
    'stitle': ('title', str),
    'salltitles': ('title_all', _list_diamond),
    # set strand as HSP property?
    'sstrand': ('strand', str),
    'qcovs': ('query_coverage', float),
}
_COLUMN_HSP = {
    'bitscore': ('bitscore', float),
    'score': ('bitscore_raw', int),
    'evalue': ('evalue', float),
    'nident': ('ident_num', int),
    'pident': ('ident_pct', float),
    'positive': ('pos_num', int),
    'ppos': ('pos_pct', float),
    'mismatch': ('mismatch_num', int),
    'gaps': ('gap_num', int),
    'gapopen': ('gapopen_num', int),
    'btop': ('btop', str),
    'qcovhsp': ('query_coverage', float),
}
_COLUMN_FRAG = {
    'length': ('aln_span', int),
    'qstart': ('query_start', int),
    'qend': ('query_end', int),
    'sstart': ('hit_start', int),
    'send': ('hit_end', int),
    'qframe': ('query_frame', int),
    'sframe': ('hit_frame', int),
    'frames': ('frames', str),
    'qseq': ('query', str),
    'sseq': ('hit', str),
}
_SUPPORTED_FIELDS = set(list(_COLUMN_QRESULT) + list(_COLUMN_HIT) +
                        list(_COLUMN_HSP) + list(_COLUMN_FRAG))

# column order in the non-commented tabular output variant
# values must be keys inside the column-attribute maps above
_DEFAULT_FIELDS = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch',
        'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
# one field from each of the following sets must exist in order for the
# parser to work
_MIN_QUERY_FIELDS = set(['qseqid', 'qacc', 'qaccver'])
_MIN_HIT_FIELDS = set(['sseqid', 'sacc', 'saccver', 'sallseqid'])

# simple function to create BLAST HSP attributes that may be computed if
# other certain attributes are present
# This was previously implemented in the HSP objects in the old model

_RE_GAPOPEN = re.compile(r'\w-')


def _compute_gapopen_num(hsp):
    """Returns the number of gap openings in the given HSP."""
    gapopen = 0
    for seq_type in ('query', 'hit'):
        seq = str(getattr(hsp, seq_type).seq)
        gapopen += len(re.findall(_RE_GAPOPEN, seq))
    return gapopen


def _augment_blast_hsp(hsp, attr):
    """Calculates the given HSP attribute, for writing."""
    if attr == 'aln_span':
        # aln_span is number of identical matches + mismatches + gaps
        func = lambda hsp: hsp.ident_num + hsp.mismatch_num + hsp.gap_num
    # ident and gap will require the num values to be computed first
    elif attr.startswith('ident'):
        func = lambda hsp: hsp.aln_span - hsp.mismatch_num - hsp.gap_num
    elif attr.startswith('gap'):
        func = lambda hsp: hsp.aln_span - hsp.ident_num - hsp.mismatch_num
    elif attr == 'mismatch_num':
        func = lambda hsp: hsp.aln_span - hsp.ident_num - hsp.gap_num
    elif attr == 'gapopen_num':
        if not hasattr(hsp, 'query') or not hasattr(hsp, 'hit'):
            # mock function so that the except clause below is triggered
            # as both the query and hit are required to compute gapopen
            def mock(hsp):
                raise AttributeError
            func = mock
        else:
            func = _compute_gapopen_num

    # set the num values
    # requires the endswith check, since we only want to set 'num' or 'span'
    # attributes here
    if not hasattr(hsp, attr) and not attr.endswith('_pct'):
        value = func(hsp)
        setattr(hsp, attr, value)

    # if the attr is a percent value, calculate it
    if attr == 'ident_pct':
        func2 = lambda hsp: hsp.ident_num / float(hsp.aln_span) * 100
    elif attr == 'pos_pct':
        func = lambda hsp: hsp.pos_num / float(hsp.aln_span) * 100
    elif attr == 'gap_pct':
        func2 = lambda hsp: hsp.gap_num / float(hsp.aln_span) * 100
    else:
        func2 = None

    # set the pct values
    if func2 is not None:
        value = func2(hsp)
        setattr(hsp, attr, value)


class BlastTabParser(object):

    """Parser for the BLAST tabular format."""

    def __init__(self, handle, comments=False, fields=_DEFAULT_FIELDS):
        self.handle = handle
        self.has_comments = comments
        self.fields = self._prep_fields(fields)
        self.line = self.handle.readline().strip()

    def __iter__(self):
        # stop iteration if file has no lines
        if not self.line:
            raise StopIteration
        # determine which iterator to use
        elif self.has_comments:
            iterfunc = self._parse_commented_qresult
        else:
            iterfunc = self._parse_qresult

        for qresult in iterfunc():
            yield qresult

    def _prep_fields(self, fields):
        """Validates and formats the given fields for use by the parser."""
        # cast into list if fields is a space-separated string
        if isinstance(fields, basestring):
            fields = fields.strip().split(' ')
        # blast allows 'std' as a proxy for the standard default lists
        # we want to transform 'std' to its proper column names
        if 'std' in fields:
            idx = fields.index('std')
            fields = fields[:idx] + _DEFAULT_FIELDS + fields[idx+1:]
        # if set(fields) has a null intersection with minimum required
        # fields for hit and query, raise an exception
        if not set(fields).intersection(_MIN_QUERY_FIELDS) or \
                not set(fields).intersection(_MIN_HIT_FIELDS):
            raise ValueError("Required query and/or hit ID field not found.")

        return fields

    def _parse_commented_qresult(self):
        """Iterator returning `QueryResult` objects from a commented file."""
        while True:
            comments = self._parse_comments()
            if comments:
                try:
                    self.fields = comments['fields']
                    # iterator for the query results
                    qres_iter = self._parse_qresult()
                except KeyError:
                    # no fields means the query has no results
                    assert 'fields' not in comments
                    # create an iterator returning one empty qresult
                    # if the query has no results
                    qres_iter = iter([QueryResult()])

                for qresult in qres_iter:
                    for key, value in comments.items():
                        setattr(qresult, key, value)
                    yield qresult

            else:
                break

    def _parse_comments(self):
        """Returns a dictionary containing tab file comments."""
        comments = {}
        while True:
            # parse program and version
            # example: # BLASTX 2.2.26+
            if 'BLAST' in self.line and 'processed' not in self.line:
                program_line = self.line[len(' #'):].split(' ')
                comments['program'] = program_line[0].lower()
                comments['version'] = program_line[1]
            # parse query id and description (if available)
            # example: # Query: gi|356995852 Mus musculus POU domain
            elif 'Query' in self.line:
                query_line = self.line[len('# Query: '):].split(' ', 1)
                comments['id'] = query_line[0]
                if len(query_line) == 2:
                    comments['description'] = query_line[1]
            # parse target database
            # example: # Database: db/minirefseq_protein
            elif 'Database' in self.line:
                comments['target'] = self.line[len('# Database: '):]
            # parse RID (from remote searches)
            elif 'RID' in self.line:
                comments['rid'] = self.line[len('# RID: '):]
            # parse column order, required for parsing the result lines
            # example: # Fields: query id, query gi, query acc., query length
            elif 'Fields' in self.line:
                comments['fields'] = self._parse_fields_line()
            # if the line has these strings, it's either the end of a comment
            # or the end of a file, so we return all the comments we've parsed
            elif ' hits found' in self.line or 'processed' in self.line:
                self.line = self.handle.readline().strip()
                return comments

            self.line = self.handle.readline()

            if not self.line:
                return comments
            else:
                self.line = self.line.strip()

    def _parse_fields_line(self):
        """Returns a list of column short names from the 'Fields'
        comment line."""
        raw_field_str = self.line[len('# Fields: '):]
        long_fields = raw_field_str.split(', ')
        fields = [_LONG_SHORT_MAP[long_name] for long_name in long_fields]
        return self._prep_fields(fields)

    def _parse_result_row(self):
        """Returns a dictionary of parsed row values."""
        fields = self.fields
        columns = self.line.strip().split('\t')
        assert len(fields) == len(columns), "Expected %i columns, found: " \
            "%i" % (len(fields), len(columns))

        qresult, hit, hsp, frag = {}, {}, {}, {}
        for idx, value in enumerate(columns):
            sname = fields[idx]
            # flag to check if any of the _COLUMNs contain sname
            in_mapping = False
            # iterate over each dict, mapping pair to determine
            # attribute name and value of each column
            for parsed_dict, mapping in (
                    (qresult, _COLUMN_QRESULT),
                    (hit, _COLUMN_HIT),
                    (hsp, _COLUMN_HSP),
                    (frag, _COLUMN_FRAG)):
                # process parsed value according to mapping
                if sname in mapping:
                    attr_name, caster = mapping[sname]
                    if caster is not str:
                        value = caster(value)
                    parsed_dict[attr_name] = value
                    in_mapping = True
            # make sure that any unhandled field is not supported
            if not in_mapping:
                assert sname not in _SUPPORTED_FIELDS

        return {'qresult': qresult, 'hit': hit, 'hsp': hsp, 'frag': frag}

    def _get_id(self, parsed):
        """Returns the value used for a QueryResult or Hit ID from a parsed row."""
        # use 'id', with 'id_all', 'accession' and 'accession_version'
        # fallbacks one of these must have a value since we've checked whether
        # they exist or not when parsing the comments
        id_cache = parsed.get('id')
        if id_cache is None and 'id_all' in parsed:
            id_cache = parsed.get('id_all')[0]
        if id_cache is None:
            id_cache = parsed.get('accession')
        if id_cache is None:
            id_cache = parsed.get('accession_version')

        return id_cache

    def _parse_qresult(self):
        """Generator function that returns QueryResult objects."""
        # state values, used to determine what to do with each line
        state_EOF = 0
        state_QRES_NEW = 1
        state_QRES_SAME = 3
        state_HIT_NEW = 2
        state_HIT_SAME = 4
        # dummies for initial states
        qres_state = None
        hit_state = None
        file_state = None
        # dummies for initial id caches
        prev_qid = None
        prev_hid = None
        # dummies for initial parsed value containers
        cur, prev = None, None
        hit_list, hsp_list = [], []

        while True:
            # store previous line's parsed values if we've past the first line
            if cur is not None:
                prev = cur
                prev_qid = cur_qid
                prev_hid = cur_hid
            # only parse the line if it's not EOF or not a comment line
            if self.line and not self.line.startswith('#'):
                cur = self._parse_result_row()
                cur_qid = self._get_id(cur['qresult'])
                cur_hid = self._get_id(cur['hit'])
            else:
                file_state = state_EOF
                # mock values for cur_qid and cur_hid since the line is empty
                cur_qid, cur_hid = None, None

            # get the state of hit and qresult
            if prev_qid != cur_qid:
                qres_state = state_QRES_NEW
            else:
                qres_state = state_QRES_SAME
            # new hits are hits with different id or hits in a new qresult
            if prev_hid != cur_hid or qres_state == state_QRES_NEW:
                hit_state = state_HIT_NEW
            else:
                hit_state = state_HIT_SAME

            # we're creating objects for the previously parsed line(s),
            # so nothing is done in the first parsed line (prev == None)
            if prev is not None:
                # every line is essentially an HSP with one fragment, so we
                # create both of these for every line
                frag = HSPFragment(prev_hid, prev_qid)
                for attr, value in prev['frag'].items():
                    # adjust coordinates to Python range
                    # NOTE: this requires both start and end coords to be
                    # present, otherwise a KeyError will be raised.
                    # Without this limitation, we might misleadingly set the
                    # start / end coords
                    for seq_type in ('query', 'hit'):
                        if attr == seq_type + '_start':
                            value = min(value,
                                    prev['frag'][seq_type + '_end']) - 1
                        elif attr == seq_type + '_end':
                            value = max(value,
                                    prev['frag'][seq_type + '_start'])
                    setattr(frag, attr, value)
                # strand and frame setattr require the full parsed values
                # to be set first
                for seq_type in ('hit', 'query'):
                    # try to set hit and query frame
                    frame = self._get_frag_frame(frag, seq_type,
                            prev['frag'])
                    setattr(frag, '%s_frame' % seq_type, frame)
                    # try to set hit and query strand
                    strand = self._get_frag_strand(frag, seq_type,
                            prev['frag'])
                    setattr(frag, '%s_strand' % seq_type, strand)

                hsp = HSP([frag])
                for attr, value in prev['hsp'].items():
                    setattr(hsp, attr, value)
                hsp_list.append(hsp)

                # create hit and append to temp hit container if hit_state
                # says we're not at the same hit or at a new query
                if hit_state == state_HIT_NEW:
                    hit = Hit(hsp_list)
                    for attr, value in prev['hit'].items():
                        if attr != 'id_all':
                            setattr(hit, attr, value)
                        else:
                            # not setting hit ID since it's already set from the
                            # prev_hid above
                            setattr(hit, '_id_alt', value[1:])
                    hit_list.append(hit)
                    hsp_list = []
                # create qresult and yield if we're at a new qresult or EOF
                if qres_state == state_QRES_NEW or file_state == state_EOF:
                    qresult = QueryResult(hit_list, prev_qid)
                    for attr, value in prev['qresult'].items():
                        setattr(qresult, attr, value)
                    yield qresult
                    # if current line is EOF, break
                    if file_state == state_EOF:
                        break
                    hit_list = []

            self.line = self.handle.readline().strip()

    def _get_frag_frame(self, frag, seq_type, parsedict):
        """Returns `HSPFragment` frame given the object, its sequence type,
        and its parsed dictionary values."""
        assert seq_type in ('query', 'hit')
        frame = getattr(frag, '%s_frame' % seq_type, None)
        if frame is not None:
            return frame
        else:
            if 'frames' in parsedict:
                # frames is 'x1/x2' string, x1 is query frame, x2 is hit frame
                idx = 0 if seq_type == 'query' else 1
                return int(parsedict['frames'].split('/')[idx])
            # else implicit None return

    def _get_frag_strand(self, frag, seq_type, parsedict):
        """Returns `HSPFragment` strand given the object, its sequence type,
        and its parsed dictionary values."""
        # NOTE: this will never set the strands as 0 for protein
        # queries / hits, since we can't detect the blast flavors
        # from the columns alone.
        assert seq_type in ('query', 'hit')
        strand = getattr(frag, '%s_strand' % seq_type, None)
        if strand is not None:
            return strand
        else:
            # using parsedict instead of the fragment object since
            # we need the unadjusted coordinated values
            start = parsedict.get('%s_start' % seq_type)
            end = parsedict.get('%s_end' % seq_type)
            if start is not None and end is not None:
                return 1 if start <= end else -1
            # else implicit None return


class BlastTabIndexer(SearchIndexer):

    """Indexer class for BLAST+ tab output."""

    _parser = BlastTabParser

    def __init__(self, filename, comments=False, fields=_DEFAULT_FIELDS):
        SearchIndexer.__init__(self, filename, comments=comments, fields=fields)

        # if the file doesn't have comments,
        # get index of column used as the key (qseqid / qacc / qaccver)
        if not self._kwargs['comments']:
            if 'qseqid' in fields:
                self._key_idx = fields.index('qseqid')
            elif 'qacc' in fields:
                self._key_idx = fields.index('qacc')
            elif 'qaccver' in fields:
                self._key_idx = fields.index('qaccver')
            else:
                raise ValueError("Custom fields is missing an ID column. "
                        "One of these must be present: 'qseqid', 'qacc', or 'qaccver'.")

    def __iter__(self):
        """Iterates over the file handle; yields key, start offset, and length."""
        handle = self._handle
        handle.seek(0)

        if not self._kwargs['comments']:
            iterfunc = self._qresult_index
        else:
            iterfunc = self._qresult_index_commented

        for key, offset, length in iterfunc():
            yield _bytes_to_string(key), offset, length

    def _qresult_index_commented(self):
        """Indexer for commented BLAST tabular files."""
        handle = self._handle
        handle.seek(0)
        start_offset = 0
        # mark of a new query
        query_mark = None
        # mark of the query's ID
        qid_mark = _as_bytes('# Query: ')
        # mark of the last line
        end_mark = _as_bytes('# BLAST processed')

        while True:
            end_offset = handle.tell()
            line = handle.readline()

            if query_mark is None:
                query_mark = line
                start_offset = end_offset
            elif line.startswith(qid_mark):
                qresult_key = line[len(qid_mark):].split()[0]
            elif line == query_mark or line.startswith(end_mark):
                yield qresult_key, start_offset, end_offset - start_offset
                start_offset = end_offset
            elif not line:
                break

    def _qresult_index(self):
        """Indexer for noncommented BLAST tabular files."""
        handle = self._handle
        handle.seek(0)
        start_offset = 0
        qresult_key = None
        key_idx = self._key_idx
        tab_char = _as_bytes('\t')

        while True:
            # get end offset here since we only know a qresult ends after
            # encountering the next one
            end_offset = handle.tell()
            # line = handle.readline()
            line = handle.readline()

            if qresult_key is None:
                qresult_key = line.split(tab_char)[key_idx]
            else:
                try:
                    curr_key = line.split(tab_char)[key_idx]
                except IndexError:
                    curr_key = _as_bytes('')

                if curr_key != qresult_key:
                    yield qresult_key, start_offset, end_offset - start_offset
                    qresult_key = curr_key
                    start_offset = end_offset

            # break if we've reached EOF
            if not line:
                break

    def get_raw(self, offset):
        """Returns the raw string of a QueryResult object from the given offset."""
        if self._kwargs['comments']:
            getfunc = self._get_raw_qresult_commented
        else:
            getfunc = self._get_raw_qresult

        return getfunc(offset)

    def _get_raw_qresult(self, offset):
        """Returns the raw string of a single QueryResult from a noncommented file."""
        handle = self._handle
        handle.seek(offset)
        qresult_raw = _as_bytes('')
        tab_char = _as_bytes('\t')
        key_idx = self._key_idx
        qresult_key = None

        while True:
            line = handle.readline()
            # get the key if the first line (qresult key)
            if qresult_key is None:
                qresult_key = line.split(tab_char)[key_idx]
            else:
                try:
                    curr_key = line.split(tab_char)[key_idx]
                except IndexError:
                    curr_key = _as_bytes('')
                # only break when qresult is finished (key is different)
                if curr_key != qresult_key:
                    break
            # append to the raw string as long as qresult is the same
            qresult_raw += line

        return qresult_raw

    def _get_raw_qresult_commented(self, offset):
        """Returns the raw string of a single QueryResult from a commented file."""
        handle = self._handle
        handle.seek(offset)
        qresult_raw = _as_bytes('')
        end_mark = _as_bytes('# BLAST processed')

        # query mark is the line marking a new query
        # something like '# TBLASTN 2.2.25+'
        query_mark = None
        line = handle.readline()
        while line:
            # since query_mark depends on the BLAST search, we need to obtain it
            # first
            if query_mark is None:
                query_mark = line
            # break when we've reached the next qresult or the search ends
            elif line == query_mark or line.startswith(end_mark):
                break

            qresult_raw += line
            line = handle.readline()

        return qresult_raw


class BlastTabWriter(object):

    """Writer for blast-tab output format."""

    def __init__(self, handle, comments=False, fields=_DEFAULT_FIELDS):
        self.handle = handle
        self.has_comments = comments
        self.fields = fields

    def write_file(self, qresults):
        """Writes to the handle, returns how many QueryResult objects are written."""
        handle = self.handle
        qresult_counter, hit_counter, hsp_counter, frag_counter = 0, 0, 0, 0

        for qresult in qresults:
            if self.has_comments:
                handle.write(self._build_comments(qresult))
            if qresult:
                handle.write(self._build_rows(qresult))
                if not self.has_comments:
                    qresult_counter += 1
                hit_counter += len(qresult)
                hsp_counter += sum(len(hit) for hit in qresult)
                frag_counter += sum(len(hit.fragments) for hit in qresult)
            # if it's commented and there are no hits in the qresult, we still
            # increment the counter
            if self.has_comments:
                qresult_counter += 1

        # commented files have a line saying how many queries were processed
        if self.has_comments:
            handle.write('# BLAST processed %i queries' % qresult_counter)

        return qresult_counter, hit_counter, hsp_counter, frag_counter

    def _build_rows(self, qresult):
        """Returns a string containing tabular rows of the QueryResult object."""
        coordinates = set(['qstart', 'qend', 'sstart', 'send'])
        qresult_lines = ''
        for hit in qresult:
            for hsp in hit:
                line = []
                for field in self.fields:
                    # get the column value ~ could either be an attribute
                    # of qresult, hit, or hsp
                    if field in _COLUMN_QRESULT:
                        value = getattr(qresult, _COLUMN_QRESULT[field][0])
                    elif field in _COLUMN_HIT:
                        if field == 'sallseqid':
                            value = getattr(hit, 'id_all')
                        else:
                            value = getattr(hit, _COLUMN_HIT[field][0])
                    # special case, since 'frames' can be determined from
                    # query frame and hit frame
                    elif field == 'frames':
                        value = '%i/%i' % (hsp.query_frame, hsp.hit_frame)
                    elif field in _COLUMN_HSP:
                        try:
                            value = getattr(hsp, _COLUMN_HSP[field][0])
                        except AttributeError:
                            attr = _COLUMN_HSP[field][0]
                            _augment_blast_hsp(hsp, attr)
                            value = getattr(hsp, attr)
                    elif field in _COLUMN_FRAG:
                        value = getattr(hsp, _COLUMN_FRAG[field][0])
                    else:
                        assert field not in _SUPPORTED_FIELDS
                        continue

                    # adjust from and to according to strand, if from and to
                    # is included in the output field
                    if field in coordinates:
                        value = self._adjust_coords(field, value, hsp)
                    # adjust output formatting
                    value = self._adjust_output(field, value)

                    line.append(value)

                hsp_line = '\t'.join(line)
                qresult_lines += hsp_line + '\n'

        return qresult_lines

    def _adjust_coords(self, field, value, hsp):
        """Adjusts start and end coordinates according to strand."""
        assert field in ('qstart', 'qend', 'sstart', 'send')
        # determine sequence type to operate on based on field's first letter
        seq_type = 'query' if field.startswith('q') else 'hit'

        strand = getattr(hsp, '%s_strand' % seq_type, None)
        if strand is None:
            raise ValueError("Required attribute %r not found." %
                    ('%s_strand' % (seq_type)))
        # switch start <--> end coordinates if strand is -1
        if strand < 0:
            if field.endswith('start'):
                value = getattr(hsp, '%s_end' % seq_type)
            elif field.endswith('end'):
                value = getattr(hsp, '%s_start' % seq_type) + 1
        elif field.endswith('start'):
            # adjust start coordinate for positive strand
            value += 1

        return value

    def _adjust_output(self, field, value):
        """Adjusts formatting of the given field and value to mimic native tab output."""
        # qseq and sseq are stored as SeqRecord, but here we only need the str
        if field in ('qseq', 'sseq'):
            value = str(value.seq)

        # evalue formatting, adapted from BLAST+ source:
        # src/objtools/align_format/align_format_util.cpp#L668
        elif field == 'evalue':
            if value < 1.0e-180:
                value = '0.0'
            elif value < 1.0e-99:
                value = '%2.0e' % value
            elif value < 0.0009:
                value = '%3.0e' % value
            elif value < 0.1:
                value = '%4.3f' % value
            elif value < 1.0:
                value = '%3.2f' % value
            elif value < 10.0:
                value = '%2.1f' % value
            else:
                value = '%5.0f' % value

        # pident and ppos formatting
        elif field in ('pident', 'ppos'):
            value = '%.2f' % value

        # evalue formatting, adapted from BLAST+ source:
        # src/objtools/align_format/align_format_util.cpp#L723
        elif field == 'bitscore':
            if value > 9999:
                value = '%4.3e' % value
            elif value > 99.9:
                value = '%4.0d' % value
            else:
                value = '%4.1f' % value

        # coverages have no comma (using floats still ~ a more proper
        # representation)
        elif field in ('qcovhsp', 'qcovs'):
            value = '%.0f' % value

        # list into '<>'-delimited string
        elif field == 'salltitles':
            value = '<>'.join(value)

        # list into ';'-delimited string
        elif field in ('sallseqid', 'sallacc', 'staxids', 'sscinames',
                'scomnames', 'sblastnames', 'sskingdoms'):
            value = ';'.join(value)

        # everything else
        else:
            value = str(value)

        return value

    def _build_comments(self, qres):
        """Returns a string of a QueryResult tabular comment."""
        comments = []
        # inverse mapping of the long-short name map, required
        # for writing comments
        inv_field_map = dict((v, k) for k, v in _LONG_SHORT_MAP.items())

        # try to anticipate qress without version
        if not hasattr(qres, 'version'):
            program_line = '# %s' % qres.program.upper()
        else:
            program_line = '# %s %s' % (qres.program.upper(), qres.version)
        comments.append(program_line)
        # description may or may not be None
        if qres.description is None:
            comments.append('# Query: %s' % qres.id)
        else:
            comments.append('# Query: %s %s' % (qres.id, qres.description))
        # try appending RID line, if present
        try:
            comments.append('# RID: %s' % qres.rid)
        except AttributeError:
            pass
        comments.append('# Database: %s' % qres.target)
        # qresults without hits don't show the Fields comment
        if qres:
            comments.append('# Fields: %s' %
                            ', '.join(inv_field_map[field] for field in self.fields))
        comments.append('# %i hits found' % len(qres))

        return '\n'.join(comments) + '\n'


# if not used as a module, run the doctest
if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()
