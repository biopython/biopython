# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO parser for BLAST+ tab output format, with and without comments."""

import warnings

from Bio._py3k import _as_bytes, _bytes_to_string
from Bio.SearchIO._objects import QueryResult, Hit, HSP
from Bio.SearchIO._index import SearchIndexer


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
    # unsupported columns
    'BTOP': 'btop',
}

# column to class attribute map
_COLUMN_QRESULT = {
    'qseqid': 'id',
    'qacc': 'acc',
    'qaccver': 'acc_ver',
    'qlen': 'seq_len',
    'qgi': 'gi',
}
_COLUMN_HIT = {
    'sseqid': 'id',
    'sallseqid': 'id_all',
    'sacc': 'acc',
    'saccver': 'acc_ver',
    'sgi': 'gi',
    'sallgi': 'gi_all',
    'slen': 'seq_len',
}
_COLUMN_HSP = {
    'length': 'ali_len',
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
_SUPPORTED_FIELDS = set(_COLUMN_QRESULT.keys() + _COLUMN_HIT.keys() + \
        _COLUMN_HSP.keys())
# ignored columns (for now) are:
# BTOP -- btop

# column order in the non-commented tabular output variant
# values must be keys inside the column-attribute maps above
_DEFAULT_FIELDS = ('qseqid', 'sseqid', 'pident', 'length', 'mismatch', \
        'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
# one field from each of the following sets must exist in order for the
# parser to work
_MIN_QUERY_FIELDS = set(['qseqid', 'qacc', 'qaccver'])
_MIN_HIT_FIELDS = set(['sseqid', 'sacc', 'saccver'])


def read_forward(handle, strip=True):
    """Reads through whitespaces, returns the first non-whitespace line."""
    while True:
        line = handle.readline()
        # return the line if it has characters
        if line and line.strip():
            if strip:
                return line.strip()
            else:
                return line
        # or if has no characters (EOF)
        elif not line:
            return line


class BlastTabIterator(object):

    """Parser for the Blast tabular format."""

    def __init__(self, handle, fields=_DEFAULT_FIELDS):
        self.handle = handle
        self.fields = fields
        self.line = read_forward(self.handle)

    def __iter__(self):
        # stop iteration if file has no lines
        if not self.line:
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
            # parse query id and description (if available)
            # example: # Query: gi|356995852 Mus musculus POU domain
            elif 'Query' in self.line:
                query_line = self.line[len('# Query: '):].split(' ')
                comments['id'] = query_line[0]
                if len(query_line) > 1:
                    comments['desc'] = ' '.join(query_line[1:])
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
                if not set(fields).intersection(_MIN_QUERY_FIELDS) or \
                        not set(fields).intersection(_MIN_HIT_FIELDS):
                    raise ValueError("Required field is not found.")
                comments['fields'] = fields
            # if the line has these strings, it's either the end of a comment
            # or the end of a file, so we return all the comments we've parsed
            elif ' hits found' in self.line or 'processed' in self.line:
                self.line = read_forward(self.handle)
                return comments

            self.line = read_forward(self.handle)

            if not self.line:
                return

    def parse_result_row(self):
        """Returns a dictionary of parsed row values."""
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
                # adjust 'from' and 'to' coordinates to 0-based ones
                if attr_name in ['qstart', 'qend', 'sstart', 'send']:
                    value = int(value) - 1
                hsp[_COLUMN_HSP[attr_name]] = value
            # make sure that any unhandled field is not supported
            else:
                assert attr_name not in _SUPPORTED_FIELDS

        return qresult, hit, hsp

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
        """Generator function that returns QueryResult objects."""
        qid_cache = ''
        hid_cache = ''
        same_query = False
        while True:
            # only parse the line if it's not EOF or not a comment line
            if self.line and not self.line.startswith('#' ):
                qres_parsed, hit_parsed, hsp_parsed = self.parse_result_row()
                qresult_id = self.get_id(qres_parsed)

            # create new qresult whenever parsed qresult id != cached qresult id
            if qid_cache != qresult_id:
                # if cache is filled (qresult # >1), we can append the latest
                # hit and yield
                if qid_cache:
                    qresult.append(hit)
                    yield qresult
                    same_query = False
                qid_cache = qresult_id
                qresult = QueryResult(qid_cache)
                for attr, value in qres_parsed.items():
                    setattr(qresult, attr, value)
            # append remaining hit and yield if it's EOF or a comment line
            elif not self.line or self.line.startswith('#' ):
                qresult.append(hit)
                yield qresult
                break
            # otherwise, we're still in the same qresult, so set the flag
            elif not same_query:
                same_query = True

            hit_id = self.get_id(hit_parsed)
            # create new hit whenever parsed hit id != cached hit id
            if hid_cache != hit_id:
                # if the qresult id is the same as the previous line,
                # append the previous line's hit before creating this one's
                if same_query:
                    qresult.append(hit)
                hid_cache = hit_id
                hit = Hit(hit_id, qresult_id)
                for attr, value in hit_parsed.items():
                    setattr(hit, attr, value)

            # every line is essentially an HSP, so we create a new one
            # for every line
            hsp = HSP(hid_cache, qid_cache)
            for attr, value in hsp_parsed.items():
                setattr(hsp, attr, value)
            # try to set hit_frame and/or query_frame if frames
            # attribute is set
            if not hasattr(hsp, 'query_frame') and hasattr(hsp, 'frames'):
                setattr(hsp, 'query_frame', hsp.frames.split('/')[0])
            if not hasattr(hsp, 'hit_frame') and hasattr(hsp, 'frames'):
                setattr(hsp, 'hit_frame', hsp.frames.split('/')[1])
            hit.append(hsp)

            self.line = read_forward(self.handle)


class BlastTabIndexer(SearchIndexer):

    """Indexer class for BLAST+ tab output."""

    _parser = BlastTabIterator

    def __init__(self, *args, **kwargs):
        SearchIndexer.__init__(self, *args, **kwargs)
        # set parser for on-the-fly parsing
        # find out if file is commented or not first
        line = read_forward(self._handle)
        if line.startswith('# '):
            self.is_commented = True
        else:
            self.is_commented = False
        # and reset handle
        self._handle.seek(0)

    def __iter__(self):
        """Iterates over the file handle; yields key, start offset, and length."""
        handle = self._handle
        handle.seek(0)
        start_offset = handle.tell()

        if not self.is_commented:
            tab_char = _as_bytes('\t')
            qresult_key = None
            while True:
                # get end offset here since we only know a qresult ends after
                # encountering the next one
                end_offset = handle.tell()
                #line = handle.readline()
                line = read_forward(handle)

                if qresult_key is None:
                    qresult_key = line.split(tab_char)[0]
                else:
                    curr_key = line.split(tab_char)[0]

                    if curr_key != qresult_key:
                        yield _bytes_to_string(qresult_key), start_offset, \
                                end_offset - start_offset
                        qresult_key = curr_key
                        start_offset = end_offset

                # break if we've reached EOF
                if not line:
                    break
        else:
            # mark of a new query
            query_mark = None
            # mark of the query's ID
            qid_mark = '# Query: '

            while True:
                end_offset = handle.tell()
                line = read_forward(handle)

                if query_mark is None:
                    query_mark = line
                    start_offset = end_offset
                elif line.startswith(qid_mark):
                    qresult_key = line[len(qid_mark):].split()[0]
                elif line == query_mark or 'BLAST processed' in line:
                    yield _bytes_to_string(qresult_key), start_offset, \
                            end_offset - start_offset
                    start_offset = end_offset
                elif not line:
                    break

    def get_raw(self, offset):
        """Returns the raw string of a QueryResult object from the given offset."""
        handle = self._handle
        handle.seek(offset)
        qresult_raw = ''

        if not self.is_commented:
            tab_char = _as_bytes('\t')
            qresult_key = None
            while True:
                line = read_forward(handle, strip=False)
                # get the key if the first line (qresult key)
                if qresult_key is None:
                    qresult_key = line.split(tab_char)[0]
                else:
                    curr_key = line.split(tab_char)[0]
                    # only break when qresult is finished (key is different)
                    if curr_key != qresult_key:
                        break
                # append to the raw string as long as qresult is the same
                qresult_raw += line
        else:
            # query mark is the line marking a new query
            # something like '# TBLASTN 2.2.25+'
            query_mark = None
            while True:
                line = read_forward(handle, strip=False)
                if query_mark is None:
                    query_mark = line
                # if we've encountered another query mark, it's the start of
                # another query
                # if 'BLAST processed' is in line, it's one line before EOF
                elif line == query_mark:
                    break
                # append to the raw string as long as qresult is the same
                qresult_raw += line

                if 'BLAST processed' in line:
                    break

        return qresult_raw


class BlastTabWriter(object):

    """Writer for blast-tab output format."""

    def __init__(self, handle, fields=_DEFAULT_FIELDS):
        self.handle = handle
        self.fields = fields

    def write_file(self, qresults):
        """Writes to the handle, returns how many QueryResult objects are written."""
        handle = self.handle
        qresult_counter, hit_counter, hsp_counter = 0, 0, 0

        for qresult in qresults:
            if qresult:
                handle.write(self.build_rows(qresult))
                qresult_counter += 1
                hit_counter += len(qresult)
                hsp_counter += sum([len(hit) for hit in qresult])

        return qresult_counter, hit_counter, hsp_counter

    def build_rows(self, qresult):
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
                        value = getattr(qresult, _COLUMN_QRESULT[field])
                    elif field in _COLUMN_HIT:
                        value = getattr(hit, _COLUMN_HIT[field])
                    # special case, since 'frames' can be determined from
                    # query frame and hit frame
                    elif field == 'frames':
                        value = '%i/%i' % (hsp.query_frame, hsp.hit_frame)
                    elif field in _COLUMN_HSP:
                        value = getattr(hsp, _COLUMN_HSP[field])
                    else:
                        assert field not in _SUPPORTED_FIELDS
                        continue

                    # adjust from and to according to strand, if from and to
                    # is included in the output field
                    if field in coordinates:
                        value = self.adjust_fromto(field, value, hsp)
                    # adjust output formatting
                    value = self.adjust_output(field, value)

                    line.append(value)

                hsp_line = '\t'.join(line)
                qresult_lines += hsp_line + '\n'

        return qresult_lines

    def adjust_fromto(self, field, value, hsp):
        """Adjusts 'from' and 'to' properties according to strand."""
        # try to determine whether strand is minus or not
        # TODO: is there a better way to do this without accessing the private
        # attributes?
        if field in ('qstart', 'qend'):
            try:
                qstrand_is_minus = hsp.query_strand < 0
            except AttributeError:
                qstrand_is_minus = hsp._query_to < hsp._query_from
            # switch from <--> to if strand is -1
            if qstrand_is_minus:
                if field == 'qstart':
                    value = hsp.query_to
                elif field == 'qend':
                    value = hsp.query_from
                else:
                   # we should not get here!
                   raise ValueError("Unexpected column name: %r" % field)
        else:
            try:
                hstrand_is_minus = hsp.hit_strand < 0
            except AttributeError:
                hstrand_is_minus = hsp._hit_to < hsp._hit_from
            # switch from <--> to if strand is -1
            if hstrand_is_minus:
                if field == 'sstart':
                    value = hsp.hit_to
                elif field == 'send':
                   value = hsp.hit_from
                else:
                   # we should not get here!
                   raise ValueError("Unexpected column name: %r" % field)

        # adjust from 0-based index to 1-based
        return value + 1

    def adjust_output(self, field, value):
        """Adjusts formatting of the given field and value to mimic native tab output."""

        # evalue formatting, adapted from BLAST+ source:
        # src/objtools/align_format/align_format_util.cpp#L668
        if field == 'evalue':
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

        # everything else
        else:
            value = str(value)

        return value


class BlastTabcWriter(BlastTabWriter):

    """Writer for blast-tabc output format."""

    def write_file(self, qresults):
        """Writes to the handle, returns how many QueryResult objects are written."""
        handle = self.handle
        qresult_counter, hit_counter, hsp_counter = 0, 0, 0

        for qresult in qresults:
            handle.write(self.build_comments(qresult))
            if qresult:
                handle.write(self.build_rows(qresult))
                hit_counter += len(qresult)
                hsp_counter += sum([len(hit) for hit in qresult])
            # if it's commented and there are no hits in the qresult, we still
            # increment the counter
            qresult_counter += 1

        # commented files have a line saying how many queries were processed
        handle.write('# BLAST processed %i queries' % qresult_counter)

        return qresult_counter, hit_counter, hsp_counter

    def build_comments(self, qresult):
        """Returns a string of a QueryResult tabular comment."""
        comments = ''
        # inverse mapping of the long-short name map, required
        # for writing comments
        inv_field_map = dict((value, key) for key, value in \
                _LONG_SHORT_MAP.items())

        # try to anticipate qresults without version
        if not hasattr(qresult, 'version'):
            program_line = '# %s\n' % qresult.program.upper()
        else:
            program_line = '# %s %s\n' % (qresult.program.upper(), qresult.version)
        comments += program_line
        # description may or may not be present, so we'll do a try here
        try:
            comments += '# Query: %s %s\n' % (qresult.id, qresult.desc)
        except AttributeError:
            comments += '# Query: %s\n' % qresult.id
        # try appending RID line, if present
        try:
            comments += '# RID: %s\n' % qresult.rid
        except AttributeError:
            pass
        comments += '# Database: %s\n' % qresult.target
        # qresults without hits don't show the Fields comment
        if qresult:
            comments += '# Fields: %s\n' % \
                    ', '.join([inv_field_map[field] for field in self.fields])
        comments += '# %i hits found\n' % len(qresult)

        return comments


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
