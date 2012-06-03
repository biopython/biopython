# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO parser for BLAST+ tab output format."""

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


class BlastTabIndexer(SearchIndexer):

    """Indexer class for BLAST+ tab output."""

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
