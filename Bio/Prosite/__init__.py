# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# Copyright 2000 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Prosite

This module provides code to work with the prosite.dat file from
Prosite.
http://www.expasy.ch/prosite/

Tested with:
Release 15.0, July 1998
Release 16.0, July 1999
Release 17.0, Dec 2001


Classes:
Record                Holds Prosite data.
PatternHit            Holds data from a hit against a Prosite pattern.
Iterator              Iterates over entries in a Prosite file.
Dictionary            Accesses a Prosite file using a dictionary interface.
ExPASyDictionary      Accesses Prosite records from ExPASy.
RecordParser          Parses a Prosite record into a Record object.

_Scanner              Scans Prosite-formatted data.
_RecordConsumer       Consumes Prosite data to a Record object.


Functions:
scan_sequence_expasy  Scan a sequence for occurrences of Prosite patterns.
index_file            Index a Prosite file for a Dictionary.
_extract_record       Extract Prosite data from a web page.
_extract_pattern_hits Extract Prosite patterns from a web page.

"""
__all__ = [
    'Pattern',
    'Prodoc',
    ]
from types import *
import string
import re
import sgmllib
from Bio import File
from Bio import Index
from Bio.ParserSupport import *
from Bio.WWW import ExPASy
from Bio.WWW import RequestLimiter

class Record:
    """Holds information from a Prosite record.

    Members:
    name           ID of the record.  e.g. ADH_ZINC
    type           Type of entry.  e.g. PATTERN, MATRIX, or RULE
    accession      e.g. PS00387
    created        Date the entry was created.  (MMM-YYYY)
    data_update    Date the 'primary' data was last updated.
    info_update    Date data other than 'primary' data was last updated.
    pdoc           ID of the PROSITE DOCumentation.
    
    description    Free-format description.
    pattern        The PROSITE pattern.  See docs.
    matrix         List of strings that describes a matrix entry.
    rules          List of rule definitions.  (strings)

    NUMERICAL RESULTS
    nr_sp_release  SwissProt release.
    nr_sp_seqs     Number of seqs in that release of Swiss-Prot. (int)
    nr_total       Number of hits in Swiss-Prot.  tuple of (hits, seqs)
    nr_positive    True positives.  tuple of (hits, seqs)
    nr_unknown     Could be positives.  tuple of (hits, seqs)
    nr_false_pos   False positives.  tuple of (hits, seqs)
    nr_false_neg   False negatives.  (int)
    nr_partial     False negatives, because they are fragments. (int)

    COMMENTS
    cc_taxo_range  Taxonomic range.  See docs for format
    cc_max_repeat  Maximum number of repetitions in a protein
    cc_site        Interesting site.  list of tuples (pattern pos, desc.)
    cc_skip_flag   Can this entry be ignored?
    cc_matrix_type
    cc_scaling_db
    cc_author
    cc_ft_key
    cc_ft_desc

    DATA BANK REFERENCES - The following are all
                           lists of tuples (swiss-prot accession,
                                            swiss-prot name)
    dr_positive
    dr_false_neg
    dr_false_pos
    dr_potential   Potential hits, but fingerprint region not yet available.
    dr_unknown     Could possibly belong

    pdb_structs    List of PDB entries.

    """
    def __init__(self):
        self.name = ''
        self.type = ''
        self.accession = ''
        self.created = ''
        self.data_update = ''
        self.info_update = ''
        self.pdoc = ''
    
        self.description = ''
        self.pattern = ''
        self.matrix = []
        self.rules = []

        self.nr_sp_release = ''
        self.nr_sp_seqs = ''
        self.nr_total = (None, None)
        self.nr_positive = (None, None)
        self.nr_unknown = (None, None)
        self.nr_false_pos = (None, None)
        self.nr_false_neg = None
        self.nr_partial = None

        self.cc_taxo_range = ''
        self.cc_max_repeat = ''
        self.cc_site = []
        self.cc_skip_flag = ''

        self.dr_positive = []
        self.dr_false_neg = []
        self.dr_false_pos = []
        self.dr_potential = []
        self.dr_unknown = []

        self.pdb_structs = []

class PatternHit:
    """Holds information from a hit against a Prosite pattern.

    Members:
    name           ID of the record.  e.g. ADH_ZINC
    accession      e.g. PS00387
    pdoc           ID of the PROSITE DOCumentation.
    description    Free-format description.
    matches        List of tuples (start, end, sequence) where
                   start and end are indexes of the match, and sequence is
                   the sequence matched.

    """
    def __init__(self):
        self.name = None
        self.accession = None
        self.pdoc = None
        self.description = None
        self.matches = []
    def __str__(self):
        lines = []
        lines.append("%s %s %s" % (self.accession, self.pdoc, self.name))
        lines.append(self.description)
        lines.append('')
        if len(self.matches) > 1:
            lines.append("Number of matches: %s" % len(self.matches))
        for i in range(len(self.matches)):
            start, end, seq = self.matches[i]
            range_str = "%d-%d" % (start, end)
            if len(self.matches) > 1:
                lines.append("%7d %10s %s" % (i+1, range_str, seq))
            else:
                lines.append("%7s %10s %s" % (' ', range_str, seq))
        return string.join(lines, '\n')
    
class Iterator:
    """Returns one record at a time from a Prosite file.

    Methods:
    next   Return the next record from the stream, or None.

    """
    def __init__(self, handle, parser=None):
        """__init__(self, handle, parser=None)

        Create a new iterator.  handle is a file-like object.  parser
        is an optional Parser object to change the results into another form.
        If set to None, then the raw contents of the file will be returned.

        """
        if type(handle) is not FileType and type(handle) is not InstanceType:
            raise ValueError, "I expected a file handle or file-like object"
        self._uhandle = File.UndoHandle(handle)
        self._parser = parser

    def next(self):
        """next(self) -> object

        Return the next Prosite record from the file.  If no more records,
        return None.

        """
        # Skip the copyright info, if it's the first record.
        line = self._uhandle.peekline()
        if line[:2] == 'CC':
            while 1:
                line = self._uhandle.readline()
                if not line:
                    break
                if line[:2] == '//':
                    break
                if line[:2] != 'CC':
                    raise SyntaxError, \
                          "Oops, where's the copyright?"
        
        lines = []
        while 1:
            line = self._uhandle.readline()
            if not line:
                break
            lines.append(line)
            if line[:2] == '//':
                break
            
        if not lines:
            return None
            
        data = string.join(lines, '')
        if self._parser is not None:
            return self._parser.parse(File.StringHandle(data))
        return data

class Dictionary:
    """Accesses a Prosite file using a dictionary interface.

    """
    __filename_key = '__filename'
    
    def __init__(self, indexname, parser=None):
        """__init__(self, indexname, parser=None)

        Open a Prosite Dictionary.  indexname is the name of the
        index for the dictionary.  The index should have been created
        using the index_file function.  parser is an optional Parser
        object to change the results into another form.  If set to None,
        then the raw contents of the file will be returned.

        """
        self._index = Index.Index(indexname)
        self._handle = open(self._index[Dictionary.__filename_key])
        self._parser = parser

    def __len__(self):
        return len(self._index)

    def __getitem__(self, key):
        start, len = self._index[key]
        self._handle.seek(start)
        data = self._handle.read(len)
        if self._parser is not None:
            return self._parser.parse(File.StringHandle(data))
        return data

    def __getattr__(self, name):
        return getattr(self._index, name)

class ExPASyDictionary:
    """Access PROSITE at ExPASy using a read-only dictionary interface.

    """
    def __init__(self, delay=5.0, parser=None):
        """__init__(self, delay=5.0, parser=None)

        Create a new Dictionary to access PROSITE.  parser is an optional
        parser (e.g. Prosite.RecordParser) object to change the results
        into another form.  If set to None, then the raw contents of the
        file will be returned.  delay is the number of seconds to wait
        between each query.

        """
        self.parser = parser
        self.limiter = RequestLimiter(delay)

    def __len__(self):
        raise NotImplementedError, "Prosite contains lots of entries"
    def clear(self):
        raise NotImplementedError, "This is a read-only dictionary"
    def __setitem__(self, key, item):
        raise NotImplementedError, "This is a read-only dictionary"
    def update(self):
        raise NotImplementedError, "This is a read-only dictionary"
    def copy(self):
        raise NotImplementedError, "You don't need to do this..."
    def keys(self):
        raise NotImplementedError, "You don't really want to do this..."
    def items(self):
        raise NotImplementedError, "You don't really want to do this..."
    def values(self):
        raise NotImplementedError, "You don't really want to do this..."
    
    def has_key(self, id):
        """has_key(self, id) -> bool"""
        try:
            self[id]
        except KeyError:
            return 0
        return 1

    def get(self, id, failobj=None):
        try:
            return self[id]
        except KeyError:
            return failobj
        raise "How did I get here?"

    def __getitem__(self, id):
        """__getitem__(self, id) -> object

        Return a Prosite entry.  id is either the id or accession
        for the entry.  Raises a KeyError if there's an error.
        
        """
        # First, check to see if enough time has passed since my
        # last query.
        self.limiter.wait()

        try:
            handle = ExPASy.get_prosite_entry(id)
        except IOError:
            raise KeyError, id
        try:
            handle = File.StringHandle(_extract_record(handle))
        except ValueError:
            raise KeyError, id
        
        if self.parser is not None:
            return self.parser.parse(handle)
        return handle.read()

class RecordParser(AbstractParser):
    """Parses Prosite data into a Record object.

    """
    def __init__(self):
        self._scanner = _Scanner()
        self._consumer = _RecordConsumer()

    def parse(self, handle):
        self._scanner.feed(handle, self._consumer)
        return self._consumer.data

class _Scanner:
    """Scans Prosite-formatted data.

    Tested with:
    Release 15.0, July 1998
    
    """
    def feed(self, handle, consumer):
        """feed(self, handle, consumer)

        Feed in Prosite data for scanning.  handle is a file-like
        object that contains prosite data.  consumer is a
        Consumer object that will receive events as the report is scanned.

        """
        if isinstance(handle, File.UndoHandle):
            uhandle = handle
        else:
            uhandle = File.UndoHandle(handle)

        while 1:
            line = uhandle.peekline()
            if not line:
                break
            elif is_blank_line(line):
                # Skip blank lines between records
                uhandle.readline()
                continue
            elif line[:2] == 'ID':
                self._scan_record(uhandle, consumer)
            elif line[:2] == 'CC':
                self._scan_copyrights(uhandle, consumer)
            else:
                raise SyntaxError, "There doesn't appear to be a record"

    def _scan_copyrights(self, uhandle, consumer):
        consumer.start_copyrights()
        self._scan_line('CC', uhandle, consumer.copyright, any_number=1)
        self._scan_terminator(uhandle, consumer)
        consumer.end_copyrights()

    def _scan_record(self, uhandle, consumer):
        consumer.start_record()
        for fn in self._scan_fns:
            fn(self, uhandle, consumer)

            # In Release 15.0, C_TYPE_LECTIN_1 has the DO line before
            # the 3D lines, instead of the other way around.
            # Thus, I'll give the 3D lines another chance after the DO lines
            # are finished.
            if fn is self._scan_do.im_func:
                self._scan_3d(uhandle, consumer)
        consumer.end_record()

    def _scan_line(self, line_type, uhandle, event_fn,
                   exactly_one=None, one_or_more=None, any_number=None,
                   up_to_one=None):
        # Callers must set exactly one of exactly_one, one_or_more, or
        # any_number to a true value.  I do not explicitly check to
        # make sure this function is called correctly.
        
        # This does not guarantee any parameter safety, but I
        # like the readability.  The other strategy I tried was have
        # parameters min_lines, max_lines.
        
        if exactly_one or one_or_more:
            read_and_call(uhandle, event_fn, start=line_type)
        if one_or_more or any_number:
            while 1:
                if not attempt_read_and_call(uhandle, event_fn,
                                             start=line_type):
                    break
        if up_to_one:
            attempt_read_and_call(uhandle, event_fn, start=line_type)

    def _scan_id(self, uhandle, consumer):
        self._scan_line('ID', uhandle, consumer.identification, exactly_one=1)

    def _scan_ac(self, uhandle, consumer):
        self._scan_line('AC', uhandle, consumer.accession, exactly_one=1)
    
    def _scan_dt(self, uhandle, consumer):
        self._scan_line('DT', uhandle, consumer.date, exactly_one=1)

    def _scan_de(self, uhandle, consumer):
        self._scan_line('DE', uhandle, consumer.description, exactly_one=1)
    
    def _scan_pa(self, uhandle, consumer):
        self._scan_line('PA', uhandle, consumer.pattern, any_number=1)
    
    def _scan_ma(self, uhandle, consumer):
        self._scan_line('MA', uhandle, consumer.matrix, any_number=1)
##        # ZN2_CY6_FUNGAL_2, DNAJ_2 in Release 15
##        # contain a CC line buried within an 'MA' line.  Need to check
##        # for that.
##        while 1:
##            if not attempt_read_and_call(uhandle, consumer.matrix, start='MA'):
##                line1 = uhandle.readline()
##                line2 = uhandle.readline()
##                uhandle.saveline(line2)
##                uhandle.saveline(line1)
##                if line1[:2] == 'CC' and line2[:2] == 'MA':
##                    read_and_call(uhandle, consumer.comment, start='CC')
##                else:
##                    break
    
    def _scan_ru(self, uhandle, consumer):
        self._scan_line('RU', uhandle, consumer.rule, any_number=1)
    
    def _scan_nr(self, uhandle, consumer):
        self._scan_line('NR', uhandle, consumer.numerical_results,
                        any_number=1)

    def _scan_cc(self, uhandle, consumer):
        self._scan_line('CC', uhandle, consumer.comment, any_number=1)
    
    def _scan_dr(self, uhandle, consumer):
        self._scan_line('DR', uhandle, consumer.database_reference,
                        any_number=1)
    
    def _scan_3d(self, uhandle, consumer):
        self._scan_line('3D', uhandle, consumer.pdb_reference,
                        any_number=1)
    
    def _scan_do(self, uhandle, consumer):
        self._scan_line('DO', uhandle, consumer.documentation, exactly_one=1)

    def _scan_terminator(self, uhandle, consumer):
        self._scan_line('//', uhandle, consumer.terminator, exactly_one=1)

    _scan_fns = [
        _scan_id,
        _scan_ac,
        _scan_dt,
        _scan_de,
        _scan_pa,
        _scan_ma,
        _scan_ru,
        _scan_nr,
        _scan_cc,

        # This is a really dirty hack, and should be fixed properly at
        # some point.  ZN2_CY6_FUNGAL_2, DNAJ_2 in Rel 15 and PS50309
        # in Rel 17 have lines out of order.  Thus, I have to rescan
        # these, which decreases performance.
        _scan_ma,   
        _scan_nr,
        _scan_cc,

        _scan_dr,
        _scan_3d,
        _scan_do,
        _scan_terminator
        ]

class _RecordConsumer(AbstractConsumer):
    """Consumer that converts a Prosite record to a Record object.

    Members:
    data    Record with Prosite data.

    """
    def __init__(self):
        self.data = None
        
    def start_record(self):
        self.data = Record()
        
    def end_record(self):
        self._clean_record(self.data)

    def identification(self, line):
        cols = string.split(line)
        if len(cols) != 3:
            raise SyntaxError, "I don't understand identification line\n%s" % \
                  line
        self.data.name = self._chomp(cols[1])    # don't want ';'
        self.data.type = self._chomp(cols[2])    # don't want '.'
    
    def accession(self, line):
        cols = string.split(line)
        if len(cols) != 2:
            raise SyntaxError, "I don't understand accession line\n%s" % line
        self.data.accession = self._chomp(cols[1])
    
    def date(self, line):
        uprline = string.upper(line)
        cols = string.split(uprline)

        # Release 15.0 contains both 'INFO UPDATE' and 'INF UPDATE'
        if cols[2] != '(CREATED);' or \
           cols[4] != '(DATA' or cols[5] != 'UPDATE);' or \
           cols[7][:4] != '(INF' or cols[8] != 'UPDATE).':
            raise SyntaxError, "I don't understand date line\n%s" % line

        self.data.created = cols[1]
        self.data.data_update = cols[3]
        self.data.info_update = cols[6]
    
    def description(self, line):
        self.data.description = self._clean(line)
    
    def pattern(self, line):
        self.data.pattern = self.data.pattern + self._clean(line)
    
    def matrix(self, line):
        self.data.matrix.append(self._clean(line))
    
    def rule(self, line):
        self.data.rules.append(self._clean(line))
    
    def numerical_results(self, line):
        cols = string.split(self._clean(line), ';')
        for col in cols:
            if not col:
                continue
            qual, data = map(string.lstrip, string.split(col, '='))
            if qual == '/RELEASE':
                release, seqs = string.split(data, ',')
                self.data.nr_sp_release = release
                self.data.nr_sp_seqs = int(seqs)
            elif qual == '/FALSE_NEG':
                self.data.nr_false_neg = int(data)
            elif qual == '/PARTIAL':
                self.data.nr_partial = int(data)
            elif qual in ['/TOTAL', '/POSITIVE', '/UNKNOWN', '/FALSE_POS']:
                m = re.match(r'(\d+)\((\d+)\)', data)
                if not m:
                    raise error, "Broken data %s in comment line\n%s" % \
                          (repr(data), line)
                hits = tuple(map(int, m.groups()))
                if(qual == "/TOTAL"):
                    self.data.nr_total = hits
                elif(qual == "/POSITIVE"):
                    self.data.nr_positive = hits
                elif(qual == "/UNKNOWN"):
                    self.data.nr_unknown = hits
                elif(qual == "/FALSE_POS"):
                    self.data.nr_false_pos = hits
            else:
                raise SyntaxError, "Unknown qual %s in comment line\n%s" % \
                      (repr(qual), line)
    
    def comment(self, line):
        cols = string.split(self._clean(line), ';')
        for col in cols:
            # DNAJ_2 in Release 15 has a non-standard comment line:
            # CC   Automatic scaling using reversed database
            # Throw it away.  (Should I keep it?)
            if not col or col[:17] == 'Automatic scaling':
                continue
            qual, data = map(string.lstrip, string.split(col, '='))
            if qual == '/TAXO-RANGE':
                self.data.cc_taxo_range = data
            elif qual == '/MAX-REPEAT':
                self.data.cc_max_repeat = data
            elif qual == '/SITE':
                pos, desc = string.split(data, ',')
                self.data.cc_site = (int(pos), desc)
            elif qual == '/SKIP-FLAG':
                self.data.cc_skip_flag = data
            elif qual == '/MATRIX_TYPE':
                self.data.cc_matrix_type = data
            elif qual == '/SCALING_DB':
                self.data.cc_scaling_db = data
            elif qual == '/AUTHOR':
                self.data.cc_author = data
            elif qual == '/FT_KEY':
                self.data.cc_ft_key = data
            elif qual == '/FT_DESC':
                self.data.cc_ft_desc = data
            else:
                raise SyntaxError, "Unknown qual %s in comment line\n%s" % \
                      (repr(qual), line)
            
    def database_reference(self, line):
        refs = string.split(self._clean(line), ';')
        for ref in refs:
            if not ref:
                continue
            acc, name, type = map(string.strip, string.split(ref, ','))
            if type == 'T':
                self.data.dr_positive.append((acc, name))
            elif type == 'F':
                self.data.dr_false_pos.append((acc, name))
            elif type == 'N':
                self.data.dr_false_neg.append((acc, name))
            elif type == 'P':
                self.data.dr_potential.append((acc, name))
            elif type == '?':
                self.data.dr_unknown.append((acc, name))
            else:
                raise SyntaxError, "I don't understand type flag %s" % type
    
    def pdb_reference(self, line):
        cols = string.split(line)
        for id in cols[1:]:  # get all but the '3D' col
            self.data.pdb_structs.append(self._chomp(id))
    
    def documentation(self, line):
        self.data.pdoc = self._chomp(self._clean(line))

    def terminator(self, line):
        pass

    def _chomp(self, word, to_chomp='.,;'):
        # Remove the punctuation at the end of a word.
        if word[-1] in to_chomp:
            return word[:-1]
        return word

    def _clean(self, line, rstrip=1):
        # Clean up a line.
        if rstrip:
            return string.rstrip(line[5:])
        return line[5:]

def scan_sequence_expasy(seq=None, id=None, exclude_frequent=None):
    """scan_sequence_expasy(seq=None, id=None, exclude_frequent=None) ->
    list of PatternHit's

    Search a sequence for occurrences of Prosite patterns.  You can
    specify either a sequence in seq or a SwissProt/trEMBL ID or accession
    in id.  Only one of those should be given.  If exclude_frequent
    is true, then the patterns with the high probability of occurring
    will be excluded.

    """
    if (seq and id) or not (seq or id):
        raise ValueError, "Please specify either a sequence or an id"
    handle = ExPASy.scanprosite1(seq, id, exclude_frequent)
    return _extract_pattern_hits(handle)

def _extract_pattern_hits(handle):
    """_extract_pattern_hits(handle) -> list of PatternHit's

    Extract hits from a web page.  Raises a ValueError if there
    was an error in the query.

    """
    class parser(sgmllib.SGMLParser):
        def __init__(self):
            sgmllib.SGMLParser.__init__(self)
            self.hits = []
            self.broken_message = 'Some error occurred'
            self._in_pre = 0
            self._current_hit = None
            self._last_found = None   # Save state of parsing
        def handle_data(self, data):
            if string.find(data, 'try again') >= 0:
                self.broken_message = data
                return
            elif data == 'illegal':
                self.broken_message = 'Sequence contains illegal characters'
                return
            if not self._in_pre:
                return
            elif not string.strip(data):
                return
            if self._last_found is None and data[:4] == 'PDOC':
                self._current_hit.pdoc = data
                self._last_found = 'pdoc'
            elif self._last_found == 'pdoc':
                if data[:2] != 'PS':
                    raise SyntaxError, "Expected accession but got:\n%s" % data
                self._current_hit.accession = data
                self._last_found = 'accession'
            elif self._last_found == 'accession':
                self._current_hit.name = data
                self._last_found = 'name'
            elif self._last_found == 'name':
                self._current_hit.description = data
                self._last_found = 'description'
            elif self._last_found == 'description':
                m = re.findall(r'(\d+)-(\d+) (\w+)', data)
                for start, end, seq in m:
                    self._current_hit.matches.append(
                        (int(start), int(end), seq))
            
        def do_hr(self, attrs):
            # <HR> inside a <PRE> section means a new hit.
            if self._in_pre:
                self._current_hit = PatternHit()
                self.hits.append(self._current_hit)
                self._last_found = None
        def start_pre(self, attrs):
            self._in_pre = 1
            self.broken_message = None   # Probably not broken
        def end_pre(self):
            self._in_pre = 0
    p = parser()
    p.feed(handle.read())
    if p.broken_message:
        raise ValueError, p.broken_message
    return p.hits


        
    
def index_file(filename, indexname, rec2key=None):
    """index_file(filename, indexname, rec2key=None)

    Index a Prosite file.  filename is the name of the file.
    indexname is the name of the dictionary.  rec2key is an
    optional callback that takes a Record and generates a unique key
    (e.g. the accession number) for the record.  If not specified,
    the id name will be used.

    """
    if not os.path.exists(filename):
        raise ValueError, "%s does not exist" % filename

    index = Index.Index(indexname, truncate=1)
    index[Dictionary._Dictionary__filename_key] = filename
    
    iter = Iterator(open(filename), parser=RecordParser())
    while 1:
        start = iter._uhandle.tell()
        rec = iter.next()
        length = iter._uhandle.tell() - start
        
        if rec is None:
            break
        if rec2key is not None:
            key = rec2key(rec)
        else:
            key = rec.name
            
        if not key:
            raise KeyError, "empty key was produced"
        elif index.has_key(key):
            raise KeyError, "duplicate key %s found" % key

        index[key] = start, length

def _extract_record(handle):
    """_extract_record(handle) -> str

    Extract PROSITE data from a web page.  Raises a ValueError if no
    data was found in the web page.

    """
    # All the data appears between tags:
    # <pre width = 80>ID   NIR_SIR; PATTERN.
    # </PRE>
    class parser(sgmllib.SGMLParser):
        def __init__(self):
            sgmllib.SGMLParser.__init__(self)
            self._in_pre = 0
            self.data = []
        def handle_data(self, data):
            if self._in_pre:
                self.data.append(data)
        def do_br(self, attrs):
            if self._in_pre:
                self.data.append('\n')
        def start_pre(self, attrs):
            self._in_pre = 1
        def end_pre(self):
            self._in_pre = 0
    p = parser()
    p.feed(handle.read())
    if not p.data:
        raise ValueError, "No data found in web page."
    return string.join(p.data, '')

