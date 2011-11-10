# Copyright 2000 by Jeffrey Chang, Brad Chapman.  All rights reserved.
# Copyright 2006-2011 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Code to work with GenBank formatted files.

Rather than using Bio.GenBank, you are now encouraged to use Bio.SeqIO with
the "genbank" or "embl" format names to parse GenBank or EMBL files into
SeqRecord and SeqFeature objects (see the Biopython tutorial for details).

Using Bio.GenBank directly to parse GenBank files is only useful if you want 
to obtain GenBank-specific Record objects, which is a much closer
representation to the raw file contents that the SeqRecord alternative from
the FeatureParser (used in Bio.SeqIO).

To use the Bio.GenBank parser, there are two helper functions:

read                  Parse a handle containing a single GenBank record
                      as Bio.GenBank specific Record objects.
parse                 Iterate over a handle containing multiple GenBank
                      records as Bio.GenBank specific Record objects.

The following internal classes are not intended for direct use and may
be deprecated in a future release.

Classes:
Iterator              Iterate through a file of GenBank entries
ErrorFeatureParser    Catch errors caused during parsing.
FeatureParser         Parse GenBank data in SeqRecord and SeqFeature objects.
RecordParser          Parse GenBank data into a Record object.

Exceptions:
ParserFailureError    Exception indicating a failure in the parser (ie.
                      scanner or consumer)
LocationParserError   Exception indiciating a problem with the spark based
                      location parser.

"""
import re

# other Biopython stuff
from Bio import SeqFeature

# other Bio.GenBank stuff
from utils import FeatureValueCleaner
from Scanner import GenBankScanner

#Constants used to parse GenBank header lines
GENBANK_INDENT = 12
GENBANK_SPACER = " " * GENBANK_INDENT

#Constants for parsing GenBank feature lines
FEATURE_KEY_INDENT = 5
FEATURE_QUALIFIER_INDENT = 21
FEATURE_KEY_SPACER = " " * FEATURE_KEY_INDENT
FEATURE_QUALIFIER_SPACER = " " * FEATURE_QUALIFIER_INDENT

#Regular expresions for location parsing
_solo_location = r"[<>]?\d+"
_pair_location = r"[<>]?\d+\.\.[<>]?\d+"
_between_location = r"\d+\^\d+"

_within_position = r"\(\d+\.\d+\)"
_re_within_position = re.compile(_within_position)
_within_location = r"([<>]?\d+|%s)\.\.([<>]?\d+|%s)" \
                   % (_within_position,_within_position)
assert _re_within_position.match("(3.9)")
assert re.compile(_within_location).match("(3.9)..10")
assert re.compile(_within_location).match("26..(30.33)")
assert re.compile(_within_location).match("(13.19)..(20.28)")

_oneof_position = r"one\-of\(\d+(,\d+)+\)"
_re_oneof_position = re.compile(_oneof_position)
_oneof_location = r"([<>]?\d+|%s)\.\.([<>]?\d+|%s)" \
                   % (_oneof_position,_oneof_position)
assert _re_oneof_position.match("one-of(6,9)")
assert re.compile(_oneof_location).match("one-of(6,9)..101")
assert re.compile(_oneof_location).match("one-of(6,9)..one-of(101,104)")
assert re.compile(_oneof_location).match("6..one-of(101,104)")

assert not _re_oneof_position.match("one-of(3)")
assert _re_oneof_position.match("one-of(3,6)")
assert _re_oneof_position.match("one-of(3,6,9)")


_simple_location = r"\d+\.\.\d+"
_re_simple_location = re.compile(_simple_location)
_re_simple_compound = re.compile(r"^(join|order|bond)\(%s(,%s)*\)$" \
                                 % (_simple_location, _simple_location))
_complex_location = r"([a-zA-z][a-zA-Z0-9]*(\.[a-zA-Z0-9]+)?\:)?(%s|%s|%s|%s|%s)" \
                    % (_pair_location, _solo_location, _between_location,
                       _within_location, _oneof_location)
_re_complex_location = re.compile(r"^%s$" % _complex_location)
_possibly_complemented_complex_location = r"(%s|complement\(%s\))" \
                                          % (_complex_location, _complex_location)
_re_complex_compound = re.compile(r"^(join|order|bond)\(%s(,%s)*\)$" \
                                 % (_possibly_complemented_complex_location,
                                    _possibly_complemented_complex_location))

assert _re_simple_location.match("104..160")
assert not _re_simple_location.match("<104..>160")
assert not _re_simple_location.match("104")
assert not _re_simple_location.match("<1")
assert not _re_simple_location.match(">99999")
assert not _re_simple_location.match("join(104..160,320..390,504..579)")
assert not _re_simple_compound.match("bond(12,63)")
assert _re_simple_compound.match("join(104..160,320..390,504..579)")
assert _re_simple_compound.match("order(1..69,1308..1465)")
assert not _re_simple_compound.match("order(1..69,1308..1465,1524)")
assert not _re_simple_compound.match("join(<1..442,992..1228,1524..>1983)")
assert not _re_simple_compound.match("join(<1..181,254..336,422..497,574..>590)")
assert not _re_simple_compound.match("join(1475..1577,2841..2986,3074..3193,3314..3481,4126..>4215)")
assert not _re_simple_compound.match("test(1..69,1308..1465)")
assert not _re_simple_compound.match("complement(1..69)")
assert not _re_simple_compound.match("(1..69)")
assert _re_complex_location.match("(3.9)..10")
assert _re_complex_location.match("26..(30.33)")
assert _re_complex_location.match("(13.19)..(20.28)")
assert _re_complex_location.match("41^42") #between
assert _re_complex_location.match("AL121804:41^42")
assert _re_complex_location.match("AL121804:41..610")
assert _re_complex_location.match("AL121804.2:41..610")
assert _re_complex_location.match("one-of(3,6)..101")
assert _re_complex_compound.match("join(153490..154269,AL121804.2:41..610,AL121804.2:672..1487)")
assert not _re_simple_compound.match("join(153490..154269,AL121804.2:41..610,AL121804.2:672..1487)")
assert _re_complex_compound.match("join(complement(69611..69724),139856..140650)")

def _pos(pos_str, offset=0):
    """Build a Position object (PRIVATE).
    
    For an end position, leave offset as zero (default):

    >>> _pos("5")
    ExactPosition(5)

    For a start position, set offset to minus one (for Python counting):

    >>> _pos("5", -1)
    ExactPosition(4)

    This also covers fuzzy positions:

    >>> p = _pos("<5")
    >>> p
    BeforePosition(5)
    >>> print p
    <5
    >>> int(p)
    5

    >>> _pos(">5")
    AfterPosition(5)

    By default assumes an end position, so note the integer behaviour:

    >>> p = _pos("one-of(5,8,11)")
    >>> p
    OneOfPosition(11, choices=[ExactPosition(5), ExactPosition(8), ExactPosition(11)])
    >>> print p
    one-of(5,8,11)
    >>> int(p)
    11

    >>> _pos("(8.10)")
    WithinPosition(10, left=8, right=10)

    Fuzzy start positions:

    >>> p = _pos("<5", -1)
    >>> p
    BeforePosition(4)
    >>> print p
    <4
    >>> int(p)
    4

    Notice how the integer behaviour changes too!

    >>> p = _pos("one-of(5,8,11)", -1)
    >>> p
    OneOfPosition(4, choices=[ExactPosition(4), ExactPosition(7), ExactPosition(10)])
    >>> print(p)
    one-of(4,7,10)
    >>> int(p)
    4

    """
    if pos_str.startswith("<"):
        return SeqFeature.BeforePosition(int(pos_str[1:])+offset)
    elif pos_str.startswith(">"):
        return SeqFeature.AfterPosition(int(pos_str[1:])+offset)
    elif _re_within_position.match(pos_str):
        s,e = pos_str[1:-1].split(".")
        s = int(s) + offset
        e = int(e) + offset
        if offset == -1:
            default = s
        else:
            default = e
        return SeqFeature.WithinPosition(default, left=s, right=e)
    elif _re_oneof_position.match(pos_str):
        assert pos_str.startswith("one-of(")
        assert pos_str[-1]==")"
        parts = [SeqFeature.ExactPosition(int(pos)+offset) \
                 for pos in pos_str[7:-1].split(",")]
        if offset == -1:
            default = min(int(pos) for pos in parts)
        else:
            default = max(int(pos) for pos in parts)
        return SeqFeature.OneOfPosition(default, choices=parts)
    else:
        return SeqFeature.ExactPosition(int(pos_str)+offset)

def _loc(loc_str, expected_seq_length, strand):
    """FeatureLocation from non-compound non-complement location (PRIVATE).
    
    Simple examples,

    >>> _loc("123..456", 1000, +1)
    FeatureLocation(ExactPosition(122), ExactPosition(456), strand=1)
    >>> _loc("<123..>456", 1000, strand = -1)
    FeatureLocation(BeforePosition(122), AfterPosition(456), strand=-1)

    A more complex location using within positions,

    >>> _loc("(9.10)..(20.25)", 1000, 1)
    FeatureLocation(WithinPosition(8, left=8, right=9), WithinPosition(25, left=20, right=25), strand=1)

    Notice how that will act as though it has overall start 8 and end 25.

    Zero length between feature,

    >>> _loc("123^124", 1000, 0)
    FeatureLocation(ExactPosition(123), ExactPosition(123), strand=0)
    
    The expected sequence length is needed for a special case, a between
    position at the start/end of a circular genome:

    >>> _loc("1000^1", 1000, 1)
    FeatureLocation(ExactPosition(1000), ExactPosition(1000), strand=1)
    
    Apart from this special case, between positions P^Q must have P+1==Q,

    >>> _loc("123^456", 1000, 1)
    Traceback (most recent call last):
       ...
    ValueError: Invalid between location '123^456'
    """
    try:
        s, e = loc_str.split("..")
    except ValueError:
        assert ".." not in loc_str
        if "^" in loc_str:
            #A between location like "67^68" (one based counting) is a
            #special case (note it has zero length). In python slice
            #notation this is 67:67, a zero length slice.  See Bug 2622
            #Further more, on a circular genome of length N you can have
            #a location N^1 meaning the junction at the origin. See Bug 3098.
            #NOTE - We can imagine between locations like "2^4", but this
            #is just "3".  Similarly, "2^5" is just "3..4"
            s, e = loc_str.split("^")
            if int(s)+1==int(e):
                pos = _pos(s)
            elif int(s)==expected_seq_length and e=="1":
                pos = _pos(s)
            else:
                raise ValueError("Invalid between location %s" % repr(loc_str))
            return SeqFeature.FeatureLocation(pos, pos, strand)
        else:
            #e.g. "123"
            s = loc_str
            e = loc_str
    return SeqFeature.FeatureLocation(_pos(s,-1), _pos(e), strand)

def _split_compound_loc(compound_loc):
    """Split a tricky compound location string (PRIVATE).
    
    >>> list(_split_compound_loc("123..145"))
    ['123..145']
    >>> list(_split_compound_loc("123..145,200..209"))
    ['123..145', '200..209']
    >>> list(_split_compound_loc("one-of(200,203)..300"))
    ['one-of(200,203)..300']
    >>> list(_split_compound_loc("complement(123..145),200..209"))
    ['complement(123..145)', '200..209']
    >>> list(_split_compound_loc("123..145,one-of(200,203)..209"))
    ['123..145', 'one-of(200,203)..209']
    >>> list(_split_compound_loc("123..145,one-of(200,203)..one-of(209,211),300"))
    ['123..145', 'one-of(200,203)..one-of(209,211)', '300']
    >>> list(_split_compound_loc("123..145,complement(one-of(200,203)..one-of(209,211)),300"))
    ['123..145', 'complement(one-of(200,203)..one-of(209,211))', '300']
    >>> list(_split_compound_loc("123..145,200..one-of(209,211),300"))
    ['123..145', '200..one-of(209,211)', '300']
    >>> list(_split_compound_loc("123..145,200..one-of(209,211)"))
    ['123..145', '200..one-of(209,211)']
    """
    if "one-of(" in compound_loc:
        #Hard case
        while "," in compound_loc:
            assert compound_loc[0] != ","
            assert compound_loc[0:2] != ".."
            i = compound_loc.find(",")
            part = compound_loc[:i]
            compound_loc = compound_loc[i:] #includes the comma
            while part.count("(") > part.count(")"):
                assert "one-of(" in part, (part, compound_loc)
                i = compound_loc.find(")")
                part += compound_loc[:i+1]
                compound_loc = compound_loc[i+1:]
            if compound_loc.startswith(".."):
                i = compound_loc.find(",")
                if i==-1:
                    part += compound_loc
                    compound_loc = ""
                else:
                    part += compound_loc[:i]
                    compound_loc = compound_loc[i:] #includes the comma
            while part.count("(") > part.count(")"):
                assert part.count("one-of(") == 2
                i = compound_loc.find(")")
                part += compound_loc[:i+1]
                compound_loc = compound_loc[i+1:]
            if compound_loc.startswith(","):
                compound_loc = compound_loc[1:]
            assert part
            yield part
        if compound_loc:
            yield compound_loc
    else:
        #Easy case
        for part in compound_loc.split(","):
            yield part

class Iterator(object):
    """Iterator interface to move over a file of GenBank entries one at a time (OBSOLETE).

    This class is likely to be deprecated in a future release of Biopython.
    Please use Bio.SeqIO.parse(..., format="gb") or Bio.GenBank.parse(...)
    for SeqRecord and GenBank specific Record objects respectively instead.
    """
    def __init__(self, handle, parser = None):
        """Initialize the iterator.

        Arguments:
        o handle - A handle with GenBank entries to iterate through.
        o parser - An optional parser to pass the entries through before
        returning them. If None, then the raw entry will be returned.
        """
        self.handle = handle
        self._parser = parser

    def next(self):
        """Return the next GenBank record from the handle.

        Will return None if we ran out of records.
        """
        if self._parser is None:
            lines = []
            while True:
                line = self.handle.readline()
                if not line : return None #Premature end of file?
                lines.append(line)
                if line.rstrip() == "//" : break
            return "".join(lines)
        try:
            return self._parser.parse(self.handle)
        except StopIteration:
            return None

    def __iter__(self):
        return iter(self.next, None)

class ParserFailureError(Exception):
    """Failure caused by some kind of problem in the parser.
    """
    pass

class LocationParserError(Exception):
    """Could not Properly parse out a location from a GenBank file.
    """
    pass
                                                          
class FeatureParser(object):
    """Parse GenBank files into Seq + Feature objects (OBSOLETE).

    Direct use of this class is discouraged, and may be deprecated in
    a future release of Biopython.

    Please use Bio.SeqIO.parse(...) or Bio.SeqIO.read(...) instead.
    """
    def __init__(self, debug_level = 0, use_fuzziness = 1, 
                 feature_cleaner = FeatureValueCleaner()):
        """Initialize a GenBank parser and Feature consumer.

        Arguments:
        o debug_level - An optional argument that species the amount of
        debugging information the parser should spit out. By default we have
        no debugging info (the fastest way to do things), but if you want
        you can set this as high as two and see exactly where a parse fails.
        o use_fuzziness - Specify whether or not to use fuzzy representations.
        The default is 1 (use fuzziness).
        o feature_cleaner - A class which will be used to clean out the
        values of features. This class must implement the function 
        clean_value. GenBank.utils has a "standard" cleaner class, which
        is used by default.
        """
        self._scanner = GenBankScanner(debug_level)
        self.use_fuzziness = use_fuzziness
        self._cleaner = feature_cleaner

    def parse(self, handle):
        """Parse the specified handle.
        """
        self._consumer = _FeatureConsumer(self.use_fuzziness, 
                                          self._cleaner)
        self._scanner.feed(handle, self._consumer)
        return self._consumer.data

class RecordParser(object):
    """Parse GenBank files into Record objects (OBSOLETE).

    Direct use of this class is discouraged, and may be deprecated in
    a future release of Biopython.

    Please use the Bio.GenBank.parse(...) or Bio.GenBank.read(...) functions
    instead.
    """
    def __init__(self, debug_level = 0):
        """Initialize the parser.

        Arguments:
        o debug_level - An optional argument that species the amount of
        debugging information the parser should spit out. By default we have
        no debugging info (the fastest way to do things), but if you want
        you can set this as high as two and see exactly where a parse fails.
        """
        self._scanner = GenBankScanner(debug_level)

    def parse(self, handle):
        """Parse the specified handle into a GenBank record.
        """
        self._consumer = _RecordConsumer()

        self._scanner.feed(handle, self._consumer)
        return self._consumer.data

class _BaseGenBankConsumer(object):
    """Abstract GenBank consumer providing useful general functions (PRIVATE).

    This just helps to eliminate some duplication in things that most
    GenBank consumers want to do.
    """
    # Special keys in GenBank records that we should remove spaces from
    # For instance, \translation keys have values which are proteins and
    # should have spaces and newlines removed from them. This class
    # attribute gives us more control over specific formatting problems.
    remove_space_keys = ["translation"]

    def __init__(self):
        pass

    def _unhandled(self, data):
        pass

    def __getattr__(self, attr):
        return self._unhandled

    def _split_keywords(self, keyword_string):
        """Split a string of keywords into a nice clean list.
        """
        # process the keywords into a python list
        if keyword_string == "" or keyword_string == ".":
            keywords = ""
        elif keyword_string[-1] == '.':
            keywords = keyword_string[:-1]
        else:
            keywords = keyword_string
        keyword_list = keywords.split(';')
        clean_keyword_list = [x.strip() for x in keyword_list]
        return clean_keyword_list

    def _split_accessions(self, accession_string):
        """Split a string of accession numbers into a list.
        """
        # first replace all line feeds with spaces
        # Also, EMBL style accessions are split with ';'
        accession = accession_string.replace("\n", " ").replace(";"," ")

        return [x.strip() for x in accession.split() if x.strip()]

    def _split_taxonomy(self, taxonomy_string):
        """Split a string with taxonomy info into a list.
        """
        if not taxonomy_string or taxonomy_string==".":
            #Missing data, no taxonomy
            return []
        
        if taxonomy_string[-1] == '.':
            tax_info = taxonomy_string[:-1]
        else:
            tax_info = taxonomy_string
        tax_list = tax_info.split(';')
        new_tax_list = []
        for tax_item in tax_list:
            new_items = tax_item.split("\n")
            new_tax_list.extend(new_items)
        while '' in new_tax_list:
            new_tax_list.remove('')
        clean_tax_list = [x.strip() for x in new_tax_list]

        return clean_tax_list

    def _clean_location(self, location_string):
        """Clean whitespace out of a location string.

        The location parser isn't a fan of whitespace, so we clean it out
        before feeding it into the parser.
        """
        #Originally this imported string.whitespace and did a replace
        #via a loop.  It's simpler to just split on whitespace and rejoin
        #the string - and this avoids importing string too.  See Bug 2684.
        return ''.join(location_string.split())

    def _remove_newlines(self, text):
        """Remove any newlines in the passed text, returning the new string.
        """
        # get rid of newlines in the qualifier value
        newlines = ["\n", "\r"]
        for ws in newlines:
            text = text.replace(ws, "")

        return text

    def _normalize_spaces(self, text):
        """Replace multiple spaces in the passed text with single spaces.
        """
        # get rid of excessive spaces
        text_parts = text.split(" ")
        text_parts = filter(None, text_parts)
        return ' '.join(text_parts)

    def _remove_spaces(self, text):
        """Remove all spaces from the passed text.
        """
        return text.replace(" ", "")

    def _convert_to_python_numbers(self, start, end):
        """Convert a start and end range to python notation.

        In GenBank, starts and ends are defined in "biological" coordinates,
        where 1 is the first base and [i, j] means to include both i and j.

        In python, 0 is the first base and [i, j] means to include i, but
        not j. 

        So, to convert "biological" to python coordinates, we need to 
        subtract 1 from the start, and leave the end and things should
        be converted happily.
        """
        new_start = start - 1
        new_end = end

        return new_start, new_end

class _FeatureConsumer(_BaseGenBankConsumer):
    """Create a SeqRecord object with Features to return (PRIVATE).

    Attributes:
    o use_fuzziness - specify whether or not to parse with fuzziness in
    feature locations.
    o feature_cleaner - a class that will be used to provide specialized
    cleaning-up of feature values.
    """
    def __init__(self, use_fuzziness, feature_cleaner = None):
        from Bio.SeqRecord import SeqRecord
        _BaseGenBankConsumer.__init__(self)
        self.data = SeqRecord(None, id = None)
        self.data.id = None
        self.data.description = ""

        self._use_fuzziness = use_fuzziness
        self._feature_cleaner = feature_cleaner

        self._seq_type = ''
        self._seq_data = []
        self._cur_reference = None
        self._cur_feature = None
        self._expected_size = None

    def locus(self, locus_name):
        """Set the locus name is set as the name of the Sequence.
        """
        self.data.name = locus_name

    def size(self, content):
        """Record the sequence length."""
        self._expected_size = int(content)

    def residue_type(self, type):
        """Record the sequence type so we can choose an appropriate alphabet.
        """
        self._seq_type = type

    def data_file_division(self, division):
        self.data.annotations['data_file_division'] = division

    def date(self, submit_date):
        self.data.annotations['date'] = submit_date 

    def definition(self, definition):
        """Set the definition as the description of the sequence.
        """
        if self.data.description:
            #Append to any existing description
            #e.g. EMBL files with two DE lines.
            self.data.description += " " + definition
        else:
            self.data.description = definition

    def accession(self, acc_num):
        """Set the accession number as the id of the sequence.

        If we have multiple accession numbers, the first one passed is
        used.
        """
        new_acc_nums = self._split_accessions(acc_num)

        #Also record them ALL in the annotations
        try:
            #On the off chance there was more than one accession line:
            for acc in new_acc_nums:
                #Prevent repeat entries
                if acc not in self.data.annotations['accessions']:
                    self.data.annotations['accessions'].append(acc)
        except KeyError:
            self.data.annotations['accessions'] = new_acc_nums

        # if we haven't set the id information yet, add the first acc num
        if self.data.id is None:
            if len(new_acc_nums) > 0:
                #self.data.id = new_acc_nums[0]
                #Use the FIRST accession as the ID, not the first on this line!
                self.data.id = self.data.annotations['accessions'][0]

    def wgs(self, content):
        self.data.annotations['wgs'] = content.split('-')

    def add_wgs_scafld(self, content):
        self.data.annotations.setdefault('wgs_scafld',[]).append(content.split('-'))

    def nid(self, content):
        self.data.annotations['nid'] = content

    def pid(self, content):
        self.data.annotations['pid'] = content

    def version(self, version_id):
        #Want to use the versioned accession as the record.id
        #This comes from the VERSION line in GenBank files, or the
        #obsolete SV line in EMBL.  For the new EMBL files we need
        #both the version suffix from the ID line and the accession
        #from the AC line.
        if version_id.count(".")==1 and version_id.split(".")[1].isdigit():
            self.accession(version_id.split(".")[0])
            self.version_suffix(version_id.split(".")[1])
        else:
            #For backwards compatibility...
            self.data.id = version_id

    def project(self, content):
        """Handle the information from the PROJECT line as a list of projects.

        e.g.
        PROJECT     GenomeProject:28471

        or:
        PROJECT     GenomeProject:13543  GenomeProject:99999

        This is stored as dbxrefs in the SeqRecord to be consistent with the
        projected switch of this line to DBLINK in future GenBank versions.
        Note the NCBI plan to replace "GenomeProject:28471" with the shorter
        "Project:28471" as part of this transition.
        """
        content = content.replace("GenomeProject:", "Project:")
        self.data.dbxrefs.extend([p for p in content.split() if p])

    def dblink(self, content):
        """Store DBLINK cross references as dbxrefs in our record object.

        This line type is expected to replace the PROJECT line in 2009. e.g.

        During transition:
        
        PROJECT     GenomeProject:28471
        DBLINK      Project:28471
                    Trace Assembly Archive:123456

        Once the project line is dropped:

        DBLINK      Project:28471
                    Trace Assembly Archive:123456

        Note GenomeProject -> Project.

        We'll have to see some real examples to be sure, but based on the
        above example we can expect one reference per line.

        Note that at some point the NCBI have included an extra space, e.g.

        DBLINK      Project: 28471
        """
        #During the transition period with both PROJECT and DBLINK lines,
        #we don't want to add the same cross reference twice.
        while ": " in content:
            content = content.replace(": ", ":")
        if content.strip() not in self.data.dbxrefs:
            self.data.dbxrefs.append(content.strip())

    def version_suffix(self, version):
        """Set the version to overwrite the id.

        Since the verison provides the same information as the accession
        number, plus some extra info, we set this as the id if we have
        a version.
        """
        #e.g. GenBank line:
        #VERSION     U49845.1  GI:1293613
        #or the obsolete EMBL line:
        #SV   U49845.1
        #Scanner calls consumer.version("U49845.1")
        #which then calls consumer.version_suffix(1)
        #
        #e.g. EMBL new line:
        #ID   X56734; SV 1; linear; mRNA; STD; PLN; 1859 BP.
        #Scanner calls consumer.version_suffix(1)
        assert version.isdigit()
        self.data.annotations['sequence_version'] = int(version)

    def db_source(self, content):
        self.data.annotations['db_source'] = content.rstrip()

    def gi(self, content):
        self.data.annotations['gi'] = content

    def keywords(self, content):
        self.data.annotations['keywords'] = self._split_keywords(content)

    def segment(self, content):
        self.data.annotations['segment'] = content

    def source(self, content):
        #Note that some software (e.g. VectorNTI) may produce an empty
        #source (rather than using a dot/period as might be expected).
        if content == "":
            source_info = ""
        elif content[-1] == '.':
            source_info = content[:-1]
        else:
            source_info = content
        self.data.annotations['source'] = source_info

    def organism(self, content):
        self.data.annotations['organism'] = content

    def taxonomy(self, content):
        """Records (another line of) the taxonomy lineage.
        """
        lineage = self._split_taxonomy(content)
        try:
            self.data.annotations['taxonomy'].extend(lineage)
        except KeyError:
            self.data.annotations['taxonomy'] = lineage
        
    def reference_num(self, content):
        """Signal the beginning of a new reference object.
        """
        # if we have a current reference that hasn't been added to
        # the list of references, add it.
        if self._cur_reference is not None:
            self.data.annotations['references'].append(self._cur_reference)
        else:
            self.data.annotations['references'] = []

        self._cur_reference = SeqFeature.Reference()

    def reference_bases(self, content):
        """Attempt to determine the sequence region the reference entails.

        Possible types of information we may have to deal with:
        
        (bases 1 to 86436)
        (sites)
        (bases 1 to 105654; 110423 to 111122)
        1  (residues 1 to 182)
        """
        # first remove the parentheses or other junk
        ref_base_info = content[1:-1]

        all_locations = []
        # parse if we've got 'bases' and 'to'
        if ref_base_info.find('bases') != -1 and \
            ref_base_info.find('to') != -1:
            # get rid of the beginning 'bases'
            ref_base_info = ref_base_info[5:]
            locations = self._split_reference_locations(ref_base_info)
            all_locations.extend(locations)
        elif (ref_base_info.find("residues") >= 0 and
              ref_base_info.find("to") >= 0):
            residues_start = ref_base_info.find("residues")
            # get only the information after "residues"
            ref_base_info = ref_base_info[(residues_start + len("residues ")):]
            locations = self._split_reference_locations(ref_base_info)
            all_locations.extend(locations)

        # make sure if we are not finding information then we have
        # the string 'sites' or the string 'bases'
        elif (ref_base_info == 'sites' or
              ref_base_info.strip() == 'bases'):
            pass
        # otherwise raise an error
        else:
            raise ValueError("Could not parse base info %s in record %s" %
                             (ref_base_info, self.data.id))

        self._cur_reference.location = all_locations

    def _split_reference_locations(self, location_string):
        """Get reference locations out of a string of reference information
        
        The passed string should be of the form:

            1 to 20; 20 to 100

        This splits the information out and returns a list of location objects
        based on the reference locations.
        """
        # split possibly multiple locations using the ';'
        all_base_info = location_string.split(';')

        new_locations = []
        for base_info in all_base_info:
            start, end = base_info.split('to')
            new_start, new_end = \
              self._convert_to_python_numbers(int(start.strip()),
                                              int(end.strip()))
            this_location = SeqFeature.FeatureLocation(new_start, new_end)
            new_locations.append(this_location)
        return new_locations

    def authors(self, content):
        if self._cur_reference.authors:
            self._cur_reference.authors += ' ' + content
        else:
            self._cur_reference.authors = content

    def consrtm(self, content):
        if self._cur_reference.consrtm:
            self._cur_reference.consrtm += ' ' + content
        else:
            self._cur_reference.consrtm = content

    def title(self, content):
        if self._cur_reference is None:
            import warnings
            from Bio import BiopythonParserWarning
            warnings.warn("GenBank TITLE line without REFERENCE line.",
                          BiopythonParserWarning)
        elif self._cur_reference.title:
            self._cur_reference.title += ' ' + content
        else:
            self._cur_reference.title = content

    def journal(self, content):
        if self._cur_reference.journal:
            self._cur_reference.journal += ' ' + content
        else:
            self._cur_reference.journal = content

    def medline_id(self, content):
        self._cur_reference.medline_id = content

    def pubmed_id(self, content):
        self._cur_reference.pubmed_id = content

    def remark(self, content):
        """Deal with a reference comment."""
        if self._cur_reference.comment:
            self._cur_reference.comment += ' ' + content
        else:
            self._cur_reference.comment = content

    def comment(self, content):
        try:
            self.data.annotations['comment'] += "\n" + "\n".join(content)
        except KeyError:
            self.data.annotations['comment'] = "\n".join(content)

    def features_line(self, content):
        """Get ready for the feature table when we reach the FEATURE line.
        """
        self.start_feature_table()

    def start_feature_table(self):
        """Indicate we've got to the start of the feature table.
        """
        # make sure we've added on our last reference object
        if self._cur_reference is not None:
            self.data.annotations['references'].append(self._cur_reference)
            self._cur_reference = None

    def feature_key(self, content):
        # start a new feature
        self._cur_feature = SeqFeature.SeqFeature()
        self._cur_feature.type = content
        self.data.features.append(self._cur_feature)

    def location(self, content):
        """Parse out location information from the location string.

        This uses simple Python code with some regular expressions to do the
        parsing, and then translates the results into appropriate objects.
        """
        # clean up newlines and other whitespace inside the location before
        # parsing - locations should have no whitespace whatsoever
        location_line = self._clean_location(content)

        # Older records have junk like replace(266,"c") in the
        # location line. Newer records just replace this with
        # the number 266 and have the information in a more reasonable
        # place. So we'll just grab out the number and feed this to the
        # parser. We shouldn't really be losing any info this way.
        if location_line.find('replace') != -1:
            comma_pos = location_line.find(',')
            location_line = location_line[8:comma_pos]
        
        cur_feature = self._cur_feature

        #Handle top level complement here for speed
        if location_line.startswith("complement("):
            assert location_line.endswith(")")
            location_line = location_line[11:-1]
            strand = -1
        elif self._seq_type.find("DNA") >= 0 \
        or self._seq_type.find("RNA") >= 0:
            #Nucleotide
            strand = 1
        else:
            #Protein
            strand = None

        #Special case handling of the most common cases for speed
        if _re_simple_location.match(location_line):
            #e.g. "123..456"
            s, e = location_line.split("..")
            cur_feature.location = SeqFeature.FeatureLocation(int(s)-1,
                                                              int(e),
                                                              strand)
            return
        if _re_simple_compound.match(location_line):
            #e.g. join(<123..456,480..>500)
            i = location_line.find("(")
            cur_feature.location_operator = location_line[:i]
            #we can split on the comma because these are simple locations
            for part in location_line[i+1:-1].split(","):
                s, e = part.split("..")
                f = SeqFeature.SeqFeature(SeqFeature.FeatureLocation(int(s)-1,
                                                                     int(e),
                                                                     strand),
                        location_operator=cur_feature.location_operator,
                        type=cur_feature.type)
                cur_feature.sub_features.append(f)
            s = cur_feature.sub_features[0].location.start
            e = cur_feature.sub_features[-1].location.end
            cur_feature.location = SeqFeature.FeatureLocation(s,e, strand)
            return
        
        #Handle the general case with more complex regular expressions
        if _re_complex_location.match(location_line):
            #e.g. "AL121804.2:41..610"
            if ":" in location_line:
                location_ref, location_line = location_line.split(":")
                cur_feature.location = _loc(location_line, self._expected_size, strand)
                cur_feature.location.ref = location_ref
            else:
                cur_feature.location = _loc(location_line, self._expected_size, strand)
            return
        if _re_complex_compound.match(location_line):
            i = location_line.find("(")
            cur_feature.location_operator = location_line[:i]
            #Can't split on the comma because of ositions like one-of(1,2,3)
            for part in _split_compound_loc(location_line[i+1:-1]):
                if part.startswith("complement("):
                    assert part[-1]==")"
                    part = part[11:-1]
                    assert strand != -1, "Double complement?"
                    part_strand = -1
                else:
                    part_strand = strand
                if ":" in part:
                    ref, part = part.split(":")
                else:
                    ref = None
                try:
                    loc = _loc(part, self._expected_size, part_strand)
                except ValueError, err:
                    print location_line
                    print part
                    raise err
                f = SeqFeature.SeqFeature(location=loc, ref=ref,
                        location_operator=cur_feature.location_operator,
                        type=cur_feature.type)
                cur_feature.sub_features.append(f)
            # Historically a join on the reverse strand has been represented
            # in Biopython with both the parent SeqFeature and its children
            # (the exons for a CDS) all given a strand of -1.  Likewise, for
            # a join feature on the forward strand they all have strand +1.
            # However, we must also consider evil mixed strand examples like
            # this, join(complement(69611..69724),139856..140087,140625..140650)
            strands = set(sf.strand for sf in cur_feature.sub_features)
            if len(strands)==1:
                strand = cur_feature.sub_features[0].strand
            else:
                strand = None # i.e. mixed strands
            s = cur_feature.sub_features[0].location.start
            e = cur_feature.sub_features[-1].location.end
            cur_feature.location = SeqFeature.FeatureLocation(s, e, strand)
            return
        #Not recognised
        if "order" in location_line and "join" in location_line:
            #See Bug 3197
            msg = 'Combinations of "join" and "order" within the same ' + \
                  'location (nested operators) are illegal:\n' + location_line
            raise LocationParserError(msg)
        raise LocationParserError(location_line)

    def feature_qualifier(self, key, value):
        """When we get a qualifier key and its value.
        
        Can receive None, since you can have valueless keys such as /pseudo
        """
        # Hack to try to preserve historical behaviour of /pseudo etc
        if value is None:
            if key not in self._cur_feature.qualifiers:
                self._cur_feature.qualifiers[key] = [""]
                return
            
        value = value.replace('"', '')
        if self._feature_cleaner is not None:
            value = self._feature_cleaner.clean_value(key, value)

        # if the qualifier name exists, append the value
        if key in self._cur_feature.qualifiers:
            self._cur_feature.qualifiers[key].append(value)
        # otherwise start a new list of the key with its values
        else:
            self._cur_feature.qualifiers[key] = [value]
       
    def feature_qualifier_name(self, content_list):
        """Use feature_qualifier instead (OBSOLETE)."""
        raise NotImplementedError("Use the feature_qualifier method instead.")

    def feature_qualifier_description(self, content):
        """Use feature_qualifier instead (OBSOLETE)."""
        raise NotImplementedError("Use the feature_qualifier method instead.")

    def contig_location(self, content):
        """Deal with CONTIG information."""
        #Historically this was stored as a SeqFeature object, but it was
        #stored under record.annotations["contig"] and not under
        #record.features with the other SeqFeature objects.
        #
        #The CONTIG location line can include additional tokens like
        #Gap(), Gap(100) or Gap(unk100) which are not used in the feature
        #location lines, so storing it using SeqFeature based location
        #objects is difficult.
        #
        #We now store this a string, which means for BioSQL we are now in
        #much better agreement with how BioPerl records the CONTIG line
        #in the database.
        #
        #NOTE - This code assumes the scanner will return all the CONTIG
        #lines already combined into one long string!
        self.data.annotations["contig"] = content

    def origin_name(self, content):
        pass

    def base_count(self, content):
        pass

    def base_number(self, content):
        pass

    def sequence(self, content):
        """Add up sequence information as we get it.

        To try and make things speedier, this puts all of the strings
        into a list of strings, and then uses string.join later to put
        them together. Supposedly, this is a big time savings
        """
        assert ' ' not in content
        self._seq_data.append(content.upper())

    def record_end(self, content):
        """Clean up when we've finished the record.
        """
        from Bio import Alphabet
        from Bio.Alphabet import IUPAC
        from Bio.Seq import Seq, UnknownSeq

        #Try and append the version number to the accession for the full id
        if self.data.id is None:
            assert 'accessions' not in self.data.annotations, \
                   self.data.annotations['accessions']
            self.data.id = self.data.name #Good fall back?
        elif self.data.id.count('.') == 0:
            try:
                self.data.id+='.%i' % self.data.annotations['sequence_version']
            except KeyError:
                pass
        
        # add the sequence information
        # first, determine the alphabet
        # we default to an generic alphabet if we don't have a
        # seq type or have strange sequence information.
        seq_alphabet = Alphabet.generic_alphabet

        # now set the sequence
        sequence = "".join(self._seq_data)

        if self._expected_size is not None \
        and len(sequence) != 0 \
        and self._expected_size != len(sequence):
            import warnings
            from Bio import BiopythonParserWarning
            warnings.warn("Expected sequence length %i, found %i (%s)." \
                          % (self._expected_size, len(sequence), self.data.id),
                          BiopythonParserWarning)

        if self._seq_type:
            # mRNA is really also DNA, since it is actually cDNA
            if self._seq_type.find('DNA') != -1 or \
               self._seq_type.find('mRNA') != -1:
                seq_alphabet = IUPAC.ambiguous_dna
            # are there ever really RNA sequences in GenBank?
            elif self._seq_type.find('RNA') != -1:
                #Even for data which was from RNA, the sequence string
                #is usually given as DNA (T not U).  Bug 2408
                if "T" in sequence and "U" not in sequence:
                    seq_alphabet = IUPAC.ambiguous_dna
                else:
                    seq_alphabet = IUPAC.ambiguous_rna
            elif self._seq_type.upper().find('PROTEIN') != -1:
                seq_alphabet = IUPAC.protein  # or extended protein?
            # work around ugly GenBank records which have circular or
            # linear but no indication of sequence type
            elif self._seq_type in ["circular", "linear", "unspecified"]:
                pass
            # we have a bug if we get here
            else:
                raise ValueError("Could not determine alphabet for seq_type %s"
                                 % self._seq_type)

        if not sequence and self.__expected_size:
            self.data.seq = UnknownSeq(self._expected_size, seq_alphabet)
        else:
            self.data.seq = Seq(sequence, seq_alphabet)

class _RecordConsumer(_BaseGenBankConsumer):
    """Create a GenBank Record object from scanner generated information (PRIVATE).
    """
    def __init__(self):
        _BaseGenBankConsumer.__init__(self)
        import Record
        self.data = Record.Record()

        self._seq_data = []
        self._cur_reference = None
        self._cur_feature = None
        self._cur_qualifier = None
        
    def wgs(self, content):
        self.data.wgs = content.split('-')

    def add_wgs_scafld(self, content):
        self.data.wgs_scafld.append(content.split('-'))

    def locus(self, content):
        self.data.locus = content

    def size(self, content):
        self.data.size = content

    def residue_type(self, content):
        self.data.residue_type = content

    def data_file_division(self, content):
        self.data.data_file_division = content

    def date(self, content):
        self.data.date = content

    def definition(self, content):
        self.data.definition = content

    def accession(self, content):
        for acc in self._split_accessions(content):
            if acc not in self.data.accession:
                self.data.accession.append(acc)

    def nid(self, content):
        self.data.nid = content

    def pid(self, content):
        self.data.pid = content

    def version(self, content):
        self.data.version = content

    def db_source(self, content):
        self.data.db_source = content.rstrip()

    def gi(self, content):
        self.data.gi = content

    def keywords(self, content):
        self.data.keywords = self._split_keywords(content)

    def project(self, content):
        self.data.projects.extend([p for p in content.split() if p])

    def dblink(self, content):
        self.data.dblinks.append(content)

    def segment(self, content):
        self.data.segment = content

    def source(self, content):
        self.data.source = content

    def organism(self, content):
        self.data.organism = content

    def taxonomy(self, content):
        self.data.taxonomy = self._split_taxonomy(content)

    def reference_num(self, content):
        """Grab the reference number and signal the start of a new reference.
        """
        # check if we have a reference to add
        if self._cur_reference is not None:
            self.data.references.append(self._cur_reference)

        import Record
        self._cur_reference = Record.Reference()
        self._cur_reference.number = content

    def reference_bases(self, content):
        self._cur_reference.bases = content

    def authors(self, content):
        self._cur_reference.authors = content

    def consrtm(self, content):
        self._cur_reference.consrtm = content

    def title(self, content):
        if self._cur_reference is None:
            import warnings
            from Bio import BiopythonParserWarning
            warnings.warn("GenBank TITLE line without REFERENCE line.",
                          BiopythonParserWarning)
            return
        self._cur_reference.title = content

    def journal(self, content):
        self._cur_reference.journal = content

    def medline_id(self, content):
        self._cur_reference.medline_id = content
        
    def pubmed_id(self, content):
        self._cur_reference.pubmed_id = content

    def remark(self, content):
        self._cur_reference.remark = content
        
    def comment(self, content):
        self.data.comment += "\n".join(content)

    def primary_ref_line(self,content):
        """Data for the PRIMARY line"""
        self.data.primary.append(content)

    def primary(self,content):
        pass
    
    def features_line(self, content):
        """Get ready for the feature table when we reach the FEATURE line.
        """
        self.start_feature_table()

    def start_feature_table(self):
        """Signal the start of the feature table.
        """
        # we need to add on the last reference
        if self._cur_reference is not None:
            self.data.references.append(self._cur_reference)

    def feature_key(self, content):
        """Grab the key of the feature and signal the start of a new feature.
        """
        # first add on feature information if we've got any
        self._add_feature()

        import Record
        self._cur_feature = Record.Feature()
        self._cur_feature.key = content

    def _add_feature(self):
        """Utility function to add a feature to the Record.

        This does all of the appropriate checking to make sure we haven't
        left any info behind, and that we are only adding info if it
        exists.
        """
        if self._cur_feature is not None:
            # if we have a left over qualifier, add it to the qualifiers
            # on the current feature
            if self._cur_qualifier is not None:
                self._cur_feature.qualifiers.append(self._cur_qualifier)

            self._cur_qualifier = None
            self.data.features.append(self._cur_feature)

    def location(self, content):
        self._cur_feature.location = self._clean_location(content)

    def feature_qualifier(self, key, value):
        self.feature_qualifier_name([key])
        if value is not None:
            self.feature_qualifier_description(value)

    def feature_qualifier_name(self, content_list):
        """Deal with qualifier names
        
        We receive a list of keys, since you can have valueless keys such as
        /pseudo which would be passed in with the next key (since no other
        tags separate them in the file)
        """
        import Record
        for content in content_list:
            # the record parser keeps the /s -- add them if we don't have 'em
            if content.find("/") != 0:
                content = "/%s" % content
            # add on a qualifier if we've got one
            if self._cur_qualifier is not None:
                self._cur_feature.qualifiers.append(self._cur_qualifier)

            self._cur_qualifier = Record.Qualifier()
            self._cur_qualifier.key = content

    def feature_qualifier_description(self, content):
        # if we have info then the qualifier key should have a ='s
        if self._cur_qualifier.key.find("=") == -1:
            self._cur_qualifier.key = "%s=" % self._cur_qualifier.key
        cur_content = self._remove_newlines(content)
        # remove all spaces from the value if it is a type where spaces
        # are not important
        for remove_space_key in self.__class__.remove_space_keys:
            if self._cur_qualifier.key.find(remove_space_key) >= 0:
                cur_content = self._remove_spaces(cur_content)
        self._cur_qualifier.value = self._normalize_spaces(cur_content)

    def base_count(self, content):
        self.data.base_counts = content

    def origin_name(self, content):
        self.data.origin = content

    def contig_location(self, content):
        """Signal that we have contig information to add to the record.
        """
        self.data.contig = self._clean_location(content) 

    def sequence(self, content):
        """Add sequence information to a list of sequence strings.

        This removes spaces in the data and uppercases the sequence, and
        then adds it to a list of sequences. Later on we'll join this
        list together to make the final sequence. This is faster than
        adding on the new string every time.
        """
        assert ' ' not in content
        self._seq_data.append(content.upper())

    def record_end(self, content):
        """Signal the end of the record and do any necessary clean-up.
        """
        # add together all of the sequence parts to create the
        # final sequence string
        self.data.sequence = "".join(self._seq_data)
        # add on the last feature
        self._add_feature()


def parse(handle):
    """Iterate over GenBank formatted entries as Record objects.

    >>> from Bio import GenBank
    >>> handle = open("GenBank/NC_000932.gb")
    >>> for record in GenBank.parse(handle):
    ...     print record.accession
    ['NC_000932']
    >>> handle.close()

    To get SeqRecord objects use Bio.SeqIO.parse(..., format="gb")
    instead.
    """
    return iter(Iterator(handle, RecordParser()))

def read(handle):
    """Read a handle containing a single GenBank entry as a Record object.

    >>> from Bio import GenBank
    >>> handle = open("GenBank/NC_000932.gb")
    >>> record = GenBank.read(handle)
    >>> print record.accession
    ['NC_000932']
    >>> handle.close()
                       
    To get a SeqRecord object use Bio.SeqIO.read(..., format="gb")
    instead.
    """
    iterator = parse(handle)
    try:
        first = iterator.next()
    except StopIteration:
        first = None
    if first is None:
        raise ValueError("No records found in handle")
    try:
        second = iterator.next()
    except StopIteration:
        second = None
    if second is not None:
        raise ValueError("More than one record found in handle")
    return first

def _test():
    """Run the Bio.GenBank module's doctests."""
    import doctest
    import os
    if os.path.isdir(os.path.join("..","..","Tests")):
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join("..","..","Tests"))
        doctest.testmod()
        os.chdir(cur_dir)
        del cur_dir
        print "Done"
    elif os.path.isdir(os.path.join("Tests")):
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join("Tests"))
        doctest.testmod()
        os.chdir(cur_dir)
        del cur_dir
        print "Done"

if __name__ == "__main__":
    _test()
