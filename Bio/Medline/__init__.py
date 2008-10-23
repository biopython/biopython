# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
This module provides code to work with Medline.

Classes:
Record           A dictionary holding Medline data.

Functions:
read             Reads one Medline record
parse            Allows you to iterate over a bunch of Medline records

Deprecated classes:
Iterator         Iterates over a file containing Medline records.
RecordParser     Parses a Medline record into a Record object.

_Scanner         Scans a Medline record.
_RecordConsumer  Consumes Medline data to a Record object.

"""

class Record(dict):
    """A dictionary holding information from a Medline record.
    All data are stored under the mnemonic appearing in the Medline
    file. These mnemonics have the following interpretations:

    Mnemonic  Description
    AB        Abstract
    CI        Copyright Information
    AD        Affiliation
    IRAD      Investigator Affiliation
    AID       Article Identifier
    AU        Author
    FAU       Full Author
    CN        Corporate Author
    DCOM      Date Completed
    DA        Date Created
    LR        Date Last Revised
    DEP       Date of Electronic Publication
    DP        Date of Publication
    EDAT      Entrez Date
    GS        Gene Symbol
    GN        General Note
    GR        Grant Number
    IR        Investigator Name
    FIR       Full Investigator Name
    IS        ISSN
    IP        Issue
    TA        Journal Title Abbreviation
    JT        Journal Title
    LA        Language
    LID       Location Identifier
    MID       Manuscript Identifier
    MHDA      MeSH Date
    MH        MeSH Terms
    JID       NLM Unique ID
    RF        Number of References
    OAB       Other Abstract
    OCI       Other Copyright Information
    OID       Other ID
    OT        Other Term
    OTO       Other Term Owner
    OWN       Owner
    PG        Pagination
    PS        Personal Name as Subject
    FPS       Full Personal Name as Subject
    PL        Place of Publication
    PHST      Publication History Status
    PST       Publication Status
    PT        Publication Type
    PUBM      Publishing Model
    PMC       PubMed Central Identifier
    PMID      PubMed Unique Identifier
    RN        Registry Number/EC Number
    NM        Substance Name
    SI        Secondary Source ID
    SO        Source
    SFM       Space Flight Mission
    STAT      Status
    SB        Subset
    TI        Title
    TT        Transliterated Title
    VI        Volume
    CON       Comment on
    CIN       Comment in
    EIN       Erratum in
    EFR       Erratum for
    CRI       Corrected and Republished in
    CRF       Corrected and Republished from
    PRIN      Partial retraction in
    PROF      Partial retraction of
    RPI       Republished in
    RPF       Republished from
    RIN       Retraction in
    ROF       Retraction of
    UIN       Update in
    UOF       Update of
    SPIN      Summary for patients in
    ORI       Original report in
    """
    def __init__(self):
        # The __init__ function can be removed when we remove the old parser
        self.id = ''
        self.pubmed_id = ''

        self.mesh_headings = []
        self.mesh_tree_numbers = []
        self.mesh_subheadings = []

        self.abstract = ''
        self.comments = []
        self.abstract_author = ''
        self.english_abstract = ''

        self.source = ''
        self.publication_types = []
        self.number_of_references = ''

        self.authors = []
        self.no_author = ''
        self.address = ''

        self.journal_title_code = ''
        self.title_abbreviation = ''
        self.issn = ''
        self.journal_subsets = []
        self.country = ''
        self.languages = []

        self.title = ''
        self.transliterated_title = ''
        self.call_number = ''
        self.issue_part_supplement = ''
        self.volume_issue = ''
        self.publication_date = ''
        self.year = ''
        self.pagination = ''

        self.special_list = ''

        self.substance_name = ''
        self.gene_symbols = []
        self.secondary_source_ids = []
        self.identifications = []
        self.registry_numbers = []

        self.personal_name_as_subjects = []

        self.record_originators = []
        self.entry_date = ''
        self.entry_month = ''
        self.class_update_date = ''
        self.last_revision_date = ''
        self.major_revision_date = ''

        self.undefined = []


def parse(handle):
    """Read Medline records one by one from the handle.

    The handle is either is a Medline file, a file-like object, or a list
    of lines describing one or more Medline records.

    Typical usage:

        from Bio import Medline
        handle = open("mymedlinefile")
        records = Medline.parse(handle)
        for record in record:
            print record['TI']

    """
    # These keys point to string values
    textkeys = ("ID", "PMID", "SO", "RF", "NI", "JC", "TA", "IS", "CY", "TT",
                "CA", "IP", "VI", "DP", "YR", "PG", "LID", "DA", "LR", "OWN",
                "STAT", "DCOM", "PUBM", "DEP", "PL", "JID", "SB", "PMC",
                "EDAT", "MHDA", "PST", "AB", "AD", "EA", "TI", "JT")
    handle = iter(handle)
    # First skip blank lines
    for line in handle:
        line = line.rstrip()
        if line:
            break
    else:
        return
    record = Record()
    finished = False
    while not finished:
        if line[:6]=="      ": # continuation line
            record[key].append(line[6:])
        elif line:
            key = line[:4].rstrip()
            if not key in record:
                record[key] = []
            record[key].append(line[6:])
        try:
            line = handle.next()
        except StopIteration:
            finished = True
        else:
            line = line.rstrip()
            if line:
                continue
        # Join each list of strings into one string.
        for key in textkeys:
            if key in record:
                record[key] = " ".join(record[key])
        if record:
            yield record
        record = Record()

def read(handle):
    """Read a single Medline records from the handle.

    The handle is either is a Medline file, a file-like object, or a list
    of lines describing a Medline record.

    Typical usage:

        from Bio import Medline
        handle = open("mymedlinefile")
        record = Medline.read(handle)
        print record['TI']

    """
    records = parse(handle)
    return records.next()

### Everything below is deprecated

from Bio import File
from Bio.ParserSupport import *

class Iterator:
    """Returns one record at a time from a file of Medline records.

    Methods:
    next   Return the next record from the stream, or None.

    """
    def __init__(self, handle, parser=None):
        """__init__(self, handle, parser=None)

        Create a new iterator.  handle is a file-like object.  parser
        is an optional Parser object to change the results into another form.
        If set to None, then the raw contents of the file will be returned.

        """
        import warnings
        warnings.warn("Bio.Medline.Iterator is deprecated. Instead of Bio.Medline.Iterator(handle, Bio.Medline.RecordParser()), please use Bio.Medline.parse(handle)", DeprecationWarning)
        self._handle = handle
        self._parser = parser

    def __iter__(self):
        return self

    def next(self):
        """next(self) -> object

        Return the next medline record from the file.  If no more records,
        return None.

        """
        lines = []
        for line in self._handle:
            lines.append(line)
            if line.strip()=='':
                break
        else:
            raise StopIteration
            
        data = ''.join(lines)

        if self._parser is not None:
            return self._parser.parse_str(data)
        return data

class RecordParser(AbstractParser):
    """Parses Medline data into a Record object.

    """
    def __init__(self):
        import warnings
        warnings.warn("Bio.Medline.RecordParser is deprecated. Instead of Bio.Medline.RecordParser().parse(handle)), please use Bio.Medline.read(handle)", DeprecationWarning)
        self._scanner = _Scanner()
        self._consumer = _RecordConsumer()

    def parse(self, handle):
        self._scanner.feed(handle, self._consumer)
        return self._consumer.data

class _Scanner:
    """Scans a Medline record.
    
    """
    # map the category qualifier to an event
    _categories = {
        "AA" : "abstract_author",
        "AB" : "abstract",
        "AD" : "address",
        "AU" : "author",
        "CA" : "call_number",
        "CM" : "comments",
        "CU" : "class_update_date",
        "CY" : "country",
        "DA" : "entry_date",
        "DP" : "publication_date",
        "EA" : "english_abstract",
        "EM" : "entry_month",
        "GS" : "gene_symbol",
        "ID" : "identification",
        "IP" : "issue_part_supplement",
        "IS" : "issn",
        "JC" : "journal_title_code",
        "LA" : "language",
        "LI" : "special_list",
        "LR" : "last_revision_date",
        "MH" : "mesh_heading",
        "MN" : "mesh_tree_number",
        "MR" : "major_revision_date",
        "NI" : "no_author",
        "NM" : "substance_name",
        "PG" : "pagination",
        "PS" : "personal_name_as_subject",
        "PT" : "publication_type",
        "RF" : "number_of_references",
        "RN" : "cas_registry_number",
        "RO" : "record_originator",
        "SB" : "journal_subset",
        "SH" : "subheadings",
        "SI" : "secondary_source_id",
        "SO" : "source",
        "TA" : "title_abbreviation",
        "TI" : "title",
        "TT" : "transliterated_title",
        "UI" : "unique_identifier",
        "VI" : "volume_issue",
        "YR" : "year",

        # Not documented.
        "PMID" : "pubmed_id",
        }

    def feed(self, handle, consumer):
        """feed(self, handle, consumer)

        Feed in a Medline unit record for scanning.  handle is a file-like
        object that contains a Medline record.  consumer is a
        Consumer object that will receive events as the report is scanned.

        """
        if isinstance(handle, File.UndoHandle):
            uhandle = handle
        else:
            uhandle = File.UndoHandle(handle)

        # Read the Entrez header information, if it exists
        if attempt_read_and_call(uhandle, consumer.noevent, start='Entrez'):
            read_and_call(uhandle, consumer.noevent, start='----------------')
        self._scan_record(uhandle, consumer)

    def _scan_record(self, uhandle, consumer):
        consumer.start_record()

        prev_qualifier = None
        while 1:
            line = uhandle.readline()
            if is_blank_line(line):
                break

            # There are 2 possible formats for a line:
            # TI  - Epidemiology of mycobacterial resistance (especially Mycoba
            #       tuberculosis).
            # 1) qualifier + '-' + data
            # 2) continuation, with just data

            # Check to see if it's a continuation line.
            qualifier = line[:4].rstrip()
            # There's a bug in some MH lines where the "isolation &
            # purification" subheading gets split across lines and
            # purification at the beginning of the line, with only 1
            # space.
            if line[0] == '\t' or qualifier == '' or \
               line[:13] == ' purification':
                if prev_qualifier is None:
                    raise ValueError("Continuation on first line\n%s" % line)
                qualifier = prev_qualifier
            else:
                # Make sure it contains a '-'
                if len(line) < 5 or line[4] != '-':
                    raise ValueError(\
                        "I don't understand the format of line %s" % line)
            prev_qualifier = qualifier
            
            try:
                fn = getattr(consumer, self._categories[qualifier])
            except KeyError:
                # call an 'undefined' function for 
                consumer.undefined(line)
            else:
                fn(line)
        
        consumer.end_record()

class _RecordConsumer(AbstractConsumer):
    """Consumer that converts a Medline record to a Record object.

    Members:
    data    Record with Medline data.

    """
    def __init__(self):
        self.data = None

    def start_record(self):
        self.data = Record()

    def end_record(self):
        self._clean_record(self.data)

    def abstract_author(self, line):
        self.data.abstract_author = self._clean(line)
    
    def abstract(self, line):
        self.data.abstract = self.data.abstract + self._clean(line, rstrip=0)
    
    def address(self, line):
        self.data.address = self.data.address + self._clean(line, rstrip=0)
    
    def author(self, line):
        self.data.authors.append(self._clean(line))
    
    def call_number(self, line):
        assert not self.data.call_number, "call numbers already defined"
        self.data.call_number = self._clean(line)
    
    def comments(self, line):
        self.data.comments.append(self._clean(line))
    
    def class_update_date(self, line):
        assert not self.data.class_update_date, \
               "class update date already defined"
        self.data.class_update_date = self._clean(line)
    
    def country(self, line):
        assert not self.data.country, "country already defined"
        self.data.country = self._clean(line)
    
    def entry_date(self, line):
        assert not self.data.entry_date, "entry date already defined"
        self.data.entry_date = self._clean(line)
    
    def publication_date(self, line):
        assert not self.data.publication_date, \
               "publication date already defined"
        self.data.publication_date = self._clean(line)
    
    def english_abstract(self, line):
        assert not self.data.english_abstract, \
               "english abstract already defined"
        self.data.english_abstract = self._clean(line)
    
    def entry_month(self, line):
        assert not self.data.entry_month, \
               "entry month already defined"
        self.data.entry_month = self._clean(line)
    
    def gene_symbol(self, line):
        self.data.gene_symbols.append(self._clean(line))
    
    def identification(self, line):
        self.data.identifications.append(self._clean(line))
    
    def issue_part_supplement(self, line):
        assert not self.data.issue_part_supplement, \
               "issue/part/supplement already defined"
        self.data.issue_part_supplement = self._clean(line)
    
    def issn(self, line):
        assert not self.data.issn, "ISSN already defined"
        self.data.issn = self._clean(line)
    
    def journal_title_code(self, line):
        assert not self.data.journal_title_code, \
               "journal title code already defined"
        self.data.journal_title_code = self._clean(line)
    
    def language(self, line):
        self.data.languages.append(self._clean(line))
    
    def special_list(self, line):
        assert not self.data.special_list, "special list already defined"
        self.data.special_list = self._clean(line)
    
    def last_revision_date(self, line):
        assert not self.data.last_revision_date, \
               "last revision date already defined"
        self.data.last_revision_date = self._clean(line)
    
    def mesh_heading(self, line):
        # Check to see whether this is a new MH line, or a
        # continuation of an old one.  If it's a continuation of an
        # old one, append it to the previous line.
        # See PMID 12107064 for an example, found by Dan Rubin.
        if line[:2] == 'MH':
            self.data.mesh_headings.append(self._clean(line))
        else:
            prev_mh = self.data.mesh_headings.pop()
            continued_mh = self._clean(line)
            self.data.mesh_headings.append("%s %s" % (prev_mh, continued_mh))
    
    def mesh_tree_number(self, line):
        self.data.mesh_tree_numbers.append(self._clean(line))
    
    def major_revision_date(self, line):
        assert not self.data.major_revision_date, \
               "major revision date already defined"
        self.data.major_revision_date = self._clean(line)
    
    def no_author(self, line):
        assert not self.data.no_author, "no author already defined"
        self.data.no_author = self._clean(line)
    
    def substance_name(self, line):
        assert not self.data.substance_name, "substance name already defined"
        self.data.substance_name = self._clean(line)
    
    def pagination(self, line):
        assert not self.data.pagination, "pagination already defined"
        self.data.pagination = self._clean(line)
    
    def personal_name_as_subject(self, line):
        self.data.personal_name_as_subjects.append(self._clean(line))
    
    def publication_type(self, line):
        self.data.publication_types.append(self._clean(line))
    
    def number_of_references(self, line):
        assert not self.data.number_of_references, \
               "num of references already defined"
        self.data.number_of_references = self._clean(line)
    
    def cas_registry_number(self, line):
        self.data.registry_numbers.append(self._clean(line))
    
    def record_originator(self, line):
        self.data.record_originators.append(self._clean(line))
    
    def journal_subset(self, line):
        self.data.journal_subsets.append(self._clean(line))
    
    def subheadings(self, line):
        self.data.mesh_subheadings.append(self._clean(line))
    
    def secondary_source_id(self, line):
        self.data.secondary_source_ids.append(self._clean(line))
    
    def source(self, line):
        self.data.source = self.data.source + self._clean(line, rstrip=0)
    
    def title_abbreviation(self, line):
        self.data.title_abbreviation = self.data.title_abbreviation + \
                                       self._clean(line, rstrip=0)
    
    def title(self, line):
        self.data.title = self.data.title + self._clean(line, rstrip=0)
    
    def transliterated_title(self, line):
        self.data.transliterated_title = self.data.transliterated_title + \
                                         self._clean(line, rstrip=0)
    
    def unique_identifier(self, line):
        assert not self.data.id, "id already defined"
        self.data.id = self._clean(line)
    
    def volume_issue(self, line):
        assert not self.data.volume_issue, "volume issue already defined"
        self.data.volume_issue = self._clean(line)
    
    def year(self, line):
        assert not self.data.year, "year already defined"
        self.data.year = self._clean(line)
    
    def pubmed_id(self, line):
        assert not self.data.pubmed_id, "PMID already defined"
        self.data.pubmed_id = self._clean(line)

    def undefined(self, line):
        # Records sometimes contain lines with qualifiers that don't match
        # any in the standard.  All these lines go into another variable.
        # Some undefined qualifiers:
        # 4098, 4099, 4100, 4101
        # 634
        # NP, PID, EDAT, MHDA
        
        self.data.undefined.append(line)

    def _clean(self, line, rstrip=1):
        tab = line.find('\t')
        if tab >= 0:
            nospace = line[tab+1:]
        elif line[:13] == ' purification':
            nospace = line[1:]
        else:
            nospace = line[6:]
        if rstrip:
            return nospace.rstrip()
        return nospace

    _needs_stripping = [
        'abstract', 'source', 'address', 'title_abbreviation',
        'title', 'transliterated_title'
        ]
    def _clean_record(self, rec):
        # Remove trailing newlines
        for m in self._needs_stripping:
            value = getattr(rec, m)
            setattr(rec, m, value.rstrip())
