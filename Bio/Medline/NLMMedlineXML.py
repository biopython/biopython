# Copyright 2001 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
This module provides code to work the NCBI's XML format for Medline.

Functions:
choose_format   Pick the right data format to use to index an XML file.
index           Index a Medline XML file.
index_many      Index multiple Medline XML files.

"""
# To Do:
# - Implement CitationParser
import os
import types
from xml.sax import handler

from Bio.ParserSupport import *
from Bio import MultiProc

import Martel

def choose_format(data):
    """choose_format(data) -> module

    Look at some data and choose the right format to parse it.  data
    should be the first 1000 characters or so of the file.  The module
    will contain 2 attributes: citation_format and format.
    citation_format is a Martel format to parse one citation.  format
    will parse the whole file.
    
    """
    formats = [
        ("nlmmedline_001211", "nlmmedline_001211_format"),
        ("nlmmedline_010319", "nlmmedline_010319_format"),
        ("nlmmedline_011101", "nlmmedline_011101_format"),
        ]
    for identifier, format_module in formats:
        if data.find(identifier) >= 0:
            break
    else:
        raise AssertionError, "I could not identify that format."
    package = '.'.join(["Bio", "Medline", format_module])
    return __import__(package, {}, {}, ["*"])

class Citation:
    """Holds information about a Medline citation.

    Members:
    medline_id         Medline ID.
    pmid               Pubmed ID.

    date_created       Tuple of (year, month, day, season, medline date).
    date_completed     Tuple of (year, month, day, season, medline date).
    date_revised       Tuple of (year, month, day, season, medline date).

    abstract           Tuple of (text, copyright info).
    journal            Tuple of (ISSN, volume, issue, date).
    article_title      Title of article.
    pagination         Tuple of (start, end, medline pagination).
    accession_numbers  List of accession numbers.
    affiliation        Affiliation.
    author_list        List of authors.
    languages          List of languages
    databank_list      List of tuples (name, accession numbers).
    grant_list         List of tuples (grant id, acronym, agency)
    publication_type_list   List of publication types.
    vernacular_title   Vernacular title.

    
    medline_journal_info   Tuple of (country, medline ta, medline code, nlm id)
    chemical_list          List of (CAS registry number, name).
    citation_subsets       List of citation subsets.
    comments_corrections   XXX not implemented
    gene_symbol_list       List of gene symbols.
    mesh_heading_list      List of (descriptor, subheadings).
    number_of_references   Number of references (int).
    personal_name_subject_list  List of personal names.
    
    """
    pass

class CitationParser(AbstractParser):
    """Parses a citation into a Record object.

    """
    def __init__(self):
        raise NotImplementedError

class _IndexerHandler(handler.ContentHandler):
    """Handles the results from the nlmmedline_format.  Saves the begin
    and end of each record as an offset from the beginning of the parse.

    """
    def __init__(self, found_citation_fn):
        """_IndexerHandler(found_citation_fn)

        found_citation_fn is called with the PMID, MedlineID, start,
        end where start and end are offsets from the beginning of the
        parse, with slice semantics.

        """
        self._citation_fn = found_citation_fn
        self._elements = []         # Open element tags.
        self._offset = 0            # Current file offset.
        self._start = None          # Offset of the start of the record.
        self._pmid = ''
        self._medline_id = ''
    def startElement(self, name, attrs):
        self._elements.append(name)
        if name == 'MedlineCitation':
            if self._start is not None:
                raise SyntaxError, "Found MedlineCitation, but already in one."
            self._start = self._offset
    def endElement(self, name):
        if not self._elements or self._elements[-1] != name:
            raise SyntaxError, "Elements not nested: %s" % name
        self._elements.pop()
        if name == 'MedlineCitation':
            if not self._pmid or not self._medline_id:  # didn't find an ID:
                raise SyntaxError, "I couldn't find an id: %s %s" % (
                    self._pmid, self._medline_id)
            self._citation_fn(
                self._pmid, self._medline_id, self._start, self._offset)
            self._start = None
            self._pmid = self._medline_id = ''
    def characters(self, content):
        self._offset += len(content)
        # Examine the tags directly under <MedlineCitation>.
        if len(self._elements)>=2 and self._elements[-2] == "MedlineCitation":
            if self._elements[-1] == "PMID":
                self._pmid = content
            elif self._elements[-1] == "MedlineID":
                self._medline_id = content

class _SavedDataHandle:
    def __init__(self, handle, saved):
        self.saved = saved
        self.handle = handle
    def read(self, length=None):
        if length is None:
            data = self.saved + self.handle.read()
            self.saved = ''
        else:
            data = self.saved[:length]
            data += self.handle.read(length-len(data))
            self.saved = self.saved[length:]
        return data

def index(handle, index_fn=None):
    """index(handle[, index_fn]) -> list of (PMID, MedlineID, start, end)

    Index a Medline XML file.  Returns where the records are, as
    offsets from the beginning of the handle.  index_fn is a callback
    function with parameters (PMID, MedlineID, start, end) and is
    called as soon as each record is indexes.

    """
    # Find the correct format to parse the data.
    data = handle.read(1000)
    format_module = choose_format(data)
    handle = _SavedDataHandle(handle, data)
    format = format_module.format
    wanted = ["MedlineCitation", "PMID", "MedlineID"]
    format = Martel.select_names(format, wanted)

    # Create an indexer that will save all the index information and
    # call index_fn if appropriate.
    indexes = []
    def citation_fn(pmid, medline_id, start, end,
                    indexes=indexes, index_fn=index_fn):
        if index_fn is not None:
            index_fn(pmid, medline_id, start, end)
        indexes.append((pmid, medline_id, start, end))
    indexer = _IndexerHandler(citation_fn)

    # Create the parser and parse the results.
    parser = format.make_parser(debug_level=0)
    parser.setContentHandler(indexer)
    parser.setErrorHandler(handler.ErrorHandler())
    parser.parseFile(handle)
    return indexes

def index_many(files_or_paths, index_fn, nprocs=1):
    """index_many(files_or_paths, index_fn[, nprocs])

    Index multiple Medline XML files.  files_or_paths can be a single
    file, a path, a list of files, or a list of paths.

    index_fn is a callback function that should take the following
    parameters:
    index_fn(file, event, data)
    
    where file is the file being indexed, event is one of "START",
    "RECORD", "END", and data is extra data dependent upon the event.
    "START" and "END" events are passed to indicate when a file is
    being indexed.  "RECORD" is passed whenever a new record has been
    indexed.  When a "RECORD" event is passed, then data is set to a
    tuple of (pmid, medline_id, start, end).  Otherwise it is None.
    start and end indicate the location of the record as offsets from
    the beginning of the file.

    """
    # This isn't a very good solution because it only allows 2 types
    # of sequences.  It's possible to use operator.isSequenceType, but
    # then we have to figure out how to exclude String types.
    if type(files_or_paths) not in [types.ListType, types.TupleType]:
        files_or_paths = [files_or_paths]
        
    files = []
    for f in files_or_paths:
        if os.path.isfile(f):
            files.append(f)
        elif os.path.isdir(f):
            names = os.listdir(f)
            for name in names:
                files.append(os.path.join(f, name))
        else:
            raise ValueError, "I can't find %s" % f

    def do_some(start, skip, files, index_fn):
        for i in range(start, len(files), skip):
            infile = files[i]
            index_fn(infile, "START", None)
            # index takes an optional index_fn with a different
            # interface than the callback for this function.  Thus, I
            # have to make an adapter to change the interface to one
            # that my client expects.
            def index_fn_adapter(pmid, medline_id, start, end,
                                 infile=infile, index_fn=index_fn):
                index_fn(infile, "RECORD", (pmid, medline_id, start, end))
            index(open(infile), index_fn_adapter)
            index_fn(infile, "END", None)
    MultiProc.run(nprocs, do_some, fn_args=(files, index_fn))
