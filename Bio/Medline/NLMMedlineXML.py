# Copyright 2001 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""NLMMedlineXML.py

This module provides code to work the NCBI's XML format for Medline.

Functions:
find_citations  Look for citations in a handle to a MedlineXML file.
index_path      Index a whole directory of MedlineXML files.

"""
# XXX To Do:
# - Implement CitationParser
# - index the indexes
import os
import time
from xml.sax import handler

from Bio.ParserSupport import *
from Bio.Tools.MultiProc import Scheduler, Task

import nlmmedline_xml_format

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
        self._pos = 0
        self._start_pos = None
        self._in_pmid = 0
        self._in_medline_id = 0
        self._pmid = ''
        self._medline_id = ''
        self._citation_fn = found_citation_fn
    def startElement(self, name, attrs):
        if name == 'MedlineCitation':
            if self._start_pos is not None:
                raise SyntaxError, "Found start MedlineCitation without end"
            self._start_pos = self._pos
        elif name == 'PMID':
            if self._in_pmid:
                raise SyntaxError, "Found start PMID without end"
            self._in_pmid = 1
        elif name == 'MedlineID':
            if self._in_medline_id:
                raise SyntaxError, "Found start MedlineID without end"
            self._in_medline_id = 1
    def endElement(self, name):
        if name == 'MedlineCitation':
            self._citation_fn(
                self._pmid, self._medline_id, self._start_pos, self._pos)
            self._start_pos = None
            self._pmid = ''
            self._medline_id = ''
        elif name == 'PMID':
            if not self._in_pmid:
                raise SyntaxError, "Found end PMID without start"
            self._in_pmid = 0
        elif name == 'MedlineID':
            if not self._in_medline_id:
                raise SyntaxError, "Found end MedlineID without start"
            self._in_medline_id = 0
    def characters(self, content):
        self._pos += len(content)
        if self._in_pmid:
            self._pmid += content
        elif self._in_medline_id:
            self._medline_id += content

def find_citations(handle):
    """find_citations(handle) -> list of (PMID, MedlineID, start, end)"""
    indexes = []
    def citation_fn(pmid, medline_id, start, end, indexes=indexes):
        indexes.append((pmid, medline_id, start, end))
        
    indexer = _IndexerHandler(citation_fn)
    wanted = ["MedlineCitation", "PMID", "MedlineID"]
    format = nlmmedline_xml_format.format
    parser = format.make_parser(debug_level=0)
    parser.setContentHandler(indexer)
    parser.setErrorHandler(handler.ErrorHandler())
    parser.parseFile(handle)
    return indexes

def index_path(inpath, outpath, noclobber=1, nprocs=1, update_fn=None):
    """index_path(inpath, outpath[, noclobber][, nprocs][, update_fn])

    Index a whole path of Medline XML files into another path of
    indexes.  update_fn is a callback that takes as a parameters the
    file being indexed.

    """
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    # Get a list of all the .xml files in the inpath.
    files = os.listdir(inpath)
    files = filter(lambda x: x.endswith(".xml"), files)

    def do_some(files, start, skip, inpath, outpath, noclobber, update_fn):
        for i in range(start, len(files), skip):
            file = files[i]
            # Create filenames for the infile and outfile.
            head, ext = os.path.splitext(file)
            infile = os.path.join(inpath, file)
            outfile = os.path.join(outpath, "%s.index" % head)

            if noclobber and os.path.exists(outfile):
                continue
            if update_fn is not None:
                update_fn(file)
                
            indexes = find_citations(open(infile))
            lines = []
            for pmid, medline_id, start, end in indexes:
                lines.append("PMID %s MedlineID %s START %d END %d\n" %
                             (pmid, medline_id, start, end))
            open(outfile, 'w').writelines(lines)

    if nprocs == 1:
        do_some(files, 0, 1, inpath, outpath, noclobber, update_fn)
    else:
        scheduler = Scheduler.Scheduler(nprocs)
        for i in range(nprocs):
            x = Task.Task(target=do_some,
                          args=(files, i, inpath, outpath,
                                noclobber, update_fn))
            scheduler.add(x)
        while scheduler.run():
            time.sleep(1.35)
