# Copyright 1999-2000 by Jeffrey Chang.  All rights reserved.
# Copyright 2008 by Michiel de Hoon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Provides code to access NCBI over the WWW.

The main Entrez web page is available at:
http://www.ncbi.nlm.nih.gov/Entrez/

A list of the Entrez utilities is available at:
http://www.ncbi.nlm.nih.gov/entrez/utils/utils_index.html


Functions:
query        Query Entrez; retrieve results in HTML format.

efetch       Retrieves records in the requested format from a list of one or
             more primary IDs or from the user's environment
epost        Posts a file containing a list of primary IDs for future use in
             the user's environment to use with subsequent search strategies
esearch      Searches and retrieves primary IDs (for use in EFetch, ELink,
             and ESummary) and term translations and optionally retains
             results for future use in the user's environment.
elink        Checks for the existence of an external or Related Articles link
             from a list of one or more primary IDs.  Retrieves primary IDs
             and relevancy scores for links to Entrez databases or Related
             Articles;  creates a hyperlink to the primary LinkOut provider
             for a specific ID and database, or lists LinkOut URLs
             and Attributes for multiple IDs.
einfo        Provides field index term counts, last update, and available
             links for each database.
esummary     Retrieves document summaries from a list of primary IDs or from
             the user's environment.
egquery      Provides Entrez database counts in XML for a single search
             using Global Query.
espell       Retrieves spelling suggestions.

_open        Internally used function.

"""
import urllib, time
from xml.sax.handler import ContentHandler, EntityResolver
from xml.sax import make_parser
import os.path
from Bio import File


def query(cmd, db, cgi='http://www.ncbi.nlm.nih.gov/sites/entrez',
          **keywds):
    """query(cmd, db, cgi='http://www.ncbi.nlm.nih.gov/sites/entrez',
    **keywds) -> handle

    Query Entrez and return a handle to the results, consisting of
    a web page in HTML format.
    See the online documentation for an explanation of the parameters:
    http://www.ncbi.nlm.nih.gov/books/bv.fcgi?rid=helplinks.chapter.linkshelp

    Raises an IOError exception if there's a network error.

    """
    variables = {'cmd' : cmd, 'db' : db}
    variables.update(keywds)
    return _open(cgi, variables)

# XXX retmode?
def epost(db, cgi='http://www.ncbi.nlm.nih.gov/entrez/eutils/epost.fcgi',
          **keywds):
    """epost(db, id|(WebEnv,query_key), [, cgi]) -> handle

    Query Entrez and return a handle to the results.

    Posts a file containing a list of UIs for future use in the user's
    environment to use with subsequent search strategies. See the online
    documentation for an explanation of the parameters:
    http://www.ncbi.nlm.nih.gov/entrez/query/static/epost_help.html

    Raises an IOError exception if there's a network error.

    """
    variables = {'db' : db}
    variables.update(keywds)
    return _open(cgi, variables)

def efetch(db, cgi='http://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi',
          **keywds):
    """efetch(db[, cgi][...]) -> handle

    Query Entrez and return a handle to the results.

    EFetch retrieves records in the requested format from a list of one or
    more UIs or from user's environment. See the online
    documentation for an explanation of the parameters:
    http://www.ncbi.nlm.nih.gov/entrez/query/static/efetch_help.html

    Raises an IOError exception if there's a network error.

    """
    variables = {'db' : db}
    variables.update(keywds)
    return _open(cgi, variables)

def esearch(db, term,
            cgi='http://www.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi',
            **keywds):
    """esearch(db, term[, cgi][...]) -> handle

    Query Entrez and return a handle to the results.

    ESearch searches and retrieves primary IDs (for use in EFetch, ELink
    and ESummary) and term translations, and optionally retains results
    for future use in the user's environment. See the online
    documentation for an explanation of the parameters:
    http://www.ncbi.nlm.nih.gov/entrez/query/static/esearch_help.html

    Raises an IOError exception if there's a network error.

    """
    variables = {'db' : db,
                 'term' : term}
    variables.update(keywds)
    return _open(cgi, variables)

def elink(cgi='http://www.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi',
          **keywds):
    """elink([, cgi][...]) -> handle

    Query Entrez and return a handle to the results.

    ELink checks for the existence of an external or Related Articles link
    from a list of one or more primary IDs;  retrieves IDs and relevancy
    scores for links to Entrez databases or Related Articles; creates a
    hyperlink to the primary LinkOut provider for a specific ID and
    database, or lists LinkOut URLs and attributes for multiple IDs. See
    the online documentation for an explanation of the parameters:
    http://www.ncbi.nlm.nih.gov/entrez/query/static/elink_help.html

    Raises an IOError exception if there's a network error.

    """
    variables = {}
    variables.update(keywds)
    return _open(cgi, variables)

def einfo(cgi='http://www.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi',
          **keywds):
    """einfo([, cgi][...]) -> handle

    Query Entrez and return a handle to the results.

    EInfo provides field names, index term counts, last update, and
    available links for each Entrez database. See the online
    documentation for an explanation of the parameters:
    http://www.ncbi.nlm.nih.gov/entrez/query/static/einfo_help.html

    Raises an IOError exception if there's a network error.

    """
    variables = {}
    variables.update(keywds)
    return _open(cgi, variables)

def esummary(cgi='http://www.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi',
          **keywds):
    """esummary([, cgi][...]) -> handle

    Query Entrez and return a handle to the results.

    ESummary retrieves document summaries from a list of primary IDs or
    from the user's environment. See the online documentation for an
    explanation of the parameters:
    http://www.ncbi.nlm.nih.gov/entrez/query/static/esummary_help.html

    Raises an IOError exception if there's a network error.

    """
    variables = {}
    variables.update(keywds)
    return _open(cgi, variables)

def egquery(cgi='http://www.ncbi.nlm.nih.gov/entrez/eutils/egquery.fcgi',
          **keywds):
    """egquery([, cgi][...]) -> handle

    Query Entrez and return a handle to the results.

    EGQuery provides Entrez database counts in XML for a single search
    using Global Query. See the online documentation for an explanation
    of the parameters:
    http://www.ncbi.nlm.nih.gov/entrez/query/static/egquery_help.html

    Raises an IOError exception if there's a network error.

    """
    variables = {}
    variables.update(keywds)
    return _open(cgi, variables)

def espell(cgi='http://www.ncbi.nlm.nih.gov/entrez/eutils/espell.fcgi',
          **keywds):
    """espell([, cgi][...]) -> handle

    Query Entrez and return a handle to the results.

    ESpell retrieves spelling suggestions, if available. See the online
    documentation for an explanation of the parameters:
    http://www.ncbi.nlm.nih.gov/entrez/query/static/espell_help.html

    Raises an IOError exception if there's a network error.

    """
    variables = {}
    variables.update(keywds)
    return _open(cgi, variables)

class AttributedInteger(int):
    def __new__(cls, value, attributes):
        self = int.__new__(cls, value)
        self.attributes = {}
        keys = attributes.keys()
        for key in keys:
            self.attributes[key] = attributes[key]
        return self

class AttributedString(str):
    def __new__(cls, value, attributes):
        self = str.__new__(cls, value)
        self.attributes = {}
        keys = attributes.keys()
        for key in keys:
            self.attributes[key] = attributes[key]
        return self

class Structure(dict):
    def __init__(self, keys):
        dict.__init__(self)
        for key in keys:
            dict.__setitem__(self, key, [])
        self.listkeys = keys
    def __setitem__(self, key, value):
        if key in self.listkeys:
            self[key].append(value)
        else:
            dict.__setitem__(self, key, value)

class DataHandler(ContentHandler, EntityResolver):
    from Bio.Entrez import EInfo, ESearch, ESummary, EPost, ELink, EGQuery, ESpell, Taxon, PubmedArticleSet, SerialSet, NCBI_Mim, Entrezgene_Set
    _NameToModule = {"eInfoResult": EInfo,
                     "eSearchResult": ESearch,
                     "eSummaryResult": ESummary,
                     "ePostResult": EPost,
                     "eLinkResult": ELink,
                     "Result": EGQuery,
                     "eSpellResult": ESpell,
                     "TaxaSet": Taxon,
                     "PubmedArticleSet": PubmedArticleSet,
                     "SerialSet": SerialSet,
                     "Mim-entries": NCBI_Mim,
                     "Entrezgene-Set": Entrezgene_Set,
                    }

    DTDs = os.path.join(__path__[0], "DTDs")

    def __init__(self):
        self.element = []
        self.path = []
	self.error = None
	self.booleans = []
	self.integers = []
	self.strings = []
	self.lists = []
        self.dictionaries = []
        self.structures = {}
        self.items = []
        self.handleStartElement = None
        self.handleEndElement = None

    def startElement(self, name, attrs):
        self.attributes = attrs
        self.element.append(name)
        self.content = ""
        if name in self.lists:
            object = []
            try:
                current = self.path[-1]
            except IndexError:
                current = None
            if isinstance(current, dict):
                current[name] = object
            elif isinstance(current, list):
                current.append(object)
            self.path.append(object)
        elif name in self.dictionaries:
            object = {}
            try:
                current = self.path[-1]
            except IndexError:
                current = None
            if isinstance(current, dict):
                current[name] = object
            elif isinstance(current, list):
                current.append(object)
            self.path.append(object)
        elif name in self.structures:
            current = self.path[-1]
            object = Structure(self.structures[name])
            if isinstance(current, dict):
                current[name] = object
            elif isinstance(current, list):
                current.append(object)
            self.path.append(object)
        elif name in self.items:	# Only appears in ESummary
            current = self.path[-1]
            itemname = str(attrs["Name"]) # convert from Unicode
            itemtype = str(attrs["Type"]) # convert from Unicode
            if itemtype=="Structure":
                object = {}
            elif itemname in ("ArticleIds", "History"):
                object = Structure(["pubmed", "medline"])
            elif itemtype=="List":
                object = []
            else:
                object = ""
            self.itemname = itemname
            self.itemtype = itemtype
            if object!="":
                if isinstance(current, dict):
                    current[itemname] = object
                elif isinstance(current, list):
                    current.append(object)
            self.path.append(object)
        elif self.strings: # Just for checking if this is the new approach
            current = self.path[-1]
            object = ""
            self.path.append(object)
        elif self.handleStartElement:
            self.handleStartElement(self, name, attrs)
        elif name in DataHandler._NameToModule:
            self.handleStartElement = DataHandler._NameToModule[name].startElement
            self.handleEndElement = DataHandler._NameToModule[name].endElement
            self.handleStartElement(self, name, attrs)
        else:
            raise RuntimeError("No parser available for " + name)

    def endElement(self, name):
        try:
            self.object = self.path[-1]
        except IndexError:
            self.object = None # This will only occur for the Taxonomy parser
        # Convert Unicode strings to plain strings
        try:
            self.content = str(self.content)
        except UnicodeEncodeError:
            pass
        value = self.content
        if name==self.error and value!="":
            raise RuntimeError(value)
        elif name in self.booleans:
            self.path = self.path[:-1]
            if self.content=='Y':
                value = True
            elif self.content=='N':
                value = False
            current = self.path[-1]
            if isinstance(current, list):
                current.append(value)
            elif isinstance(current, dict):
                current[name] = value
        elif name in self.integers:
            self.path = self.path[:-1]
            if self.attributes:
                value = AttributedInteger(self.content, self.attributes)
            else:
                value = int(self.content)
            current = self.path[-1]
            if isinstance(current, list):
                current.append(value)
            elif isinstance(current, dict):
                current[name] = value
        elif name in self.strings:
            self.path = self.path[:-1]
            if self.attributes:
                value = AttributedString(self.content, self.attributes)
            else:
                value = self.content
            current = self.path[-1]
            if isinstance(current, list):
                current.append(value)
            elif isinstance(current, dict):
                current[name] = value
        elif name in self.items:
            object = self.path.pop()
            if object=="":
                current = self.path[-1]
                if isinstance(current, list):
                    current.append(self.content)
                elif isinstance(current, dict):
                    itemname = self.itemname
                    itemtype = self.itemtype
                    value = self.content
                    if itemtype=="Integer": value = int(value)
                    current[itemname] = value
        elif name in DataHandler._NameToModule:
            self.handleStartElement = None
            self.handleEndElement = None
        elif self.handleEndElement:
            self.handleEndElement(self, name)
        self.element = self.element[:-1]

    def characters(self, content):
        self.content += content

    def resolveEntity(self, publicId, systemId):
        location, filename = os.path.split(systemId)
        self.load_dtd_definitions(filename)
        path = os.path.join(DataHandler.DTDs, filename)
        try:
            handle = open(path)
        except IOError:
            import warnings
            warnings.warn("DTD file %s not found in Biopython installation; trying to retrieve it from NCBI" % filename)
            handle = EntityResolver.resolveEntity(self, publicId, systemId)
        return handle

    def load_dtd_definitions(self, filename):
        if filename=="eInfo_020511.dtd":
            import EInfo as module
        elif filename=="eSearch_020511.dtd":
            import ESearch as module
        elif filename=="ePost_020511.dtd":
            import EPost as module
        elif filename=="eSummary_041029.dtd":
            import ESummary as module
        elif filename=="eLink_020511.dtd":
            import ELink as module
        elif filename=="eSpell.dtd":
            import ESpell as module
        elif filename=="egquery.dtd":
            import EGQuery as module
        elif filename=="NCBI_BioSource.mod.dtd":
            import NCBI_BioSource as module
        elif filename=="NCBI_Entrezgene.mod.dtd":
            import NCBI_Entrezgene as module
        elif filename=="NCBI_Seqloc.mod.dtd":
            import NCBI_Seqloc as module
        else:
            return
        self.error = module.error
        self.booleans.extend(module.booleans)
        self.integers.extend(module.integers)
        self.strings.extend(module.strings)
        self.lists.extend(module.lists)
        self.dictionaries.extend(module.dictionaries)
        self.structures.update(module.structures)
        self.items.extend(module.items)
        self.handleStartElement = module.startElement
        self.handleEndElement = module.endElement

def read(handle):
    """read(hande) -> record
    
    This function parses an XML from Entrez, typically returning the
    data as a list of dictionaries.  An appropriate parser will be
    automatically be selected if available.
    """
    saxparser = make_parser()
    handler = DataHandler()
    saxparser.setContentHandler(handler)
    saxparser.setEntityResolver(handler)
    saxparser.parse(handle)
    if handler.object:
        record = handler.object
    else:
        record = handler.record # Only for the Taxonompy parser
    return record

def _open(cgi, params={}):
    """_open(cgi, params={}) -> UndoHandle

    Open a handle to Entrez.  cgi is the URL for the cgi script to access.
    params is a dictionary with the options to pass to it.  Does some
    simple error checking, and will raise an IOError if it encounters one.

    """
    # NCBI requirement: At least three seconds between queries
    delay = 3.0
    current = time.time()
    wait = _open.previous + delay - current
    if wait > 0:
        time.sleep(wait)
        _open.previous = current + wait
    else:
        _open.previous = current
    # Tell Entrez that we are Biopython
    if not "tool" in params:
        params["tool"] = "biopython"
    # Open a handle to Entrez.
    options = urllib.urlencode(params, doseq=True)
    cgi += "?" + options
    handle = urllib.urlopen(cgi)

    # Wrap the handle inside an UndoHandle.
    uhandle = File.UndoHandle(handle)

    # Check for errors in the first 5 lines.
    # This is kind of ugly.
    lines = []
    for i in range(5):
        lines.append(uhandle.readline())
    for i in range(4, -1, -1):
        uhandle.saveline(lines[i])
    data = ''.join(lines)
                   
    if "500 Proxy Error" in data:
        # Sometimes Entrez returns a Proxy Error instead of results
        raise IOError, "500 Proxy Error (NCBI busy?)"
    elif "502 Proxy Error" in data:
        raise IOError, "502 Proxy Error (NCBI busy?)"
    elif "WWW Error 500 Diagnostic" in data:
        raise IOError, "WWW Error 500 Diagnostic (NCBI busy?)"
    elif data[:5] == "ERROR":
        # XXX Possible bug here, because I don't know whether this really
        # occurs on the first line.  I need to check this!
        raise IOError, "ERROR, possibly because id not available?"
    # Should I check for 404?  timeout?  etc?
    return uhandle

_open.previous = 0
