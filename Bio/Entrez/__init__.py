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

# The following four classes are used to add a member .attributes to integers,
# strings, lists, and dictionaries, respectively.

class AttributedInteger(int): pass

class AttributedString(str): pass

class AttributedList(list): pass

class AttributedDictionary(dict): pass

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

    DTDs = os.path.join(__path__[0], "DTDs")

    def __init__(self):
        self.path = []
	self.error = None
	self.booleans = []
	self.integers = []
	self.strings = []
	self.lists = []
        self.dictionaries = []
        self.structures = {}
        self.items = []
        self.initialized = False

    def startElement(self, name, attrs):
        if not self.initialized:
            # This XML file does not have a DTD; load its definitions here
            # using the first element in the XML.
            self.load_dtd_definitions(name)
        self.content = ""
        if name in self.lists:
            if attrs:
                object = AttributedList()
                object.attributes = dict(attrs)
            else:
                object = []
        elif name in self.dictionaries:
            if attrs:
                object = AttributedDictionary()
                object.attributes = dict(attrs)
            else:
                object = {}
        elif name in self.structures:
            object = Structure(self.structures[name])
            if attrs:
                object.attributes = dict(attrs)
        elif name in self.items:	# Only appears in ESummary
            name = str(attrs["Name"]) # convert from Unicode
            itemtype = str(attrs["Type"]) # convert from Unicode
            if itemtype=="Structure":
                object = {}
            elif name in ("ArticleIds", "History"):
                object = Structure(["pubmed", "medline"])
            elif itemtype=="List":
                object = []
            else:
                object = ""
                self.itemname = name
                self.itemtype = itemtype
        else:
            object = ""
            self.attributes = attrs
        if object!="" and len(self.path)!=0:
            current = self.path[-1]
            try:
                current.append(object)
            except AttributeError:
                current[name] = object
        self.path.append(object)

    def endElement(self, name):
        self.object = self.path.pop()
        # Convert Unicode strings to plain strings
        try:
            self.content = str(self.content)
        except UnicodeEncodeError:
            pass
        value = self.content
        if name==self.error and value!="":
            raise RuntimeError(value)
        elif name in self.booleans:
            if value=='Y':
                value = True
            elif value=='N':
                value = False
        elif name in self.integers:
            if self.attributes:
                value = AttributedInteger(self.content)
                value.attributes = dict(self.attributes)
                del self.attributes
            else:
                value = int(self.content)
        elif name in self.strings:
            if self.attributes:
                value = AttributedString(self.content)
                value.attributes = dict(self.attributes)
                del self.attributes
            else:
                value = self.content
        elif name in self.items:
            if self.object!="": return
            name = self.itemname
            value = self.content
            if self.itemtype=="Integer": value = int(value)
        else:
            return
        current = self.path[-1]
        try:
            current.append(value)
        except AttributeError:
            current[name] = value

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
        elif filename=="pubmed_080101.dtd":
            import PubmedArticleSet as module
        elif filename=="NCBI_BioSource.mod.dtd":
            import NCBI_BioSource as module
        elif filename=="NCBI_Entity.mod.dtd":
            # Nothing of interest in this DTD
            return
        elif filename=="NCBI_Entrezgene.dtd":
            # Nothing of interest in this DTD
            return
        elif filename=="NCBI_Entrezgene.mod.dtd":
            import NCBI_Entrezgene as module
        elif filename=="NCBI_Gene.mod.dtd":
            import NCBI_Gene as module
        elif filename=="NCBI_General.mod.dtd":
            import NCBI_General as module
        elif filename=="NCBI_Seqloc.mod.dtd":
            import NCBI_Seqloc as module
        elif filename=="NCBI_Mim.dtd":
            # Nothing of interest in this DTD
            return
        elif filename=="NCBI_Mim.mod.dtd":
            import NCBI_Mim as module
        elif filename=="NCBI_Organism.mod.dtd":
            import NCBI_Organism as module
        elif filename=="NCBI_Protein.mod.dtd":
            import NCBI_Protein as module
        elif filename=="nlmmedline_080101.dtd":
            import NLMMedline as module
        elif filename=="nlmmedlinecitation_080101.dtd":
            import NLMMedlineCitation as module
        elif filename=="nlmsharedcatcit_080101.dtd":
            import NLMSharedCatCit as module
        elif filename=="nlmcommon_080101.dtd":
            import NLMCommon as module
        elif filename=="taxon.dtd":
            import Taxon as module
        elif filename=="SerialSet":
            import  SerialSet as module
        else:
            import warnings
            warnings.warn("No parser available for %s; skipping its elements" % filename)
            return
        self.error = module.error
        self.booleans.extend(module.booleans)
        self.integers.extend(module.integers)
        self.strings.extend(module.strings)
        self.lists.extend(module.lists)
        self.dictionaries.extend(module.dictionaries)
        self.structures.update(module.structures)
        self.items.extend(module.items)
        self.initialized = True

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
    record = handler.object
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
