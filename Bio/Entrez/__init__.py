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

read         Parses the XML results returned by any of the above functions.
             Typical usage is:
             >>> handle = Entrez.einfo() # or esearch, efetch, ...
             >>> record = Entrez.read(handle)
             where record is now a Python dictionary or list.

_open        Internally used function.

"""
import urllib, time, warnings
import os.path
from Bio import File


email = None

def query(cmd, db, cgi='http://www.ncbi.nlm.nih.gov/sites/entrez',
          **keywds):
    """Query Entrez and return a handle to the HTML results (DEPRECATED).

    See the online documentation for an explanation of the parameters:
    http://www.ncbi.nlm.nih.gov/books/bv.fcgi?rid=helplinks.chapter.linkshelp

    Return a handle to the results.

    Raises an IOError exception if there's a network error.
    """
    import warnings
    warnings.warn("Bio.Entrez.query is deprecated, since it breaks NCBI's rule to only use the E-Utilities URL.", DeprecationWarning)

# XXX retmode?
def epost(db, **keywds):
    """Post a file of identifiers for future use.

    Posts a file containing a list of UIs for future use in the user's
    environment to use with subsequent search strategies.

    See the online documentation for an explanation of the parameters:
    http://www.ncbi.nlm.nih.gov/entrez/query/static/epost_help.html

    Return a handle to the results.

    Raises an IOError exception if there's a network error.
    """
    cgi='http://eutils.ncbi.nlm.nih.gov/entrez/eutils/epost.fcgi'
    variables = {'db' : db}
    variables.update(keywds)
    return _open(cgi, variables, post=True)

def efetch(db, **keywds):
    """Fetches Entrez results which are returned as a handle.

    EFetch retrieves records in the requested format from a list of one or
    more UIs or from user's environment.

    See the online documentation for an explanation of the parameters:
    http://www.ncbi.nlm.nih.gov/entrez/query/static/efetch_help.html

    Return a handle to the results.

    Raises an IOError exception if there's a network error.

    Short example:

    from Bio import Entrez
    handle = Entrez.efetch(db="nucleotide", id="57240072", rettype="gb")
    print handle.read()
    """
    for key in keywds :
        if key.lower()=="rettype" and keywds[key].lower()=="genbank" :
            import warnings
            warnings.warn('As of Easter 2009, Entrez EFtech no longer '
                          'supports the unofficial return type "genbank", '
                          'use "gb" or "gp" instead.', DeprecationWarning)
            if db.lower()=="protein" :
                keywds[key] = "gp" #GenPept
            else :
                keywds[key] = "gb" #GenBank
    cgi='http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    variables = {'db' : db}
    variables.update(keywds)
    return _open(cgi, variables)

def esearch(db, term, **keywds):
    """ESearch runs an Entrez search and returns a handle to the results.

    ESearch searches and retrieves primary IDs (for use in EFetch, ELink
    and ESummary) and term translations, and optionally retains results
    for future use in the user's environment.

    See the online documentation for an explanation of the parameters:
    http://www.ncbi.nlm.nih.gov/entrez/query/static/esearch_help.html

    Return a handle to the results which are always in XML format.

    Raises an IOError exception if there's a network error.

    Short example:

    from Bio import Entez
    handle = Entrez.esearch(db="nucleotide", retmax=10, term="Opuntia")
    record = Entrez.read(handle)
    print record["Count"]
    print record["IdList"]
    """
    cgi='http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
    variables = {'db' : db,
                 'term' : term}
    variables.update(keywds)
    return _open(cgi, variables)

def elink(**keywds):
    """ELink checks for linked external articles and returns a handle.

    ELink checks for the existence of an external or Related Articles link
    from a list of one or more primary IDs;  retrieves IDs and relevancy
    scores for links to Entrez databases or Related Articles; creates a
    hyperlink to the primary LinkOut provider for a specific ID and
    database, or lists LinkOut URLs and attributes for multiple IDs.

    See the online documentation for an explanation of the parameters:
    http://www.ncbi.nlm.nih.gov/entrez/query/static/elink_help.html

    Return a handle to the results, by default in XML format.

    Raises an IOError exception if there's a network error.
    """
    cgi='http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi'
    variables = {}
    variables.update(keywds)
    return _open(cgi, variables)

def einfo(**keywds):
    """EInfo returns a summary of the Entez databases as a results handle.

    EInfo provides field names, index term counts, last update, and
    available links for each Entrez database.

    See the online documentation for an explanation of the parameters:
    http://www.ncbi.nlm.nih.gov/entrez/query/static/einfo_help.html

    Return a handle to the results, by default in XML format.

    Raises an IOError exception if there's a network error.

    Short example:

    from Bio import Entrez 
    record = Entrez.read(Entrez.einfo())
    print record['DbList']
    """
    cgi='http://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi'
    variables = {}
    variables.update(keywds)
    return _open(cgi, variables)

def esummary(**keywds):
    """ESummary retrieves document summaries as a results handle.

    ESummary retrieves document summaries from a list of primary IDs or
    from the user's environment.

    See the online documentation for an explanation of the parameters:
    http://www.ncbi.nlm.nih.gov/entrez/query/static/esummary_help.html

    Return a handle to the results, by default in XML format.

    Raises an IOError exception if there's a network error.
    """
    cgi='http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi'
    variables = {}
    variables.update(keywds)
    return _open(cgi, variables)

def egquery(**keywds):
    """EGQuery provides Entrez database counts for a global search.

    EGQuery provides Entrez database counts in XML for a single search
    using Global Query.

    See the online documentation for an explanation of the parameters:
    http://www.ncbi.nlm.nih.gov/entrez/query/static/egquery_help.html

    Return a handle to the results in XML format.

    Raises an IOError exception if there's a network error.
    """
    cgi='http://eutils.ncbi.nlm.nih.gov/entrez/eutils/egquery.fcgi'
    variables = {}
    variables.update(keywds)
    return _open(cgi, variables)

def espell(**keywds):
    """ESpell retrieves spelling suggestions, returned in a results handle.

    ESpell retrieves spelling suggestions, if available.

    See the online documentation for an explanation of the parameters:
    http://www.ncbi.nlm.nih.gov/entrez/query/static/espell_help.html

    Return a handle to the results, by default in XML format.

    Raises an IOError exception if there's a network error.

    Short example:

    from Bio import Entrez 
    record = Entrez.read(Entrez.espell(term="biopythooon")) 
    print record["Query"] 
    print record["CorrectedQuery"] 
    """
    cgi='http://eutils.ncbi.nlm.nih.gov/entrez/eutils/espell.fcgi'
    variables = {}
    variables.update(keywds)
    return _open(cgi, variables)

def read(handle):
    """Parses an XML file from the NCBI Entrez Utilities into python objects.
    
    This function parses an XML file created by NCBI's Entrez Utilities,
    returning a multilevel data structure of Python lists and dictionaries.
    Most XML files returned by NCBI's Entrez Utilities can be parsed by
    this function, provided its DTD is available. Biopython includes the
    DTDs for most commonly used Entrez Utilities.

    Whereas the data structure seems to consist of generic Python lists,
    dictionaries, strings, and so on, each of these is actually a class
    derived from the base type. This allows us to store the attributes
    (if any) of each element in a dictionary my_element.attributes, and
    the tag name in my_element.tag.
    """
    from Parser import DataHandler
    DTDs = os.path.join(__path__[0], "DTDs")
    handler = DataHandler(DTDs)
    record = handler.run(handle)
    return record

def parse(handle):
    from Parser import DataHandler
    DTDs = os.path.join(__path__[0], "DTDs")
    handler = DataHandler(DTDs)
    records = handler.parse(handle)
    return records

def _open(cgi, params={}, post=False):
    """Helper function to build the URL and open a handle to it (PRIVATE).

    Open a handle to Entrez.  cgi is the URL for the cgi script to access.
    params is a dictionary with the options to pass to it.  Does some
    simple error checking, and will raise an IOError if it encounters one.

    This function also enforces the "up to three queries per second rule"
    to avoid abusing the NCBI servers.
    """
    # NCBI requirement: At most three queries per second.
    # Equivalently, at least a third of second between queries
    delay = 0.333333334
    current = time.time()
    wait = _open.previous + delay - current
    if wait > 0:
        time.sleep(wait)
        _open.previous = current + wait
    else:
        _open.previous = current
    # Remove None values from the parameters
    for key, value in params.items():
        if value is None:
            del params[key]
    # Tell Entrez that we are using Biopython
    if not "tool" in params:
        params["tool"] = "biopython"
    # Tell Entrez who we are
    if not "email" in params:
        if email!=None:
            params["email"] = email
    # Open a handle to Entrez.
    options = urllib.urlencode(params, doseq=True)
    if post :
        #HTTP POST
        handle = urllib.urlopen(cgi, data=options)
    else :
        #HTTP GET
        cgi += "?" + options
        handle = urllib.urlopen(cgi)

    # Wrap the handle inside an UndoHandle.
    uhandle = File.UndoHandle(handle)

    # Check for errors in the first 7 lines.
    # This is kind of ugly.
    lines = []
    for i in range(7):
        lines.append(uhandle.readline())
    for i in range(6, -1, -1):
        uhandle.saveline(lines[i])
    data = ''.join(lines)
                   
    if "500 Proxy Error" in data:
        # Sometimes Entrez returns a Proxy Error instead of results
        raise IOError("500 Proxy Error (NCBI busy?)")
    elif "502 Proxy Error" in data:
        raise IOError("502 Proxy Error (NCBI busy?)")
    elif "WWW Error 500 Diagnostic" in data:
        raise IOError("WWW Error 500 Diagnostic (NCBI busy?)")
    elif "<title>Service unavailable!</title>" in data :
        #Probably later in the file it will say "Error 503"
        raise IOError("Service unavailable!")
    elif "<title>Bad Gateway!</title>" in data :
        #Probably later in the file it will say:
        #  "The proxy server received an invalid
        #   response from an upstream server."
        raise IOError("Bad Gateway!")
    elif "<title>414 Request-URI Too Large</title>" in data \
    or "<h1>Request-URI Too Large</h1>" in data :
        raise IOError("Requested URL too long (try using EPost?)")
    elif data.startswith("Error:") :
        #e.g. 'Error: Your session has expired. Please repeat your search.\n'
        raise IOError(data.strip())
    elif data.startswith("The resource is temporarily unavailable") :
        #This can occur with an invalid query_key
        #Perhaps this should be a ValueError?
        raise IOError("The resource is temporarily unavailable")
    elif data.startswith("download dataset is empty") :
        #This can occur when omit the identifier, or the WebEnv and query_key
        #Perhaps this should be a ValueError?
        raise IOError("download dataset is empty")
    elif data[:5] == "ERROR":
        # XXX Possible bug here, because I don't know whether this really
        # occurs on the first line.  I need to check this!
        raise IOError("ERROR, possibly because id not available?")
    # Should I check for 404?  timeout?  etc?
    return uhandle

_open.previous = 0
