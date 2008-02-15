# Copyright 1999-2000 by Jeffrey Chang.  All rights reserved.
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

_open

"""
import urllib

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

def pmfetch(db, id, report=None, mode=None,
            cgi="http://www.ncbi.nlm.nih.gov/entrez/utils/pmfetch.fcgi"):
    """pmfetch(db, id, report=None, mode=None,
    cgi="http://www.ncbi.nlm.nih.gov/entrez/utils/pmfetch.fcgi")

    Query PmFetch and return a handle to the results.  See the
    online documentation for an explanation of the parameters:
    http://www.ncbi.nlm.nih.gov/entrez/utils/pmfetch_help.html
    
    Raises an IOError exception if there's a network error.
    
    """
    # NCBI has retired PmFetch!!!
    import warnings
    warnings.warn("pmfetch is deprecated, as NCBI has retired PmFetch. Please let the Biopython developers know (biopython-dev@biopython.org) if you still use this function", DeprecationWarning)

    variables = {'db' : db, 'id' : id}
    if report is not None:
        variables['report'] = report
    if mode is not None:
        variables['mode'] = mode
    return _open(cgi, variables)

def pmqty(db, term, dopt=None, 
          cgi='http://www.ncbi.nlm.nih.gov/entrez/utils/pmqty.fcgi',
          **keywds):
    """pmqty(db, term, dopt=None,
    cgi='http://www.ncbi.nlm.nih.gov/entrez/utils/pmqty.fcgi') -> handle

    Query PmQty and return a handle to the results.  See the
    online documentation for an explanation of the parameters:
    http://www.ncbi.nlm.nih.gov/entrez/utils/pmqty_help.html
    
    Raises an IOError exception if there's a network error.
    
    """
    # NCBI has retired PmQty!!!
    import warnings
    warnings.warn("pmqty is deprecated, as NCBI has retired PmQty. Please let the Biopython developers know (biopython-dev@biopython.org) if you still use this function", DeprecationWarning)
    variables = {'db' : db, 'term' : term}
    if dopt is not None:
        variables['dopt'] = dopt
    variables.update(keywds)
    return _open(cgi, variables)

def pmneighbor(pmid, display,
               cgi='http://www.ncbi.nlm.nih.gov/entrez/utils/pmneighbor.fcgi'):
    """pmneighbor(pmid, display,
    cgi='http://www.ncbi.nlm.nih.gov/entrez/utils/pmneighbor.fcgi') -> handle

    Query PMNeighbor and return a handle to the results.  See the
    online documentation for an explanation of the parameters:
    http://www.ncbi.nlm.nih.gov/entrez/utils/pmneighbor_help.html
    
    Raises an IOError exception if there's a network error.
    
    """
    # NCBI has retired PmNeighbor!!!
    import warnings
    warnings.warn("pmneighbor is deprecated, as NCBI has retired PmNeighbor. Please let the Biopython developers know (biopython-dev@biopython.org) if you still use this function", DeprecationWarning)
    #
    # Warning: HUGE HACK HERE!  pmneighbor expects the display
    # parameter to be passed as just a tag, with no value.
    # Unfortunately, _open doesn't support these types of parameters,
    # so I'm building my own cgi string.  This is really due to the
    # limitations of urllib.urlencode.  We'll have to figure out a
    # good workaround.
    fullcgi = "%s?pmid=%s&%s" % (cgi, pmid, display)
    return _open(fullcgi)

# XXX retmode?
def epost(db, id, cgi='http://www.ncbi.nlm.nih.gov/entrez/eutils/epost.fcgi',
          **keywds):
    """epost(db, id[, cgi]) -> handle

    Query Entrez and return a handle to the results.

    Posts a file containing a list of UIs for future use in the user's
    environment to use with subsequent search strategies. See the online
    documentation for an explanation of the parameters:
    http://www.ncbi.nlm.nih.gov/entrez/query/static/epost_help.html

    Raises an IOError exception if there's a network error.

    """
    variables = {'db' : db, 'id' : id}
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

def _open(cgi, params={}, get=1):
    """_open(cgi, params={}, get=1) -> UndoHandle

    Open a handle to Entrez.  cgi is the URL for the cgi script to access.
    params is a dictionary with the options to pass to it.  get is a boolean
    that describes whether a GET should be used.  Does some
    simple error checking, and will raise an IOError if it encounters one.

    """
    # Open a handle to Entrez.
    options = urllib.urlencode(params, doseq=True)
    if get:  # do a GET
        fullcgi = cgi
        if options:
            fullcgi = "%s?%s" % (cgi, options)
        # print fullcgi
        handle = urllib.urlopen(fullcgi)
    else:    # do a POST
        handle = urllib.urlopen(cgi, options)

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
