# Copyright 1999-2000 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Provides code to access NCBI over the WWW.

The main Entrez web page is available at:
http://www.ncbi.nlm.nih.gov/Entrez/

A list of the Entrez utilities (will go away Dec 2002) is available
at:
http://www.ncbi.nlm.nih.gov/entrez/utils/utils_index.html

Documentation for the e-utilies are available at:
http://www.ncbi.nlm.nih.gov/entrez/query/static/eutils_help.html

The main Blast web page is available at:
http://www.ncbi.nlm.nih.gov/BLAST/


Functions:
query        Query Entrez.
pmfetch      Retrieve results using a unique identifier.
pmqty        Search PubMed.
pmneighbor   Return a list of related articles for a PubMed entry.

efetch       Access the efetch script.
_open

"""
import string
import urllib

from Bio import File

def query(cmd, db, cgi='http://www.ncbi.nlm.nih.gov/entrez/query.fcgi',
          **keywds):
    """query(cmd, db, cgi='http://www.ncbi.nlm.nih.gov/entrez/query.fcgi',
    **keywds) -> handle

    Query Entrez and return a handle to the results.  See the online
    documentation for an explanation of the parameters:
    http://www.ncbi.nlm.nih.gov/entrez/query/static/linking.html

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

    Query Entrez and return a handle to the results.  See the online
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

    Query Entrez and return a handle to the results.  See the online
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

    Query Entrez and return a handle to the results.  See the online
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

    Query Entrez and return a handle to the results.  See the online
    documentation for an explanation of the parameters:
    http://www.ncbi.nlm.nih.gov/entrez/query/static/elink_help.html

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
    options = urllib.urlencode(params)
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
    data = string.join(lines, '')
                   
    if string.find(data, "500 Proxy Error") >= 0:
        # Sometimes Entrez returns a Proxy Error instead of results
        raise IOError, "500 Proxy Error (NCBI busy?)"
    elif string.find(data, "502 Proxy Error") >= 0:
        raise IOError, "502 Proxy Error (NCBI busy?)"
    elif string.find(data, "WWW Error 500 Diagnostic") >= 0:
        raise IOError, "WWW Error 500 Diagnostic (NCBI busy?)"
    elif data[:5] == "ERROR":
        # XXX Possible bug here, because I don't know whether this really
        # occurs on the first line.  I need to check this!
        raise IOError, "ERROR, possibly because id not available?"
    # Should I check for 404?  timeout?  etc?
    return uhandle
    
