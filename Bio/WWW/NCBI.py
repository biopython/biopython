# Copyright 1999-2000 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""NCBI.py

Provides code to access NCBI over the WWW.

The main Entrez web page is available at:
http://www.ncbi.nlm.nih.gov/Entrez/

A list of the Entrez utilities is available at:
http://www.ncbi.nlm.nih.gov/entrez/utils/utils_index.html

The main Blast web page is available at:
http://www.ncbi.nlm.nih.gov/BLAST/


Functions:
query        Query Entrez.
pmfetch      Retrieve results using a unique identifier.
pmqty        Search PubMed.
pmneighbor   Return a list of related articles for a PubMed entry.
_open

"""
import string
import urllib
import sgmllib
import urlparse
import time

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
    variables = {'pmid' : pmid, 'display' : display}
    return _open(cgi, variables)

def blast(program, datalib, sequence,
          input_type='Sequence in FASTA format',
          double_window=None, gi_list='(None)', expect='10',
          filter='L', genetic_code='Standard (1)',
          mat_param='PAM30     9       1',
          other_advanced=None, ncbi_gi=None, overview=None,
          alignment_view='0', descriptions=None, alignments=None,
          email=None, path=None, html=None, 
          cgi='http://www.ncbi.nlm.nih.gov/blast/blast.cgi'
          ):
    """

    """
    params = {'PROGRAM' : program,
              'DATALIB' : datalib,
              'SEQUENCE' : sequence,
              'DOUBLE_WINDOW' : double_window,
              'INPUT_TYPE' : input_type,
              'EXPECT' : expect,
              'FILTER' : filter,
              'GENETIC_CODE' : genetic_code,
              'MAT_PARAM' : mat_param,
              'OTHER_ADVANCED' : other_advanced,
              'NCBI_GI' : ncbi_gi,
              'OVERVIEW' : overview,
              'ALIGNMENT_VIEW' : alignment_view,
              'DESCRIPTIONS' : descriptions,
              'ALIGNMENTS' : alignments,
              'EMAIL' : email,
              'PATH' : path,
              'HTML' : html
              }
    variables = {}
    for k in params.keys():
        if params[k] is not None:
            variables[k] = str(params[k])
    # This returns a handle to the HTML file that points to the results.
    handle = _open(cgi, variables, get=0)
    # Now parse the HTML from the handle and figure out how to retrieve
    # the results.
    class RefPageParser(sgmllib.SGMLParser):
        def __init__(self, cgi):
            sgmllib.SGMLParser.__init__(self)
            self.cgi = cgi
            self.params = {}
        def do_form(self, attributes):
            for attr, value in attributes:
                attr = string.upper(attr)
                if attr == 'ACTION':
                    self.cgi = urlparse.urljoin(self.cgi, value)
        def do_input(self, attributes):
            is_rid = 0
            rid = None
            for attr, value in attributes:
                attr, value = string.upper(attr), string.upper(value)
                if attr == 'NAME' and value == 'RID':
                    is_rid = 1
                elif attr == 'VALUE':
                    rid = value
            if is_rid and rid:
                self.params['RID'] = rid
    parser = RefPageParser(cgi)
    parser.feed(handle.read())
    if not parser.params.has_key('RID'):
        raise SyntaxError, "Error getting BLAST results: RID not found"
    # This gets the CGI page and the params to retrieve the BLAST results.
    cgi, params = parser.cgi, parser.params
    class ResultsParser(sgmllib.SGMLParser):
        def __init__(self):
            sgmllib.SGMLParser.__init__(self)
            self.ready = 0
        def handle_comment(self, comment):
            comment = string.lower(comment)
            if string.find(comment, 'status=ready') >= 0:
                self.ready = 1
    start = time.time()
    while 1:
        # Sometimes the BLAST results aren't done yet.  Look at the page
        # to see if the results are there.  If not, then try again later.
        results = _open(cgi, params, get=0).read()
        parser = ResultsParser()
        parser.feed(results)
        if parser.ready:
            break
        # Time out if it's not done after 15 minutes.
        if time.time() - start > 15*60:
            raise IOError, "timed out"
        # pause for 5 seconds and try again.
        time.sleep(5)
    return File.UndoHandle(File.StringHandle(results))

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
        raise IOError, "500 Proxy Error (Entrez busy?)"
    elif string.find(data, "WWW Error 500 Diagnostic") >= 0:
        raise IOError, "WWW Error 500 Diagnostic (Entrez busy?)"
    elif data[:5] == "ERROR":
        # XXX Possible bug here, because I don't know whether this really
        # occurs on the first line.  I need to check this!
        raise IOError, "ERROR, possibly because id not available?"
    # Should I check for 404?  timeout?  etc?
    return uhandle
    
