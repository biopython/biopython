# Copyright 1999-2000 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""PubMed.py

This module provides code to work with PubMed from the NCBI.
http://www.ncbi.nlm.nih.gov/PubMed/

Online documentation for linking to PubMed is available at:
http://www.ncbi.nlm.nih.gov/PubMed/linking.html


Classes:
Iterator

Functions:
retrieve_entries
query_entrez

"""

import string
import urllib
import time

from Bio import File

class Iterator:
    """Returns one record at a time from a list of Medline ids.
    Queries Entrez for records, sleeping for some length of time
    between queries so as not to use too much bandwidth.

    Methods:
    next   Return the next record from the stream, or None.
    
    """
    def __init__(self, ids, delay=5.0, parser=None):
        """__init__(self, ids, delay=5.0, parser=None)

        Create a new Iterator over PubMed entries.  ids is a list of
        Medline Unique ID's or PubMed ID's.  parser is an optional
        parser (e.g. Medline.RecordParser) object to change the results
        into another form.  If set to None, then the raw contents of the
        file will be returned.  delay is the number of seconds to wait
        between each query.

        """
        self.ids = ids[:]
        self.delay = delay
        self.parser = parser
        self.last_query = None
        
    def next():
        """next() -> object

        Return the next Medline entry.  If no more entries, return None.
        Raises an IOError if there's an error.
        
        """
        if not self.ids:
            return None
        
        if self.last_query is not None:
            delay = self.last_query + self.delay - time.time()
            if delay > 0.0:
                time.sleep(delay)
        self.last_query = time.time()
        
        id = self.ids.pop(0)
        handle = query_entrez(db='m', uid=id, dopt='l').read()
        if self._parser is not None:
            return self._parser.parse(handle)
        return handle.read()

def query_entrez(db, uid, dopt, base_cgi="http://www.ncbi.nlm.nih.gov/" +
                 "htbin-post/Entrez/query"):
    """query_entrez(db, uid, dopt,
    base_cgi="http://www.ncbi.nlm.nih.gov/htbin-post/Entrez/query") -> handle

    Query Entrez and return a handle to the results.  db is the database
    to query.  uid is the ID of the thing to get.  dopt, Display Options,
    specifies the format of the return value.  base_cgi should point
    to the CGI script for querying the database.

    http://www.ncbi.nlm.nih.gov/Entrez/

    Raises an IOError exception if there's a network error.

    """
    # Format the query variables.  For some reason, the Entrez server
    # will still occasionally return HTML tags.
    options = urllib.urlencode({'db' : db,
                                'uid' : uid,
                                'Dopt' : dopt,
                                'form' : 6,      # required for some reason
                                'title' : 'no',  # no title buttons, etc.
                                'html' : 'no'    # no HTML tags
                                })
    # to do a "GET" instead of a "POST"
    #handle = urllib.urlopen("%s?%s" % (base_cgi, options))
    handle = urllib.urlopen(base_cgi, options)

    # Wrap the handle inside an UndoHandle.
    # I need to peek at the first line to see if there are errors.
    # XXX Possible bug here: what is the error message is not on the
    # first line?  Better check this.
    uhandle = File.UndoHandle(handle)

    # Check for errors
    line = uhandle.peekline()
    # Sometimes PubMed returns a Proxy Error instead of results
    if string.find(line, "500 Proxy Error") >= 0:
        raise IOError, "Proxy Error (PubMed busy?)"
    elif line[:5] == "ERROR":
        # XXX argh, need to check this
        raise IOError, "ERROR, possibly because uid not available"
    # Should I check for 404?  timeout?  etc?
    return uhandle

