# Copyright 1999-2000 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""WWW.py

Provides code to access the Entrez over the WWW.


Functions:
query

"""

def query(db, uid, dopt, base_cgi="http://www.ncbi.nlm.nih.gov/" +
                 "htbin-post/Entrez/query"):
    """query(db, uid, dopt,
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
    # Sometimes Entrez returns a Proxy Error instead of results
    if string.find(line, "500 Proxy Error") >= 0:
        raise IOError, "Proxy Error (Entrez busy?)"
    elif line[:5] == "ERROR":
        # XXX argh, need to check this
        raise IOError, "ERROR, possibly because uid not available"
    # Should I check for 404?  timeout?  etc?
    return uhandle
