# Copyright 2000 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
Provides code to access SCOP over the WWW.  The main SCOP web page
is available at:
http://scop.mrc-lmb.cam.ac.uk/scop/

Functions:
search       Access the main CGI script.
_open

"""
import string
import urllib

from Bio import File

def search(pdb=None, key=None, sid=None, disp=None, dir=None, loc=None,
           cgi='http://scop.mrc-lmb.cam.ac.uk/scop/search.cgi', **keywds):
    """search(pdb=None, key=None, sid=None, disp=None, dir=None, loc=None,
    cgi='http://scop.mrc-lmb.cam.ac.uk/scop/search.cgi', **keywds)

    Access search.cgi and return a handle to the results.  See the
    online help file for an explanation of the parameters:
    http://scop.mrc-lmb.cam.ac.uk/scop/help.html

    Raises an IOError if there's a network error.
    
    """
    params = {'pdb' : pdb, 'key' : key, 'sid' : sid, 'disp' : disp,
              'dir' : dir, 'loc' : loc}
    variables = {}
    for k in params.keys():
        if params[k] is not None:
            variables[k] = params[k]
    variables.update(keywds)
    return _open(cgi, variables)

def _open(cgi, params={}, get=1):
    """_open(cgi, params={}, get=1) -> UndoHandle

    Open a handle to SCOP.  cgi is the URL for the cgi script to access.
    params is a dictionary with the options to pass to it.  get is a boolean
    that describes whether a GET should be used.  Does some
    simple error checking, and will raise an IOError if it encounters one.

    """
    # Open a handle to SCOP.
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
    # Should I check for 404?  timeout?  etc?
    return uhandle
