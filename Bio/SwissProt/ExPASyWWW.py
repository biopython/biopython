# Copyright 2000 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""ExPASyWWW.py

This module provides code to access SwissProt at ExPASy over the WWW.
http://www.expasy.ch/sprot/


Classes:
Dictionary        Accesses a SwissProt entry using a dictionary interface.

Functions:
get_sprot_raw     Interface to the get-sprot-raw.pl CGI script.
sprot_search_ful  Interface to the sprot-search-ful CGI script.
sprot_search_de   Interface to the sprot-search-de CGI script.
_open

"""
import time
import urllib

from Bio import File

class Dictionary:
    """Access SwissProt at ExPASy using a read-only dictionary interface.

    """
    _cgi = 'http://www.expasy.ch/cgi-bin/get-sprot-raw.pl'
    
    def __init__(self, delay=5.0, parser=None):
        """__init__(self, delay=5.0, parser=None)

        Create a new Dictionary to access SwissProt.  parser is an optional
        parser (e.g. SProt.RecordParser) object to change the results
        into another form.  If set to None, then the raw contents of the
        file will be returned.  delay is the number of seconds to wait
        between each query.

        """
        self.delay = delay
        self.parser = parser
        self.last_query_time = None

    def __len__(self):
        raise NotImplementedError, "SwissProt contains lots of entries"
    def clear(self):
        raise NotImplementedError, "This is a read-only dictionary"
    def __setitem__(self, key, item):
        raise NotImplementedError, "This is a read-only dictionary"
    def update(self):
        raise NotImplementedError, "This is a read-only dictionary"
    def copy(self):
        raise NotImplementedError, "You don't need to do this..."
    def keys(self):
        raise NotImplementedError, "You don't really want to do this..."
    def items(self):
        raise NotImplementedError, "You don't really want to do this..."
    def values(self):
        raise NotImplementedError, "You don't really want to do this..."
    
    def has_key(self, id):
        """has_key(self, id) -> bool"""
        try:
            self[id]
        except KeyError:
            return 0
        return 1

    def get(self, id, failobj=None):
        try:
            return self[id]
        except KeyError:
            return failobj
        raise "How did I get here?"

    def __getitem__(self, id):
        """__getitem__(self, id) -> object

        Return a SwissProt entry.  id is either the id or accession
        for the entry.  Raises a KeyError if there's an error.
        
        """
        # First, check to see if enough time has passed since my
        # last query.
        if self.last_query_time is not None:
            delay = self.last_query_time + self.delay - time.time()
            if delay > 0.0:
                time.sleep(delay)
        self.last_query_time = time.time()

        try:
            handle = get_sprot_raw(id)
        except IOError:
            raise KeyError, id
        
        if self.parser is not None:
            return self.parser.parse(handle)
        return handle.read()

def get_sprot_raw(id, cgi='http://www.expasy.ch/cgi-bin/get-sprot-raw.pl'):
    """get_sprot_raw(id, cgi='http://www.expasy.ch/cgi-bin/get-sprot-raw.pl')
    -> handle

    Get a handle to a raw SwissProt entry at ExPASy.

    """
    return _open("%s?%s" % (cgi, id))

def sprot_search_ful(text, make_wild=None, swissprot=1, trembl=None,
                     cgi='http://www.expasy.ch/cgi-bin/sprot-search-ful'):
    """sprot_search_ful(text, make_wild=None, swissprot=1, trembl=None,
    cgi='http://www.expasy.ch/cgi-bin/sprot-search-ful') -> handle

    Search SwissProt by full text.

    """
    variables = {'SEARCH' : text}
    if make_wild:
        variables['makeWild'] = 'on'
    if swissprot:
        variables['S'] = 'on'
    if trembl:
        variables['T'] = 'on'
    return _open(cgi, variables)

def sprot_search_de(text, swissprot=1, trembl=None,
                    cgi='http://www.expasy.ch/cgi-bin/sprot-search-de'):
    """sprot_search_de(text, swissprot=1, trembl=None,
    cgi='http://www.expasy.ch/cgi-bin/sprot-search-de') -> handle

    Search SwissProt by name, description, gene name, species, or
    organelle.

    """
    variables = {'SEARCH' : text}
    if swissprot:
        variables['S'] = 'on'
    if trembl:
        variables['T'] = 'on'
    return _open(cgi, variables)

def _open(cgi, params={}, get=1):
    """_open(cgi, params={}, get=1) -> UndoHandle

    Open a handle to ExPASy.  cgi is the URL for the cgi script to access.
    params is a dictionary with the options to pass to it.  get is a boolean
    that describes whether a GET should be used.  Does some
    simple error checking, and will raise an IOError if it encounters one.

    """
    # Open a handle to ExPASy.
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

    # If the key doesn't exist, ExPASy returns nothing.
    if not uhandle.peekline():
        raise IOError, "no results"
    
    return uhandle
