# Copyright 2000 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Genbank/__init__.py

This module provides code to work with GenBank.
http://www.ncbi.nlm.nih.gov/

Classes:
NCBIDictionary    Access GenBank using a dictionary interface.

Functions:
search_for        Do a query against GenBank.

"""

import re
import sgmllib
import urlparse

from Bio.WWW import NCBI
from Bio.WWW import RequestLimiter

class NCBIDictionary:
    """Access GenBank using a read-only dictionary interface.

    Methods:
    
    """
    def __init__(self, database='Nucleotide', delay=5.0, parser=None):
        """NCBIDictionary([database][, delay][, parser])

        Create a new Dictionary to access GenBank.  database should be
        either 'Nucleotide' or 'Protein'.  delay is the number of
        seconds to wait between each query (5 default).  parser is an
        optional parser object to change the results into another
        form.  If unspecified, then the raw contents of the file will
        be returned.

        """
        self.parser = parser
        self.limiter = RequestLimiter(delay)
        if database == 'Nucleotide':
            self.format = 'GenBank'
        elif database == 'Protein':
            self.format = 'GenPept'
        else:
            raise ValueError, "database should be 'Nucleotide' or 'Protein'."
        self.database = database

    def __len__(self):
        raise NotImplementedError, "GenBank contains lots of entries"
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
        """S.has_key(id) -> bool"""
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
        """S.__getitem__(id) -> object

        Return the GenBank entry.  id is the GenBank ID (gi) of the
        entry.  Raises a KeyError if there's an error.
        
        """
        # First, check to see if enough time has passed since my
        # last query.
        self.limiter.wait()
        
        try:
            handle = NCBI.query(
                'Text', self.database, dopt=self.format, uid=id)
        except IOError, x:
            # raise a KeyError instead of an IOError
            # XXX I really should distinguish between a real IOError and
            # if the id is not in the database.
            raise KeyError, x
        # If the id is not in the database, I get a message like:
        # 'GenPept does not exist for GI "433174"\012\012' 
        line = handle.peekline()
        if line.find('does not exist') >= 0:
            raise KeyError, line
        elif line.lower().find('html') >= 0:
            raise KeyError, "I unexpectedly got back html-formatted data."
        # Parse the record if a parser was passed in.
        if self.parser is not None:
            return self.parser.parse(handle)
        return handle.read()

def search_for(search, database='Nucleotide', max_ids=500):
    """search_for(search[, database][, max_ids])

    Search GenBank and return a list of GenBank identifiers (gi's).
    search is the search string used to search the database.  database
    should be either 'Nucleotide' or 'Protein'.  max_ids is the maximum
    number of ids to retrieve (default 500).
    
    """
    if database not in ['Nucleotide', 'Protein']:
        raise ValueError, "database must be 'Nucleotide' or 'Protein'"
    
    class ResultParser(sgmllib.SGMLParser):
        # GenBank returns an HTML formatted page with lots of pretty stuff.
        # I want to rip out all the genbank id's from the page.  I'm going
        # to do this by looking for links that retrieve records using
        # the query.fcgi script.
        # <a href="http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?
        # cmd=Retrieve&amp;db=Nucleotide&amp;list_uids=5174616&amp;
        # dopt=GenBank">NM_006092</a>
        def __init__(self):
            sgmllib.SGMLParser.__init__(self)
            self.ids = []
        def start_a(self, attrs):
            # If I see a href back to the Nucleotide database from
            # Entrez, then keep it.
            href = None
            for name, value in attrs:
                if name == 'href':
                    href = value
            if not href:
                return
            scheme, netloc, path, params, query, frag = urlparse.urlparse(href)
            if path[-10:] != 'query.fcgi':   # only want links to query.fcgi
                return
            # Valid queries look like (truncated):
            # cmd=Retrieve&amp;db=Nucleotide&amp;list_uids=4503354&amp;dopt=Gen
            # I have to first split on '&amp;' and then on '='.
            params = query.split('&amp;')
            params = [x.split('=') for x in params]
            list_uids = None
            db = None
            for name, value in params:
                if name == 'list_uids':
                    list_uids = value
            if list_uids is not None:
                self.ids.append(list_uids)

    parser = ResultParser()
    handle = NCBI.query("Search", database, term=search, doptcmdl='DocSum',
                        dispmax=max_ids)
    parser.feed(handle.read())
    return parser.ids
