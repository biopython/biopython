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
Dictionary    Access PubMed articles using a dictionary interface.

Functions:
search_for    Search PubMed.

"""

import time
import string
import re
import sgmllib

from Bio.WWW import NCBI

class Dictionary:
    """Access PubMed using a read-only dictionary interface.

    Methods:
    
    """
    def __init__(self, delay=5.0, parser=None):
        """Dictionary(delay=5.0, parser=None)

        Create a new Dictionary to access PubMed.  parser is an optional
        parser (e.g. Medline.RecordParser) object to change the results
        into another form.  If set to None, then the raw contents of the
        file will be returned.  delay is the number of seconds to wait
        between each query.

        """
        self.delay = delay
        self.parser = parser
        self.last_query_time = None

    def __len__(self):
        raise NotImplementedError, "PubMed contains lots of entries"
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

        Return the Medline entry.  id is either the Medline Unique ID
        or the Pubmed ID of the article.  Raises a KeyError if there's an
        error.
        
        """
        # First, check to see if enough time has passed since my
        # last query.
        if self.last_query_time is not None:
            delay = self.last_query_time + self.delay - time.time()
            if delay > 0.0:
                time.sleep(delay)
        self.last_query_time = time.time()
        
        try:
            handle = NCBI.pmfetch(
                db='PubMed', id=id, report='medlars', mode='text')
        except IOError, x:
            # raise a KeyError instead of an IOError
            # XXX I really should distinguish between a real IOError and
            # if the id is not in the database.
            raise KeyError, x
        if self.parser is not None:
            return self.parser.parse(handle)
        return handle.read()

def search_for(search, batchsize=10000, delay=1, callback_fn=None):
    """search_for(search, batchsize=10000, delay=1, callback_fn=None) -> ids

    Search PubMed and return a list of the PMID's that match the criteria.
    search is the search string used to search the database.  PubMed only
    allows users to retrieve the search results in batches of up to 10000
    ID's at a time.  batchsize is the size of the batch to use.  delay
    is the number of seconds to wait between queries.  callback_fn is
    an optional callback function that will be called as results are
    retrieved.  It should take the PMID as an argument.

    """
    class ResultParser(sgmllib.SGMLParser):
        # Parse the ID's out of the HTML-formatted page that PubMed
        # returns.  The format of the page is:
        # <Title>QueryResult</Title>
        # <Body>
        # 10807727<Br>
        # [...]
        # </Body>
        def __init__(self):
            sgmllib.SGMLParser.__init__(self)
            self.ids = []
            self.in_body = 0
        def start_body(self, attributes):
            self.in_body = 1
        def end_body(self):
            self.in_body = 0
        _not_pmid_re = re.compile(r'\D')
        def handle_data(self, data):
            # The ID's only appear in the body.  If I'm not in the body,
            # then don't do anything.
            if not self.in_body:
                return
            # If data is just whitespace, then ignore it.
            data = string.strip(data)
            if not data:
                return
            # Everything here should be a PMID.  Check and make sure
            # data really is one.  A PMID should be a string consisting
            # of only integers.  Should I check to make sure it
            # meets a certain minimum length?
            if self._not_pmid_re.search(data):
                raise SyntaxError, \
                      "I expected an ID, but '%s' doesn't look like one." % \
                      repr(data)
            self.ids.append(data)

    last_search = None
    ids = []
    while 1:
        parser = ResultParser()
        
        # Check to make sure enough time has passed before my
        # last search.  If not, then wait.
        if last_search is not None:
            time.sleep(time.time() - (last_search + delay))
        last_search = time.time()
        
        # Do a query against PmQty.  Search medline, using the
        # search string, and get only the ID's in the results.
        h = NCBI.pmqty('m', search, dopt='d',
                       dispmax=batchsize, dispstart=len(ids))
        parser.feed(h.read())
        if not parser.ids:  # no more id's to read
            break
        if callback_fn is not None:
            # Call the callback function with each of the new ID's.
            for id in parser.ids:
                callback_fn(id)
        ids.extend(parser.ids)
    return ids
