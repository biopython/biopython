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
Dictionary

"""

import time

from Bio.WWW import NCBI

class Dictionary:
    """Access PubMed using a read-only dictionary interface.

    Methods:
    
    """
    def __init__(self, delay=5.0, parser=None):
        """__init__(self, delay=5.0, parser=None)

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
            raise KeyError, x
        if self.parser is not None:
            return self.parser.parse(handle)
        return handle.read()
