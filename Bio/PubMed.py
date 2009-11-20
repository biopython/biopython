# Copyright 1999-2000 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
This module provides code to work with PubMed from the NCBI (DEPRECATED).

This module has been deprecated and is likely to be removed in a future
release of Biopython.  Please use Bio.Entrez instead, which is described
in the Biopython Tutorial.

See also:
http://www.ncbi.nlm.nih.gov/PubMed/

Online documentation for linking to PubMed is available at:
http://www.ncbi.nlm.nih.gov/PubMed/linking.html


Classes:
Dictionary     Access PubMed articles using a dictionary interface.

Functions:
search_for     Search PubMed.
find_related   Find related articles in PubMed.
download_many  Download many articles from PubMed in batch mode.

"""

import warnings
warnings.warn("Bio.PubMed has been deprecated, and we intend to remove it in" \
              +" a future release of Biopython.  Please use Bio.Entrez"\
              +" instead as described in the Tutorial.  If you need help" \
              +" with this transition, or wish to continue to use this code,"\
              +" please get in contact via the mailing lists.", \
              DeprecationWarning)

import re
import sgmllib

from Bio import File
from Bio import Entrez
from Bio import Medline

class Dictionary:
    """Access PubMed using a read-only dictionary interface (DEPRECATED).

    Please use the Bio.Entrez.efetch(...) function instead as described in the
    Biopython Tutorial.
    """
    def __init__(self, parser=None):
        """Dictionary(parser=None)

        Create a new Dictionary to access PubMed.  parser is an optional
        parser (e.g. Medline.RecordParser) object to change the results
        into another form.  If set to None, then the raw contents of the
        file will be returned.

        """
        self.parser = parser

    def __len__(self):
        raise NotImplementedError("PubMed contains lots of entries")
    def clear(self):
        raise NotImplementedError("This is a read-only dictionary")
    def __setitem__(self, key, item):
        raise NotImplementedError("This is a read-only dictionary")
    def update(self):
        raise NotImplementedError("This is a read-only dictionary")
    def copy(self):
        raise NotImplementedError("You don't need to do this...")
    def keys(self):
        raise NotImplementedError("You don't really want to do this...")
    def items(self):
        raise NotImplementedError("You don't really want to do this...")
    def values(self):
        raise NotImplementedError("You don't really want to do this...")
    
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

    def __getitem__(self, id):
        """S.__getitem__(id) -> object

        Return the Medline entry.  id is either the Medline Unique ID
        or the Pubmed ID of the article.  Raises a KeyError if there's an
        error.
        
        """
        try:
            handle = Entrez.efetch(
                db="pubmed", id=id, retmode='text', rettype='medlars')
        except IOError, x:
            # raise a KeyError instead of an IOError
            # XXX I really should distinguish between a real IOError and
            # if the id is not in the database.
            raise KeyError(x)
        if self.parser is not None:
            return self.parser.parse(handle)
        return handle.read()

def search_for(search, reldate=None, mindate=None, maxdate=None,
               batchsize=100, callback_fn=None, start_id=0, max_ids=None):
    """Search PubMed, returns a list of IDs (DEPRECATED).

    Please use Bio.Entrez instead as described in the Biopython Tutorial.

    Search PubMed and return a list of the PMID's that match the
    criteria.  search is the search string used to search the
    database.  reldate is the number of dates prior to the current
    date to restrict the search.  mindate and maxdate are the dates to
    restrict the search, e.g. 2002/01/01.  batchsize specifies the
    number of ids to return at one time.  By default, it is set to
    10000, the maximum.  callback_fn is an optional callback function
    that will be called as passed a PMID as results are retrieved.
    start_id specifies the index of the first id to retrieve and
    max_ids specifies the maximum number of id's to retrieve.

    XXX The date parameters don't seem to be working with NCBI's
    script.  Please let me know if you can get it to work.
    
    """
    params = {
        'db' : 'pubmed',
        'term' : search,
        'reldate' : reldate,
        'mindate' : mindate,
        'maxdate' : maxdate
        }
    #Note that Bio.Entrez can now cope with None arguments (it ignores them)

    ids = []
    while max_ids is None or len(ids) < max_ids:
        start = start_id + len(ids)
        max = batchsize
        if max_ids is not None and max > max_ids - len(ids):
            max = max_ids - len(ids)

        params['retstart'] = start
        params['retmax'] = max
        h = Entrez.esearch(**params)
        record = Entrez.read(h)
        idlist = record["IdList"]
        ids.extend(idlist)
        if callback_fn is not None:
            # Call the callback function with each of the new ID's.
            for id in idlist:
                callback_fn(id)
        if len(idlist) < max:  # no more id's to read
            break
    return ids

def find_related(pmid):
    """Find related articles in PubMed, returns an ID list (DEPRECATED).

    Search PubMed for a list of citations related to pmid.  pmid can
    be a PubMed ID, a MEDLINE UID, or a list of those.

    Please use Bio.Entrez instead as described in the Biopython Tutorial.
    e.g.

    >>> from Bio import Entrez
    >>> Entrez.email = "A.N.Other@example.com"
    >>> pmid = "12230038"
    >>> handle = Entrez.elink(dbfrom='pubmed', id=pmid)
    >>> result = Entrez.read(handle)
    >>> for link in result[0]["LinkSetDb"][0]['Link']:
    ...     print link

    (Output ommitted)

    """
    class ResultParser(sgmllib.SGMLParser):
        # Parse the ID's out of the HTML-formatted page that PubMed
        # returns.  The format of the page is:
        # [...]
        #   <Link>
        #      <Id>######</Id>
        #      <Score>######</Score>
        #      [...]
        #   </Link>
        # [...]
        def __init__(self):
            sgmllib.SGMLParser.__init__(self)
            self.ids = []
            self.in_link = 0
            self.in_id = 0
        def start_id(self, attributes):
            self.in_id = 1
        def end_id(self):
            self.in_id = 0
        def start_link(self, attributes):
            self.in_link = 1
        def end_link(self):
            self.in_link = 0
        _not_pmid_re = re.compile(r'\D')
        def handle_data(self, data):
            if not self.in_link or not self.in_id:
                return
            # Everything here should be a PMID.  Check and make sure
            # data really is one.  A PMID should be a string consisting
            # of only integers.  Should I check to make sure it
            # meets a certain minimum length?
            if self._not_pmid_re.search(data):
                raise ValueError(\
                      "I expected an ID, but '%s' doesn't look like one." % \
                      repr(data))
            self.ids.append(data)

    parser = ResultParser()
    if type(pmid) is type([]):
        pmid = ','.join(pmid)
    h = Entrez.elink(dbfrom='pubmed', id=pmid)
    parser.feed(h.read())
    return parser.ids

def download_many(ids, callback_fn, broken_fn=None, 
                  batchsize=500, parser=None):
    """Download multiple PubMed records, no return value (DEPRECATED).

    Please use Bio.Entrez instead as described in the Biopython Tutorial.

    Download many records from PubMed.  ids is a list of either the
    Medline Unique ID or the PubMed ID's of the articles.  Each time a
    record is downloaded, callback_fn is called with the text of the
    record.  broken_fn is an optional function that is called with the
    id of records that were not able to be downloaded.  batchsize is the
    number of records to request each time.

    """
    # parser is an undocumented parameter that allows people to
    # specify an optional parser to handle each record.  This is
    # dangerous because the results may be malformed, and exceptions
    # in the parser may disrupt the whole download process.
    if batchsize > 500 or batchsize < 1:
        raise ValueError("batchsize must be between 1 and 500")
    current_batchsize = batchsize
    
    # Loop until all the ids are processed.  We want to process as
    # many as possible with each request.  Unfortunately, errors can
    # occur.  Some id may be incorrect, or the server may be
    # unresponsive.  In addition, one broken id out of a list of id's
    # can cause a non-specific error.  Thus, the strategy I'm going to
    # take, is to start by downloading as many as I can.  If the
    # request fails, I'm going to half the number of records I try to
    # get.  If there's only one more record, then I'll report it as
    # broken and move on.  If the request succeeds, I'll double the
    # number of records until I get back up to the batchsize.
    nsuccesses = 0
    while ids:
        if current_batchsize > len(ids):
            current_batchsize = len(ids)
        
        id_str = ','.join(ids[:current_batchsize])

        try:
            # Query PubMed.  If one or more of the id's are broken,
            # this will raise an IOError.
            handle = Entrez.efetch(
                db="pubmed", id=id_str, retmode='text', rettype='medlars')

            # I'm going to check to make sure PubMed returned the same
            # number of id's as I requested.  If it didn't then I'm going
            # to raise an exception.  This could take a lot of memory if
            # the batchsize is large.
            results = handle.read()
            num_ids = 0
            for x in Medline.Iterator(File.StringHandle(results)):
                num_ids = num_ids + 1
            if num_ids != current_batchsize:
                raise IOError
            handle = File.StringHandle(results)
        except IOError:   # Query did not work.
            if current_batchsize == 1:
                # There was only 1 id in the query.  Report it as
                # broken and move on.
                id = ids.pop(0)
                if broken_fn is not None:
                    broken_fn(id)
            else:
                # I don't know which one is broken.  Try again with
                # fewer id's.
                current_batchsize = current_batchsize / 2
            nsuccesses = 0
            continue
        nsuccesses = nsuccesses + 1

        # Iterate through the results and pass the records to the
        # callback.
        idnum = 0
        for rec in Medline.Iterator(handle, parser):
            callback_fn(ids[idnum], rec)
            idnum = idnum + 1

        ids = ids[current_batchsize:]

        # If I'm not downloading the maximum number of articles,
        # double the number for next time.
        if nsuccesses >= 2 and current_batchsize < batchsize:
            current_batchsize = current_batchsize * 2
            if current_batchsize > batchsize:
                current_batchsize = batchsize
