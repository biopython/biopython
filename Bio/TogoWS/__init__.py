# Copyright 2010-2011 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Provides code to access the TogoWS integrated websevices of DBCLS, Japan.

This module aims to make the TogoWS (from DBCLS, Japan) easier to use. See:
http://togows.dbcls.jp/

The TogoWS REST service provides simple access to a range of databases, acting
as a proxy to shield you from all the different provider APIs. This works using
simple URLs (which this module will construct for you). For more details, see
http://togows.dbcls.jp/site/en/rest.html

The functionality is somewhat similar to Biopython's Bio.Entrez module which
provides access to the NCBI's Entrez Utilities (E-Utils) which also covers a
wide range of databases.

Currently TogoWS does not provide any usage guidelines (unlike the NCBI whose
requirements are reasonably clear). To avoid risking overloading the service,
Biopython will only allow three calls per second.

The TogoWS SOAP service offers a more complex API for calling web services
(essentially calling remote functions) provided by DDBJ, KEGG and PDBj. For
example, this allows you to run a remote BLAST search at the DDBJ. This is
not yet covered by this module, however there are lots of Python examples
on the TogoWS website using the SOAPpy python library. See:
http://togows.dbcls.jp/site/en/soap.html
http://soapy.sourceforge.net/
"""

import urllib
import urllib2
import time
from Bio._py3k import _binary_to_string_handle, _as_bytes

#Constant
_BASE_URL = "http://togows.dbcls.jp"

#Caches:
_search_db_names = None
_entry_db_names = None
_entry_db_fields = {}
_entry_db_formats = {}
_convert_formats = []

def _get_fields(url):
    """Queries a TogoWS URL for a plain text list of values (PRIVATE)."""
    handle = _open(url)
    fields = handle.read().strip().split()
    handle.close()
    return fields

def _get_entry_dbs():
    return _get_fields(_BASE_URL + "/entry")

def _get_entry_fields(db):
    return _get_fields(_BASE_URL + "/entry/%s?fields" % db)

def _get_entry_formats(db):
    return _get_fields(_BASE_URL + "/entry/%s?formats" % db)

def _get_convert_formats():
    return [pair.split(".") for pair in \
            _get_fields(_BASE_URL + "/convert/")]

def entry(db, id, format=None, field=None):
    """TogoWS fetch entry (returns a handle).

    db - database (string), see list below.
    id - identier (string) or a list of identifiers (either as a list of
         strings or a single string with comma separators).
    format - return data file format (string), options depend on the database
             e.g. "xml", "json", "gff", "fasta", "ttl" (RDF Turtle)
    field - specific field from within the database record (string)
            e.g. "au" or "authors" for pubmed.

    At the time of writing, this includes the following:

    KEGG: compound, drug, enzyme, genes, glycan, orthology, reaction,
          module, pathway
    DDBj: ddbj, dad, pdb
    NCBI: nuccore, nucest, nucgss, nucleotide, protein, gene, onim,
          homologue, snp, mesh, pubmed
    EBI:  embl, uniprot, uniparc, uniref100, uniref90, uniref50

    For the current list, please see http://togows.dbcls.jp/entry/

    This function is essentially equivalent to the NCBI Entrez service
    EFetch, available in Biopython as Bio.Entrez.efetch(...), but that
    does not offer field extraction.
    """
    global _entry_db_names, _entry_db_fields, fetch_db_formats
    if _entry_db_names is None:
        _entry_db_names = _get_entry_dbs()
    if db not in _entry_db_names:
        raise ValueError("TogoWS entry fetch does not officially support "
                         "database '%s'." % db)
    if field:
        try:
            fields = _entry_db_fields[db]
        except KeyError:
            fields = _get_entry_fields(db)
            _entry_db_fields[db] = fields
        if field not in fields:
            raise ValueError("TogoWS entry fetch does not explicitly support "
                             "field '%s' for database '%s'. Only: %s" \
                             % (field, db, ", ".join(sorted(fields))))
    if format:
        try:
            formats = _entry_db_formats[db]
        except KeyError:
            formats = _get_entry_formats(db)
            _entry_db_formats[db] = formats
        if format not in formats:
            raise ValueError("TogoWS entry fetch does not explicitly support "
                             "format '%s' for database '%s'. Only: %s" \
                             % (format, db, ", ".join(sorted(formats))))

    if isinstance(id, list):
        id = ",".join(id)
    url = _BASE_URL + "/entry/%s/%s" % (db, urllib.quote(id))
    if field:
        url += "/" + field
    if format:
        url += "." + format
    return _open(url)

def search_count(db, query):
    """TogoWS search count (returns an integer).

    db - database (string), see http://togows.dbcls.jp/search
    query - search term (string)

    You could then use the count to download a large set of search results in
    batches using the offset and limit options to Bio.TogoWS.search(). In
    general however the Bio.TogoWS.search_iter() function is simpler to use.
    """
    global _search_db_names
    if _search_db_names is None:
        _search_db_names = _get_fields(_BASE_URL + "/search")
    if db not in _search_db_names:
        #TODO - Make this a ValueError? Right now despite the HTML website
        #claiming to, the "gene" or "ncbi-gene" don't work and are not listed.
        import warnings
        warnings.warn("TogoWS search does not officially support database '%s'. "
                      "See %s/search/ for options." % (db, _BASE_URL))
    handle = _open(_BASE_URL + "/search/%s/%s/count" \
                   % (db, urllib.quote(query)))
    count = int(handle.read().strip())
    handle.close()
    return count

def search_iter(db, query, limit=None, batch=100):
    """TogoWS search iteratating over the results (generator function).

    db - database (string), see http://togows.dbcls.jp/search
    query - search term (string)
    limit - optional upper bound on number of search results
    batch - number of search results to pull back each time talk to
            TogoWS (currently limited to 100).

    You would use this function within a for loop, e.g.

    >>> for id in search_iter("pubmed", "lung+cancer+drug", limit=10):
    ...     print id #maybe fetch data with entry?

    Internally this first calls the Bio.TogoWS.search_count() and then
    uses Bio.TogoWS.search() to get the results in batches.
    """
    count = search_count(db, query)
    if not count:
        raise StopIteration
    #NOTE - We leave it to TogoWS to enforce any upper bound on each
    #batch, they currently return an HTTP 400 Bad Request if above 100.
    remain = count
    if limit is not None:
        remain = min(remain, limit)
    offset = 1 #They don't use zero based counting
    prev_ids = [] #Just cache the last batch for error checking
    while remain:
        batch = min(batch, remain)
        #print "%r left, asking for %r" % (remain, batch)
        ids = search(db, query, offset, batch).read().strip().split()
        assert len(ids)==batch, "Got %i, expected %i" % (len(ids), batch)
        #print "offset %i, %s ... %s" % (offset, ids[0], ids[-1])
        if ids == prev_ids:
            raise RuntimeError("Same search results for previous offset")
        for identifier in ids:
            if identifier in prev_ids:
                raise RuntimeError("Result %s was in previous batch" \
                                   % identifier)
            yield identifier
        offset += batch
        remain -= batch
        prev_ids = ids

def search(db, query, offset=None, limit=None, format=None):
    """TogoWS search (returns a handle).

    This is a low level wrapper for the TogoWS search function, which
    can return results in a several formats. In general, the search_iter
    function is more suitable for end users.

    db - database (string), see http://togows.dbcls.jp/search/
    query - search term (string)
    offset, limit - optional integers specifying which result to start from
            (1 based) and the number of results to return.
    format - return data file format (string), e.g. "json", "ttl" (RDF)
             By default plain text is returned, one result per line.

    At the time of writing, TogoWS applies a default count limit of 100
    search results, and this is an upper bound. To access more results,
    use the offset argument or the search_iter(...) function.

    TogoWS supports a long list of databases, including many from the NCBI
    (e.g. "ncbi-pubmed" or "pubmed", "ncbi-genbank" or "genbank", and
    "ncbi-taxonomy"), EBI (e.g. "ebi-ebml" or "embl", "ebi-uniprot" or
    "uniprot, "ebi-go"), and KEGG (e.g. "kegg-compound" or "compound").
    For the current list, see http://togows.dbcls.jp/search/

    The NCBI provide the Entrez Search service (ESearch) which is similar,
    available in Biopython as the Bio.Entrez.esearch() function.

    See also the function Bio.TogoWS.search_count() which returns the number
    of matches found, and the Bio.TogoWS.search_iter() function which allows
    you to iterate over the search results (taking care of batching for you).
    """
    global _search_db_names
    if _search_db_names is None:
        _search_db_names = _get_fields(_BASE_URL + "/search")
    if db not in _search_db_names:
        #TODO - Make this a ValueError? Right now despite the HTML website
        #claiming to, the "gene" or "ncbi-gene" don't work and are not listed.
        import warnings
        warnings.warn("TogoWS search does not explicitly support database '%s'. "
                      "See %s/search/ for options." % (db, _BASE_URL))
    url = _BASE_URL + "/search/%s/%s" % (db, urllib.quote(query))
    if offset is not None and limit is not None:
        try:
            offset = int(offset)
        except:
            raise ValueError("Offset should be an integer (at least one), not %r" % offset)
        try:
            limit = int(limit)
        except:
            raise ValueError("Limit should be an integer (at least one), not %r" % limit)
        if offset <= 0:
            raise ValueError("Offset should be at least one, not %i" % offset)
        if limit <= 0:
            raise ValueError("Count should be at least one, not %i" % limit)
        url += "/%i,%i" % (offset, limit)
    elif offset is not None or limit is not None:
        raise ValueError("Expect BOTH offset AND limit to be provided (or neither)")
    if format:
        url += "." + format
    #print url
    return _open(url)

def convert(data, in_format, out_format):
    """TogoWS convert (returns a handle).
    
    data - string or handle containing input record(s)
    in_format - string describing the input file format (e.g. "genbank")
    out_format - string describing the requested output format (e.g. "fasta")

    For a list of supported conversions (e.g. "genbank" to "fasta"), see
    http://togows.dbcls.jp/convert/
    
    Note that Biopython has built in support for conversion of sequence and
    alignnent file formats (functions Bio.SeqIO.convert and Bio.AlignIO.convert)
    """
    global _convert_formats
    if not _convert_formats:
        _convert_formats = _get_convert_formats()
    if [in_format, out_format] not in _convert_formats:
        msg = "\n".join("%s -> %s" % tuple(pair) for pair in _convert_formats)
        raise ValueError("Unsupported conversion. Choose from:\n%s" % msg)
    url = _BASE_URL + "/convert/%s.%s" % (in_format, out_format)
    #TODO - Should we just accept a string not a handle? What about a filename?
    if hasattr(data, "read"):
        #Handle
        return _open(url, post={"data":data.read()})
    else:
        #String
        return _open(url, post={"data":data})

def _open(url, post=None):
    """Helper function to build the URL and open a handle to it (PRIVATE).

    Open a handle to TogoWS, will raise an IOError if it encounters an error.

    In the absense of clear guidelines, this function enforces a limit of
    "up to three queries per second" to avoid abusing the TogoWS servers.
    """
    delay = 0.333333333 #one third of a second
    current = time.time()
    wait = _open.previous + delay - current
    if wait > 0:
        time.sleep(wait)
        _open.previous = current + wait
    else:
        _open.previous = current

    #print url
    try:
        if post:
            handle = urllib2.urlopen(url, _as_bytes(urllib.urlencode(post)))
        else:
            handle = urllib2.urlopen(url)
    except urllib2.HTTPError, exception:
        raise exception

    #We now trust TogoWS to have set an HTTP error code, that
    #suffices for my current unit tests. Previously we would
    #examine the start of the data returned back.
    return _binary_to_string_handle(handle)

_open.previous = 0
