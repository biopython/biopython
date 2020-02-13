# Copyright 2010-2011, 2013-2014, 2016-2018 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
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


import io
import time

from urllib.request import urlopen
from urllib.parse import quote


# Constant
_BASE_URL = "http://togows.dbcls.jp"

# Caches:
_search_db_names = None
_entry_db_names = None
_entry_db_fields = {}
_entry_db_formats = {}
_convert_formats = []


def _get_fields(url):
    """Query a TogoWS URL for a plain text list of values (PRIVATE)."""
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
    return [pair.split(".") for pair in _get_fields(_BASE_URL + "/convert/")]


def entry(db, id, format=None, field=None):
    """Call TogoWS 'entry' to fetch a record.

    Arguments:
     - db - database (string), see list below.
     - id - identier (string) or a list of identifiers (either as a list of
       strings or a single string with comma separators).
     - format - return data file format (string), options depend on the database
       e.g. "xml", "json", "gff", "fasta", "ttl" (RDF Turtle)
     - field - specific field from within the database record (string)
       e.g. "au" or "authors" for pubmed.

    At the time of writing, this includes the following::

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
        raise ValueError(
            "TogoWS entry fetch does not officially support database '%s'." % db
        )
    if field:
        try:
            fields = _entry_db_fields[db]
        except KeyError:
            fields = _get_entry_fields(db)
            _entry_db_fields[db] = fields
        if db == "pubmed" and field == "ti" and "title" in fields:
            # Backwards compatibility fix for TogoWS change Nov/Dec 2013
            field = "title"
            import warnings

            warnings.warn(
                "TogoWS dropped 'pubmed' field alias 'ti', please use 'title' instead."
            )
        if field not in fields:
            raise ValueError(
                "TogoWS entry fetch does not explicitly support "
                "field '%s' for database '%s'. Only: %s"
                % (field, db, ", ".join(sorted(fields)))
            )
    if format:
        try:
            formats = _entry_db_formats[db]
        except KeyError:
            formats = _get_entry_formats(db)
            _entry_db_formats[db] = formats
        if format not in formats:
            raise ValueError(
                "TogoWS entry fetch does not explicitly support "
                "format '%s' for database '%s'. Only: %s"
                % (format, db, ", ".join(sorted(formats)))
            )

    if isinstance(id, list):
        id = ",".join(id)
    url = _BASE_URL + "/entry/%s/%s" % (db, quote(id))
    if field:
        url += "/" + field
    if format:
        url += "." + format
    return _open(url)


def search_count(db, query):
    """Call TogoWS search count to see how many matches a search gives.

    Arguments:
     - db - database (string), see http://togows.dbcls.jp/search
     - query - search term (string)

    You could then use the count to download a large set of search results in
    batches using the offset and limit options to Bio.TogoWS.search(). In
    general however the Bio.TogoWS.search_iter() function is simpler to use.
    """
    global _search_db_names
    if _search_db_names is None:
        _search_db_names = _get_fields(_BASE_URL + "/search")
    if db not in _search_db_names:
        # TODO - Make this a ValueError? Right now despite the HTML website
        # claiming to, the "gene" or "ncbi-gene" don't work and are not listed.
        import warnings

        warnings.warn(
            "TogoWS search does not officially support database '%s'. "
            "See %s/search/ for options." % (db, _BASE_URL)
        )
    url = _BASE_URL + "/search/%s/%s/count" % (db, quote(query))
    handle = _open(url)
    data = handle.read()
    handle.close()
    if not data:
        raise ValueError("TogoWS returned no data from URL %s" % url)
    try:
        return int(data.strip())
    except ValueError:
        raise ValueError(
            "Expected an integer from URL %s, got: %r" % (url, data)
        ) from None


def search_iter(db, query, limit=None, batch=100):
    """Call TogoWS search iterating over the results (generator function).

    Arguments:
     - db - database (string), see http://togows.dbcls.jp/search
     - query - search term (string)
     - limit - optional upper bound on number of search results
     - batch - number of search results to pull back each time talk to
       TogoWS (currently limited to 100).

    You would use this function within a for loop, e.g.

    >>> from Bio import TogoWS
    >>> for id in TogoWS.search_iter("pubmed", "diabetes+human", limit=10):
    ...     print("PubMed ID: %s" %id) # maybe fetch data with entry?
    PubMed ID: ...

    Internally this first calls the Bio.TogoWS.search_count() and then
    uses Bio.TogoWS.search() to get the results in batches.
    """
    count = search_count(db, query)
    if not count:
        return
    # NOTE - We leave it to TogoWS to enforce any upper bound on each
    # batch, they currently return an HTTP 400 Bad Request if above 100.
    remain = count
    if limit is not None:
        remain = min(remain, limit)
    offset = 1  # They don't use zero based counting
    prev_ids = []  # Just cache the last batch for error checking
    while remain:
        batch = min(batch, remain)
        # print("%r left, asking for %r" % (remain, batch))
        ids = search(db, query, offset, batch).read().strip().split()
        assert len(ids) == batch, "Got %i, expected %i" % (len(ids), batch)
        # print("offset %i, %s ... %s" % (offset, ids[0], ids[-1]))
        if ids == prev_ids:
            raise RuntimeError("Same search results for previous offset")
        for identifier in ids:
            if identifier in prev_ids:
                raise RuntimeError("Result %s was in previous batch" % identifier)
            yield identifier
        offset += batch
        remain -= batch
        prev_ids = ids


def search(db, query, offset=None, limit=None, format=None):
    """Call TogoWS search.

    This is a low level wrapper for the TogoWS search function, which
    can return results in a several formats. In general, the search_iter
    function is more suitable for end users.

    Arguments:
     - db - database (string), see http://togows.dbcls.jp/search/
     - query - search term (string)
     - offset, limit - optional integers specifying which result to start from
       (1 based) and the number of results to return.
     - format - return data file format (string), e.g. "json", "ttl" (RDF)
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
        # TODO - Make this a ValueError? Right now despite the HTML website
        # claiming to, the "gene" or "ncbi-gene" don't work and are not listed.
        import warnings

        warnings.warn(
            "TogoWS search does not explicitly support database '%s'. "
            "See %s/search/ for options." % (db, _BASE_URL)
        )
    url = _BASE_URL + "/search/%s/%s" % (db, quote(query))
    if offset is not None and limit is not None:
        try:
            offset = int(offset)
        except ValueError:
            raise ValueError(
                "Offset should be an integer (at least one), not %r" % offset
            ) from None
        try:
            limit = int(limit)
        except ValueError:
            raise ValueError(
                "Limit should be an integer (at least one), not %r" % limit
            ) from None
        if offset <= 0:
            raise ValueError("Offset should be at least one, not %i" % offset)
        if limit <= 0:
            raise ValueError("Count should be at least one, not %i" % limit)
        url += "/%i,%i" % (offset, limit)
    elif offset is not None or limit is not None:
        raise ValueError("Expect BOTH offset AND limit to be provided (or neither)")
    if format:
        url += "." + format
    # print(url)
    return _open(url)


def convert(data, in_format, out_format):
    """Call TogoWS for file format conversion.

    Arguments:
     - data - string or handle containing input record(s)
     - in_format - string describing the input file format (e.g. "genbank")
     - out_format - string describing the requested output format (e.g. "fasta")

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
    # TODO - Should we just accept a string not a handle? What about a filename?
    try:
        # Handle
        data = data.read()
    except AttributeError:
        # String
        pass
    return _open(url, post=data)


def _open(url, post=None):
    """Build the URL and open a handle to it (PRIVATE).

    Open a handle to TogoWS, will raise an IOError if it encounters an error.

    In the absence of clear guidelines, this function enforces a limit of
    "up to three queries per second" to avoid abusing the TogoWS servers.
    """
    delay = 0.333333333  # one third of a second
    current = time.time()
    wait = _open.previous + delay - current
    if wait > 0:
        time.sleep(wait)
        _open.previous = current + wait
    else:
        _open.previous = current

    if post:
        handle = urlopen(url, post.encode())
    else:
        handle = urlopen(url)

    # We now trust TogoWS to have set an HTTP error code, that
    # suffices for my current unit tests. Previously we would
    # examine the start of the data returned back.
    text_handle = io.TextIOWrapper(handle, encoding="UTF-8")
    text_handle.url = handle.url
    return text_handle


_open.previous = 0


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest(verbose=0)
