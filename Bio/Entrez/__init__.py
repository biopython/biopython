# Copyright 1999-2000 by Jeffrey Chang.  All rights reserved.
# Copyright 2008-2013 by Michiel de Hoon.  All rights reserved.
# Revisions copyright 2011-2016 by Peter Cock. All rights reserved.
# Revisions copyright 2015 by Eric Rasche. All rights reserved.
# Revisions copyright 2015 by Carlos Pena. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Provides code to access NCBI over the WWW.

The main Entrez web page is available at:
http://www.ncbi.nlm.nih.gov/Entrez/

Entrez Programming Utilities web page is available at:
http://www.ncbi.nlm.nih.gov/books/NBK25501/

This module provides a number of functions like ``efetch`` (short for
Entrez Fetch) which will return the data as a handle object. This is
a standard interface used in Python for reading data from a file, or
in this case a remote network connection, and provides methods like
``.read()`` or offers iteration over the contents line by line. See
also "What the heck is a handle?" in the Biopython Tutorial and
Cookbook: http://biopython.org/DIST/docs/tutorial/Tutorial.html
http://biopython.org/DIST/docs/tutorial/Tutorial.pdf

Unlike a handle to a file on disk from the ``open(filename)`` function,
which has a ``.name`` attribute giving the filename, the handles from
``Bio.Entrez`` all have a ``.url`` attribute instead giving the URL
used to connect to the NCBI Entrez API.

The Entrez module also provides an XML parser which takes a handle
as input.

Variables:

    - email        Set the Entrez email parameter (default is not set).
    - tool         Set the Entrez tool parameter (default is ``biopython``).

Functions:

    - efetch       Retrieves records in the requested format from a list of one or
      more primary IDs or from the user's environment
    - epost        Posts a file containing a list of primary IDs for future use in
      the user's environment to use with subsequent search strategies
    - esearch      Searches and retrieves primary IDs (for use in EFetch, ELink,
      and ESummary) and term translations and optionally retains
      results for future use in the user's environment.
    - elink        Checks for the existence of an external or Related Articles link
      from a list of one or more primary IDs.  Retrieves primary IDs
      and relevancy scores for links to Entrez databases or Related
      Articles;  creates a hyperlink to the primary LinkOut provider
      for a specific ID and database, or lists LinkOut URLs
      and Attributes for multiple IDs.
    - einfo        Provides field index term counts, last update, and available
      links for each database.
    - esummary     Retrieves document summaries from a list of primary IDs or from
      the user's environment.
    - egquery      Provides Entrez database counts in XML for a single search
      using Global Query.
    - espell       Retrieves spelling suggestions.
    - ecitmatch    Retrieves PubMed IDs (PMIDs) that correspond to a set of
      input citation strings.

    - read         Parses the XML results returned by any of the above functions.
      Typical usage is:

          >>> from Bio import Entrez
          >>> Entrez.email = "Your.Name.Here@example.org"
          >>> handle = Entrez.einfo() # or esearch, efetch, ...
          >>> record = Entrez.read(handle)
          >>> handle.close()

       where record is now a Python dictionary or list.

    - parse        Parses the XML results returned by those of the above functions
      which can return multiple records - such as efetch, esummary
      and elink. Typical usage is:

          >>> handle = Entrez.esummary(db="pubmed", id="19304878,14630660", retmode="xml")
          >>> records = Entrez.parse(handle)
          >>> for record in records:
          ...     # each record is a Python dictionary or list.
          ...     print(record['Title'])
          Biopython: freely available Python tools for computational molecular biology and bioinformatics.
          PDB file parser and structure class implemented in Python.
          >>> handle.close()

      This function is appropriate only if the XML file contains
      multiple records, and is particular useful for large files.

    - _open        Internally used function.

"""
from __future__ import print_function

import time
import warnings

# Importing these functions with leading underscore as not intended for reuse
from Bio._py3k import urlopen as _urlopen
from Bio._py3k import urlencode as _urlencode
from Bio._py3k import HTTPError as _HTTPError

from Bio._py3k import _binary_to_string_handle, _as_bytes


email = None
tool = "biopython"


# XXX retmode?
def epost(db, **keywds):
    """Post a file of identifiers for future use.

    Posts a file containing a list of UIs for future use in the user's
    environment to use with subsequent search strategies.

    See the online documentation for an explanation of the parameters:
    http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EPost

    Return a handle to the results.

    Raises an IOError exception if there's a network error.
    """
    cgi = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/epost.fcgi'
    variables = {'db': db}
    variables.update(keywds)
    return _open(cgi, variables, post=True)


def efetch(db, **keywords):
    """Fetches Entrez results which are returned as a handle.

    EFetch retrieves records in the requested format from a list of one or
    more UIs or from user's environment.

    See the online documentation for an explanation of the parameters:
    http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch

    Return a handle to the results.

    Raises an IOError exception if there's a network error.

    Short example:

    >>> from Bio import Entrez
    >>> Entrez.email = "Your.Name.Here@example.org"
    >>> handle = Entrez.efetch(db="nucleotide", id="AY851612", rettype="gb", retmode="text")
    >>> print(handle.readline().strip())
    LOCUS       AY851612                 892 bp    DNA     linear   PLN 10-APR-2007
    >>> handle.close()

    This will automatically use an HTTP POST rather than HTTP GET if there
    are over 200 identifiers as recommended by the NCBI.

    **Warning:** The NCBI changed the default retmode in Feb 2012, so many
    databases which previously returned text output now give XML.
    """
    cgi = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    variables = {'db': db}
    variables.update(keywords)
    post = False
    try:
        ids = variables["id"]
    except KeyError:
        pass
    else:
        if isinstance(ids, list):
            ids = ",".join(ids)
            variables["id"] = ids
        elif isinstance(ids, int):
            ids = str(ids)
            variables["id"] = ids

        if ids.count(",") >= 200:
            # NCBI prefers an HTTP POST instead of an HTTP GET if there are
            # more than about 200 IDs
            post = True
    return _open(cgi, variables, post=post)


def esearch(db, term, **keywds):
    """ESearch runs an Entrez search and returns a handle to the results.

    ESearch searches and retrieves primary IDs (for use in EFetch, ELink
    and ESummary) and term translations, and optionally retains results
    for future use in the user's environment.

    See the online documentation for an explanation of the parameters:
    http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch

    Return a handle to the results which are always in XML format.

    Raises an IOError exception if there's a network error.

    Short example:

    >>> from Bio import Entrez
    >>> Entrez.email = "Your.Name.Here@example.org"
    >>> handle = Entrez.esearch(db="nucleotide", retmax=10, term="opuntia[ORGN] accD", idtype="acc")
    >>> record = Entrez.read(handle)
    >>> handle.close()
    >>> int(record["Count"]) >= 2
    True
    >>> "EF590893.1" in record["IdList"]
    True
    >>> "EF590892.1" in record["IdList"]
    True

    """
    cgi = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
    variables = {'db': db,
                 'term': term}
    variables.update(keywds)
    return _open(cgi, variables)


def elink(**keywds):
    """ELink checks for linked external articles and returns a handle.

    ELink checks for the existence of an external or Related Articles link
    from a list of one or more primary IDs;  retrieves IDs and relevancy
    scores for links to Entrez databases or Related Articles; creates a
    hyperlink to the primary LinkOut provider for a specific ID and
    database, or lists LinkOut URLs and attributes for multiple IDs.

    See the online documentation for an explanation of the parameters:
    http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ELink

    Return a handle to the results, by default in XML format.

    Raises an IOError exception if there's a network error.

    This example finds articles related to the Biopython application
    note's entry in the PubMed database:

    >>> from Bio import Entrez
    >>> Entrez.email = "Your.Name.Here@example.org"
    >>> pmid = "19304878"
    >>> handle = Entrez.elink(dbfrom="pubmed", id=pmid, linkname="pubmed_pubmed")
    >>> record = Entrez.read(handle)
    >>> handle.close()
    >>> print(record[0]["LinkSetDb"][0]["LinkName"])
    pubmed_pubmed
    >>> linked = [link["Id"] for link in record[0]["LinkSetDb"][0]["Link"]]
    >>> "17121776" in linked
    True

    This is explained in much more detail in the Biopython Tutorial.
    """
    cgi = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi'
    variables = {}
    variables.update(keywds)
    return _open(cgi, variables)


def einfo(**keywds):
    """EInfo returns a summary of the Entez databases as a results handle.

    EInfo provides field names, index term counts, last update, and
    available links for each Entrez database.

    See the online documentation for an explanation of the parameters:
    http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EInfo

    Return a handle to the results, by default in XML format.

    Raises an IOError exception if there's a network error.

    Short example:

    >>> from Bio import Entrez
    >>> Entrez.email = "Your.Name.Here@example.org"
    >>> record = Entrez.read(Entrez.einfo())
    >>> 'pubmed' in record['DbList']
    True

    """
    cgi = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi'
    variables = {}
    variables.update(keywds)
    return _open(cgi, variables)


def esummary(**keywds):
    """ESummary retrieves document summaries as a results handle.

    ESummary retrieves document summaries from a list of primary IDs or
    from the user's environment.

    See the online documentation for an explanation of the parameters:
    http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESummary

    Return a handle to the results, by default in XML format.

    Raises an IOError exception if there's a network error.

    This example discovers more about entry 19923 in the structure
    database:

    >>> from Bio import Entrez
    >>> Entrez.email = "Your.Name.Here@example.org"
    >>> handle = Entrez.esummary(db="structure", id="19923")
    >>> record = Entrez.read(handle)
    >>> handle.close()
    >>> print(record[0]["Id"])
    19923
    >>> print(record[0]["PdbDescr"])
    Crystal Structure Of E. Coli Aconitase B

    """
    cgi = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi'
    variables = {}
    variables.update(keywds)
    return _open(cgi, variables)


def egquery(**keywds):
    """EGQuery provides Entrez database counts for a global search.

    EGQuery provides Entrez database counts in XML for a single search
    using Global Query.

    See the online documentation for an explanation of the parameters:
    http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EGQuery

    Return a handle to the results in XML format.

    Raises an IOError exception if there's a network error.

    This quick example based on a longer version from the Biopython
    Tutorial just checks there are over 60 matches for 'Biopython'
    in PubMedCentral:

    >>> from Bio import Entrez
    >>> Entrez.email = "Your.Name.Here@example.org"
    >>> handle = Entrez.egquery(term="biopython")
    >>> record = Entrez.read(handle)
    >>> handle.close()
    >>> for row in record["eGQueryResult"]:
    ...     if "pmc" in row["DbName"]:
    ...         print(int(row["Count"]) > 60)
    True

    """
    cgi = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/egquery.fcgi'
    variables = {}
    variables.update(keywds)
    return _open(cgi, variables)


def espell(**keywds):
    """ESpell retrieves spelling suggestions, returned in a results handle.

    ESpell retrieves spelling suggestions, if available.

    See the online documentation for an explanation of the parameters:
    http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESpell

    Return a handle to the results, by default in XML format.

    Raises an IOError exception if there's a network error.

    Short example:

    >>> from Bio import Entrez
    >>> Entrez.email = "Your.Name.Here@example.org"
    >>> record = Entrez.read(Entrez.espell(term="biopythooon"))
    >>> print(record["Query"])
    biopythooon
    >>> print(record["CorrectedQuery"])
    biopython

    """
    cgi = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/espell.fcgi'
    variables = {}
    variables.update(keywds)
    return _open(cgi, variables)


def _update_ecitmatch_variables(keywds):
    # XML is the only supported value, and it actually returns TXT.
    variables = {'retmode': 'xml'}
    citation_keys = ('journal_title', 'year', 'volume', 'first_page', 'author_name', 'key')

    # Accept pre-formatted strings
    if isinstance(keywds['bdata'], str):
        variables.update(keywds)
    else:
        # Alternatively accept a nicer interface
        variables['db'] = keywds['db']
        bdata = []
        for citation in keywds['bdata']:
            formatted_citation = '|'.join([citation.get(key, "") for key in citation_keys])
            bdata.append(formatted_citation)
        variables['bdata'] = '\r'.join(bdata)
    return variables


def ecitmatch(**keywds):
    """ECitMatch retrieves PMIDs-Citation linking

    ECitMatch retrieves PubMed IDs (PMIDs) that correspond to a set of input citation strings.

    See the online documentation for an explanation of the parameters:
    http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ECitMatch

    Return a handle to the results, by default in plain text

    Raises an IOError exception if there's a network error.

    Short example:

    >>> from Bio import Entrez
    >>> Entrez.email = "Your.Name.Here@example.org"
    >>> citation_1 = {"journal_title": "proc natl acad sci u s a",
    ...               "year": "1991", "volume": "88", "first_page": "3248",
    ...               "author_name": "mann bj", "key": "citation_1"}
    >>> handle = Entrez.ecitmatch(db="pubmed", bdata=[citation_1])
    >>> print(handle.read().strip().split("|"))
    ['proc natl acad sci u s a', '1991', '88', '3248', 'mann bj', 'citation_1', '2014248']
    >>> handle.close()

    """
    cgi = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ecitmatch.cgi'
    variables = _update_ecitmatch_variables(keywds)
    return _open(cgi, variables, ecitmatch=True)


def read(handle, validate=True):
    """Parses an XML file from the NCBI Entrez Utilities into python objects.

    This function parses an XML file created by NCBI's Entrez Utilities,
    returning a multilevel data structure of Python lists and dictionaries.
    Most XML files returned by NCBI's Entrez Utilities can be parsed by
    this function, provided its DTD is available. Biopython includes the
    DTDs for most commonly used Entrez Utilities.

    If validate is True (default), the parser will validate the XML file
    against the DTD, and raise an error if the XML file contains tags that
    are not represented in the DTD. If validate is False, the parser will
    simply skip such tags.

    Whereas the data structure seems to consist of generic Python lists,
    dictionaries, strings, and so on, each of these is actually a class
    derived from the base type. This allows us to store the attributes
    (if any) of each element in a dictionary my_element.attributes, and
    the tag name in my_element.tag.
    """
    from .Parser import DataHandler
    handler = DataHandler(validate)
    record = handler.read(handle)
    return record


def parse(handle, validate=True):
    """Parses an XML file from the NCBI Entrez Utilities into python objects.

    This function parses an XML file created by NCBI's Entrez Utilities,
    returning a multilevel data structure of Python lists and dictionaries.
    This function is suitable for XML files that (in Python) can be represented
    as a list of individual records. Whereas 'read' reads the complete file
    and returns a single Python list, 'parse' is a generator function that
    returns the records one by one. This function is therefore particularly
    useful for parsing large files.

    Most XML files returned by NCBI's Entrez Utilities can be parsed by
    this function, provided its DTD is available. Biopython includes the
    DTDs for most commonly used Entrez Utilities.

    If validate is True (default), the parser will validate the XML file
    against the DTD, and raise an error if the XML file contains tags that
    are not represented in the DTD. If validate is False, the parser will
    simply skip such tags.

    Whereas the data structure seems to consist of generic Python lists,
    dictionaries, strings, and so on, each of these is actually a class
    derived from the base type. This allows us to store the attributes
    (if any) of each element in a dictionary my_element.attributes, and
    the tag name in my_element.tag.
    """
    from .Parser import DataHandler
    handler = DataHandler(validate)
    records = handler.parse(handle)
    return records


def _open(cgi, params=None, post=None, ecitmatch=False):
    """Helper function to build the URL and open a handle to it (PRIVATE).

    Open a handle to Entrez.  cgi is the URL for the cgi script to access.
    params is a dictionary with the options to pass to it.  Does some
    simple error checking, and will raise an IOError if it encounters one.

    The arugment post should be a boolean to explicitly control if an HTTP
    POST should be used rather an HTTP GET based on the query length.
    By default (post=None), POST is used if the URL encoded paramters would
    be over 1000 characters long.

    This function also enforces the "up to three queries per second rule"
    to avoid abusing the NCBI servers.
    """
    # NCBI requirement: At most three queries per second.
    # Equivalently, at least a third of second between queries
    delay = 0.333333334
    current = time.time()
    wait = _open.previous + delay - current
    if wait > 0:
        time.sleep(wait)
        _open.previous = current + wait
    else:
        _open.previous = current

    params = _construct_params(params)
    options = _encode_options(ecitmatch, params)

    # By default, post is None. Set to a boolean to over-ride length choice:
    if post is None and len(options) > 1000:
        post = True
    cgi = _construct_cgi(cgi, post, options)

    try:
        if post:
            handle = _urlopen(cgi, data=_as_bytes(options))
        else:
            handle = _urlopen(cgi)
    except _HTTPError as exception:
        raise exception

    return _binary_to_string_handle(handle)


_open.previous = 0


def _construct_params(params):
    if params is None:
        params = {}

    # Remove None values from the parameters
    for key, value in list(params.items()):
        if value is None:
            del params[key]
    # Tell Entrez that we are using Biopython (or whatever the user has
    # specified explicitly in the parameters or by changing the default)
    if "tool" not in params:
        params["tool"] = tool
    # Tell Entrez who we are
    if "email" not in params:
        if email is not None:
            params["email"] = email
        else:
            warnings.warn("""
Email address is not specified.

To make use of NCBI's E-utilities, NCBI requires you to specify your
email address with each request.  As an example, if your email address
is A.N.Other@example.com, you can specify it as follows:
   from Bio import Entrez
   Entrez.email = 'A.N.Other@example.com'
In case of excessive usage of the E-utilities, NCBI will attempt to contact
a user at the email address provided before blocking access to the
E-utilities.""", UserWarning)
    return params


def _encode_options(ecitmatch, params):
    # Open a handle to Entrez.
    options = _urlencode(params, doseq=True)
    # _urlencode encodes pipes, which NCBI expects in ECitMatch
    if ecitmatch:
        options = options.replace('%7C', '|')
    return options


def _construct_cgi(cgi, post, options):
    if not post:
        # HTTP GET
        cgi += "?" + options
    return cgi


if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()
