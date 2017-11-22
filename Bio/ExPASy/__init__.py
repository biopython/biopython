# Copyright 2000 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Code to access resources at ExPASy over the WWW.

See https://www.expasy.org/


Functions:
 - get_prodoc_entry  Interface to the get-prodoc-entry CGI script.
 - get_prosite_entry Interface to the get-prosite-entry CGI script.
 - get_prosite_raw   Interface to the get-prosite-raw CGI script.
 - get_sprot_raw     Interface to the get-sprot-raw CGI script.
 - sprot_search_ful  Interface to the sprot-search-ful CGI script.
 - sprot_search_de   Interface to the sprot-search-de CGI script.

"""

# Importing these functions with leading underscore as not intended for reuse
from Bio._py3k import urlopen as _urlopen
from Bio._py3k import urlencode as _urlencode
from Bio._py3k import _binary_to_string_handle


def get_prodoc_entry(id,
                     cgi='https://prosite.expasy.org/cgi-bin/prosite/get-prodoc-entry'):
    """Get a text handle to a PRODOC entry at ExPASy in HTML format.

    >>> from Bio import ExPASy
    >>> with ExPASy.get_prodoc_entry('PDOC00001') as in_handle:
    ...     html = in_handle.read()
    ...
    >>> with open("myprodocrecord.html", "w") as out_handle:
    ...     out_handle.write(html)
    ...

    For a non-existing key XXX, ExPASy returns an HTML-formatted page
    containing this text: 'There is currently no PROSITE entry for'
    """
    return _binary_to_string_handle(_urlopen("%s?%s" % (cgi, id)))


def get_prosite_entry(id,
                      cgi='https://prosite.expasy.org/cgi-bin/prosite/get-prosite-entry'):
    """Get a text handle to a PROSITE entry at ExPASy in HTML format.

    >>> from Bio import ExPASy
    >>> with ExPASy.get_prosite_entry('PS00001') as in_handle:
    ...     html = in_handle.read()
    ...
    >>> with open("myprositerecord.html", "w") as out_handle:
    ...     out_handle.write(html)
    ...

    For a non-existing key XXX, ExPASy returns an HTML-formatted page
    containing this text: 'There is currently no PROSITE entry for'
    """
    return _binary_to_string_handle(_urlopen("%s?%s" % (cgi, id)))


def get_prosite_raw(id, cgi=None):
    """Get a text handle to a raw PROSITE or PRODOC record at ExPASy.

    The cgi argument is deprecated due to changes in the ExPASy
    website.

    For a non-existing key, ExPASy returns nothing.

    >>> from Bio import ExPASy
    >>> from Bio.ExPASy import Prosite
    >>> with ExPASy.get_prosite_raw('PS00001') as handle:
    ...    record = Prosite.read(handle)
    ...
    >>> print(record.accession)
    PS00001


    For a non-existing key, ExPASy returns an error:

    >>> handle = get_prosite_raw("does_not_exist")
    Traceback (most recent call last):
    ...
    urllib.error.HTTPError: HTTP Error 404: Not Found

    """
    url = "https://prosite.expasy.org/%s.txt" % id
    return _binary_to_string_handle(_urlopen(url))


def get_sprot_raw(id):
    """Get a text handle to a raw SwissProt entry at ExPASy.

    For an ID of XXX, fetches http://www.uniprot.org/uniprot/XXX.txt
    (as per the https://www.expasy.org/expasy_urls.html documentation).

    >>> from Bio import ExPASy
    >>> from Bio import SwissProt
    >>> with ExPASy.get_sprot_raw("O23729") as handle:
    ...     record = SwissProt.read(handle)
    ...
    >>> print(record.entry_name)
    CHS3_BROFI

    For a non-existing identifier, UniProt returns an error:

    >>> ExPASy.get_sprot_raw("DOES_NOT_EXIST")
    Traceback (most recent call last):
    ...
    urllib.error.HTTPError: HTTP Error 404: Not Found

    """
    url = "http://www.uniprot.org/uniprot/%s.txt" % id
    return _binary_to_string_handle(_urlopen(url))


def sprot_search_ful(text, make_wild=None, swissprot=1, trembl=None,
                     cgi='https://www.expasy.org/cgi-bin/sprot-search-ful'):
    """Search SwissProt by full text (BROKEN)."""
    variables = {'SEARCH': text}
    if make_wild:
        variables['makeWild'] = 'on'
    if swissprot:
        variables['S'] = 'on'
    if trembl:
        variables['T'] = 'on'
    options = _urlencode(variables)
    fullcgi = "%s?%s" % (cgi, options)
    handle = _binary_to_string_handle(_urlopen(fullcgi))
    return handle


def sprot_search_de(text, swissprot=1, trembl=None,
                    cgi='https://www.expasy.org/cgi-bin/sprot-search-de'):
    """Search SwissProt (BROKEN).

    Search by name, description, gene name, species, or organelle.
    """
    variables = {'SEARCH': text}
    if swissprot:
        variables['S'] = 'on'
    if trembl:
        variables['T'] = 'on'
    options = _urlencode(variables)
    fullcgi = "%s?%s" % (cgi, options)
    return _binary_to_string_handle(_urlopen(fullcgi))
