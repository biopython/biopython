# Copyright 2000 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Code to access resources at ExPASy over the WWW.

See http://www.expasy.ch/


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

__docformat__ = "restructuredtext en"


def get_prodoc_entry(id, cgi='http://www.expasy.ch/cgi-bin/get-prodoc-entry'):
    """get_prodoc_entry(id,
    cgi='http://www.expasy.ch/cgi-bin/get-prodoc-entry') -> handle

    Get a handle to a PRODOC entry at ExPASy in HTML format.

    For a non-existing key XXX, ExPASy returns an HTML-formatted page
    containing this line:
    'There is no PROSITE documentation entry XXX. Please try again.'
    """
    # Open a handle to ExPASy.
    return _urlopen("%s?%s" % (cgi, id))


def get_prosite_entry(id,
                      cgi='http://www.expasy.ch/cgi-bin/get-prosite-entry'):
    """get_prosite_entry(id,
    cgi='http://www.expasy.ch/cgi-bin/get-prosite-entry') -> handle

    Get a handle to a PROSITE entry at ExPASy in HTML format.

    For a non-existing key XXX, ExPASy returns an HTML-formatted page
    containing this line:
    'There is currently no PROSITE entry for XXX. Please try again.'
    """
    return _urlopen("%s?%s" % (cgi, id))


def get_prosite_raw(id, cgi='http://www.expasy.ch/cgi-bin/get-prosite-raw.pl'):
    """get_prosite_raw(id,
                       cgi='http://www.expasy.ch/cgi-bin/get-prosite-raw.pl')
    -> handle

    Get a handle to a raw PROSITE or PRODOC entry at ExPASy.

    For a non-existing key, ExPASy returns nothing.
    """
    return _urlopen("%s?%s" % (cgi, id))


def get_sprot_raw(id):
    """Get a handle to a raw SwissProt entry at ExPASy.

    For an ID of XXX, fetches http://www.uniprot.org/uniprot/XXX.txt
    (as per the http://www.expasy.ch/expasy_urls.html documentation).
    """
    return _urlopen("http://www.uniprot.org/uniprot/%s.txt" % id)


def sprot_search_ful(text, make_wild=None, swissprot=1, trembl=None,
                     cgi='http://www.expasy.ch/cgi-bin/sprot-search-ful'):
    """sprot_search_ful(text, make_wild=None, swissprot=1, trembl=None,
    cgi='http://www.expasy.ch/cgi-bin/sprot-search-ful') -> handle

    Search SwissProt by full text.

    """
    variables = {'SEARCH': text}
    if make_wild:
        variables['makeWild'] = 'on'
    if swissprot:
        variables['S'] = 'on'
    if trembl:
        variables['T'] = 'on'
    options = _urlencode(variables)
    fullcgi = "%s?%s" % (cgi, options)
    handle = _urlopen(fullcgi)
    return handle


def sprot_search_de(text, swissprot=1, trembl=None,
                    cgi='http://www.expasy.ch/cgi-bin/sprot-search-de'):
    """sprot_search_de(text, swissprot=1, trembl=None,
    cgi='http://www.expasy.ch/cgi-bin/sprot-search-de') -> handle

    Search SwissProt by name, description, gene name, species, or
    organelle.

    """
    variables = {'SEARCH': text}
    if swissprot:
        variables['S'] = 'on'
    if trembl:
        variables['T'] = 'on'
    options = _urlencode(variables)
    fullcgi = "%s?%s" % (cgi, options)
    handle = _urlopen(fullcgi)
    return handle
