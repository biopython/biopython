# Copyright 2000 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
This module provides code to access resources at ExPASy over the WWW.
http://www.expasy.ch/


Functions:
get_prodoc_entry  Interface to the get-prodoc-entry CGI script.
get_prosite_entry Interface to the get-prosite-entry CGI script.
get_prosite_raw   Interface to the get-prosite-raw CGI script.
get_sprot_raw     Interface to the get-sprot-raw CGI script.
sprot_search_ful  Interface to the sprot-search-ful CGI script.
sprot_search_de   Interface to the sprot-search-de CGI script.

The function scanprosite1 is OBSOLETE; please see the
Bio.ExPASy.ScanProsite module for this functionality.
"""

import urllib


def get_prodoc_entry(id, cgi='http://www.expasy.ch/cgi-bin/get-prodoc-entry'):
    """get_prodoc_entry(id,
    cgi='http://www.expasy.ch/cgi-bin/get-prodoc-entry') -> handle

    Get a handle to a PRODOC entry at ExPASy in HTML format. 

    For a non-existing key XXX, ExPASy returns an HTML-formatted page
    containing this line:
    'There is no PROSITE documentation entry XXX. Please try again.'
    """
    # Open a handle to ExPASy.
    handle = urllib.urlopen("%s?%s" % (cgi, id))
    return handle

def get_prosite_entry(id,
                      cgi='http://www.expasy.ch/cgi-bin/get-prosite-entry'):
    """get_prosite_entry(id,
    cgi='http://www.expasy.ch/cgi-bin/get-prosite-entry') -> handle

    Get a handle to a PROSITE entry at ExPASy in HTML format.

    For a non-existing key XXX, ExPASy returns an HTML-formatted page
    containing this line:
    'There is currently no PROSITE entry for XXX. Please try again.'
    """
    handle = urllib.urlopen("%s?%s" % (cgi, id))
    return handle

def get_prosite_raw(id, cgi='http://www.expasy.ch/cgi-bin/get-prosite-raw.pl'):
    """get_prosite_raw(id,
                       cgi='http://www.expasy.ch/cgi-bin/get-prosite-raw.pl')
    -> handle

    Get a handle to a raw PROSITE or PRODOC entry at ExPASy.

    For a non-existing key, ExPASy returns nothing.
    """
    handle = urllib.urlopen("%s?%s" % (cgi, id))
    return handle

def get_sprot_raw(id, cgi=None):
    """Get a handle to a raw SwissProt entry at ExPASy.

    For an ID of XXX, fetches http://www.uniprot.org/uniprot/XXX.txt
    (as per the http://www.expasy.ch/expasy_urls.html documentation).


    For a non-existing key XXX, ExPASy returns an HTML Error 404 page.

    This function used to take a cgi option to specify the URL, but that
    is no longer supported. This is because prior to November 2009 we
    used to use http://www.expasy.ch/cgi-bin/get-sprot-raw.pl?XXX
    However, at the time of writting this returns FASTA format instead
    (probably an ExPASy/UniProt oversight). Under the new URL scheme,
    we cannot just append "?XXX" to the cgi argument.
    """
    if cgi :
        import warnings
        warnings.warn("The cgi argument in get_sprot_raw is not "
                      "supported anymore", DeprecationWarning)
    return urllib.urlopen("http://www.uniprot.org/uniprot/%s.txt" % id)

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
    options = urllib.urlencode(variables)
    fullcgi = "%s?%s" % (cgi, options)
    handle = urllib.urlopen(fullcgi)
    return handle

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
    options = urllib.urlencode(variables)
    fullcgi = "%s?%s" % (cgi, options)
    handle = urllib.urlopen(fullcgi)
    return handle

def scanprosite1(seq=None, id=None, exclude_frequent=None,
                 cgi='http://www.expasy.org/cgi-bin/scanprosite/scanprosite?1'):
    """scanprosite1(seq=None, id=None, exclude_frequent=None,
    cgi='http://www.expasy.org/cgi-bin/scanprosite/scanprosite?1') -> handle

    Scan a sequence for a Prosite pattern.  Either a sequence or a SwissProt/
    trEMBL sequence can be passed.  exclude_frequent specifies whether to
    exclude patterns with high probability.

    """
    import warnings
    warnings.warn("Bio.ExPASy.scanprosite1() has been deprecated, and we" \
                  +" intend to remove it in a future release of Biopython."\
                  +" Please use the Bio.ExPASy.ScanProsite module instead,"\
                  +" as described in the Tutorial.",
                  DeprecationWarning)
    variables = {}
    if seq:
        variables['SEQ'] = seq
    if id:
        variables['ID'] = id
    if exclude_frequent:
        variables['box'] = 'ok'
    options = urllib.urlencode(variables)
    handle = urllib.urlopen(cgi, options)
    return handle
