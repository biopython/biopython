# Copyright 2002 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""ExPASy web page: http://www.expasy.ch/

"""
from Bio.sources import CGI

get_prodoc_entry = CGI(
    name="get-prodoc-entry",
    delay=5.0,
    cgi='http://www.expasy.ch/cgi-bin/get-prodoc-entry',
    doc="Retrieve a prodoc entry by ID",
    )

get_prosite_entry = CGI(
    name="get-prosite-entry",
    cgi='http://www.expasy.ch/cgi-bin/get-prosite-entry',
    delay=5.0,
    doc="Retrieve a prosite entry by ID",
    )

get_sprot_raw = CGI(
    name="get-sprot-raw",
    cgi='http://www.expasy.ch/cgi-bin/get-sprot-raw.pl',
    delay=5.0,
    doc="Retrieve a swiss-prot entry by ID",
    )

sprot_search_ful = CGI(
    name="sprot-search-ful",
    cgi='http://www.expasy.ch/cgi-bin/sprot-search-ful',
    delay=5.0,
    doc="Search swiss-prot",
    )
# 'makeWild' : 'on'  if make_wild
# 'S' : 'on'    FOR SWISSPROT
# 'T' : 'on'    FOR TREMBL
    
sprot_search_de = CGI(
    name="sprot-search-de",
    cgi='http://www.expasy.ch/cgi-bin/sprot-search-de',
    delay=5.0,
    doc="Search swissprot by name, description, gene name, species, or organelle.",
    )
##    variables = {'SEARCH' : text}
##    if swissprot:
##        variables['S'] = 'on'
##    if trembl:
##        variables['T'] = 'on'
##    return _open(cgi, variables)

scanprosite1 = CGI(
    name="scanprosite1",
    cgi='http://expasy.cbr.nrc.ca/cgi-bin/scanprosite/scanprosite?1',
    delay=5.0,
    doc="Scan a sequence for a Prosite pattern."
    )

##    cgi='http://expasy.cbr.nrc.ca/cgi-bin/scanprosite/scanprosite?1') -> handle
    
##    Scan a sequence for a Prosite pattern.  Either a sequence or a SwissProt/
##    trEMBL sequence can be passed.  exclude_frequent specifies whether to
##    exclude patterns with high probability.
    
##    """
##    variables = {}
##    if seq:
##        variables['SEQ'] = seq
##    if id:
##        variables['ID'] = id
##    if exclude_frequent:
##        variables['box'] = 'ok'
##    return _open(cgi, variables, get=0)
