# Copyright 2002 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""NCBI Entrez: http://www.ncbi.nlm.nih.gov/Entrez/

A list of the Entrez utilities is available at:
http://www.ncbi.nlm.nih.gov/entrez/utils/utils_index.html

"""
from Bio.sources import CGI
from Martel import *

from Bio.iodefs._support import *

proxy_error_expr = has_expr(Alt(Str("500"), Str("502")) + Str(" Proxy Error"))
diagnostic_error_expr = has_str("WWW Error 500 Diagnostic")
error_expr = Str("ERROR")

ncbi_failure_cases = [(proxy_error_expr, "proxy error"),
                      (diagnostic_error_expr, "diagnostic error"),
                      (error_expr, "ERROR")
                      ]

query = CGI(
    name="query",
    delay=5.0,
    cgi="http://www.ncbi.nlm.nih.gov/entrez/query.fcgi",
    url="http://www.ncbi.nlm.nih.gov/entrez/query/static/linking.html",
    doc="Query Entrez",
    failure_cases=ncbi_failure_cases,
    )

pmfetch = CGI(
    name="pmfetch",
    delay=5.0,
    cgi='http://www.ncbi.nlm.nih.gov/entrez/utils/pmfetch.fcgi',
    url="http://www.ncbi.nlm.nih.gov/entrez/utils/pmfetch_help.html",
    doc='Query PmFetch',
    failure_cases=ncbi_failure_cases,
    )

pmqty = CGI(
    name="pmqty",
    delay=5.0,
    cgi='http://www.ncbi.nlm.nih.gov/entrez/utils/pmqty.fcgi',
    url='http://www.ncbi.nlm.nih.gov/entrez/utils/pmqty_help.html',
    doc="Query PmQty",
    failure_cases=ncbi_failure_cases,
##            required_params=["db", "term"]
##            optional_params=[("dopt", None),
##                             ("dispmax", None),
##                             ("dispstart", None)
##                             ]
    )

# CHUNKED ACCESSOR
pmneighbor = CGI(
    name="pmneighbor",
    delay=5.0,
    cgi='http://www.ncbi.nlm.nih.gov/entrez/utils/pmneighbor.fcgi',
    url='http://www.ncbi.nlm.nih.gov/entrez/utils/pmneighbor_help.html',
    doc="Query PMNeighbor",
    failure_cases=ncbi_failure_cases,
##            required_params=["pmid", "display"],
    )
    # Warning: HUGE HACK HERE!  pmneighbor expects the display
    # parameter to be passed as just a tag, with no value.
    # Unfortunately, _open doesn't support these types of parameters,
    # so I'm building my own cgi string.  This is really due to the
    # limitations of urllib.urlencode.  We'll have to figure out a
    # good workaround.


