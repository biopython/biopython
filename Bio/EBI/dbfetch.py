# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Provides code to access the EBI dbfetch tool.

This module aims to make dbfetch easy to use. See:
http://www.ebi.ac.uk/Tools/dbfetch/

The dbfetch tool from EBI provides access to different databases via a simple
URL that the module constructs. Errors are handled regarding databse access and
HTTP level.

The functionality is somewhat similar to Biopython's Bio.TogoWS, Bio.KEGG.REST
and Bio.Entrez modules.

EBI provides extensive guidelines for the databases, their format and the
syntax necessary to use the tool. See:
http://www.ebi.ac.uk/Tools/dbfetch/syntax.jsp

This module does not support the SOAP service for dbfetch. See:
http://www.ebi.ac.uk/Tools/webservices/services/dbfetch

References:
Lopez R, Cowley A, Li W, McWilliam H.;
(2013) Using EMBL-EBI Services via Web Interface and Programmatically via
Web Services. Curr Protoc Bioinformatics Volume 48 (2014) p.3.12.1-3.12.50.

McWilliam H., Li W., Uludag M., Squizzato S., Park Y.M., Buso N., Cowley A.P.,
Lopez R.;(2013) Analysis Tool Web Services from the EMBL-EBI Nucleic Acids
Research 41: W597-W600.

"""

from __future__ import print_function

from Bio._py3k import urlopen as _urlopen
from Bio._py3k import _binary_to_string_handle

import requires_internet
requires_internet.check()


_BASE_URL = "http://www.ebi.ac.uk/Tools/dbfetch/dbfetch%s"


def _query_by_url(url):
    """Use a URL to perform the request to the service (PRIVATE)."""
    resp = _urlopen(_BASE_URL % url)
    return _binary_to_string_handle(resp)


def db_info():
    """Display the current information for all databases in JSON format.

    A valid list of databases with and relative information is available
    here http://www.ebi.ac.uk/Tools/dbfetch/dbfetch/dbfetch.databases

    """
    return _query_by_url("dbfetch.databases?style=json")


def dbfetch(db, ids, format=None, style=None):
    """Create a URL to use the tool via RESTful method.

    - db    - database to query.
    - ids   - single id or list of ids separated by commas. Can be provided as
              list or string object.
    - format - output format. Depend on the database queried.
    - style - raw by default, otherwise it can be specified as HTML.

    """
    if isinstance(ids, list):
        ids = ",".join(ids)
    elif isinstance(ids, int):
        ids = str(ids)
    else:
        ids = ids.replace(" ", "")
    return _query_by_db(db, ids, format, style)


def _query_by_db(db, ids=None, format=None, style=None):
    """Prepare the URL for the database request to the service (PRIVATE)."""
    url = '/%s' % db
    if ids:
        url += '/%s' % ids
    if format:
        url += '/%s' % format
    if style:
        url += '?%s' % style
    return _query_by_url(url)
