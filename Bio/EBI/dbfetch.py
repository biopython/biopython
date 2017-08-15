# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Provides code to access the EBI dbfetch tool.

This module aims to make dbfetch easy to use. See:
http://www.ebi.ac.uk/Tools/dbfetch/

The dbfetch tool from EBI provvides access to different databases via a simple
URL that the module constructs. Errors are handled regarding databse access and
HTPP level.

The functionality is somewhat similar to Biopython's Bio.TogoWS, Bio.KEGG.REST
and Bio.Entrez modules.

EBI provvides extensive guidelines for the databases, their format and the
syntax necessary to use the tool. See:
http://www.ebi.ac.uk/Tools/dbfetch/syntax.jsp

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


def _q(arg1, arg2=None, arg3=None, arg4=None):
    URL = "http://www.ebi.ac.uk/Tools/dbfetch/dbfetch%s"
    if arg3 and arg4:
        args = "/%s/%s/%s?%s" % (arg1, arg2, arg3, arg4)
    elif arg3:
        args = "/%s/%s/%s" % (arg1, arg2, arg3)
    elif arg2:
        args = "/%s/%s" % (arg1, arg2)
    else:
        args = "/%s" % (arg1)
    resp = _urlopen(URL % (args))

    return _binary_to_string_handle(resp)


def db_info():
    """Display the current information for all databases in JSON format.

    A valid list of databases with and relative information is available
    here http://www.ebi.ac.uk/Tools/dbfetch/dbfetch/dbfetch.databases

    """
    return _q("dbfetch.databases?style=json")


def dbfetch(db, ids, format=None, style=None):
    """Create a URL to use the tool via RESTful method.

    db = database to query.
    ids = single id or list of ids separated by commas. Can be provvided as
          list or string object.
    format = output format. Depend on the database queried.
    style = raw by default, otherwise it can be specified as HTML.

    """
    if isinstance(ids, list):
        ids = ",".join(ids)
    elif isinstance(ids, int):
        ids = str(ids)
    else:
        ids = ids.replace(" ", "")
    return _q(db, ids, format, style)
