# Copyright 2002 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio.config.DBRegistry import CGIDB, DBGroup
from _support import *

# This returns HTML-formatted data.  Is there a way to get raw text?
interpro_ebi_cgi = CGIDB(
    name="interpro-ebi-cgi",
    cgi='http://www.ebi.ac.uk/interpro/IEntry',
    doc="Retrieve an InterPro entry",
    delay=5.0,
    key="ac",
    failure_cases=[(has_str("No InterPro entry"), "No InterPro entry")]
    )
