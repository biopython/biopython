# Copyright 2002 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio.config.DBRegistry import CGIDB, DBGroup
from _support import *

prosite_expasy_cgi = CGIDB(
    name="prosite-expasy-cgi",
    doc="Retrieve a prosite entry by ID",
    cgi='http://www.expasy.ch/cgi-bin/get-prosite-entry',
    delay=5.0,
    params=[],
    key="",
    failure_cases=[(has_str("There is currently no PROSITE entry"),
              "No PROSITE entry")],
    )
