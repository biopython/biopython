# Copyright 2002 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio import register_db
from Bio.sources import ExPASy

from _support import *

register_db(
    name="prosite-expasy-cgi",
    source=ExPASy.get_prosite_entry,
    key="",
    failure=[(has_str("There is currently no PROSITE entry"),
              "No PROSITE entry")],
    )
