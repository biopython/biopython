# Copyright 2002 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio import register_db, group_db
from Bio.sources import ExPASy

from _support import *

register_db(
    name="swissprot-expasy-cgi",
    source=ExPASy.get_sprot_raw,
    key="",
    failure=[(blank_expr, "no results")]
    )

register_db(
    name="swissprot",
    behavior="serial",
    )
group_db("swissprot", "swissprot-expasy-cgi")
