# Copyright 2002 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio.config.DBRegistry import CGIDB, DBGroup
from _support import *

swissprot_expasy_cgi = CGIDB(
    name="swissprot-expasy-cgi",
    doc="Retrieve a swiss-prot entry by ID from ExPASy.",
    cgi="http://www.expasy.ch/cgi-bin/get-sprot-raw.pl",
    delay=5.0,
    timeout=10,
    params=[],
    key="",
    failure_cases=[(blank_expr, "no results")]
    )

swissprot_usmirror_cgi = CGIDB(
    name="swissprot-usmirror-cgi",
    doc="Retrieve a swiss-prot entry by ID from the US mirror.",
    cgi="http://us.expasy.org/cgi-bin/get-sprot-raw.pl",
    delay=5.0,
    timeout=10,
    params=[],
    key="",
    failure_cases=[(blank_expr, "no results")]
    )

swissprot = DBGroup(
    name="swissprot",
    behavior="serial",
    doc="Retrieve a swiss-prot entry by ID.  Will try different servers until one works.",
    )
swissprot.add(swissprot_expasy_cgi)
swissprot.add(swissprot_usmirror_cgi)
