# Copyright 2002 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio.config.DBRegistry import CGIDB, DBGroup

prodoc_expasy_cgi = CGIDB(
    name="prodoc-expasy-cgi",
    doc="Retrieve a prodoc entry by ID",
    cgi='http://us.expasy.org/cgi-bin/get-prosite-raw.pl',
    delay=5.0,
    params=[],
    key="",
    )

prodoc = DBGroup(
        name = "prodoc",
        behavior = "serial"
    )
prodoc.add(prodoc_expasy_cgi)
