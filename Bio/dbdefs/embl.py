# Copyright 2002 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio.config.DBRegistry import CGIDB, DBGroup
from _support import *

embl_xembl_cgi = CGIDB(
    name="embl-xembl-cgi",
    doc="Query XEMBL for EMBL sequence data in XML format.",
    cgi="http://www.ebi.ac.uk/cgi-bin/xembl/XEMBL.pl",
    url="http://www.ebi.ac.uk/xembl/",
    delay=5.0,
    params=[("format", "Bsml")],
    key="id",
    failure_cases=[(has_str("NOT EXIST"), "id does not exist")],
    )

embl_dbfetch_cgi = CGIDB(
    name="embl-dbfetch-cgi",
    cgi="http://www.ebi.ac.uk/cgi-bin/dbfetch",
    url="http://www.ebi.ac.uk/cgi-bin/dbfetch",
    doc="dbfetch provides EMBL, Genbank, and SWALL sequences",
    delay=5.0,
    params=[("db", "embl"),
            ("style", "raw"),
            ("format", "embl"),
            ],
    key="id",
    failure_cases=[(has_str("not found in database"), "id does not exist")]
    )

embl_ebi_cgi = CGIDB(
    name="embl-ebi-cgi",
    cgi="http://www.ebi.ac.uk/cgi-bin/emblfetch",
    url="http://www.ebi.ac.uk/cgi-bin/emblfetch",
    doc="Retrieve many kinds of sequences from EBI",
    delay=5.0,
    params=[("db", "EMBL"),
            ("format", "default"),   # also Fasta, bsml, agave available
            ("style", "raw")
            ],
    key="id",
    failure_cases=[(blank_expr, "No results returned")]
    )

embl = DBGroup(
    name="embl",
    behavior="serial",
##    cache="XXX"
    )
embl.add(embl_dbfetch_cgi)
embl.add(embl_xembl_cgi)
embl.add(embl_ebi_cgi)


embl_fast = DBGroup(
    name="embl-fast",
    behavior="concurrent",
    )
embl_fast.add(embl_dbfetch_cgi)
embl_fast.add(embl_xembl_cgi)
embl_fast.add(embl_ebi_cgi)
