# Copyright 2002 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio.config.DBRegistry import CGIDB, DBGroup
from _support import *

from Martel import *

not_header_expr = AssertNot(Str("HEADER"))
pdb_rcsb_cgi = CGIDB(
    name="pdb-rcsb-cgi",
    cgi="http://www.rcsb.org/pdb/cgi/export.cgi",
    url="XXX PLEASE FILL THIS IN XXX",
    delay=5.0,
    params=[("format", "PDB"),
            ("compression", "None")
            ],
    key="pdbId",
    # failure cases for file not found are making retrieval freeze up 
    # while Martel checks for them, for some reason I can't figure
    # so we go with checking to make sure results look like PDB
    # failure_cases=[(has_str("File not found"), "ID does not exist")],
    failure_cases=[(not_header_expr, "results do not look like PDB format")]
    )

pdb_ebi_cgi = CGIDB(
    name="pdb-ebi-cgi",
    cgi="http://www.ebi.ac.uk/cgi-bin/emblfetch",
    url="http://www.ebi.ac.uk/cgi-bin/emblfetch",
    delay=5.0,
    params=[("db", "PDB"),
            ("format", "default"),   # also Fasta, bsml, agave available
            ("style", "raw"),
            ],
    key="id",
    failure_cases=[(not_header_expr, "results do not look like PDB format")]
    )

pdb = DBGroup(
    name="pdb",
    behavior="serial"
    )
pdb.add(pdb_rcsb_cgi)
pdb.add(pdb_ebi_cgi)
