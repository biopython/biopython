# Copyright 2002 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio import register_io
from Bio.sources import RCSB, EBI

from Martel import *
from _support import *

register_io(
    name="pdb-rcsb-cgi",
    source=RCSB.export,
    params=[("format", "PDB"),
            ("compression", "None")
            ],
    key="pdbId",
    failure=[(has_str("File not found"), "ID does not exist")]
    )
    

not_header_expr = AssertNot(Str("HEADER"))
register_io(
    name="pdb-ebi-cgi",
    source=EBI.emblfetch,
    params=[("db", "PDB"),
            ("format", "default"),   # also Fasta, bsml, agave available
            ("style", "raw"),
            ],
    key="id",
    failure=[(not_header_expr, "results do not look like PDB format")]
    )

