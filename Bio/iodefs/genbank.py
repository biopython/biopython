# Copyright 2002 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio import register_io
from Bio.sources import NCBI

from Martel import *

from _support import *

ncbi_failures=[
    (html_expr, "I got HTML and shouldn't have"),
    (Str("Please try again later"), "Please try again later"),
    (Str("The sequence has been intentionally withdrawn"),
     "Sequence withdrawn"),
    (blank_expr, "No data returned")
    ]

    
register_io(
    name="nucleotide-genbank-cgi",
    source=NCBI.query,
    params=[("cmd", "Text"),
            ("db", "Nucleotide"),
            ("dopt", "GenBank")
            ],
    key="uid",
    failure=ncbi_failures
    )


# If the id is not in the database, I get a message like:
# ERROR : GenPept does not exist for gi = 433174
not_exist_expr = Str("ERROR") + Re("[^d]*") + Str("does not exist for gi")

register_io(
    name="protein-genbank-cgi",
    source=NCBI.query,
    params=[("cmd", "Text"),
            ("db", "Protein"),
            ("dopt", "GenPept")
            ],
    key="uid",
    failure=ncbi_failures+[(not_exist_expr, "GI does not exist")]
    )

