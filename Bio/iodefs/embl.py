# Copyright 2002 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio import register_io, group_io
from Bio.sources import EBI

from _support import *

register_io(
    name="embl-xembl-cgi",
    source=EBI.XEMBL,
    params=[("format", "Bsml")],
    key="id",
    failure=[(has_str("NOT EXIST"), "id does not exist")],
    )

register_io(
    name="embl-dbfetch-cgi",
    source=EBI.dbfetch,
    params=[("db", "embl"),
            ("style", "raw"),
            ("format", "embl"),
            ],
    key="id",
    failure=[(has_str("not found in database"), "id does not exist")]
    )

register_io(
    name="embl-ebi-cgi",
    source=EBI.emblfetch,
    params=[("db", "EMBL"),
            ("format", "default"),   # also Fasta, bsml, agave available
            ("style", "raw")
            ],
    key="id",
    failure=[(blank_expr, "No results returned")]
    )

register_io(
    name="embl",
    behavior="serial",
##    cache="XXX"
    )
group_io("embl", "embl-dbfetch-cgi")
group_io("embl", "embl-xembl-cgi")
group_io("embl", "embl-ebi-cgi")


register_io(
    name="embl-fast",
    behavior="concurrent",
    )
group_io("embl-fast", "embl-dbfetch-cgi")
group_io("embl-fast", "embl-xembl-cgi")
group_io("embl-fast", "embl-ebi-cgi")
