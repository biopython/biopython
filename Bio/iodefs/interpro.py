# Copyright 2002 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio import register_io
from Bio.sources import EBI

from _support import *


register_io(
    name="interpro-ebi-cgi",
    source=EBI.IEntry,
    key="ac",
    failure=[(has_str("No InterPro entry"), "No InterPro entry")]
    )
