# Copyright 2002 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio import register_io
from Bio.sources import ExPASy

from Martel import *

register_io(
    name="prodoc-expasy-cgi",
    source=ExPASy.get_prodoc_entry,
    key="",
    )
