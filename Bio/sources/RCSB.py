# Copyright 2002 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""http://www.rcsb.org/

"""
from Bio.sources import CGI
from Martel import *


export = CGI(
    name="export",
    delay=5.0,
    cgi="http://www.rcsb.org/pdb/cgi/export.cgi",
    )

