# Copyright 2002 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""XEMBL at http://www.ebi.ac.uk/xembl/

"""
from Bio.sources import CGI

XEMBL = CGI(
    name="XEMBL",
    delay=5.0,
    cgi="http://www.ebi.ac.uk/cgi-bin/xembl/XEMBL.pl",
    url="http://www.ebi.ac.uk/xembl/",
    doc="Query XEMBL for EMBL sequence data in XML format.",
    )
