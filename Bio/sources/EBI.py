# Copyright 2002 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""EBI web page: http://www.ebi.ac.uk/

"""
from Bio.sources import CGI

# emblfetch can retrieve EBI sequences, Ensemble sequences, Medline
# records, PDB structures, and SWALL.
emblfetch = CGI(
    name="emblfetch",
    delay=5.0,
    cgi="http://www.ebi.ac.uk/cgi-bin/emblfetch",
    url="http://www.ebi.ac.uk/cgi-bin/emblfetch",
    doc="Retrieve many kinds of sequences from EBI"
    )

XEMBL = CGI(
    name="XEMBL",
    delay=5.0,
    cgi="http://www.ebi.ac.uk/cgi-bin/xembl/XEMBL.pl",
    url="http://www.ebi.ac.uk/xembl/",
    doc="Query XEMBL for EMBL sequence data in XML format.",
    )

IEntry = CGI(
    name="IEntry",
    delay=5.0,
    cgi='http://www.ebi.ac.uk/interpro/IEntry',
    doc="Retrieve an InterPro entry",
    )

dbfetch = CGI(
    name="dbfetch",
    delay=5.0,
    cgi="http://www.ebi.ac.uk/cgi-bin/dbfetch",
    doc="dbfetch provides EMBL, Genbank, and SWALL sequences",
    )
