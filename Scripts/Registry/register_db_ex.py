"""Example of registering a database with the Biopython-specific system.
"""
from Bio.sources import CGI
local_cgi = CGI(name = "local_cgi",
                delay = 0.0,
                cgi = "http://www.myserver.org/cgi-bin/my_local.cgi",
                url = "http://www.myserver.org/cgi_documentation.html",
                doc = "Query a local databases",
                failure_cases = [])

import Martel
my_failures = [
     (Martel.Str("Sequence not available"), "No sequence found")]

from Bio import register_db
register_db(name = "nucleotide-genbank-local",
            key = "uid",
            source = local_cgi,
            failure = my_failures)

register_db(name = "genbank", behavior = "concurrent")
from Bio import group_db
group_db("genbank", "nucleotide-genbank-local")
group_db("genbank", "nucleotide-genbank-cgi")

from Bio import db
print db.keys()
