"""Example of registering a BioSQL database using Biopython-specific system.
"""
# register our BioSQL database

from Bio.sources import BioSQL
fred_db = BioSQL(name = "fred",
                 doc = "Test MySQL server on fred",
                 db_host = "192.168.0.192",
                 db_user = "root",
                 db_passwd = "",
                 sql_db = "test_biosql",
                 namespace_db = "embl_rod")

from Bio import register_db
register_db(name = "biosql-embl-fred",
            key = "accession",
            source = fred_db,
            failure = [])

# deal with the database
from Bio import db
biosql_db = db["biosql-embl-fred"]


rec = biosql_db["AB030760"]
print rec
print rec.description
print rec.id
print rec.name
print rec.seq.data[0:20]

