"""Example of registering a BioCorba database using Biopython-specific system.
"""
# register our BioCorba database
from Bio.sources import BioCorba
fred_corba = BioCorba(name = "fred",
                      doc = "Test CORBA server on fred",
                      ior_ref = "http://fred/CORBA/py_collection.ior")

from Bio import register_db
register_db(name = "biocorba-test-fred",
            key = "accession",
            source = fred_corba,
            failure = [])

# deal with the database
from Bio import db
biocorba_db = db["biocorba-test-fred"]


rec = biocorba_db["X55053"]
print rec
print rec.description
print rec.id
print rec.name
print rec.seq.data[0:20]

