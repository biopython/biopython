"""Example of registering an Indexed file using the Biopython system.
"""
# register our indexed file
from Bio.sources import IndexedFile
index_file = IndexedFile(name = "test",
                         doc = "Test indexed file",
                         dbname = "/home/chapmanb/bioppjx/tmp/test")

from Bio import register_db
register_db(name = "indexfile-test-swissprot",
            key = "id",
            source = index_file,
            failure = [])

# deal with the database
from Bio import db
indexed_db = db["indexfile-test-swissprot"]


rec = indexed_db["N33_HUMAN"]
print rec
print rec.description
print rec.id
print rec.name
print rec.seq.data[0:20]

