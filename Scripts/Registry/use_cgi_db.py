"""Example of fetching sequence information from a CGI database.
"""
# get any swissprot database
from Bio import db
sp = db["swissprot"]

rec = sp['O23731']
# rec = sp['P15697']
if hasattr(rec, "read"): # it is a handle
    out = open("test.sp", "w")
    out.write(rec.read())
    out.close()
    # print rec.read()[:200]
else: # it is a real SeqRecord
    print rec
    print rec.description
    print rec.id
    print rec.name
    print rec.seq.data[0:20]

