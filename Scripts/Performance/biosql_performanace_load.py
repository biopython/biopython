#/usr/bin/env python
"""Small script to test timing of loading records into a BioSQL database.
"""
import time
# set up the connection
from Bio import GenBank
from BioSQL import BioSeqDatabase
server = BioSeqDatabase.open_database (host = "192.168.0.192", user = "root", 
                                       passwd = "", db = "pythonloadtest")

# remove the database if it already exists
db_name = "testload"
try:
    server[db_name]
    server.remove_database(db_name)
except KeyError:
    pass
db = server.new_database(db_name)

input_file = "/home/hack/install/biopython/Tests/GenBank/cor6_6.gb"
handle = open(input_file, "r")
parser = GenBank.FeatureParser()
iterator = GenBank.Iterator(handle, parser)

# -- do the timing part
start_time = time.time()
num_records = db.load(iterator)
end_time = time.time()
elapsed_time = end_time - start_time
print "Loading"
print "\tDid %s records in %s seconds for\n\t%f records per second" % \
      (num_records, elapsed_time, float(num_records) / float(elapsed_time))

