#!/usr/bin/env python
# Copyright 2002 Brad Chapman.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Test timing of getting records from a BioSQL database."""
from __future__ import print_function

import time
# set up the connection
from BioSQL import BioSeqDatabase


server = BioSeqDatabase.open_database(host="192.168.0.192", user="root",
                                      passwd="", db="test_biosql")
db = server["embl_rod"]

# -- do the fasta-only timing part
start_time = time.time()
num_records = 0
for junk_id, record in db.items():
    num_records += 1
    sequence = record.seq.data
    d = record.description
    i = record.id
    n = record.name

end_time = time.time()
elapsed_time = end_time - start_time
print("Fasta")
print("\tDid %s records in %s seconds for\n\t%f records per second" %
      (num_records, elapsed_time, float(num_records) / float(elapsed_time)))

# -- do the "EMBL" timing part
start_time = time.time()
num_records = 0
for junk_id, record in db.items():
    num_records += 1
    sequence = record.seq.data
    d = record.description
    i = record.id
    n = record.name
    features = record.features
    anns = record.annotations
    dates = record.dates
    species = record.species
    keywords = record.keywords
end_time = time.time()
elapsed_time = end_time - start_time
print("EMBL")
print("\tDid %s records in %s seconds for\n\t%f records per second" %
      (num_records, elapsed_time, float(num_records) / float(elapsed_time)))
