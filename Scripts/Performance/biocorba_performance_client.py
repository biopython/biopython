#/usr/bin/env python
"""Test timing of getting records from BioCorba BioSequenceCollection servers.

Usage:
    python biocorba_performance_client.py <url of ior> 
"""
import time
import sys

assert len(sys.argv) == 2, "Need to specify IOR location on the commandline"
ior_url = sys.argv[1]

# connect to the CORBA server
from BioCorba.Client.BiocorbaConnect import GenericCorbaClient
from BioCorba.Client.Seqcore.CorbaCollection import BioSequenceCollection
from BioCorba.Client.Bsane.Base import StopIteration

from Bio.SeqFeature import FeatureLocation

client_connector = GenericCorbaClient(BioSequenceCollection)
corba_db = client_connector.from_url_ior(ior_url)

test_ids = ["AB030644", "AB030700", "AB030729", "AB030734", "AB030737",
            "AB030738", "AB030755", "AB030756", "AB030757", "AB030758",
            "AB024565", "AB024566", "AB024567", "AB024573", "AB024614",
            "AB024688", "AB024689", "AB024713", "AB024717", "AB024719"]

# test_ids = ["HUMBDNF", "NT_010368", "HUMBETGLOA"]

# -- do the sequence-only timing part
start_time = time.time()

for test_id in test_ids:
    bio_seq = corba_db.resolve(test_id)
    sequence = bio_seq.seq()

end_time = time.time()
num_records = len(test_ids)
elapsed_time = end_time - start_time
print "Sequence"
print "\tDid %s records in %s seconds for\n\t%f records per second" % \
      (num_records, elapsed_time, float(num_records) / float(elapsed_time))

# -- do the "EMBL" timing part
start_time = time.time()

for test_id in test_ids:
    bio_seq = corba_db.resolve(test_id)
    
    # primary seq stuff
    sequence = bio_seq.seq()
    t = bio_seq.get_type()
    length = bio_seq.get_length()

    # features
    f_col = bio_seq.get_seq_features()
    nothing, f_it = f_col.get_features_on_region(0, FeatureLocation(1, length))
    while 1:
        try:
            f = f_it.next()
        except StopIteration:
            break
        
        # do stuff with the features
        f.get_start()
        f.get_end()
        f.get_locations()

        ann_col = f.get_annotations()
        nothing, ann_it = ann_col.get_annotations(0)
        while 1:
            try:
                ann = ann_it.next()
            except StopIteration:
                break
            ann.get_name()
            ann.get_value()

end_time = time.time()
num_records = len(test_ids)
elapsed_time = end_time - start_time
print "Full"
print "\tDid %s records in %s seconds for\n\t%f records per second" % \
      (num_records, elapsed_time, float(num_records) / float(elapsed_time))
