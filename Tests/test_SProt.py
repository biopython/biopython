#!/usr/bin/env python
"""Test for the SProt parser on SwissProt files.
"""
import os
from Bio.SwissProt import SProt

test_files = ['sp001', 'sp002', 'sp003', 'sp004', 'sp005',
              'sp006', 'sp007', 'sp008', 'sp009', 'sp010',
              'sp011', 'sp012', 'sp013', 'sp014']

record_parser = SProt.RecordParser()
sequence_parser = SProt.SequenceParser()

# test the record parser
for test_file in test_files:
    print "\ntesting %s..." % test_file
    datafile = os.path.join('SwissProt', test_file)
    test_handle = open(datafile)

    record = record_parser.parse(test_handle)

    # test a couple of things on the record -- this is not exhaustive
    print record.entry_name
    print record.accessions
    print record.organism_classification
    print record.seqinfo
    
    print "***Features:"
    for feature in record.features:
        print feature

    print "***References:"
    for ref in record.references:
        print "authors:", ref.authors
        print "title:", ref.title
        print "references:", ref.references

    test_handle.close()

# test the sequence parser
for test_file in test_files:
    print "\ntesting %s..." % test_file
    datafile = os.path.join('SwissProt', test_file)
    test_handle = open(datafile)

    seq_record = sequence_parser.parse(test_handle)
    print seq_record.id
    print seq_record.name
    print seq_record.description
    print seq_record.seq

    test_handle.close()
    

    
