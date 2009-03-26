#!/usr/bin/env python
"""Test for the SProt parser on SwissProt files.
"""
import os
from Bio import SeqIO
from Bio import SwissProt
from Bio.SeqRecord import SeqRecord
from Bio.SwissProt import SProt

test_files = ['sp001', 'sp002', 'sp003', 'sp004', 'sp005',
              'sp006', 'sp007', 'sp008', 'sp009', 'sp010',
              'sp011', 'sp012', 'sp013', 'sp014', 'sp015']

# test the record parser
for test_file in test_files:
    print "\ntesting %s..." % test_file
    datafile = os.path.join('SwissProt', test_file)

    print "*Using SequenceParser"
    test_handle = open(datafile)
    seq_record = SeqIO.read(test_handle, "swiss")
    test_handle.close()

    assert isinstance(seq_record, SeqRecord)

    print seq_record.id
    print seq_record.name
    print seq_record.description
    print repr(seq_record.seq)

    print "*Using RecordParser"
    test_handle = open(datafile)
    record = SProt.read(test_handle)
    test_handle.close()

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

    #Check the two parsers agree on the essentials
    assert seq_record.seq.tostring() == record.sequence    
    assert seq_record.description == record.description
    assert seq_record.name == record.entry_name
    assert seq_record.id in record.accessions

    #Now try using the Iterator - note that all these
    #test cases have only one record.

    # With the SequenceParser
    test_handle = open(datafile)
    records = list(SeqIO.parse(test_handle, "swiss"))
    test_handle.close()

    assert len(records) == 1
    assert isinstance(records[0], SeqRecord)

    #Check matches what we got earlier without the iterator:
    assert records[0].seq.tostring() == seq_record.seq.tostring()
    assert records[0].description == seq_record.description
    assert records[0].name == seq_record.name
    assert records[0].id == seq_record.id
    
    # With the RecordParser
    test_handle = open(datafile)
    records = list(SProt.parse(test_handle))
    test_handle.close()

    assert len(records) == 1
    assert isinstance(records[0], SProt.Record)
    
    #Check matches what we got earlier without the iterator:
    assert records[0].sequence == record.sequence
    assert records[0].description == record.description
    assert records[0].entry_name == record.entry_name
    assert records[0].accessions == record.accessions
    
