#!/usr/bin/env python
"""Test the GenBank parser and make sure everything is working smoothly.
"""
# standard library
import os
import copy

# Biopython
from Bio import File

# GenBank stuff to test
from Bio import GenBank

gb_file_dir = os.path.join(os.getcwd(), 'GenBank')

#test_files = ['cor6_6.gb', 'arab1.gb', 'iro.gb', 'arab2.gb', 'pri1.gb',
#              'arab3.gb', 'arab4.gb', 'pri2.gb', 'bct1.gb', 'bct2.gb']

test_files = ['noref.gb', 'cor6_6.gb', 'iro.gb', 'pri1.gb', 'arab1.gb']

files_to_parse = []
for file in test_files:
    files_to_parse.append(os.path.join(gb_file_dir, file))

# parse the bioperl test files
# comment this out for now -- there are a bunch of junky records in here
# that no longer exist in GenBank -- do we really need to support those?
# files_to_parse = [os.path.join(os.getcwd(), 'GenBank', 'bioperl_test.gb')]

# parse the biojava test files
# files_to_parse += [os.path.join(os.getcwd(), 'GenBank', 'biojava_test.gb')]

# test the parsers
feature_parser = GenBank.FeatureParser(debug_level = 0)
record_parser = GenBank.RecordParser(debug_level = 0)

all_parsers = [feature_parser, record_parser]
print "Testing parsers..."
for parser in all_parsers:
    for file in files_to_parse:
        handle = open(file, 'r')
        iterator = GenBank.Iterator(handle, parser)
        
        while 1:
            cur_record = iterator.next()

            if not(cur_record):
                break

            if isinstance(parser, GenBank.FeatureParser):
                print "***Record from the FeatureParser"
                print "Seq:", cur_record.seq
                print "Id:", cur_record.id
                print "Name:", cur_record.name
                print "Description", cur_record.description
                print "Annotations***"
                for ann_key in cur_record.annotations.keys():
                    if ann_key != 'references':
                        print "Key: %s" % ann_key
                        print "Value: %s" % \
                              cur_record.annotations[ann_key]
                    else:
                        print "References*"
                        for reference in cur_record.annotations[ann_key]:
                            print str(reference) 
                print "Feaures"
                for feature in cur_record.features:
                    print feature
            elif isinstance(parser, GenBank.RecordParser):
                print "***Record from the RecordParser"
                print "locus:", cur_record.locus
                print "definition:", cur_record.definition
                print "accession:", cur_record.accession
                for reference in cur_record.references:
                    print "reference title:", reference.title

                for feature in cur_record.features:
                    print "feature key:", feature.key
                    print "location:", feature.location
                    print "num qualifiers:", len(feature.qualifiers)

        handle.close()

# test GenBank dictionary
print "Testing dictionaries..."
dict_file = os.path.join(gb_file_dir, 'cor6_6.gb')
index_file = os.path.join(gb_file_dir, 'cor6_6.idx')

print "Indexing file to serve as a dictionary..."
GenBank.index_file(dict_file, index_file)
gb_dict = GenBank.Dictionary(index_file, GenBank.FeatureParser())

print "len:", len(gb_dict)
print "keys:", gb_dict.keys()

# pick out some keys and make sure we are getting back decent records
for key in gb_dict.keys()[:3]:
    print "Retrieving record with key %s" % key
    cur_seqrecord = gb_dict[key]
    print "description:", cur_seqrecord.description
    print "id:", cur_seqrecord.id


