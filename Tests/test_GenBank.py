#!/usr/bin/env python
"""Test the GenBank parser and make sure everything is working smoothly.
"""
# standard library
import os
import copy
import cStringIO

# Biopython
from Bio import File

# GenBank stuff to test
from Bio import GenBank
from Bio.GenBank import utils

gb_file_dir = os.path.join(os.getcwd(), 'GenBank')

#test_files = ['cor6_6.gb', 'arab1.gb', 'iro.gb', 'arab2.gb', 'pri1.gb',
#              'arab3.gb', 'arab4.gb', 'pri2.gb', 'bct1.gb', 'bct2.gb']

test_files = ['noref.gb', 'cor6_6.gb', 'iro.gb', 'pri1.gb', 'arab1.gb',
              'protein_refseq.gb', 'extra_keywords.gb']
test_files += ['one_of.gb', 'NT_019265.gb', 'origin_line.gb']

write_format_files = test_files[:]
# don't test writing on protein_refseq, since it is horribly nasty
write_format_files.remove("protein_refseq.gb")
# don't test writing on the CONTIG refseq, because the wrapping of
# locations won't work exactly
write_format_files.remove("NT_019265.gb")

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
                ann_keys = cur_record.annotations.keys()
                ann_keys.sort()
                for ann_key in ann_keys:
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
                    for qualifier in feature.qualifiers:
                        print "key:", qualifier.key, "value:", qualifier.value
                         
        handle.close()
        
# test GenBank dictionary
print "Testing dictionaries..."
dict_file = os.path.join(gb_file_dir, 'cor6_6.gb')
index_file = os.path.join(gb_file_dir, 'cor6_6.idx')

print "Indexing file to serve as a dictionary..."
GenBank.index_file(dict_file, index_file)
gb_dict = GenBank.Dictionary(index_file, GenBank.FeatureParser())

print "len:", len(gb_dict)
k = gb_dict.keys()
k.sort()
print "keys:", k

# pick out some keys and make sure we are getting back decent records
for key in k[:3]:
    print "Retrieving record with key %s" % key
    cur_seqrecord = gb_dict[key]
    print "description:", cur_seqrecord.description
    print "id:", cur_seqrecord.id

# test writing GenBank format
print "Testing writing GenBank format..."

def do_comparison(good_record, test_record):
    """Compare two records to see if they are the same.

    Ths compares the two GenBank record, and will raise an AssertionError
    if two lines do not match, showing the non-matching lines.
    """
    good_handle = cStringIO.StringIO(good_record)
    test_handle = cStringIO.StringIO(test_record)

    while 1:
        good_line = good_handle.readline()
        test_line = test_handle.readline()

        if not(good_line) and not(test_line):
            break

        if not(good_line):
            raise AssertionError("Extra info in Test: `%s`" % test_line)
        if not(test_line):
            raise AssertionError("Extra info in Expected: `%s`" % good_line)

        assert test_line == good_line, \
               "Expected does not match Test.\nExpect:`%s`\nTest  :`%s`\n" % \
               (good_line, test_line)
    
def t_write_format():
    record_parser = GenBank.RecordParser(debug_level = 0)

    for file in write_format_files:
        print "Testing GenBank writing for %s..." % os.path.basename(file)
        cur_handle = open(os.path.join("GenBank", file), "r")
        compare_handle = open(os.path.join("GenBank", file), "r")
        
        iterator = GenBank.Iterator(cur_handle, record_parser)
        compare_iterator = GenBank.Iterator(compare_handle)
        
        while 1:
            cur_record = iterator.next()
            compare_record = compare_iterator.next()
            
            if cur_record is None or compare_record is None:
                break

            print "\tTesting for %s" % cur_record.version

            output_record = str(cur_record) + "\n"
            do_comparison(compare_record, output_record)

        cur_handle.close()

t_write_format()    

def t_cleaning_features():
    """Test the ability to clean up feature values.
    """
    parser = GenBank.FeatureParser(feature_cleaner = \
                                   utils.FeatureValueCleaner())
    handle = open(os.path.join("GenBank", "arab1.gb"))
    iterator = GenBank.Iterator(handle, parser)

    first_record = iterator.next()
    
    # test for cleaning of translation
    translation_feature = first_record.features[1]
    test_trans = translation_feature.qualifiers["translation"][0]
    assert test_trans.find(" ") == -1, \
      "Did not clean spaces out of the translation"
    assert test_trans.find("\012") == -1, \
      "Did not clean newlines out of the translation"

print "Testing feature cleaning..."
t_cleaning_features()
