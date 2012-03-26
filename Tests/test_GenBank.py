#!/usr/bin/env python
"""Test the GenBank parser and make sure everything is working smoothly.
"""
# standard library
import os
import cStringIO

# GenBank stuff to test
from Bio import GenBank
from Bio.GenBank import utils

gb_file_dir = os.path.join(os.getcwd(), 'GenBank')

test_files = ['noref.gb', 'cor6_6.gb', 'iro.gb', 'pri1.gb', 'arab1.gb',
              'protein_refseq.gb', 'extra_keywords.gb', 'one_of.gb',
              'NT_019265.gb', 'origin_line.gb', 'blank_seq.gb',
              'dbsource_wrap.gb', 'gbvrl1_start.seq', 'NC_005816.gb']

# We only test writing on a subset of the examples:
write_format_files = ['noref.gb', 'cor6_6.gb', 'iro.gb', 'pri1.gb', 'arab1.gb',
                      'extra_keywords.gb', 'one_of.gb', 'origin_line.gb']
# don't test writing on protein_refseq, since it is horribly nasty
# don't test writing on the CONTIG refseq, because the wrapping of
# locations won't work exactly
# don't test writing on blank_seq because it lacks a sequence type
# don't test dbsource_wrap because it is a junky RefSeq file

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
    for filename in files_to_parse:
        if not os.path.isfile(filename):
            print "Missing test input file: %s" % filename
            continue
        
        handle = open(filename, 'r')
        iterator = GenBank.Iterator(handle, parser)
        
        while 1:
            cur_record = iterator.next()

            if cur_record is None:
                break

            if isinstance(parser, GenBank.FeatureParser):
                print "***Record from %s with the FeatureParser" \
                      % filename.split(os.path.sep)[-1]
                print "Seq:", repr(cur_record.seq)
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
                print "DB cross refs", cur_record.dbxrefs
            elif isinstance(parser, GenBank.RecordParser):
                print "***Record from %s with the RecordParser" \
                      % filename.split(os.path.sep)[-1]
                print "sequence length: %i" % len(cur_record.sequence)
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
        
#The dictionaries code has been deprecated
#print "Testing dictionaries..."
#...

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
        test_normalized = ' '.join([x for x in test_line.split() if x])
        good_normalized = ' '.join([x for x in good_line.split() if x])
        assert test_normalized == good_normalized, \
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
        compare_handle.close()

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

    handle.close()

print "Testing feature cleaning..."
t_cleaning_features()

def t_ensembl_locus():
    line = "LOCUS       HG531_PATCH 1000000 bp DNA HTG 18-JUN-2011\n"
    s = GenBank.Scanner.GenBankScanner()
    c = GenBank._FeatureConsumer(True)
    s._feed_first_line(c, line)
    assert c.data.name == "HG531_PATCH", c.data.name
    assert c._expected_size == 1000000, c._expected_size

    line = "LOCUS       HG531_PATCH 759984 bp DNA HTG 18-JUN-2011\n"
    s = GenBank.Scanner.GenBankScanner()
    c = GenBank._FeatureConsumer(True)
    s._feed_first_line(c, line)
    assert c.data.name == "HG531_PATCH", c.data.name
    assert c._expected_size == 759984, c._expected_size

    line = "LOCUS       HG506_HG1000_1_PATCH 814959 bp DNA HTG 18-JUN-2011\n"
    s = GenBank.Scanner.GenBankScanner()
    c = GenBank._FeatureConsumer(True)
    s._feed_first_line(c, line)
    assert c.data.name == "HG506_HG1000_1_PATCH", c.data.name
    assert c._expected_size == 814959, c._expected_size

    line = "LOCUS       HG506_HG1000_1_PATCH 1219964 bp DNA HTG 18-JUN-2011\n"
    s = GenBank.Scanner.GenBankScanner()
    c = GenBank._FeatureConsumer(True)
    s._feed_first_line(c, line)
    assert c.data.name == "HG506_HG1000_1_PATCH", c.data.name
    assert c._expected_size == 1219964, c._expected_size

    print "Done"
print "Testing EnsEMBL LOCUS lines..."
t_ensembl_locus()
