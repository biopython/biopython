#!/usr/bin/env python
"""Check for the ability to read and write identical GenBank records.

This script takes as input a single file to be tested for reading records.
It will run through this file and check that the output of the parsed and
printed record matches the original record.

Usage:
python check_output.py <name of file to parse>
"""
# standard modules
import sys
import os
import cStringIO
import gzip

# biopython
from Bio import GenBank

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
            if good_line.strip():
                raise AssertionError("Extra info in Test: `%s`" % test_line)
        if not(test_line):
            if test_line.strip():
                raise AssertionError("Extra info in Expected: `%s`"
                                     % good_line)

        assert test_line == good_line, \
               "Expected does not match Test.\nExpect:`%s`\nTest  :`%s`\n" % \
               (good_line, test_line)
    
def write_format(file):
    record_parser = GenBank.RecordParser(debug_level = 2)

    print "Testing GenBank writing for %s..." % os.path.basename(file)
    # be able to handle gzipped files
    if file.find(".gz") >= 0:
        cur_handle = gzip.open(file, "r")
        compare_handle = gzip.open(file, "r")
    else:
        cur_handle = open(file, "r")
        compare_handle = open(file, "r")

    iterator = GenBank.Iterator(cur_handle, record_parser)
    compare_iterator = GenBank.Iterator(compare_handle)
        
    while 1:
        cur_record = iterator.next()
        compare_record = compare_iterator.next()
            
        if cur_record is None or compare_record is None:
            break

        # print "\tTesting for %s" % cur_record.version

        output_record = str(cur_record) + "\n"
        try:
            do_comparison(compare_record, output_record)
        except AssertionError, msg:
            print "\tTesting for %s" % cur_record.version
            print msg

    cur_handle.close()
    compare_handle.close()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print __doc__
        sys.exit()

    write_format(sys.argv[1])  
