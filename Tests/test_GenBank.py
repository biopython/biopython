#!/usr/bin/env python
# Copyright 2001-2004 by Brad Chapman.  All rights reserved.
# Revisions copyright 2007-2016 by Peter Cock. All rights reserved.
# Revisions copyright 2013 by Kai Blin. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Test the GenBank parser and make sure everything is working smoothly.
"""
# standard library
from __future__ import print_function

import os
from Bio._py3k import StringIO
import warnings

from Bio import BiopythonParserWarning

# GenBank stuff to test
from Bio import GenBank
from Bio.GenBank import utils

from Bio.Alphabet import _get_base_alphabet, ProteinAlphabet

gb_file_dir = os.path.join(os.getcwd(), "GenBank")

test_files = ["noref.gb", "cor6_6.gb", "iro.gb", "pri1.gb", "arab1.gb",
              "protein_refseq.gb", "extra_keywords.gb", "one_of.gb",
              "NT_019265.gb", "origin_line.gb", "blank_seq.gb",
              "dbsource_wrap.gb", "gbvrl1_start.seq", "NC_005816.gb",
              "no_end_marker.gb", "wrong_sequence_indent.gb",
              "invalid_locus_line_spacing.gb", "empty_feature_qualifier.gb",
              "invalid_misc_feature.gb", "1MRR_A.gp"]

# We only test writing on a subset of the examples:
write_format_files = ["noref.gb", "cor6_6.gb", "iro.gb", "pri1.gb", "arab1.gb",
                      "extra_keywords.gb", "one_of.gb", "origin_line.gb"]
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
# files_to_parse = [os.path.join(os.getcwd(), "GenBank", "bioperl_test.gb")]

# parse the biojava test files
# files_to_parse += [os.path.join(os.getcwd(), "GenBank", "biojava_test.gb")]

# test the parsers
feat_parser = GenBank.FeatureParser(debug_level=0)
rec_parser = GenBank.RecordParser(debug_level=0)

all_parsers = [feat_parser, rec_parser]
print("Testing parsers...")
for parser in all_parsers:
    for filename in files_to_parse:
        if not os.path.isfile(filename):
            print("Missing test input file: %s" % filename)
            continue

        handle = open(filename, "r")
        gb_iterator = GenBank.Iterator(handle, parser)

        while True:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", BiopythonParserWarning)
                # e.g. BiopythonParserWarning: Premature end of file in sequence data
                cur_record = next(gb_iterator)

            if cur_record is None:
                break

            if isinstance(parser, GenBank.FeatureParser):
                print("***Record from %s with the FeatureParser"
                      % filename.split(os.path.sep)[-1])
                print("Seq: %r" % cur_record.seq)
                print("Id: %s" % cur_record.id)
                print("Name: %s" % cur_record.name)
                print("Description: %s" % cur_record.description)
                print("Annotations***")
                ann_keys = sorted(cur_record.annotations)
                for ann_key in ann_keys:
                    if ann_key != "references":
                        print("Key: %s" % ann_key)
                        print("Value: %s" %
                              cur_record.annotations[ann_key])
                    else:
                        print("References*")
                        for reference in cur_record.annotations[ann_key]:
                            print(str(reference))
                print("Features")
                for feature in cur_record.features:
                    print(feature)
                    if isinstance(_get_base_alphabet(cur_record.seq.alphabet),
                                  ProteinAlphabet):
                        assert feature.strand is None
                    else:
                        # Assuming no mixed strand examples...
                        assert feature.strand is not None
                print("DB cross refs %s" % cur_record.dbxrefs)
            elif isinstance(parser, GenBank.RecordParser):
                print("***Record from %s with the RecordParser"
                      % filename.split(os.path.sep)[-1])
                print("sequence length: %i" % len(cur_record.sequence))
                print("locus: %s" % cur_record.locus)
                print("definition: %s" % cur_record.definition)
                print("accession: %s" % cur_record.accession)
                for reference in cur_record.references:
                    print("reference title: %s" % reference.title)

                for feature in cur_record.features:
                    print("feature key: %s" % feature.key)
                    print("location: %s" % feature.location)
                    print("num qualifiers: %i" % len(feature.qualifiers))
                    for qualifier in feature.qualifiers:
                        print("key: %s value: %s" % (qualifier.key, qualifier.value))

        handle.close()

# test writing GenBank format
print("Testing writing GenBank format...")


def do_comparison(good_record, test_record):
    """Compare two records to see if they are the same.

    Ths compares the two GenBank record, and will raise an AssertionError
    if two lines do not match, showing the non-matching lines.
    """
    good_handle = StringIO(good_record)
    test_handle = StringIO(test_record)

    while True:
        good_line = good_handle.readline()
        test_line = test_handle.readline()

        if not good_line and not test_line:
            break
        if not good_line:
            raise AssertionError("Extra info in Test: %r" % test_line)
        if not test_line:
            raise AssertionError("Extra info in Expected: %r" % good_line)
        test_normalized = " ".join(x for x in test_line.split() if x)
        good_normalized = " ".join(x for x in good_line.split() if x)
        assert test_normalized == good_normalized, \
            "Expected does not match Test.\nExpect: %r\nTest:   %r\n" % (good_line, test_line)


def t_write_format():
    record_parser = GenBank.RecordParser(debug_level=0)

    for next_file in write_format_files:
        print("Testing GenBank writing for %s..." % os.path.basename(next_file))
        cur_handle = open(os.path.join("GenBank", next_file), "r")
        compare_handle = open(os.path.join("GenBank", next_file), "r")

        iterator = GenBank.Iterator(cur_handle, record_parser)
        compare_iterator = GenBank.Iterator(compare_handle)

        while True:
            cur_rec = next(iterator)
            compare_record = next(compare_iterator)

            if cur_rec is None or compare_record is None:
                break

            print("\tTesting for %s" % cur_rec.version)

            output_record = str(cur_rec) + "\n"
            do_comparison(compare_record, output_record)

        cur_handle.close()
        compare_handle.close()


t_write_format()


def t_cleaning_features():
    """Test the ability to clean up feature values.
    """
    gb_parser = GenBank.FeatureParser(feature_cleaner=utils.FeatureValueCleaner())
    handle = open(os.path.join("GenBank", "arab1.gb"))
    iterator = GenBank.Iterator(handle, gb_parser)

    first_record = next(iterator)

    # test for cleaning of translation
    translation_feature = first_record.features[1]
    test_trans = translation_feature.qualifiers["translation"][0]
    assert " " not in test_trans, "Did not clean spaces out of the translation"
    assert "\012" not in test_trans, "Did not clean newlines out of the translation"

    handle.close()


print("Testing feature cleaning...")
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

    print("Done")


print("Testing EnsEMBL LOCUS lines...")
t_ensembl_locus()
