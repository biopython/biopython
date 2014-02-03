#!/usr/bin/env python
"""Find GenBank records that the parser has problems with within a big file.

This is meant to make it easy to get accession numbers for records that
don't parse properly.

Usage:
find_parser_problems.py <GenBank file to parse>
"""
# standard library
from __future__ import print_function

import sys

# GenBank
from Bio import GenBank

verbose = 0

if len(sys.argv) != 2:
    print("Usage ./find_parser_problems <GenBank file to parse>")
    sys.exit()

feature_parser = GenBank.FeatureParser(debug_level=0)
parser = GenBank.ErrorParser(feature_parser)

handle = open(sys.argv[1], 'r')
iterator = GenBank.Iterator(handle, parser, has_header=1)

while True:
    have_record = 0

    while have_record == 0:
        try:
            cur_record = next(iterator)
            have_record = 1
        except GenBank.ParserFailureError as msg:
            print("Parsing Problem: %s" % msg)
            sys.exit()

    if cur_record is None:
        break

    print("Successfully parsed record %s" % cur_record.id)

    if verbose:
        print("***Record")
        print("Seq: %s" % cur_record.seq)
        print("Id: %s" % cur_record.id)
        print("Name: %s" % cur_record.name)
        print("Description: %s" % cur_record.description)
        print("Annotations: %s" % cur_record.annotations)
        print("Feaures")
        for feature in cur_record.features:
            print(feature)

handle.close()
