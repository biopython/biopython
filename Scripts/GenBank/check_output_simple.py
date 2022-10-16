#!/usr/bin/env python
# Copyright 2000 Brad Chapman.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Display the SeqFeatures produced by the parser.

This produces a ton of output so it is possible to hand check what is
produced by the parser with the original GenBank file to make sure
everything is being parsed and output properly.

Usage:
./check_output.py <GenBank file to parse>

"""
# standard library

import sys

# GenBank stuff to test
from Bio import GenBank


if len(sys.argv) != 2:
    print("Usage ./check_output.py <GenBank file to parse>")
    sys.exit()

parser = GenBank.FeatureParser(debug_level=2)

with open(sys.argv[1]) as handle:

    iterator = GenBank.Iterator(handle, parser)

    while True:
        cur_record = next(iterator)

        if not cur_record:
            break

        print("***Record")
        print(f"Seq: {cur_record.seq}")
        print(f"Id: {cur_record.id}")
        print(f"Name: {cur_record.name}")
        print(f"Description: {cur_record.description}")
        print("Annotations****")
        for annotation_key in cur_record.annotations:
            if annotation_key != "references":
                print(f"Key: {annotation_key}")
                print(f"Value: {cur_record.annotations[annotation_key]}")
            else:
                print("References*")
                for reference in cur_record.annotations[annotation_key]:
                    print(reference)
        print("Features")
        for feature in cur_record.features:
            print(feature)
