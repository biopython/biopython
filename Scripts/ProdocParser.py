#!/usr/bin/env python

# Copyright 2000 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Sample script for accessing the Prodoc Parser.

"""

import sys
import getopt

def scan(file, format):
    from Bio import ParserSupport
    consumer = ParserSupport.TaggingConsumer()
    
    if format == 'text':
        from Bio.Prosite import Prodoc
        scanner = Prodoc._Scanner()
    else:
        sys.stderr.write("I don't understand the format %s\n" % format)
        sys.exit(-1)
    scanner.feed(open(file), consumer)

def parse(file, format):
    from Bio.Prosite import Prodoc
    if format == 'text':
        parser = Prodoc.RecordParser()
    else:
        sys.stderr.write("I don't understand the format %s\n" % format)
        sys.exit(-1)

    iter = Prodoc.Iterator(open(file), parser=parser)
    while 1:
        rec = iter.next()
        if not rec:
            break
        print "Successfully parsed %s" % rec.accession

def print_instructions(handle=sys.stdout):
    instructions = """
ProdocParser.py [OPTIONS] <filename>

XXX describe script here

OPTIONS
    -f FORMAT
        The format of the input file.  FORMAT can be only be
        text (default).

    -s
        Only scan the file, printing out the scanned results

"""
    handle.write(instructions)


if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print_instructions()
        sys.exit(0)

    options, files = getopt.getopt(sys.argv[1:], 'f:s')
    
    if len(files) != 1:
        print_instructions()
        sys.exit(0)
    file = files[0]

    format = 'text'
    scan_only = 0
    
    for opt, value in options:
        if opt == '-s':
            scan_only = 1
        elif opt == '-f':
            format = value

    if scan_only:
        scan(file, format)
    else:
        parse(file, format)
