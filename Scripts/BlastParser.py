#!/usr/bin/env python

# Copyright 1999-2000 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Sample script for accessing the Blast Parser.

"""

import sys
import getopt

def scan(file, format):
    from Bio import ParserSupport
    consumer = ParserSupport.TaggingConsumer()
    
    if format == 'ncbi-blast' or format == 'ncbi-psiblast':
        from Bio.Blast import NCBIStandalone
        scanner = NCBIStandalone._Scanner()
    elif format == 'ncbi-www':
        from Bio.Blast import NCBIWWW
        scanner = NCBIWWW._Scanner()
    else:
        sys.stderr.write("I don't understand the format %s\n" % format)
        sys.exit(-1)
    scanner.feed(open(file), consumer)

def parse(file, format):
    if format == 'ncbi-blast':
        from Bio.Blast import NCBIStandalone
        parser = NCBIStandalone.BlastParser()
    elif format == 'ncbi-psiblast':
        from Bio.Blast import NCBIStandalone
        parser = NCBIStandalone.PSIBlastParser()
    elif format == 'ncbi-www':
        from Bio.Blast import NCBIWWW
        sys.stderr.write("The parser for NCBIWWW is not written.  Sorry.\n")
        sys.exit(-1)
    else:
        sys.stderr.write("I don't understand the format %s\n" % format)
        sys.exit(-1)
    obj = parser.parse(open(file))
    print "Success!"

def print_instructions(handle=sys.stdout):
    instructions = """
BlastParser.py [OPTIONS] <filename>

XXX describe script here

OPTIONS
    -f FORMAT
        The format of the input file.  FORMAT can be either
        ncbi-blast, ncbi-psiblast, or ncbi-www (default).

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

    format = 'ncbi-www'
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
