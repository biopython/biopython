#!/usr/bin/env python

# Copyright 2000 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import sys
import getopt

from Bio import PubMed

DEFAULT_DELAY = 5

def print_usage():
    print """query_pubmed.py [-h] [-c] [-d delay] query

This script sends a query to PubMed* and prints the MEDLARS-
formatted results to the screen.

Arguments:
    -h           Print out this help message.
    -c           Count the hits, and don't print them out.
    -d delay     Set the delay between queries.  This script limits
                 the query rate to prevent saturating NCBI's bandwidth.
                 By default one query is sent every %g seconds.

* http://www.ncbi.nlm.nih.gov/PubMed/
""" % DEFAULT_DELAY

if __name__ == '__main__':
    try:
        optlist, args = getopt.getopt(sys.argv[1:], "hcd:")
    except getopt.error, x:
        print x
        sys.exit(0)
    if len(args) != 1:     # If they gave extraneous arguments,
        print_usage()      # print the instructions and quit.
        sys.exit(0)
    query = args[0]

    help = 0
    delay = DEFAULT_DELAY
    count_only = 0
    for opt, arg in optlist:
        if opt == '-h':
            help = 1
        elif opt == '-c':
            count_only = 1
        elif opt == '-d':
            try:
                delay = float(arg)
            except ValueError:
                print "Delay must be a floating point value"
                sys.exit(0)
            if delay < 0:
                print "Delay cannot be negative"
                sys.exit(0)
    if help:
        print_usage()
        sys.exit(0)

    print "Doing a PubMed search for %s..." % repr(query)
    
    ids = PubMed.search_for(query)
    print "Found %d citations" % len(ids)

    if count_only:
        sys.exit(0)

    pm = PubMed.Dictionary(delay=delay)
    for id in ids:
        try:
            print pm[id]
        except KeyError, x:
            print "Couldn't download %s, %s" % (id, x)
        sys.stdout.flush()
