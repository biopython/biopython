#!/usr/bin/env python

import sys
import getopt

import Bio

USAGE = """%s [-h] [dbname]

This script tests the databases in the db Registry.  If dbname is not
specified, then tests all the databases.


OPTIONS:

-h      Show this help file.

dbname  An optional name of the database to test.

""" % sys.argv[0]


# List of test cases as dbname, id.
tests = [
    ('prosite-expasy-cgi', 'PS50108'),
    ('prodoc-expasy-cgi', 'pdoc50108'),
    ('interpro-ebi-cgi', 'ipr000999'),
    ('embl-dbfetch-cgi', 'ab050095'),
    ('embl-ebi-cgi', 'ab050095'),
    ('embl-xembl-cgi', 'ab050095'),
    ('embl', 'ab050095'),
    ('embl-fast', 'ab050095'),
#    ('nucleotide-dbfetch-cgi', '19923388'),   # Disappeared?
    ('nucleotide-genbank-cgi', '19923388'),
    ('protein-genbank-cgi', '6323655'),
    ('swissprot', 'p50105'),
    ('swissprot-expasy-cgi', 'p50105'),
    ('pdb', '1hlb'),
    ('pdb-ebi-cgi', '1hlb'),
    ('pdb-rcsb-cgi', '1hlb'),
    ]

def test_db(dbname, id):
    linelen = (60 - len(dbname)+2)/2
    line = "=" * linelen
    print "%s %s %s" % (line, dbname, line)

    handle = Bio.db[dbname][id]
    print handle.read()[:500]
    

if __name__ == '__main__':
    try:
        optlist, args = getopt.getopt(sys.argv[1:], "h")
    except getopt.error, x:
        print >>sys.stderr, x
        sys.exit(-1)
    if len(args) > 1:
        print >>sys.stderr, USAGE
        sys.exit(-1)

    if args:
        DBNAME, = args
    else:
        DBNAME = None

    for opt, arg in optlist:
        if opt == '-h':
            print USAGE
            sys.exit(0)

    for dbname, id in tests:
        if DBNAME is not None and dbname != DBNAME:
            continue
        test_db(dbname, id)
