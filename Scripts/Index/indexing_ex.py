#!/usr/bin/env python
"""Example code for indexing a file using the new Flat and Berkeley indexes.

Usage:
    python indexing_ex.py --type=<index_type> --dbname=<name> 
        <files> <to> <index>

Where:
    --type : Specifies the index type to use, either "berkeley" or "flat"
             Defaults to "flat"
    --dbname : The name of the database to create. Defaults to "test"
"""
import sys
import getopt

from Bio.Mindy import SimpleSeqRecord

def main(index_type, files, db_name):
    if index_type == "berkeley":
        create_index = SimpleSeqRecord.create_berkeleydb
    elif index_type == "flat":
        create_index = SimpleSeqRecord.create_flatdb
    else:
        raise ValueError("Bad index type: %s" % index_type)

    create_index(files, db_name)

if __name__ == "__main__":
    try:
        options, args = getopt.getopt(sys.argv[1:], '', 
                                      ["type=", "dbname="])
    except getopt.error, msg:
        print msg
        print __doc__
        sys.exit()
    
    index_type = "flat"
    db_name = "test"
    for opt, value in options:
        if opt == "--type":
            index_type = value
        elif opt == "--dbname":
            db_name = value

    sys.exit(main(index_type, args, db_name))
