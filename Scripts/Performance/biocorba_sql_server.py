#!/usr/bin/env python
"""SQL-backed BioCorba server for testing purposes.

Usage:
    python biocorba_sql_server.py <file to write ior to>
"""
import sys

from BioSQL import BioSeqDatabase
from BioCorba.Server.Seqcore.CorbaCollection import BioSequenceCollection

def main(ior_file):
    db = get_database() 
    db_dict = BioSQLAdapterDictionary(db)
    corba_server = BioSequenceCollection(db_dict, "embl_rod", 1.0, "")
    corba_server.string_ior_to_file(ior_file)
    
    print "Server up and running with IOR at %s" % ior_file
    corba_server.run()

def get_database():
    """Perform a connection with the database.
    
    XXX The info here shouldn't be hard coded and should be specified
    on the commandline.
    """
    server = BioSeqDatabase.open_database(host = "192.168.0.192",
                user = "root", passwd = "", db = "biosql_new")
    return server["embl_rod"]

class BioSQLAdapterDictionary:
    """Simple dictionary-like class to plug BioSQL into BioCorba.
    """
    def __init__(self, sql_db):
        """Initialize with a connected BioSeqDatabase.
        """
        self.sql_db = sql_db

    def keys(self):
        """Just return all of the primary ids in the database.
        """
        return self.sql_db.keys()

    def __getitem__(self, key):
        """Fetch items by either primary id or accession information.
        """
        # primary id keys -- assumed to be either ints or longs
        if type(key) == type(1) or type(key) == type(1L):
            return self.sql_db[key]
        # strings are the actual accession or locus-type names
        elif type(key) == type(""):
            for fetch_type in ["display_id", "accession"]:
                try: # do the actual lookups
                    return apply(self.sql_db.lookup, (), {fetch_type : key})
                except IndexError:
                    pass
            # if we got here we didn't fetch it
            raise KeyError("Did not find record for key %s" % key)
        else:
            raise ValueError("Not sure how to handle key: %s" % key)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print "Incorrect arguments -- need to specify IOR file"
        print __doc__
        sys.exit()
    else:
        sys.exit(main(sys.argv[1]))
