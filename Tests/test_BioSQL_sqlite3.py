#!/usr/bin/env python
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Run BioSQL tests using SQLite"""
import os
import tempfile
from Bio import MissingExternalDependencyError
from Bio import SeqIO
from BioSQL import BioSeqDatabase

from common_BioSQL import *

# Constants for the database driver
DBHOST = 'localhost'
DBUSER = 'root'
DBPASSWD = ''

DBDRIVER = 'sqlite3'
DBTYPE = 'sqlite'

# In memory SQLite does not work with current test structure since the tests
# expect databases to be retained between individual tests.
#TESTDB = ':memory:'
# Instead, we use (if we can) /dev/shm
try:
    test_db_fname = tempfile.mkstemp(dir='/dev/shm')[1]
except OSError:
    # We can't use /dev/shm
    h, test_db_fname = tempfile.mkstemp("_BioSQL.db")
    os.close(h)

TESTDB = test_db_fname


#This will abort if driver not installed etc:
check_config(DBDRIVER, DBTYPE, DBHOST, DBUSER, DBPASSWD, TESTDB)

#Some of the unit tests don't create their own database,
#so just in case there is no database already:
create_database()


if False:
    #This is how I generated test file Tests/BioSQL/cor6_6.db
    #which is test cross-checked with the latest bindings to
    #catch any regressions in how we map GenBank entries to
    #the database.
    assert not os.path.isfile("BioSQL/cor6_6.db")
    server = BioSeqDatabase.open_database(driver=DBDRIVER,
                                          db="BioSQL/cor6_6.db")
    DBSCHEMA = "biosqldb-" + DBTYPE + ".sql"
    SQL_FILE = os.path.join(os.getcwd(), "BioSQL", DBSCHEMA)
    assert os.path.isfile(SQL_FILE), SQL_FILE
    server.load_database_sql(SQL_FILE)
    server.commit()
    db = server.new_database("OLD")
    count = db.load(SeqIO.parse("GenBank/cor6_6.gb", "gb"))
    assert count == 6
    server.commit()
    assert len(db) == 6
    server.close()


class BackwardsCompatibilityTest(unittest.TestCase):
    def test_backwards_compatibility(self):
        """Check can re-use an old BioSQL SQLite3 database."""
        original_records = list(SeqIO.parse("GenBank/cor6_6.gb", "gb"))
        # now open a connection to load the database
        server = BioSeqDatabase.open_database(driver=DBDRIVER,
                                              db="BioSQL/cor6_6.db")
        db = server["OLD"]
        self.assertEqual(len(db), len(original_records))
        #Now read them back...
        biosql_records = [db.lookup(name=rec.name)
                          for rec in original_records]
        #And check they agree
        self.assertTrue(compare_records(original_records, biosql_records))

if __name__ == "__main__":
    #Run the test cases
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
