#!/usr/bin/env python
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Run BioSQL tests using SQLite"""
import os

from Bio import MissingExternalDependencyError
from Bio import SeqIO
from BioSQL import BioSeqDatabase

from common_BioSQL import *
import BioSQL_settings

# Constants for the database driver
BioSQL_settings.DBHOST = 'localhost'
BioSQL_settings.DBUSER = 'root'
BioSQL_settings.DBPASSWD = ''

BioSQL_settings.DBDRIVER = 'sqlite3'
BioSQL_settings.DBTYPE = 'sqlite'

BioSQL_settings.TESTDB = temp_db_filename()


# This will abort if driver not installed etc:
check_config(BioSQL_settings.DBDRIVER,
             BioSQL_settings.DBTYPE,
             BioSQL_settings.DBHOST,
             BioSQL_settings.DBUSER,
             BioSQL_settings.DBPASSWD,
             BioSQL_settings.TESTDB)

# Some of the unit tests don't create their own database,
# so just in case there is no database already:
create_database()


if False:
    # This is how I generated test file Tests/BioSQL/cor6_6.db
    # which is test cross-checked with the latest bindings to
    # catch any regressions in how we map GenBank entries to
    # the database.
    assert not os.path.isfile("BioSQL/cor6_6.db")
    server = BioSeqDatabase.open_database(driver=BioSQL_settings.DBDRIVER,
                                          db="BioSQL/cor6_6.db")
    BioSQL_settings.DBSCHEMA = "biosqldb-" + BioSQL_settings.DBTYPE + ".sql"
    BioSQL_settings.SQL_FILE = os.path.join(os.getcwd(), "BioSQL", BioSQL_settings.DBSCHEMA)
    assert os.path.isfile(BioSQL_settings.SQL_FILE), BioSQL_settings.SQL_FILE
    server.load_database_sql(BioSQL_settings.SQL_FILE)
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
        server = BioSeqDatabase.open_database(driver=BioSQL_settings.DBDRIVER,
                                              db="BioSQL/cor6_6.db")
        db = server["OLD"]
        self.assertEqual(len(db), len(original_records))
        # Now read them back...
        biosql_records = [db.lookup(name=rec.name)
                          for rec in original_records]
        # And check they agree
        self.assertTrue(compare_records(original_records, biosql_records))

if __name__ == "__main__":
    # Run the test cases
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
