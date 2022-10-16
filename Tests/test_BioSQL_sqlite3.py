# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Run BioSQL tests using SQLite."""

import os
import unittest

from Bio import SeqIO
from BioSQL import BioSeqDatabase

from seq_tests_common import SeqRecordTestBaseClass

# Really do want "import *" to get all the test classes:
from common_BioSQL import *  # noqa: F403

# Import these explicitly to avoid flake8 F405 below:
from common_BioSQL import load_biosql_ini, check_config, temp_db_filename

# Constants for the database driver
DBDRIVER = "sqlite3"
DBTYPE = "sqlite"

DBHOST = None
DBUSER = "root"
DBPASSWD = None
TESTDB = temp_db_filename()

# This will abort if driver not installed etc:
check_config(DBDRIVER, DBTYPE, DBHOST, DBUSER, DBPASSWD, TESTDB)


if False:
    # This is how I generated test file Tests/BioSQL/cor6_6.db
    # which is test cross-checked with the latest bindings to
    # catch any regressions in how we map GenBank entries to
    # the database.
    assert not os.path.isfile("BioSQL/cor6_6.db")
    server = BioSeqDatabase.open_database(driver=DBDRIVER, db="BioSQL/cor6_6.db")
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


class BackwardsCompatibilityTest(SeqRecordTestBaseClass):
    def test_backwards_compatibility(self):
        """Check can re-use an old BioSQL SQLite3 database."""
        original_records = []
        for record in SeqIO.parse("GenBank/cor6_6.gb", "gb"):
            if record.annotations["molecule_type"] == "mRNA":
                record.annotations["molecule_type"] = "DNA"
            original_records.append(record)
        # now open a connection to load the database
        server = BioSeqDatabase.open_database(driver=DBDRIVER, db="BioSQL/cor6_6.db")
        db = server["OLD"]
        self.assertEqual(len(db), len(original_records))
        # Now read them back...
        biosql_records = [db.lookup(name=rec.name) for rec in original_records]
        # And check they agree
        self.compare_records(original_records, biosql_records)
        server.close()


if __name__ == "__main__":
    # Run the test cases
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
