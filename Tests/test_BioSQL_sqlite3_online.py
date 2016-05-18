#!/usr/bin/env python
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Run BioSQL tests using SQLite"""
import os

from Bio import MissingExternalDependencyError
from Bio import SeqIO
from BioSQL import BioSeqDatabase

from common_BioSQL_online import *
import BioSQL_settings
from common_BioSQL import check_config, create_database

# Constants for the database driver
BioSQL_settings.DBHOST = 'localhost'
BioSQL_settings.DBUSER = 'root'
BioSQL_settings.DBPASSWD = ''

BioSQL_settings.DBDRIVER = 'sqlite3'
BioSQL_settings.DBTYPE = 'sqlite'



BioSQL_settings.TESTDB = temp_db_filename()


# This will abort if driver not installed etc:
check_config(BioSQL_settings.DBDRIVER, BioSQL_settings.DBTYPE, BioSQL_settings.DBHOST, BioSQL_settings.DBUSER, BioSQL_settings.DBPASSWD, BioSQL_settings.TESTDB)

# Some of the unit tests don't create their own database,
# so just in case there is no database already:

create_database()

if __name__ == "__main__":
    # Run the test cases
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
