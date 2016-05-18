#!/usr/bin/env python
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Run BioSQL tests using PostgreSQL"""
from Bio import MissingExternalDependencyError
from BioSQL import BioSeqDatabase

from common_BioSQL_online import *
import BioSQL_settings

##################################
# Start of user-editable section #
##################################

# Constants for the database driver
BioSQL_settings.DBHOST = 'localhost'
BioSQL_settings.DBUSER = 'postgres'
BioSQL_settings.DBPASSWD = ''
BioSQL_settings.TESTDB = 'biosql_test'

################################
# End of user-editable section #
################################

DBDRIVER = 'psycopg2'
DBTYPE = 'pg'

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

if __name__ == "__main__":
    # Run the test cases
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
