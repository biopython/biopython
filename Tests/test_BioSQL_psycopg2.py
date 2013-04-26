#!/usr/bin/env python
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Run BioSQL tests using PostgreSQL"""
from Bio import MissingExternalDependencyError
from BioSQL import BioSeqDatabase

from common_BioSQL import *

##################################
# Start of user-editable section #
##################################

# Constants for the database driver
DBHOST = 'localhost'
DBUSER = 'postgres'
DBPASSWD = ''
TESTDB = 'biosql_test'

################################
# End of user-editable section #
################################

DBDRIVER = 'psycopg2'
DBTYPE = 'pg'

#This will abort if driver not installed etc:
check_config(DBDRIVER, DBTYPE, DBHOST, DBUSER, DBPASSWD, TESTDB)

#Some of the unit tests don't create their own database,
#so just in case there is no database already:
create_database()

if __name__ == "__main__":
    #Run the test cases
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
