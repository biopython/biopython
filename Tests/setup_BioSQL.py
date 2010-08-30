#!/usr/bin/env python
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Preparation for BioSQL tests, setting passwords etc
"""
import os
from Bio import MissingExternalDependencyError
from BioSQL import BioSeqDatabase

##################################
# Start of user-editable section #
##################################

# You are expected to edit the following lines to match your system.
# The BioSQL unit tests will call this code, and will only run if it works.

# Constants for the database driver
DBHOST = 'localhost'
DBUSER = 'root'
DBPASSWD = ''
TESTDB = 'biosql_test'

# -- MySQL
#DBDRIVER = 'MySQLdb'
#DBTYPE = 'mysql'
# -- PostgreSQL
#DBDRIVER = 'psycopg2'
#DBTYPE = 'pg'
# -- SQLite
DBDRIVER = 'sqlite3'
DBTYPE = 'sqlite'
TESTDB = os.path.join(os.getcwd(), "BioSQL", "temp_sqlite.db")
# In memory SQLite does not work with current test structure since the tests
# expect databases to be retained between individual tests. 
#TESTDB = ':memory:'

################################
# End of user-editable section #
################################

# Works for mysql and postgresql, not oracle
try:
    DBSCHEMA = "biosqldb-" + DBTYPE + ".sql"
except NameError:
    #This happens if the lines above are commented out
    message = "Enter your settings in Tests/setup_BioSQL.py " \
              "(not important if you do not plan to use BioSQL)."
    raise MissingExternalDependencyError(message)

# Uses the SQL file in the Tests/BioSQL directory -- try to keep this current
# with what is going on with BioSQL
SQL_FILE = os.path.join(os.getcwd(), "BioSQL", DBSCHEMA)
assert os.path.isfile(SQL_FILE), "Missing %s" % SQL_FILE

#Check the database driver is installed:
try:
    __import__(DBDRIVER)
except ImportError:
    message = "Install %s or correct Tests/setup_BioSQL.py "\
              "(not important if you do not plan to use BioSQL)." % DBDRIVER
    raise MissingExternalDependencyError(message)

#Could check the username, password and host work here,
#but this only seems to work for the first unit test
#that tries to import this file.
