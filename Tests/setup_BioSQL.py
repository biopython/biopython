#!/usr/bin/env python
"""Preparation for BioSQL tests, setting passwords etc
"""
import os
import Bio
##################################
# Start of user-editable section #
##################################

# You are expected to edit the following lines to match your system.
# The BioSQL unit tests will call this code, and will only run if it works.

# -- MySQL
DBDRIVER = 'MySQLdb'
DBTYPE = 'mysql'
# -- PostgreSQL
#DBDRIVER = 'psycopg'
#DBTYPE = 'pg'

# Constants for the database driver
DBHOST = 'localhost'
DBUSER = 'root'
DBPASSWD = ''
TESTDB = 'biosql_test'

################################
# End of user-editable section #
################################

# Works for mysql and postgresql, not oracle
try:
    DBSCHEMA = "biosqldb-" + DBTYPE + ".sql"
# don't run the tests unless a valid DBTYPE has been set. This
# should be done if you have a MySQL or PostgreSQL database set up and want
# to run the tests. You will also need to set the constants for the database
# driver below.
except NameError:
    message = "Enable tests in Tests/setup_BioSQL.py (not important if you do not plan to use BioSQL)."
    raise Bio.MissingExternalDependencyError(message)

# Uses the SQL file in the Tests/BioSQL directory -- try to keep this current
# with what is going on with BioSQL
SQL_FILE = os.path.join(os.getcwd(), "BioSQL", DBSCHEMA)
assert os.path.isfile(SQL_FILE), "Missing %s" % SQL_FILE
