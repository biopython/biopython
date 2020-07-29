# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Run BioSQL tests using SQLite."""

import unittest

# Really do want "import *" to get all the test clases:
from common_BioSQL import *  # noqa: F403
from common_BioSQL_online import *  # noqa: F403

# Import these explicitly to avoid flake8 F405 below
from common_BioSQL import load_biosql_ini, check_config, temp_db_filename
from common_BioSQL_online import share_config

import requires_internet

requires_internet.check()

# Constants for the database driver
DBDRIVER = "sqlite3"
DBTYPE = "sqlite"

DBHOST = None
DBUSER = "root"
DBPASSWD = None
TESTDB = temp_db_filename()

# This will abort if driver not installed etc:
check_config(DBDRIVER, DBTYPE, DBHOST, DBUSER, DBPASSWD, TESTDB)
share_config(DBDRIVER, DBTYPE, DBHOST, DBUSER, DBPASSWD, TESTDB)

if __name__ == "__main__":
    # Run the test cases
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
