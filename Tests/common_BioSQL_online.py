# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for dealing with storage of biopython objects in a relational db.
"""
from __future__ import print_function

import os
import platform
import unittest
import tempfile
import time

from Bio._py3k import StringIO
from Bio._py3k import zip
from Bio._py3k import basestring

# Hide annoying warnings from things like bonds in GenBank features,
# or PostgreSQL schema rules. TODO - test these warnings are raised!
import warnings
from Bio import BiopythonWarning

# local stuff
from Bio import MissingExternalDependencyError
from Bio.Seq import Seq, MutableSeq
from Bio.SeqFeature import SeqFeature
from Bio import Alphabet
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from BioSQL import BioSeqDatabase
from BioSQL import BioSeq
from Bio import Entrez

from seq_tests_common import compare_record, compare_records

import requires_internet

if __name__ == "__main__":
    raise RuntimeError("Call this via test_BioSQL_*.py not directly")

# Exporting these to the test_BioSQL_XXX.py files which import this file:
# DBDRIVER, DBTYPE, DBHOST, DBUSER, DBPASSWD, TESTDB, DBSCHEMA, SQL_FILE, SYSTEM

SYSTEM = platform.system()


def temp_db_filename():
    # In memory SQLite does not work with current test structure since the tests
    # expect databases to be retained between individual tests.
    # TESTDB = ':memory:'
    # Instead, we use (if we can) /dev/shm
    try:
        h, test_db_fname = tempfile.mkstemp("_BioSQL.db", dir='/dev/shm')
    except OSError:
        # We can't use /dev/shm
        h, test_db_fname = tempfile.mkstemp("_BioSQL.db")
    os.close(h)
    return test_db_fname


def check_config(dbdriver, dbtype, dbhost, dbuser, dbpasswd, testdb):
    global DBDRIVER, DBTYPE, DBHOST, DBUSER, DBPASSWD, TESTDB, DBSCHEMA
    global SYSTEM, SQL_FILE
    DBDRIVER = dbdriver
    DBTYPE = dbtype
    DBHOST = dbhost
    DBUSER = dbuser
    DBPASSWD = dbpasswd
    TESTDB = testdb

    # Check the database driver is installed:
    if SYSTEM == "Java":
        try:
            if DBDRIVER in ["MySQLdb"]:
                import com.mysql.jdbc.Driver
            elif DBDRIVER in ["psycopg2"]:
                import org.postgresql.Driver
        except ImportError:
            message = "Install the JDBC driver for %s to use BioSQL " % DBTYPE
            raise MissingExternalDependencyError(message)
    else:
        try:
            __import__(DBDRIVER)
        except ImportError:
            message = "Install %s if you want to use %s with BioSQL " % (DBDRIVER, DBTYPE)
            raise MissingExternalDependencyError(message)

    try:
        if DBDRIVER in ["sqlite3"]:
            server = BioSeqDatabase.open_database(driver=DBDRIVER, db=TESTDB)
        else:
            server = BioSeqDatabase.open_database(driver=DBDRIVER, host=DBHOST,
                                                  user=DBUSER, passwd=DBPASSWD)
        server.close()
        del server
    except Exception as e:
        message = "Connection failed, check settings if you plan to use BioSQL: %s" % e
        raise MissingExternalDependencyError(message)

    DBSCHEMA = "biosqldb-" + DBTYPE + ".sql"
    SQL_FILE = os.path.join(os.getcwd(), "BioSQL", DBSCHEMA)

    if not os.path.isfile(SQL_FILE):
        message = "Missing SQL schema file: %s" % SQL_FILE
        raise MissingExternalDependencyError(message)


def _do_db_create():
    """Do the actual work of database creation.

    Relevant for MySQL and PostgreSQL.
    """
    # first open a connection to create the database
    server = BioSeqDatabase.open_database(driver=DBDRIVER, host=DBHOST,
                                          user=DBUSER, passwd=DBPASSWD)

    if DBDRIVER == "pgdb":
        # The pgdb postgres driver does not support autocommit, so here we
        # commit the current transaction so that 'drop database' query will
        # be outside a transaction block
        server.adaptor.cursor.execute("COMMIT")
    else:
        # Auto-commit: postgresql cannot drop database in a transaction
        try:
            server.adaptor.autocommit()
        except AttributeError:
            pass

    # drop anything in the database
    try:
        # with Postgres, can get errors about database still being used and
        # not able to be dropped. Wait briefly to be sure previous tests are
        # done with it.
        time.sleep(1)
        sql = r"DROP DATABASE " + TESTDB
        server.adaptor.cursor.execute(sql, ())
    except (server.module.OperationalError,
            server.module.Error,
            server.module.DatabaseError) as e:  # the database doesn't exist
        pass
    except (server.module.IntegrityError,
            server.module.ProgrammingError) as e:  # ditto--perhaps
        if str(e).find('database "%s" does not exist' % TESTDB) == -1:
            server.close()
            raise
    # create a new database
    sql = r"CREATE DATABASE " + TESTDB
    server.adaptor.execute(sql, ())
    server.close()


def create_database():
    """Delete any existing BioSQL test database, then (re)create an empty BioSQL database."""
    if DBDRIVER in ["sqlite3"]:
        global TESTDB
        if os.path.exists(TESTDB):
            try:
                os.remove(TESTDB)
            except:
                time.sleep(1)
                try:
                    os.remove(TESTDB)
                except Exception:
                    # Seen this with PyPy 2.1 (and older) on Windows -
                    # which suggests an open handle still exists?
                    print("Could not remove %r" % TESTDB)
                    pass
        # Now pick a new filename - just in case there is a stale handle
        # (which might be happening under Windows...)
        TESTDB = temp_db_filename()
    else:
        _do_db_create()

    # now open a connection to load the database
    server = BioSeqDatabase.open_database(driver=DBDRIVER,
                                          user=DBUSER, passwd=DBPASSWD,
                                          host=DBHOST, db=TESTDB)
    try:
        server.load_database_sql(SQL_FILE)
        server.commit()
        server.close()
    except:
        # Failed, but must close the handle...
        server.close()
        raise


def destroy_database():
    """Delete any temporary BioSQL sqlite3 database files."""
    if DBDRIVER in ["sqlite3"]:
        if os.path.exists(TESTDB):
            os.remove(TESTDB)


def load_database(gb_filename_or_handle):
    """Load a GenBank file into a new BioSQL database.

    This is useful for running tests against a newly created database.
    """

    create_database()
    # now open a connection to load the database
    db_name = "biosql-test"
    server = BioSeqDatabase.open_database(driver=DBDRIVER,
                                          user=DBUSER, passwd=DBPASSWD,
                                          host=DBHOST, db=TESTDB)
    db = server.new_database(db_name)

    # get the GenBank file we are going to put into it
    iterator = SeqIO.parse(gb_filename_or_handle, "gb")
    # finally put it in the database
    count = db.load(iterator)
    server.commit()
    server.close()
    return count


class TaxonomyTest(unittest.TestCase):
    """Test proper insertion and retrieval of taxonomy data
    """
    def setUp(self):
        Entrez.email = "biopython-dev@biopython.org"
        # create TESTDB
        create_database()

        # load the database
        db_name = "biosql-test"
        self.server = BioSeqDatabase.open_database(driver=DBDRIVER,
                                                   user=DBUSER, passwd=DBPASSWD,
                                                   host=DBHOST, db=TESTDB)

        # remove the database if it already exists
        try:
            self.server[db_name]
            self.server.remove_database(db_name)
        except KeyError:
            pass

        self.db = self.server.new_database(db_name)

        # get the GenBank file we are going to put into it
        self.iterator = SeqIO.parse("GenBank/cor6_6.gb", "gb")

    def tearDown(self):
        self.server.close()
        destroy_database()
        del self.db
        del self.server

    def test_taxon_left_right_values(self):
        self.db.load(self.iterator, True)
        sql = """SELECT DISTINCT include.ncbi_taxon_id FROM taxon
                  INNER JOIN taxon AS include ON
                      (include.left_value BETWEEN taxon.left_value
                                  AND taxon.right_value)
                  WHERE taxon.taxon_id IN
                      (SELECT taxon_id FROM taxon_name
                                  WHERE name = 'Brassicales')
                      AND include.right_value - include.left_value = 1"""

        rows = self.db.adaptor.execute_and_fetchall(sql)
        self.assertEqual(4, len(rows))
        values = set()
        for row in rows:
            values.add(row[0])
        self.assertEqual(set([3704, 3711, 3708, 3702]), set(values))


