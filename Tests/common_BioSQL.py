# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Dealing with storage of biopython objects in a BioSQL relational db."""

import configparser
import os
import platform
import tempfile
import time
import unittest

from io import StringIO

# Hide annoying warnings from things like bonds in GenBank features,
# or PostgreSQL schema rules. TODO - test these warnings are raised!
import warnings
from Bio import BiopythonWarning

# local stuff
from Bio import MissingExternalDependencyError
from Bio.Seq import Seq, MutableSeq
from Bio.SeqFeature import SeqFeature, UnknownPosition, ExactPosition
from Bio import Alphabet
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from BioSQL import BioSeqDatabase
from BioSQL import BioSeq

from seq_tests_common import compare_record, compare_records

if __name__ == "__main__":
    raise RuntimeError("Call this via test_BioSQL_*.py not directly")

# Exporting these to the test_BioSQL_XXX.py files which import this file:
# DBDRIVER, DBTYPE, DBHOST, DBUSER, DBPASSWD, TESTDB, DBSCHEMA, SQL_FILE, SYSTEM

SYSTEM = platform.system()


def load_biosql_ini(DBTYPE):
    """Load the database settings from INI file."""
    if not os.path.isfile("biosql.ini"):
        raise MissingExternalDependencyError(
            "BioSQL test configuration file biosql.ini missing (see biosql.ini.sample)"
        )

    config = configparser.ConfigParser()
    config.read("biosql.ini")
    DBHOST = config.get(DBTYPE, "dbhost")
    DBUSER = config.get(DBTYPE, "dbuser")
    DBPASSWD = config.get(DBTYPE, "dbpasswd")
    TESTDB = config.get(DBTYPE, "testdb")
    return DBHOST, DBUSER, DBPASSWD, TESTDB


def temp_db_filename():
    """Generate a temporary filename for SQLite database."""
    # In memory SQLite does not work with current test structure since the tests
    # expect databases to be retained between individual tests.
    # TESTDB = ':memory:'
    # Instead, we use (if we can) /dev/shm
    try:
        h, test_db_fname = tempfile.mkstemp("_BioSQL.db", dir="/dev/shm")
    except OSError:
        # We can't use /dev/shm
        h, test_db_fname = tempfile.mkstemp("_BioSQL.db")
    os.close(h)
    return test_db_fname


def check_config(dbdriver, dbtype, dbhost, dbuser, dbpasswd, testdb):
    """Verify the database settings work for connecting."""
    global DBDRIVER, DBTYPE, DBHOST, DBUSER, DBPASSWD, TESTDB, DBSCHEMA
    global SYSTEM, SQL_FILE
    DBDRIVER = dbdriver
    DBTYPE = dbtype
    DBHOST = dbhost
    DBUSER = dbuser
    DBPASSWD = dbpasswd
    TESTDB = testdb

    if not DBDRIVER or not DBTYPE or not DBUSER:
        # No point going any further...
        raise MissingExternalDependencyError("Incomplete BioSQL test settings")

    # Check the database driver is installed:
    if SYSTEM == "Java":
        try:
            if DBDRIVER in ["MySQLdb"]:
                import com.mysql.jdbc.Driver
            elif DBDRIVER in ["psycopg2", "pgdb"]:
                import org.postgresql.Driver
        except ImportError:
            message = "Install the JDBC driver for %s to use BioSQL " % DBTYPE
            raise MissingExternalDependencyError(message) from None
    else:
        try:
            __import__(DBDRIVER)
        except ImportError:
            if DBDRIVER in ["MySQLdb"]:
                message = (
                    "Install MySQLdb or mysqlclient if you want to use %s with BioSQL "
                    % (DBTYPE)
                )
            else:
                message = "Install %s if you want to use %s with BioSQL " % (
                    DBDRIVER,
                    DBTYPE,
                )
            raise MissingExternalDependencyError(message) from None

    try:
        if DBDRIVER in ["sqlite3"]:
            server = BioSeqDatabase.open_database(driver=DBDRIVER, db=TESTDB)
        else:
            server = BioSeqDatabase.open_database(
                driver=DBDRIVER, host=DBHOST, user=DBUSER, passwd=DBPASSWD
            )
        server.close()
        del server
    except Exception as e:
        message = "Connection failed, check settings if you plan to use BioSQL: %s" % e
        raise MissingExternalDependencyError(message) from None

    DBSCHEMA = "biosqldb-" + DBTYPE + ".sql"
    SQL_FILE = os.path.join(os.getcwd(), "BioSQL", DBSCHEMA)

    if not os.path.isfile(SQL_FILE):
        message = "Missing SQL schema file: %s" % SQL_FILE
        raise MissingExternalDependencyError(message)


def _do_db_cleanup():
    """Cleanup everything from TESTDB.

    Relevant for MySQL and PostgreSQL.
    """
    if DBDRIVER in ["psycopg2", "pgdb"]:
        # first open a connection the database
        # notice that postgres doesn't have createdb privileges, so
        # the TESTDB must exist
        server = BioSeqDatabase.open_database(
            driver=DBDRIVER, host=DBHOST, user=DBUSER, passwd=DBPASSWD, db=TESTDB
        )

        # The pgdb postgres driver does not support autocommit, so here we
        # commit the current transaction so that 'drop database' query will
        # be outside a transaction block
        server.adaptor.cursor.execute("COMMIT")
        # drop anything in the database
        # with Postgres, can get errors about database still being used.
        # Wait briefly to be sure previous tests are done with it.
        time.sleep(1)
        # drop anything in the database
        sql = r"DROP OWNED BY " + DBUSER
        server.adaptor.cursor.execute(sql, ())
        server.close()
    else:
        # first open a connection to create the database
        server = BioSeqDatabase.open_database(
            driver=DBDRIVER, host=DBHOST, user=DBUSER, passwd=DBPASSWD
        )
        # Auto-commit
        try:
            server.adaptor.autocommit()
        except AttributeError:
            pass
        # drop the database
        try:
            sql = r"DROP DATABASE " + TESTDB
            server.adaptor.cursor.execute(sql, ())
        except (
            server.module.OperationalError,
            server.module.Error,
            server.module.DatabaseError,
        ) as e:  # the database doesn't exist
            pass
        except (
            server.module.IntegrityError,
            server.module.ProgrammingError,
        ) as e:  # ditto--perhaps
            if str(e).find('database "%s" does not exist' % TESTDB) == -1:
                server.close()
                raise
        # create a new database
        sql = r"CREATE DATABASE " + TESTDB
        server.adaptor.execute(sql, ())
        server.close()


def create_database():
    """Delete any existing BioSQL test DB, then (re)create an empty BioSQL DB.

    Returns TESTDB name which will change for for SQLite.
    """
    if DBDRIVER in ["sqlite3"]:
        global TESTDB
        if os.path.exists(TESTDB):
            try:
                os.remove(TESTDB)
            except Exception:
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
        _do_db_cleanup()

    # now open a connection to load the database
    server = BioSeqDatabase.open_database(
        driver=DBDRIVER, user=DBUSER, passwd=DBPASSWD, host=DBHOST, db=TESTDB
    )
    try:
        server.load_database_sql(SQL_FILE)
        server.commit()
        server.close()
    except Exception:
        # Failed, but must close the handle...
        server.close()
        raise

    return TESTDB


def destroy_database():
    """Delete any temporary BioSQL sqlite3 database files."""
    if DBDRIVER in ["sqlite3"]:
        if os.path.exists(TESTDB):
            os.remove(TESTDB)


def load_database(gb_filename_or_handle):
    """Load a GenBank file into a new BioSQL database.

    This is useful for running tests against a newly created database.
    """
    TESTDB = create_database()
    # now open a connection to load the database
    db_name = "biosql-test"
    server = BioSeqDatabase.open_database(
        driver=DBDRIVER, user=DBUSER, passwd=DBPASSWD, host=DBHOST, db=TESTDB
    )
    db = server.new_database(db_name)

    # get the GenBank file we are going to put into it
    iterator = SeqIO.parse(gb_filename_or_handle, "gb")
    # finally put it in the database
    count = db.load(iterator)
    server.commit()
    server.close()
    return count


def load_multi_database(gb_filename_or_handle, gb_filename_or_handle2):
    """Load two GenBank files into a new BioSQL database as different subdatabases.

    This is useful for running tests against a newly created database.
    """
    TESTDB = create_database()
    # now open a connection to load the database
    db_name = "biosql-test"
    db_name2 = "biosql-test2"
    server = BioSeqDatabase.open_database(
        driver=DBDRIVER, user=DBUSER, passwd=DBPASSWD, host=DBHOST, db=TESTDB
    )
    db = server.new_database(db_name)

    # get the GenBank file we are going to put into it
    iterator = SeqIO.parse(gb_filename_or_handle, "gb")
    count = db.load(iterator)

    db = server.new_database(db_name2)

    # get the GenBank file we are going to put into it
    iterator = SeqIO.parse(gb_filename_or_handle2, "gb")
    # finally put it in the database
    count2 = db.load(iterator)
    server.commit()

    server.close()
    return count + count2


class MultiReadTest(unittest.TestCase):
    """Test reading a database with multiple namespaces."""

    loaded_db = 0

    def setUp(self):
        """Connect to and load up the database."""
        load_multi_database("GenBank/cor6_6.gb", "GenBank/NC_000932.gb")

        self.server = BioSeqDatabase.open_database(
            driver=DBDRIVER, user=DBUSER, passwd=DBPASSWD, host=DBHOST, db=TESTDB
        )

        self.db = self.server["biosql-test"]
        self.db2 = self.server["biosql-test2"]

    def tearDown(self):
        self.server.close()
        destroy_database()
        del self.db
        del self.db2
        del self.server

    def test_server(self):
        """Check BioSeqDatabase methods."""
        server = self.server
        self.assertIn("biosql-test", server)
        self.assertIn("biosql-test2", server)
        self.assertEqual(2, len(server))
        self.assertEqual(["biosql-test", "biosql-test2"], list(server.keys()))
        # Check we can delete the namespace...
        del server["biosql-test"]
        del server["biosql-test2"]
        self.assertEqual(0, len(server))
        with self.assertRaises(KeyError):
            del server["non-existant-name"]

    def test_get_db_items(self):
        """Check list, keys, length etc."""
        db = self.db
        items = list(db.values())
        keys = list(db)
        length = len(items)
        self.assertEqual(length, len(db))
        self.assertEqual(length, len(list(db)))
        self.assertEqual(length, len(list(db.items())))
        self.assertEqual(length, len(list(db.keys())))
        self.assertEqual(length, len(list(db.values())))
        for (k1, r1), (k2, r2) in zip(zip(keys, items), db.items()):
            self.assertEqual(k1, k2)
            self.assertEqual(r1.id, r2.id)
        for k in keys:
            del db[k]
        self.assertEqual(0, len(db))
        with self.assertRaises(KeyError):
            del db["non-existant-name"]

    def test_cross_retrieval_of_items(self):
        """Test that valid ids can't be retrieved between namespaces."""
        db = self.db
        db2 = self.db2
        for db2_id in db2.keys():
            with self.assertRaises(KeyError):
                rec = db[db2_id]


class ReadTest(unittest.TestCase):
    """Test reading a database from an already built database."""

    loaded_db = 0

    def setUp(self):
        """Connect to and load up the database."""
        load_database("GenBank/cor6_6.gb")

        self.server = BioSeqDatabase.open_database(
            driver=DBDRIVER, user=DBUSER, passwd=DBPASSWD, host=DBHOST, db=TESTDB
        )

        self.db = self.server["biosql-test"]

    def tearDown(self):
        self.server.close()
        destroy_database()
        del self.db
        del self.server

    def test_server(self):
        """Check BioSeqDatabase methods."""
        server = self.server
        self.assertIn("biosql-test", server)
        self.assertEqual(1, len(server))
        self.assertEqual(["biosql-test"], list(server.keys()))
        # Check we can delete the namespace...
        del server["biosql-test"]
        self.assertEqual(0, len(server))
        with self.assertRaises(KeyError):
            del server["non-existant-name"]

    def test_get_db_items(self):
        """Check list, keys, length etc."""
        db = self.db
        items = list(db.values())
        keys = list(db)
        length = len(items)
        self.assertEqual(length, len(db))
        self.assertEqual(length, len(list(db.items())))
        self.assertEqual(length, len(list(db)))
        self.assertEqual(length, len(list(db.values())))
        for (k1, r1), (k2, r2) in zip(zip(keys, items), db.items()):
            self.assertEqual(k1, k2)
            self.assertEqual(r1.id, r2.id)
        for k in keys:
            del db[k]
        self.assertEqual(0, len(db))
        with self.assertRaises(KeyError):
            del db["non-existant-name"]

    def test_lookup_items(self):
        """Test retrieval of items using various ids."""
        self.db.lookup(accession="X62281")
        try:
            self.db.lookup(accession="Not real")
            raise AssertionError("No problem on fake id retrieval")
        except IndexError:
            pass
        self.db.lookup(display_id="ATKIN2")
        try:
            self.db.lookup(display_id="Not real")
            raise AssertionError("No problem on fake id retrieval")
        except IndexError:
            pass

        # primary id retrieval
        self.db.lookup(primary_id="16353")
        try:
            self.db.lookup(primary_id="Not Real")
            raise AssertionError("No problem on fake primary id retrieval")
        except IndexError:
            pass


class SeqInterfaceTest(unittest.TestCase):
    """Make sure the BioSQL objects implement the expected biopython interface."""

    def setUp(self):
        """Load a database."""
        load_database("GenBank/cor6_6.gb")

        self.server = BioSeqDatabase.open_database(
            driver=DBDRIVER, user=DBUSER, passwd=DBPASSWD, host=DBHOST, db=TESTDB
        )
        self.db = self.server["biosql-test"]
        self.item = self.db.lookup(accession="X62281")

    def tearDown(self):
        self.server.close()
        destroy_database()
        del self.db
        del self.item
        del self.server

    def test_seq_record(self):
        """Make sure SeqRecords from BioSQL implement the right interface."""
        test_record = self.item
        self.assertIsInstance(test_record.seq, BioSeq.DBSeq)
        self.assertEqual(test_record.id, "X62281.1", test_record.id)
        self.assertEqual(test_record.name, "ATKIN2")
        self.assertEqual(test_record.description, "A.thaliana kin2 gene")
        self.assertTrue(hasattr(test_record, "annotations"))
        # XXX should do something with annotations once they are like
        # a dictionary
        for feature in test_record.features:
            self.assertIsInstance(feature, SeqFeature)
        # shouldn't cause any errors!
        self.assertIsInstance(str(test_record), str)
        # Confirm can delete annotations etc to test these properties
        del test_record.annotations
        del test_record.dbxrefs
        del test_record.features
        del test_record.seq

    def test_seq(self):
        """Make sure Seqs from BioSQL implement the right interface."""
        test_seq = self.item.seq
        alphabet = test_seq.alphabet
        self.assertIsInstance(alphabet, Alphabet.Alphabet)
        data = test_seq.data
        self.assertEqual(type(data), type(""))
        string_rep = str(test_seq)
        self.assertEqual(string_rep, str(test_seq))  # check __str__ too
        self.assertEqual(type(string_rep), type(""))
        self.assertEqual(len(test_seq), 880)
        self.assertEqual(test_seq[879], "A")
        self.assertEqual(test_seq[-1], "A")
        self.assertEqual(test_seq[0], "A")
        self.assertEqual(test_seq[-880], "A")
        self.assertRaises(IndexError, test_seq.__getitem__, 880)
        self.assertRaises(IndexError, test_seq.__getitem__, -881)
        self.assertRaises(TypeError, test_seq.__getitem__, None)

    def test_convert(self):
        """Check can turn a DBSeq object into a Seq or MutableSeq."""
        test_seq = self.item.seq

        other = test_seq.toseq()
        self.assertEqual(str(test_seq), str(other))
        self.assertIsInstance(other, Seq)

        other = test_seq.tomutable()
        self.assertEqual(str(test_seq), str(other))
        self.assertIsInstance(other, MutableSeq)

    def test_addition(self):
        """Check can add DBSeq objects together."""
        test_seq = self.item.seq
        for other in [
            Seq("ACGT", test_seq.alphabet),
            MutableSeq("ACGT", test_seq.alphabet),
            "ACGT",
            test_seq,
        ]:
            test = test_seq + other
            self.assertEqual(str(test), str(test_seq) + str(other))
            self.assertIsInstance(test, Seq)
            test = other + test_seq
            self.assertEqual(str(test), str(other) + str(test_seq))

    def test_multiplication(self):
        """Check can multiply DBSeq objects by integers."""
        test_seq = self.item.seq
        alphabet = test_seq.alphabet
        tripled = test_seq * 3
        # Test DBSeq.__mul__
        self.assertIsInstance(tripled, Seq)
        self.assertNotIsInstance(tripled, BioSeq.DBSeq)
        self.assertEqual(tripled, str(test_seq) * 3)
        # Test DBSeq.__rmul__
        tripled = 3 * test_seq
        self.assertIsInstance(tripled, Seq)
        self.assertNotIsInstance(tripled, BioSeq.DBSeq)
        self.assertEqual(tripled, str(test_seq) * 3)
        # Test DBSeq.__imul__
        original = self.item.seq
        tripled = test_seq
        tripled *= 3
        self.assertIsInstance(tripled, Seq)
        self.assertNotIsInstance(tripled, BioSeq.DBSeq)
        self.assertEqual(tripled, str(original) * 3)

    def test_seq_slicing(self):
        """Check that slices of sequences are retrieved properly."""
        test_seq = self.item.seq
        new_seq = test_seq[:10]
        self.assertIsInstance(new_seq, BioSeq.DBSeq)
        # simple slicing
        self.assertEqual(str(test_seq[:5]), "ATTTG")
        self.assertEqual(str(test_seq[0:5]), "ATTTG")
        self.assertEqual(str(test_seq[2:3]), "T")
        self.assertEqual(str(test_seq[2:4]), "TT")
        self.assertEqual(str(test_seq[870:]), "TTGAATTATA")
        # getting more fancy
        self.assertEqual(test_seq[-1], "A")
        self.assertEqual(test_seq[1], "T")
        self.assertEqual(str(test_seq[-10:][5:]), "TTATA")
        self.assertEqual(str(test_seq[-10:][5:]), "TTATA")

    def test_record_slicing(self):
        """Check that slices of DBSeqRecord are retrieved properly."""
        new_rec = self.item[400:]
        self.assertIsInstance(new_rec, SeqRecord)
        self.assertEqual(len(new_rec), 480)
        self.assertEqual(len(new_rec.features), 5)

    def test_seq_features(self):
        """Check SeqFeatures of a sequence."""
        test_features = self.item.features
        cds_feature = test_features[6]
        self.assertEqual(cds_feature.type, "CDS")
        self.assertEqual(
            str(cds_feature.location), "join{[103:160](+), [319:390](+), [503:579](+)}"
        )

        try:
            self.assertEqual(cds_feature.qualifiers["gene"], ["kin2"])
            self.assertEqual(cds_feature.qualifiers["protein_id"], ["CAA44171.1"])
            self.assertEqual(cds_feature.qualifiers["codon_start"], ["1"])
        except KeyError:
            raise KeyError(
                "Missing expected entries, have %s" % repr(cds_feature.qualifiers)
            ) from None

        self.assertIn("db_xref", cds_feature.qualifiers)
        multi_ann = cds_feature.qualifiers["db_xref"]
        self.assertEqual(len(multi_ann), 2)
        self.assertIn("GI:16354", multi_ann)
        self.assertIn("SWISS-PROT:P31169", multi_ann)


class LoaderTest(unittest.TestCase):
    """Load a database from a GenBank file."""

    def setUp(self):
        # create TESTDB
        TESTDB = create_database()

        # load the database
        db_name = "biosql-test"
        self.server = BioSeqDatabase.open_database(
            driver=DBDRIVER, user=DBUSER, passwd=DBPASSWD, host=DBHOST, db=TESTDB
        )

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

    def test_load_database(self):
        """Load SeqRecord objects into a BioSQL database."""
        self.db.load(self.iterator)

        # do some simple tests to make sure we actually loaded the right
        # thing. More advanced tests in a different module.
        items = list(self.db.values())
        self.assertEqual(len(items), 6)
        self.assertEqual(len(self.db), 6)
        item_names = []
        item_ids = []
        for item in items:
            item_names.append(item.name)
            item_ids.append(item.id)
        item_names.sort()
        item_ids.sort()
        self.assertEqual(
            item_names,
            ["AF297471", "ARU237582", "ATCOR66M", "ATKIN2", "BNAKINI", "BRRBIF72"],
        )
        self.assertEqual(
            item_ids,
            [
                "AF297471.1",
                "AJ237582.1",
                "L31939.1",
                "M81224.1",
                "X55053.1",
                "X62281.1",
            ],
        )


class DeleteTest(unittest.TestCase):
    """Test proper deletion of entries from a database."""

    loaded_db = 0

    def setUp(self):
        """Connect to and load up the database."""
        load_database("GenBank/cor6_6.gb")

        self.server = BioSeqDatabase.open_database(
            driver=DBDRIVER, user=DBUSER, passwd=DBPASSWD, host=DBHOST, db=TESTDB
        )

        self.db = self.server["biosql-test"]

    def tearDown(self):
        self.server.close()
        destroy_database()
        del self.db
        del self.server

    def test_server(self):
        """Check BioSeqDatabase methods."""
        server = self.server
        self.assertIn("biosql-test", server)
        self.assertEqual(1, len(server))
        self.assertEqual(["biosql-test"], list(server.keys()))
        # Check we can delete the namespace...
        del server["biosql-test"]
        self.assertEqual(0, len(server))
        with self.assertRaises(KeyError):
            del server["non-existant-name"]

    def test_del_db_items(self):
        """Check all associated data is deleted from an item."""
        db = self.db
        items = list(db.values())
        keys = list(db)
        length = len(items)

        for seq_id in keys:
            sql = "SELECT seqfeature_id from seqfeature where bioentry_id = '%s'"
            # get the original number of seqfeatures associated with the bioentry
            seqfeatures = self.db.adaptor.execute_and_fetchall(sql % (seq_id))

            del db[seq_id]
            # check to see that the entry in the bioentry table is removed
            self.assertEqual(seq_id in db, False)

            # no need to check seqfeature presence if it had none to begin with
            if len(seqfeatures):
                rows_d = self.db.adaptor.execute_and_fetchall(sql % (seq_id))
                # check to see that associated data is removed
                self.assertEqual(len(rows_d), 0)

        self.assertEqual(0, len(list(db.values())))


class DupLoadTest(unittest.TestCase):
    """Check a few duplicate conditions fail."""

    def setUp(self):
        # drop any old database and create a new one:
        TESTDB = create_database()
        # connect to new database:
        self.server = BioSeqDatabase.open_database(
            driver=DBDRIVER, user=DBUSER, passwd=DBPASSWD, host=DBHOST, db=TESTDB
        )
        # Create new namespace within new empty database:
        self.db = self.server.new_database("biosql-test")

    def tearDown(self):
        self.server.rollback()
        self.server.close()
        destroy_database()
        del self.db
        del self.server

    def test_duplicate_load(self):
        """Make sure can't import a single record twice (in one go)."""
        record = SeqRecord(Seq("ATGCTATGACTAT", Alphabet.generic_dna), id="Test1")
        try:
            count = self.db.load([record, record])
        except Exception as err:
            # Good!
            # Note we don't do a specific exception handler because the
            # exception class will depend on which DB back end is in use.
            self.assertTrue(
                err.__class__.__name__
                in [
                    "IntegrityError",
                    "UniqueViolation",
                    "AttributeError",
                    "OperationalError",
                ],
                err.__class__.__name__,
            )
            return
        raise Exception("Should have failed! Loaded %i records" % count)

    def test_duplicate_load2(self):
        """Make sure can't import a single record twice (in steps)."""
        record = SeqRecord(Seq("ATGCTATGACTAT", Alphabet.generic_dna), id="Test2")
        count = self.db.load([record])
        self.assertEqual(count, 1)
        try:
            count = self.db.load([record])
        except Exception as err:
            # Good!
            self.assertTrue(
                err.__class__.__name__
                in ["IntegrityError", "UniqueViolation", "AttributeError"],
                err.__class__.__name__,
            )
            return
        raise Exception("Should have failed! Loaded %i records" % count)

    def test_duplicate_id_load(self):
        """Make sure can't import records with same ID (in one go)."""
        record1 = SeqRecord(Seq("ATGCTATGACTAT", Alphabet.generic_dna), id="TestA")
        record2 = SeqRecord(Seq("GGGATGCGACTAT", Alphabet.generic_dna), id="TestA")
        try:
            count = self.db.load([record1, record2])
        except Exception as err:
            # Good!
            self.assertTrue(
                err.__class__.__name__
                in ["IntegrityError", "UniqueViolation", "AttributeError"],
                err.__class__.__name__,
            )
            return
        raise Exception("Should have failed! Loaded %i records" % count)


class ClosedLoopTest(unittest.TestCase):
    """Test file -> BioSQL -> file."""

    @classmethod
    def setUpClass(cls):
        # NOTE - For speed I don't bother to create a new database each time,
        # simply a new unique namespace is used for each test.
        TESTDB = create_database()

    def test_NC_005816(self):
        """From GenBank file to BioSQL and back to a GenBank file, NC_005816."""
        with warnings.catch_warnings():
            # BiopythonWarning: order location operators are not fully supported
            warnings.simplefilter("ignore", BiopythonWarning)
            self.loop("GenBank/NC_005816.gb", "gb")

    def test_NC_000932(self):
        """From GenBank file to BioSQL and back to a GenBank file, NC_000932."""
        self.loop("GenBank/NC_000932.gb", "gb")

    def test_NT_019265(self):
        """From GenBank file to BioSQL and back to a GenBank file, NT_019265."""
        self.loop("GenBank/NT_019265.gb", "gb")

    def test_protein_refseq2(self):
        """From GenBank file to BioSQL and back to a GenBank file, protein_refseq2."""
        with warnings.catch_warnings():
            # BiopythonWarning: order location operators are not fully supported
            warnings.simplefilter("ignore", BiopythonWarning)
            self.loop("GenBank/protein_refseq2.gb", "gb")

    def test_no_ref(self):
        """From GenBank file to BioSQL and back to a GenBank file, noref."""
        self.loop("GenBank/noref.gb", "gb")

    def test_one_of(self):
        """From GenBank file to BioSQL and back to a GenBank file, one_of."""
        self.loop("GenBank/one_of.gb", "gb")

    def test_cor6_6(self):
        """From GenBank file to BioSQL and back to a GenBank file, cor6_6."""
        self.loop("GenBank/cor6_6.gb", "gb")

    def test_arab1(self):
        """From GenBank file to BioSQL and back to a GenBank file, arab1."""
        self.loop("GenBank/arab1.gb", "gb")

    def loop(self, filename, format):
        original_records = list(SeqIO.parse(filename, format))
        # now open a connection to load the database
        server = BioSeqDatabase.open_database(
            driver=DBDRIVER, user=DBUSER, passwd=DBPASSWD, host=DBHOST, db=TESTDB
        )
        db_name = "test_loop_%s" % filename  # new namespace!
        db = server.new_database(db_name)
        count = db.load(original_records)
        self.assertEqual(count, len(original_records))
        server.commit()
        # Now read them back...
        biosql_records = [db.lookup(name=rec.name) for rec in original_records]
        # And check they agree
        self.assertTrue(compare_records(original_records, biosql_records))
        # Now write to a handle...
        handle = StringIO()
        SeqIO.write(biosql_records, handle, "gb")
        # Now read them back...
        handle.seek(0)
        new_records = list(SeqIO.parse(handle, "gb"))
        # And check they still agree
        self.assertEqual(len(new_records), len(original_records))
        for old, new in zip(original_records, new_records):
            # TODO - remove this hack because we don't yet write these (yet):
            for key in ["comment", "references", "db_source"]:
                if key in old.annotations and key not in new.annotations:
                    del old.annotations[key]
            self.assertTrue(compare_record(old, new))
        # Done
        handle.close()
        server.close()


class TransferTest(unittest.TestCase):
    """Test file -> BioSQL, BioSQL -> BioSQL."""

    # NOTE - For speed I don't bother to create a new database each time,
    # simply a new unique namespace is used for each test.

    def setUp(self):
        TESTDB = create_database()

    def test_NC_005816(self):
        """From GenBank file to BioSQL, then again to a new namespace, NC_005816."""
        with warnings.catch_warnings():
            # BiopythonWarning: order location operators are not fully supported
            warnings.simplefilter("ignore", BiopythonWarning)
            self.trans("GenBank/NC_005816.gb", "gb")

    def test_NC_000932(self):
        """From GenBank file to BioSQL, then again to a new namespace, NC_000932."""
        self.trans("GenBank/NC_000932.gb", "gb")

    def test_NT_019265(self):
        """From GenBank file to BioSQL, then again to a new namespace, NT_019265."""
        self.trans("GenBank/NT_019265.gb", "gb")

    def test_protein_refseq2(self):
        """From GenBank file to BioSQL, then again to a new namespace, protein_refseq2."""
        with warnings.catch_warnings():
            # BiopythonWarning: order location operators are not fully supported
            warnings.simplefilter("ignore", BiopythonWarning)
            self.trans("GenBank/protein_refseq2.gb", "gb")

    def test_no_ref(self):
        """From GenBank file to BioSQL, then again to a new namespace, noref."""
        self.trans("GenBank/noref.gb", "gb")

    def test_one_of(self):
        """From GenBank file to BioSQL, then again to a new namespace, one_of."""
        self.trans("GenBank/one_of.gb", "gb")

    def test_cor6_6(self):
        """From GenBank file to BioSQL, then again to a new namespace, cor6_6."""
        self.trans("GenBank/cor6_6.gb", "gb")

    def test_arab1(self):
        """From GenBank file to BioSQL, then again to a new namespace, arab1."""
        self.trans("GenBank/arab1.gb", "gb")

    def trans(self, filename, format):
        original_records = list(SeqIO.parse(filename, format))
        # now open a connection to load the database
        server = BioSeqDatabase.open_database(
            driver=DBDRIVER, user=DBUSER, passwd=DBPASSWD, host=DBHOST, db=TESTDB
        )
        db_name = "test_trans1_%s" % filename  # new namespace!
        db = server.new_database(db_name)
        count = db.load(original_records)
        self.assertEqual(count, len(original_records))
        server.commit()
        # Now read them back...
        biosql_records = [db.lookup(name=rec.name) for rec in original_records]
        # And check they agree
        self.assertTrue(compare_records(original_records, biosql_records))
        # Now write to a second name space...
        db_name = "test_trans2_%s" % filename  # new namespace!
        db = server.new_database(db_name)
        count = db.load(biosql_records)
        self.assertEqual(count, len(original_records))
        # Now read them back again,
        biosql_records2 = [db.lookup(name=rec.name) for rec in original_records]
        # And check they also agree
        self.assertTrue(compare_records(original_records, biosql_records2))
        # Done
        server.close()

    def tearDown(self):
        destroy_database()


class InDepthLoadTest(unittest.TestCase):
    """Make sure we are loading and retreiving in a semi-lossless fashion."""

    def setUp(self):
        gb_file = os.path.join(os.getcwd(), "GenBank", "cor6_6.gb")
        load_database(gb_file)

        self.server = BioSeqDatabase.open_database(
            driver=DBDRIVER, user=DBUSER, passwd=DBPASSWD, host=DBHOST, db=TESTDB
        )
        self.db = self.server["biosql-test"]

    def tearDown(self):
        self.server.close()
        destroy_database()
        del self.db
        del self.server

    def test_transfer(self):
        """Make sure can load record into another namespace."""
        # Should be in database already...
        db_record = self.db.lookup(accession="X55053")
        # Make a new namespace
        db2 = self.server.new_database("biosql-test-alt")
        # Should be able to load this DBSeqRecord there...
        count = db2.load([db_record])
        self.assertEqual(count, 1)

    def test_reload(self):
        """Make sure can't reimport existing records."""
        gb_file = os.path.join(os.getcwd(), "GenBank", "cor6_6.gb")
        with open(gb_file) as gb_handle:
            record = next(SeqIO.parse(gb_handle, "gb"))
        # Should be in database already...
        db_record = self.db.lookup(accession="X55053")
        self.assertEqual(db_record.id, record.id)
        self.assertEqual(db_record.name, record.name)
        self.assertEqual(db_record.description, record.description)
        self.assertEqual(str(db_record.seq), str(record.seq))
        # Good... now try reloading it!
        try:
            count = self.db.load([record])
        except Exception as err:
            # Good!
            self.assertTrue(
                err.__class__.__name__
                in ["IntegrityError", "UniqueViolation", "AttributeError"],
                err.__class__.__name__,
            )
            return
        raise Exception("Should have failed! Loaded %i records" % count)

    def test_record_loading(self):
        """Make sure all records are correctly loaded."""
        test_record = self.db.lookup(accession="X55053")
        self.assertEqual(test_record.name, "ATCOR66M")
        self.assertEqual(test_record.id, "X55053.1")
        self.assertEqual(test_record.description, "A.thaliana cor6.6 mRNA")
        self.assertIsInstance(test_record.seq.alphabet, Alphabet.DNAAlphabet)
        self.assertEqual(str(test_record.seq[:10]), "AACAAAACAC")

        test_record = self.db.lookup(accession="X62281")
        self.assertEqual(test_record.name, "ATKIN2")
        self.assertEqual(test_record.id, "X62281.1")
        self.assertEqual(test_record.description, "A.thaliana kin2 gene")
        self.assertIsInstance(test_record.seq.alphabet, Alphabet.DNAAlphabet)
        self.assertEqual(str(test_record.seq[:10]), "ATTTGGCCTA")

    def test_seq_feature(self):
        """In depth check that SeqFeatures are transmitted through the db."""
        test_record = self.db.lookup(accession="AJ237582")
        features = test_record.features
        self.assertEqual(len(features), 7)

        # test single locations
        test_feature = features[0]
        self.assertEqual(test_feature.type, "source")
        self.assertEqual(str(test_feature.location), "[0:206](+)")
        self.assertEqual(len(test_feature.qualifiers), 3)
        self.assertEqual(test_feature.qualifiers["country"], ["Russia:Bashkortostan"])
        self.assertEqual(test_feature.qualifiers["organism"], ["Armoracia rusticana"])
        self.assertEqual(test_feature.qualifiers["db_xref"], ["taxon:3704"])

        # test split locations
        test_feature = features[4]
        self.assertEqual(test_feature.type, "CDS")
        self.assertEqual(str(test_feature.location), "join{[0:48](+), [142:206](+)}")
        self.assertEqual(len(test_feature.location.parts), 2)
        self.assertEqual(str(test_feature.location.parts[0]), "[0:48](+)")
        self.assertEqual(str(test_feature.location.parts[1]), "[142:206](+)")
        self.assertEqual(test_feature.location.operator, "join")
        self.assertEqual(len(test_feature.qualifiers), 6)
        self.assertEqual(test_feature.qualifiers["gene"], ["csp14"])
        self.assertEqual(test_feature.qualifiers["codon_start"], ["2"])
        self.assertEqual(test_feature.qualifiers["product"], ["cold shock protein"])
        self.assertEqual(test_feature.qualifiers["protein_id"], ["CAB39890.1"])
        self.assertEqual(test_feature.qualifiers["db_xref"], ["GI:4538893"])
        self.assertEqual(
            test_feature.qualifiers["translation"],
            ["DKAKDAAAAAGASAQQAGKNISDAAAGGVNFVKEKTG"],
        )

        # test passing strand information
        # XXX We should be testing complement as well
        test_record = self.db.lookup(accession="AJ237582")
        test_feature = test_record.features[4]  # DNA, no complement
        self.assertEqual(test_feature.strand, 1)
        for loc in test_feature.location.parts:
            self.assertEqual(loc.strand, 1)

        test_record = self.db.lookup(accession="X55053")
        test_feature = test_record.features[0]
        # mRNA, so really cDNA, so the strand should be 1 (not complemented)
        self.assertEqual(test_feature.strand, 1)


#####################################################################


class AutoSeqIOTests(unittest.TestCase):
    """Test SeqIO and BioSQL together."""

    server = None
    db = None

    @classmethod
    def setUpClass(cls):
        # Create and reuse on database for all tests in this class
        TESTDB = create_database()

    def setUp(self):
        """Connect to the database."""
        db_name = "biosql-test-seqio"
        server = BioSeqDatabase.open_database(
            driver=DBDRIVER, user=DBUSER, passwd=DBPASSWD, host=DBHOST, db=TESTDB
        )
        self.server = server
        if db_name not in server:
            self.db = server.new_database(db_name)
            server.commit()
        self.db = self.server[db_name]

    def tearDown(self):
        if self.db:
            del self.db
        if self.server:
            self.server.close()
            del self.server

    def check(self, t_format, t_filename, t_count=1):
        db = self.db

        iterator = SeqIO.parse(t_filename, t_format)
        count = db.load(iterator)
        assert count == t_count
        self.server.commit()

        iterator = SeqIO.parse(t_filename, t_format)
        for record in iterator:
            # print(" - %s, %s" % (checksum_summary(record), record.id))
            key = record.name
            # print(" - Retrieving by name/display_id '%s'," % key)
            db_rec = db.lookup(name=key)
            compare_record(record, db_rec)
            db_rec = db.lookup(display_id=key)
            compare_record(record, db_rec)

            key = record.id
            if key.count(".") == 1 and key.split(".")[1].isdigit():
                # print(" - Retrieving by version '%s'," % key)
                db_rec = db.lookup(version=key)
                compare_record(record, db_rec)

            if "accessions" in record.annotations:
                # Only expect FIRST accession to work!
                key = record.annotations["accessions"][0]
                assert key, "Blank accession in annotation %s" % repr(
                    record.annotations
                )
                if key != record.id:
                    # print(" - Retrieving by accession '%s'," % key)
                    db_rec = db.lookup(accession=key)
                    compare_record(record, db_rec)

            if "gi" in record.annotations:
                key = record.annotations["gi"]
                if key != record.id:
                    # print(" - Retrieving by GI '%s'," % key)
                    db_rec = db.lookup(primary_id=key)
                    compare_record(record, db_rec)

    def test_SeqIO_loading(self):
        self.check("fasta", "Fasta/lupine.nu")
        self.check("fasta", "Fasta/elderberry.nu")
        self.check("fasta", "Fasta/phlox.nu")
        self.check("fasta", "Fasta/centaurea.nu")
        self.check("fasta", "Fasta/wisteria.nu")
        self.check("fasta", "Fasta/sweetpea.nu")
        self.check("fasta", "Fasta/lavender.nu")
        self.check("fasta", "Fasta/aster.pro")
        self.check("fasta", "Fasta/loveliesbleeding.pro")
        self.check("fasta", "Fasta/rose.pro")
        self.check("fasta", "Fasta/rosemary.pro")
        self.check("fasta", "Fasta/f001")
        self.check("fasta", "Fasta/f002", 3)
        self.check("fasta", "Fasta/fa01", 2)
        self.check("fasta", "GFF/NC_001802.fna")
        self.check("fasta", "GFF/multi.fna", 3)
        self.check("fasta", "Registry/seqs.fasta", 2)
        self.check("swiss", "SwissProt/sp001")
        self.check("swiss", "SwissProt/sp002")
        self.check("swiss", "SwissProt/sp003")
        self.check("swiss", "SwissProt/P0A186.txt")
        self.check("swiss", "SwissProt/sp005")
        self.check("swiss", "SwissProt/sp006")
        self.check("swiss", "SwissProt/sp007")
        self.check("swiss", "SwissProt/sp008")
        self.check("swiss", "SwissProt/sp009")
        self.check("swiss", "SwissProt/sp010")
        self.check("swiss", "SwissProt/sp011")
        self.check("swiss", "SwissProt/sp012")
        self.check("swiss", "SwissProt/sp013")
        self.check("swiss", "SwissProt/P60137.txt")
        self.check("swiss", "SwissProt/sp015")
        self.check("swiss", "SwissProt/sp016")
        self.check("swiss", "Registry/EDD_RAT.dat")
        self.check("genbank", "GenBank/noref.gb")
        self.check("genbank", "GenBank/cor6_6.gb", 6)
        self.check("genbank", "GenBank/iro.gb")
        self.check("genbank", "GenBank/pri1.gb")
        self.check("genbank", "GenBank/arab1.gb")
        with warnings.catch_warnings():
            # BiopythonWarning: order location operators are not fully
            # supported
            warnings.simplefilter("ignore", BiopythonWarning)
            self.check("genbank", "GenBank/protein_refseq2.gb")
        self.check("genbank", "GenBank/extra_keywords.gb")
        self.check("genbank", "GenBank/one_of.gb")
        self.check("genbank", "GenBank/NT_019265.gb")
        self.check("genbank", "GenBank/origin_line.gb")
        self.check("genbank", "GenBank/blank_seq.gb")
        with warnings.catch_warnings():
            # BiopythonWarning: bond location operators are not fully supported
            warnings.simplefilter("ignore", BiopythonWarning)
            self.check("genbank", "GenBank/dbsource_wrap.gb")
            # BiopythonWarning: order location operators are not fully
            # supported
            self.check("genbank", "GenBank/NC_005816.gb")
        self.check("genbank", "GenBank/gbvrl1_start.seq", 3)
        self.check("genbank", "GFF/NC_001422.gbk")
        self.check("embl", "EMBL/TRBG361.embl")
        self.check("embl", "EMBL/DD231055_edited.embl")
        self.check("embl", "EMBL/SC10H5.embl")
        self.check("embl", "EMBL/U87107.embl")
        self.assertEqual(len(self.db), 66)


class SwissProtUnknownPositionTest(unittest.TestCase):
    """Handle SwissProt unknown position by setting value to null in database."""

    def setUp(self):
        # drop any old database and create a new one:
        TESTDB = create_database()
        # connect to new database:
        self.server = BioSeqDatabase.open_database(
            driver=DBDRIVER, user=DBUSER, passwd=DBPASSWD, host=DBHOST, db=TESTDB
        )
        # Create new namespace within new empty database:
        self.db = self.server.new_database("biosql-test")

    def tearDown(self):
        self.server.rollback()
        self.server.close()
        destroy_database()
        del self.db
        del self.server

    def test_ambiguous_location(self):
        """Loaded uniprot-xml with ambiguous location in BioSQL."""
        id = "P97881"
        seqiter = SeqIO.parse("SwissProt/%s.xml" % id, "uniprot-xml")
        self.assertEqual(self.db.load(seqiter), 1)

        dbrecord = self.db.lookup(primary_id=id)
        for feature in dbrecord.features:
            if feature.type == "signal peptide":
                self.assertIsInstance(feature.location.end, UnknownPosition)
            elif feature.type == "chain":
                self.assertIsInstance(feature.location.start, UnknownPosition)
            else:
                self.assertIsInstance(feature.location.start, ExactPosition)
