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

from common_BioSQL import create_database, destroy_database, check_config

from seq_tests_common import compare_record, compare_records

import requires_internet

if __name__ == "__main__":
    raise RuntimeError("Call this via test_BioSQL_*online.py not directly")

# Sharing these with test_BioSQL_XXX_online.py files which import this file:
# DBDRIVER, DBTYPE, DBHOST, DBUSER, DBPASSWD, TESTDB, DBSCHEMA, SQL_FILE, SYSTEM
SYSTEM = platform.system()


def share_config(dbdriver, dbtype, dbhost, dbuser, dbpasswd, testdb):
    """Make sure we can access the DB settings from this file."""
    global DBDRIVER, DBTYPE, DBHOST, DBUSER, DBPASSWD, TESTDB, DBSCHEMA
    global SYSTEM, SQL_FILE
    DBDRIVER = dbdriver
    DBTYPE = dbtype
    DBHOST = dbhost
    DBUSER = dbuser
    DBPASSWD = dbpasswd
    TESTDB = testdb


class TaxonomyTest(unittest.TestCase):
    """Test proper insertion and retrieval of taxonomy data."""
    def setUp(self):
        global DBDRIVER, DBTYPE, DBHOST, DBUSER, DBPASSWD, TESTDB, DBSCHEMA
        global SYSTEM, SQL_FILE

        Entrez.email = "biopython-dev@biopython.org"
        # create TESTDB
        TESTDB = create_database()

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
