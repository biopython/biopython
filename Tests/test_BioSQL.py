#!/usr/bin/env python
"""Tests for dealing with storage of biopython objects in a relational db.

Currently these tests require a MySQL db loaded with the GenBank info
in GenBank/cor6_6.gb. This loading can be done with bioperl-db.
"""
# standard library
import sys
import os

# PyUnit
import unittest

# local stuff
import Bio
from BioSQL import BioSeqDatabase
from BioSQL import BioSeq

# Constants for the MySQL database
DBHOST = 'localhost'
DBUSER = 'root'
DBPASSWD = 'howareyou'
TESTDB = 'biosql'

def run_tests(argv):
    test_suite = testing_suite()
    runner = unittest.TextTestRunner(sys.stdout, verbosity = 2)
    runner.run(test_suite)

def testing_suite():
    """Generate the suite of tests.
    """
    test_suite = unittest.TestSuite()

    test_loader = unittest.TestLoader()
    test_loader.testMethodPrefix = 't_'
    tests = [LoadTest, SeqInterfaceTest]
    
    for test in tests:
        cur_suite = test_loader.loadTestsFromTestCase(test)
        test_suite.addTest(cur_suite)

    return test_suite

class LoadTest(unittest.TestCase):
    """Test loading a database from an already built database.

    XXX Once we have loading ability, this should use that instead
    of insisting on an existing database.
    """
    def setUp(self):
        """Connect to and load up the database.
        """
        server = BioSeqDatabase.open_database(user = DBUSER, passwd = DBPASSWD,
                                              host = DBHOST, db = TESTDB)
        self.db = server["biosql-test"]

    def t_get_db_items(self):
        """Get a list of all items in the database.
        """
        items = self.db.values()

    def t_lookup_items(self):
        """Test retrieval of itmes using various ids.
        """
        item = self.db.lookup(accession = "X62281")
        try:
            item = self.db.lookup(accession = "Not real")
            raise Assertionerror("No problem on fake id retrieval")
        except IndexError:
            pass
        item = self.db.lookup(display_id = "ATKIN2")
        try:
            item = self.db.lookup(display_id = "Not real")
            raise AssertionError("No problem on fake id retrieval")
        except IndexError:
            pass
        
        # primary id doesn't work right now
        try:
            item = self.db.lookup(primary_id = "16353")
            raise AssertionError("Need to write tests for primary_id fetch")
        except NotImplementedError:
            pass

class SeqInterfaceTest(unittest.TestCase):
    """Make sure the BioSQL objects implement the expected biopython interfaces
    """
    def setUp(self):
        """Load a database.

        XXX This is an already created database. We should actually
        build our own for testing when possible.
        """
        server = BioSeqDatabase.open_database(user = DBUSER, passwd = DBPASSWD,
                                              host = DBHOST, db = TESTDB)
        db = server["biosql-test"]
        self.item = db.lookup(accession = "X62281")
    
    def t_seq_record(self):
        """Make sure SeqRecords from BioSQL implement the right interface.
        """
        test_record = self.item
        assert isinstance(test_record.seq, BioSeq.DBSeq), \
          "Seq retrieval is not correct"
        assert test_record.id == "X62281"
        assert test_record.name == "ATKIN2"
        assert test_record.description == "" # should have a real description

        annotations = test_record.annotations
        # XXX should do something with annotations once they are like
        # a dictionary
        for feature in test_record.features:
            assert isinstance(feature, Bio.SeqFeature.SeqFeature)

    def t_seq(self):
        """Make sure Seqs from BioSQL implement the right interface.
        """
        test_seq = self.item.seq

    def t_seq_slicing(self):
        """Check that slices of sequences are retrieved properly.
        """
        test_seq = self.item.seq

    def t_seq_features(self):
        """Check SeqFeatures of a sequence.
        """
        test_features = self.item.features

if __name__ == "__main__":
    sys.exit(run_tests(sys.argv))
