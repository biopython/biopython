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
from Bio.Seq import Seq
from Bio import Alphabet
from Bio import GenBank

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
    tests = [ReadTest, SeqInterfaceTest, LoaderTest]
    tests = [LoaderTest]
    
    for test in tests:
        cur_suite = test_loader.loadTestsFromTestCase(test)
        test_suite.addTest(cur_suite)

    return test_suite

class ReadTest(unittest.TestCase):
    """Test reading a database from an already built database.

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
        alphabet = test_seq.alphabet
        assert isinstance(alphabet, Alphabet.Alphabet)

        data = test_seq.data
        assert type(data) == type("")
    
        string_rep = test_seq.tostring()
        assert type(string_rep) == type("")
    
        assert len(test_seq) == 880, len(test_seq)

    def t_seq_slicing(self):
        """Check that slices of sequences are retrieved properly.
        """
        test_seq = self.item.seq
        new_seq = test_seq[:10]
        assert isinstance(new_seq, BioSeq.DBSeq)

        # simple slicing
        assert test_seq[:5].tostring() == 'ATTTG'
        assert test_seq[0:5].tostring() == 'ATTTG'
        assert test_seq[2:3].tostring() == 'T'
        assert test_seq[2:4].tostring() == 'TT'
        assert test_seq[870:].tostring() == 'TTGAATTATA'

        # getting more fancy
        assert test_seq[-1] == 'A'
        assert test_seq[1] == 'T'
        assert test_seq[-10:][5:].tostring() == "TTATA"

    def t_seq_features(self):
        """Check SeqFeatures of a sequence.
        """
        test_features = self.item.features
        cds_feature = test_features[6]
        assert cds_feature.type == "CDS", cds_feature.type
        assert str(cds_feature.location) == "(103..579)"
        for sub_feature in cds_feature.sub_features:
            assert sub_feature.type == "CDS"
            assert sub_feature.location_operator == "join"
       
        ann = cds_feature.qualifiers["gene"]
        assert ann == ["kin2"]
        multi_ann = cds_feature.qualifiers["db_xref"]
        assert len(multi_ann) == 2
        assert "GI:16354" in multi_ann
        assert "SWISS-PROT:P31169" in multi_ann

class LoaderTest(unittest.TestCase):
    """Load a database from a GenBank file.
    """
    def setUp(self):
        # load the database
        db_name = "biosql-loadertest"
        server = BioSeqDatabase.open_database(user = DBUSER, passwd = DBPASSWD,
                                              host = DBHOST, db = TESTDB)
        
        # remove the database if it already exists
        try:
            server[db_name]
            server.remove_database(db_name)
        except KeyError:
            pass
        
        self.db = server.new_database(db_name)

        # get the GenBank file we are going to put into it
        input_file = os.path.join(os.getcwd(), "GenBank", "cor6_6.gb")
        handle = open(input_file, "r")
        parser = GenBank.FeatureParser()
        self.iterator = GenBank.Iterator(handle, parser)

    def t_load_database(self):
        """Load SeqRecord objects into a BioSQL database.
        """
        self.db.load(self.iterator)

        # do some simple tests to make sure we actually loaded the right
        # thing. More advanced tests in a different module.
        items = self.db.values()
        assert len(items) == 6
        item_names = []
        item_ids = []
        for item in items:
            item_names.append(item.name)
            item_ids.append(item.id)
        item_names.sort()
        item_ids.sort()
        assert item_names == ['AF297471', 'ARU237582', 'ATCOR66M', 
                              'ATKIN2', 'BNAKINI', 'BRRBIF72']
        assert item_ids == ['AF297471', 'AJ237582', 'L31939', 'M81224', 
                            'X55053', 'X62281']

if __name__ == "__main__":
    sys.exit(run_tests(sys.argv))
