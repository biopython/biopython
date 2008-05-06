#!/usr/bin/env python
"""Tests for dealing with storage of biopython objects in a relational db.
"""
# standard library
import sys
import os

# PyUnit
import unittest

# local stuff
from Bio import MissingExternalDependencyError
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio import Alphabet
from Bio import GenBank

from BioSQL import BioSeqDatabase
from BioSQL import BioSeq

# This testing suite should try to detect whether a valid database
# installation exists on this computer.  Only run the tests if it
# does.
try :
    from setup_BioSQL import DBDRIVER, DBTYPE
    from setup_BioSQL import DBHOST, DBUSER, DBPASSWD, TESTDB
    from setup_BioSQL import DBSCHEMA, SQL_FILE
except NameError :
    message = "Check settings in Tests/setup_BioSQL.py "\
              "if you plan to use BioSQL."
    raise MissingExternalDependencyError(message)

try :
    server = BioSeqDatabase.open_database(driver = DBDRIVER,
                                          user = DBUSER, passwd = DBPASSWD,
                                          host = DBHOST)
    del server
except Exception, e :
    message = "Connection failed, check settings in Tests/setup_BioSQL.py "\
              "if you plan to use BioSQL: %s" % str(e)
    raise MissingExternalDependencyError(message)
  
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
    tests = [LoaderTest, ReadTest, SeqInterfaceTest, InDepthLoadTest]
    
    for test in tests:
        cur_suite = test_loader.loadTestsFromTestCase(test)
        test_suite.addTest(cur_suite)

    return test_suite

def create_database():
    """Create an empty BioSQL database."""
    # first open a connection to create the database
    server = BioSeqDatabase.open_database(driver = DBDRIVER,
                                          user = DBUSER, passwd = DBPASSWD,
                                          host = DBHOST)

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
        import time
        time.sleep(1)

        sql = r"DROP DATABASE " + TESTDB
        server.adaptor.cursor.execute(sql, ())
    except server.module.OperationalError: # the database doesn't exist
        pass
    except (server.module.IntegrityError,
            server.module.ProgrammingError), e: # ditto--perhaps
        if str(e).find('database "%s" does not exist' % TESTDB) > 0:
            pass
        else:
            raise
    # create a new database
    sql = r"CREATE DATABASE " + TESTDB
    server.adaptor.execute(sql, ())

    server.adaptor.conn.close()

    # now open a connection to load the database
    server = BioSeqDatabase.open_database(driver = DBDRIVER,
                                          user = DBUSER, passwd = DBPASSWD,
                                          host = DBHOST, db = TESTDB)
    server.load_database_sql(SQL_FILE)
    server.adaptor.conn.commit()
    server.adaptor.conn.close()

def load_database(gb_handle):
    """Load a GenBank file into a BioSQL database.
    
    This is useful for running tests against a newly created database.
    """

    create_database()
    # now open a connection to load the database
    db_name = "biosql-test"
    server = BioSeqDatabase.open_database(driver = DBDRIVER,
                                          user = DBUSER, passwd = DBPASSWD,
                                          host = DBHOST, db = TESTDB)
    db = server.new_database(db_name)
    
    # get the GenBank file we are going to put into it
    parser = GenBank.FeatureParser()
    iterator = GenBank.Iterator(gb_handle, parser)
    # finally put it in the database
    db.load(iterator)
    server.adaptor.conn.commit()
    server.adaptor.conn.close()

class ReadTest(unittest.TestCase):
    """Test reading a database from an already built database.
    """
    loaded_db = 0
    
    def setUp(self):
        """Connect to and load up the database.
        """
        gb_file = os.path.join(os.getcwd(), "GenBank", "cor6_6.gb")
        gb_handle = open(gb_file, "r")
        load_database(gb_handle)
        gb_handle.close()
            
        server = BioSeqDatabase.open_database(driver = DBDRIVER,
                                              user = DBUSER, 
                                              passwd = DBPASSWD,
                                              host = DBHOST, db = TESTDB)
            
        self.db = server["biosql-test"]

    def tearDown(self):
        self.db.adaptor.conn.close()
        del self.db

    def t_get_db_items(self):
        """Get a list of all items in the database.
        """
        items = self.db.values()

    def t_lookup_items(self):
        """Test retrieval of items using various ids.
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
        
        # primary id retrieval
        item = self.db.lookup(primary_id = "16353")
        try:
            item = self.db.lookup(primary_id = "Not Real")
            raise AssertionError("No problem on fake primary id retrieval")
        except IndexError:
            pass

class SeqInterfaceTest(unittest.TestCase):
    """Make sure the BioSQL objects implement the expected biopython interfaces
    """
    def setUp(self):
        """Load a database.
        """
        gb_file = os.path.join(os.getcwd(), "GenBank", "cor6_6.gb")
        gb_handle = open(gb_file, "r")
        load_database(gb_handle)
        gb_handle.close()

        server = BioSeqDatabase.open_database(driver = DBDRIVER,
                                              user = DBUSER, passwd = DBPASSWD,
                                              host = DBHOST, db = TESTDB)
        self.db = server["biosql-test"]
        self.item = self.db.lookup(accession = "X62281")

    def tearDown(self):
        self.db.adaptor.conn.close()
        del self.db
        del self.item
    
    def t_seq_record(self):
        """Make sure SeqRecords from BioSQL implement the right interface.
        """
        test_record = self.item
        assert isinstance(test_record.seq, BioSeq.DBSeq), \
          "Seq retrieval is not correct"
        assert test_record.id == "X62281.1", test_record.id
        assert test_record.name == "ATKIN2"
        assert test_record.description == "A.thaliana kin2 gene."

        annotations = test_record.annotations
        # XXX should do something with annotations once they are like
        # a dictionary
        for feature in test_record.features:
            assert isinstance(feature, SeqFeature)

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
        assert str(cds_feature.location) == "[103:579]", \
            str(cds_feature.location)
        for sub_feature in cds_feature.sub_features:
            assert sub_feature.type == "CDS"
            assert sub_feature.location_operator == "join"

        try :
            assert cds_feature.qualifiers["gene"] == ["kin2"]
            assert cds_feature.qualifiers["protein_id"] == ["CAA44171.1"]
            assert cds_feature.qualifiers["codon_start"] == ["1"]
        except KeyError :
            assert False, \
                   "Missing expected entries, have %s" % repr(cds_feature.qualifiers)
        
        assert "db_xref" in cds_feature.qualifiers, \
               cds_feature.qualifiers.keys()
        multi_ann = cds_feature.qualifiers["db_xref"]
        assert len(multi_ann) == 2
        assert "GI:16354" in multi_ann
        assert "SWISS-PROT:P31169" in multi_ann

class LoaderTest(unittest.TestCase):
    """Load a database from a GenBank file.
    """
    def setUp(self):
        # create TESTDB
        create_database()
        
        # load the database
        db_name = "biosql-test"
        server = BioSeqDatabase.open_database(driver = DBDRIVER,
                                              user = DBUSER, passwd = DBPASSWD,
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

    def tearDown(self):
        self.db.adaptor.conn.close()
        del self.db

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
        assert item_ids == ['AF297471.1', 'AJ237582.1', 'L31939.1', 'M81224.1',
                            'X55053.1', 'X62281.1'], item_ids

class InDepthLoadTest(unittest.TestCase):
    """Make sure we are loading and retreiving in a semi-lossless fashion.
    """
    def setUp(self):
        gb_file = os.path.join(os.getcwd(), "GenBank", "cor6_6.gb")
        gb_handle = open(gb_file, "r")
        load_database(gb_handle)
        gb_handle.close()

        server = BioSeqDatabase.open_database(driver = DBDRIVER,
                                              user = DBUSER, passwd = DBPASSWD,
                                              host = DBHOST, db = TESTDB)
        self.db = server["biosql-test"]

    def tearDown(self):
        self.db.adaptor.conn.close()
        del self.db

    def t_record_loading(self):
        """Make sure all records are correctly loaded.
        """
        test_record = self.db.lookup(accession = "X55053")
        assert test_record.name == "ATCOR66M"
        assert test_record.id == "X55053.1", test_record.id
        assert test_record.description == "A.thaliana cor6.6 mRNA."
        assert isinstance(test_record.seq.alphabet, Alphabet.DNAAlphabet)
        assert test_record.seq[:10].tostring() == 'AACAAAACAC'

        test_record = self.db.lookup(accession = "X62281")
        assert test_record.name == "ATKIN2"
        assert test_record.id == "X62281.1", test_record.id
        assert test_record.description == "A.thaliana kin2 gene."
        assert isinstance(test_record.seq.alphabet, Alphabet.DNAAlphabet)
        assert test_record.seq[:10].tostring() == 'ATTTGGCCTA'

    def t_seq_feature(self):
        """Indepth check that SeqFeatures are transmitted through the db.
        """
        test_record = self.db.lookup(accession = "AJ237582")
        features = test_record.features
        assert len(features) == 7
       
        # test single locations
        test_feature = features[0]
        assert test_feature.type == "source"
        assert str(test_feature.location) == "[0:206]"
        assert len(test_feature.qualifiers.keys()) == 3, \
               "Expected three keys, have %s" % repr(test_feature.qualifiers.keys())
        assert test_feature.qualifiers.has_key("country")
        assert test_feature.qualifiers["country"] == ["Russia:Bashkortostan"]
        assert test_feature.qualifiers.has_key("organism")
        assert test_feature.qualifiers["organism"] == ["Armoracia rusticana"]
        assert test_feature.qualifiers.has_key("db_xref")
        assert test_feature.qualifiers["db_xref"] == ["taxon:3704"], \
               "%s <> ['taxon:3704']" % test_feature.qualifiers["db_xref"]

        # test split locations
        test_feature = features[4]
        assert test_feature.type == "CDS", test_feature.type
        assert str(test_feature.location) == "[0:206]"
        assert len(test_feature.sub_features) == 2
        assert str(test_feature.sub_features[0].location) == "[0:48]"
        assert test_feature.sub_features[0].type == "CDS"
        assert test_feature.sub_features[0].location_operator == "join"
        assert str(test_feature.sub_features[1].location) == "[142:206]"
        assert test_feature.sub_features[1].type == "CDS"
        assert test_feature.sub_features[1].location_operator == "join"
        assert len(test_feature.qualifiers.keys()) == 6
        assert test_feature.qualifiers.has_key("product")
        assert test_feature.qualifiers["gene"] == ["csp14"]
        assert test_feature.qualifiers["codon_start"] == ["2"]
        assert test_feature.qualifiers["product"] == ["cold shock protein"]
        assert test_feature.qualifiers["protein_id"] == ["CAB39890.1"]
        assert test_feature.qualifiers["db_xref"] == ["GI:4538893"]
        assert test_feature.qualifiers["translation"] \
               == ["DKAKDAAAAAGASAQQAGKNISDAAAGGVNFVKEKTG"]

        # test passing strand information
        # XXX We should be testing complement as well
        test_record = self.db.lookup(accession = "AJ237582")
        test_feature = test_record.features[4] # DNA, no complement
        assert test_feature.strand == 1
        for sub_feature in test_feature.sub_features:
            assert sub_feature.strand == 1

        test_record = self.db.lookup(accession = "X55053")
        test_feature = test_record.features[0]
        # mRNA, so really cDNA, so the strand should be 1 (not complemented)
        assert test_feature.strand == 1, test_feature.strand

if __name__ == "__main__":
    sys.exit(run_tests(sys.argv))

