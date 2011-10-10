#!/usr/bin/env python
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for dealing with storage of biopython objects in a relational db.
"""
# standard library
import os
import unittest
from StringIO import StringIO

# local stuff
from Bio import MissingExternalDependencyError
from Bio.Seq import Seq, MutableSeq
from Bio.SeqFeature import SeqFeature
from Bio import Alphabet
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from BioSQL import BioSeqDatabase
from BioSQL import BioSeq

# This testing suite should try to detect whether a valid database
# installation exists on this computer.  Only run the tests if it
# does.
try:
    from setup_BioSQL import DBDRIVER, DBTYPE
    from setup_BioSQL import DBHOST, DBUSER, DBPASSWD, TESTDB
    from setup_BioSQL import DBSCHEMA, SQL_FILE
except (NameError, ImportError):
    message = "Check settings in Tests/setup_BioSQL.py "\
              "if you plan to use BioSQL."
    raise MissingExternalDependencyError(message)

try:
    if DBDRIVER in ["sqlite3"]:
        server = BioSeqDatabase.open_database(driver = DBDRIVER, db = TESTDB)
    else:
        server = BioSeqDatabase.open_database(driver = DBDRIVER,
                                              user = DBUSER, passwd = DBPASSWD,
                                              host = DBHOST)
    del server
except Exception, e:
    message = "Connection failed, check settings in Tests/setup_BioSQL.py "\
              "if you plan to use BioSQL: %s" % str(e)
    raise MissingExternalDependencyError(message)

from seq_tests_common import compare_record, compare_records

def _do_db_create():
    """Do the actual work of database creation. Relevant for MySQL and PostgreSQL
    """
    # first open a connection to create the database
    server = BioSeqDatabase.open_database(driver = DBDRIVER,
                                          user = DBUSER, passwd = DBPASSWD,
                                          host = DBHOST)

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
        import time
        time.sleep(1)

        sql = r"DROP DATABASE " + TESTDB
        server.adaptor.cursor.execute(sql, ())
    except (server.module.OperationalError,
            server.module.DatabaseError), e: # the database doesn't exist
        pass
    except (server.module.IntegrityError,
            server.module.ProgrammingError), e: # ditto--perhaps
        if str(e).find('database "%s" does not exist' % TESTDB) == -1:
            raise
    # create a new database
    sql = r"CREATE DATABASE " + TESTDB
    server.adaptor.execute(sql, ())
    server.close()

def create_database():
    """Delete any existing BioSQL test database, then (re)create an empty BioSQL database."""
    if DBDRIVER in ["sqlite3"]: 
        if os.path.exists(TESTDB):
            os.remove(TESTDB)
    else:
        _do_db_create()

    # now open a connection to load the database
    server = BioSeqDatabase.open_database(driver = DBDRIVER,
                                          user = DBUSER, passwd = DBPASSWD,
                                          host = DBHOST, db = TESTDB)
    server.load_database_sql(SQL_FILE)
    server.commit()
    server.close()

def load_database(gb_handle):
    """Load a GenBank file into a new BioSQL database.
    
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
    iterator = SeqIO.parse(gb_handle, "gb")
    # finally put it in the database
    count = db.load(iterator)
    server.commit()
    server.close()
    return count

class ReadTest(unittest.TestCase):
    """Test reading a database from an already built database.
    """
    loaded_db = 0
    
    def setUp(self):
        """Connect to and load up the database.
        """
        gb_file = "GenBank/cor6_6.gb"
        gb_handle = open(gb_file, "r")
        load_database(gb_handle)
        gb_handle.close()
            
        self.server = BioSeqDatabase.open_database(driver = DBDRIVER,
                                              user = DBUSER, 
                                              passwd = DBPASSWD,
                                              host = DBHOST, db = TESTDB)
            
        self.db = self.server["biosql-test"]

    def tearDown(self):
        self.server.close()
        del self.db
        del self.server

    def test_server(self):
        """Check BioSeqDatabase methods"""
        server = self.server
        self.assertTrue("biosql-test" in server)
        self.assertEqual(1, len(server))
        self.assertEqual(["biosql-test"], server.keys())
        #Check we can delete the namespace...
        del server["biosql-test"]
        self.assertEqual(0, len(server))
        try:
            del server["non-existant-name"]
            assert False, "Should have raised KeyError"
        except KeyError:
            pass

    def test_get_db_items(self):
        """Check list, keys, length etc"""
        db = self.db
        items = db.values()
        keys = db.keys()
        l = len(items)
        self.assertEqual(l, len(db))
        self.assertEqual(l, len(list(db.iteritems())))
        self.assertEqual(l, len(list(db.iterkeys())))
        self.assertEqual(l, len(list(db.itervalues())))
        for (k1,r1), (k2,r2) in zip(zip(keys, items), db.iteritems()):
            self.assertEqual(k1, k2)
            self.assertEqual(r1.id, r2.id)
        for k in keys:
            del db[k]
        self.assertEqual(0, len(db))
        try:
            del db["non-existant-name"]
            assert False, "Should have raised KeyError"
        except KeyError:
            pass

    def test_lookup_items(self):
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

        self.server = BioSeqDatabase.open_database(driver = DBDRIVER,
                                              user = DBUSER, passwd = DBPASSWD,
                                              host = DBHOST, db = TESTDB)
        self.db = self.server["biosql-test"]
        self.item = self.db.lookup(accession = "X62281")

    def tearDown(self):
        self.server.close()
        del self.db
        del self.item
        del self.server
    
    def test_seq_record(self):
        """Make sure SeqRecords from BioSQL implement the right interface.
        """
        test_record = self.item
        self.assertTrue(isinstance(test_record.seq, BioSeq.DBSeq))
        self.assertEqual(test_record.id, "X62281.1", test_record.id)
        self.assertEqual(test_record.name, "ATKIN2")
        self.assertEqual(test_record.description, "A.thaliana kin2 gene.")
        annotations = test_record.annotations
        # XXX should do something with annotations once they are like
        # a dictionary
        for feature in test_record.features:
            self.assertTrue(isinstance(feature, SeqFeature))
        s = str(test_record) #shouldn't cause any errors!

    def test_seq(self):
        """Make sure Seqs from BioSQL implement the right interface.
        """
        test_seq = self.item.seq
        alphabet = test_seq.alphabet
        self.assertTrue(isinstance(alphabet, Alphabet.Alphabet))
        data = test_seq.data
        self.assertEqual(type(data), type(""))
        string_rep = test_seq.tostring()
        self.assertEqual(string_rep, str(test_seq)) #check __str__ too
        self.assertEqual(type(string_rep), type(""))
        self.assertEqual(len(test_seq), 880)
        
    def test_convert(self):
        """Check can turn a DBSeq object into a Seq or MutableSeq."""
        test_seq = self.item.seq

        other = test_seq.toseq()
        self.assertEqual(str(test_seq), str(other))
        self.assertEqual(test_seq.alphabet, other.alphabet)
        self.assertTrue(isinstance(other, Seq))

        other = test_seq.tomutable()
        self.assertEqual(str(test_seq), str(other))
        self.assertEqual(test_seq.alphabet, other.alphabet)
        self.assertTrue(isinstance(other, MutableSeq))

    def test_addition(self):
        """Check can add DBSeq objects together."""
        test_seq = self.item.seq
        for other in [Seq("ACGT",test_seq.alphabet),
                      MutableSeq("ACGT",test_seq.alphabet),
                      "ACGT",
                      test_seq]:
            test = test_seq + other
            self.assertEqual(str(test), str(test_seq) + str(other))
            self.assertTrue(isinstance(test, Seq))
            test = other + test_seq
            self.assertEqual(str(test), str(other) + str(test_seq))

    def test_seq_slicing(self):
        """Check that slices of sequences are retrieved properly.
        """
        test_seq = self.item.seq
        new_seq = test_seq[:10]
        self.assertTrue(isinstance(new_seq, BioSeq.DBSeq))
        # simple slicing
        self.assertEqual(test_seq[:5].tostring(), 'ATTTG')
        self.assertEqual(test_seq[0:5].tostring(), 'ATTTG')
        self.assertEqual(test_seq[2:3].tostring(), 'T')
        self.assertEqual(test_seq[2:4].tostring(), 'TT')
        self.assertEqual(test_seq[870:].tostring(), 'TTGAATTATA')
        # getting more fancy
        self.assertEqual(test_seq[-1], 'A')
        self.assertEqual(test_seq[1], 'T')
        self.assertEqual(test_seq[-10:][5:].tostring(), "TTATA")
        self.assertEqual(str(test_seq[-10:][5:]), "TTATA")

    def test_seq_features(self):
        """Check SeqFeatures of a sequence.
        """
        test_features = self.item.features
        cds_feature = test_features[6]
        self.assertEqual(cds_feature.type, "CDS")
        self.assertEqual(str(cds_feature.location), "[103:579](+)")
        for sub_feature in cds_feature.sub_features:
            self.assertEqual(sub_feature.type, "CDS")
            self.assertEqual(sub_feature.location_operator, "join")

        try:
            self.assertEqual(cds_feature.qualifiers["gene"], ["kin2"])
            self.assertEqual(cds_feature.qualifiers["protein_id"], ["CAA44171.1"])
            self.assertEqual(cds_feature.qualifiers["codon_start"], ["1"])
        except KeyError:
            raise KeyError("Missing expected entries, have %s" \
                           % repr(cds_feature.qualifiers))
        
        self.assertTrue("db_xref" in cds_feature.qualifiers)
        multi_ann = cds_feature.qualifiers["db_xref"]
        self.assertEqual(len(multi_ann), 2)
        self.assertTrue("GI:16354" in multi_ann)
        self.assertTrue("SWISS-PROT:P31169" in multi_ann)

class LoaderTest(unittest.TestCase):
    """Load a database from a GenBank file.
    """
    def setUp(self):
        # create TESTDB
        create_database()
        
        # load the database
        db_name = "biosql-test"
        self.server = BioSeqDatabase.open_database(driver = DBDRIVER,
                                              user = DBUSER, passwd = DBPASSWD,
                                              host = DBHOST, db = TESTDB)
        
        # remove the database if it already exists
        try:
            self.server[db_name]
            self.server.remove_database(db_name)
        except KeyError:
            pass
        
        self.db = self.server.new_database(db_name)

        # get the GenBank file we are going to put into it
        input_file = os.path.join(os.getcwd(), "GenBank", "cor6_6.gb")
        handle = open(input_file, "r")
        self.iterator = SeqIO.parse(handle, "gb")

    def tearDown(self):
        self.server.close()
        del self.db
        del self.server

    def test_load_database(self):
        """Load SeqRecord objects into a BioSQL database.
        """
        self.db.load(self.iterator)

        # do some simple tests to make sure we actually loaded the right
        # thing. More advanced tests in a different module.
        items = self.db.values()
        self.assertEqual(len(items), 6)
        self.assertEqual(len(self.db), 6)
        item_names = []
        item_ids = []
        for item in items:
            item_names.append(item.name)
            item_ids.append(item.id)
        item_names.sort()
        item_ids.sort()
        self.assertEqual(item_names, ['AF297471', 'ARU237582', 'ATCOR66M',
                                      'ATKIN2', 'BNAKINI', 'BRRBIF72'])
        self.assertEqual(item_ids, ['AF297471.1', 'AJ237582.1', 'L31939.1',
                                    'M81224.1', 'X55053.1', 'X62281.1'])

class DupLoadTest(unittest.TestCase):
    """Check a few duplicate conditions fail."""
    def setUp(self):
        #drop any old database and create a new one:
        create_database()
        #connect to new database:
        self.server = BioSeqDatabase.open_database(driver = DBDRIVER,
                                              user = DBUSER, passwd = DBPASSWD,
                                              host = DBHOST, db = TESTDB)
        #Create new namespace within new empty database:
        self.db = self.server.new_database("biosql-test")

    def tearDown(self):
        self.server.rollback()
        self.server.close()
        del self.db
        del self.server

    def test_duplicate_load(self):
        """Make sure can't import a single record twice (in one go)."""
        record = SeqRecord(Seq("ATGCTATGACTAT", Alphabet.generic_dna),id="Test1")
        try:
            count = self.db.load([record,record])
        except Exception, err:
            #Good!
            #Note we don't do a specific exception handler because the
            #exception class will depend on which DB back end is in use.            
            self.assertTrue(err.__class__.__name__ in ["IntegrityError",
                                                       "OperationalError"],
                            err.__class__.__name__)
            return
        raise Exception("Should have failed! Loaded %i records" % count)

    def test_duplicate_load2(self):
        """Make sure can't import a single record twice (in steps)."""
        record = SeqRecord(Seq("ATGCTATGACTAT", Alphabet.generic_dna),id="Test2")
        count = self.db.load([record])
        self.assertEqual(count,1)
        try:
            count = self.db.load([record])
        except Exception, err:
            #Good!
            self.assertEqual("IntegrityError", err.__class__.__name__)
            return
        raise Exception("Should have failed! Loaded %i records" % count)

    def test_duplicate_id_load(self):
        """Make sure can't import records with same ID (in one go)."""
        record1 = SeqRecord(Seq("ATGCTATGACTAT", Alphabet.generic_dna),id="TestA")
        record2 = SeqRecord(Seq("GGGATGCGACTAT", Alphabet.generic_dna),id="TestA")
        try:
            count = self.db.load([record1,record2])
        except Exception, err:
            #Good!
            self.assertEqual("IntegrityError", err.__class__.__name__)
            return
        raise Exception("Should have failed! Loaded %i records" % count)

class ClosedLoopTest(unittest.TestCase):
    """Test file -> BioSQL -> file."""
    #NOTE - For speed I don't bother to create a new database each time,
    #simple a new unique namespace is used for each test.
    
    def test_NC_005816(self):
        """GenBank file to BioSQL and back to a GenBank file, NC_005816."""
        self.loop("GenBank/NC_005816.gb", "gb")

    def test_NC_000932(self):
        """GenBank file to BioSQL and back to a GenBank file, NC_000932."""
        self.loop("GenBank/NC_000932.gb", "gb")

    def test_NT_019265(self):
        """GenBank file to BioSQL and back to a GenBank file, NT_019265."""
        self.loop("GenBank/NT_019265.gb", "gb")

    def test_protein_refseq2(self):
        """GenBank file to BioSQL and back to a GenBank file, protein_refseq2."""
        self.loop("GenBank/protein_refseq2.gb", "gb")

    def test_no_ref(self):
        """GenBank file to BioSQL and back to a GenBank file, noref."""
        self.loop("GenBank/noref.gb", "gb")

    def test_one_of(self):
        """GenBank file to BioSQL and back to a GenBank file, one_of."""
        self.loop("GenBank/one_of.gb", "gb")

    def test_cor6_6(self):
        """GenBank file to BioSQL and back to a GenBank file, cor6_6."""
        self.loop("GenBank/cor6_6.gb", "gb")

    def test_arab1(self):
        """GenBank file to BioSQL and back to a GenBank file, arab1."""
        self.loop("GenBank/arab1.gb", "gb")

    def loop(self, filename, format):
        original_records = list(SeqIO.parse(open(filename, "rU"), format))
        # now open a connection to load the database
        server = BioSeqDatabase.open_database(driver = DBDRIVER,
                                              user = DBUSER, passwd = DBPASSWD,
                                              host = DBHOST, db = TESTDB)
        db_name = "test_loop_%s" % filename #new namespace!
        db = server.new_database(db_name)
        count = db.load(original_records)
        self.assertEqual(count, len(original_records))
        server.commit()
        #Now read them back...
        biosql_records = [db.lookup(name=rec.name) \
                          for rec in original_records]
        #And check they agree
        self.assertTrue(compare_records(original_records, biosql_records))
        #Now write to a handle...
        handle = StringIO()
        SeqIO.write(biosql_records, handle, "gb")
        #Now read them back...
        handle.seek(0)
        new_records = list(SeqIO.parse(handle, "gb"))
        #And check they still agree
        self.assertEqual(len(new_records), len(original_records))
        for old, new in zip(original_records, new_records):
            #TODO - remove this hack because we don't yet write these (yet):
            for key in ["comment", "references", "db_source"]:
                if key in old.annotations and key not in new.annotations:
                    del old.annotations[key]
            self.assertTrue(compare_record(old, new))
        #Done
        server.close()

class TransferTest(unittest.TestCase):
    """Test file -> BioSQL, BioSQL -> BioSQL."""
    #NOTE - For speed I don't bother to create a new database each time,
    #simple a new unique namespace is used for each test.
    
    def test_NC_005816(self):
        """GenBank file to BioSQL, then again to a new namespace, NC_005816."""
        self.trans("GenBank/NC_005816.gb", "gb")

    def test_NC_000932(self):
        """GenBank file to BioSQL, then again to a new namespace, NC_000932."""
        self.trans("GenBank/NC_000932.gb", "gb")

    def test_NT_019265(self):
        """GenBank file to BioSQL, then again to a new namespace, NT_019265."""
        self.trans("GenBank/NT_019265.gb", "gb")

    def test_protein_refseq2(self):
        """GenBank file to BioSQL, then again to a new namespace, protein_refseq2."""
        self.trans("GenBank/protein_refseq2.gb", "gb")

    def test_no_ref(self):
        """GenBank file to BioSQL, then again to a new namespace, noref."""
        self.trans("GenBank/noref.gb", "gb")

    def test_one_of(self):
        """GenBank file to BioSQL, then again to a new namespace, one_of."""
        self.trans("GenBank/one_of.gb", "gb")

    def test_cor6_6(self):
        """GenBank file to BioSQL, then again to a new namespace, cor6_6."""
        self.trans("GenBank/cor6_6.gb", "gb")

    def test_arab1(self):
        """GenBank file to BioSQL, then again to a new namespace, arab1."""
        self.trans("GenBank/arab1.gb", "gb")

    def trans(self, filename, format):
        original_records = list(SeqIO.parse(open(filename, "rU"), format))
        # now open a connection to load the database
        server = BioSeqDatabase.open_database(driver = DBDRIVER,
                                              user = DBUSER, passwd = DBPASSWD,
                                              host = DBHOST, db = TESTDB)
        db_name = "test_trans1_%s" % filename #new namespace!
        db = server.new_database(db_name)
        count = db.load(original_records)
        self.assertEqual(count, len(original_records))
        server.commit()
        #Now read them back...
        biosql_records = [db.lookup(name=rec.name) \
                          for rec in original_records]
        #And check they agree
        self.assertTrue(compare_records(original_records, biosql_records))
        #Now write to a second name space...
        db_name = "test_trans2_%s" % filename #new namespace!
        db = server.new_database(db_name)
        count = db.load(biosql_records)
        self.assertEqual(count, len(original_records))
        #Now read them back again,
        biosql_records2 = [db.lookup(name=rec.name) \
                          for rec in original_records]
        #And check they also agree
        self.assertTrue(compare_records(original_records, biosql_records2))
        #Done
        server.close()


class InDepthLoadTest(unittest.TestCase):
    """Make sure we are loading and retreiving in a semi-lossless fashion.
    """
    def setUp(self):
        gb_file = os.path.join(os.getcwd(), "GenBank", "cor6_6.gb")
        gb_handle = open(gb_file, "r")
        load_database(gb_handle)
        gb_handle.close()

        self.server = BioSeqDatabase.open_database(driver = DBDRIVER,
                                              user = DBUSER, passwd = DBPASSWD,
                                              host = DBHOST, db = TESTDB)
        self.db = self.server["biosql-test"]

    def tearDown(self):
        self.server.close()
        del self.db
        del self.server

    def test_transfer(self):
        """Make sure can load record into another namespace."""
        #Should be in database already...
        db_record = self.db.lookup(accession = "X55053")
        #Make a new namespace
        db2 = self.server.new_database("biosql-test-alt")
        #Should be able to load this DBSeqRecord there...
        count = db2.load([db_record])
        self.assertEqual(count,1)

    def test_reload(self):
        """Make sure can't reimport existing records."""
        gb_file = os.path.join(os.getcwd(), "GenBank", "cor6_6.gb")
        gb_handle = open(gb_file, "r")
        record = SeqIO.parse(gb_handle, "gb").next()
        gb_handle.close()
        #Should be in database already...
        db_record = self.db.lookup(accession = "X55053")
        self.assertEqual(db_record.id, record.id)
        self.assertEqual(db_record.name, record.name)
        self.assertEqual(db_record.description, record.description)
        self.assertEqual(str(db_record.seq), str(record.seq))
        #Good... now try reloading it!
        try:
            count = self.db.load([record])
        except Exception, err:
            #Good!
            self.assertEqual("IntegrityError", err.__class__.__name__)
            return
        raise Exception("Should have failed! Loaded %i records" % count)

    def test_record_loading(self):
        """Make sure all records are correctly loaded.
        """
        test_record = self.db.lookup(accession = "X55053")
        self.assertEqual(test_record.name, "ATCOR66M")
        self.assertEqual(test_record.id, "X55053.1")
        self.assertEqual(test_record.description, "A.thaliana cor6.6 mRNA.")
        self.assertTrue(isinstance(test_record.seq.alphabet, Alphabet.DNAAlphabet))
        self.assertEqual(test_record.seq[:10].tostring(), 'AACAAAACAC')

        test_record = self.db.lookup(accession = "X62281")
        self.assertEqual(test_record.name, "ATKIN2")
        self.assertEqual(test_record.id, "X62281.1")
        self.assertEqual(test_record.description, "A.thaliana kin2 gene.")
        self.assertTrue(isinstance(test_record.seq.alphabet, Alphabet.DNAAlphabet))
        self.assertEqual(test_record.seq[:10].tostring(), 'ATTTGGCCTA')

    def test_seq_feature(self):
        """Indepth check that SeqFeatures are transmitted through the db.
        """
        test_record = self.db.lookup(accession = "AJ237582")
        features = test_record.features
        self.assertEqual(len(features), 7)
       
        # test single locations
        test_feature = features[0]
        self.assertEqual(test_feature.type, "source")
        self.assertEqual(str(test_feature.location), "[0:206](+)")
        self.assertEqual(len(test_feature.qualifiers.keys()), 3)
        self.assertEqual(test_feature.qualifiers["country"], ["Russia:Bashkortostan"])
        self.assertEqual(test_feature.qualifiers["organism"], ["Armoracia rusticana"])
        self.assertEqual(test_feature.qualifiers["db_xref"], ["taxon:3704"])

        # test split locations
        test_feature = features[4]
        self.assertEqual(test_feature.type, "CDS")
        self.assertEqual(str(test_feature.location), "[0:206](+)")
        self.assertEqual(len(test_feature.sub_features), 2)
        self.assertEqual(str(test_feature.sub_features[0].location), "[0:48](+)")
        self.assertEqual(test_feature.sub_features[0].type, "CDS")
        self.assertEqual(test_feature.sub_features[0].location_operator, "join")
        self.assertEqual(str(test_feature.sub_features[1].location), "[142:206](+)")
        self.assertEqual(test_feature.sub_features[1].type, "CDS")
        self.assertEqual(test_feature.sub_features[1].location_operator, "join")
        self.assertEqual(len(test_feature.qualifiers.keys()), 6)
        self.assertEqual(test_feature.qualifiers["gene"], ["csp14"])
        self.assertEqual(test_feature.qualifiers["codon_start"], ["2"])
        self.assertEqual(test_feature.qualifiers["product"],
                         ["cold shock protein"])
        self.assertEqual(test_feature.qualifiers["protein_id"], ["CAB39890.1"])
        self.assertEqual(test_feature.qualifiers["db_xref"], ["GI:4538893"])
        self.assertEqual(test_feature.qualifiers["translation"],
                         ["DKAKDAAAAAGASAQQAGKNISDAAAGGVNFVKEKTG"])

        # test passing strand information
        # XXX We should be testing complement as well
        test_record = self.db.lookup(accession = "AJ237582")
        test_feature = test_record.features[4] # DNA, no complement
        self.assertEqual(test_feature.strand, 1)
        for sub_feature in test_feature.sub_features:
            self.assertEqual(sub_feature.strand, 1)

        test_record = self.db.lookup(accession = "X55053")
        test_feature = test_record.features[0]
        # mRNA, so really cDNA, so the strand should be 1 (not complemented)
        self.assertEqual(test_feature.strand, 1)

#Some of the unit tests don't create their own database,
#so just in case there is no database already:
create_database()

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
