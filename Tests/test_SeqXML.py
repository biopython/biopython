# Copyright 2010 by Thomas Schmitt.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
import unittest

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from StringIO import StringIO

test_files = {
    "dna" : ["SeqXML/dna_example.xml",4,None],
    "rna" : ["SeqXML/rna_example.xml",5,None],
    "protein" : ["SeqXML/protein_example.xml",5,None],
}

corrupt_files = {"corrupt1" : ["SeqXML/corrupt_example1.xml",None],
                 "corrupt2" : ["SeqXML/corrupt_example2.xml",None],}

def open_test_files():
    for _file in test_files:
        test_files[_file][2] = open(test_files[_file][0],"r")
    for _file in corrupt_files:
        corrupt_files[_file][1] = open(corrupt_files[_file][0],"r")

def close_test_files():
    for _file in test_files:
        test_files[_file][2].close()
    for _file in corrupt_files:
        corrupt_files[_file][1].close()
        
def assert_equal_records(testCase,record_a,record_b):
    
    testCase.assertEqual(record_a.id,record_b.id)
    testCase.assertEqual(record_a.name,record_b.name)
    testCase.assertEqual(record_a.description,record_b.description)
    testCase.assertEqual(record_a.seq.tostring(),record_b.seq.tostring())
    testCase.assertEqual(record_a.dbxrefs,record_b.dbxrefs)
    testCase.assertEqual(record_a.annotations,record_a.annotations)
    
    

class TestSimpleRead(unittest.TestCase):
    
    def setUp(self):
        open_test_files()

    def tearDown(self):
        close_test_files()
        
    def test_check_SeqIO(self):
        """Files readable using parser via SeqIO."""
        
        for _file in test_files:
            records = list(SeqIO.parse(test_files[_file][2],"seqxml"))
            
            self.assertEqual(len(records),test_files[_file][1])
            
class TestDetalidRead(unittest.TestCase):
    
    records = {}
    
    def setUp(self):
        open_test_files()
        
        for _file in test_files:
            self.records[_file] = list(SeqIO.parse(test_files[_file][2],"seqxml"))

    def tearDown(self):
        close_test_files()
        
        
    def test_special_characters_desc(self):
        """Read special XML characters in description."""
        
        self.assertEqual(self.records["dna"][2].description, u'some special characters in the description\n<tag> "quoted string"')
        
    def test_unicode_characters_desc(self):
        """Test special unicode characters in the description."""
        self.assertEqual(self.records["rna"][2].description, u"\u00E5\u00C5\u00FC\u00F6\u00D6\u00DF\u00F8\u00E4\u00A2\u00A3$\u20AC\u9999\u80A0")
        
    def test_full_characters_set_read(self):
        """Read full characters set for each type"""
        
        self.assertEqual(self.records["dna"][1].seq.tostring(),"ACGTMRWSYKVHDBXN.-" )
        self.assertEqual(self.records["rna"][1].seq.tostring(),"ACGUMRWSYKVHDBXN.-" )
        self.assertEqual(self.records["protein"][1].seq.tostring(),"ABCDEFGHIJKLMNOPQRSTUVWXYZ.-*")
    
    def test_duplicated_property(self):
        """Read property with multiple values"""
        
        self.assertEqual(self.records["protein"][2].annotations["test"],[u"1",u"2",u"3"])
        
    def test_duplicated_alternativeID(self):
        """Read multiple alternative identifier form single source"""
        
        self.assertEqual(self.records["protein"][2].dbxrefs,[u"someDB:G001",u"someDB:G002"])
    
    def test_read_minial_required(self):
        
        minimalRecord = SeqRecord(id="test",seq=Seq("abc"))
        minimalRecord.annotations["source"] = u"Ensembl"
        
        self.assertEqual(self.records["rna"][3].name,minimalRecord.name)
        self.assertEqual(self.records["dna"][3].annotations,minimalRecord.annotations)
        self.assertEqual(self.records["rna"][3].dbxrefs,minimalRecord.dbxrefs)
        self.assertEqual(self.records["protein"][3].description,minimalRecord.description)
        
    def test_species(self):
        
        self.assertEqual(self.records["rna"][1].annotations["organism"],"Mus musculus")
        self.assertEqual(self.records["rna"][1].annotations["ncbi_taxid"],"10090")
        
        self.assertEqual(self.records["rna"][0].annotations["organism"],"Gallus gallus")
        self.assertEqual(self.records["rna"][0].annotations["ncbi_taxid"],"9031")
        
        
    def test_local_source_definition(self):
        
        self.assertEqual(self.records["protein"][4].annotations["source"],u"Uniprot")

    def test_empty_description(self):
         
        self.assertEqual(self.records["rna"][4].description,SeqRecord(id="",seq=Seq("")).description)


class TestReadAndWrite(unittest.TestCase):
    
    def setUp(self):
        open_test_files()

    def tearDown(self):
        close_test_files()
    
    def test_read_write_rna(self):
        
        read1_records = list(SeqIO.parse(test_files["rna"][2],"seqxml"))
        self._write_parse_and_compare(read1_records)
    
    def test_read_write_dna(self):
        
        read1_records = list(SeqIO.parse(test_files["dna"][2],"seqxml"))
        self._write_parse_and_compare(read1_records)
    
    def test_read_write_protein(self):
        
        read1_records = list(SeqIO.parse(test_files["protein"][2],"seqxml"))
        self._write_parse_and_compare(read1_records)
        
    
    def _write_parse_and_compare(self,read1_records):
        
        handle = StringIO()
        
        SeqIO.write(read1_records,handle,"seqxml")
        
        handle.seek(0)
        read2_records = list(SeqIO.parse(handle,"seqxml"))
        
        self.assertEquals(len(read1_records),len(read2_records))
        
        for record1,record2 in zip(read1_records,read2_records):
            assert_equal_records(self,record1,record2)
        
    
class TestReadCorruptFiles(unittest.TestCase):
    
    def setUp(self):
        open_test_files()

    def tearDown(self):
        close_test_files()
        
    
    def test_for_errors(self):
        
        for _file in corrupt_files:
            seqIO = SeqIO.parse(corrupt_files[_file][1],"seqxml")
            
            self.assertRaises(ValueError,seqIO.next);

    
class TestWriteCorruptRecords(unittest.TestCase):
    pass
        
if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)