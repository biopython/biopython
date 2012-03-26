# Copyright 2010 by Thomas Schmitt.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
import unittest
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from StringIO import StringIO

test_files = {
    "dna" : ["SeqXML/dna_example.xml",4],
    "rna" : ["SeqXML/rna_example.xml",5],
    "protein" : ["SeqXML/protein_example.xml",5],
    "globalSpecies" : ["SeqXML/global_species_example.xml",2],
}

corrupt_files = ["SeqXML/corrupt_example1.xml",
                 "SeqXML/corrupt_example2.xml",
                ]


def assert_equal_records(testCase,record_a,record_b):
    testCase.assertEqual(record_a.id,record_b.id)
    testCase.assertEqual(record_a.name,record_b.name)
    testCase.assertEqual(record_a.description,record_b.description)
    testCase.assertEqual(record_a.seq.tostring(),record_b.seq.tostring())
    testCase.assertEqual(record_a.dbxrefs,record_b.dbxrefs)
    testCase.assertEqual(record_a.annotations,record_a.annotations)
    

class TestSimpleRead(unittest.TestCase):
        
    def test_check_SeqIO(self):
        """Files readable using parser via SeqIO."""
        for key in test_files:
            records = list(SeqIO.parse(test_files[key][0],"seqxml"))
            self.assertEqual(len(records),test_files[key][1])
            
class TestDetailedRead(unittest.TestCase):
    
    records = {}
    
    def setUp(self):
        for key in test_files:
            self.records[key] = list(SeqIO.parse(test_files[key][0],"seqxml"))

    def test_special_characters_desc(self):
        """Read special XML characters in description."""
        self.assertEqual(self.records["dna"][2].description, u'some special characters in the description\n<tag> "quoted string"')

    #TODO - Fix this failure under Windows with Python 3.1 and 3.2   
    if not (sys.platform=="win32" and sys.version_info[0] >= 3):
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
        
    def test_duplicated_dbxref(self):
        """Read multiple cross references to a single source"""
        self.assertEqual(self.records["protein"][2].dbxrefs,[u"someDB:G001",u"someDB:G002"])
    
    def test_read_minimal_required(self):
        """Check minimal record."""
        minimalRecord = SeqRecord(id="test",seq=Seq("abc"))
        minimalRecord.annotations["source"] = u"Ensembl"
        
        self.assertEqual(self.records["rna"][3].name,minimalRecord.name)
        self.assertEqual(self.records["dna"][3].annotations,minimalRecord.annotations)
        self.assertEqual(self.records["rna"][3].dbxrefs,minimalRecord.dbxrefs)
        self.assertEqual(self.records["protein"][3].description,minimalRecord.description)
        
    def test_local_species(self):
        """Check local species."""
        self.assertEqual(self.records["rna"][1].annotations["organism"],"Mus musculus")
        self.assertEqual(self.records["rna"][1].annotations["ncbi_taxid"],"10090")
        
        self.assertEqual(self.records["rna"][0].annotations["organism"],"Gallus gallus")
        self.assertEqual(self.records["rna"][0].annotations["ncbi_taxid"],"9031")
        
    def test_global_species(self):
        """Check global species."""
        self.assertEqual(self.records["globalSpecies"][0].annotations["organism"],"Mus musculus")
        self.assertEqual(self.records["globalSpecies"][0].annotations["ncbi_taxid"],"10090")
        
        self.assertEqual(self.records["globalSpecies"][1].annotations["organism"],"Homo sapiens")
        self.assertEqual(self.records["globalSpecies"][1].annotations["ncbi_taxid"],"9606")

    def test_local_source_definition(self):
        """Check local source."""
        self.assertEqual(self.records["protein"][4].annotations["source"],u"Uniprot")

    def test_empty_description(self):
        """Check empty description."""
        self.assertEqual(self.records["rna"][4].description,SeqRecord(id="",seq=Seq("")).description)


class TestReadAndWrite(unittest.TestCase):
    
    def test_read_write_rna(self):
        """Read and write RNA."""
        read1_records = list(SeqIO.parse(test_files["rna"][0],"seqxml"))
        self._write_parse_and_compare(read1_records)
    
    def test_read_write_dna(self):
        """Read and write DNA."""
        read1_records = list(SeqIO.parse(test_files["dna"][0],"seqxml"))
        self._write_parse_and_compare(read1_records)
    
    def test_read_write_protein(self):
        """Read and write protein."""
        read1_records = list(SeqIO.parse(test_files["protein"][0],"seqxml"))
        self._write_parse_and_compare(read1_records)
        
    def test_read_write_globalSpecies(self):
        """Read and write global species."""
        read1_records = list(SeqIO.parse(test_files["globalSpecies"][0],"seqxml"))
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
    
    def test_for_errors(self):
        """Handling of corrupt files."""        
        for filename in corrupt_files:
            iterator = SeqIO.parse(filename,"seqxml")
            self.assertRaises(ValueError,iterator.next);

        
if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
