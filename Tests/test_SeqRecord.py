# Copyright 2009 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""SeqFeature related tests for SeqRecord objects from Bio.SeqIO.

Initially this takes matched tests of GenBank and FASTA files from the NCBI
and confirms they are consistent using our different parsers.
"""
import unittest
from Bio import SeqIO
from Bio.Alphabet import generic_dna, generic_rna, generic_protein
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
from Bio.SeqFeature import WithinPosition, BeforePosition, AfterPosition, OneOfPosition

class SeqRecordCreation(unittest.TestCase):
    """Test basic creation of SeqRecords."""
    def test_annotations(self):
        """Pass in annotations to SeqRecords"""
        rec = SeqRecord(Seq("ACGT", generic_dna),
                        id="Test", name="Test", description="Test")
        self.assertEqual(rec.annotations, {})
        rec = SeqRecord(Seq("ACGT", generic_dna),
                        id="Test", name="Test", description="Test",
                        annotations={"test" : ["a test"]})
        self.assertEqual(rec.annotations["test"], ["a test"])

    def test_letter_annotations(self):
        """Pass in letter annotations to SeqRecords"""
        rec = SeqRecord(Seq("ACGT", generic_dna),
                        id="Test", name="Test", description="Test")
        self.assertEqual(rec.annotations, {})
        rec = SeqRecord(Seq("ACGT", generic_dna),
                        id="Test", name="Test", description="Test",
                        letter_annotations={"test" : [1, 2, 3, 4]})
        self.assertEqual(rec.letter_annotations["test"], [1, 2, 3, 4])
        #Now try modifying it to a bad value...
        try:
            rec.letter_annotations["bad"] = "abc"
            self.assertTrue(False, "Adding a bad letter_annotation should fail!")
        except (TypeError, ValueError), e:
            pass
        #Now try setting it afterwards to a bad value...
        rec = SeqRecord(Seq("ACGT", generic_dna),
                        id="Test", name="Test", description="Test")
        try:
            rec.letter_annotations={"test" : [1, 2, 3]}
            self.assertTrue(False, "Changing to bad letter_annotations should fail!")
        except (TypeError, ValueError), e:
            pass
        #Now try setting it at creation time to a bad value...
        try:
            rec = SeqRecord(Seq("ACGT", generic_dna),
                            id="Test", name="Test", description="Test",
                            letter_annotations={"test" : [1, 2, 3]})
            self.assertTrue(False, "Wrong length letter_annotations should fail!")
        except (TypeError, ValueError), e:
            pass

class SeqRecordMethods(unittest.TestCase):
    """Test SeqRecord methods."""

    def setUp(self) :
        f0 = SeqFeature(FeatureLocation(0,26), type="source",
                        qualifiers={"mol_type":["fake protein"]})
        f1 = SeqFeature(FeatureLocation(0,ExactPosition(10)))
        f2 = SeqFeature(FeatureLocation(WithinPosition(12, left=12,right=15),BeforePosition(22)))
        f3 = SeqFeature(FeatureLocation(AfterPosition(16),
                                        OneOfPosition(26, [ExactPosition(25),AfterPosition(26)])))
        self.record = SeqRecord(Seq("ABCDEFGHIJKLMNOPQRSTUVWZYX", generic_protein),
                                id="TestID", name="TestName", description="TestDescr",
                                dbxrefs=["TestXRef"], annotations={"k":"v"},
                                letter_annotations = {"fake":"X"*26},
                                features = [f0,f1,f2,f3])

    def test_slice_variantes(self):
        """Simple slices using different start/end values"""
        for start in range(-30,30)+[None] :
            for end in range(-30,30)+[None] :
                if start is None and end is None : continue
                rec = self.record[start:end]
                seq = self.record.seq[start:end]
                seq_str = str(self.record.seq)[start:end]
                self.assertEqual(seq_str, str(seq))
                self.assertEqual(seq_str, str(rec.seq))
                self.assertEqual("X"*len(seq_str), rec.letter_annotations["fake"])

    def test_slice_simple(self):
        """Simple slice"""
        rec = self.record
        self.assertEqual(len(rec), 26)
        left = rec[:10]
        self.assertEqual(str(left.seq), str(rec.seq[:10]))
        right = rec[-10:]
        self.assertEqual(str(right.seq), str(rec.seq[-10:]))
        mid = rec[12:22]
        self.assertEqual(str(mid.seq), str(rec.seq[12:22]))
        for sub in [left, right, mid] :
            self.assertEqual(len(sub), 10)
            self.assertEqual(sub.id, "TestID")
            self.assertEqual(sub.name, "TestName")
            self.assertEqual(sub.description, "TestDescr")
            self.assertEqual(sub.letter_annotations, {"fake":"X"*10})
            self.assertEqual(sub.dbxrefs, []) # May change this...
            self.assertEqual(sub.annotations, {}) # May change this...
            self.assertEqual(len(sub.features), 1)
            #By construction, each feature matches the full sliced region:
            self.assertEqual(str(sub.features[0].extract(sub.seq)), str(sub.seq))
            self.assertEqual(sub.features[0].extract(str(sub.seq)), str(sub.seq))

    def test_slice_zero(self):
        """Zero slice"""
        rec = self.record
        self.assertEqual(len(rec), 26)
        self.assertEqual(len(rec[2:-2]), 22)
        self.assertEqual(len(rec[5:2]), 0)
        self.assertEqual(len(rec[5:2][2:-2]), 0)

    def test_add_simple(self):
        """Simple addition"""
        rec = self.record + self.record
        self.assertEqual(len(rec), 52)
        self.assertEqual(rec.id, "TestID")
        self.assertEqual(rec.name, "TestName")
        self.assertEqual(rec.description, "TestDescr")
        self.assertEqual(rec.dbxrefs, ["TestXRef"])
        self.assertEqual(rec.annotations, {"k":"v"})
        self.assertEqual(rec.letter_annotations, {"fake":"X"*52})
        self.assertEqual(len(rec.features), 2*len(self.record.features))

    def test_add_seq(self):
        """Simple addition of Seq or string"""
        for other in [Seq("BIO"), "BIO"] :
            rec = self.record + other # will use SeqRecord's __add__ method
            self.assertEqual(len(rec), 26+3)
            self.assertEqual(str(rec.seq), str(self.record.seq)+"BIO")
            self.assertEqual(rec.id, "TestID")
            self.assertEqual(rec.name, "TestName")
            self.assertEqual(rec.description, "TestDescr")
            self.assertEqual(rec.dbxrefs, ["TestXRef"])
            self.assertEqual(rec.annotations, {"k":"v"})
            self.assertEqual(rec.letter_annotations, {})
            self.assertEqual(len(rec.features), len(self.record.features))
            self.assertEqual(rec.features[0].type, "source")
            self.assertEqual(rec.features[0].location.nofuzzy_start, 0)
            self.assertEqual(rec.features[0].location.nofuzzy_end, 26) #not +3

    def test_add_seqrecord(self):
        """Simple left addition of SeqRecord from genbank file."""
        other = SeqIO.read("GenBank/dbsource_wrap.gb", "gb")
        other.dbxrefs = ["dummy"]
        rec = self.record + other
        self.assertEqual(len(rec), len(self.record)+len(other))
        self.assertEqual(str(rec.seq), str(self.record.seq)+str(other.seq))
        self.assertEqual(rec.id, "<unknown id>")
        self.assertEqual(rec.name, "<unknown name>")
        self.assertEqual(rec.description, "<unknown description>")
        self.assertEqual(rec.dbxrefs, ["TestXRef", "dummy"])
        self.assertEqual(len(rec.annotations), 0)
        self.assertEqual(len(rec.letter_annotations),0)
        self.assertEqual(len(rec.features),
                         len(self.record.features) + len(other.features))
        self.assertEqual(rec.features[0].type, "source")
        self.assertEqual(rec.features[0].location.nofuzzy_start, 0)
        self.assertEqual(rec.features[0].location.nofuzzy_end, len(self.record)) #not +3
        i = len(self.record.features)
        self.assertEqual(rec.features[i].type, "source")
        self.assertEqual(rec.features[i].location.nofuzzy_start, len(self.record))
        self.assertEqual(rec.features[i].location.nofuzzy_end, len(rec))

    def test_add_seq_left(self):
        """Simple left addition of Seq or string"""
        for other in [Seq("BIO"), "BIO"] :
            rec = other + self.record # will use SeqRecord's __radd__ method
            self.assertEqual(len(rec), 26+3)
            self.assertEqual(str(rec.seq), "BIO"+str(self.record.seq))
            self.assertEqual(rec.id, "TestID")
            self.assertEqual(rec.name, "TestName")
            self.assertEqual(rec.description, "TestDescr")
            self.assertEqual(rec.dbxrefs, ["TestXRef"])
            self.assertEqual(rec.annotations, {"k":"v"})
            self.assertEqual(rec.letter_annotations, {})
            self.assertEqual(len(rec.features), len(self.record.features))
            self.assertEqual(rec.features[0].type, "source")
            self.assertEqual(rec.features[0].location.nofuzzy_start, 3)
            self.assertEqual(rec.features[0].location.nofuzzy_end, 26+3)
            
    def test_slice_add_simple(self):
        """Simple slice and add"""
        for cut in range(27) :
            rec = self.record[:cut] + self.record[cut:]
            self.assertEqual(str(rec.seq), str(self.record.seq))
            self.assertEqual(len(rec), 26)
            self.assertEqual(rec.id, "TestID")
            self.assertEqual(rec.name, "TestName")
            self.assertEqual(rec.description, "TestDescr")
            self.assertEqual(rec.dbxrefs, []) # May change this...
            self.assertEqual(rec.annotations, {}) # May change this...
            self.assertEqual(rec.letter_annotations, {"fake":"X"*26})
            self.assertTrue(len(rec.features) <= len(self.record.features))

    def test_slice_add_shift(self):
        """Simple slice and add to shift"""
        for cut in range(27) :
            rec = self.record[cut:] + self.record[:cut]
            self.assertEqual(str(rec.seq), str(self.record.seq[cut:] + self.record.seq[:cut]))
            self.assertEqual(len(rec), 26)
            self.assertEqual(rec.id, "TestID")
            self.assertEqual(rec.name, "TestName")
            self.assertEqual(rec.description, "TestDescr")
            self.assertEqual(rec.dbxrefs, []) # May change this...
            self.assertEqual(rec.annotations, {}) # May change this...
            self.assertEqual(rec.letter_annotations, {"fake":"X"*26})
            self.assertTrue(len(rec.features) <= len(self.record.features))
            
if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
