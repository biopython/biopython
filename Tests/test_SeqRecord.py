# Copyright 2009 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""SeqFeature related tests for SeqRecord objects from Bio.SeqIO.

Initially this takes matched tests of GenBank and FASTA files from the NCBI
and confirms they are consistent using our different parsers.
"""
import unittest
from Bio.Alphabet import generic_dna, generic_rna, generic_protein
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

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
            self.assert_(False, "Adding a bad letter_annotation should fail!")
        except (TypeError, ValueError), e:
            pass
        #Now try setting it afterwards to a bad value...
        rec = SeqRecord(Seq("ACGT", generic_dna),
                        id="Test", name="Test", description="Test")
        try:
            rec.letter_annotations={"test" : [1, 2, 3]}
            self.assert_(False, "Changing to bad letter_annotations should fail!")
        except (TypeError, ValueError), e:
            pass
        #Now try setting it at creation time to a bad value...
        try:
            rec = SeqRecord(Seq("ACGT", generic_dna),
                            id="Test", name="Test", description="Test",
                            letter_annotations={"test" : [1, 2, 3]})
            self.assert_(False, "Wrong length letter_annotations should fail!")
        except (TypeError, ValueError), e:
            pass

class SeqRecordSlicing(unittest.TestCase):
    """Test SeqRecord slicing."""

    def test_simple(self):
        """Simple slice"""
        rec = SeqRecord(Seq("ABCDEFGHIJKLMNOPQRSTUVWZYX", generic_protein),
                        id="TestID", name="TestName", description="TestDescr",
                        dbxrefs=["TestXRef"], annotations={"k":"v"},
                        letter_annotations = {"fake":"X"*26})
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

class SeqRecordAddition(unittest.TestCase):
    """Test SeqRecord addition."""

    def test_simple(self):
        """Simple addition"""
        rec1 = SeqRecord(Seq("ABCDEFGHIJKLMNOPQRSTUVWZYX", generic_protein),
                         id="TestID", name="TestName", description="TestDescr",
                         dbxrefs=["TestXRef"], annotations={"k":"v"},
                         letter_annotations = {"fake":[0,1]*13})
        rec = rec1 + rec1
        self.assertEqual(len(rec), 52)
        self.assertEqual(rec.id, "TestID")
        self.assertEqual(rec.name, "TestName")
        self.assertEqual(rec.description, "TestDescr")
        self.assertEqual(rec.dbxrefs, ["TestXRef"])
        self.assertEqual(rec.annotations, {"k":"v"})
        self.assertEqual(rec.letter_annotations, {"fake":[0,1]*26})
        

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
