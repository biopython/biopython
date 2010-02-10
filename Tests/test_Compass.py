"""Tests for parsing Compass output.
"""
import os
import unittest

from Bio import Compass


class CompassTest(unittest.TestCase):
    def setUp(self):
        file_dir = os.path.join("Compass")
        self.test_files = [
          os.path.join(file_dir, "comtest1"),
          os.path.join(file_dir, "comtest2"),
          os.path.join(file_dir, "comtest3")]

    def testCompassScanAndConsume(self):
        com_record = Compass.read(open(self.test_files[0]))

        self.assertEquals("60456.blo.gz.aln", com_record.query)
        self.assertEquals("60456.blo.gz.aln", com_record.hit)
        self.assertEquals(0.5, com_record.gap_threshold)

        self.assertEquals(388, com_record.query_length)
        self.assertEquals(386, com_record.query_filtered_length)
        self.assertEquals(388, com_record.hit_length)
        self.assertEquals(386, com_record.hit_filtered_length)

        self.assertEquals(399, com_record.query_nseqs)
        self.assertEquals(12.972, com_record.query_neffseqs)
        self.assertEquals(399, com_record.hit_nseqs)
        self.assertEquals(12.972, com_record.hit_neffseqs)
                                                      
        self.assertEquals(2759, com_record.sw_score)
        self.assertEquals(float("0.00e+00"), com_record.evalue)

    def testCompassParser(self):
        com_record = Compass.read(open(self.test_files[0]))

        self.assertEquals("60456.blo.gz.aln", com_record.query)

    def testCompassIteratorEasy(self):
        records = Compass.parse(open(self.test_files[0]))

        com_record = records.next()
        self.assertEquals("60456.blo.gz.aln", com_record.query)

        self.assertRaises(StopIteration, records.next)
        
    def testCompassIteratorHard(self):
        records = Compass.parse(open(self.test_files[1]))

        com_record = records.next()
        self.assertEquals("allscop//14982.blo.gz.aln", com_record.hit)
        self.assertEquals(float('1.01e+03'), com_record.evalue)
        
        com_record = records.next()
        self.assertEquals("allscop//14983.blo.gz.aln", com_record.hit)
        self.assertEquals(float('1.01e+03'), com_record.evalue)
                                      
        com_record = records.next()
        self.assertEquals("allscop//14984.blo.gz.aln", com_record.hit)
        self.assertEquals(float('5.75e+02'), com_record.evalue)
                                                                                    
    def testAlignmentParsingOne(self):
        records = Compass.parse(open(self.test_files[1]))

        com_record = records.next()
        self.assertEquals(178, com_record.query_start)
        self.assertEquals("KKDLEEIAD", com_record.query_aln)
        self.assertEquals(9, com_record.hit_start)
        self.assertEquals("QAAVQAVTA", com_record.hit_aln)
        self.assertEquals("++ ++++++", com_record.positives)
        
        com_record = records.next()
        com_record = records.next()
        self.assertEquals(371, com_record.query_start)
        self.assertEquals("LEEAMDRMER~~~V", com_record.query_aln)
        self.assertEquals(76, com_record.hit_start)
        self.assertEquals("LQNFIDQLDNpddL", com_record.hit_aln)
        self.assertEquals("+ ++++ + +   +", com_record.positives)

    def testAlignmentParsingTwo(self):
        records = Compass.parse(open(self.test_files[0]))
        
        com_record = records.next()
        self.assertEquals(2, com_record.query_start)
        self.assertEquals(2, com_record.hit_start)
        self.assertEquals("LKERKL", com_record.hit_aln[-6:])

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
