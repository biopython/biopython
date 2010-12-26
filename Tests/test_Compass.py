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
        handle = open(self.test_files[0])
        com_record = Compass.read(handle)
        handle.close()

        self.assertEqual("60456.blo.gz.aln", com_record.query)
        self.assertEqual("60456.blo.gz.aln", com_record.hit)
        self.assertEqual(0.5, com_record.gap_threshold)

        self.assertEqual(388, com_record.query_length)
        self.assertEqual(386, com_record.query_filtered_length)
        self.assertEqual(388, com_record.hit_length)
        self.assertEqual(386, com_record.hit_filtered_length)

        self.assertEqual(399, com_record.query_nseqs)
        self.assertEqual(12.972, com_record.query_neffseqs)
        self.assertEqual(399, com_record.hit_nseqs)
        self.assertEqual(12.972, com_record.hit_neffseqs)
                                                      
        self.assertEqual(2759, com_record.sw_score)
        self.assertEqual(float("0.00e+00"), com_record.evalue)

    def testCompassParser(self):
        handle = open(self.test_files[0])
        com_record = Compass.read(handle)
        handle.close()

        self.assertEqual("60456.blo.gz.aln", com_record.query)

    def testCompassIteratorEasy(self):
        handle = open(self.test_files[0])
        records = Compass.parse(handle)
        com_record = records.next()
        self.assertEqual("60456.blo.gz.aln", com_record.query)
        self.assertRaises(StopIteration, records.next)
        handle.close()
        
    def testCompassIteratorHard(self):
        handle = open(self.test_files[1])
        records = Compass.parse(handle)

        com_record = records.next()
        self.assertEqual("allscop//14982.blo.gz.aln", com_record.hit)
        self.assertEqual(float('1.01e+03'), com_record.evalue)
        
        com_record = records.next()
        self.assertEqual("allscop//14983.blo.gz.aln", com_record.hit)
        self.assertEqual(float('1.01e+03'), com_record.evalue)
                                      
        com_record = records.next()
        self.assertEqual("allscop//14984.blo.gz.aln", com_record.hit)
        self.assertEqual(float('5.75e+02'), com_record.evalue)

        handle.close()
                                                                                    
    def testAlignmentParsingOne(self):
        handle = open(self.test_files[1])
        records = Compass.parse(handle)

        com_record = records.next()
        self.assertEqual(178, com_record.query_start)
        self.assertEqual("KKDLEEIAD", com_record.query_aln)
        self.assertEqual(9, com_record.hit_start)
        self.assertEqual("QAAVQAVTA", com_record.hit_aln)
        self.assertEqual("++ ++++++", com_record.positives)
        
        com_record = records.next()
        com_record = records.next()
        self.assertEqual(371, com_record.query_start)
        self.assertEqual("LEEAMDRMER~~~V", com_record.query_aln)
        self.assertEqual(76, com_record.hit_start)
        self.assertEqual("LQNFIDQLDNpddL", com_record.hit_aln)
        self.assertEqual("+ ++++ + +   +", com_record.positives)

        handle.close()

    def testAlignmentParsingTwo(self):
        handle = open(self.test_files[0])
        records = Compass.parse(handle)
        com_record = records.next()
        self.assertEqual(2, com_record.query_start)
        self.assertEqual(2, com_record.hit_start)
        self.assertEqual("LKERKL", com_record.hit_aln[-6:])
        handle.close()

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
