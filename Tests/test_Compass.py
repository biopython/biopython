"""Tests for parsing Compass output.
"""
import os
import sys
import unittest

from Bio import Compass

def run_tests(argv):
    test_suite = testing_suite()
    runner = unittest.TextTestRunner(sys.stdout, verbosity = 2)
    runner.run(test_suite)

def testing_suite():
    """Generate the suite of tests.
    """
    test_suite = unittest.TestSuite()

    test_loader = unittest.TestLoader()
    test_loader.testMethodPrefix = 'test'
    tests = [CompassTest]
    
    for test in tests:
        cur_suite = test_loader.loadTestsFromTestCase(test)
        test_suite.addTest(cur_suite)

    return test_suite

class CompassTest(unittest.TestCase):
    def setUp(self):
        file_dir = os.path.join("Compass")
        self.test_files = [
          os.path.join(file_dir, "comtest1"),
          os.path.join(file_dir, "comtest2")]

    def testCompassScanAndConsume(self):
        cons = Compass._Consumer()
        scan = Compass._Scanner()
        scan.feed(open(self.test_files[0]), cons)

        com_record = cons.data

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
        parser = Compass.RecordParser()
        com_record = parser.parse(open(self.test_files[0]))

        self.assertEquals("60456.blo.gz.aln", com_record.query)

    def testCompassIteratorEasy(self):
        it = Compass.Iterator(open(self.test_files[0]))

        com_record = it.next()
        self.assertEquals("60456.blo.gz.aln", com_record.query)

        com_record = it.next()
        self.assertEquals(None, com_record)
        pass
        
    def testCompassIteratorHard(self):
        it = Compass.Iterator(open(self.test_files[1]))

        com_record = it.next()
        self.assertEquals("allscop//14982.blo.gz.aln", com_record.hit)
        self.assertEquals(float('1.01e+03'), com_record.evalue)
        
        com_record = it.next()
        self.assertEquals("allscop//14983.blo.gz.aln", com_record.hit)
        self.assertEquals(float('1.01e+03'), com_record.evalue)
                                      
        com_record = it.next()
        self.assertEquals("allscop//14984.blo.gz.aln", com_record.hit)
        self.assertEquals(float('5.75e+02'), com_record.evalue)
                                                                                    
    def testAlignmentParsingOne(self):
        it = Compass.Iterator(open(self.test_files[1]))

        com_record = it.next()
        self.assertEquals(178, com_record.query_start)
        self.assertEquals("KKDLEEIAD", com_record.query_aln)
        self.assertEquals(9, com_record.hit_start)
        self.assertEquals("QAAVQAVTA", com_record.hit_aln)
        self.assertEquals("++ ++++++", com_record.positives)
        
        com_record = it.next()
        com_record = it.next()
        self.assertEquals(371, com_record.query_start)
        self.assertEquals("LEEAMDRMER~~~V", com_record.query_aln)
        self.assertEquals(76, com_record.hit_start)
        self.assertEquals("LQNFIDQLDNpddL", com_record.hit_aln)
        self.assertEquals("+ ++++ + +   +", com_record.positives)

    def testAlignmentParsingTwo(self):
        it = Compass.Iterator(open(self.test_files[0]))
        
        com_record = it.next()
        self.assertEquals(2, com_record.query_start)
        self.assertEquals(2, com_record.hit_start)
        self.assertEquals("LKERKL", com_record.hit_aln[-6:])

if __name__ == "__main__":
    sys.exit(run_tests(sys.argv))
