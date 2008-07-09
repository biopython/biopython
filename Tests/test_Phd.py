import os
import sys
import unittest

from Bio.Sequencing import Phd

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
    tests = [PhdTestOne]
    
    for test in tests:
        cur_suite = test_loader.loadTestsFromTestCase(test)
        test_suite.addTest(cur_suite)

    return test_suite

class PhdTestOne(unittest.TestCase):
    def setUp(self):
        self.handle = open("Phd/phd1")

    def tearDown(self):
        self.handle.close()

    def t_check_record_parser(self):
        """Test to check that record parser parses all records of a contig.
        """
        rs = Phd.parse(self.handle)
        x=0
        print
        for r in rs:
            x+=1
            print '---- Record',x,'----'
            print 'Filename:',r.file_name
            comments=r.comments.items()
            comments.sort()
            print '\n'.join([c[0]+': '+str(c[1]) for c in comments])
            allsites=['-'.join(s) for s in r.sites]
            print allsites[:10],'...',allsites[len(allsites)/2-5:len(allsites)/2+5],'...',allsites[-10:]
            print r.seq.tostring()[:10],r.seq.tostring()[-10:]
            print r.seq_trimmed.tostring()[:10],r.seq_trimmed.tostring()[-10:]
        
if __name__ == "__main__":
    sys.exit(run_tests(sys.argv))
