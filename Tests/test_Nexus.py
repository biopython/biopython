import os
import sys
import StringIO
import unittest

from Bio.Nexus import Nexus

def run_tests(argv):
    test_suite = testing_suite()
    runner = unittest.TextTestRunner(sys.stdout, verbosity = 2)
    runner.run(test_suite)

def testing_suite():
    """Generate the suite of tests."""
    test_suite = unittest.TestSuite()

    test_loader = unittest.TestLoader()
    test_loader.testMethodPrefix = 't_'
    tests = [NexusTest1]
    
    for test in tests:
        cur_suite = test_loader.loadTestsFromTestCase(test)
        test_suite.addTest(cur_suite)

    return test_suite

class NexusTest1(unittest.TestCase):
    def setUp(self):
        self.handle = open("Nexus/test_Nexus_input.nex")
        #Remove any output files from a previous test run
        f1=os.path.join('Nexus','f1.nex')
        f2=os.path.join('Nexus','f2.nex')
        if os.path.isfile(f1) : os.remove(f1)
        if os.path.isfile(f2) : os.remove(f2)

    def tearDown(self):
        self.handle.close()
        #Remove our test output files
        f1=os.path.join('Nexus','f1.nex')
        f2=os.path.join('Nexus','f2.nex')
        if os.path.isfile(f1) : os.remove(f1)
        if os.path.isfile(f2) : os.remove(f2)

    def print_dictnlist(self,name,d):
        """Output a list or a dict alphabetically by keys()"""
        if type(d)==dict:
            dk=d.keys()
            dk.sort()
            print '\n%s:\n%s' % (name, '\n'.join([str(k)+'\t'+str(d[k]) for k in dk]))
        else:
            print '\n%s:\n%s' % (name, '\n'.join([k for k in d]))

    def output_basics(self,n):
        """Output basic info of a nexus file."""
        # When run on Windows, want to match the expected test output
        # (in file output/test_Nexus) which uses "Nexus/f1.nex", not "Nexus\f1.nex"
        print '\nName:',n.filename.replace(os.sep,"/")
        print '\n',n.ntax, n.nchar, n.datatype, n.interleave, n.missing, n.gap
        self.print_dictnlist('Taxa',n.taxlabels)
        self.print_dictnlist('Charlabels',n.charlabels)
        self.print_dictnlist('Charsets',n.charsets)
        self.print_dictnlist('Taxsets',n.taxsets)
        ps=n.charpartitions.keys()
        ps.sort()
        print '\nCharpartitions:'
        for pk in ps:
                self.print_dictnlist(pk,n.charpartitions[pk])
        ps=n.taxpartitions.keys()
        ps.sort()
        print '\nTaxpartitions:'
        for pk in ps:
                self.print_dictnlist(pk,n.taxpartitions[pk])
        
    def t_NexusTest1(self):
        """Test Nexus module"""
        n=Nexus.Nexus(self.handle)
        # check data of main nexus file
        self.output_basics(n)
        # now we check excluding characters, deleting taxa, and exporting adjusted sets
        f1=os.path.join('Nexus','f1.nex')
        f2=os.path.join('Nexus','f2.nex')
        n.write_nexus_data(f1,delete=['t1','t7'],exclude=n.invert(n.charsets['big']))
        n.write_nexus_data(f2,delete=['t2_the_name'],exclude=range(3,40,4))
        nf1=Nexus.Nexus(f1)
        self.output_basics(nf1)
        nf2=Nexus.Nexus(f2)
        self.output_basics(nf2)
        # check the stepmatrix
        print n.weighted_stepmatrix(name='matrix_test')

    def t_TreeTest1(self):
        """Test Tree module."""
        n=Nexus.Nexus(self.handle)
        t3=n.trees[2]
        t2=n.trees[2]
        t3.root_with_outgroup(['t1','t5'])
        print t3
        print 'Monophyletic t6,t7,t8,t9:',t3.is_monophyletic(['t8','t9','t6','t7'])
        print 'Monophyletic t1,t5:',t3.is_monophyletic(['t1','t5'])
        t3.split(parent_id=t3.search_taxon('t9'))
        print 'The tree looks as follows:'
        t3.display()
        print 'self incompatibility...:',t3.is_compatible(t2,threshold=0.3)
        

if __name__ == "__main__":
    sys.exit(run_tests(sys.argv))

