# Copyright 2006 by Tiago Antao <tiagoantao@gmail.com>.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import commands
import os
import shutil
import sys
import tempfile
import unittest
from Bio.PopGen import FDist
from Bio.PopGen.FDist import Controller
from Bio import MissingExternalDependencyError

#Tests fdist related code. Note: this case requires fdist
#test_PopGen_FDist_nodepend tests code that does not require fdist

not_found_types = ["command not found", ": not found"]
fdist_output = commands.getoutput("fdist2")
for not_found in not_found_types:
    if not_found in fdist_output:
        raise MissingExternalDependencyError("Fdist not found (not a problem if you do not intend to use it).")

def run_tests(argv):
    test_suite = testing_suite()
    runner = unittest.TextTestRunner(sys.stdout, verbosity = 2)
    runner.run(test_suite)

def testing_suite():
    """Generate the suite of tests.
    """
    print "Running fdist tests, which might take some time, please wait"
    test_suite = unittest.TestSuite()

    test_loader = unittest.TestLoader()
    test_loader.testMethodPrefix = 't_'
    tests = [AppTest]
    
    for test in tests:
        cur_suite = test_loader.loadTestsFromTestCase(test)
        test_suite.addTest(cur_suite)

    return test_suite

class AppTest(unittest.TestCase):
    """Tests the fdist suite of applications.
    """
    def _copyfile(self, inname, outname):
        shutil.copyfile(
            'PopGen' + os.sep + inname,
            self.dirname + os.sep + outname)

    def setUp(self):
        self.ctrl = Controller.FDistController()
        self.dirname = tempfile.mkdtemp()
        self._copyfile('data_fst_outfile', 'data_fst_outfile')
        self._copyfile('fdist1', 'infile')
        self._copyfile('out.dat', 'out.dat')
        self._copyfile('out.cpl', 'out.cpl')

    def tearDown(self):
        for file in os.listdir(self.dirname):
            os.remove(self.dirname + os.sep + file)
        os.rmdir(self.dirname)

    def t_datacal(self):
        """Test datacal execution.
        """
        fst, samp_size = self.ctrl.run_datacal(data_dir = self.dirname)
        assert (fst - 0.44 < 0.01)
        assert (samp_size == 11)

    def t_fdist(self):
        """Test fdist execution.
        """
        fst = self.ctrl.run_fdist(npops = 15, nsamples = 10, fst = 0.1,
                sample_size = 20, mut = 0, num_sims = 10000,
                data_dir = self.dirname)
        assert(abs(fst - 0.1) < 0.02) #Stochastic result...

    def t_fdist_force_fst(self):
        """Test fdist execution approximating Fst.
        """
        fst = self.ctrl.run_fdist_force_fst(npops = 15, nsamples = 10,
                fst = 0.1,
                sample_size = 20, mut = 0, num_sims = 10000,
                data_dir = self.dirname)
        assert(abs(fst - 0.09) < 0.05) #Stochastic result...

    def t_cplot(self):
        """Test cplot execution.
        """
        cpl_interval =self.ctrl.run_cplot(data_dir = self.dirname)
        assert(len(cpl_interval) == 8)

    def t_pv(self):
        """Test pv execution.
        """
        pv_data = self.ctrl.run_pv(data_dir = self.dirname)
        assert(len(pv_data) == 4)


if __name__ == "__main__":
    sys.exit(run_tests(sys.argv))
