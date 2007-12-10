# Copyright 2007 by Tiago Antao <tiagoantao@gmail.com>.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import commands
import os
import shutil
import sys
import tempfile
import unittest
from Bio.PopGen import SimCoal
from Bio.PopGen.SimCoal.Controller import SimCoalController
from Bio import MissingExternalDependencyError

#Tests simcoal related code. Note: this case requires simcoal
#test_PopGen_SimCoal_nodepend tests code that does not require simcoal

found = False
for path in os.environ['PATH'].split(os.pathsep):
    try:
        list = os.listdir(path)
        for file in os.listdir(path):
            if file.startswith('simcoal2'):
                found = True
                simcoal_dir = path
    except os.error:
        pass #Path doesn't exist - correct to pass
if not found:
    raise MissingExternalDependencyError("SimCoal not found (not a problem if you do not intend to use it).")

def run_tests(argv):
    test_suite = testing_suite()
    runner = unittest.TextTestRunner(sys.stdout, verbosity = 2)
    runner.run(test_suite)

def testing_suite():
    """Generate the suite of tests.
    """
    print "Running simcoal tests, which might take some time, please wait"
    test_suite = unittest.TestSuite()

    test_loader = unittest.TestLoader()
    test_loader.testMethodPrefix = 't_'
    tests = [AppTest]
    
    for test in tests:
        cur_suite = test_loader.loadTestsFromTestCase(test)
        test_suite.addTest(cur_suite)

    return test_suite

class AppTest(unittest.TestCase):
    """Tests simcoal execution via biopython.
    """
    def tearDown(self):
        for file in os.listdir('PopGen' + os.sep + 'simple'):
            os.remove(os.sep.join(['PopGen', 'simple', file]))
        os.rmdir('PopGen' + os.sep + 'simple')

    def t_simcoal(self):
        """Test simcoal execution.
        """
        ctrl = SimCoalController(simcoal_dir)
        ctrl.run_simcoal('simple.par', 50, par_dir = 'PopGen')
        assert( len(os.listdir('PopGen' + os.sep + 'simple')) == 52)


if __name__ == "__main__":
    sys.exit(run_tests(sys.argv))
