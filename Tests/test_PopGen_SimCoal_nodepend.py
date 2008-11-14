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
from Bio.PopGen import SimCoal
from Bio.PopGen.SimCoal.Template import generate_simcoal_from_template

#Tests simcoal related code. Note: this case doesn't require simcoal
#test_PopGen_SimCoal tests code that requires simcoal

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
    tests = [TemplateTest]
    
    for test in tests:
        cur_suite = test_loader.loadTestsFromTestCase(test)
        test_suite.addTest(cur_suite)

    return test_suite

class TemplateTest(unittest.TestCase):
    def t_template_full(self):
        """Full template creation test
        """
        generate_simcoal_from_template('simple',
            [(1, [('SNP', [24, 0.0005, 0.0])])],
            [('sample_size', [30]),
             ('pop_size', [100])],
            'PopGen')
        #Confirm the files match (ignoring any switch of line endings
        #possible if the input file used a different OS convention)
        old = open('PopGen' + os.sep + 'simple.par', "rU").readlines()
        new = open('PopGen' + os.sep + 'simple_100_30.par').readlines()
        assert old==new, "Error - Old:\n%s\n\nNew:\n%s\n" % (old, new)
        #assert(os.stat('PopGen' + os.sep + 'simple.par').st_size ==
        #       os.stat('PopGen' + os.sep + 'simple_100_30.par').st_size)

    def tearDown(self):
        os.remove('PopGen' + os.sep + 'tmp.par')
        os.remove('PopGen' + os.sep + 'simple_100_30.par')


if __name__ == "__main__":
    sys.exit(run_tests(sys.argv))
