# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
import sys
import unittest

from Bio import MyModule

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
    tests = [MyModuleTestOne]
    
    for test in tests:
        cur_suite = test_loader.loadTestsFromTestCase(test)
        test_suite.addTest(cur_suite)

    return test_suite

class MyModuleTestOne(unittest.TestCase):
    def setUp(self):
        self.handle = open("MyModule/input_file.txt")
        self.output_file = "MyModule/output.txt"

    def tearDown(self):
        self.handle.close()
        if os.path.exists(self.output_file):
            os.remove(self.output_file)

    def t_simple_parsing(self):
        """Test to be sure that MyModule can parse input files.
        """
        parser = MyModule.RecordParser()
        rec = parser.parse(self.handle)
        assert rec.id = "TheExpectedID"

    def t_output(self):
        """Ensure that we can write proper output files.
        """
        parser = MyModule.RecordParser()
        rec = parser.parse(self.handle)
        output_handle = open(self.output_file, "w")
        rec.write_to_file(output_handle)
        output_handle.close()

if __name__ == "__main__":
    sys.exit(run_tests(sys.argv))
