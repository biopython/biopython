#!/usr/bin/env python
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import doctest, unittest, sys

from Bio import Seq, SeqRecord, SeqIO, AlignIO
test_modules = [Seq, SeqRecord, SeqIO, AlignIO]

test_suite = unittest.TestSuite((doctest.DocTestSuite(module) \
                                for module in test_modules))

#Use stdout so that run_tests.py can capture the output.
#Even verbosity=0 outputs something, e.g.
"""
----------------------------------------------------------------------
Ran 15 tests in 0.456s

OK
"""
#However, the only bits that change here are the number of tests and the
#time, and run_tests.py knows to ignore these lines.  This means we don't
#have to update output/test_docstrings when more doctests are added :)
runner = unittest.TextTestRunner(sys.stdout, verbosity = 0)
runner.run(test_suite)
