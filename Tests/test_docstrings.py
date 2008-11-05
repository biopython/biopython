#!/usr/bin/env python
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import doctest, unittest

from Bio import Seq, SeqRecord, SeqIO, AlignIO
test_modules = [Seq, SeqRecord, SeqIO, AlignIO]

test_suite = unittest.TestSuite((doctest.DocTestSuite(module) \
                                for module in test_modules))

#Using sys.stdout prevent this working nicely when run from idle:
#runner = unittest.TextTestRunner(sys.stdout, verbosity = 0)

#Using verbosity = 0 means we won't have to regenerate the unit
#test output file used by the run_tests.py framework whenever a
#new module or doctest is added.
runner = unittest.TextTestRunner(verbosity = 0)
runner.run(test_suite)
