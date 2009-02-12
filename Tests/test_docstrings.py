#!/usr/bin/env python
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import sys
from Bio import MissingExternalDependencyError
if sys.version_info[:2] < (2, 4):
    #On python 2.3, doctest uses slightly different formatting
    #which would be a problem as the expected output won't match.
    #Also, it can't cope with <BLANKLINE> in a doctest string.
    raise MissingExternalDependencyError(\
          "This unit test requires Python 2.4 or later")
import doctest, unittest

# This is the list of modules containing docstring tests.
# If you develop docstring tests for other modules, please add
# those modules here.
DOCTEST_MODULES = ["Bio.Seq",
                   "Bio.SeqRecord",
                   "Bio.SeqIO",
                   "Bio.Align.Generic",
                   "Bio.AlignIO",
                   "Bio.KEGG.Compound",
                   "Bio.KEGG.Enzyme",
                   "Bio.Wise",
                   "Bio.Wise.psw",
                  ]
#Silently ignore any doctests for modules requiring numpy!
try:
    import numpy
    DOCTEST_MODULES.extend(["Bio.Statistics.lowess"])
except ImportError:
    pass

#Assuming we keep adding doctests, we'll need to have this module
#list populated automatically - or otherwise integrated into the
#Biopython test framework.


#Can't use fromlist=name.split(".") until python 2.5+
modules = [doctest.DocTestSuite(module = __import__(name,
                                                    None,
                                                    None,
                                                    name.split("."))) \
           for name in DOCTEST_MODULES]

test_suite = unittest.TestSuite(modules)

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
