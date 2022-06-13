# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for pairwise2 module using the default C functions.

This test file imports the TestCases from ``pairwise2_testCases.py``.
If you want to add more tests, do this over there.

"""

import unittest

# Import all test classes deliberately:
from pairwise2_testCases import *  # noqa: F401, F403

# Implicitly using functions from C extension:
from Bio import pairwise2

if pairwise2.rint == pairwise2._python_rint:
    from Bio import MissingExternalDependencyError

    raise MissingExternalDependencyError("Missing or non-compiled file: 'cpairwise2'")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
