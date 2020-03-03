# Copyright 2008-2016 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Testing online code for fetching NCBI qblast.

This test file imports the TestCases from ``NCBI_qblast_testCases.py``.
If you want to add more tests, do this over there.

These tests will be running online and may take a long time to complete!

"""
import unittest
from io import BytesIO

# Import all test classes deliberately:
from NCBI_qblast_testCases import *  # noqa: F401, F403

import requires_internet


requires_internet.check()

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
