# Copyright 2010 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Testing online code in Bio.SCOP
"""

import unittest

import requires_internet
requires_internet.check()

from Bio import SCOP

class ScopSearch(unittest.TestCase):
    """SCOP search tests."""
    def test_search(self):
        """Bio.SCOP.search(...)"""
        handle = SCOP.search("1JOY")

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
