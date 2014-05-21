# Copyright 2014 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import unittest
from Bio import SeqRecord
from Bio.SeqIO import _lazy

class SeqRecordProxyBaseClassTests(unittest.TestCase):

    def setUp(self):
        pass

    def test_nothing(self):
        """An addition test"""
        a = _lazy.SeqRecordProxyBase("sequencefake", "fakeid")
        self.assertEqual(5, 5)

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)