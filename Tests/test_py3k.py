# Copyright 2009 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for our Python 2/3 compatibility layer, Bio._py3k"""

import unittest


class ODTest(unittest.TestCase):
    def test_od(self):
        """Quick test OrderedDict works."""
        from Bio._py3k import OrderedDict
        d = OrderedDict()
        d[5] = "five"
        d[1] = "one"
        d[3] = "three"
        self.assertEqual(list(d.keys()), [5, 1, 3])

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
