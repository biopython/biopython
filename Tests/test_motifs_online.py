# Copyright 2015 by Peter Cock.  All rights reserved.
# Copyright 2015 by Antony Lee. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Testing online code for Bio.motifs (weblogo etc)."""
import os
import unittest

import requires_internet
requires_internet.check()

# We want to test these:
from Bio import motifs

# In order to check any sequences returned
from Bio.Seq import Seq


class TestMotifWeblogo(unittest.TestCase):
    def setUp(self):
        self.m = motifs.create([
            Seq("TACAA"), Seq("TACGC"), Seq("TACAC"), Seq("TACCC"),
            Seq("AACCC"), Seq("AATGC"), Seq("AATGC")])

    def test_weblogo(self):
        self.m.weblogo(os.devnull)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
