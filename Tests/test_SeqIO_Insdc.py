# Copyright 2013 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import unittest

from Bio import SeqIO

class TestEmbl(unittest.TestCase):
    def test_annotation1(self):
        """Check parsing of annotation from EMBL files (1)."""
        record = SeqIO.read("EMBL/TRBG361.embl", "embl")
        self.assertEqual(len(record), 1859)
        #Single keyword:
        self.assertEqual(record.annotations["keywords"], ["beta-glucosidase"])

    def test_annotation2(self):
        """Check parsing of annotation from EMBL files (2)."""
        record = SeqIO.read("EMBL/DD231055_edited.embl", "embl")
        self.assertEqual(len(record), 315)
        #Multiple keywords:
        self.assertEqual(record.annotations["keywords"],
                         ['JP 2005522996-A/12', 'test-data',
                          'lot and lots of keywords for this example',
                          'multi-line keywords'])


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
