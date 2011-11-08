# Copyright 2010-2011 by Peter Cock.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Test code for working with BGZF files (used in BAM files).

See also the doctests in bgzf.py which are called via run_tests.py
"""

import unittest
import gzip
import os

from Bio.bgzf import BgzfWriter, BgzfBlocks

class BgzfTests(unittest.TestCase):
    def setUp(self):
        self.temp_file = "temp.bgzf"
        if os.path.isfile(self.temp_file):
            os.remove(self.temp_file)

    def tearDown(self):
        if os.path.isfile(self.temp_file):
            os.remove(self.temp_file)

    def test_one(self):
        """Reproduce BGZF compression for BAM file"""
        temp_file = self.temp_file

        #Note this example is from an old version of samtools
        #and all the blocks are full (except the last one)
        h = gzip.open("SamBam/ex1.bam", "rb")
        data = h.read()
        h.close()

        h = BgzfWriter(temp_file, "wb")
        h.write(data)
        h.flush()
        h.flush() #Second flush gives empty BGZF block as BAM EOF marker
        h.close()

        h = gzip.open(temp_file)
        new_data = h.read()
        h.close()

        #Check the decompressed files agree
        self.assert_(new_data, "Empty BGZF file?")
        self.assertEqual(len(data), len(new_data))
        self.assertEqual(data, new_data)

        #Now check the blocks agree (using the fact that
        #this example BAM file has simple block usage)
        h = open("SamBam/ex1.bam", "rb")
        old = list(BgzfBlocks(h))
        h.close()
        h = open(temp_file, "rb")
        new = list(BgzfBlocks(h))
        h.close()
        self.assertEqual(len(old), len(new))
        self.assertEqual(old, new)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
