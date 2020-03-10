# Copyright 2001 by Gavin E. Crooks.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Unit test for Des."""

import unittest

from Bio.SCOP import Des


class DesTests(unittest.TestCase):
    def setUp(self):
        self.filename = "./SCOP/dir.des.scop.txt_test"

    def testParse(self):
        """Test if all records in a DES file are being read."""
        count = 0
        with open(self.filename) as f:
            records = Des.parse(f)
            for record in records:
                count += 1
        self.assertEqual(count, 20)

    def testStr(self):
        """Test if we can convert each record to a string correctly."""
        with open(self.filename) as f:
            for line in f:
                record = Des.Record(line)
                # End of line is platform dependent. Strip it off
                self.assertEqual(str(record).rstrip(), line.rstrip())

    def testError(self):
        """Test if a corrupt record raises the appropriate exception."""
        corruptRec = "49268\tsp\tb.1.2.1\t-\n"
        self.assertRaises(ValueError, Des.Record, corruptRec)

    def testRecord(self):
        """Test one record in detail."""
        recLine = "49268\tsp\tb.1.2.1\t-\tHuman (Homo sapiens)    \n"
        recFields = (49268, "sp", "b.1.2.1", "", "Human (Homo sapiens)")

        record = Des.Record(recLine)
        self.assertEqual(record.sunid, recFields[0])
        self.assertEqual(record.nodetype, recFields[1])
        self.assertEqual(record.sccs, recFields[2])
        self.assertEqual(record.name, recFields[3])
        self.assertEqual(record.description, recFields[4])


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
