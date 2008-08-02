#!/usr/bin/env python
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import sys
import unittest

from Bio import Medline

class TestMedline(unittest.TestCase):

    def test_read(self):
        handle = open("Medline/pubmed_result1.txt")
        parser = Medline.RecordParser()
        record = parser.parse(handle)

    def test_parse(self):
        handle = open("Medline/pubmed_result2.txt")
        parser = Medline.RecordParser()
        records = Medline.Iterator(handle, parser)
        for record in records:
            pass


def run_tests(argv):
    suite = unittest.TestLoader().loadTestsFromTestCase(TestMedline)
    runner = unittest.TextTestRunner(sys.stdout, verbosity = 2)
    runner.run(suite)


if __name__ == "__main__":
    sys.exit(run_tests(sys.argv))
