# Copyright 2009 by James Casbon.  All rights reserved.
# Revisions copyright 2009-2010 by Michiel de Hoon. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for parsing Compass output."""

import os
import unittest

from Bio import Compass


class CompassTest(unittest.TestCase):
    def setUp(self):
        file_dir = os.path.join("Compass")
        self.test_files = [
            os.path.join(file_dir, "comtest1"),
            os.path.join(file_dir, "comtest2"),
        ]

    def testCompassScanAndConsume(self):
        with open(self.test_files[0]) as handle:
            com_record = Compass.read(handle)

        self.assertEqual("60456.blo.gz.aln", com_record.query)
        self.assertEqual("60456.blo.gz.aln", com_record.hit)
        self.assertAlmostEqual(0.5, com_record.gap_threshold)

        self.assertEqual(388, com_record.query_length)
        self.assertEqual(386, com_record.query_filtered_length)
        self.assertEqual(388, com_record.hit_length)
        self.assertEqual(386, com_record.hit_filtered_length)

        self.assertEqual(399, com_record.query_nseqs)
        self.assertAlmostEqual(12.972, com_record.query_neffseqs)
        self.assertEqual(399, com_record.hit_nseqs)
        self.assertAlmostEqual(12.972, com_record.hit_neffseqs)

        self.assertEqual(2759, com_record.sw_score)
        self.assertAlmostEqual(0.0, com_record.evalue)

    def testCompassParser(self):
        with open(self.test_files[0]) as handle:
            com_record = Compass.read(handle)

        self.assertEqual("60456.blo.gz.aln", com_record.query)

    def testCompassIteratorEasy(self):
        with open(self.test_files[0]) as handle:
            records = Compass.parse(handle)
            com_record = next(records)
        self.assertEqual("60456.blo.gz.aln", com_record.query)
        self.assertRaises(StopIteration, next, records)

    def testCompassIteratorHard(self):
        with open(self.test_files[1]) as handle:
            records = Compass.parse(handle)

            com_record = next(records)
            self.assertEqual("allscop//14982.blo.gz.aln", com_record.hit)
            self.assertAlmostEqual(1.01e03, com_record.evalue)

            com_record = next(records)
            self.assertEqual("allscop//14983.blo.gz.aln", com_record.hit)
            self.assertAlmostEqual(1.01e03, com_record.evalue)

            com_record = next(records)
            self.assertEqual("allscop//14984.blo.gz.aln", com_record.hit)
            self.assertAlmostEqual(5.75e02, com_record.evalue)

    def testAlignmentParsingOne(self):
        with open(self.test_files[1]) as handle:
            records = Compass.parse(handle)

            com_record = next(records)
            self.assertEqual(178, com_record.query_start)
            self.assertEqual("KKDLEEIAD", com_record.query_aln)
            self.assertEqual(9, com_record.hit_start)
            self.assertEqual("QAAVQAVTA", com_record.hit_aln)
            self.assertEqual("++ ++++++", com_record.positives)

            com_record = next(records)
            com_record = next(records)
            self.assertEqual(371, com_record.query_start)
            self.assertEqual("LEEAMDRMER~~~V", com_record.query_aln)
            self.assertEqual(76, com_record.hit_start)
            self.assertEqual("LQNFIDQLDNpddL", com_record.hit_aln)
            self.assertEqual("+ ++++ + +   +", com_record.positives)

    def testAlignmentParsingTwo(self):
        with open(self.test_files[0]) as handle:
            records = Compass.parse(handle)
            com_record = next(records)
        self.assertEqual(2, com_record.query_start)
        self.assertEqual(2, com_record.hit_start)
        self.assertEqual("LKERKL", com_record.hit_aln[-6:])


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
