#!/usr/bin/env python
# Copyright 2001 by Brad Chapman.  All rights reserved.
# Revisions copyright 2011-2013 by Peter Cock. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Test the Location code in the Bio.SeqFeature module.

This checks to be sure fuzzy and non-fuzzy representations of locations
are working properly.
"""

import unittest

from Bio import SeqFeature


class TestLocations(unittest.TestCase):

    def test_fuzzy(self):
        """Test fuzzy representations."""
        # check the positions alone
        exact_pos = SeqFeature.ExactPosition(5)
        within_pos_s = SeqFeature.WithinPosition(10, left=10, right=13)
        within_pos_e = SeqFeature.WithinPosition(13, left=10, right=13)
        between_pos_e = SeqFeature.BetweenPosition(24, left=20, right=24)
        before_pos = SeqFeature.BeforePosition(15)
        after_pos = SeqFeature.AfterPosition(40)
        self.assertEqual(int(within_pos_s), 10)
        self.assertEqual(str(within_pos_s), "(10.13)")
        self.assertEqual(int(within_pos_e), 13)
        self.assertEqual(str(within_pos_e), "(10.13)")
        self.assertEqual(int(between_pos_e), 24)
        self.assertEqual(str(between_pos_e), "(20^24)")
        self.assertEqual(str(before_pos), "<15")
        self.assertEqual(str(after_pos), ">40")
        # put these into Locations
        location1 = SeqFeature.FeatureLocation(exact_pos, within_pos_e)
        location2 = SeqFeature.FeatureLocation(before_pos, between_pos_e)
        location3 = SeqFeature.FeatureLocation(within_pos_s, after_pos)
        self.assertEqual(str(location1), "[5:(10.13)]")
        self.assertEqual(str(location1.start), "5")
        self.assertEqual(str(location1.end), "(10.13)")
        self.assertEqual(str(location2), "[<15:(20^24)]")
        self.assertEqual(str(location2.start), "<15")
        self.assertEqual(str(location2.end), "(20^24)")
        self.assertEqual(str(location3), "[(10.13):>40]")
        self.assertEqual(str(location3.start), "(10.13)")
        self.assertEqual(str(location3.end), ">40")
        # --- test non-fuzzy representations
        self.assertEqual(location1.nofuzzy_start, 5)
        self.assertEqual(location1.nofuzzy_end, 13)
        self.assertEqual(location2.nofuzzy_start, 15)
        self.assertEqual(location2.nofuzzy_end, 24)
        self.assertEqual(location3.nofuzzy_start, 10)
        self.assertEqual(location3.nofuzzy_end, 40)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
