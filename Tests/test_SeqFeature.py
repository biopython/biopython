# Copyright 2001 by Brad Chapman.  All rights reserved.
# Revisions copyright 2011-2013 by Peter Cock. All rights reserved.
# Copyright 2015-2017 by Kai Blin.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests Bio.SeqFeature."""

import unittest

from os import path

from Bio import Seq
from Bio import SeqIO
from Bio import SeqRecord
from Bio.Data.CodonTable import TranslationError
from Bio.SeqFeature import FeatureLocation, AfterPosition, BeforePosition
from Bio.SeqFeature import CompoundLocation, UnknownPosition, SeqFeature
from Bio.SeqFeature import ExactPosition, WithinPosition, BetweenPosition
from Bio.SeqFeature import OneOfPosition


class TestReference(unittest.TestCase):
    """Tests for the SeqFeature.Reference class."""

    def test_eq_identical(self):
        """Test two identical references eq() to True."""
        testfile = path.join("GenBank", "origin_line.gb")
        rec1 = SeqIO.read(testfile, "genbank")
        rec2 = SeqIO.read(testfile, "genbank")

        self.assertEqual(
            rec1.annotations["references"][0], rec1.annotations["references"][0]
        )
        self.assertEqual(
            rec1.annotations["references"][0], rec2.annotations["references"][0]
        )
        self.assertNotEqual(
            rec1.annotations["references"][0], rec1.annotations["references"][1]
        )
        self.assertNotEqual(
            rec1.annotations["references"][0], rec2.annotations["references"][1]
        )
        self.assertEqual(
            rec1.annotations["references"][1], rec1.annotations["references"][1]
        )
        self.assertEqual(
            rec1.annotations["references"][1], rec2.annotations["references"][1]
        )


class TestFeatureLocation(unittest.TestCase):
    """Tests for the SeqFeature.FeatureLocation class."""

    def test_eq_identical(self):
        """Test two identical locations are equal."""
        loc1 = FeatureLocation(23, 42, 1)
        loc2 = FeatureLocation(23, 42, 1)
        self.assertEqual(loc1, loc2)

        loc1 = FeatureLocation(23, 42, -1)
        loc2 = FeatureLocation(23, 42, -1)
        self.assertEqual(loc1, loc2)

        loc1 = FeatureLocation(BeforePosition(23), AfterPosition(42), 1)
        loc2 = FeatureLocation(23, 42, 1)
        self.assertEqual(loc1, loc2)

        loc1 = FeatureLocation(23, 42, 1, "foo", "bar")
        loc2 = FeatureLocation(23, 42, 1, "foo", "bar")
        self.assertEqual(loc1, loc2)

    def test_eq_not_identical(self):
        """Test two different locations are not equal."""
        loc1 = FeatureLocation(22, 42, 1)
        loc2 = FeatureLocation(23, 42, 1)
        self.assertNotEqual(loc1, loc2)

        loc1 = FeatureLocation(23, 42, 1)
        loc2 = FeatureLocation(23, 43, 1)
        self.assertNotEqual(loc1, loc2)

        loc1 = FeatureLocation(23, 42, 1)
        loc2 = FeatureLocation(23, 42, -1)
        self.assertNotEqual(loc1, loc2)

        loc1 = FeatureLocation(23, 42, 1)
        loc2 = (23, 42, 1)
        self.assertNotEqual(loc1, loc2)

        loc1 = FeatureLocation(23, 42, 1, "foo")
        loc2 = FeatureLocation(23, 42, 1, "bar")
        self.assertNotEqual(loc1, loc2)

        loc1 = FeatureLocation(23, 42, 1, "foo", "bar")
        loc2 = FeatureLocation(23, 42, 1, "foo", "baz")
        self.assertNotEqual(loc1, loc2)

    def test_start_before_end(self):
        expected = "must be greater than or equal to start location"
        with self.assertRaises(ValueError) as err:
            FeatureLocation(42, 23, 1)
        self.assertIn(expected, str(err.exception))

        with self.assertRaises(ValueError) as err:
            FeatureLocation(42, 0, 1)
        self.assertIn(expected, str(err.exception))

        with self.assertRaises(ValueError) as err:
            FeatureLocation(BeforePosition(42), AfterPosition(23), -1)
        self.assertIn(expected, str(err.exception))

        with self.assertRaises(ValueError) as err:
            FeatureLocation(42, AfterPosition(0), 1)
        self.assertIn(expected, str(err.exception))

        # Features with UnknownPositions should pass check
        FeatureLocation(42, UnknownPosition())
        FeatureLocation(UnknownPosition(), 42)

        # Same start and end should pass check
        FeatureLocation(42, 42)


class TestCompoundLocation(unittest.TestCase):
    """Tests for the SeqFeature.CompoundLocation class."""

    def test_eq_identical(self):
        """Test two identical locations are equal."""
        loc1 = FeatureLocation(12, 17, 1) + FeatureLocation(23, 42, 1)
        loc2 = FeatureLocation(12, 17, 1) + FeatureLocation(23, 42, 1)
        self.assertEqual(loc1, loc2)

        loc1 = FeatureLocation(12, 17, 1) + FeatureLocation(23, 42, 1)
        loc2 = CompoundLocation(
            [FeatureLocation(12, 17, 1), FeatureLocation(23, 42, 1)]
        )
        self.assertEqual(loc1, loc2)

    def test_eq_not_identical(self):
        """Test two different locations are not equal."""
        loc1 = FeatureLocation(12, 17, 1) + FeatureLocation(23, 42, 1)
        loc2 = (
            FeatureLocation(12, 17, 1)
            + FeatureLocation(23, 42, 1)
            + FeatureLocation(50, 60, 1)
        )
        self.assertNotEqual(loc1, loc2)

        loc1 = FeatureLocation(12, 17, 1) + FeatureLocation(23, 42, 1)
        loc2 = FeatureLocation(12, 17, -1) + FeatureLocation(23, 42, -1)
        self.assertNotEqual(loc1, loc2)

        loc1 = CompoundLocation(
            [FeatureLocation(12, 17, 1), FeatureLocation(23, 42, 1)]
        )
        loc2 = CompoundLocation(
            [FeatureLocation(12, 17, 1), FeatureLocation(23, 42, 1)], "order"
        )
        self.assertNotEqual(loc1, loc2)

        loc1 = FeatureLocation(12, 17, 1) + FeatureLocation(23, 42, 1)
        loc2 = 5
        self.assertNotEqual(loc1, loc2)


class TestSeqFeature(unittest.TestCase):
    """Tests for the SeqFeature.SeqFeature class."""

    def test_translation_checks_cds(self):
        """Test that a CDS feature is subject to respective checks."""
        seq = Seq.Seq("GGTTACACTTACCGATAATGTCTCTGATGA")
        f = SeqFeature(FeatureLocation(0, 30), type="CDS")
        f.qualifiers["transl_table"] = [11]
        with self.assertRaises(TranslationError):
            f.translate(seq)


class TestLocations(unittest.TestCase):
    def test_fuzzy(self):
        """Test fuzzy representations."""
        # check the positions alone
        exact_pos = ExactPosition(5)
        within_pos_s = WithinPosition(10, left=10, right=13)
        within_pos_e = WithinPosition(13, left=10, right=13)
        between_pos_e = BetweenPosition(24, left=20, right=24)
        before_pos = BeforePosition(15)
        after_pos = AfterPosition(40)
        self.assertEqual(int(within_pos_s), 10)
        self.assertEqual(str(within_pos_s), "(10.13)")
        self.assertEqual(int(within_pos_e), 13)
        self.assertEqual(str(within_pos_e), "(10.13)")
        self.assertEqual(int(between_pos_e), 24)
        self.assertEqual(str(between_pos_e), "(20^24)")
        self.assertEqual(str(before_pos), "<15")
        self.assertEqual(str(after_pos), ">40")
        # put these into Locations
        location1 = FeatureLocation(exact_pos, within_pos_e)
        location2 = FeatureLocation(before_pos, between_pos_e)
        location3 = FeatureLocation(within_pos_s, after_pos)
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


class TestPositions(unittest.TestCase):
    def test_pickle(self):
        """Test pickle behaviour of position instances."""
        # setup
        import pickle

        within_pos = WithinPosition(10, left=10, right=13)
        between_pos = BetweenPosition(24, left=20, right=24)
        oneof_pos = OneOfPosition(1888, [ExactPosition(1888), ExactPosition(1901)])
        # test __getnewargs__
        self.assertEqual(within_pos.__getnewargs__(), (10, 10, 13))
        self.assertEqual(between_pos.__getnewargs__(), (24, 20, 24))
        self.assertEqual(
            oneof_pos.__getnewargs__(),
            (1888, [ExactPosition(1888), ExactPosition(1901)]),
        )
        # test pickle behaviour
        within_pos2 = pickle.loads(pickle.dumps(within_pos))
        between_pos2 = pickle.loads(pickle.dumps(between_pos))
        oneof_pos2 = pickle.loads(pickle.dumps(oneof_pos))
        self.assertEqual(within_pos, within_pos2)
        self.assertEqual(between_pos, between_pos2)
        self.assertEqual(oneof_pos, oneof_pos2)
        self.assertEqual(within_pos._left, within_pos2._left)
        self.assertEqual(within_pos._right, within_pos2._right)
        self.assertEqual(between_pos._left, between_pos2._left)
        self.assertEqual(between_pos._right, between_pos2._right)
        self.assertEqual(oneof_pos.position_choices, oneof_pos2.position_choices)


class TestExtract(unittest.TestCase):
    def test_reference_in_location_record(self):
        """Test location with reference to another record."""
        parent_record = SeqRecord.SeqRecord(seq=Seq.Seq("actg"))
        another_record = SeqRecord.SeqRecord(seq=Seq.Seq("gtcagctac"))
        location = FeatureLocation(5, 8, ref="ANOTHER.7")
        with self.assertRaisesRegexp(
            ValueError,
            r"Feature references another sequence \(ANOTHER\.7\), references mandatory",
        ):
            location.extract(parent_record)
        with self.assertRaisesRegexp(
            ValueError,
            r"Feature references another sequence \(ANOTHER\.7\), not found in references",
        ):
            location.extract(parent_record, references={"SOMEOTHER.2": another_record})
        self.assertEqual(
            location.extract(parent_record, references={"ANOTHER.7": another_record}),
            "cta",
        )

    def test_reference_in_location_sequence(self):
        """Test location with reference to another sequence."""
        parent_sequence = Seq.Seq("actg")
        another_sequence = Seq.Seq("gtcagctac")
        location = FeatureLocation(5, 8, ref="ANOTHER.7")
        self.assertEqual(
            location.extract(
                parent_sequence, references={"ANOTHER.7": another_sequence}
            ),
            "cta",
        )

    def test_reference_in_compound_location_record(self):
        """Test compound location with reference to another record."""
        parent_record = SeqRecord.SeqRecord(Seq.Seq("aaccaaccaaccaaccaa"))
        another_record = SeqRecord.SeqRecord(Seq.Seq("ttggttggttggttggtt"))
        location = FeatureLocation(2, 6) + FeatureLocation(5, 8, ref="ANOTHER.7")
        with self.assertRaisesRegexp(
            ValueError,
            r"Feature references another sequence \(ANOTHER\.7\), references mandatory",
        ):
            location.extract(parent_record)
        with self.assertRaisesRegexp(
            ValueError,
            r"Feature references another sequence \(ANOTHER\.7\), not found in references",
        ):
            location.extract(parent_record, references={"SOMEOTHER.2": another_record})
        self.assertEqual(
            location.extract(
                parent_record, references={"ANOTHER.7": another_record}
            ).seq,
            "ccaatgg",
        )

    def test_reference_in_compound_location_sequence(self):
        """Test compound location with reference to another sequence."""
        parent_sequence = Seq.Seq("aaccaaccaaccaaccaa")
        another_sequence = Seq.Seq("ttggttggttggttggtt")
        location = FeatureLocation(2, 6) + FeatureLocation(5, 8, ref="ANOTHER.7")
        self.assertEqual(
            location.extract(
                parent_sequence, references={"ANOTHER.7": another_sequence}
            ),
            "ccaatgg",
        )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
