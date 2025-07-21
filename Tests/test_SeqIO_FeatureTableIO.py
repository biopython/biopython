"""Tests for SeqIO FeatureTable module."""

from io import StringIO
import unittest

from Bio import BiopythonWarning, BiopythonParserWarning
from Bio import SeqIO


class TestRead(unittest.TestCase):
    def test_read_tbl1(self):
        """Test parsing valid .tbl files."""
        # Documentation gives example1.tbl as the equivalent of example1.gb, so
        # make sure they share every feature they should.
        records = list(SeqIO.parse("FeatureTable/example1.tbl", "feature-table"))
        self.assertEqual(len(records), 1)
        tbl_record = records[0]
        tbl_record_i = 0
        with self.assertWarns(BiopythonParserWarning):
            gb_record = next(SeqIO.parse("FeatureTable/example1.gb", "genbank"))
        gb_record_i = 0
        self.assertEqual(tbl_record.id, gb_record.id)
        while tbl_record_i < len(tbl_record.features):
            tbl_feature = tbl_record.features[tbl_record_i]
            gb_feature = gb_record.features[gb_record_i]

            if gb_feature.type == "source":
                # Feature tables don't have source features
                gb_record_i += 1
                continue
            self.assertEqual(tbl_feature.type, gb_feature.type)
            self.assertEqual(tbl_feature.location, gb_feature.location)
            for qualifier in gb_feature.qualifiers:
                if qualifier == "translation":
                    # Feature tables don't have translation qualifiers
                    continue
                self.assertEqual(
                    tbl_feature.qualifiers[qualifier], gb_feature.qualifiers[qualifier]
                )
            tbl_record_i += 1
            gb_record_i += 1

    def test_read_tbl2(self):
        """Test parsing multiple records."""
        records = list(SeqIO.parse("FeatureTable/example2.tbl", "feature-table"))
        self.assertEqual(len(records), 4)
        self.assertEqual(len(records[0].features), 1)
        self.assertEqual(len(records[1].features), 2)
        self.assertEqual(len(records[2].features), 4)
        self.assertEqual(len(records[3].features), 2)

    def test_read_tbl3(self):
        """Test parsing .tbl files not from documentation."""
        with self.assertWarns(BiopythonParserWarning):
            records = list(SeqIO.parse("FeatureTable/example3.tbl", "feature-table"))
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].id, "Chr_1")
        self.assertNotIn("references", records[0].annotations)
        self.assertEqual(len(records[0].features), 6)


class TestWrite(unittest.TestCase):
    def test_write_tbl1(self):
        """Test writing valid .tbl files."""
        handle = StringIO()
        records = list(SeqIO.parse("FeatureTable/example1.tbl", "feature-table"))
        self.assertEqual(len(records), 1)
        SeqIO.write(records, handle, "feature-table")
        handle.seek(0)
        records2 = list(SeqIO.parse(handle, "feature-table"))
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].id, records2[0].id)
        self.assertEqual(records[0].annotations, records2[0].annotations)
        self.assertEqual(records[0].features, records2[0].features)

    def test_write_tbl2(self):
        """Test converting GenBank to .tbl."""
        handle = StringIO()
        with self.assertWarns(BiopythonParserWarning):
            gb_records = list(SeqIO.parse("FeatureTable/example1.gb", "genbank"))
        with self.assertWarns(BiopythonParserWarning):
            SeqIO.write(gb_records, handle, "feature-table")
        handle.seek(0)
        records = list(SeqIO.parse(handle, "feature-table"))
        self.assertEqual(len(records), 1)
        records2 = list(SeqIO.parse("FeatureTable/example1.tbl", "feature-table"))
        self.assertEqual(len(records2), 1)
        self.assertEqual(records[0].id, records2[0].id)
        self.assertEqual(records[0].annotations, records2[0].annotations)
        self.assertEqual(records[0].features, records2[0].features)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
