# Copyright 2019 Damien Goutte-Gattat.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Tests for the SeqIO Xdna module."""
import unittest

from io import BytesIO

from Bio import BiopythonWarning
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import BeforePosition
from Bio.SeqFeature import SimpleLocation
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord


class TestXdna(unittest.TestCase):

    sample_data = {
        "sample-a": {
            "file": "Xdna/sample-a.xdna",
            "name": "Sample",
            "id": "Sample",
            "description": "Sample sequence A",
            "length": 1000,
            "molecule_type": "DNA",
            "topology": "linear",
            "features": [
                {
                    "type": "promoter",
                    "start": 49,
                    "end": 150,
                    "strand": 1,
                    "label": "FeatureA",
                },
                {
                    "type": "misc_binding",
                    "start": 499,
                    "end": 700,
                    "strand": -1,
                    "label": "FeatureB",
                },
            ],
        },
        "sample-b": {
            "file": "Xdna/sample-b.xdna",
            "name": "Sample",
            "id": "Sample",
            "description": "Sample sequence B",
            "length": 1000,
            "molecule_type": "DNA",
            "topology": "circular",
            "features": [
                {
                    "type": "rep_origin",
                    "start": 160,
                    "end": 241,
                    "strand": 1,
                    "label": "FeatureA",
                },
                {
                    "type": "terminator",
                    "start": 399,
                    "end": 750,
                    "strand": -1,
                    "label": "FeatureB",
                },
            ],
        },
        "sample-c": {
            "file": "Xdna/sample-c.xprt",
            "name": "Sample",
            "id": "Sample",
            "description": "Sample Sequence C",
            "length": 1000,
            "molecule_type": "protein",
            "topology": "linear",
            "features": [
                {
                    "type": "misc_feature",
                    "start": 10,
                    "end": 11,
                    "strand": 1,
                    "label": "S11",
                },
                {
                    "type": "misc_binding",
                    "start": 164,
                    "end": 195,
                    "strand": 1,
                    "label": "RIP1",
                },
            ],
        },
    }

    def test_read(self):
        """Read sample files."""
        for sample in self.sample_data.values():
            record = SeqIO.read(sample["file"], "xdna")
            self.assertEqual(sample["name"], record.name)
            self.assertEqual(sample["id"], record.id)
            self.assertEqual(sample["description"], record.description)
            self.assertEqual(sample["length"], len(record))
            self.assertEqual(
                sample["molecule_type"], record.annotations["molecule_type"]
            )
            self.assertEqual(sample["topology"], record.annotations["topology"])

            self.assertEqual(len(sample["features"]), len(record.features))
            for i in range(len(sample["features"])):
                exp_feat = sample["features"][i]
                read_feat = record.features[i]
                self.assertEqual(exp_feat["type"], read_feat.type)
                self.assertEqual(exp_feat["start"], read_feat.location.start)
                self.assertEqual(exp_feat["end"], read_feat.location.end)
                self.assertEqual(exp_feat["strand"], read_feat.location.strand)
                self.assertEqual(exp_feat["label"], read_feat.qualifiers["label"][0])


class TestInvalidXdna(unittest.TestCase):
    def setUp(self):
        with open("Xdna/sample-a.xdna", "rb") as f:
            self.buffer = f.read()

    def munge_buffer(self, position, value):
        mod_buffer = bytearray(self.buffer)
        if isinstance(value, list):
            mod_buffer[position : position + len(value) - 1] = value
        else:
            mod_buffer[position] = value
        return BytesIO(mod_buffer)

    def test_unsupported_version(self):
        """Read a file with unexpected version number."""
        h = self.munge_buffer(0, 0x01)  # Change version byte
        with self.assertRaisesRegex(ValueError, "Unsupported XDNA version"):
            SeqIO.read(h, "xdna")
        h.close()

    def test_invalid_sequence_type(self):
        """Read a file with an unknown sequence type."""
        h = self.munge_buffer(1, 0x0A)  # Change type byte
        with self.assertRaisesRegex(ValueError, "Unknown sequence type"):
            SeqIO.read(h, "xdna")
        h.close()

    def test_corrupted_length(self):
        """Read a file with incorrect length."""
        # Set a length shorter than the actual length of the sequence
        h = self.munge_buffer(29, [0x00, 0x00, 0x00, 0x80])
        with self.assertRaisesRegex(ValueError, "invalid literal"):
            SeqIO.read(h, "xdna")
        h.close()

        # Set a length larger than the actual length of the sequence
        h = self.munge_buffer(29, [0x00, 0x08, 0x00, 0x00])
        with self.assertRaisesRegex(ValueError, "Cannot read 2048 bytes from handle"):
            SeqIO.read(h, "xdna")
        h.close()

    def test_missing_features(self):
        """Read a file with an incorrect number of features."""
        # Set a larger number of features than the file actually contains
        # Offset of the features number byte:
        # header + length of sequence + length of comment + overhangs
        feature_byte = 112 + 1000 + len("Sample sequence A") + 5
        h = self.munge_buffer(feature_byte, 3)
        with self.assertRaisesRegex(ValueError, "Cannot read 1 bytes from handle"):
            SeqIO.read(h, "xdna")
        h.close()


class TestXdnaWriter(unittest.TestCase):
    def test_write_sequence_type(self):
        """Write correct sequence type."""
        h = BytesIO()

        record = SeqRecord(Seq("ACGT"))

        for molecule_type, expected_byte in [
            (None, 0),
            ("DNA", 1),
            ("RNA", 3),
            ("protein", 4),
        ]:
            record.annotations["molecule_type"] = molecule_type
            h.seek(0, 0)
            SeqIO.write([record], h, "xdna")
            buf = bytearray(h.getvalue())
            self.assertEqual(expected_byte, buf[1])

        h.close()

    def test_warnings_on_data_loss(self):
        """Emit warnings when dropping data on write."""
        h = BytesIO()

        # Fabricate a record with > 255 features
        record = SeqRecord(Seq("ACGT"))
        for i in range(260):
            feature = SeqFeature(SimpleLocation(1, 2), type="misc_feature")
            record.features.append(feature)
        with self.assertWarnsRegex(BiopythonWarning, "Too many features"):
            SeqIO.write([record], h, "xdna")

        # Now a record with a fuzzy-located feature
        feature = SeqFeature(SimpleLocation(BeforePosition(2), 3), type="misc_feature")
        record.features = [feature]
        with self.assertWarnsRegex(
            BiopythonWarning, r"Dropping \d+ features with fuzzy locations"
        ):
            SeqIO.write([record], h, "xdna")

        # Now a record with a feature with a qualifier too long
        qualifiers = {"note": ["x" * 260]}
        feature = SeqFeature(
            SimpleLocation(2, 3), type="misc_feature", qualifiers=qualifiers
        )
        record.features = [feature]
        with self.assertWarnsRegex(
            BiopythonWarning, "Some annotations were truncated to 255 characters"
        ):
            SeqIO.write([record], h, "xdna")

        h.close()


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
