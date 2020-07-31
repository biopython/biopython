# Copyright 2019 Damien Goutte-Gattat.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Tests for the SeqIO SnapGene module."""

import datetime
from io import BytesIO
import unittest

from Bio import SeqIO


class TestSnapGene(unittest.TestCase):

    sample_data = {
        "sample-d": {
            "file": "SnapGene/sample-d.dna",
            "name": "Sample",
            "id": "Sample",
            "description": "Sample Sequence D",
            "length": 1000,
            "topology": "linear",
            "date": datetime.datetime(2019, 8, 3, 0, 0),
            "features": [
                {
                    "type": "misc_binding",
                    "start": 499,
                    "end": 700,
                    "strand": 1,
                    "label": "FeatureB",
                },
                {
                    "type": "promoter",
                    "start": 49,
                    "end": 150,
                    "strand": 1,
                    "label": "FeatureA",
                },
            ],
        },
        "sample-e": {
            "file": "SnapGene/sample-e.dna",
            "name": "Sample",
            "id": "Sample",
            "description": "Sample Sequence E",
            "length": 1000,
            "date": datetime.datetime(2019, 8, 3, 0, 0),
            "topology": "circular",
            "features": [
                {
                    "type": "terminator",
                    "start": 399,
                    "end": 750,
                    "strand": -1,
                    "label": "FeatureB",
                },
                {
                    "type": "rep_origin",
                    "start": 160,
                    "end": 241,
                    "strand": 1,
                    "label": "FeatureA",
                },
            ],
        },
        "pFA-KanMX4": {
            "file": "SnapGene/pFA-KanMX4.dna",
            "name": "<unknown name>",
            "id": "<unknown id>",
            "description": "<unknown description>",
            "length": 3941,
            "date": datetime.datetime(2020, 7, 30, 0, 0),
            "topology": "circular",
            "features": [
                {
                    "type": "promoter",
                    "start": 0,
                    "end": 3941,
                    "strand": 1,
                    "label": "SP6 promoter",
                },
                {
                    "type": "promoter",
                    "start": 1578,
                    "end": 1597,
                    "strand": -1,
                    "label": "T7 promoter",
                },
                {
                    "type": "promoter",
                    "start": 3474,
                    "end": 3579,
                    "strand": -1,
                    "label": "AmpR promoter",
                },
                {
                    "type": "terminator",
                    "start": 1273,
                    "end": 1471,
                    "strand": 1,
                    "label": "TEF terminator",
                },
                {
                    "type": "promoter",
                    "start": 114,
                    "end": 458,
                    "strand": 1,
                    "label": "TEF promoter",
                },
                {
                    "type": "rep_origin",
                    "start": 1854,
                    "end": 2443,
                    "strand": -1,
                    "label": "ori",
                },
                {
                    "type": "CDS",
                    "start": 458,
                    "end": 1268,
                    "strand": 1,
                    "label": "KanR",
                },
                {
                    "type": "CDS",
                    "start": 2613,
                    "end": 3474,
                    "strand": -1,
                    "label": "AmpR",
                },
                {
                    "type": "gene",
                    "start": 114,
                    "end": 1471,
                    "strand": 1,
                    "label": "kanMX",
                },
            ],
        },
    }

    def test_read(self):
        """Read sample files."""
        for sample in self.sample_data.values():
            record = SeqIO.read(sample["file"], "snapgene")
            self.assertEqual(sample["name"], record.name)
            self.assertEqual(sample["id"], record.id)
            self.assertEqual(sample["description"], record.description)
            self.assertEqual(sample["length"], len(record))
            self.assertEqual(sample["date"], record.annotations["date"])
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


class TestCorruptedSnapGene(unittest.TestCase):
    def setUp(self):
        with open("SnapGene/sample-d.dna", "rb") as f:
            self.buffer = f.read()

    def munge_buffer(self, position, value):
        mod_buffer = bytearray(self.buffer)
        if isinstance(value, list):
            mod_buffer[position : position + len(value) - 1] = value
        else:
            mod_buffer[position] = value
        return BytesIO(mod_buffer)

    def test_invalid_cookie(self):
        """Read a file with missing or invalid cookie packet."""
        # Remove the first packet
        h = BytesIO(self.buffer[19:])
        with self.assertRaisesRegex(
            ValueError, "The file does not start with a SnapGene cookie packet"
        ):
            SeqIO.read(h, "snapgene")
        h.close()

        # Keep the first packet but destroy the magic cookie
        h = self.munge_buffer(5, [0x4B, 0x41, 0x42, 0x4F, 0x4F, 0x4D])
        with self.assertRaisesRegex(
            ValueError, "The file is not a valid SnapGene file"
        ):
            SeqIO.read(h, "snapgene")
        h.close()

    def test_missing_dna(self):
        """Read a file without a DNA packet."""
        # Simulate a missing DNA packet by changing the tag byte to an
        # unknown packet type, so that the parser will skip the packet.
        h = self.munge_buffer(19, 0x80)
        with self.assertRaisesRegex(ValueError, "No DNA packet in file"):
            SeqIO.read(h, "snapgene")
        h.close()

    def test_extra_dna(self):
        """Read a file with supernumerary DNA packet."""
        # Fabricate a file with a duplicated DNA packet
        buf = bytearray(self.buffer)
        buf.extend(self.buffer[19:1025])  # Append duplicated DNA packet
        h = BytesIO(buf)
        with self.assertRaisesRegex(
            ValueError, "The file contains more than one DNA packet"
        ):
            SeqIO.read(h, "snapgene")
        h.close()

    def test_truncated_packet(self):
        """Read a file with incomplete packet."""
        # Truncate before the end of the length bytes
        h = BytesIO(self.buffer[3:])
        with self.assertRaisesRegex(ValueError, "Unexpected end of packet"):
            SeqIO.read(h, "snapgene")
        h.close()

        # Truncate before the end of the data
        h = BytesIO(self.buffer[10:])
        with self.assertRaisesRegex(ValueError, "Unexpected end of packet"):
            SeqIO.read(h, "snapgene")
        h.close()


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
