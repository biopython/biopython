# Copyright 2019 Damien Goutte-Gattat.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Tests for the SeqIO Gck module."""
import os.path
import shutil
import unittest

from io import BytesIO
from zipfile import ZipFile

from Bio import MissingExternalDependencyError
from Bio import SeqIO


class TestGckWithArtificialData(unittest.TestCase):
    def setUp(self):
        with open("Gck/artificial.gck", "rb") as f:
            self.buffer = f.read()

    def test_read(self):
        """Read an artificial sample file."""
        h = BytesIO(self.buffer)
        record = SeqIO.read(h, "gck")
        self.assertEqual("ACGTACGTACGT", record.seq)
        self.assertEqual("Sample construct", record.description)
        self.assertEqual("linear", record.annotations["topology"])
        self.assertEqual(2, len(record.features))

        self.assertEqual(2, record.features[0].location.start)
        self.assertEqual(6, record.features[0].location.end)
        self.assertEqual(1, record.features[0].location.strand)
        self.assertEqual("misc_feature", record.features[0].type)
        self.assertEqual("FeatureA", record.features[0].qualifiers["label"][0])

        self.assertEqual(7, record.features[1].location.start)
        self.assertEqual(11, record.features[1].location.end)
        self.assertEqual(-1, record.features[1].location.strand)
        self.assertEqual("CDS", record.features[1].type)
        self.assertEqual("FeatureB", record.features[1].qualifiers["label"][0])

        h.close()

    def munge_buffer(self, position, value):
        mod_buffer = bytearray(self.buffer)
        if isinstance(value, list):
            mod_buffer[position : position + len(value) - 1] = value
        else:
            mod_buffer[position] = value
        return BytesIO(mod_buffer)

    def test_conflicting_lengths(self):
        """Read a file with incorrect length."""
        # Change the sequence length as indicated in the sequence packet
        h = self.munge_buffer(0x1C, [0x00, 0x00, 0x20, 0x15])
        with self.assertRaisesRegex(ValueError, "Conflicting sequence length values"):
            SeqIO.read(h, "gck")
        h.close()

        # Change the sequence length as indicated in the features packet
        h = self.munge_buffer(0x36, [0x00, 0x00, 0x20, 0x15])
        with self.assertRaisesRegex(ValueError, "Conflicting sequence length values"):
            SeqIO.read(h, "gck")
        h.close()

        # Change the number of features
        h = self.munge_buffer(0x3B, 0x30)
        with self.assertRaisesRegex(
            ValueError, "Features packet size inconsistent with number of features"
        ):
            SeqIO.read(h, "gck")
        h.close()

        # Change the number of restriction sites
        h = self.munge_buffer(0x137, 0x30)
        with self.assertRaisesRegex(
            ValueError, "Sites packet size inconsistent with number of sites"
        ):
            SeqIO.read(h, "gck")
        h.close()


class TestGckWithImproperHeader(unittest.TestCase):
    def test_read(self):
        """Read a file with an incomplete header."""
        handle = BytesIO(b"tiny")
        with self.assertRaisesRegex(
            ValueError, "Improper header, cannot read 24 bytes from handle"
        ):
            SeqIO.read(handle, "gck")
        handle.close()


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
