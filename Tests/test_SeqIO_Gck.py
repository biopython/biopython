# Copyright 2019 Damien Goutte-Gattat.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Tests for the SeqIO Gck module."""

from io import BytesIO
import os.path
import shutil
import unittest
from zipfile import ZipFile

from Bio import SeqIO, MissingExternalDependencyError
from urllib.request import urlopen, HTTPError
import requires_internet


class TestGckWithDGVC(unittest.TestCase):

    sample_data = {
        "pACW": {
            "file": "Drosophila Gateway Vectors GCK/pACW",
            "name": "Construct:",
            "id": "Construct:",
            "description": "Construct:  pACTIN-RW-SV",
            "length": 7957,
            "topology": "circular",
            "features": [
                {
                    "type": "CDS",
                    "start": 6155,
                    "end": 7013,
                    "strand": 1,
                    "label": "ampR",
                },
                {
                    "type": "misc_feature",
                    "start": 5216,
                    "end": 6071,
                    "strand": 1,
                    "label": "SV40 sti/polyA",
                },
                {
                    "type": "misc_feature",
                    "start": 89,
                    "end": 2662,
                    "strand": 1,
                    "label": "actin5C promoter",
                },
                {
                    "type": "CDS",
                    "start": 3722,
                    "end": 4400,
                    "strand": 1,
                    "label": "chlR",
                },
                {
                    "type": "CDS",
                    "start": 4722,
                    "end": 5025,
                    "strand": 1,
                    "label": "ccdB",
                },
                {
                    "type": "misc_feature",
                    "start": 3489,
                    "end": 3507,
                    "strand": 1,
                    "label": "attR1",
                },
                {
                    "type": "misc_feature",
                    "start": 5175,
                    "end": 5193,
                    "strand": -1,
                    "label": "attR2",
                },
                {
                    "type": "misc_feature",
                    "start": 3489,
                    "end": 5193,
                    "strand": 1,
                    "label": "Gateway cassette",
                },
                {
                    "type": "misc_feature",
                    "start": 5192,
                    "end": 5205,
                    "strand": 1,
                    "label": "triple STOP",
                },
                {
                    "type": "CDS",
                    "start": 2763,
                    "end": 3480,
                    "strand": 1,
                    "label": "ECFP",
                },
                {
                    "type": "misc_feature",
                    "start": 2755,
                    "end": 3482,
                    "strand": 1,
                    "label": "pACTIN-SV",
                },
                {
                    "type": "misc_feature",
                    "start": 2755,
                    "end": 3482,
                    "strand": 1,
                    "label": "Construct:  pACTIN-RW-SV",
                },
            ],
        },
        "pPWF": {
            "file": "Drosophila Gateway Vectors GCK/pPWG",
            "name": "Construct:",
            "id": "Construct:",
            "description": "Construct:  pPWF",
            "length": 12320,
            "topology": "circular",
            "features": [
                {
                    "type": "misc_feature",
                    "start": 0,
                    "end": 587,
                    "strand": 1,
                    "label": "P 5' end",
                },
                {
                    "type": "misc_feature",
                    "start": 9327,
                    "end": 9560,
                    "strand": -1,
                    "label": "P 3' end",
                },
                {
                    "type": "misc_feature",
                    "start": 1363,
                    "end": 4244,
                    "strand": -1,
                    "label": "mini-white",
                },
                {
                    "type": "CDS",
                    "start": 10466,
                    "end": 11324,
                    "strand": 1,
                    "label": "ampR",
                },
                {
                    "type": "misc_feature",
                    "start": 7930,
                    "end": 9314,
                    "strand": 1,
                    "label": "K10 terminator",
                },
                {
                    "type": "misc_feature",
                    "start": 4762,
                    "end": 4829,
                    "strand": 1,
                    "label": "GAGA repeats",
                },
                {
                    "type": "misc_feature",
                    "start": 4855,
                    "end": 5177,
                    "strand": 1,
                    "label": "GAL4 sites",
                },
                {
                    "type": "misc_feature",
                    "start": 5279,
                    "end": 5415,
                    "strand": 1,
                    "label": "P intron",
                },
                {
                    "type": "misc_feature",
                    "start": 5184,
                    "end": 5279,
                    "strand": 1,
                    "label": "P promoter",
                },
                {
                    "type": "misc_feature",
                    "start": 4762,
                    "end": 5416,
                    "strand": 1,
                    "label": "UASp promoter",
                },
                {
                    "type": "misc_feature",
                    "start": 10060,
                    "end": 12092,
                    "strand": 1,
                    "label": "pUC8",
                },
                {
                    "type": "misc_feature",
                    "start": 7106,
                    "end": 7124,
                    "strand": -1,
                    "label": "attR2",
                },
                {
                    "type": "CDS",
                    "start": 6653,
                    "end": 6956,
                    "strand": 1,
                    "label": "ccdB",
                },
                {
                    "type": "CDS",
                    "start": 5653,
                    "end": 6331,
                    "strand": 1,
                    "label": "chlR",
                },
                {
                    "type": "misc_feature",
                    "start": 5420,
                    "end": 7124,
                    "strand": 1,
                    "label": "Gateway Cassette",
                },
                {
                    "type": "misc_feature",
                    "start": 5420,
                    "end": 5438,
                    "strand": 1,
                    "label": "attR1",
                },
                {
                    "type": "CDS",
                    "start": 7137,
                    "end": 7854,
                    "strand": 1,
                    "label": "EGFP",
                },
                {
                    "type": "misc_feature",
                    "start": 7129,
                    "end": 7856,
                    "strand": 1,
                    "label": "pACTIN-SV",
                },
                {
                    "type": "misc_feature",
                    "start": 7129,
                    "end": 7856,
                    "strand": 1,
                    "label": "Construct:  pACTIN-WC-SV",
                },
                {
                    "type": "misc_feature",
                    "start": 5416,
                    "end": 7875,
                    "strand": 1,
                    "label": "Construct:  pPWF",
                },
            ],
        },
    }

    def setUp(self):
        # We are using the files of the Drosophila Gateway Vector Collection
        # (<https://emb.carnegiescience.edu/drosophila-gateway-vector-collection>)
        # as sample Gck files. We cannot redistribute those files along with
        # Biopython, so we need to download them now for the tests to run.
        if not os.path.exists("Gck/DGVC_GCK.zip"):
            try:
                requires_internet.check()
            except MissingExternalDependencyError:
                self.skipTest("Sample files missing and no Internet access")
                return

            try:
                with urlopen("https://emb.carnegiescience.edu/sites/default/files/DGVC_GCK.zip") as src, open("Gck/DGVC_GCK.zip", "wb") as dst:
                    shutil.copyfileobj(src, dst)
            except HTTPError:
                self.skipTest("Cannot download the sample files")
                return

        self.zipdata = ZipFile("Gck/DGVC_GCK.zip")

    def tearDown(self):
        self.zipdata.close()

    def test_read(self):
        """Read sample files."""
        for sample in self.sample_data.values():
            with self.zipdata.open(sample["file"]) as f:
                record = SeqIO.read(f, "gck")
            self.assertEqual(sample["name"], record.name)
            self.assertEqual(sample["id"], record.id)
            self.assertEqual(sample["description"], record.description)
            self.assertEqual(sample["length"], len(record))
            self.assertEqual(sample["topology"], record.annotations["topology"])

            self.assertEqual(len(sample["features"]), len(record.features))
            for i, exp_feat in enumerate(sample["features"]):
                read_feat = record.features[i]
                self.assertEqual(exp_feat["type"], read_feat.type)
                self.assertEqual(exp_feat["start"], read_feat.location.start)
                self.assertEqual(exp_feat["end"], read_feat.location.end)
                self.assertEqual(exp_feat["strand"], read_feat.location.strand)
                self.assertEqual(exp_feat["label"], read_feat.qualifiers["label"][0])


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
