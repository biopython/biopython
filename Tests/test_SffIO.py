# Copyright 2012 by Jeff Hussmann.  All rights reserved.
# Revisions copyright 2013 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import re
import unittest
from io import BytesIO

from Bio import SeqIO

# sffinfo E3MFGYR02_random_10_reads.sff | sed -n '/>\|Run Prefix\|Region\|XY/p'
test_data = """
>E3MFGYR02JWQ7T
  Run Prefix:   R_2008_01_09_16_16_00_
  Region #:     2
  XY Location:  3946_2103
>E3MFGYR02JA6IL
  Run Prefix:   R_2008_01_09_16_16_00_
  Region #:     2
  XY Location:  3700_3115
>E3MFGYR02JHD4H
  Run Prefix:   R_2008_01_09_16_16_00_
  Region #:     2
  XY Location:  3771_2095
>E3MFGYR02GFKUC
  Run Prefix:   R_2008_01_09_16_16_00_
  Region #:     2
  XY Location:  2520_2738
>E3MFGYR02FTGED
  Run Prefix:   R_2008_01_09_16_16_00_
  Region #:     2
  XY Location:  2268_2739
>E3MFGYR02FR9G7
  Run Prefix:   R_2008_01_09_16_16_00_
  Region #:     2
  XY Location:  2255_0361
>E3MFGYR02GAZMS
  Run Prefix:   R_2008_01_09_16_16_00_
  Region #:     2
  XY Location:  2468_1618
>E3MFGYR02HHZ8O
  Run Prefix:   R_2008_01_09_16_16_00_
  Region #:     2
  XY Location:  2958_1574
>E3MFGYR02GPGB1
  Run Prefix:   R_2008_01_09_16_16_00_
  Region #:     2
  XY Location:  2633_0607
>E3MFGYR02F7Z7G
  Run Prefix:   R_2008_01_09_16_16_00_
  Region #:     2
  XY Location:  2434_1658"""


class TestUAN(unittest.TestCase):
    def setUp(self):
        self.records = [record for record in SeqIO.parse('Roche/E3MFGYR02_random_10_reads.sff', 'sff')]
        self.test_annotations = {}
        for line in test_data.splitlines():
            fields = re.split(r"\s+", line.strip())
            if '>' in line:
                current_name = fields[0].lstrip('>')
                self.test_annotations[current_name] = {}
            elif 'Prefix' in line:
                time_list = [int(v) for v in fields[2].split('_')[1:-1]]
                self.test_annotations[current_name]["time"] = time_list
            elif 'Region' in line:
                region = int(fields[-1])
                self.test_annotations[current_name]["region"] = region
            elif 'XY' in line:
                x, y = [int(v) for v in fields[-1].split('_')]
                self.test_annotations[current_name]["coords"] = (x, y)

    def test_time(self):
        for record in self.records:
            self.assertEqual(record.annotations["time"], self.test_annotations[record.name]["time"])

    def test_region(self):
        for record in self.records:
            self.assertEqual(record.annotations["region"], self.test_annotations[record.name]["region"])

    def test_coords(self):
        for record in self.records:
            self.assertEqual(record.annotations["coords"], self.test_annotations[record.name]["coords"])


class TestErrors(unittest.TestCase):
    def test_empty(self):
        fh = BytesIO()
        try:
            records = list(SeqIO.parse(fh, "sff"))
        except ValueError as err:
            self.assertEqual(str(err), "Empty file.")
        else:
            self.assertTrue(False, "Empty file did not raise exception")

    def check_bad_header(self, header, msg):
        try:
            records = list(SeqIO.parse(BytesIO(header), "sff"))
        except ValueError as err:
            if isinstance(msg, (tuple, list)):
                self.assertTrue(str(err) in msg, "Unexpected error: %s" % err)
            else:
                self.assertEqual(str(err), msg)
        else:
            self.assertTrue(False, "Test SFF header only did not raise exception")

    def test_30bytes(self):
        self.check_bad_header(b"x" * 30,
                              "File too small to hold a valid SFF header.")

    def test_31bytes(self):
        self.check_bad_header(b"x" * 31,
                              ("SFF file did not start '.sff', but 'xxxx'",
                               "SFF file did not start '.sff', but b'xxxx'"))

    def test_31bytes_bad_ver(self):
        self.check_bad_header(b".sff1.00" + b"x" * 23,
                              "Unsupported SFF version in header, 49.46.48.48")

    def test_31bytes_bad_flowgram(self):
        self.check_bad_header(b".sff\x00\x00\x00\x01" + b"x" * 23,
                              "Flowgram format code 120 not supported")


class TestConcatenated(unittest.TestCase):
    def test_parses_gzipped_stream(self):
        import gzip
        count = 0
        fh = gzip.open("Roche/E3MFGYR02_random_10_reads.sff.gz", 'rb')
        for record in SeqIO.parse(fh, 'sff'):
            count += 1
        self.assertEqual(10, count)

    def test_parse1(self):
        count = 0
        caught = False
        try:
            for record in SeqIO.parse("Roche/invalid_greek_E3MFGYR02.sff", "sff"):
                count += 1
        except ValueError as err:
            self.assertTrue("Additional data at end of SFF file, perhaps "
                            "multiple SFF files concatenated? "
                            "See offset 65296" in str(err), err)
            caught = True
        self.assertTrue(caught, "Didn't spot concatenation")
        self.assertEqual(count, 24)

    def test_index1(self):
        try:
            d = SeqIO.index("Roche/invalid_greek_E3MFGYR02.sff", "sff")
        except ValueError as err:
            self.assertTrue("Additional data at end of SFF file, perhaps "
                            "multiple SFF files concatenated? "
                            "See offset 65296" in str(err), err)
        else:
            raise ValueError("Indxing Roche/invalid_greek_E3MFGYR02.sff should fail")

    def test_parse2(self):
        count = 0
        caught = False
        try:
            for record in SeqIO.parse("Roche/invalid_paired_E3MFGYR02.sff", "sff"):
                count += 1
        except ValueError as err:
            self.assertTrue("Your SFF file is invalid, post index 5 byte "
                            "null padding region ended '.sff' which could "
                            "be the start of a concatenated SFF file? "
                            "See offset 54371" in str(err), err)
            caught = True
        self.assertTrue(caught, "Didn't spot concatenation")
        self.assertEqual(count, 20)

    def test_index2(self):
        try:
            d = SeqIO.index("Roche/invalid_paired_E3MFGYR02.sff", "sff")
        except ValueError as err:
            self.assertTrue("Your SFF file is invalid, post index 5 byte "
                            "null padding region ended '.sff' which could "
                            "be the start of a concatenated SFF file? "
                            "See offset 54371" in str(err), err)
        else:
            raise ValueError("Indxing Roche/invalid_paired_E3MFGYR02.sff should fail")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
