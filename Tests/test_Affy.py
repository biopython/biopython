# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for Affy module."""

import unittest

import struct
import os
import sys

try:
    from numpy import array
    import numpy.testing
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.Affy.CelFile") from None

from Bio.Affy import CelFile


class AffyTest(unittest.TestCase):
    def setUp(self):
        self.affy3 = "Affy/affy_v3_example.CEL"
        self.affy4 = "Affy/affy_v4_example.CEL"
        self.affy4Bad = "Affy/affy_v4_bad_example.CEL"
        with open(self.affy4Bad, "wb") as f:
            self.writeExampleV4(f, bad=True)

    def tearDown(self):
        os.remove(self.affy4Bad)

    # tests the old Affymetrix v3 parser
    def testAffy3(self):
        with open(self.affy3, "rb") as f:
            record = CelFile.read(f)
            self.assertTrue(len(record.DatHeader) > 0)
            self.assertEqual(record.intensities.shape, (5, 5))
            self.assertEqual(record.intensities.shape, record.stdevs.shape)
            self.assertEqual(record.intensities.shape, record.npix.shape)
            self.assertEqual(record.ncols, 5)
            self.assertEqual(record.nrows, 5)

    def testAffy3Backwards(self):
        # tests the old Affymetrix v3 parser
        with open(self.affy3, "rb") as f:
            record = CelFile.read(f, version=3)

        self.assertTrue(len(record.DatHeader) > 0)
        self.assertEqual(record.intensities.shape, (5, 5))
        self.assertEqual(record.intensities.shape, record.stdevs.shape)
        self.assertEqual(record.intensities.shape, record.npix.shape)
        self.assertEqual(record.ncols, 5)
        self.assertEqual(record.nrows, 5)

    # tests the new Affymetrix v4 parser
    def testAffy4(self):
        with open(self.affy4, "rb") as f:
            record = CelFile.read(f)
        self.assertEqual(record.intensities.shape, (5, 5))
        self.assertEqual(record.intensities.shape, record.stdevs.shape)
        self.assertEqual(record.intensities.shape, record.npix.shape)
        self.assertEqual(record.ncols, 5)
        self.assertEqual(record.nrows, 5)
        global message
        try:
            numpy.testing.assert_allclose(record.intensities,
                                          [[0., 1., 2., 3., 4.],
                                           [5., 6., 7., 8., 9.],
                                           [10., 11., 12., 13., 14.],
                                           [15., 16., 17., 18., 19.],
                                           [20., 21., 22., 23., 24.]])
            message = None
        except AssertionError as err:
            message = str(err)
        if message is not None:
            self.fail(message)
        try:
            numpy.testing.assert_allclose(record.stdevs,
                                          [[0., -1., -2., -3., -4.],
                                           [-5., -6., -7., -8., -9.],
                                           [-10., -11., -12., -13., -14.],
                                           [-15., -16., -17., -18., -19.],
                                           [-20., -21., -22., -23., -24.]])
            message = None
        except AssertionError as err:
            message = str(err)
        if message is not None:
            self.fail(message)
        try:
            numpy.testing.assert_allclose(record.npix,
                                          [[9, 9, 9, 9, 9],
                                           [9, 9, 9, 9, 9],
                                           [9, 9, 9, 9, 9],
                                           [9, 9, 9, 9, 9],
                                           [9, 9, 9, 9, 9]])
            message = None
        except AssertionError as err:
            message = str(err)
        if message is not None:
            self.fail(message)
        self.assertEqual(len(record.AlgorithmParameters), 329)
        self.assertEqual(len(record.GridCornerUL), 7)
        self.assertEqual(record.AlgorithmParameters[-3:], "169")

    def testAffyBadHeader(self):
        with self.assertRaises(CelFile.ParserError):
            with open(self.affy4Bad, "rb") as f:
                record = CelFile.read(f)

    def testAffyWrongModeReadV3(self):
        with self.assertRaises(ValueError):
            with open(self.affy4, "rt") as f:
                record = CelFile.read(f, version=4)

    def testAffyWrongModeReadV4(self):
        with self.assertRaises(ValueError):
            with open(self.affy4, "rt") as f:
                record = CelFile.read(f, version=4)

    # Writes a small example Affymetrix V4 CEL File
    def writeExampleV4(self, f, bad=False):
        preHeaders = {"cellNo": 25,
                      "columns": 5,
                      "headerLen": 752,
                      "magic": 64,
                      "rows": 5,
                      "version": 4}
        goodH = {"Axis-invertX": b"0"}
        badH = {"Axis-invertX": b"1"}

        headers = {"Algorithm": b"Percentile",
                   "AlgorithmParameters": b"Percentile:75;CellMargin:4;Outlie"
                   b"rHigh:1.500;OutlierLow:1.004;AlgVersion:6.0;FixedCellSize"
                   b":TRUE;FullFeatureWidth:7;FullFeatureHeight:7;IgnoreOutlie"
                   b"rsInShiftRows:FALSE;FeatureExtraction:TRUE;PoolWidthExten"
                   b"stion:1;PoolHeightExtension:1;UseSubgrids:FALSE;Randomize"
                   b"Pixels:FALSE;ErrorBasis:StdvMean;StdMult:1.000000;NumDATS"
                   b"ubgrids:169",
                   "AxisInvertY": b"0",
                   "Cols": b"5",
                   "DatHeader": b"[0..65534]  20_10N:CLS=19420RWS=19420XIN=0"
                   b" YIN=0  VE=30        2.0 05/25/05 23:19:07 50102310  M10 "
                   b"  \x14  \x14 HuEx-1_0-st-v2.1sq \x14  \x14  \x14  \x14  "
                   b"\x14570 \x14 25540.671875 \x14 3.500000 \x14 0.7000 \x14"
                   b" 3",
                   "GridCornerLL": b"518 18668",
                   "GridCornerLR": b"18800 18825",
                   "GridCornerUL": b"659 469",
                   "GridCornerUR": b"18942 623",
                   "OffsetX": b"0",
                   "OffsetY": b"0",
                   "Rows": b"5",
                   "TotalX": b"2560",
                   "TotalY": b"2560",
                   "swapXY": b"0"}
        if not bad:
            headers.update(goodH)
        else:
            headers.update(badH)
        prePadding = b"this text doesn't matter and is ignored\x04"
        preHeadersOrder = ["magic",
                           "version",
                           "columns",
                           "rows",
                           "cellNo",
                           "headerLen"]
        headersEncoded = struct.pack("<" + "i" * len(preHeadersOrder),
                                     *[preHeaders[header] for header in
                                     preHeadersOrder])

        def packData(intensity, sdev, pixel):
            return struct.pack("< f f h", intensity, sdev, pixel)
        f.write(headersEncoded)
        for header in headers:
            try:
                f.write(bytes(header, encoding="utf-8") +
                        b"=" +
                        headers[header] +
                        b"\n")
            except TypeError:
                f.write(header + b"=" + headers[header] + b"\n")
        f.write(prePadding)
        f.write(b"\x00" * 15)
        for i in range(25):
            f.write(packData(float(i), float(-i), 9))


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=0)
    unittest.main(testRunner=runner)
