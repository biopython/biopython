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
        "Install NumPy if you want to use Bio.Affy.CelFile"
    ) from None

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

    # tests the Affymetrix v3 parser
    def testAffy3(self):
        with open(self.affy3) as f:
            record = CelFile.read(f)
            self.assertGreater(len(record.DatHeader), 0)
            self.assertEqual(record.intensities.shape, (5, 5))
            self.assertEqual(record.intensities.shape, record.stdevs.shape)
            self.assertEqual(record.intensities.shape, record.npix.shape)
            self.assertEqual(record.ncols, 5)
            self.assertEqual(record.nrows, 5)
            self.assertEqual(record.version, 3)
            self.assertEqual(record.GridCornerUL, (206, 129))
            self.assertEqual(record.GridCornerUR, (3570, 107))
            self.assertEqual(record.GridCornerLR, (3597, 3470))
            self.assertEqual(record.GridCornerLL, (234, 3492))
            self.assertEqual(record.DatHeader["filename"], "1g_A9AF")
            self.assertEqual(record.DatHeader["CLS"], 3684)
            self.assertEqual(record.DatHeader["RWS"], 3684)
            self.assertEqual(record.DatHeader["XIN"], 1)
            self.assertEqual(record.DatHeader["YIN"], 1)
            self.assertEqual(record.DatHeader["VE"], 30)
            self.assertAlmostEqual(record.DatHeader["laser-power"], 2.0)
            self.assertEqual(record.DatHeader["scan-date"], "08/23/07")
            self.assertEqual(record.DatHeader["scan-time"], "11:23:24")
            self.assertEqual(record.DatHeader["scanner-id"], "50205880")
            self.assertEqual(record.DatHeader["scanner-type"], "M10")
            self.assertEqual(record.DatHeader["array-type"], "Tgondii_SNP1.1sq")
            self.assertEqual(record.DatHeader["filter-wavelength"], 570)
            self.assertAlmostEqual(record.DatHeader["arc-radius"], 25356.509766)
            self.assertAlmostEqual(record.DatHeader["laser-spotsize"], 3.5)
            self.assertAlmostEqual(record.DatHeader["pixel-size"], 1.56)
            self.assertEqual(record.DatHeader["image-orientation"], 6)
            self.assertEqual(record.Algorithm, "Percentile")
            self.assertEqual(len(record.AlgorithmParameters), 16)
            self.assertEqual(record.AlgorithmParameters["Percentile"], 75)
            self.assertEqual(record.AlgorithmParameters["CellMargin"], 2)
            self.assertAlmostEqual(record.AlgorithmParameters["OutlierHigh"], 1.500)
            self.assertAlmostEqual(record.AlgorithmParameters["OutlierLow"], 1.004)
            self.assertEqual(record.AlgorithmParameters["AlgVersion"], "6.0")
            self.assertEqual(
                record.AlgorithmParameters["FixedCellSize"], True
            )  # noqa: A502
            self.assertEqual(record.AlgorithmParameters["FullFeatureWidth"], 7)
            self.assertEqual(record.AlgorithmParameters["FullFeatureHeight"], 7)
            self.assertEqual(
                record.AlgorithmParameters["IgnoreOutliersInShiftRows"], False
            )  # noqa: A502
            self.assertEqual(
                record.AlgorithmParameters["FeatureExtraction"], True
            )  # noqa: A502
            self.assertEqual(record.AlgorithmParameters["PoolWidthExtenstion"], 2)
            self.assertEqual(record.AlgorithmParameters["PoolHeightExtension"], 2)
            self.assertEqual(
                record.AlgorithmParameters["UseSubgrids"], False
            )  # noqa: A502
            self.assertEqual(
                record.AlgorithmParameters["RandomizePixels"], False
            )  # noqa: A502
            self.assertEqual(record.AlgorithmParameters["ErrorBasis"], "StdvMean")
            self.assertAlmostEqual(record.AlgorithmParameters["StdMult"], 1.0)
            self.assertEqual(record.NumberCells, 25)

            global message
            try:
                numpy.testing.assert_allclose(
                    record.intensities,
                    [
                        [234.0, 170.0, 22177.0, 164.0, 22104.0],
                        [188.0, 188.0, 21871.0, 168.0, 21883.0],
                        [188.0, 193.0, 21455.0, 198.0, 21300.0],
                        [188.0, 182.0, 21438.0, 188.0, 20945.0],
                        [193.0, 20370.0, 174.0, 20605.0, 168.0],
                    ],
                )
                message = None
            except AssertionError as err:
                message = str(err)
            if message is not None:
                self.fail(message)
            try:
                numpy.testing.assert_allclose(
                    record.stdevs,
                    [
                        [24.0, 34.5, 2669.0, 19.7, 3661.2],
                        [29.8, 29.8, 2795.9, 67.9, 2792.4],
                        [29.8, 88.7, 2976.5, 62.0, 2914.5],
                        [29.8, 76.2, 2759.5, 49.2, 2762.0],
                        [38.8, 2611.8, 26.6, 2810.7, 24.1],
                    ],
                )
                message = None
            except AssertionError as err:
                message = str(err)
            if message is not None:
                self.fail(message)
            try:
                numpy.testing.assert_array_equal(
                    record.npix,
                    [
                        [25, 25, 25, 25, 25],
                        [25, 25, 25, 25, 25],
                        [25, 25, 25, 25, 25],
                        [25, 25, 25, 25, 25],
                        [25, 25, 25, 25, 25],
                    ],
                )
                message = None
            except AssertionError as err:
                message = str(err)
            if message is not None:
                self.fail(message)
            self.assertEqual(record.nmask, 3)
            try:
                numpy.testing.assert_array_equal(
                    record.mask,
                    [
                        [False, False, False, False, False],
                        [False, False, False, True, True],
                        [False, False, False, False, True],
                        [False, False, False, False, False],
                        [False, False, False, False, False],
                    ],
                )
                message = None
            except AssertionError as err:
                message = str(err)
            if message is not None:
                self.fail(message)
            self.assertEqual(record.noutliers, 3)
            try:
                numpy.testing.assert_array_equal(
                    record.outliers,
                    [
                        [False, False, False, False, False],
                        [False, True, True, False, False],
                        [False, False, False, False, False],
                        [False, True, False, False, False],
                        [False, False, False, False, False],
                    ],
                )
                message = None
            except AssertionError as err:
                message = str(err)
            if message is not None:
                self.fail(message)
            self.assertEqual(record.nmodified, 3)
            try:
                numpy.testing.assert_allclose(
                    record.modified,
                    [
                        [0.0, 0.0, 0.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0, 189.0, 220.0],
                        [0.0, 0.0, 0.0, 21775.0, 0.0],
                        [0.0, 0.0, 0.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0, 0.0, 0.0],
                    ],
                )
                message = None
            except AssertionError as err:
                message = str(err)
            if message is not None:
                self.fail(message)

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
            numpy.testing.assert_allclose(
                record.intensities,
                [
                    [0.0, 1.0, 2.0, 3.0, 4.0],
                    [5.0, 6.0, 7.0, 8.0, 9.0],
                    [10.0, 11.0, 12.0, 13.0, 14.0],
                    [15.0, 16.0, 17.0, 18.0, 19.0],
                    [20.0, 21.0, 22.0, 23.0, 24.0],
                ],
            )
            message = None
        except AssertionError as err:
            message = str(err)
        if message is not None:
            self.fail(message)
        try:
            numpy.testing.assert_allclose(
                record.stdevs,
                [
                    [0.0, -1.0, -2.0, -3.0, -4.0],
                    [-5.0, -6.0, -7.0, -8.0, -9.0],
                    [-10.0, -11.0, -12.0, -13.0, -14.0],
                    [-15.0, -16.0, -17.0, -18.0, -19.0],
                    [-20.0, -21.0, -22.0, -23.0, -24.0],
                ],
            )
            message = None
        except AssertionError as err:
            message = str(err)
        if message is not None:
            self.fail(message)
        try:
            numpy.testing.assert_allclose(
                record.npix,
                [
                    [9, 9, 9, 9, 9],
                    [9, 9, 9, 9, 9],
                    [9, 9, 9, 9, 9],
                    [9, 9, 9, 9, 9],
                    [9, 9, 9, 9, 9],
                ],
            )
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
            with open(self.affy3, "rb") as f:
                record = CelFile.read(f, version=3)

    def testAffyWrongModeReadV4(self):
        with self.assertRaises(ValueError):
            with open(self.affy4) as f:
                record = CelFile.read(f, version=4)

    # Writes a small example Affymetrix V4 CEL File
    def writeExampleV4(self, f, bad=False):
        preHeaders = {
            "cellNo": 25,
            "columns": 5,
            "headerLen": 752,
            "magic": 64,
            "rows": 5,
            "version": 4,
        }
        goodH = {"Axis-invertX": b"0"}
        badH = {"Axis-invertX": b"1"}

        headers = {
            "Algorithm": b"Percentile",
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
            "swapXY": b"0",
        }
        if not bad:
            headers.update(goodH)
        else:
            headers.update(badH)
        prePadding = b"this text doesn't matter and is ignored\x04"
        preHeadersOrder = ["magic", "version", "columns", "rows", "cellNo", "headerLen"]
        headersEncoded = struct.pack(
            "<" + "i" * len(preHeadersOrder),
            *(preHeaders[header] for header in preHeadersOrder),
        )

        f.write(headersEncoded)
        for header in headers:
            try:
                f.write(
                    bytes(header, encoding="utf-8") + b"=" + headers[header] + b"\n"
                )
            except TypeError:
                f.write(header + b"=" + headers[header] + b"\n")
        f.write(prePadding)
        f.write(b"\x00" * 15)
        for i in range(25):
            f.write(struct.pack("< f f h", i, -i, 9))  # intensity, sdev, pixel


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
