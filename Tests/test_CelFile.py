# Copyright 2014 by Vincent Davis. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

import unittest

try:
    import numpy
    from numpy.testing import assert_array_equal
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.Affy.CelFile")

from Bio.Affy import CelFile


class testCelFile(unittest.TestCase):
    """Test reading of CEL file."""

    def test_read(self):
        version = 3
        GridCornerUL = (206, 129)
        GridCornerUR = (3570, 107)
        GridCornerLR = (3597, 3470)
        GridCornerLL = (234, 3492)
        DatHeader = '[11..65533]  1g_A9AF:CLS=3684 RWS=3684 XIN=1  YIN=1  VE=30        2.0 08/23/07 11:23:24 50205880  M10      Tgondii_SNP1.1sq          570  25356.509766  3.500000  1.5600  6'
        Algorithm = 'Percentile'
        AlgorithmParameters = 'Percentile:75;CellMargin:2;OutlierHigh:1.500;OutlierLow:1.004;AlgVersion:6.0;FixedCellSize:TRUE;FullFeatureWidth:7;FullFeatureHeight:7;IgnoreOutliersInShiftRows:FALSE;FeatureExtraction:TRUE;PoolWidthExtenstion:2;PoolHeightExtension:2;UseSubgrids:FALSE;RandomizePixels:FALSE;ErrorBasis:StdvMean;StdMult:1.000000'
        NumberCells = 25
        intensities = numpy.array([[234.0, 170.0, 22177.0, 164.0, 22104.0],
                                   [188.0, 188.0, 21871.0, 168.0, 21883.0],
                                   [188.0, 193.0, 21455.0, 198.0, 21300.0],
                                   [188.0, 182.0, 21438.0, 188.0, 20945.0],
                                   [193.0, 20370.0, 174.0, 20605.0, 168.0]])
        stdevs = numpy.array([[24.0, 34.5, 2669.0, 19.7, 3661.2],
                              [29.8, 29.8, 2795.9, 67.9, 2792.4],
                              [29.8, 88.7, 2976.5, 62.0, 2914.5],
                              [29.8, 76.2, 2759.5, 49.2, 2762.0],
                              [38.8, 2611.8, 26.6, 2810.7, 24.1]])
        npix = numpy.array([[25, 25, 25, 25, 25],
                            [25, 25, 25, 25, 25],
                            [25, 25, 25, 25, 25],
                            [25, 25, 25, 25, 25],
                            [25, 25, 25, 25, 25]])
        nrows = 5
        ncols = 5
        nmask = 3
        mask = numpy.array([[0., 0., 0., 0., 0.],
                            [0., 0., 0., 1., 1.],
                            [0., 0., 0., 0., 1.],
                            [0., 0., 0., 0., 0.],
                            [0., 0., 0., 0., 0.]])
        noutliers = 3
        outliers = [[0., 0., 0., 0., 0.],
                    [0., 1., 1., 0., 0.],
                    [0., 0., 0., 0., 0.],
                    [0., 1., 0., 0., 0.],
                    [0., 0., 0., 0., 0.]]
        nmodified = 3
        modified = [[0., 0., 0., 0., 0.],
                    [0., 0., 0., 189., 220.],
                    [0., 0., 0., 21775., 0.],
                    [0., 0., 0., 0., 0.],
                    [0., 0., 0., 0., 0.]]

        with open('./Affy/affy_v3_example.CEL') as handle:
            cel = CelFile.read(handle)
            self.assertEqual(version, cel.version)
            self.assertEqual(GridCornerUL, cel.GridCornerUL)
            self.assertEqual(GridCornerUR, cel.GridCornerUR)
            self.assertEqual(GridCornerLR, cel.GridCornerLR)
            self.assertEqual(GridCornerLL, cel.GridCornerLL)
            self.assertEqual(DatHeader, cel.DatHeader)
            self.assertEqual(Algorithm, cel.Algorithm)
            self.assertEqual(AlgorithmParameters, cel.AlgorithmParameters)
            self.assertEqual(NumberCells, cel.NumberCells)
            assert_array_equal(intensities, cel.intensities)
            assert_array_equal(stdevs, cel.stdevs)
            assert_array_equal(npix, cel.npix)
            self.assertEqual(nrows, cel.nrows)
            self.assertEqual(ncols, cel.ncols)
            self.assertEqual(nmask, cel.nmask)
            assert_array_equal(mask, cel.mask)
            self.assertEqual(noutliers, cel.noutliers)
            assert_array_equal(outliers, cel.outliers)
            self.assertEqual(nmodified, cel.nmodified)
            assert_array_equal(modified, cel.modified)


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
