
import unittest
from Bio.Affy import CelFile
try:
    import numpy
    from numpy.testing import assert_array_equal
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.Affy.CelFile")

class testCelFile(unittest.TestCase):
    """
    Test reading of CEL file
    """
    def setUp(self):
        self.version = 3
        self.GridCornerUL = (206, 129)
        self.GridCornerUR = (3570, 107)
        self.GridCornerLR = (3597, 3470)
        self.GridCornerLL = (234, 3492)
        self.DatHeader = '[11..65533]  1g_A9AF:CLS=3684 RWS=3684 XIN=1  YIN=1  VE=30        2.0 08/23/07 11:23:24 50205880  M10      Tgondii_SNP1.1sq          570  25356.509766  3.500000  1.5600  6'
        self.Algorithm = 'Percentile'
        self.AlgorithmParameters = 'Percentile:75;CellMargin:2;OutlierHigh:1.500;OutlierLow:1.004;AlgVersion:6.0;FixedCellSize:TRUE;FullFeatureWidth:7;FullFeatureHeight:7;IgnoreOutliersInShiftRows:FALSE;FeatureExtraction:TRUE;PoolWidthExtenstion:2;PoolHeightExtension:2;UseSubgrids:FALSE;RandomizePixels:FALSE;ErrorBasis:StdvMean;StdMult:1.000000'
        self.NumberCells = 25
        self.intensities = numpy.array([[   234.0,    170.0,  22177.0,    164.0,  22104.0],
                                        [   188.0,    188.0,  21871.0,    168.0,  21883.0],
                                        [   188.0,    193.0,  21455.0,    198.0,  21300.0],
                                        [   188.0,    182.0,  21438.0,    188.0,  20945.0],
                                        [   193.0,  20370.0,    174.0,  20605.0,    168.0]])
        self.stdevs = numpy.array([[   24.0, 34.5, 2669.0, 19.7, 3661.2],
                                  [   29.8, 29.8, 2795.9, 67.9, 2792.4],
                                  [   29.8, 88.7, 2976.5, 62.0, 2914.5],
                                  [   29.8, 76.2, 2759.5, 49.2, 2762.0],
                                  [   38.8, 2611.8, 26.6, 2810.7, 24.1]])
        self.npix = numpy.array([[25, 25, 25, 25, 25],
                                  [25, 25, 25, 25, 25],
                                  [25, 25, 25, 25, 25],
                                  [25, 25, 25, 25, 25],
                                  [25, 25, 25, 25, 25]])
        self.nrows = 5
        self.ncols = 5
        self.nmask = 3
        self.mask = numpy.array([[ 0.,  0.,  0.,  0.,  0.],
                                  [ 0.,  0.,  0.,  1.,  1.],
                                  [ 0.,  0.,  0.,  0.,  1.],
                                  [ 0.,  0.,  0.,  0.,  0.],
                                  [ 0.,  0.,  0.,  0.,  0.]])
        self.noutliers = 3
        self.outliers = [[ 0.,  0.,  0.,  0.,  0.],
                          [ 0.,  1.,  1.,  0.,  0.],
                          [ 0.,  0.,  0.,  0.,  0.],
                          [ 0.,  1.,  0.,  0.,  0.],
                          [ 0.,  0.,  0.,  0.,  0.]]
        self.nmodified = 3
        self.modified = [[ 0., 0., 0., 0.,     0.],
                         [ 0., 0., 0., 189.,   220.],
                         [ 0., 0., 0., 21775., 0.],
                         [ 0., 0., 0., 0.,     0.],
                         [ 0., 0., 0., 0.,     0.]]
    
    def test_read(self):
        with open('./Affy/affy_v3_example.CEL') as handle:
            cel = CelFile.read(handle)
            self.assertEqual(self.version, cel.version)
            self.assertEqual(self.GridCornerUL, cel.GridCornerUL)
            self.assertEqual(self.GridCornerUR, cel.GridCornerUR)
            self.assertEqual(self.GridCornerLR, cel.GridCornerLR)
            self.assertEqual(self.GridCornerLL, cel.GridCornerLL)
            self.assertEqual(self.DatHeader, cel.DatHeader)
            self.assertEqual(self.Algorithm, cel.Algorithm)
            self.assertEqual(self.AlgorithmParameters, cel.AlgorithmParameters)
            self.assertEqual(self.NumberCells, cel.NumberCells)
            assert_array_equal(self.intensities, cel.intensities)
            assert_array_equal(self.stdevs, cel.stdevs)
            assert_array_equal(self.npix, cel.npix)
            self.assertEqual(self.nrows, cel.nrows)
            self.assertEqual(self.ncols, cel.ncols)
            self.assertEqual(self.nmask, cel.nmask)
            assert_array_equal(self.mask, cel.mask)
            self.assertEqual(self.noutliers, cel.noutliers)
            assert_array_equal(self.outliers, cel.outliers)
            self.assertEqual(self.nmodified, cel.nmodified)
            assert_array_equal(self.modified, cel.modified)

if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)