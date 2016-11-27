import unittest

from Bio.Affy import CelFile
import struct
from numpy import array
import numpy.testing

class AffyTest(unittest.TestCase):
    def setUp(self):
        self.affy3 = "Affy/affy_v3_example.CEL"
        self.affy4 = "Affy/affy_v4_example.CEL"
    def tearDown(self):
        pass
    
    # tests if the new strict mode complains about passing a non-file
    def testAffyStrict(self):
        with self.assertRaises(IOError) as context:
            with open(self.affy3, "rb") as f:
                record = CelFile.read("hello", strict=True)

    # tests if the new strict mode is backwards compatible
    def testAffyStrict(self):
        with open(self.affy3, "rb") as f:
            record = CelFile.read("hello")
        assert record.DatHeader is None

    # tests the old Affymetrix v3 parser
    def testAffy3(self):
        with open(self.affy3, "r") as f:
            record = CelFile.read(f)
            assert(len(record.DatHeader) > 0)
            assert(record.intensities.shape == (5, 5))
            assert(record.intensities.shape == record.stdevs.shape)
            assert(record.intensities.shape == record.npix.shape)
            assert(record.ncols == 5)
            assert(record.nrows == 5)
    # tests the new Affymetrix v4 parser
    def testAffy4(self):
        with open(self.affy4, "rb") as f:
            record = CelFile.read(f)
            assert(record.intensities.shape == (5, 5))
            assert(record.intensities.shape == record.stdevs.shape)
            assert(record.intensities.shape == record.npix.shape)
            assert(record.ncols == 5)
            assert(record.nrows == 5)
            numpy.allclose(record.intensities, [[  0.,   1.,   2.,   3.,   4.],
                                                [  5.,   6.,   7.,   8.,   9.],
                                                [ 10.,  11.,  12.,  13.,  14.],
                                                [ 15.,  16.,  17.,  18.,  19.],
                                                [ 20.,  21.,  22.,  23.,  24.]])
            numpy.allclose(record.stdevs, [[  0.,  -1.,  -2.,  -3.,  -4.],
                                           [ -5.,  -6.,  -7.,  -8.,  -9.],
                                           [-10., -11., -12., -13., -14.],
                                           [-15., -16., -17., -18., -19.],
                                           [-20., -21., -22., -23., -24.]])
            numpy.allclose(record.npix, [[9, 9, 9, 9, 9],
                                         [9, 9, 9, 9, 9],
                                         [9, 9, 9, 9, 9],
                                         [9, 9, 9, 9, 9],
                                         [9, 9, 9, 9, 9]])

# Writes a small example Affymetrix V4 CEL File
def writeExampleV4():
    preHeaders = {   'cellNo': 25,
        'columns': 5,
        'headerLen': 752,
        'magic': 64,
        'rows': 5,
        'version': 4}
    headers = {   u'Algorithm': u'Percentile',
        u'AlgorithmParameters': u'Percentile:75;CellMargin:4;OutlierHigh:1.500;OutlierLow:1.004;AlgVersion:6.0;FixedCellSize:TRUE;FullFeatureWidth:7;FullFeatureHeight:7;IgnoreOutliersInShiftRows:FALSE;FeatureExtraction:TRUE;PoolWidthExtenstion:1;PoolHeightExtension:1;UseSubgrids:FALSE;RandomizePixels:FALSE;ErrorBasis:StdvMean;StdMult:1.000000;NumDATSubgrids:169',
        u'Axis-invertX': u'0',
        u'AxisInvertY': u'0',
        u'Cols': u'5',
        u'DatHeader': u'[0..65534]  20_10N:CLS=19420RWS=19420XIN=0  YIN=0  VE=30        2.0 05/25/05 23:19:07 50102310  M10   \x14  \x14 HuEx-1_0-st-v2.1sq \x14  \x14  \x14  \x14  \x14 570 \x14 25540.671875 \x14 3.500000 \x14 0.7000 \x14 3',
        u'GridCornerLL': u'518 18668',
        u'GridCornerLR': u'18800 18825',
        u'GridCornerUL': u'659 469',
        u'GridCornerUR': u'18942 623',
        u'OffsetX': u'0',
        u'OffsetY': u'0',
        u'Rows': u'5',
        u'TotalX': u'2560',
        u'TotalY': u'2560',
        u'swapXY': u'0'}
    prePadding = b"this text doesn't matter and is ignored\x04"
    preHeadersOrder = ["magic", "version", "columns", "rows", "cellNo", "headerLen"]
    headersEncoded = struct.pack('<' + 'i' * len(preHeadersOrder), *[preHeaders[header] for header in preHeadersOrder])
    def packData(intensity, sdev, pixel):
        return struct.pack("< f f h", intensity, sdev, pixel)
    with open("exampleFile", "wb") as f:
        f.write(headersEncoded)
        for header in headers:
            f.write(header + "=" + headers[header] + "\n")
        f.write(prePadding)
        f.write("\x00"*15)
        for i in range(25):
            f.write(packData(float(i), float(-i), 9))
            
if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=0)
    unittest.main(testRunner=runner)
