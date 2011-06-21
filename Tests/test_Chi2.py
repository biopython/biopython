# Copyright (C) 2011 by Brandon Invergo (b.invergo@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

import unittest
import os
import os.path
import sys
from Bio.Phylo.PAML import chi2

class ModTest(unittest.TestCase):

    def testCdfChi2(self):
        self.assertRaises(ValueError, chi2.cdf_chi2, df = 0, stat = 3.84)
        self.assertRaises(ValueError, chi2.cdf_chi2, df = 1, stat = -3.84)
        self.assertRaises(TypeError, chi2.cdf_chi2, df = "d", stat ="stat")
        self.assertAlmostEqual(chi2.cdf_chi2(2, 3.84), 0.1466070, places=5)

    def testLnGamma(self):
        self.assertRaises(ValueError, chi2._ln_gamma_function, -1)
        self.assertAlmostEqual(chi2._ln_gamma_function(10), 12.80183,
                               places=5)

    def testIncompleteGamma(self):
        self.assertRaises(ValueError, chi2._incomplete_gamma, x = 0.5, 
                          alpha = -1)
        self.assertAlmostEqual(chi2._incomplete_gamma(0.5, 0.5), 0.6826895,
                               places=5)

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
