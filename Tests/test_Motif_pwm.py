# Copyright 2010 by Bartek Wilczynski.  All rights reserved.
# Copyright 2010 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
import unittest

from Bio import Seq
from Bio import Motif

try:
    from Bio.Motif import _pwm
except ImportError:
    #Failed to import the C code extension module
    from Bio import MissingExternalDependencyError
    raise MissingExternalDependencyError(\
        "C module in Bio.Motif not compiled")

class MotifTestsBasic(unittest.TestCase):
    def test_simple(self):
        """Test if Motif C code works."""
        m=Motif.read(open("Motif/SRF.pfm"),"jaspar-pfm")
        result = m.scanPWM(Seq.Seq("ACGTGTGCGTAGTGCGT",m.alphabet))
        self.assertEqual(6, len(result))
        self.assertAlmostEqual(result[0], -29.18363571)
        self.assertAlmostEqual(result[1], -38.3365097)
        self.assertAlmostEqual(result[2], -29.17756271)
        self.assertAlmostEqual(result[3], -38.04542542)
        self.assertAlmostEqual(result[4], -20.3014183)
        self.assertAlmostEqual(result[5], -25.18009186)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
