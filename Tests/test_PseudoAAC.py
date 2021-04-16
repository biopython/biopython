# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Bio.SeqUtils.PseudoAAC and related code."""

import unittest
from Bio.SeqUtils import PseudoAAC
from Bio.SeqUtils.ProtParamData import Flex, ja, kd


class PseudoAACTest(unittest.TestCase):
    """Tests for PseudoAAC class."""

    def setUp(self):
        """Initialise objects."""
        self.seq_text = "MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQLSAESVGEVYIKSTETGQYLAMDTSGLLYGSQTPSEECLFLERLEENHYNTYTSKKHAEKNWFVGLKKNGSCKRGPRTHYGQKAILFLPLPV"
        self.paac = PseudoAAC.PseudoAAC(self.seq_text)

    def test_standard_params(self):
        """Test method with standard parameters."""
        # Turn black code style off
        # fmt: off
        expected = [0.0377, 0.0188, 0.0314, 0.0753, 0.0377, 0.0879, 0.0314,
                    0.0314, 0.0753, 0.1130, 0.0126, 0.0439, 0.0502, 0.0377,
                    0.0377, 0.0628, 0.0816, 0.0314, 0.0063, 0.0502, 0.0141,
                    0.0156, 0.0161]
        # Turn black code style on
        # fmt: on

        for i, e in zip(self.paac.pseudoAAC(), expected):
            self.assertAlmostEqual(i, e, places=4)

    def test_different_numeric_params(self):
        """Test different values for l_param and weight."""
        # Turn black code style off
        # fmt: off
        expected = [0.0275, 0.0137, 0.0229, 0.0549, 0.0275, 0.0641, 0.0229,
                    0.0229, 0.0549, 0.0824, 0.0092, 0.0320, 0.0366, 0.0275,
                    0.0275, 0.0458, 0.0595, 0.0229, 0.0046, 0.0366, 0.0206,
                    0.0228, 0.0235, 0.0246, 0.0237, 0.0238, 0.0257, 0.0226,
                    0.0245, 0.0216, 0.0246, 0.0222, 0.0244]
        # Turn black code style on
        # fmt: on

        for i, e in zip(self.paac.pseudoAAC(l_param=13, weight=0.1), expected):
            self.assertAlmostEqual(i, e, places=4)

    def test_different_scales(self):
        """Test using user given scales."""
        # Turn black code style off
        # fmt: off
        expected = [0.0387, 0.0194, 0.0323, 0.0775, 0.0387, 0.0904, 0.0323,
                    0.0323, 0.0775, 0.1162, 0.0129, 0.0452, 0.0516, 0.0387,
                    0.0387, 0.0645, 0.0839, 0.0323, 0.0065, 0.0516, 0.0061,
                    0.0062, 0.0065]
        # Turn black code style on
        # fmt: on

        for i, e in zip(self.paac.pseudoAAC(scales=[Flex, ja, kd]), expected):
            self.assertAlmostEqual(i, e, places=4)

        with self.assertRaises(KeyError) as cm:
            del kd["C"]
            self.paac.pseudoAAC(scales=[kd])
        self.assertEqual("scale 0 is missing value for aa: C", cm.exception.args[0])


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
