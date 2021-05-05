# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Bio.SeqUtils.ctd and related code."""

import unittest
from Bio.SeqUtils import ctd


class CTD_PropertyTest(unittest.TestCase):
    """Tests for CTD_Property class."""

    def test_incomplete_groups(self):
        """Test method with incomplete group definition."""
        with self.assertRaises(ValueError) as cm:
            mock = ctd.CTD_Property(["K", "R", "ANCQGHILMFPSTWYVD"])
        self.assertEqual(
            "Given groups do not contain all Amino Acids.", cm.exception.args[0]
        )


class CTDTest(unittest.TestCase):
    """Tests for CTD class."""

    def setUp(self):
        """Initialise objects."""
        self.seq_text = (
            "MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRD"
            "RSDQHIQLQLSAESVGEVYIKSTETGQYLAMDTSGLLYGSQTPSEECLFLE"
            "RLEENHYNTYTSKKHAEKNWFVGLKKNGSCKRGPRTHYGQKAILFLPLPV"
        )
        self.desc = ctd.CTD(self.seq_text)

    def test_standard_properties(self):
        """Test method with standard parameters."""
        # Turn black code style off
        # fmt: off
        expected = [0.3158, 0.4211, 0.2632, 0.3882, 0.3487, 0.2632, 0.1184,
                    0.7697, 0.1118, 0.3158, 0.3355, 0.3487, 0.3158, 0.4211,
                    0.2632, 0.4408, 0.2697, 0.2895, 0.3816, 0.3158, 0.3026,
                    0.2781, 0.1457, 0.2119, 0.2781, 0.1722, 0.2053, 0.1523,
                    0.0331, 0.1656, 0.1854, 0.1987, 0.2517, 0.2715, 0.1391,
                    0.2384, 0.2318, 0.2252, 0.1457, 0.2384, 0.2583, 0.1854,
                    0.0132, 0.2105, 0.3882, 0.6316, 0.7829, 0.0066, 0.1908,
                    0.4079, 0.5789, 0.7697, 0.0000, 0.1974, 0.3947, 0.6382,
                    0.8224, 0.0066, 0.1974, 0.3421, 0.5461, 0.7697, 0.0132,
                    0.2105, 0.4145, 0.5987, 0.7763, 0.0000, 0.1908, 0.4671,
                    0.7303, 0.8355, 0.0921, 0.1711, 0.4671, 0.7566, 0.8684,
                    0.0000, 0.1974, 0.3947, 0.5855, 0.8026, 0.0132, 0.2763,
                    0.3487, 0.5395, 0.6645, 0.0000, 0.1908, 0.4276, 0.5789,
                    0.8092, 0.0066, 0.1645, 0.4013, 0.5526, 0.7697, 0.0132,
                    0.2303, 0.3882, 0.6645, 0.7829, 0.0066, 0.2171, 0.3487,
                    0.5263, 0.7368, 0.0132, 0.1842, 0.3947, 0.6250, 0.8092,
                    0.0000, 0.1908, 0.4671, 0.7303, 0.8355, 0.0000, 0.1842,
                    0.4079, 0.6447, 0.7829, 0.0329, 0.1908, 0.4276, 0.5789,
                    0.7961, 0.0197, 0.2105, 0.3421, 0.5592, 0.8158, 0.0066,
                    0.1974, 0.3947, 0.5724, 0.8224, 0.0132, 0.2105, 0.3882,
                    0.6316, 0.7829, 0.0000, 0.1908, 0.4211, 0.5789, 0.7368]
        # Turn black code style on
        # fmt: on

        for i, e in zip(self.desc.all_descriptors(), expected):
            self.assertAlmostEqual(i, e, places=4)

        # test partial methods
        for i, e in zip(self.desc.composition_descriptors(), expected[:21]):
            self.assertAlmostEqual(i, e, places=4)

        for i, e in zip(self.desc.transition_descriptors(), expected[21:42]):
            self.assertAlmostEqual(i, e, places=4)

        for i, e in zip(self.desc.distribution_descriptors(), expected[-105:]):
            self.assertAlmostEqual(i, e, places=4)

    def test_different_properties(self):
        """Test using user defined properties."""
        # Turn black code style off
        # fmt: off
        expected = [0.3289, 0.5461, 0.1250, 0.4013, 0.5987, 0.0329, 0.0789,
                    0.0789, 0.1250, 0.4211, 0.2632, 0.3642, 0.0861, 0.1258,
                    0.4901, 0.0000, 0.0000, 0.0199, 0.0331, 0.0132, 0.0132,
                    0.0132, 0.0728, 0.0331, 0.0265, 0.0464, 0.0331, 0.1258,
                    0.0662, 0.2119, 0.0132, 0.1842, 0.4145, 0.6382, 0.7829,
                    0.0066, 0.2105, 0.4013, 0.5855, 0.7895, 0.0000, 0.0987,
                    0.3684, 0.5329, 0.8092, 0.0132, 0.2105, 0.3750, 0.6645,
                    0.7895, 0.0000, 0.1974, 0.4079, 0.5658, 0.7961, 0.2763,
                    0.3026, 0.3289, 0.3487, 0.5395, 0.0132, 0.0855, 0.4408,
                    0.6316, 0.6842, 0.0921, 0.1579, 0.4671, 0.7829, 0.8355,
                    0.1053, 0.2500, 0.3750, 0.6711, 0.8421, 0.0066, 0.1908,
                    0.4079, 0.5789, 0.7697, 0.0000, 0.1974, 0.3947, 0.6382,
                    0.8224]
        # Turn black code style on
        # fmt: on

        # random scales
        new_props = [
            ctd.CTD_Property(["DEKCL", "AGHPSTYNQR", "VIMFW"]),
            ctd.CTD_Property(["DEKNQRHP", "AGSTYCLVIMFW"]),
            ctd.CTD_Property(["D", "E", "K", "NQR", "AGHPSTY", "CLVIMFW"]),
        ]

        for i, e in zip(self.desc.all_descriptors(properties=new_props), expected):
            self.assertAlmostEqual(i, e, places=4)

        # test partial methods
        for i, e in zip(
            self.desc.composition_descriptors(properties=new_props), expected[:11]
        ):
            self.assertAlmostEqual(i, e, places=4)

        for i, e in zip(
            self.desc.transition_descriptors(properties=new_props), expected[11:30]
        ):
            self.assertAlmostEqual(i, e, places=4)

        for i, e in zip(
            self.desc.distribution_descriptors(properties=new_props), expected[-55:]
        ):
            self.assertAlmostEqual(i, e, places=4)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
