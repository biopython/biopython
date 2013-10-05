# Copyright 2003-2004 by Iddo Friedberg.  All rights reserved.
# Revisions copyright 2008-2010 by Peter Cock. All rights reserved.
# Revisions copyright 2012 by Matt Fenwick. All rights reserved.
# Revisions copyright 2012 by Kai Blin. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import unittest
from Bio.SeqUtils import ProtParam, ProtParamData


class ProtParamTest(unittest.TestCase):
    def setUp(self):
        self.seq_text = "MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQLSAESVGEVYIKSTETGQYLAMDTSGLLYGSQTPSEECLFLERLEENHYNTYTSKKHAEKNWFVGLKKNGSCKRGPRTHYGQKAILFLPLPV"
        self.analysis = ProtParam.ProteinAnalysis(self.seq_text)

    def test_count_amino_acids(self):
        "Test getting amino acid counts"
        count_dict = self.analysis.count_amino_acids()
        for i in sorted(count_dict):
            self.assertEqual(count_dict[i], self.seq_text.count(i))

    def test_get_amino_acids_percent(self):
        "Test getting amino acid percentage"
        percent_dict = self.analysis.get_amino_acids_percent()
        seq_len = len(self.seq_text)
        for i in sorted(percent_dict):
            self.assertAlmostEqual(percent_dict[i], self.seq_text.count(i)/float(seq_len))

    def test_get_molecular_weight(self):
        "Test calculating protein molecular weight"
        self.assertAlmostEqual(self.analysis.molecular_weight(), 17102.76)

    def test_get_monoisotopic_molecular_weight(self):
        "Test calculating the monoisotopic molecular weight"
        self.analysis = ProtParam.ProteinAnalysis(self.seq_text, monoisotopic=True)
        self.assertAlmostEqual(self.analysis.molecular_weight(), 17092.53)

    def test_aromaticity(self):
        "Test calculating protein aromaticity"
        # Old test used a number rounded to two digits, so use the same
        self.assertEqual(round(self.analysis.aromaticity(), 2), 0.10)

    def test_instability_index(self):
        "Test calculating protein instability index"
        # Old test used a number rounded to two digits, so use the same
        self.assertEqual(round(self.analysis.instability_index(), 2), 41.98)

    def test_flexibility(self):
        "Test calculating protein flexibility"
        flexibility = self.analysis.flexibility()
        expected_flexibility = [
            0.9825119047619049, 1.0166904761904763, 0.9947857142857144,
            0.9660238095238095, 0.9890714285714285, 0.9737261904761906,
            0.9789166666666669, 1.004547619047619, 1.0235357142857144,
            1.0163214285714286, 0.981297619047619, 1.0388809523809523,
            0.9956309523809524, 1.0379047619047619, 1.014654761904762,
            1.015154761904762, 1.0317619047619049, 1.0100833333333334,
            1.0738333333333334, 1.0460952380952382, 1.0333571428571429,
            1.0429761904761905, 0.9842738095238095, 0.9984404761904762,
            0.9814404761904763, 0.9715357142857144, 1.0063690476190477,
            0.988952380952381, 0.9930952380952381, 0.9962619047619047,
            0.9774523809523811, 0.9747857142857144, 0.9701547619047618,
            0.9759404761904762, 0.9515119047619047, 0.9745714285714286,
            1.007642857142857, 1.0024523809523809, 1.0019761904761904,
            1.0053571428571428, 1.000595238095238, 1.0385238095238094,
            1.0090357142857143, 1.0095, 1.0207142857142857, 1.0371071428571428,
            1.0223690476190477, 1.0373809523809523, 1.030095238095238,
            1.0166190476190475, 0.9939404761904762, 0.9935833333333335,
            1.009547619047619, 0.9678095238095237, 1.0027380952380953,
            0.9720595238095241, 1.0215952380952382, 0.996904761904762,
            1.0330238095238098, 1.0140714285714285, 0.9977976190476192,
            1.0319285714285715, 1.016714285714286, 0.9713928571428571,
            0.9921666666666669, 0.9904404761904763, 1.0350238095238096,
            1.0021904761904763, 1.0092738095238096, 1.0462380952380952,
            1.012392857142857, 1.0289761904761903, 1.0108095238095236,
            0.9762619047619049, 0.9885357142857144, 0.9901428571428571,
            0.9795119047619048, 1.014059523809524, 0.9857976190476189,
            1.0114404761904763, 0.9970357142857144, 0.9755238095238097,
            0.9871428571428573, 0.9820952380952382, 1.006095238095238,
            0.9997023809523811, 1.0065238095238094, 1.0149047619047618,
            1.0451071428571428, 1.0396785714285717, 1.0452142857142857,
            1.0262619047619048, 0.9735952380952381, 0.998404761904762,
            0.9777976190476191, 0.9773809523809524, 1.008214285714286,
            0.9772619047619048, 0.9955000000000002, 1.0482261904761907,
            1.0262499999999999, 1.0189404761904763, 0.9988452380952382,
            1.0009285714285714, 1.0204404761904762, 0.9839880952380953,
            0.9809166666666669, 1.01475, 1.0234166666666669, 1.0355476190476192,
            1.033107142857143, 1.011154761904762, 1.0480714285714288,
            1.0591190476190477, 1.033809523809524, 1.008107142857143,
            0.9757261904761906, 0.9885714285714287, 0.9831428571428572,
            1.000595238095238, 0.9739761904761908, 1.0377619047619049,
            1.0295357142857142, 1.0269642857142858, 1.0371904761904762,
            1.0385476190476193, 1.0057023809523808, 1.06075, 1.006714285714286,
            1.0269642857142858, 1.0229761904761905, 1.0046309523809522,
            1.0053690476190476, 0.9886666666666667, 0.9931666666666669, 1.01775,
            1.003297619047619, 1.0161666666666667, 0.977440476190476,
            0.9762738095238096, 0.9785833333333332, 0.9609642857142857,
            0.9650833333333334]

        self.assertEqual(len(flexibility), len(expected_flexibility), "Output length differs")
        for f, e in zip(flexibility, expected_flexibility):
            self.assertAlmostEqual(f, e)

    def test_isoelectric_point(self):
        "Test calculating the isoelectric point"
        # Old test used a number rounded to two digits, so use the same
        self.assertAlmostEqual(round(self.analysis.isoelectric_point(), 2), 7.72)

    def test_secondary_structure_fraction(self):
        "Test calculating secondary structure fractions"
        helix, turn, sheet = self.analysis.secondary_structure_fraction()
        # Old test used numbers rounded to two digits, so use the same
        self.assertAlmostEqual(round(helix, 2), 0.28)
        self.assertAlmostEqual(round(turn, 2), 0.26)
        self.assertAlmostEqual(round(sheet, 2), 0.25)

    def test_protein_scale(self):
        "Test calculating the Kite Doolittle scale"
        expected = [-0.0783, +0.0358, +0.1258, +0.6950, +0.8775, +0.8350, +0.2925, +0.3383,
                    -0.1733, -0.4142, -0.5292, -0.6108, -0.8308, -0.8100, -0.8208, -1.0283,
                    -1.6300, -1.8233, -2.4267, -2.2292, -1.7817, -1.4742, -0.7467, -0.1608,
                    +0.1108, +0.2142, +0.1792, -0.1217, -0.4808, -0.4333, -0.5167, -0.2833,
                    +0.3758, +0.7225, +0.4958, +0.6033, +0.5625, +0.3108, -0.2408, -0.0575,
                    -0.3717, -0.7800, -1.1242, -1.4083, -1.7550, -2.2642, -2.8575, -2.9175,
                    -2.5358, -2.5325, -1.8142, -1.4667, -0.6058, -0.4483, +0.1300, +0.1225,
                    +0.2825, +0.1650, +0.3317, -0.2000, +0.2683, +0.1233, +0.4092, +0.1392,
                    +0.4192, +0.2758, -0.2350, -0.5750, -0.5983, -1.2067, -1.3867, -1.3583,
                    -0.8708, -0.5383, -0.3675, +0.0667, +0.0825, -0.0150, +0.1817, +0.4692,
                    +0.3017, +0.3800, +0.4825, +0.4675, +0.1575, -0.1783, -0.5175, -1.2017,
                    -1.7033, -1.5500, -1.2375, -0.8500, -0.0583, +0.3125, +0.4242, +0.7133,
                    +0.5633, +0.0483, -0.7167, -1.3158, -1.9217, -2.5033, -2.4117, -2.2483,
                    -2.3758, -2.0633, -1.8900, -1.8667, -1.9292, -1.8625, -2.0050, -2.2708,
                    -2.4050, -2.3508, -2.1758, -1.5533, -1.0350, -0.1983, -0.0233, +0.1800,
                    +0.0317, -0.0917, -0.6375, -0.9650, -1.4500, -1.6008, -1.7558, -1.5450,
                    -1.7900, -1.8133, -2.0125, -2.1383, -2.3142, -2.1525, -2.1425, -1.9733,
                    -1.4742, -0.8083, -0.2100, +0.8067, +1.3092, +1.8367, +2.0283, +2.3558]
        for i, e in zip(self.analysis.protein_scale(ProtParamData.kd, 9, 0.4), expected):
            # Expected values have 4 decimal places, so restrict to that exactness
            self.assertAlmostEqual(i, e, places=4)

    def test_gravy(self):
        "Test calculating gravy"
        self.assertAlmostEqual(self.analysis.gravy(), -0.5974, places=4)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
