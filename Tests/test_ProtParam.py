# Copyright 2003-2004 by Iddo Friedberg.  All rights reserved.
# Revisions copyright 2008-2010 by Peter Cock. All rights reserved.
# Revisions copyright 2012 by Matt Fenwick. All rights reserved.
# Revisions copyright 2012 by Kai Blin. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Bio.SeqUtils.ProtParam and related code."""

import unittest

from Bio import BiopythonDeprecationWarning
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import molecular_weight
from Bio.SeqUtils import ProtParam
from Bio.SeqUtils import ProtParamData


class ProtParamTest(unittest.TestCase):
    """Tests for ProtParam."""

    def setUp(self):
        """Initialise objects."""
        text = "MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQLSAESVGEVYIKSTETGQYLAMDTSGLLYGSQTPSEECLFLERLEENHYNTYTSKKHAEKNWFVGLKKNGSCKRGPRTHYGQKAILFLPLPV"
        seq = Seq(text)
        record = SeqRecord(seq)
        analysis_text = ProtParam.ProteinAnalysis(text)
        analysis_seq = ProtParam.ProteinAnalysis(seq)
        analysis_record = ProtParam.ProteinAnalysis(record)
        self.text = text
        self.sequences = (text, seq, record)
        self.analyses = (analysis_text, analysis_seq, analysis_record)

    def test_count_amino_acids(self):
        """Calculate amino acid counts."""
        for analysis in self.analyses:
            count_dict = analysis.count_amino_acids()
            for i in count_dict:
                self.assertEqual(count_dict[i], self.text.count(i))

    def test_get_amino_acids_percent(self):
        """Calculate amino acid percentages (DEPRECATED)."""
        with self.assertWarns(BiopythonDeprecationWarning):
            for analysis in self.analyses:
                percent_dict = analysis.get_amino_acids_percent()
                seq_len = len(self.text)
                for i in percent_dict:
                    self.assertAlmostEqual(
                        percent_dict[i], self.text.count(i) / seq_len
                    )

    def test_amino_acids_percent(self):
        """Calculate amino acid percentages."""
        for analysis in self.analyses:
            seq_len = len(self.text)
            for i in analysis.amino_acids_percent:
                self.assertAlmostEqual(
                    analysis.amino_acids_percent[i],
                    (self.text.count(i) * 100 / seq_len),
                )

    def test_get_molecular_weight(self):
        """Calculate protein molecular weight."""
        for analysis in self.analyses:
            self.assertAlmostEqual(analysis.molecular_weight(), 17103.16, 2)

    def test_get_monoisotopic_molecular_weight(self):
        """Calculate monoisotopic molecular weight."""
        for sequence in self.sequences:
            analysis = ProtParam.ProteinAnalysis(sequence, monoisotopic=True)
            self.assertAlmostEqual(analysis.molecular_weight(), 17092.61, 2)

    def test_get_molecular_weight_identical(self):
        """Confirm protein molecular weight agrees with calculation from Bio.SeqUtils."""
        # This test is somehow useless, since ProteinAnalysis.molecular_weight
        # is internally calling SeqUtils.molecular_weight.
        mw_2 = molecular_weight(self.text, seq_type="protein")
        for analysis in self.analyses:
            mw_1 = analysis.molecular_weight()
            self.assertAlmostEqual(mw_1, mw_2)

    def test_get_monoisotopic_molecular_weight_identical(self):
        """Confirm protein molecular weight agrees with calculation from Bio.SeqUtils."""
        # This test is somehow useless, since ProteinAnalysis.molecular_weight
        # is internally calling SeqUtils.molecular_weight.
        mw_2 = molecular_weight(self.text, seq_type="protein", monoisotopic=True)
        for sequence in self.sequences:
            analysis = ProtParam.ProteinAnalysis(sequence, monoisotopic=True)
            mw_1 = analysis.molecular_weight()
            self.assertAlmostEqual(mw_1, mw_2)

    def test_aromaticity(self):
        """Calculate protein aromaticity."""
        for analysis in self.analyses:
            # Old test used a number rounded to two digits, so use the same
            self.assertAlmostEqual(analysis.aromaticity(), 0.10, 2)

    def test_instability_index(self):
        """Calculate protein instability index."""
        for analysis in self.analyses:
            # Old test used a number rounded to two digits, so use the same
            self.assertAlmostEqual(analysis.instability_index(), 41.98, 2)

    def test_flexibility(self):
        """Calculate protein flexibility."""
        # Turn black code style off
        # fmt: off
        expected_flexibility = [
            0.9684642857142857, 0.9118452380952381, 0.8612142857142857,
            0.8256428571428571, 0.8262857142857143, 0.8450833333333334,
            0.8921428571428572, 0.9348690476190478, 0.984202380952381,
            1.021095238095238, 1.0644285714285713, 1.0964880952380953,
            1.1290952380952382, 1.1534880952380953, 1.175452380952381,
            1.188559523809524, 1.2004523809523808, 1.1957261904761904,
            1.1861309523809522, 1.1626547619047618, 1.117190476190476,
            1.0648452380952382, 1.0110238095238095, 0.9854166666666667,
            0.9859880952380952, 1.0114166666666666, 1.0553571428571429,
            1.0887380952380952, 1.0985357142857142, 1.0642619047619046,
            1.0205357142857143, 0.9753809523809525, 0.9456309523809525,
            0.9478333333333334, 0.9782261904761905, 1.011202380952381,
            1.0450833333333334, 1.0662738095238093, 1.0629285714285714,
            1.0645952380952381, 1.062904761904762, 1.079452380952381,
            1.1027142857142858, 1.156642857142857, 1.2014880952380953,
            1.2422857142857144, 1.2526666666666666, 1.2410357142857145,
            1.19075, 1.1423333333333332, 1.076297619047619,
            1.0399642857142857, 1.0037142857142858, 1.0114999999999998,
            1.0185238095238096, 1.0415357142857142, 1.0574285714285714,
            1.0859642857142857, 1.0827619047619048, 1.0870238095238096,
            1.0682857142857143, 1.0475833333333333, 1.0145,
            0.9929285714285715, 0.990952380952381, 0.9903690476190476,
            1.0199047619047619, 1.044607142857143, 1.0721190476190476,
            1.0760952380952378, 1.0725714285714285, 1.0367142857142857,
            1.0014761904761904, 0.9526071428571429, 0.9401428571428573,
            0.9283690476190478, 0.9478690476190476, 0.9693095238095237,
            0.9955238095238095, 1.0019523809523807, 0.9999285714285713,
            1.0046547619047619, 1.0084404761904762, 1.0303333333333333,
            1.067238095238095, 1.1213809523809524, 1.16275,
            1.2009166666666669, 1.2023095238095236, 1.189404761904762,
            1.1374285714285712, 1.0780833333333333, 1.015809523809524,
            0.9733452380952382, 0.9550595238095239, 0.9662738095238097,
            1.013345238095238, 1.0564761904761906, 1.095595238095238,
            1.1098452380952382, 1.1004880952380953, 1.065095238095238,
            1.0202261904761905, 0.9785952380952381, 0.9500595238095239,
            0.951059523809524, 0.966857142857143, 0.9973214285714287,
            1.026690476190476, 1.0583095238095237, 1.0710238095238094,
            1.0807738095238093, 1.0795952380952383, 1.0684285714285713,
            1.0515357142857142, 1.0156547619047618, 0.9936904761904761,
            0.97325, 0.9851785714285713, 1.0112261904761906,
            1.0835833333333333, 1.1376547619047617, 1.1984285714285714,
            1.2216547619047617, 1.2332023809523807, 1.2188928571428572,
            1.2253333333333334, 1.211952380952381, 1.2090476190476191,
            1.186952380952381, 1.1572142857142855, 1.1029880952380953,
            1.0719642857142857, 1.0360833333333335, 1.0257142857142856,
            1.0202261904761905, 1.0217619047619049, 0.9919523809523811,
            0.9644642857142858, 0.9181309523809524, 0.9000952380952382,
            0.8824642857142858
        ]
        # Turn black code style on
        # fmt: on

        for analysis in self.analyses:
            flexibility = analysis.flexibility()
            self.assertEqual(
                len(flexibility), len(expected_flexibility), "Output length differs"
            )
            for f, e in zip(flexibility, expected_flexibility):
                self.assertAlmostEqual(f, e)

    def test_isoelectric_point(self):
        """Calculate the isoelectric point."""
        for analysis in self.analyses:
            # Old test used a number rounded to two digits, so use the same
            self.assertAlmostEqual(analysis.isoelectric_point(), 7.72, 2)

    def test_charge_at_pH(self):
        """Test charge_at_pH function."""
        for analysis in self.analyses:
            self.assertAlmostEqual(analysis.charge_at_pH(7.72), 0.00, 2)

    def test_secondary_structure_fraction(self):
        """Calculate secondary structure fractions."""
        for analysis in self.analyses:
            helix, turn, sheet = analysis.secondary_structure_fraction()
            # Old test used numbers rounded to two digits, so use the same
            self.assertAlmostEqual(helix, 0.33, 2)
            self.assertAlmostEqual(turn, 0.29, 2)
            self.assertAlmostEqual(sheet, 0.37, 2)

    def test_protein_scale(self):
        """Calculate the Kite Doolittle scale."""
        # Turn black code style off
        # fmt: off
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
        # Turn black code style on
        # fmt: on
        for analysis in self.analyses:
            for i, e in zip(analysis.protein_scale(ProtParamData.kd, 9, 0.4), expected):
                # Expected values have 4 decimal places, so restrict to that exactness
                self.assertAlmostEqual(i, e, places=4)

    def test_gravy(self):
        """Calculate gravy. Tests all pre-defined scales."""
        expected_values = {
            "KyteDoolitle": -0.5974,
            "Aboderin": 4.5671,
            "AbrahamLeo": 0.2378,
            "Argos": 0.8607,
            "BlackMould": 0.5074,
            "BullBreese": -0.0445,
            "Casari": -0.2414,
            "Cid": -0.0678,
            "Cowan3.4": 0.0234,
            "Cowan7.5": -0.0733,
            "Eisenberg": -0.0435,
            "Engelman": 1.600,
            "Fasman": -0.3614,
            "Fauchere": 0.327,
            "GoldSack": 1.1564,
            "Guy": 0.0675,
            "Jones": 1.223,
            "Juretic": -0.6672,
            "Kidera": 0.1383,
            "Miyazawa": 5.3109,
            "Parker": 1.7487,
            "Ponnuswamy": 0.3491,
            "Rose": 0.7147,
            "Roseman": -0.4729,
            "Sweet": -0.0791,
            "Tanford": 0.0625,
            "Wilson": 1.5493,
            "Zimmerman": 1.2841,
        }

        for analysis in self.analyses:
            for scale, exp_v in expected_values.items():
                self.assertAlmostEqual(analysis.gravy(scale=scale), exp_v, places=4)

            with self.assertRaises(ValueError) as cm:
                analysis.gravy("Wrong Scale")
            self.assertEqual("scale: Wrong Scale not known", str(cm.exception))

    def test_molar_extinction_coefficient(self):
        """Molar extinction coefficient."""
        for analysis in self.analyses:
            self.assertAlmostEqual(
                analysis.molar_extinction_coefficient()[0], 17420, places=5
            )
            self.assertAlmostEqual(
                analysis.molar_extinction_coefficient()[1], 17545, places=5
            )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
