# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file included with this package.

"""Regression tests for twofold-site weighting in the LWL85 estimator."""

import unittest

from Bio.Align import Alignment
from Bio.Align.analysis import calculate_dn_ds
from Bio.codonalign.codonseq import CodonSeq
from Bio.codonalign.codonseq import cal_dn_ds


class TestLWL85SiteWeights(unittest.TestCase):
    def check_rates(self, phe, ala, substitutions, expected):
        before = ["TTT"] * phe + ["GCT"] * ala
        after = before.copy()
        for index, codon in substitutions.items():
            after[index] = codon
        sequence1, sequence2 = "".join(before), "".join(after)
        for reverse in (False, True):
            if reverse:
                sequence1, sequence2 = sequence2, sequence1
            for api in ("Align", "codonalign"):
                with self.subTest(api=api, reverse=reverse):
                    if api == "Align":
                        actual = calculate_dn_ds(
                            Alignment([sequence1, sequence2]), method="LWL85"
                        )
                    else:
                        actual = cal_dn_ds(
                            CodonSeq(sequence1), CodonSeq(sequence2), method="LWL85"
                        )
                    self.assertAlmostEqual(actual[0], expected[0], places=12)
                    self.assertAlmostEqual(actual[1], expected[1], places=12)

    def test_unequal_twofold_and_fourfold_counts(self):
        # Phe and Ala have third positions of fold 2 and 4, respectively.
        # TTT -> TTC is a twofold transition; TTT -> TAT a zero-fold
        # transversion. Hence L0=50, L2=20, L4=5, P2=2/20, Q0=1/50.
        # Independent LWL85 formulas give the following dN and dS. In
        # particular dS = -6*log(4/5)/7, not -3*log(4/5)/14.
        self.check_rates(
            20,
            5,
            {0: "TTC", 1: "TTC", 2: "TAT"},
            (0.01603172549117645, 0.19126590112646552),
        )

    def test_twofold_transversion_weight(self):
        # The additional TTT -> TTG is a twofold transversion. The Leu
        # endpoint adds half a twofold site: L0=99/2, L2=41/2, L4=5.
        # Thus P2=4/41, Q2=2/41 and Q0=2/99, testing the dN weight too.
        self.check_rates(
            20,
            5,
            {0: "TTC", 1: "TTC", 2: "TTG", 3: "TAT"},
            (0.03273413610312013, 0.19771625055113906),
        )

    def test_equal_site_counts_control(self):
        # L0=40, L2=L4=10. Equal class counts conceal the weighting bug.
        self.check_rates(
            10,
            10,
            {0: "TTC", 1: "TTC", 2: "TAT"},
            (0.021841909362027918, 0.1915596089122465),
        )


if __name__ == "__main__":
    unittest.main()
