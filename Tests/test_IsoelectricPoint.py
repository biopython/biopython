# Copyright 2026 by Martin Mokrejs.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for IsoelectricPoint multi-method support.

Tests the 18 pKa scale methods added to Bio.SeqUtils.IsoelectricPoint
and the corresponding parameter forwarding in Bio.SeqUtils.ProtParam.
"""

import unittest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint as IP
from Bio.SeqUtils.IsoelectricPoint import pKa_scales
from Bio.SeqUtils.ProtParam import ProteinAnalysis as PA

# Reference protein from the original Biopython ProtParam test suite
REF_PROTEIN = (
    "MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGT"
    "RDRSDQHIQLQLSAESVGEVYIKSTETGQYLAMDTSGLLYGSQTPSEEC"
    "LFLERLEENHYNTYTSKKHAEKNWFVGLKKNGSCKRGPRTHYGQKAILF"
    "LPLPV"
)

# Shorter test sequences
BASIC_PEPTIDE = "INGAR"
ACIDIC_PEPTIDE = "EFYTVDGQK"
LONG_ACIDIC = "ADDKNPLEECFRETDYEEFLEIARNGLKATSNPKRVV"


class TestPKaScalesCatalog(unittest.TestCase):
    """Verify the pKa_scales dictionary structure."""

    def test_all_18_scales_present(self):
        """All 18 pKa scales must be present."""
        expected = {
            "IPC2_protein",
            "IPC2_peptide",
            "IPC_protein",
            "IPC_peptide",
            "Bjellqvist",
            "EMBOSS",
            "DTASelect",
            "Solomon",
            "Sillero",
            "Rodwell",
            "Patrickios",
            "Wikipedia",
            "Grimsley",
            "Lehninger",
            "Toseland",
            "Thurlkill",
            "Nozaki",
            "Dawson",
        }
        self.assertEqual(set(pKa_scales.keys()), expected)

    def test_scale_structure(self):
        """Each scale has [neg_pKs, pos_pKs] with correct keys."""
        neg_keys = {"Cterm", "D", "E", "C", "Y"}
        pos_keys = {"H", "Nterm", "K", "R"}
        for name, (neg, pos) in pKa_scales.items():
            with self.subTest(scale=name):
                self.assertEqual(set(neg.keys()), neg_keys, f"{name} neg keys")
                self.assertEqual(set(pos.keys()), pos_keys, f"{name} pos keys")

    def test_pka_values_are_float(self):
        """All pKa values must be numeric."""
        for name, (neg, pos) in pKa_scales.items():
            for key, val in {**neg, **pos}.items():
                with self.subTest(scale=name, key=key):
                    self.assertIsInstance(val, (int, float))


class TestBackwardCompatibility(unittest.TestCase):
    """Ensure the default Bjellqvist behavior is unchanged."""

    def test_default_is_bjellqvist(self):
        """Default pKa_scale should be Bjellqvist."""
        p = IP(REF_PROTEIN)
        self.assertEqual(p.pKa_scale, "Bjellqvist")

    def test_default_pi_unchanged(self):
        """Default pI matches the historic Biopython value."""
        p = IP(REF_PROTEIN)
        self.assertAlmostEqual(p.pi(), 7.72, places=2)

    def test_default_charge_at_pi_near_zero(self):
        """Charge at the isoelectric point should be near zero."""
        p = IP(REF_PROTEIN)
        pi_val = p.pi()
        self.assertAlmostEqual(p.charge_at_pH(pi_val), 0.0, places=2)

    def test_legacy_constants_still_exported(self):
        """Module-level positive_pKs, negative_pKs etc. still available."""
        from Bio.SeqUtils.IsoelectricPoint import (
            charged_aas,
            negative_pKs,
            pKcterminal,
            pKnterminal,
            positive_pKs,
        )

        self.assertEqual(positive_pKs["K"], 10.0)
        self.assertEqual(negative_pKs["D"], 4.05)
        self.assertEqual(pKcterminal["D"], 4.55)
        self.assertEqual(pKnterminal["A"], 7.59)
        self.assertIn("K", charged_aas)

    def test_ingar_pi_unchanged(self):
        """The INGAR peptide pI from the original docstring."""
        p = IP("INGAR")
        self.assertAlmostEqual(p.pi(), 9.75, places=2)

    def test_ingar_charge_unchanged(self):
        """The INGAR charge at pH 7 from the original docstring."""
        p = IP("INGAR")
        self.assertAlmostEqual(p.charge_at_pH(7.0), 0.76, places=2)


class TestBjellqvistTerminalCorrections(unittest.TestCase):
    """Test Bjellqvist-specific N- and C-terminal pKa corrections."""

    def test_nterm_methionine(self):
        """M at N-terminus gets adjusted pKa (7.0)."""
        p = IP("MAAAA")
        self.assertAlmostEqual(p.pos_pKs["Nterm"], 7.0)

    def test_nterm_alanine(self):
        """A at N-terminus gets adjusted pKa (7.59)."""
        p = IP("AAAAA")
        self.assertAlmostEqual(p.pos_pKs["Nterm"], 7.59)

    def test_nterm_serine(self):
        """S at N-terminus gets adjusted pKa (6.93)."""
        p = IP("SAAAA")
        self.assertAlmostEqual(p.pos_pKs["Nterm"], 6.93)

    def test_nterm_proline(self):
        """P at N-terminus gets adjusted pKa (8.36)."""
        p = IP("PAAAA")
        self.assertAlmostEqual(p.pos_pKs["Nterm"], 8.36)

    def test_nterm_threonine(self):
        """T at N-terminus gets adjusted pKa (6.82)."""
        p = IP("TAAAA")
        self.assertAlmostEqual(p.pos_pKs["Nterm"], 6.82)

    def test_nterm_valine(self):
        """V at N-terminus gets adjusted pKa (7.44)."""
        p = IP("VAAAA")
        self.assertAlmostEqual(p.pos_pKs["Nterm"], 7.44)

    def test_nterm_glutamate(self):
        """E at N-terminus gets adjusted pKa (7.7)."""
        p = IP("EAAAA")
        self.assertAlmostEqual(p.pos_pKs["Nterm"], 7.7)

    def test_nterm_no_correction(self):
        """G at N-terminus gets no correction (default 7.5)."""
        p = IP("GAAAA")
        self.assertAlmostEqual(p.pos_pKs["Nterm"], 7.5)

    def test_cterm_aspartate(self):
        """D at C-terminus gets adjusted pKa (4.55)."""
        p = IP("AAAAD")
        self.assertAlmostEqual(p.neg_pKs["Cterm"], 4.55)

    def test_cterm_glutamate(self):
        """E at C-terminus gets adjusted pKa (4.75)."""
        p = IP("AAAAE")
        self.assertAlmostEqual(p.neg_pKs["Cterm"], 4.75)

    def test_cterm_no_correction(self):
        """K at C-terminus gets no correction (default 3.55)."""
        p = IP("AAAAK")
        self.assertAlmostEqual(p.neg_pKs["Cterm"], 3.55)

    def test_no_terminal_corrections_for_other_scales(self):
        """Non-Bjellqvist scales should NOT apply terminal corrections."""
        # M at N-terminus - only Bjellqvist should adjust
        bjellqvist = IP("MAAAA", "Bjellqvist")
        ipc2 = IP("MAAAA", "IPC2_protein")
        # Bjellqvist adjusts Nterm for M: 7.5 -> 7.0
        self.assertAlmostEqual(bjellqvist.pos_pKs["Nterm"], 7.0)
        # IPC2_protein uses its own fixed value (5.779)
        self.assertAlmostEqual(ipc2.pos_pKs["Nterm"], 5.779)


class TestAllScalesProduceDifferentResults(unittest.TestCase):
    """Verify that different scales produce different pI values."""

    def test_all_scales_give_valid_pi(self):
        """pI from every scale should be between 1 and 14."""
        for scale in pKa_scales:
            with self.subTest(scale=scale):
                p = IP(REF_PROTEIN, scale)
                pi_val = p.pi()
                self.assertGreater(pi_val, 1.0)
                self.assertLess(pi_val, 14.0)

    def test_all_scales_charge_near_zero_at_pi(self):
        """Charge at pI should be near zero for every scale."""
        for scale in pKa_scales:
            with self.subTest(scale=scale):
                p = IP(REF_PROTEIN, scale)
                pi_val = p.pi()
                self.assertAlmostEqual(p.charge_at_pH(pi_val), 0.0, places=2)

    def test_different_scales_give_different_pi(self):
        """At least some scales should produce different pI values."""
        pis = set()
        for scale in pKa_scales:
            p = IP(REF_PROTEIN, scale)
            pis.add(round(p.pi(), 2))
        # With 18 scales, we should have at least 5 distinct values
        self.assertGreater(len(pis), 5)


class TestSpecificScaleValues(unittest.TestCase):
    """Test specific known pI values for selected scales and sequences."""

    def test_ipc2_protein(self):
        """IPC2_protein pI for the long acidic peptide."""
        p = IP(LONG_ACIDIC, "IPC2_protein")
        self.assertAlmostEqual(p.pi(), 4.83, places=2)

    def test_ipc2_peptide(self):
        """IPC2_peptide pI for EFYTVDGQK."""
        p = IP(ACIDIC_PEPTIDE, "IPC2_peptide")
        self.assertAlmostEqual(p.pi(), 4.28, places=2)

    def test_ipc2_protein_charge(self):
        """IPC2_protein charge at pH 7.4 for the long acidic peptide."""
        p = IP(LONG_ACIDIC, "IPC2_protein")
        self.assertAlmostEqual(p.charge_at_pH(7.4), -4.22, places=2)

    def test_emboss_differs_from_bjellqvist(self):
        """EMBOSS and Bjellqvist should give different results."""
        bj = IP(REF_PROTEIN, "Bjellqvist")
        em = IP(REF_PROTEIN, "EMBOSS")
        self.assertNotAlmostEqual(bj.pi(), em.pi(), places=1)

    def test_ref_protein_ipc2(self):
        """Reference protein pI with IPC2_protein."""
        p = IP(REF_PROTEIN, "IPC2_protein")
        self.assertAlmostEqual(p.pi(), 6.96, places=2)


class TestInvalidInput(unittest.TestCase):
    """Test error handling for bad inputs."""

    def test_invalid_scale_raises(self):
        """An invalid pKa scale name should raise ValueError."""
        with self.assertRaises(ValueError) as cm:
            IP("INGAR", "NonexistentScale")
        self.assertIn("NonexistentScale", str(cm.exception))
        self.assertIn("Available scales", str(cm.exception))

    def test_empty_string_scale_raises(self):
        """An empty string pKa scale should raise ValueError."""
        with self.assertRaises(ValueError):
            IP("INGAR", "")

    def test_none_scale_defaults_to_bjellqvist(self):
        """pKa_scale=None should default to Bjellqvist."""
        p = IP("INGAR", None)
        self.assertEqual(p.pKa_scale, "Bjellqvist")
        self.assertAlmostEqual(p.pi(), 9.75, places=2)


class TestEdgeCases(unittest.TestCase):
    """Test edge cases and unusual sequences."""

    def test_single_amino_acid(self):
        """A single amino acid should produce a valid pI."""
        p = IP("K")
        pi = p.pi()
        self.assertGreater(pi, 1.0)
        self.assertLess(pi, 14.0)

    def test_all_aspartate(self):
        """All-Asp protein should have very low pI."""
        p = IP("DDDDDDDDD")
        self.assertLess(p.pi(), 4.0)

    def test_all_arginine(self):
        """All-Arg protein should have very high pI."""
        p = IP("RRRRRRRRR")
        self.assertGreater(p.pi(), 11.0)

    def test_all_alanine(self):
        """All-Ala protein (no charged side chains) pI from termini only."""
        p = IP("AAAAAAAAAA")
        pi = p.pi()
        # pI should be determined solely by N- and C-terminal groups
        self.assertGreater(pi, 3.0)
        self.assertLess(pi, 10.0)

    def test_histidine_rich(self):
        """A His-rich sequence should have pI near 6-8."""
        p = IP("HHHHHKKR")
        pi = p.pi()
        self.assertGreater(pi, 6.0)
        self.assertLess(pi, 12.0)

    def test_cysteine_rich(self):
        """A Cys-rich sequence."""
        p = IP("CCCCCCCK")
        pi = p.pi()
        self.assertGreater(pi, 5.0)
        self.assertLess(pi, 10.0)

    def test_seq_object_input(self):
        """Bio.Seq input should work like string."""
        pi_str = IP("INGAR").pi()
        pi_seq = IP(Seq("INGAR")).pi()
        self.assertAlmostEqual(pi_str, pi_seq)

    def test_lowercase_input(self):
        """Lowercase sequences should be handled (uppercased)."""
        pi_upper = IP("INGAR").pi()
        pi_lower = IP("ingar").pi()
        self.assertAlmostEqual(pi_upper, pi_lower)

    def test_patrickios_zero_pkas(self):
        """Patrickios scale has pKa=0 for C, Y, H — should still converge."""
        p = IP(REF_PROTEIN, "Patrickios")
        pi = p.pi()
        self.assertGreater(pi, 1.0)
        self.assertLess(pi, 14.0)
        self.assertAlmostEqual(p.charge_at_pH(pi), 0.0, places=2)


class TestProtParamIntegration(unittest.TestCase):
    """Test that ProtParam properly forwards pKa_scale parameter."""

    def test_default_pi(self):
        """Default ProtParam.isoelectric_point() matches historic value."""
        analysis = PA(REF_PROTEIN)
        self.assertAlmostEqual(analysis.isoelectric_point(), 7.72, places=2)

    def test_default_charge_at_ph(self):
        """Default charge_at_pH at the isoelectric point is near zero."""
        analysis = PA(REF_PROTEIN)
        self.assertAlmostEqual(analysis.charge_at_pH(7.72), 0.00, places=2)

    def test_ipc2_protein_via_protparam(self):
        """IPC2_protein pI via ProtParam should match direct calculation."""
        analysis = PA(REF_PROTEIN)
        pi_protparam = analysis.isoelectric_point("IPC2_protein")
        pi_direct = IP(REF_PROTEIN, "IPC2_protein").pi()
        self.assertAlmostEqual(pi_protparam, pi_direct)

    def test_charge_at_ph_with_scale(self):
        """charge_at_pH with a specified scale works through ProtParam."""
        analysis = PA(REF_PROTEIN)
        charge_protparam = analysis.charge_at_pH(7.4, "IPC2_protein")
        charge_direct = IP(REF_PROTEIN, "IPC2_protein").charge_at_pH(7.4)
        self.assertAlmostEqual(charge_protparam, charge_direct)

    def test_invalid_scale_via_protparam(self):
        """Invalid scale via ProtParam should raise ValueError."""
        analysis = PA(REF_PROTEIN)
        with self.assertRaises(ValueError):
            analysis.isoelectric_point("BadScale")

    def test_protparam_seq_object(self):
        """ProtParam with Seq object and scale parameter."""
        analysis = PA(Seq(REF_PROTEIN))
        pi = analysis.isoelectric_point("Toseland")
        self.assertGreater(pi, 1.0)
        self.assertLess(pi, 14.0)

    def test_protparam_seqrecord(self):
        """ProtParam with SeqRecord and scale parameter."""
        analysis = PA(SeqRecord(Seq(REF_PROTEIN)))
        pi = analysis.isoelectric_point()
        self.assertAlmostEqual(pi, 7.72, places=2)

    def test_all_scales_via_protparam(self):
        """All 18 scales work through ProtParam."""
        analysis = PA(REF_PROTEIN)
        for scale in pKa_scales:
            with self.subTest(scale=scale):
                pi = analysis.isoelectric_point(scale)
                self.assertGreater(pi, 1.0)
                self.assertLess(pi, 14.0)

    def test_nozaki_via_protparam(self):
        """Specific known value: Nozaki pI via ProtParam."""
        analysis = PA(REF_PROTEIN)
        pi = analysis.isoelectric_point("Nozaki")
        self.assertAlmostEqual(pi, 8.00, places=2)


class TestNumericalStability(unittest.TestCase):
    """Test that the bisection converges reliably."""

    def test_convergence_precision(self):
        """pI should be reproducible to 4 decimal places."""
        p1 = IP(REF_PROTEIN)
        p2 = IP(REF_PROTEIN)
        self.assertAlmostEqual(p1.pi(), p2.pi(), places=4)

    def test_convergence_all_scales(self):
        """Bisection should converge for every scale, even extreme pKa values."""
        for scale in pKa_scales:
            with self.subTest(scale=scale):
                p = IP("DDDDDRRRRR", scale)
                pi = p.pi()
                charge = p.charge_at_pH(pi)
                self.assertAlmostEqual(charge, 0.0, places=2)

    def test_very_acidic_protein(self):
        """A very acidic protein should have pI < 3 with wide search range."""
        # All Asp + Glu, no basic residues except termini
        p = IP("DDDDDEEEEEDDDDDEEEEE")
        pi = p.pi()
        self.assertLess(pi, 4.0)
        self.assertGreater(pi, 1.0)

    def test_very_basic_protein(self):
        """A very basic protein should have pI > 11 with wide search range."""
        p = IP("RRRRRKKKKKRRRRR")
        pi = p.pi()
        self.assertGreater(pi, 10.0)
        self.assertLess(pi, 14.0)


class TestPeptideVsProtein(unittest.TestCase):
    """IPC2_peptide vs IPC2_protein should give different results."""

    def test_peptide_vs_protein_differ(self):
        """IPC2_peptide and IPC2_protein give different pI for the same seq."""
        p_prot = IP(ACIDIC_PEPTIDE, "IPC2_protein")
        p_pep = IP(ACIDIC_PEPTIDE, "IPC2_peptide")
        # These should differ noticeably
        self.assertNotAlmostEqual(p_prot.pi(), p_pep.pi(), places=0)

    def test_ipc_peptide_vs_ipc_protein(self):
        """IPC_peptide and IPC_protein also differ."""
        p_prot = IP(ACIDIC_PEPTIDE, "IPC_protein")
        p_pep = IP(ACIDIC_PEPTIDE, "IPC_peptide")
        self.assertNotAlmostEqual(p_prot.pi(), p_pep.pi(), places=1)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
