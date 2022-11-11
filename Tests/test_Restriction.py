# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#

"""Testing code for Restriction enzyme classes of Biopython."""
import unittest

from Bio.Restriction import Analysis, Restriction, RestrictionBatch
from Bio.Restriction import CommOnly, NonComm, AllEnzymes
from Bio.Restriction import (
    Acc65I,
    Asp718I,
    BamHI,
    EcoRI,
    EcoRV,
    KpnI,
    SmaI,
    MluCI,
    McrI,
    NdeI,
    BsmBI,
    AanI,
    EarI,
    SnaI,
    SphI,
)
from Bio.Restriction import FormattedSeq
from Bio.Seq import Seq, MutableSeq
from Bio import BiopythonWarning


class SequenceTesting(unittest.TestCase):
    """Tests for dealing with input."""

    def test_sequence_object(self):
        """Test if sequence must be a Seq or MutableSeq object."""
        with self.assertRaises(TypeError):
            seq = FormattedSeq("GATC")
        seq = FormattedSeq(Seq("TAGC"))
        seq = FormattedSeq(MutableSeq("AGTC"))
        seq = FormattedSeq(seq)
        with self.assertRaises(TypeError):
            EcoRI.search("GATC")
        EcoRI.search(Seq("ATGC"))
        EcoRI.search(MutableSeq("TCAG"))

    def test_non_iupac_letters(self):
        """Test if non-IUPAC letters raise a TypeError."""
        with self.assertRaises(TypeError):
            seq = FormattedSeq(Seq("GATCZ"))

    def test_formatted_seq(self):
        """Test several methods of FormattedSeq."""
        self.assertEqual(
            str(FormattedSeq(Seq("GATC"))), "FormattedSeq(Seq('GATC'), linear=True)"
        )
        self.assertNotEqual(FormattedSeq(Seq("GATC")), FormattedSeq(Seq("TAGC")))
        self.assertNotEqual(FormattedSeq(Seq("TAGC")), Seq("TAGC"))
        self.assertEqual(FormattedSeq(Seq("ATGC")), FormattedSeq(Seq("ATGC")))
        linear_seq = FormattedSeq(Seq("T"))
        self.assertTrue(linear_seq.is_linear())
        linear_seq.circularise()
        self.assertFalse(linear_seq.is_linear())
        linear_seq.linearise()
        circular_seq = linear_seq.to_circular()
        self.assertFalse(circular_seq.is_linear())
        linear_seq = circular_seq.to_linear()
        self.assertTrue(linear_seq.is_linear())


class SimpleEnzyme(unittest.TestCase):
    """Tests for dealing with basic enzymes using the Restriction package."""

    def test_init(self):
        """Check for error during __init__."""
        with self.assertRaises(ValueError) as ve:
            Restriction.OneCut("bla-me", (Restriction.RestrictionType,), {})
            self.assertIn("hyphen", str(ve.exception))

    def setUp(self):
        """Set up some sequences for later use."""
        base_seq = Seq("AAAA")
        self.ecosite_seq = base_seq + Seq(EcoRI.site) + base_seq
        self.smasite_seq = base_seq + Seq(SmaI.site) + base_seq
        self.kpnsite_seq = base_seq + Seq(KpnI.site) + base_seq

    def test_eco_cutting(self):
        """Test basic cutting with EcoRI (5'overhang)."""
        self.assertEqual(EcoRI.site, "GAATTC")
        self.assertTrue(EcoRI.cut_once())
        self.assertFalse(EcoRI.is_blunt())
        self.assertTrue(EcoRI.is_5overhang())
        self.assertFalse(EcoRI.is_3overhang())
        self.assertEqual(EcoRI.overhang(), "5' overhang")
        self.assertTrue(EcoRI.is_defined())
        self.assertFalse(EcoRI.is_ambiguous())
        self.assertFalse(EcoRI.is_unknown())
        self.assertTrue(EcoRI.is_palindromic())
        self.assertTrue(EcoRI.is_comm())
        self.assertIn("Thermo Fisher Scientific", EcoRI.supplier_list())
        self.assertEqual(EcoRI.elucidate(), "G^AATT_C")
        self.assertEqual(EcoRI.search(self.ecosite_seq), [6])
        self.assertEqual(EcoRI.characteristic(), (1, -1, None, None, "GAATTC"))

        parts = EcoRI.catalyse(self.ecosite_seq)
        self.assertEqual(len(parts), 2)
        self.assertEqual(str(parts[1]), "AATTCAAAA")
        parts = EcoRI.catalyze(self.ecosite_seq)
        self.assertEqual(len(parts), 2)

    def test_kpn_cutting(self):
        """Test basic cutting with KpnI (3'overhang)."""
        self.assertTrue(KpnI.is_3overhang())
        self.assertFalse(KpnI.is_5overhang())
        self.assertFalse(KpnI.is_blunt())
        self.assertEqual(KpnI.overhang(), "3' overhang")
        parts = KpnI.catalyse(self.kpnsite_seq)
        self.assertEqual(len(parts), 2)
        self.assertEqual(
            KpnI.catalyse(self.kpnsite_seq), KpnI.catalyze(self.kpnsite_seq)
        )

    def test_sma_cutting(self):
        """Test basic cutting with SmaI (blunt cutter)."""
        self.assertTrue(SmaI.is_blunt())
        self.assertFalse(SmaI.is_3overhang())
        self.assertFalse(SmaI.is_5overhang())
        self.assertEqual(SmaI.overhang(), "blunt")
        parts = SmaI.catalyse(self.smasite_seq)
        self.assertEqual(len(parts), 2)
        self.assertEqual(str(parts[1]), "GGGAAAA")
        parts = SmaI.catalyze(self.smasite_seq)
        self.assertEqual(len(parts), 2)

    def test_ear_cutting(self):
        """Test basic cutting with EarI (ambiguous overhang)."""
        self.assertFalse(EarI.is_palindromic())
        self.assertFalse(EarI.is_defined())
        self.assertTrue(EarI.is_ambiguous())
        self.assertFalse(EarI.is_unknown())
        self.assertEqual(EarI.elucidate(), "CTCTTCN^NNN_N")

    def test_sna_cutting(self):
        """Test basic cutting with SnaI (unknown)."""
        self.assertEqual(SnaI.elucidate(), "? GTATAC ?")
        self.assertFalse(SnaI.is_defined())
        self.assertFalse(SnaI.is_ambiguous())
        self.assertTrue(SnaI.is_unknown())
        self.assertFalse(SnaI.is_comm())
        self.assertIsNone(SnaI.suppliers())
        self.assertEqual(SnaI.supplier_list(), [])
        with self.assertRaises(TypeError):
            SnaI.buffers("no company")

    def test_circular_sequences(self):
        """Deal with cutting circular sequences."""
        parts = EcoRI.catalyse(self.ecosite_seq, linear=False)
        self.assertEqual(len(parts), 1)
        locations = EcoRI.search(parts[0], linear=False)
        self.assertEqual(locations, [1])

        parts = KpnI.catalyse(self.kpnsite_seq, linear=False)
        self.assertEqual(len(parts), 1)
        locations = KpnI.search(parts[0], linear=False)
        self.assertEqual(locations, [1])

        parts = SmaI.catalyse(self.smasite_seq, linear=False)
        self.assertEqual(len(parts), 1)
        locations = SmaI.search(parts[0], linear=False)
        self.assertEqual(locations, [1])

        self.assertEqual(
            EarI.search(FormattedSeq(Seq("CTCTTCAAAAA")), linear=False), [8]
        )
        self.assertEqual(
            SnaI.search(FormattedSeq(Seq("GTATACAAAAA")), linear=False), [1]
        )

    def test_shortcuts(self):
        """Check if '/' and '//' work as '.search' and '.catalyse'."""
        self.assertEqual(EcoRI / self.ecosite_seq, [6])
        self.assertEqual(self.ecosite_seq / EcoRI, [6])
        self.assertEqual(len(EcoRI // self.ecosite_seq), 2)
        self.assertEqual(len(self.ecosite_seq // EcoRI), 2)

    def test_cutting_border_positions(self):
        """Check if cutting after first and penultimate position works."""
        # Use EarI, cuts as follows: CTCTTCN^NNN_N
        seq = Seq("CTCTTCA")
        self.assertEqual(EarI.search(seq), [])
        seq += "A"
        self.assertEqual(EarI.search(seq), [8])
        # Recognition site on reverse-complement strand
        seq = Seq("AAAAGAAGAG")
        self.assertEqual(EarI.search(seq), [])
        seq = "A" + seq
        self.assertEqual(EarI.search(seq), [2])

    def test_recognition_site_on_both_strands(self):
        """Check if recognition sites on both strands are properly handled."""
        seq = Seq("CTCTTCGAAGAG")
        self.assertEqual(EarI.search(seq), [3, 8])

    def test_overlapping_cut_sites(self):
        """Check if overlapping recognition sites are properly handled."""
        seq = Seq("CATGCACGCATGCATGCACGC")
        self.assertEqual(SphI.search(seq), [13, 17])


class EnzymeComparison(unittest.TestCase):
    """Tests for comparing various enzymes."""

    def test_basic_isochizomers(self):
        """Test to be sure isochizomer and neoschizomers are as expected."""
        self.assertEqual(Acc65I.isoschizomers(), [Asp718I, KpnI])
        self.assertEqual(Acc65I.elucidate(), "G^GTAC_C")
        self.assertEqual(Asp718I.elucidate(), "G^GTAC_C")
        self.assertEqual(KpnI.elucidate(), "G_GTAC^C")
        self.assertTrue(Acc65I.is_isoschizomer(KpnI))
        self.assertFalse(Acc65I.is_equischizomer(KpnI))
        self.assertTrue(Acc65I.is_neoschizomer(KpnI))
        self.assertIn(Acc65I, Asp718I.equischizomers())
        self.assertIn(KpnI, Asp718I.neoschizomers())
        self.assertIn(KpnI, Acc65I.isoschizomers())

    def test_comparisons(self):
        """Test comparison operators between different enzymes."""
        # Comparison of iso- and neoschizomers
        self.assertEqual(Acc65I, Acc65I)
        self.assertNotEqual(Acc65I, KpnI)
        self.assertFalse(Acc65I == Asp718I)  # noqa: A500
        # self.assertNotEqual(Acc65I, Asp718I) it doesn't work as expected
        self.assertFalse(Acc65I != Asp718I)  # noqa: A500
        self.assertNotEqual(Acc65I, EcoRI)
        self.assertTrue(Acc65I >> KpnI)
        self.assertFalse(Acc65I >> Asp718I)

        # Compare length of recognition sites
        self.assertFalse(EcoRI >= EcoRV)
        self.assertGreaterEqual(EcoRV, EcoRI)
        with self.assertRaises(NotImplementedError):
            EcoRV >= 3
        self.assertFalse(EcoRI > EcoRV)
        self.assertGreater(EcoRV, EcoRI)
        with self.assertRaises(NotImplementedError):
            EcoRV > 3
        self.assertLessEqual(EcoRI, EcoRV)
        self.assertFalse(EcoRV <= EcoRI)
        with self.assertRaises(NotImplementedError):
            EcoRV <= 3
        self.assertLess(EcoRI, EcoRV)
        self.assertFalse(EcoRV < EcoRI)
        with self.assertRaises(NotImplementedError):
            EcoRV < 3

        # Compare compatible overhangs
        self.assertTrue(Acc65I % Asp718I)
        self.assertTrue(Acc65I % Acc65I)
        self.assertFalse(Acc65I % KpnI)
        with self.assertRaises(TypeError):
            Acc65I % "KpnI"
        self.assertTrue(SmaI % EcoRV)
        self.assertTrue(EarI % EarI)
        self.assertIn(EcoRV, SmaI.compatible_end())
        self.assertIn(Acc65I, Asp718I.compatible_end())


class RestrictionBatchPrintTest(unittest.TestCase):
    """Tests Restriction.Analysis printing functionality."""

    def createAnalysis(self, seq_str, batch_ary):
        """Restriction.Analysis creation helper method."""
        rb = Restriction.RestrictionBatch(batch_ary)
        seq = Seq(seq_str)
        return Restriction.Analysis(rb, seq)

    def assertAnalysisFormat(self, analysis, expected):
        """Test make_format.

        Test that the Restriction.Analysis make_format(print_that) matches
        some string.
        """
        dct = analysis.mapping
        ls, nc = [], []
        for k, v in dct.items():
            if v:
                ls.append((k, v))
            else:
                nc.append(k)
        result = analysis.make_format(ls, "", [], "")
        self.assertEqual(result.replace(" ", ""), expected.replace(" ", ""))

    def test_make_format_map1(self):
        """Test that print_as('map'); print_that() correctly wraps round.

        1. With no marker.
        """
        analysis = self.createAnalysis(
            "CCAGTCTATAATTCG"
            + Restriction.BamHI.site
            + "GCGGCATCATACTCGAATATCGCGTGATGATACGTAGTAATTACGCATG",
            ["BamHI"],
        )
        analysis.print_as("map")
        expected = [
            "                17 BamHI",
            "                |                                           ",
            "CCAGTCTATAATTCGGGATCCGCGGCATCATACTCGAATATCGCGTGATGATACGTAGTA",
            "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||",
            "GGTCAGATATTAAGCCCTAGGCGCCGTAGTATGAGCTTATAGCGCACTACTATGCATCAT",
            "1                                                         60",
            "",
            "ATTACGCATG",
            "||||||||||",
            "TAATGCGTAC",
            "61                          70",
            "",
            "",
        ]
        self.assertAnalysisFormat(analysis, "\n".join(expected))

    def test_make_format_map2(self):
        """Test that print_as('map'); print_that() correctly wraps round.

        2. With marker.
        """
        analysis = self.createAnalysis(
            "CCAGTCTATAATTCG"
            + Restriction.BamHI.site
            + "GCGGCATCATACTCGA"
            + Restriction.BamHI.site
            + "ATATCGCGTGATGATA"
            + Restriction.NdeI.site
            + "CGTAGTAATTACGCATG",
            ["NdeI", "EcoRI", "BamHI", "BsmBI"],
        )
        analysis.print_as("map")
        expected = [
            "                17 BamHI",
            "                |                                           ",
            "                |                     39 BamHI",
            "                |                     |                     ",
            "CCAGTCTATAATTCGGGATCCGCGGCATCATACTCGAGGATCCATATCGCGTGATGATAC",
            "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||",
            "GGTCAGATATTAAGCCCTAGGCGCCGTAGTATGAGCTCCTAGGTATAGCGCACTACTATG",
            "1                                                         60",
            "",
            " 62 NdeI",
            " |                                                          ",
            "ATATGCGTAGTAATTACGCATG",
            "||||||||||||||||||||||",
            "TATACGCATCATTAATGCGTAC",
            "61                          82",
            "",
            "",
        ]
        self.assertAnalysisFormat(analysis, "\n".join(expected))

    def test_make_format_map3(self):
        """Test that print_as('map'); print_that() correctly wraps round.

        3. With marker restricted.
        """
        analysis = self.createAnalysis(
            "CCAGTCTATAATTCG"
            + Restriction.BamHI.site
            + "GCGGCATCATACTCGA"
            + Restriction.BamHI.site
            + "ATATCGCGTGATGATA"
            + Restriction.EcoRV.site
            + "CGTAGTAATTACGCATG",
            ["NdeI", "EcoRI", "BamHI", "BsmBI"],
        )
        analysis.print_as("map")
        expected = [
            "                17 BamHI",
            "                |                                           ",
            "                |                     39 BamHI",
            "                |                     |                     ",
            "CCAGTCTATAATTCGGGATCCGCGGCATCATACTCGAGGATCCATATCGCGTGATGATAG",
            "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||",
            "GGTCAGATATTAAGCCCTAGGCGCCGTAGTATGAGCTCCTAGGTATAGCGCACTACTATC",
            "1                                                         60",
            "",
            "ATATCCGTAGTAATTACGCATG",
            "||||||||||||||||||||||",
            "TATAGGCATCATTAATGCGTAC",
            "61                          82",
            "",
            "",
        ]
        self.assertAnalysisFormat(analysis, "\n".join(expected))

    def test_change(self):
        """Test that change() changes something."""
        seq = Seq(
            "CCAGTCTATAATTCG"
            + BamHI.site
            + "GCGGCATCATACTCGA"
            + BamHI.site
            + "ATATCGCGTGATGATA"
            + EcoRV.site
            + "CGTAGTAATTACGCATG"
        )
        batch = NdeI + EcoRI + BamHI + BsmBI
        analysis = Analysis(batch, seq)
        self.assertEqual(analysis.full()[BamHI], [17, 39])
        batch = NdeI + EcoRI + BsmBI
        seq += NdeI.site
        analysis.change(sequence=seq)
        analysis.change(rb=batch)
        self.assertEqual(len(analysis.full()), 3)
        self.assertEqual(analysis.full()[NdeI], [85])
        with self.assertRaises(AttributeError):
            analysis.change(**{"NameWidth": 3, "KonsoleWidth": 40})  # Console


class RestrictionBatches(unittest.TestCase):
    """Tests for dealing with batches of restriction enzymes."""

    def test_creating_batch(self):
        """Creating and modifying a restriction batch."""
        batch = RestrictionBatch()
        self.assertEqual(batch.suppl_codes()["N"], "New England Biolabs")
        self.assertTrue(batch.is_restriction(EcoRI))
        batch = RestrictionBatch([EcoRI])
        batch.add(KpnI)
        batch += EcoRV
        self.assertEqual(len(batch), 3)
        self.assertEqual(batch.elements(), ["EcoRI", "EcoRV", "KpnI"])
        # Problem with Python 3, as sequence of list may be different:
        # self.assertEqual(batch.as_string(), ['EcoRI', 'KpnI', 'EcoRV'])
        self.assertIn("EcoRI", batch.as_string())

        # The usual way to test batch membership
        self.assertIn(EcoRV, batch)
        self.assertIn(EcoRI, batch)
        self.assertIn(KpnI, batch)
        self.assertNotIn(SmaI, batch)
        # Syntax sugar for the above
        self.assertIn("EcoRV", batch)
        self.assertNotIn("SmaI", batch)

        batch.get(EcoRV)
        self.assertRaises(ValueError, batch.get, SmaI)
        batch.get(SmaI, add=True)
        self.assertEqual(len(batch), 4)
        batch.remove(SmaI)
        batch.remove(EcoRV)
        self.assertEqual(len(batch), 2)

        self.assertNotIn(EcoRV, batch)
        self.assertNotIn("EcoRV", batch)

        # Creating a batch by addition of restriction enzymes
        new_batch = EcoRI + KpnI
        self.assertEqual(batch, new_batch)
        # or by addition of a batch with an enzyme
        another_new_batch = new_batch + EcoRV
        new_batch += EcoRV
        self.assertEqual(another_new_batch, new_batch)
        self.assertRaises(TypeError, EcoRI.__add__, 1)

        # Create a batch with suppliers and other supplier related methods
        # These tests may be 'update sensitive' since company names and
        # products may change often...
        batch = RestrictionBatch((), ("S"))  # Sigma
        self.assertEqual(batch.current_suppliers(), ["Sigma Chemical Corporation"])
        self.assertIn(EcoRI, batch)
        self.assertNotIn(AanI, batch)
        batch.add_supplier("B")  # Thermo Fisher Scientific
        self.assertIn(AanI, batch)

    def test_batch_analysis(self):
        """Sequence analysis with a restriction batch."""
        seq = Seq("AAAA" + EcoRV.site + "AAAA" + EcoRI.site + "AAAA")
        batch = RestrictionBatch([EcoRV, EcoRI])

        hits = batch.search(seq)
        self.assertEqual(hits[EcoRV], [8])
        self.assertEqual(hits[EcoRI], [16])

    def test_premade_batches(self):
        """Test content of premade batches CommOnly, NoComm, AllEnzymes."""
        self.assertEqual(len(AllEnzymes), (len(CommOnly) + len(NonComm)))
        self.assertTrue(len(AllEnzymes) > len(CommOnly) > len(NonComm))

    def test_search_premade_batches(self):
        """Test search with pre-made batches CommOnly, NoComm, AllEnzymes."""
        seq = Seq("ACCCGAATTCAAAACTGACTGATCGATCGTCGACTG")
        search = AllEnzymes.search(seq)
        self.assertEqual(search[MluCI], [6])
        # Check if '/' operator works as 'search':
        search = CommOnly / seq
        self.assertEqual(search[MluCI], [6])
        # Also in reverse order:
        search = seq / NonComm
        self.assertEqual(search[McrI], [28])

    def test_analysis_restrictions(self):
        """Test Fancier restriction analysis."""
        new_seq = Seq("TTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAA")
        rb = RestrictionBatch([EcoRI, KpnI, EcoRV])
        ana = Analysis(rb, new_seq, linear=False)
        # Output only the result for enzymes which cut blunt:
        self.assertEqual(ana.blunt(), {EcoRV: []})
        self.assertEqual(ana.full(), {KpnI: [], EcoRV: [], EcoRI: [33]})
        # Output only the result for enzymes which have a site:
        self.assertEqual(ana.with_sites(), {EcoRI: [33]})
        # Output only the enzymes which have no site:
        self.assertEqual(ana.without_site(), {KpnI: [], EcoRV: []})
        self.assertEqual(ana.with_site_size([32]), {})
        # Output only enzymes which produce 5' overhangs
        self.assertEqual(ana.overhang5(), {EcoRI: [33]})
        # Output only enzymes which produce 3' overhangs
        self.assertEqual(ana.overhang3(), {KpnI: []})
        # Output only enzymes which produce defined ends
        self.assertEqual(ana.defined(), {KpnI: [], EcoRV: [], EcoRI: [33]})
        # Output only enzymes hich cut N times
        self.assertEqual(ana.with_N_sites(2), {})
        # The enzymes which cut between position x and y:
        with self.assertRaises(TypeError):
            ana.only_between("t", 20)
        with self.assertRaises(TypeError):
            ana.only_between(1, "t")
        self.assertEqual(ana.only_between(1, 20), {})
        self.assertEqual(ana.only_between(20, 34), {EcoRI: [33]})
        # Mix start/end order:
        self.assertEqual(ana.only_between(34, 20), {EcoRI: [33]})
        self.assertEqual(ana.only_outside(20, 34), {})
        with self.assertWarns(BiopythonWarning):
            ana.with_name(["fake"])
        self.assertEqual(ana.with_name([EcoRI]), {EcoRI: [33]})
        self.assertEqual((ana._boundaries(1, 20)[:2]), (1, 20))
        # Reverse order:
        self.assertEqual((ana._boundaries(20, 1)[:2]), (1, 20))
        # Fix negative start:
        self.assertEqual((ana._boundaries(-1, 20)[:2]), (20, 33))
        # Fix negative end:
        self.assertEqual((ana._boundaries(1, -1)[:2]), (1, 33))
        # Sites in- and outside of boundaries
        new_seq = Seq("GAATTCAAAAAAGAATTC")
        rb = RestrictionBatch([EcoRI])
        ana = Analysis(rb, new_seq)
        # Cut at least inside
        self.assertEqual(ana.between(1, 7), {EcoRI: [2, 14]})
        # Cut at least inside and report only inside site
        self.assertEqual(ana.show_only_between(1, 7), {EcoRI: [2]})
        # Cut at least outside
        self.assertEqual(ana.outside(1, 7), {EcoRI: [2, 14]})
        # Don't cut within
        self.assertEqual(ana.do_not_cut(7, 12), {EcoRI: [2, 14]})


class TestPrintOutputs(unittest.TestCase):
    """Class to test various print outputs."""

    import sys
    from io import StringIO

    def test_supplier(self):
        """Test output of supplier list for different enzyme types."""
        out = self.StringIO()
        self.sys.stdout = out
        EcoRI.suppliers()
        self.assertIn("Thermo Fisher Scientific", out.getvalue())
        self.assertIsNone(SnaI.suppliers())
        EcoRI.all_suppliers()  # Independent of enzyme, list of all suppliers
        self.assertIn("Agilent Technologies", out.getvalue())
        batch = EcoRI + SnaI
        batch.show_codes()
        self.assertIn("N = New England Biolabs", out.getvalue())
        self.sys.stdout = self.sys.__stdout__

    def test_print_that(self):
        """Test print_that function."""
        out = self.StringIO()
        self.sys.stdout = out
        my_batch = EcoRI + SmaI + KpnI
        my_seq = Seq("GAATTCCCGGGATATA")  # EcoRI and SmaI sites
        analysis = Analysis(my_batch, my_seq)
        analysis.print_that(None, title="My sequence\n\n", s1="Non Cutters\n\n")
        self.assertIn("My sequence", out.getvalue())
        self.assertIn("Non Cutters", out.getvalue())
        self.assertIn("2.", out.getvalue())
        self.sys.stdout = self.sys.__stdout__

    def test_str_method(self):
        """Test __str__ and __repr__ outputs."""
        batch = EcoRI + SmaI + KpnI
        self.assertEqual(str(batch), "EcoRI+KpnI+SmaI")
        batch += Asp718I
        batch += SnaI
        self.assertEqual(str(batch), "Asp718I+EcoRI...SmaI+SnaI")
        self.assertEqual(
            repr(batch),
            "RestrictionBatch(['Asp718I', 'EcoRI', 'KpnI', 'SmaI', 'SnaI'])",
        )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
