# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#

"""Testing code for Restriction enzyme classes of Biopython.
"""

from Bio.Restriction import Analysis, Restriction, RestrictionBatch
from Bio.Restriction import Acc65I, Asp718I, EcoRI, EcoRV, KpnI, SmaI
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio import BiopythonWarning

from sys import version_info
if version_info[0] < 3:
    try:
        import unittest2 as unittest
    except ImportError:
        from Bio import MissingPythonDependencyError
        raise MissingPythonDependencyError("Under Python 2 this test needs the unittest2 library")
else:
    import unittest


class SimpleEnzyme(unittest.TestCase):
    """Tests for dealing with basic enzymes using the Restriction package.
    """

    def setUp(self):
        base_seq = Seq("AAAA", IUPACAmbiguousDNA())
        self.ecosite_seq = base_seq + Seq(EcoRI.site,
                                          IUPACAmbiguousDNA()) + base_seq

    def test_eco_cutting(self):
        """Test basic cutting with EcoRI.
        """
        self.assertEqual(EcoRI.site, 'GAATTC')
        self.assertFalse(EcoRI.is_blunt())
        self.assertTrue(EcoRI.is_5overhang())
        self.assertFalse(EcoRI.is_3overhang())
        self.assertEqual(EcoRI.elucidate(), "G^AATT_C")
        self.assertEqual(EcoRI.search(self.ecosite_seq), [6])

        parts = EcoRI.catalyse(self.ecosite_seq)
        self.assertEqual(len(parts), 2)
        self.assertEqual(str(parts[1]), "AATTCAAAA")
        parts = EcoRI.catalyze(self.ecosite_seq)
        self.assertEqual(len(parts), 2)

    def test_circular_sequences(self):
        """Deal with cutting circular sequences.
        """
        parts = EcoRI.catalyse(self.ecosite_seq, linear=False)
        self.assertEqual(len(parts), 1)
        locations = EcoRI.search(parts[0], linear=False)
        self.assertEqual(locations, [1])


class EnzymeComparison(unittest.TestCase):
    """Tests for comparing various enzymes.
    """

    def test_basic_isochizomers(self):
        """Test to be sure isochizomer and neoschizomers are as expected.
        """
        self.assertEqual(Acc65I.isoschizomers(), [Asp718I, KpnI])
        self.assertEqual(Acc65I.elucidate(), 'G^GTAC_C')
        self.assertEqual(Asp718I.elucidate(), 'G^GTAC_C')
        self.assertEqual(KpnI.elucidate(), 'G_GTAC^C')

    def test_comparisons(self):
        """Comparison operators between iso and neoschizomers.
        """
        self.assertEqual(Acc65I, Acc65I)
        self.assertNotEqual(Acc65I, KpnI)
        self.assertFalse(Acc65I == Asp718I)
        self.assertFalse(Acc65I != Asp718I)
        self.assertNotEqual(Acc65I, EcoRI)

        self.assertTrue(Acc65I >> KpnI)
        self.assertFalse(Acc65I >> Asp718I)

        self.assertTrue(Acc65I % Asp718I)
        self.assertTrue(Acc65I % Acc65I)
        self.assertFalse(Acc65I % KpnI)


class RestrictionBatchPrintTest(unittest.TestCase):
    """Tests Restriction.Analysis printing functionality.
    """

    def createAnalysis(self, seq_str, batch_ary):
        """Restriction.Analysis creation helper method."""
        rb = Restriction.RestrictionBatch(batch_ary)
        seq = Seq(seq_str)
        return Restriction.Analysis(rb, seq)

    def assertAnalysisFormat(self, analysis, expected):
        """Asserts that the Restriction.Analysis make_format(print_that) matches some string."""
        dct = analysis.mapping
        ls, nc = [], []
        for k, v in dct.items():
            if v:
                ls.append((k, v))
            else:
                nc.append(k)
        result = analysis.make_format(ls, '', [], '')
        self.assertEqual(result.replace(' ', ''), expected.replace(' ', ''))

    def test_make_format_map1(self):
        """Make sure print_as('map'); print_that() does not error on wrap round with no markers.
        """
        analysis = self.createAnalysis(
            'CCAGTCTATAATTCG' +
            Restriction.BamHI.site +
            'GCGGCATCATACTCGAATATCGCGTGATGATACGTAGTAATTACGCATG',
            ["BamHI"])
        analysis.print_as('map')
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
            "", ""]
        self.assertAnalysisFormat(analysis, '\n'.join(expected))

    def test_make_format_map2(self):
        """Make sure print_as('map'); print_that() does not error on wrap round with marker.
        """
        analysis = self.createAnalysis(
            'CCAGTCTATAATTCG' +
            Restriction.BamHI.site +
            'GCGGCATCATACTCGA' +
            Restriction.BamHI.site +
            'ATATCGCGTGATGATA' +
            Restriction.NdeI.site +
            'CGTAGTAATTACGCATG',
            ["NdeI", "EcoRI", "BamHI", "BsmBI"])
        analysis.print_as('map')
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
            "", ""]
        self.assertAnalysisFormat(analysis, '\n'.join(expected))

    def test_make_format_map3(self):
        """Make sure print_as('map'); print_that() does not error on wrap round with marker restricted.
        """
        analysis = self.createAnalysis(
            'CCAGTCTATAATTCG' +
            Restriction.BamHI.site +
            'GCGGCATCATACTCGA' +
            Restriction.BamHI.site +
            'ATATCGCGTGATGATA' +
            Restriction.EcoRV.site +
            'CGTAGTAATTACGCATG',
            ["NdeI", "EcoRI", "BamHI", "BsmBI"])
        analysis.print_as('map')
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
            "", ""]
        self.assertAnalysisFormat(analysis, '\n'.join(expected))


class RestrictionBatches(unittest.TestCase):
    """Tests for dealing with batches of restriction enzymes.
    """

    def test_creating_batch(self):
        """Creating and modifying a restriction batch.
        """
        batch = RestrictionBatch([EcoRI])
        batch.add(KpnI)
        batch += EcoRV
        self.assertEqual(len(batch), 3)

        # The usual way to test batch membership
        self.assertIn(EcoRV, batch)
        self.assertIn(EcoRI, batch)
        self.assertIn(KpnI, batch)
        self.assertNotIn(SmaI, batch)
        # Syntax sugar for the above
        self.assertIn('EcoRV', batch)
        self.assertNotIn('SmaI', batch)

        batch.get(EcoRV)
        self.assertRaises(ValueError, batch.get, SmaI)

        batch.remove(EcoRV)
        self.assertEqual(len(batch), 2)

        self.assertNotIn(EcoRV, batch)
        self.assertNotIn('EcoRV', batch)

    def test_batch_analysis(self):
        """Sequence analysis with a restriction batch.
        """
        seq = Seq("AAAA" + EcoRV.site + "AAAA" + EcoRI.site + "AAAA",
                  IUPACAmbiguousDNA())
        batch = RestrictionBatch([EcoRV, EcoRI])

        hits = batch.search(seq)
        self.assertEqual(hits[EcoRV], [8])
        self.assertEqual(hits[EcoRI], [16])

    def test_analysis_restrictions(self):
        """Test Fancier restriction analysis
        """
        new_seq = Seq('TTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAA', IUPACAmbiguousDNA())
        rb = RestrictionBatch([EcoRI, KpnI, EcoRV])
        ana = Analysis(rb, new_seq, linear=False)
        self.assertEqual(ana.blunt(), {EcoRV: []})  # output only the result for enzymes which cut blunt
        self.assertEqual(ana.full(), {KpnI: [], EcoRV: [], EcoRI: [33]})
        self.assertEqual(ana.with_sites(), {EcoRI: [33]})  # output only the result for enzymes which have a site
        self.assertEqual(ana.without_site(), {KpnI: [], EcoRV: []})  # output only the enzymes which have no site
        self.assertEqual(ana.with_site_size([32]), {})
        self.assertEqual(ana.only_between(1, 20), {})  # the enzymes which cut between position 1 and 20
        self.assertEqual(ana.only_between(20, 34), {EcoRI: [33]})  # etc...
        self.assertEqual(ana.only_between(34, 20), {EcoRI: [33]})  # mix start end order
        self.assertEqual(ana.only_outside(20, 34), {})
        with self.assertWarns(BiopythonWarning):
            ana.with_name(['fake'])
        self.assertEqual(ana.with_name([EcoRI]), {EcoRI: [33]})
        self.assertEqual((ana._boundaries(1, 20)[:2]), (1, 20))
        self.assertEqual((ana._boundaries(20, 1)[:2]), (1, 20))  # reverse order
        self.assertEqual((ana._boundaries(-1, 20)[:2]), (20, 33))  # fix negative start


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
