# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#

"""Testing code for Restriction enzyme classes of Biopython.
"""

import unittest

from Bio.Restriction import *
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA


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
        self.assertTrue(EcoRV in batch)
        self.assertTrue(EcoRI in batch)
        self.assertTrue(KpnI in batch)
        self.assertTrue(SmaI not in batch)
        # Syntax sugar for the above
        self.assertTrue('EcoRV' in batch)
        self.assertFalse('SmaI' in batch)

        batch.get(EcoRV)
        self.assertRaises(ValueError, batch.get, SmaI)

        batch.remove(EcoRV)
        self.assertEqual(len(batch), 2)

        self.assertTrue(EcoRV not in batch)
        self.assertTrue('EcoRV' not in batch)

    def test_batch_analysis(self):
        """Sequence analysis with a restriction batch.
        """
        seq = Seq("AAAA" + EcoRV.site + "AAAA" + EcoRI.site + "AAAA",
                IUPACAmbiguousDNA())
        batch = RestrictionBatch([EcoRV, EcoRI])

        hits = batch.search(seq)
        self.assertEqual(hits[EcoRV], [8])
        self.assertEqual(hits[EcoRI], [16])


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
