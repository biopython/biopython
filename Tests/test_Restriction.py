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
        self.assertEqual(parts[1].tostring(), "AATTCAAAA")
        parts = EcoRI.catalyze(self.ecosite_seq)
        self.assertEqual(len(parts), 2)

    def test_circular_sequences(self):
        """Deal with cutting circular sequences.
        """
        parts = EcoRI.catalyse(self.ecosite_seq, linear = False)
        self.assertEqual(len(parts), 1)
        locations = EcoRI.search(parts[0], linear = False)
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
        self.assert_(Acc65I == Acc65I)
        self.assert_(Acc65I != KpnI)
        self.assert_(not (Acc65I == Asp718I))
        self.assert_(not (Acc65I != Asp718I))
        self.assert_(Acc65I != EcoRI)

        self.assert_(Acc65I >> KpnI)
        self.assert_(not (Acc65I >> Asp718I))

        self.assert_(Acc65I % Asp718I)
        self.assert_(Acc65I % Acc65I)
        self.assert_(not (Acc65I % KpnI))


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
        self.assert_(EcoRV in batch)
        self.assert_(EcoRI in batch)
        self.assert_(KpnI in batch)
        self.assert_(SmaI not in batch)
        # Syntax sugar for the above
        self.assert_('EcoRV' in batch)
        self.assert_('SmaI' not in batch)

        batch.get(EcoRV)
        self.assertRaises(ValueError, batch.get, SmaI)

        batch.remove(EcoRV)
        self.assertEqual(len(batch), 2)

        self.assert_(EcoRV not in batch)
        self.assert_('EcoRV' not in batch)

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
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
