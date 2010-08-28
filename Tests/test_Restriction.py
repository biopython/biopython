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

    def tearDown(self):
        pass

    def test_eco_cutting(self):
        """Test basic cutting with EcoRI.
        """
        assert EcoRI.site == 'GAATTC'
        assert EcoRI.is_blunt() == False
        assert EcoRI.is_5overhang() == True
        assert EcoRI.is_3overhang() == False
        assert EcoRI.elucidate() == "G^AATT_C"
        assert EcoRI.search(self.ecosite_seq) == [6]

        parts = EcoRI.catalyse(self.ecosite_seq)
        assert len(parts) == 2
        assert parts[1].tostring() == "AATTCAAAA"
        parts = EcoRI.catalyze(self.ecosite_seq)
        assert len(parts) == 2

    def test_circular_sequences(self):
        """Deal with cutting circular sequences.
        """
        parts = EcoRI.catalyse(self.ecosite_seq, linear = False)
        assert len(parts) == 1
        locations = EcoRI.search(parts[0], linear = False)
        assert locations == [1]

class EnzymeComparison(unittest.TestCase):
    """Tests for comparing various enzymes.
    """
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_basic_isochizomers(self):
        """Test to be sure isochizomer and neoschizomers are as expected.
        """
        assert Acc65I.isoschizomers() == [Asp718I, KpnI]
        assert Acc65I.elucidate() == 'G^GTAC_C'
        assert Asp718I.elucidate() == 'G^GTAC_C'
        assert KpnI.elucidate() == 'G_GTAC^C'

    def test_comparisons(self):
        """Comparison operators between iso and neoschizomers.
        """
        assert Acc65I == Acc65I
        assert Acc65I != KpnI
        assert not(Acc65I == Asp718I)
        assert not(Acc65I != Asp718I)
        assert Acc65I != EcoRI

        assert Acc65I >> KpnI
        assert not(Acc65I >> Asp718I)

        assert Acc65I % Asp718I
        assert Acc65I % Acc65I
        assert not(Acc65I % KpnI)

class RestrictionBatches(unittest.TestCase):
    """Tests for dealing with batches of restriction enzymes.
    """
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_creating_batch(self):
        """Creating and modifying a restriction batch.
        """
        batch = RestrictionBatch([EcoRI])
        batch.add(KpnI)
        batch += EcoRV
        assert len(batch) == 3

        batch.get(EcoRV)
        try:
            batch.get(SmaI)
            raise AssertionError("No error with non-existent enzyme.")
        except ValueError:
            pass

        batch.remove(EcoRV)
        assert len(batch) == 2

    def test_batch_analysis(self):
        """Sequence analysis with a restriction batch.
        """
        seq = Seq("AAAA" + EcoRV.site + "AAAA" + EcoRI.site + "AAAA",
                IUPACAmbiguousDNA())
        batch = RestrictionBatch([EcoRV, EcoRI])

        hits = batch.search(seq)
        assert hits[EcoRV] == [8] and hits[EcoRI] == [16]


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
