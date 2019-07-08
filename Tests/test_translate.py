# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


"""Tests of the transcription and translation methods of Seq objects."""


import unittest

from Bio import Seq
from Bio import Alphabet
from Bio.Alphabet import IUPAC


class TestTranscriptionTranslation(unittest.TestCase):

    def test_transcription(self):
        s = "ATA"
        dna = Seq.Seq(s, IUPAC.unambiguous_dna)
        rna = dna.transcribe()
        self.assertEqual(str(rna), "AUA")
        s = "GAAAATTCATTTTCTTTGGACTTTCTCTGAAATCCGAGTCCTAGGAAAGATGCGTGAGATTCTTCATATT"
        dna = Seq.Seq(s, IUPAC.unambiguous_dna)
        rna = dna.transcribe()
        self.assertEqual(str(rna), 'GAAAAUUCAUUUUCUUUGGACUUUCUCUGAAAUCCGAGUCCUAGGAAAGAUGCGUGAGAUUCUUCAUAUU')
        s = "GAAAAUUCAUUUUCUUUGGACUUUCUCUGAAAUCCGAGUCCUAGGAAAGAUGCGUGAGAUUCUUCAUAUU"
        rna = Seq.Seq(s, IUPAC.unambiguous_rna)
        dna = rna.back_transcribe()
        self.assertEqual(str(dna), 'GAAAATTCATTTTCTTTGGACTTTCTCTGAAATCCGAGTCCTAGGAAAGATGCGTGAGATTCTTCATATT')

    def test_translation(self):
        s = ""
        dna = Seq.Seq(s, IUPAC.unambiguous_dna)
        protein = dna.translate(to_stop=True)
        self.assertEqual(str(protein), "")
        s = "TAA"
        dna = Seq.Seq(s, IUPAC.unambiguous_dna)
        protein = dna.translate(to_stop=True)
        self.assertEqual(str(protein), "")
        s = "GAAAATTCATTTTCTTTGGACTTTCTCTGAAATCCGAGTCCTAGGAAAGATGCGTGAGATTCTTCA"
        dna = Seq.Seq(s, IUPAC.unambiguous_dna)
        protein = dna.translate(to_stop=True)
        self.assertEqual(str(protein), 'ENSFSLDFL')
        s = "GAA"
        dna = Seq.Seq(s, IUPAC.unambiguous_dna)
        protein = dna.translate(15, to_stop=True)
        self.assertEqual(str(protein), "E")
        s = "ATA"
        dna = Seq.Seq(s, IUPAC.unambiguous_dna)
        protein = dna.translate('Vertebrate Mitochondrial', to_stop=True)
        self.assertEqual(str(protein), 'M')
        s = "GAAAATTCATTTTCTTTGGACTTTCTCTGAAATCCGAGTCCTAGGAAAGATGCGTGAGATTCTTCATAT"
        dna = Seq.Seq(s, IUPAC.unambiguous_dna)
        protein = dna.translate('SGC8', to_stop=True)
        self.assertEqual(str(protein), 'ENSFSLDFLWNPSPSNDAWDSSY')

    def test_dna_rna_translation(self):
        s = "TCAAAAAGGTGCATCTAGATG"
        dna = Seq.Seq(s, IUPAC.unambiguous_dna)
        protein = dna.translate(to_stop=True)
        self.assertIsInstance(protein.alphabet, IUPAC.IUPACProtein)
        self.assertEqual(str(protein), "SKRCI")
        gapped_protein = dna.translate()
        self.assertIsInstance(gapped_protein.alphabet, Alphabet.HasStopCodon)
        self.assertEqual(str(gapped_protein), "SKRCI*M")
        # The table used here has "AGG" as a stop codon:
        p2 = dna.translate(table=2, to_stop=True)
        self.assertEqual(str(p2), "SK")
        p2 = dna.translate(table=2)
        self.assertEqual(str(p2), "SK*CI*M")
        p2 = dna.translate(table=2, stop_symbol="+")
        self.assertEqual(str(p2), "SK+CI+M")
        r = s.replace("T", "U")
        rna = Seq.Seq(r, IUPAC.unambiguous_rna)
        protein = rna.translate(to_stop=True)
        self.assertIsInstance(protein.alphabet, IUPAC.IUPACProtein)
        self.assertEqual(str(protein), "SKRCI")
        gapped_protein = rna.translate()
        self.assertEqual(str(gapped_protein), "SKRCI*M")

    def test_ambiguous(self):
        s = "RATGATTARAATYTA"
        dna = Seq.Seq(s, IUPAC.ambiguous_dna)
        protein = dna.translate('Vertebrate Mitochondrial')
        self.assertEqual(str(protein), "BD*NL")
        stop_protein = dna.translate('SGC1', to_stop=True)
        self.assertEqual(str(stop_protein), "BD")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
