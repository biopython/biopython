# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


"""Tests of the transcription and translation methods of Seq objects."""


import unittest

from Bio import Seq


class TestTranscriptionTranslation(unittest.TestCase):
    def test_transcription(self):
        s = "ATA"
        dna = Seq.Seq(s)
        rna = dna.transcribe()
        self.assertEqual(rna, "AUA")
        s = "GAAAATTCATTTTCTTTGGACTTTCTCTGAAATCCGAGTCCTAGGAAAGATGCGTGAGATTCTTCATATT"
        dna = Seq.Seq(s)
        rna = dna.transcribe()
        self.assertEqual(
            rna,
            "GAAAAUUCAUUUUCUUUGGACUUUCUCUGAAAUCCGAGUCCUAGGAAAGAUGCGUGAGAUUCUUCAUAUU",
        )
        s = "GAAAAUUCAUUUUCUUUGGACUUUCUCUGAAAUCCGAGUCCUAGGAAAGAUGCGUGAGAUUCUUCAUAUU"
        rna = Seq.Seq(s)
        dna = rna.back_transcribe()
        self.assertEqual(
            dna,
            "GAAAATTCATTTTCTTTGGACTTTCTCTGAAATCCGAGTCCTAGGAAAGATGCGTGAGATTCTTCATATT",
        )

    def test_translation(self):
        s = ""
        dna = Seq.Seq(s)
        protein = dna.translate(to_stop=True)
        self.assertEqual(protein, "")
        s = "TAA"
        dna = Seq.Seq(s)
        protein = dna.translate(to_stop=True)
        self.assertEqual(protein, "")
        s = "GAAAATTCATTTTCTTTGGACTTTCTCTGAAATCCGAGTCCTAGGAAAGATGCGTGAGATTCTTCA"
        dna = Seq.Seq(s)
        protein = dna.translate(to_stop=True)
        self.assertEqual(protein, "ENSFSLDFL")
        s = "GAA"
        dna = Seq.Seq(s)
        protein = dna.translate(15, to_stop=True)
        self.assertEqual(protein, "E")
        s = "ATA"
        dna = Seq.Seq(s)
        protein = dna.translate("Vertebrate Mitochondrial", to_stop=True)
        self.assertEqual(protein, "M")
        s = "GAAAATTCATTTTCTTTGGACTTTCTCTGAAATCCGAGTCCTAGGAAAGATGCGTGAGATTCTTCATAT"
        dna = Seq.Seq(s)
        protein = dna.translate("SGC8", to_stop=True)
        self.assertEqual(protein, "ENSFSLDFLWNPSPSNDAWDSSY")

    def test_dna_rna_translation(self):
        s = "TCAAAAAGGTGCATCTAGATG"
        dna = Seq.Seq(s)
        protein = dna.translate(to_stop=True)
        self.assertEqual(protein, "SKRCI")
        gapped_protein = dna.translate()
        self.assertEqual(gapped_protein, "SKRCI*M")
        # The table used here has "AGG" as a stop codon:
        p2 = dna.translate(table=2, to_stop=True)
        self.assertEqual(p2, "SK")
        p2 = dna.translate(table=2)
        self.assertEqual(p2, "SK*CI*M")
        p2 = dna.translate(table=2, stop_symbol="+")
        self.assertEqual(p2, "SK+CI+M")
        r = s.replace("T", "U")
        rna = Seq.Seq(r)
        protein = rna.translate(to_stop=True)
        self.assertEqual(protein, "SKRCI")
        gapped_protein = rna.translate()
        self.assertEqual(gapped_protein, "SKRCI*M")

    def test_ambiguous(self):
        s = "RATGATTARAATYTA"
        dna = Seq.Seq(s)
        protein = dna.translate("Vertebrate Mitochondrial")
        self.assertEqual(protein, "BD*NL")
        stop_protein = dna.translate("SGC1", to_stop=True)
        self.assertEqual(stop_protein, "BD")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
