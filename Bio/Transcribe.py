import string

from Bio import Alphabet, Seq
from Bio.Alphabet import IUPAC

class Transcribe:
    def __init__(self, dna_alphabet, rna_alphabet):
        self.dna_alphabet = dna_alphabet
        self.rna_alphabet = rna_alphabet
        
    def transcribe(self, dna):
        assert dna.alphabet == self.dna_alphabet, \
               "transcribe has the wrong DNA alphabet"
        s = dna.data
        return Seq.Seq(string.replace(s, "T", "U"), self.rna_alphabet)
    def back_transcribe(self, rna):
        assert rna.alphabet == self.rna_alphabet, \
               "back transcribe has the wrong RNA alphabet"
        s = rna.data
        return Seq.Seq(string.replace(s, "U", "T"), self.dna_alphabet)

generic_transcriber = Transcribe(Alphabet.generic_dna,
                                 Alphabet.generic_rna)
ambiguous_transcriber = Transcribe(IUPAC.ambiguous_dna,
                                   IUPAC.ambiguous_rna)
unambiguous_transcriber = Transcribe(IUPAC.unambiguous_dna,
                                     IUPAC.unambiguous_rna)
