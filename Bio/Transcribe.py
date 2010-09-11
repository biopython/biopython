"""Code to transcribe DNA into RNA or back (DEPRECATED).

You are now encouraged to use the Seq object methods or the functions
in Bio.Seq instead.

This module is now deprecated, and will be removed in a future release of
Biopython.
"""
import warnings
import Bio
warnings.warn("Bio.Translate and Bio.Transcribe are deprecated, and will be "\
              "removed in a future release of Biopython. Please use the "\
              "functions or object methods defined in Bio.Seq instead "\
              "(described in the tutorial). If you want to continue to use "\
              "this code, please get in contact with the Biopython developers "\
              "via the mailing lists to avoid its permanent removal from "
              +"Biopython.", \
              Bio.BiopythonDeprecationWarning)


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
        return Seq.Seq(s.replace("T", "U"), self.rna_alphabet)
    def back_transcribe(self, rna):
        assert rna.alphabet == self.rna_alphabet, \
               "back transcribe has the wrong RNA alphabet"
        s = rna.data
        return Seq.Seq(s.replace("U", "T"), self.dna_alphabet)

generic_transcriber = Transcribe(Alphabet.generic_dna,
                                 Alphabet.generic_rna)
ambiguous_transcriber = Transcribe(IUPAC.ambiguous_dna,
                                   IUPAC.ambiguous_rna)
unambiguous_transcriber = Transcribe(IUPAC.unambiguous_dna,
                                     IUPAC.unambiguous_rna)
