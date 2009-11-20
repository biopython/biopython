"""Code to translate DNA or RNA into proteins (DEPRECATED).

Instead of Bio.Translate, for translation you are now encouraged to use the
Seq object's translate method, or the translate function in the Bio.Seq
module.  Translate-to-stop functionality is via an optional argument.

Bio.Seq does not offer any back-translation function like the one here. It
was concluded that a since a simple back-translation giving a Seq or python
string could only capture some of the possible back translations, there were
no practical uses for such a method/function.

This module is now deprecated, and will be removed in a future release of
Biopython.
"""
import warnings
warnings.warn("Bio.Translate and Bio.Transcribe are deprecated, and will be "\
              "removed in a future release of Biopython. Please use the "\
              "functions or object methods defined in Bio.Seq instead "\
              "(described in the tutorial). If you want to continue to use "\
              "this code, please get in contact with the Biopython developers "\
              "via the mailing lists to avoid its permanent removal from "
              +"Biopython.", \
              DeprecationWarning)

from Bio import Alphabet, Seq
from Bio.Data import CodonTable

class Translator:
    def __init__(self, table):
        self.table = table
        self._encoded = {}

    def __str__(self):
        return "Translator object\n" + str(self.table)

    def translate(self, seq, stop_symbol = "*"):
        #Allow different instances of the same class to be used:
        assert seq.alphabet.__class__ == \
               self.table.nucleotide_alphabet.__class__, \
               "cannot translate from given alphabet (have %s, need %s)" %\
               (seq.alphabet, self.table.nucleotide_alphabet)
        s = seq.data
        letters = []
        append = letters.append
        table = self.table
        get = table.forward_table.get
        n = len(seq)
        for i in range(0, n-n%3, 3):
            append(get(s[i:i+3], stop_symbol))

        # return with the correct alphabet encoding (cache the encoding)
        try:
            alphabet = self._encoded[stop_symbol]
        except KeyError:
            alphabet = Alphabet.HasStopCodon(table.protein_alphabet,
                                             stop_symbol)
            self._encoded[stop_symbol] = alphabet

        return Seq.Seq("".join(letters), alphabet)
                           
    def translate_to_stop(self, seq):
        # This doesn't have a stop encoding

        #Allow different instances of the same class to be used:
        assert seq.alphabet.__class__ == \
               self.table.nucleotide_alphabet.__class__, \
               "cannot translate from given alphabet (have %s, need %s)" %\
               (seq.alphabet, self.table.nucleotide_alphabet)
        s = seq.data
        letters = []
        append = letters.append
        table = self.table.forward_table
        n = len(seq)
        try:
            for i in range(0, n-n%3, 3):
                append(table[s[i:i+3]])
        except KeyError:
            # Stop at the first codon failure
            pass
        return Seq.Seq("".join(letters), self.table.protein_alphabet)

    def back_translate(self, seq):
        # includes the stop codon
        if not isinstance(seq.alphabet, Alphabet.HasStopCodon):
            return self._back_translate_no_stop(seq)
        assert seq.alphabet.alphabet == self.table.protein_alphabet, \
               "cannot back translate from the given alphabet (%s)" % \
               seq.alphabet.alphabet
        s = seq.data
        letter = seq.alphabet.stop_symbol
        letters = []
        append = letters.append
        table = self.table.back_table
        for c in seq.data:
            if c == letter:
                append(table[None])
            else:
                append(table[c])
        return Seq.Seq("".join(letters),
                       self.table.nucleotide_alphabet)

    def _back_translate_no_stop(self, seq):
        # does not allow a stop codon
        assert seq.alphabet == self.table.protein_alphabet, \
               "cannot back translate from the given alphabet (%s)" % \
               seq.alphabet
        s = seq.data
        letters = []
        append = letters.append
        table = self.table.back_table
        for c in seq.data:
            append(table[c])
        return Seq.Seq("".join(letters),
                       self.table.nucleotide_alphabet)

unambiguous_dna_by_name = {}
for key, value in CodonTable.unambiguous_dna_by_name.items():
    unambiguous_dna_by_name[key] = Translator(value)
unambiguous_dna_by_id = {}
for key, value in CodonTable.unambiguous_dna_by_id.items():
    unambiguous_dna_by_id[key] = Translator(value)

unambiguous_rna_by_name = {}
for key, value in CodonTable.unambiguous_rna_by_name.items():
    unambiguous_rna_by_name[key] = Translator(value)
unambiguous_rna_by_id = {}
for key, value in CodonTable.unambiguous_rna_by_id.items():
    unambiguous_rna_by_id[key] = Translator(value)

# XXX Ambiguous - can be done the same except for stop codons!
ambiguous_dna_by_name = {}
for key, value in CodonTable.ambiguous_dna_by_name.items():
    ambiguous_dna_by_name[key] = Translator(value)
ambiguous_dna_by_id = {}
for key, value in CodonTable.ambiguous_dna_by_id.items():
    ambiguous_dna_by_id[key] = Translator(value)

ambiguous_rna_by_name = {}
for key, value in CodonTable.ambiguous_rna_by_name.items():
    ambiguous_rna_by_name[key] = Translator(value)
ambiguous_rna_by_id = {}
for key, value in CodonTable.ambiguous_rna_by_id.items():
    ambiguous_rna_by_id[key] = Translator(value)
