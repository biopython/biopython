import string, re

__all__ = [
    'IUPAC',
    ]

# This is used by sequences which contain a finite number of similar
# words.

class Alphabet:
    size = None     # no fixed size for words
    letters = None  # no fixed alphabet; implement as a list-like
                    # interface,
    def __repr__(self):
        return self.__class__.__name__ + "()"

    def contains(self, other):
        return isinstance(other, self.__class__)

generic_alphabet = Alphabet()

class SingleLetterAlphabet(Alphabet):
    size = 1
    letters = None   # string of all letters in the alphabet

single_letter_alphabet = SingleLetterAlphabet()

########### Protein

class ProteinAlphabet(SingleLetterAlphabet):
    pass

generic_protein = ProteinAlphabet()

########### DNA
class NucleotideAlphabet(SingleLetterAlphabet):
    pass

generic_nucleotide = NucleotideAlphabet()

class DNAAlphabet(NucleotideAlphabet):
    pass

generic_dna = DNAAlphabet()


########### RNA

class RNAAlphabet(NucleotideAlphabet):
    pass

generic_rna = RNAAlphabet()



########### Other per-sequence encodings

class SecondaryStructure(SingleLetterAlphabet):
    letters = "HSTC"

class ThreeLetterProtein(Alphabet):
    size = 3
    letters = [
        "Ala", "Asx", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile",
        "Lys", "Leu", "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr",
        "Sec", "Val", "Trp", "Xaa", "Tyr", "Glx",
        ]
        
###### Non per-sequence modifications

# (These are Decorator classes)

class AlphabetEncoder:
    def __init__(self, alphabet, new_letters):
        self.alphabet = alphabet
        self.new_letters = new_letters
        if alphabet.letters is not None:
            self.letters = alphabet.letters + new_letters
        else:
            self.letters = None
    def __getattr__(self, key):
        if key[:2] == "__" and key[-2:] == "__":
            raise AttributeError(key)
        return getattr(self.alphabet, key)

    def __repr__(self):
        return "%s(%r, %r)" % (self.__class__.__name__, self.alphabet,
                               self.new_letters)

    def contains(self, other):
        return 0
    
class Gapped(AlphabetEncoder):
    gap_char = '-'
    def __init__(self, alphabet, gap_char = gap_char):
        AlphabetEncoder.__init__(self, alphabet, gap_char)

    def contains(self, other):
        return other.gap_char == self.gap_char and \
               self.alphabet.contains(other.alphabet)
               
class HasStopCodon(AlphabetEncoder):
    stop_symbol = "*"
    def __init__(self, alphabet, stop_symbol = stop_symbol):
        AlphabetEncoder.__init__(self, alphabet, stop_symbol)
    def __cmp__(self, other):
        x = cmp(self.alphabet, other.alphabet)
        if x == 0:
            return cmp(self.stop_symbol, other.stop_symbol)
        return x

    def contains(self, other):
        return other.stop_symbol == self.stop_symbol and \
               self.alphabet.contains(other.alphabet)
