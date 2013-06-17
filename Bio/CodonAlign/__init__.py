# Copyright 2013 by Zheng Ruan.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Code for dealing with Codon Alignment.

CodonAlignment class is interited from MultipleSeqAlignment class. This is
the core class to deal with codon alignment in biopython.

"""
__docformat__ = "epytext en"  # Don't just use plain text in epydoc API pages!

#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord
from Bio import Alphabet
from Bio.Alphabet import IUPAC
from Bio.Align import MultipleSeqAlignment

def get_codon_alphabet(alphabet, gap="-", stop="*"):
    """function to get alignment alphabet for codon alignment. Only nucleotide
    alphabet is accepted. Raise an error with the type of alphabet is
    incompatible.
    """
    from Bio.Alphabet import NucleotideAlphabet
    if isinstance(alphabet, NucleotideAlphabet):
        alpha = alphabet
        if gap:
            alpha = Alphabet.Gapped(alpha, gap_char=gap)
        if stop:
            alpha = Alphabet.HasStopCodon(alpha, stop_symbol=stop)
    else:
        raise TypeError("Only Nuclteotide Alphabet is accepted!")
    return alpha

default_alphabet = get_codon_alphabet(IUPAC.unambiguous_dna)


class CodonAlignment(MultipleSeqAlignment):
    """Codon Alignment class that inherits from MultipleSeqAlignment.
    """
    def __init__(self, records, name=None, alphabet=default_alphabet):

        MultipleSeqAlignment.__init__(self, records, alphabet=alphabet)

        # check the type of the alignment to be nucleotide
        for rec in self:
            rec.seq.alphabet = alphabet
            assert Alphabet._verify_alphabet(rec.seq), \
                "%s is incompatible with the %s\n" % (rec.id, str(alphabet))
    
        # check the length of the alignment to be a triple
        assert self.get_alignment_length() % 3 == 0, \
            "Alignment length is not a triple number"


if __name__ == "__main__":
#    from Bio._utils import run_doctest
#    run_doctest()

     # These command should run correctly
     from Bio.Alphabet import generic_dna
     from Bio.Seq import Seq
     from Bio.SeqRecord import SeqRecord
     a = SeqRecord(Seq("AAAACGTC", generic_dna), id="Alpha")
     b = SeqRecord(Seq("AAA-CGTC", generic_dna), id="Beta")
     c = SeqRecord(Seq("AAAAGGTG", generic_dna), id="Gamma")
     MultipleSeqAlignment([a, b, c])
     align = CodonAlignment([a, b, c])

