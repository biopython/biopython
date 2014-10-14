# Copyright 2013 by Zheng Ruan (zruan1991@gmail.com).
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Code for Codon Alphabet.

CodonAlphabet class is interited from Alphabet class. It is an
alphabet for CodonSeq class.
"""
__docformat__ = "restructuredtext en"  # Don't just use plain text in epydoc API pages!


import copy
try:
    from itertools import izip
except ImportError:
    izip = zip

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC, Gapped, HasStopCodon, Alphabet, generic_dna
from Bio.Data.CodonTable import generic_by_id

default_codon_table = copy.deepcopy(generic_by_id[1])


def get_codon_alphabet(alphabet, gap="-", stop="*"):
    """Gets alignment alphabet for codon alignment.
    
    Only nucleotide alphabet is accepted. Raise an error when the type of 
    alphabet is incompatible.
    """
    from Bio.Alphabet import NucleotideAlphabet
    if isinstance(alphabet, NucleotideAlphabet):
        alpha = alphabet
        if gap:
            alpha = Gapped(alpha, gap_char=gap)
        if stop:
            alpha = HasStopCodon(alpha, stop_symbol=stop)
    else:
        raise TypeError("Only Nuclteotide Alphabet is accepted!")
    return alpha

default_alphabet = get_codon_alphabet(IUPAC.unambiguous_dna)


class CodonAlphabet(Alphabet):
    """Generic Codon Alphabet with a size of three"""
    size = 3
    letters = None
    name = ''

    def __repr__(self):
        return "%s(%s)" % (self.__class__.__name__, self.names[0])


def get_codon_alphabet(codon_table, gap_char="-"):
    letters = list(codon_table.forward_table.keys())
    letters.extend(codon_table.stop_codons)
    letters.extend(codon_table.start_codons)
    if gap_char:
        letters.append(gap_char*3)
    generic_codon_alphabet = CodonAlphabet()
    generic_codon_alphabet.letters = letters
    generic_codon_alphabet.gap_char = '-'
    generic_codon_alphabet.names = codon_table.names
    return generic_codon_alphabet

default_codon_alphabet = get_codon_alphabet(default_codon_table)


if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()
