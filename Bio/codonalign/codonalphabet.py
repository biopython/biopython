# Copyright 2013 by Zheng Ruan (zruan1991@gmail.com).
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Code for Codon Alphabet.

CodonAlphabet class is inherited from Alphabet class. It is an
alphabet for CodonSeq class.
"""

import copy
try:
    from itertools import izip
except ImportError:
    izip = zip

from Bio.Alphabet import IUPAC, Gapped, HasStopCodon, Alphabet
from Bio.Data.CodonTable import generic_by_id

default_codon_table = copy.deepcopy(generic_by_id[1])


class CodonAlphabet(Alphabet):
    """Generic Codon Alphabet with a size of three."""

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
        letters.append(gap_char * 3)
    generic_codon_alphabet = CodonAlphabet()
    generic_codon_alphabet.letters = letters
    generic_codon_alphabet.gap_char = '-'
    generic_codon_alphabet.names = codon_table.names
    return generic_codon_alphabet


default_codon_alphabet = get_codon_alphabet(default_codon_table)


if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()
