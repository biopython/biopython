# Copyright 2000 by Andrew Dalke.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Miscellaneous functions for dealing with sequences (DEPRECATED).

This module is deprecated, and is expected to be removed in the next release.
If you use this module, please contact the Biopython developers via the
mailing lists.
"""

import warnings
import Bio
warnings.warn("Bio.utils has been deprecated, and we intend to remove it in "
              "the next release of Biopython.", Bio.BiopythonDeprecationWarning)
raise NotImplementedError("Error to force a traceback")

import Seq
import Alphabet

from Bio.Alphabet import _verify_alphabet as verify_alphabet

from PropertyManager import default_manager

def ungap(seq):
    """given a sequence with gap encoding, return the ungapped sequence"""
    #TODO - Fix this?  It currently assumes the outmost AlphabetEncoder
    #is for the gap.  Consider HasStopCodon(Gapped(Protein())) as a test case.
    warnings.warn("Bio.utils has been deprecated, and we intend to remove it "
                  "in the next release of Biopython. Instead of function "
                  "Bio.utils.ungap please use the ungap method of the Seq "
                  "object (added in Biopython 1.53).",
                  Bio.BiopythonDeprecationWarning)
    gap = seq.gap_char
    letters = []
    for c in seq:
        if c != gap:
            letters.append(c)
    return Seq.Seq("".join(letters), seq.alphabet.alphabet)

def count_monomers(seq):
    dict = {}
    for c in seq.alphabet.letters:
        dict[c] = seq.count(c)
    return dict

def percent_monomers(seq):
    dict2 = {}
    seq_len = len(seq)
    dict = count_monomers(seq)
    for m in dict:
        dict2[m] = dict[m] * 100. / seq_len
    return dict2

def sum(seq, table, zero = 0.0):
    total = zero
    for c in getattr(seq, "data", seq):
        total = total + table[c]
    return total

# For ranged addition
def sum_2ple(seq, table, zero = (0.0, 0.0)):
    x, y = zero
    data = getattr(seq, "data", seq)
    for c in data:
        x2, y2 = table[c]
        x = x + x2
        y = y + y2
    return (x, y)

def total_weight(seq, weight_table = None):
    if weight_table is None:
        weight_table = default_manager.resolve(seq.alphabet, "weight_table")
    return sum(seq, weight_table)

def total_weight_range(seq, weight_table = None):
    if weight_table is None:
        weight_table = default_manager.resolve(seq.alphabet, "weight_range_table")
    return sum_2ple(seq, weight_table)

def reduce_sequence(seq, reduction_table,new_alphabet=None):
   """ given an amino-acid sequence, return it in reduced alphabet form based
       on the letter-translation table passed. Some "standard" tables are in
       Alphabet.Reduced.
       seq: a Seq.Seq type sequence
       reduction_table: a dictionary whose keys are the "from" alphabet, and values
       are the "to" alphabet"""
   if new_alphabet is None:
      new_alphabet = Alphabet.single_letter_alphabet
      new_alphabet.letters = ''
      for letter in reduction_table:
         new_alphabet.letters += letter
      new_alphabet.size = len(new_alphabet.letters)
   new_seq = Seq.Seq('',new_alphabet)
   for letter in seq:
      new_seq += reduction_table[letter]
   return new_seq


