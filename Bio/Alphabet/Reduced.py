# Copyright 2004 by Iddo Friedberg.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Reduced alphabets which lump together several amino-acids into one letter.

Reduced (redundant or simplified) alphabets are used to represent protein sequences using an
alternative alphabet which lumps together several amino-acids into one letter, based
on physico-chemical traits. For example, all the aliphatics (I,L,V) are usually
quite interchangeable, so many sequence studies group them into one letter

Examples of reduced alphabets are available in:

http://viscose.ifg.uni-muenster.de/html/alphabets.html

The Murphy tables are from here:

Murphy L.R., Wallqvist A, Levy RM. (2000) Simplified amino acid
alphabets for protein fold recognition and implications for folding.
Protein Eng. 13(3):149-152

Bio.utils.reduce_sequence is used to take a Protein alphabet, and reduce it using one of
the tables here, or a user-defined table.
"""

from Bio import Alphabet

__docformat__ = "restructuredtext en"

murphy_15_tab = {"L": "L",
                 "V": "L",
                 "I": "L",
                 "M": "L",
                 "C": "C",
                 "A": "A",
                 "G": "G",
                 "S": "S",
                 "T": "T",
                 "P": "P",
                 "F": "F",
                 "Y": "F",
                 "W": "W",
                 "E": "E",
                 "D": "D",
                 "N": "N",
                 "Q": "Q",
                 "K": "K",
                 "R": "K",
                 "H": "H",
}


class Murphy15(Alphabet.ProteinAlphabet):
    letters = "LCAGSTPFWEDNQKH"
    size = 15
murphy_15 = Murphy15()

murphy_10_tab = {"L": "L",
                 "V": "L",
                 "I": "L",
                 "M": "L",
                 "C": "C",
                 "A": "A",
                 "G": "G",
                 "S": "S",
                 "T": "S",
                 "P": "P",
                 "F": "F",
                 "Y": "F",
                 "W": "F",
                 "E": "E",
                 "D": "E",
                 "N": "E",
                 "Q": "E",
                 "K": "K",
                 "R": "K",
                 "H": "H"}


class Murphy10(Alphabet.ProteinAlphabet):
    letters = "LCAGSPFEKH"
    size = 10
murphy_10 = Murphy10()

murphy_8_tab = {"L": "L",
                "V": "L",
                "I": "L",
                "M": "L",
                "C": "L",
                "A": "A",
                "G": "A",
                "S": "S",
                "T": "S",
                "P": "P",
                "F": "F",
                "Y": "F",
                "W": "F",
                "E": "E",
                "D": "E",
                "N": "E",
                "Q": "E",
                "K": "K",
                "R": "K",
                "H": "H"}


class Murphy8(Alphabet.ProteinAlphabet):
    letters = "LASPFEKH"
    size = 8
murphy_8 = Murphy8()

murphy_4_tab = {"L": "L",
                "V": "L",
                "I": "L",
                "M": "L",
                "C": "L",
                "A": "A",
                "G": "A",
                "S": "A",
                "T": "A",
                "P": "A",
                "F": "F",
                "Y": "F",
                "W": "F",
                "E": "E",
                "D": "E",
                "N": "E",
                "Q": "E",
                "K": "E",
                "R": "E",
                "H": "E"}


class Murphy4(Alphabet.ProteinAlphabet):
    letters = "LAFE"
    size = 4
murphy_4 = Murphy4()

hp_model_tab = {"A": "P",   # Hydrophilic
                "G": "P",
                "T": "P",
                "S": "P",
                "N": "P",
                "Q": "P",
                "D": "P",
                "E": "P",
                "H": "P",
                "R": "P",
                "K": "P",
                "P": "P",
                "C": "H",  # Hydrophobic
                "M": "H",
                "F": "H",
                "I": "H",
                "L": "H",
                "V": "H",
                "W": "H",
                "Y": "H"}


class HPModel(Alphabet.ProteinAlphabet):
    letters = "HP"
    size = 2
hp_model = HPModel()

pc_5_table = {"I": "A",  # Aliphatic
              "V": "A",
              "L": "A",
              "F": "R",  # Aromatic
              "Y": "R",
              "W": "R",
              "H": "R",
              "K": "C",  # Charged
              "R": "C",
              "D": "C",
              "E": "C",
              "G": "T",  # Tiny
              "A": "T",
              "C": "T",
              "S": "T",
              "T": "D",  # Diverse
              "M": "D",
              "Q": "D",
              "N": "D",
              "P": "D"}


class PC5(Alphabet.ProteinAlphabet):
    letters = "ARCTD"
    size = 5
hp_model = HPModel()
