# Copyright 2004 by Iddo Friedberg.
# All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Reduced alphabets which lump together several amino-acids into one letter.

Reduced (redundant or simplified) alphabets are used to represent protein
sequences using an alternative alphabet which lumps together several
amino-acids into one letter, based on physico-chemical traits. For example,
all the aliphatics (I,L,V) are usually quite interchangeable, so many sequence
studies group them into one letter

Examples of reduced alphabets are available in:

http://viscose.herokuapp.com/html/alphabets.html

The Murphy tables are from here:

Murphy L.R., Wallqvist A, Levy RM. (2000) Simplified amino acid
alphabets for protein fold recognition and implications for folding.
Protein Eng. 13(3):149-152

These alphabets have been used with Bio.utils.reduce_sequence, which has been
removed from Biopython. You can use this is alphabets and tables like this:

    >>> from Bio.Seq import Seq
    >>> from Bio import Alphabet
    >>> from Bio.Alphabet import Reduced
    >>> my_protein = Seq('MAGSKEWKRFCELTINEA', Alphabet.ProteinAlphabet())

Now, we convert this sequence into a sequence which only recognizes polar (P)
or hydrophobic (H) residues:

    >>> new_protein = Seq('', Alphabet.Reduced.HPModel())
    >>> for aa in my_protein:
    ...     new_protein += Alphabet.Reduced.hp_model_tab[aa]
    >>> new_protein
    Seq('HPPPPPHPPHHPHPHPPP', HPModel())

The following Alphabet classes are available:

 - Murphy15: Maps 20 amino acids to 15, use murphy_15_tab for conversion,
             ambiguous letters: L: LVIM, F: FY, K: KR
 - Murphy10: Maps 20 amino acids to 10, use murphy_10_tab for conversion,
             ambiguous letters: L: LVIM, S: ST, F: FYW, E: EDNQ, K: KR
 - Murphy8: Maps 20 amino acids to 8, use murphy_8_tab for conversion,
            ambiguous letters: L: LVIMC, A: AG, S: ST, F: FYW, E: EDNQ,
            K: KR
 - Murphy4: Maps 20 amino acids to 4, use murphy_4_tab for conversion,
            ambiguous letters: L: LVIMC, A: AGSTP, F: FYW, E: EDNQKRH
 - HPModel: Groups amino acids as polar (hydrophilic) or hydrophobic
            (non-polar), use hp_model_tab for conversion,
            P: AGTSNQDEHRKP, H: CMFILVWY
 - PC5: Amino acids grouped according to 5 physico-chemical properties,
        use pc_5_table for conversion,
        A (Aliphatic): IVL, R (aRomatic): FYWH, C (Charged): KRDE, T (Tiny):
        GACS, D (Diverse): TMQNP
"""

from Bio import Alphabet


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
                 "H": "H"}


class Murphy15(Alphabet.ProteinAlphabet):
    """Reduced protein alphabet with 15 letters.

    Letters: A, C, D, E, G, H, N, P, Q, S, T, W,
             L(LVIM), F(FY), K(KR)
    """

    letters = "LCAGSTPFWEDNQKH"
    size = 1


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
    """Reduced protein alphabet with 10 letters.

    Letters: A, C, G, H, P, L(LVIM), S(ST), F(FYW),
             E(EDNQ), K(KR)
    """

    letters = "LCAGSPFEKH"
    size = 1


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
    """Reduced protein alphabet with 8 letters.

    Letters: H, P, L(LVIMC), A(AG), S(ST), F(FYW),
             E(EDNQ), K(KR)
    """

    letters = "LASPFEKH"
    size = 1


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
    """Reduced protein alphabet with 4 letters.

    Letters: L(LVIMC), A(AGSTP), F(FYW), E(EDNQKRH)
    """

    letters = "LAFE"
    size = 1


murphy_4 = Murphy4()

hp_model_tab = {"A": "P",  # Hydrophilic
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
    """Reduced protein alphabet with only two letters for polar or hydophobic.

    Letters: P (polar: AGTSNQDEHRKP), H (hydrophobic: CMFILVWY)
    """

    letters = "HP"
    size = 1


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
    """Reduced protein alphabet with 5 letters for physico-chemical properties.

    Letters: A (Aliphatic: IVL), R (aRomatic: FYWH), C (Charged: KRDE),
             T (Tiny: GACS), D (Diverse: TMQNP)
    """

    letters = "ARCTD"
    size = 1


pc5 = PC5()
