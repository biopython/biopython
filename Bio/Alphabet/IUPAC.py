# Copyright 2000-2001 by Andrew Dalke.
# Revisions copyright 2008 by Peter Cock.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Standard nucleotide and protein alphabets defined by IUPAC."""

from Bio import Alphabet
from Bio.Data import IUPACData

##################### Protein

# From the IUPAC definition at:
#   http://www.chem.qmw.ac.uk/iupac/AminoAcid/A2021.html#AA21

assert IUPACData.extended_protein_letters == IUPACData.extended_protein_letters.upper()
class ExtendedIUPACProtein(Alphabet.ProteinAlphabet):
    """Extended uppercase IUPAC protein single letter alphabet including X etc.

    In addition to the standard 20 single letter protein codes, this includes:
    
    B = "Asx";  Aspartic acid (R) or Asparagine (N)
    X = "Xxx";  Unknown or 'other' amino acid
    Z = "Glx";  Glutamic acid (E) or Glutamine (Q)
    J = "Xle";  Leucine (L) or Isoleucine (I), used in mass-spec (NMR)
    U = "Sec";  Selenocysteine
    O = "Pyl";  Pyrrolysine

    This alphabet is not intended to be used with X for Selenocysteine
    (an ad-hoc standard prior to the IUPAC adoption of U instead).
    """
    letters = IUPACData.extended_protein_letters

extended_protein = ExtendedIUPACProtein()

assert IUPACData.protein_letters == IUPACData.protein_letters.upper()
class IUPACProtein(ExtendedIUPACProtein):
    """Uppercase IUPAC protein single letter alphabet of the 20 standard amino acids."""
    letters = IUPACData.protein_letters

protein = IUPACProtein()

##################### DNA

# The next two are the IUPAC definitions, from:
#   http://www.chem.qmw.ac.uk/iubmb/misc/naseq.html
class IUPACAmbiguousDNA(Alphabet.DNAAlphabet):
    """Uppercase IUPAC ambiguous DNA."""
    letters = IUPACData.ambiguous_dna_letters

ambiguous_dna = IUPACAmbiguousDNA()

class IUPACUnambiguousDNA(IUPACAmbiguousDNA):
    """Uppercase IUPAC unambiguous DNA (letters GATC only)."""
    letters = IUPACData.unambiguous_dna_letters

unambiguous_dna = IUPACUnambiguousDNA()


# Also from the URL, but not part of the standard
class ExtendedIUPACDNA(Alphabet.DNAAlphabet):
    """Extended IUPAC DNA alphabet.

    In addition to the standard letter codes GATC, this includes:

    B = 5-bromouridine
    D = 5,6-dihydrouridine
    S = thiouridine
    W = wyosine
    """
    letters = IUPACData.extended_dna_letters

extended_dna = ExtendedIUPACDNA()

##################### RNA

class IUPACAmbiguousRNA(Alphabet.RNAAlphabet):
    """Uppercase IUPAC ambiguous RNA."""
    letters = IUPACData.ambiguous_rna_letters

ambiguous_rna = IUPACAmbiguousRNA()

class IUPACUnambiguousRNA(IUPACAmbiguousRNA):
    """Uppercase IUPAC unambiguous RNA (letters GAUC only)."""
    letters = IUPACData.unambiguous_rna_letters

unambiguous_rna = IUPACUnambiguousRNA()

# are there extended forms?
#class ExtendedIUPACRNA(Alphabet.RNAAlphabet):
#    letters = extended_rna_letters
#    #   B == 5-bromouridine
#    #   D == 5,6-dihydrouridine
#    #   S == thiouridine
#    #   W == wyosine
