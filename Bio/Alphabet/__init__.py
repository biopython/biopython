# Copyright 2000-2002 by Andrew Dalke.
# Revisions copyright 2007-2010 by Peter Cock.
# All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Alphabets used in Seq objects etc to declare sequence type and letters (OBSOLETE).

This is used by sequences which contain a finite number of similar words.

The design of Bio.Aphabet included a number of historic design choices
which, with the benefit of hindsight, were regretable. While the details
remain to be agreed, we intend to remove or replace Bio.Alphabet in 2020.
Please avoid using this module explicitly in your code. See also:
https://github.com/biopython/biopython/issues/2046
"""

import warnings

warnings.warn("We intend to remove or replace Bio.Alphabet in 2020, "
              "ideally avoid using it explicitly in your code. Please "
              "get in touch if you will be adversely affected by this. "
              "https://github.com/biopython/biopython/issues/2046",
              PendingDeprecationWarning)


class Alphabet(object):
    """Generic alphabet base class.

    This class is used as a base class for other types of alphabets.

    Attributes:
        - letters - list-like object containing the letters of the alphabet.
               Usually it is a string when letters are single characters.
        - size    - size of the alphabet's letters (e.g. 1 when letters are
               single characters).

    """

    size = None     # default to no fixed size for words
    letters = None  # default to no fixed alphabet
    # In general, a list-like object. However,
    # assuming letters are single characters, use a
    # string. This is expected for use with Seq like
    # objects.

    def __repr__(self):
        """Represent the alphabet class as a string for debugging."""
        return self.__class__.__name__ + "()"

    def contains(self, other):
        """Test if the other alphabet is contained in this one (OBSOLETE?).

        Returns a boolean.  This relies on the Alphabet subclassing
        hierarchy only, and does not check the letters property.
        This isn't ideal, and doesn't seem to work as intended
        with the AlphabetEncoder classes.
        """
        return isinstance(other, self.__class__)

    def _case_less(self):
        """Return a case-less variant of the current alphabet (PRIVATE)."""
        # TODO - remove this method by dealing with things in subclasses?
        if isinstance(self, ProteinAlphabet):
            return generic_protein
        elif isinstance(self, DNAAlphabet):
            return generic_dna
        elif isinstance(self, RNAAlphabet):
            return generic_rna
        elif isinstance(self, NucleotideAlphabet):
            return generic_nucleotide
        elif isinstance(self, SingleLetterAlphabet):
            return single_letter_alphabet
        else:
            return generic_alphabet

    def _upper(self):
        """Return an upper case variant of the current alphabet (PRIVATE)."""
        if not self.letters or self.letters == self.letters.upper():
            # Easy case, no letters or already upper case!
            return self
        else:
            # TODO - Raise NotImplementedError and handle via subclass?
            return self._case_less()

    def _lower(self):
        """Return a lower case variant of the current alphabet (PRIVATE)."""
        if not self.letters or self.letters == self.letters.lower():
            # Easy case, no letters or already lower case!
            return self
        else:
            # TODO - Raise NotImplementedError and handle via subclass?
            return self._case_less()


generic_alphabet = Alphabet()


class SingleLetterAlphabet(Alphabet):
    """Generic alphabet with letters of size one."""

    size = 1
    letters = None   # string of all letters in the alphabet


single_letter_alphabet = SingleLetterAlphabet()

# ########## Protein


class ProteinAlphabet(SingleLetterAlphabet):
    """Generic single letter protein alphabet."""

    pass


generic_protein = ProteinAlphabet()

# ########## DNA


class NucleotideAlphabet(SingleLetterAlphabet):
    """Generic single letter nucleotide alphabet."""

    pass


generic_nucleotide = NucleotideAlphabet()


class DNAAlphabet(NucleotideAlphabet):
    """Generic single letter DNA alphabet."""

    pass


generic_dna = DNAAlphabet()


# ########## RNA


class RNAAlphabet(NucleotideAlphabet):
    """Generic single letter RNA alphabet."""

    pass


generic_rna = RNAAlphabet()

# ########## Other per-sequence encodings


class SecondaryStructure(SingleLetterAlphabet):
    """Alphabet used to describe secondary structure.

    Letters are 'H' (helix), 'S' (strand), 'T' (turn) and 'C' (coil).
    """

    letters = "HSTC"


class ThreeLetterProtein(Alphabet):
    """Three letter protein alphabet."""

    size = 3
    letters = [
        "Ala", "Asx", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile",
        "Lys", "Leu", "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr",
        "Sec", "Val", "Trp", "Xaa", "Tyr", "Glx",
        ]

    def _upper(self):
        raise NotImplementedError("We don't have an uppercase three letter "
                                  "protein alphabet.")

    def _lower(self):
        raise NotImplementedError("We don't have a lowercase three letter "
                                  "protein alphabet.")

# ##### Non per-sequence modifications

# (These are Decorator classes)


class AlphabetEncoder(object):
    """A class to construct a new, extended alphabet from an existing one."""

    def __init__(self, alphabet, new_letters):
        """Initialize the class."""
        self.alphabet = alphabet
        self.new_letters = new_letters
        if alphabet.letters is not None:
            self.letters = alphabet.letters + new_letters
        else:
            self.letters = None

    def __getattr__(self, key):
        """Proxy method for accessing attributes of the wrapped alphabet."""
        if key[:2] == "__" and key[-2:] == "__":
            raise AttributeError(key)
        return getattr(self.alphabet, key)

    def __repr__(self):
        """Represent the alphabet encoder class as a string for debugging."""
        return "%s(%r, %r)" % (self.__class__.__name__, self.alphabet,
                               self.new_letters)

    def contains(self, other):
        """Test if the other alphabet is contained in this one (OBSOLETE?).

        This is isn't implemented for the base AlphabetEncoder,
        which will always return 0 (False).
        """
        return 0

    def _upper(self):
        """Return an upper case variant of the current alphabet (PRIVATE)."""
        return AlphabetEncoder(self.alphabet._upper(),
                               self.new_letters.upper())

    def _lower(self):
        """Return a lower case variant of the current alphabet (PRIVATE)."""
        return AlphabetEncoder(self.alphabet._lower(),
                               self.new_letters.lower())


class Gapped(AlphabetEncoder):
    """Alphabets which contain a gap character."""

    def __init__(self, alphabet, gap_char="-"):
        """Initialize the class."""
        AlphabetEncoder.__init__(self, alphabet, gap_char)
        self.gap_char = gap_char

    def contains(self, other):
        """Test if the other alphabet is contained in this one (OBSOLETE?).

        Returns a boolean.  This relies on the Alphabet subclassing
        hierarchy, and attempts to check the gap character.  This fails
        if the other alphabet does not have a gap character!
        """
        return (other.gap_char == self.gap_char and
                self.alphabet.contains(other.alphabet))

    def _upper(self):
        """Return an upper case variant of the current alphabet (PRIVATE)."""
        return Gapped(self.alphabet._upper(), self.gap_char.upper())

    def _lower(self):
        """Return a lower case variant of the current alphabet (PRIVATE)."""
        return Gapped(self.alphabet._lower(), self.gap_char.lower())


class HasStopCodon(AlphabetEncoder):
    """Alphabets which contain a stop symbol."""

    def __init__(self, alphabet, stop_symbol="*"):
        """Initialize the class."""
        AlphabetEncoder.__init__(self, alphabet, stop_symbol)
        self.stop_symbol = stop_symbol

    def contains(self, other):
        """Test if the other alphabet is contained in this one (OBSOLETE?).

        Returns a boolean.  This relies on the Alphabet subclassing
        hierarchy, and attempts to check the stop symbol.  This fails
        if the other alphabet does not have a stop symbol!
        """
        return (other.stop_symbol == self.stop_symbol and
                self.alphabet.contains(other.alphabet))

    def _upper(self):
        """Return an upper case variant of the current alphabet (PRIVATE)."""
        return HasStopCodon(self.alphabet._upper(), self.stop_symbol.upper())

    def _lower(self):
        """Return a lower case variant of the current alphabet (PRIVATE)."""
        return HasStopCodon(self.alphabet._lower(), self.stop_symbol.lower())


def _get_base_alphabet(alphabet):
    """Return the non-gapped non-stop-codon Alphabet object (PRIVATE)."""
    a = alphabet
    while isinstance(a, AlphabetEncoder):
        a = a.alphabet
    if not isinstance(a, Alphabet):
        raise TypeError("Invalid alphabet found, %s." % repr(a))
    return a


def _ungap(alphabet):
    """Return the alphabet without any gap encoder (PRIVATE)."""
    # TODO - Handle via method of the objects?
    if not hasattr(alphabet, "gap_char"):
        return alphabet
    elif isinstance(alphabet, Gapped):
        return alphabet.alphabet
    elif isinstance(alphabet, HasStopCodon):
        return HasStopCodon(_ungap(alphabet.alphabet),
                            stop_symbol=alphabet.stop_symbol)
    elif isinstance(alphabet, AlphabetEncoder):
        return AlphabetEncoder(_ungap(alphabet.alphabet),
                               letters=alphabet.letters)
    else:
        raise NotImplementedError


def _consensus_base_alphabet(alphabets):
    """Return a common but often generic base alphabet object (PRIVATE).

    This throws away any AlphabetEncoder information, e.g. Gapped alphabets.

    Note that DNA+RNA -> Nucleotide, and Nucleotide+Protein-> generic single
    letter.  These DO NOT raise an exception!
    """
    common = None
    for alpha in alphabets:
        a = _get_base_alphabet(alpha)
        if common is None:
            common = a
        elif common == a:
            pass
        elif isinstance(a, common.__class__):
            pass
        elif isinstance(common, a.__class__):
            common = a
        elif (isinstance(a, NucleotideAlphabet) and
              isinstance(common, NucleotideAlphabet)):
            # e.g. Give a mix of RNA and DNA alphabets
            common = generic_nucleotide
        elif (isinstance(a, SingleLetterAlphabet) and
              isinstance(common, SingleLetterAlphabet)):
            # This is a pretty big mis-match!
            common = single_letter_alphabet
        else:
            # We have a major mis-match... take the easy way out!
            return generic_alphabet
    if common is None:
        # Given NO alphabets!
        return generic_alphabet
    return common


def _consensus_alphabet(alphabets):
    """Return a common but often generic alphabet object (PRIVATE).

        >>> from Bio.Alphabet import IUPAC
        >>> _consensus_alphabet([IUPAC.extended_protein, IUPAC.protein])
        ExtendedIUPACProtein()
        >>> _consensus_alphabet([generic_protein, IUPAC.protein])
        ProteinAlphabet()

    Note that DNA+RNA -> Nucleotide, and Nucleotide+Protein-> generic single
    letter.  These DO NOT raise an exception!

        >>> _consensus_alphabet([generic_dna, generic_nucleotide])
        NucleotideAlphabet()
        >>> _consensus_alphabet([generic_dna, generic_rna])
        NucleotideAlphabet()
        >>> _consensus_alphabet([generic_dna, generic_protein])
        SingleLetterAlphabet()
        >>> _consensus_alphabet([single_letter_alphabet, generic_protein])
        SingleLetterAlphabet()

    This is aware of Gapped and HasStopCodon and new letters added by
    other AlphabetEncoders.  This WILL raise an exception if more than
    one gap character or stop symbol is present.

        >>> from Bio.Alphabet import IUPAC
        >>> _consensus_alphabet([Gapped(IUPAC.extended_protein),
        ...                     HasStopCodon(IUPAC.protein)])
        HasStopCodon(Gapped(ExtendedIUPACProtein(), '-'), '*')
        >>> _consensus_alphabet([Gapped(IUPAC.protein, "-"),
        ...                     Gapped(IUPAC.protein, "=")])
        Traceback (most recent call last):
            ...
        ValueError: More than one gap character present
        >>> _consensus_alphabet([HasStopCodon(IUPAC.protein, "*"),
        ...                     HasStopCodon(IUPAC.protein, "+")])
        Traceback (most recent call last):
            ...
        ValueError: More than one stop symbol present
    """
    base = _consensus_base_alphabet(alphabets)
    gap = None
    stop = None
    new_letters = ""
    for alpha in alphabets:
        # Gaps...
        if not hasattr(alpha, "gap_char"):
            pass
        elif gap is None:
            gap = alpha.gap_char
        elif gap == alpha.gap_char:
            pass
        else:
            raise ValueError("More than one gap character present")
        # Stops...
        if not hasattr(alpha, "stop_symbol"):
            pass
        elif stop is None:
            stop = alpha.stop_symbol
        elif stop == alpha.stop_symbol:
            pass
        else:
            raise ValueError("More than one stop symbol present")
        # New letters...
        if hasattr(alpha, "new_letters"):
            for letter in alpha.new_letters:
                if letter not in new_letters and letter != gap \
                   and letter != stop:
                    new_letters += letter

    alpha = base
    if new_letters:
        alpha = AlphabetEncoder(alpha, new_letters)
    if gap:
        alpha = Gapped(alpha, gap_char=gap)
    if stop:
        alpha = HasStopCodon(alpha, stop_symbol=stop)
    return alpha


def _check_type_compatible(alphabets):
    """Return True except for DNA+RNA or Nucleotide+Protein (PRIVATE).

        >>> _check_type_compatible([generic_dna, generic_nucleotide])
        True
        >>> _check_type_compatible([generic_dna, generic_rna])
        False
        >>> _check_type_compatible([generic_dna, generic_protein])
        False
        >>> _check_type_compatible([single_letter_alphabet, generic_protein])
        True

    This relies on the Alphabet subclassing hierarchy.  It does not
    check things like gap characters or stop symbols.
    """
    dna, rna, nucl, protein = False, False, False, False
    for alpha in alphabets:
        a = _get_base_alphabet(alpha)
        if isinstance(a, DNAAlphabet):
            dna = True
            nucl = True
            if rna or protein:
                return False
        elif isinstance(a, RNAAlphabet):
            rna = True
            nucl = True
            if dna or protein:
                return False
        elif isinstance(a, NucleotideAlphabet):
            nucl = True
            if protein:
                return False
        elif isinstance(a, ProteinAlphabet):
            protein = True
            if nucl:
                return False
    return True


def _verify_alphabet(sequence):
    """Check all letters in sequence are in the alphabet (PRIVATE).

        >>> from Bio.Seq import Seq
        >>> from Bio.Alphabet import IUPAC
        >>> my_seq = Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF",
        ...              IUPAC.protein)
        >>> _verify_alphabet(my_seq)
        True

        This example has an X, which is not in the IUPAC protein alphabet
        (you should be using the IUPAC extended protein alphabet):

        >>> bad_seq = Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVFX",
        ...                IUPAC.protein)
        >>> _verify_alphabet(bad_seq)
        False

    This replaces Bio.utils.verify_alphabet() since we are deprecating
    that. Potentially this could be added to the Alphabet object, and
    I would like it to be an option when creating a Seq object... but
    that might slow things down.
    """
    letters = sequence.alphabet.letters
    if not letters:
        raise ValueError("Alphabet does not define letters.")
    for letter in sequence:
        if letter not in letters:
            return False
    return True
