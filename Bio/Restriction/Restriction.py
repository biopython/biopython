#!/usr/bin/env python
#
#      Restriction Analysis Libraries.
#      Copyright (C) 2004. Frederic Sohm.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#

"""Restriction Enzyme classes.

Notes about the diverses class of the restriction enzyme implementation::

            RestrictionType is the type of all restriction enzymes.
        -----------------------------------------------------------------------
            AbstractCut implements some methods that are common to all enzymes.
        -----------------------------------------------------------------------
            NoCut, OneCut,TwoCuts   represent the number of double strand cuts
                                    produced by the enzyme.
                                    they correspond to the 4th field of the
                                    rebase record emboss_e.NNN.
                    0->NoCut    : the enzyme is not characterised.
                    2->OneCut   : the enzyme produce one double strand cut.
                    4->TwoCuts  : two double strand cuts.
        -----------------------------------------------------------------------
            Meth_Dep, Meth_Undep    represent the methylation susceptibility to
                                    the enzyme.
                                    Not implemented yet.
        -----------------------------------------------------------------------
            Palindromic,            if the site is palindromic or not.
            NotPalindromic          allow some optimisations of the code.
                                    No need to check the reverse strand
                                    with palindromic sites.
        -----------------------------------------------------------------------
            Unknown, Blunt,         represent the overhang.
            Ov5, Ov3                Unknown is here for symmetry reasons and
                                    correspond to enzymes that are not
                                    characterised in rebase.
        -----------------------------------------------------------------------
            Defined, Ambiguous,     represent the sequence of the overhang.
            NotDefined
                                    NotDefined is for enzymes not characterised
                                    in rebase.

                                    Defined correspond to enzymes that display
                                    a constant overhang whatever the sequence.
                                    ex : EcoRI. G^AATTC -> overhang :AATT
                                                CTTAA^G

                                    Ambiguous : the overhang varies with the
                                    sequence restricted.
                                    Typically enzymes which cut outside their
                                    restriction site or (but not always)
                                    inside an ambiguous site.
                                    ex:
                                    AcuI CTGAAG(22/20)  -> overhang : NN
                                    AasI GACNNN^NNNGTC  -> overhang : NN
                                         CTGN^NNNNNCAG

                note : these 3 classes refers to the overhang not the site.
                   So the enzyme ApoI (RAATTY) is defined even if its
                   restriction site is ambiguous.

                        ApoI R^AATTY -> overhang : AATT -> Defined
                             YTTAA^R
                   Accordingly, blunt enzymes are always Defined even
                   when they cut outside their restriction site.
        -----------------------------------------------------------------------
            Not_available,          as found in rebase file emboss_r.NNN files.
            Commercially_available
                                    allow the selection of the enzymes
                                    according to their suppliers to reduce the
                                    quantity of results.
                                    Also will allow the implementation of
                                    buffer compatibility tables. Not
                                    implemented yet.

                                    the list of suppliers is extracted from
                                    emboss_s.NNN
        -----------------------------------------------------------------------

"""


import warnings

import re
import string
import itertools

from Bio.Seq import Seq, MutableSeq
from Bio.Restriction.Restriction_Dictionary import rest_dict as enzymedict
from Bio.Restriction.Restriction_Dictionary import typedict
from Bio.Restriction.Restriction_Dictionary import suppliers as suppliers_dict
from Bio.Restriction.PrintFormat import PrintFormat
from Bio import BiopythonWarning


matching = {
    "A": "ARWMHVDN",
    "C": "CYSMHBVN",
    "G": "GRSKBVDN",
    "T": "TYWKHBDN",
    "R": "ABDGHKMNSRWV",
    "Y": "CBDHKMNSTWVY",
    "W": "ABDHKMNRTWVY",
    "S": "CBDGHKMNSRVY",
    "M": "ACBDHMNSRWVY",
    "K": "BDGHKNSRTWVY",
    "H": "ACBDHKMNSRTWVY",
    "B": "CBDGHKMNSRTWVY",
    "V": "ACBDGHKMNSRWVY",
    "D": "ABDGHKMNSRTWVY",
    "N": "ACBDGHKMNSRTWVY",
}

DNA = Seq


class FormattedSeq:
    """A linear or circular sequence object for restriction analysis.

    Translates a Bio.Seq into a formatted sequence to be used with Restriction.

    Roughly: remove anything which is not IUPAC alphabet and then add a space
             in front of the sequence to get a biological index instead of a
             python index (i.e. index of the first base is 1 not 0).

    Retains information about the shape of the molecule linear (default) or
    circular. Restriction sites are search over the edges of circular sequence.
    """

    _remove_chars = string.whitespace.encode() + string.digits.encode()
    _table = bytearray(256)
    upper_to_lower = ord("A") - ord("a")
    for c in b"ABCDGHKMNRSTVWY":  # Only allow IUPAC letters
        _table[c] = c  # map uppercase to uppercase
        _table[c - upper_to_lower] = c  # map lowercase to uppercase
    del upper_to_lower
    _table = bytes(_table)

    def __init__(self, seq, linear=True):
        """Initialize ``FormattedSeq`` with sequence and topology (optional).

        ``seq`` is either a ``Bio.Seq``, ``Bio.MutableSeq`` or a
        ``FormattedSeq``. If ``seq`` is a ``FormattedSeq``, ``linear``
        will have no effect on the shape of the sequence.
        """
        if isinstance(seq, (Seq, MutableSeq)):
            self.lower = seq.islower()
            data = bytes(seq)
            self.data = data.translate(self._table, delete=self._remove_chars)
            if 0 in self.data:  # Check if all letters were IUPAC
                raise TypeError(f"Invalid character found in {data.decode()}")
            # Note this adds a leading space to the sequence (!)
            self.data = " " + self.data.decode("ASCII")
            self.linear = linear
            self.klass = seq.__class__
        elif isinstance(seq, FormattedSeq):
            self.lower = seq.lower
            self.data = seq.data
            self.linear = seq.linear
            self.klass = seq.klass
        else:
            raise TypeError(f"expected Seq or MutableSeq, got {type(seq)}")

    def __len__(self):
        """Return length of ``FormattedSeq``.

        ``FormattedSeq`` has a leading space, thus subtract 1.
        """
        return len(self.data) - 1

    def __repr__(self):
        """Represent ``FormattedSeq`` class as a string."""
        return f"FormattedSeq({self[1:]!r}, linear={self.linear!r})"

    def __eq__(self, other):
        """Implement equality operator for ``FormattedSeq`` object."""
        if isinstance(other, FormattedSeq):
            if repr(self) == repr(other):
                return True
            else:
                return False
        return False

    def circularise(self):
        """Circularise sequence in place."""
        self.linear = False

    def linearise(self):
        """Linearise sequence in place."""
        self.linear = True

    def to_linear(self):
        """Make a new instance of sequence as linear."""
        new = self.__class__(self)
        new.linear = True
        return new

    def to_circular(self):
        """Make a new instance of sequence as circular."""
        new = self.__class__(self)
        new.linear = False
        return new

    def is_linear(self):
        """Return if sequence is linear (True) or circular (False)."""
        return self.linear

    def finditer(self, pattern, size):
        """Return a list of a given pattern which occurs in the sequence.

        The list is made of tuple (location, pattern.group).
        The latter is used with non palindromic sites.
        Pattern is the regular expression pattern corresponding to the
        enzyme restriction site.
        Size is the size of the restriction enzyme recognition-site size.
        """
        if self.is_linear():
            data = self.data
        else:
            data = self.data + self.data[1:size]
        return [(i.start(), i.group) for i in re.finditer(pattern, data)]

    def __getitem__(self, i):
        """Return substring of ``FormattedSeq``.

        The class of the returned object is the class of the respective
        sequence. Note that due to the leading space, indexing is 1-based:

        >>> from Bio.Seq import Seq
        >>> from Bio.Restriction.Restriction import FormattedSeq
        >>> f_seq = FormattedSeq(Seq('ATGCATGC'))
        >>> f_seq[1]
        Seq('A')

        """
        if self.lower:
            return self.klass(self.data[i].lower())
        return self.klass(self.data[i])


class RestrictionType(type):
    """RestrictionType. Type from which all enzyme classes are derived.

    Implement the operator methods.
    """

    def __init__(cls, name="", bases=(), dct=None):
        """Initialize RestrictionType instance.

        Not intended to be used in normal operation. The enzymes are
        instantiated when importing the module.
        See below.
        """
        if "-" in name:
            raise ValueError(f"Problem with hyphen in {name!r} as enzyme name")
        # 2011/11/26 - Nobody knows what this call was supposed to accomplish,
        # but all unit tests seem to pass without it.
        # super().__init__(cls, name, bases, dct)
        try:
            cls.compsite = re.compile(cls.compsite)
        except AttributeError:
            # Can happen if initialised wrongly.
            # (This was seen when Sphinx api-doc imports the classes, and
            # tried to automatically general documentation for them)
            pass
        except Exception:
            raise ValueError(
                f"Problem with regular expression, re.compiled({cls.compsite!r})"
            ) from None

    def __add__(cls, other):
        """Add restriction enzyme to a RestrictionBatch().

        If other is an enzyme returns a batch of the two enzymes.
        If other is already a RestrictionBatch add enzyme to it.
        """
        if isinstance(other, RestrictionType):
            return RestrictionBatch([cls, other])
        elif isinstance(other, RestrictionBatch):
            return other.add_nocheck(cls)
        else:
            raise TypeError

    def __truediv__(cls, other):
        """Override '/' operator to use as search method.

        >>> from Bio.Restriction import EcoRI
        >>> EcoRI/Seq('GAATTC')
        [2]

        Returns RE.search(other).
        """
        return cls.search(other)

    def __rtruediv__(cls, other):
        """Override division with reversed operands to use as search method.

        >>> from Bio.Restriction import EcoRI
        >>> Seq('GAATTC')/EcoRI
        [2]

        Returns RE.search(other).
        """
        return cls.search(other)

    def __floordiv__(cls, other):
        """Override '//' operator to use as catalyse method.

        >>> from Bio.Restriction import EcoRI
        >>> EcoRI//Seq('GAATTC')
        (Seq('G'), Seq('AATTC'))

        Returns RE.catalyse(other).
        """
        return cls.catalyse(other)

    def __rfloordiv__(cls, other):
        """As __floordiv__, with reversed operands.

        >>> from Bio.Restriction import EcoRI
        >>> Seq('GAATTC')//EcoRI
        (Seq('G'), Seq('AATTC'))

        Returns RE.catalyse(other).
        """
        return cls.catalyse(other)

    def __str__(cls):
        """Return the name of the enzyme as string."""
        return cls.__name__

    def __repr__(cls):
        """Implement repr method.

        Used with eval or exec will instantiate the enzyme.
        """
        return f"{cls.__name__}"

    def __len__(cls):
        """Return length of recognition site of enzyme as int."""
        try:
            return cls.size
        except AttributeError:
            # Happens if the instance was not initialised as expected.
            # e.g. if instance created by a documentation framework
            # like Sphinx trying to inspect the class automatically,
            # Also seen within IPython.
            return 0

    def __hash__(cls):
        """Implement ``hash()`` method for ``RestrictionType``.

        Python default is to use ``id(...)``
        This is consistent with the ``__eq__`` implementation
        """
        return id(cls)

    def __eq__(cls, other):
        """Override '==' operator.

        True if RE and other are the same enzyme.

        Specifically this checks they are the same Python object.
        """
        # assert (id(cls)==id(other)) == (other is cls) == (cls is other)
        return id(cls) == id(other)

    def __ne__(cls, other):
        """Override '!=' operator.

        Isoschizomer strict (same recognition site, same restriction) -> False
        All the other-> True

        WARNING - This is not the inverse of the __eq__ method

        >>> from Bio.Restriction import SacI, SstI
        >>> SacI != SstI  # true isoschizomers
        False
        >>> SacI == SstI
        False
        """
        if not isinstance(other, RestrictionType):
            return True
        elif cls.charac == other.charac:
            return False
        else:
            return True

    def __rshift__(cls, other):
        """Override '>>' operator to test for neoschizomers.

        neoschizomer : same recognition site, different restriction. -> True
        all the others :                                             -> False

        >>> from Bio.Restriction import SmaI, XmaI
        >>> SmaI >> XmaI
        True
        """
        if not isinstance(other, RestrictionType):
            return False
        elif cls.site == other.site and cls.charac != other.charac:
            return True
        else:
            return False

    def __mod__(cls, other):
        """Override '%' operator to test for compatible overhangs.

        True if a and b have compatible overhang.

        >>> from Bio.Restriction import XhoI, SalI
        >>> XhoI % SalI
        True
        """
        if not isinstance(other, RestrictionType):
            raise TypeError(f"expected RestrictionType, got {type(other)} instead")
        return cls._mod1(other)

    def __ge__(cls, other):
        """Compare length of recognition site of two enzymes.

        Override '>='. a is greater or equal than b if the a site is longer
        than b site. If their site have the same length sort by alphabetical
        order of their names.

        >>> from Bio.Restriction import EcoRI, EcoRV
        >>> EcoRI.size
        6
        >>> EcoRV.size
        6
        >>> EcoRI >= EcoRV
        False
        """
        if not isinstance(other, RestrictionType):
            raise NotImplementedError
        if len(cls) > len(other):
            return True
        elif cls.size == len(other) and cls.__name__ >= other.__name__:
            return True
        else:
            return False

    def __gt__(cls, other):
        """Compare length of recognition site of two enzymes.

        Override '>'. Sorting order:

        1. size of the recognition site.
        2. if equal size, alphabetical order of the names.

        """
        if not isinstance(other, RestrictionType):
            raise NotImplementedError
        if len(cls) > len(other):
            return True
        elif cls.size == len(other) and cls.__name__ > other.__name__:
            return True
        else:
            return False

    def __le__(cls, other):
        """Compare length of recognition site of two enzymes.

        Override '<='. Sorting order:

        1. size of the recognition site.
        2. if equal size, alphabetical order of the names.

        """
        if not isinstance(other, RestrictionType):
            raise NotImplementedError
        elif len(cls) < len(other):
            return True
        elif len(cls) == len(other) and cls.__name__ <= other.__name__:
            return True
        else:
            return False

    def __lt__(cls, other):
        """Compare length of recognition site of two enzymes.

        Override '<'. Sorting order:

        1. size of the recognition site.
        2. if equal size, alphabetical order of the names.

        """
        if not isinstance(other, RestrictionType):
            raise NotImplementedError
        elif len(cls) < len(other):
            return True
        elif len(cls) == len(other) and cls.__name__ < other.__name__:
            return True
        else:
            return False


class AbstractCut(RestrictionType):
    """Implement the methods that are common to all restriction enzymes.

    All the methods are classmethod.

    For internal use only. Not meant to be instantiated.
    """

    @classmethod
    def search(cls, dna, linear=True):
        """Return a list of cutting sites of the enzyme in the sequence.

        Compensate for circular sequences and so on.

        dna must be a Bio.Seq.Seq instance or a Bio.Seq.MutableSeq instance.

        If linear is False, the restriction sites that span over the boundaries
        will be included.

        The positions are the first base of the 3' fragment,
        i.e. the first base after the position the enzyme will cut.
        """
        #
        #   Separating search from _search allow a (very limited) optimisation
        #   of the search when using a batch of restriction enzymes.
        #   in this case the DNA is tested once by the class which implements
        #   the batch instead of being tested by each enzyme single.
        #   see RestrictionBatch.search() for example.
        #
        if isinstance(dna, FormattedSeq):
            cls.dna = dna
            return cls._search()
        else:
            cls.dna = FormattedSeq(dna, linear)
            return cls._search()

    @classmethod
    def all_suppliers(cls):
        """Print all the suppliers of restriction enzyme."""
        supply = sorted(x[0] for x in suppliers_dict.values())
        print(",\n".join(supply))

    @classmethod
    def is_equischizomer(cls, other):
        """Test for real isoschizomer.

        True if other is an isoschizomer of RE, but not an neoschizomer,
        else False.

        Equischizomer: same site, same position of restriction.

        >>> from Bio.Restriction import SacI, SstI, SmaI, XmaI
        >>> SacI.is_equischizomer(SstI)
        True
        >>> SmaI.is_equischizomer(XmaI)
        False

        """
        return not cls != other

    @classmethod
    def is_neoschizomer(cls, other):
        """Test for neoschizomer.

        True if other is an isoschizomer of RE, else False.
        Neoschizomer: same site, different position of restriction.
        """
        return cls >> other

    @classmethod
    def is_isoschizomer(cls, other):
        """Test for same recognition site.

        True if other has the same recognition site, else False.

        Isoschizomer: same site.

        >>> from Bio.Restriction import SacI, SstI, SmaI, XmaI
        >>> SacI.is_isoschizomer(SstI)
        True
        >>> SmaI.is_isoschizomer(XmaI)
        True

        """
        return (not cls != other) or cls >> other

    @classmethod
    def equischizomers(cls, batch=None):
        """List equischizomers of the enzyme.

        Return a tuple of all the isoschizomers of RE.
        If batch is supplied it is used instead of the default AllEnzymes.

        Equischizomer: same site, same position of restriction.
        """
        if not batch:
            batch = AllEnzymes
        r = [x for x in batch if not cls != x]
        i = r.index(cls)
        del r[i]
        r.sort()
        return r

    @classmethod
    def neoschizomers(cls, batch=None):
        """List neoschizomers of the enzyme.

        Return a tuple of all the neoschizomers of RE.
        If batch is supplied it is used instead of the default AllEnzymes.

        Neoschizomer: same site, different position of restriction.
        """
        if not batch:
            batch = AllEnzymes
        r = sorted(x for x in batch if cls >> x)
        return r

    @classmethod
    def isoschizomers(cls, batch=None):
        """List all isoschizomers of the enzyme.

        Return a tuple of all the equischizomers and neoschizomers of RE.
        If batch is supplied it is used instead of the default AllEnzymes.
        """
        if not batch:
            batch = AllEnzymes
        r = [x for x in batch if (cls >> x) or (not cls != x)]
        i = r.index(cls)
        del r[i]
        r.sort()
        return r

    @classmethod
    def frequency(cls):
        """Return the theoretically cutting frequency of the enzyme.

        Frequency of the site, given as 'one cut per x bases' (int).
        """
        return cls.freq


class NoCut(AbstractCut):
    """Implement the methods specific to the enzymes that do not cut.

    These enzymes are generally enzymes that have been only partially
    characterised and the way they cut the DNA is unknown or enzymes for
    which the pattern of cut is to complex to be recorded in Rebase
    (ncuts values of 0 in emboss_e.###).

    When using search() with these enzymes the values returned are at the start
    of the restriction site.

    Their catalyse() method returns a TypeError.

    Unknown and NotDefined are also part of the base classes of these enzymes.

    Internal use only. Not meant to be instantiated.
    """

    @classmethod
    def cut_once(cls):
        """Return if the cutting pattern has one cut.

        True if the enzyme cut the sequence one time on each strand.
        """
        return False

    @classmethod
    def cut_twice(cls):
        """Return if the cutting pattern has two cuts.

        True if the enzyme cut the sequence twice on each strand.
        """
        return False

    @classmethod
    def _modify(cls, location):
        """Return a generator that moves the cutting position by 1 (PRIVATE).

        For internal use only.

        location is an integer corresponding to the location of the match for
        the enzyme pattern in the sequence.
        _modify returns the real place where the enzyme will cut.

        Example::

            EcoRI pattern : GAATTC
            EcoRI will cut after the G.
            so in the sequence:
                     ______
            GAATACACGGAATTCGA
                     |
                     10
            dna.finditer(GAATTC, 6) will return 10 as G is the 10th base
            EcoRI cut after the G so:
            EcoRI._modify(10) -> 11.

        If the enzyme cut twice _modify will returns two integer corresponding
        to each cutting site.
        """
        yield location

    @classmethod
    def _rev_modify(cls, location):
        """Return a generator that moves the cutting position by 1 (PRIVATE).

        For internal use only.

        As _modify for site situated on the antiparallel strand when the
        enzyme is not palindromic.
        """
        yield location

    @classmethod
    def characteristic(cls):
        """Return a list of the enzyme's characteristics as tuple.

        the tuple contains the attributes:

        - fst5 -> first 5' cut ((current strand) or None
        - fst3 -> first 3' cut (complementary strand) or None
        - scd5 -> second 5' cut (current strand) or None
        - scd5 -> second 3' cut (complementary strand) or None
        - site -> recognition site.

        """
        return None, None, None, None, cls.site


class OneCut(AbstractCut):
    """Implement the methods for enzymes that cut the DNA only once.

    Correspond to ncuts values of 2 in emboss_e.###

    Internal use only. Not meant to be instantiated.
    """

    @classmethod
    def cut_once(cls):
        """Return if the cutting pattern has one cut.

        True if the enzyme cut the sequence one time on each strand.
        """
        return True

    @classmethod
    def cut_twice(cls):
        """Return if the cutting pattern has two cuts.

        True if the enzyme cut the sequence twice on each strand.
        """
        return False

    @classmethod
    def _modify(cls, location):
        """Return a generator that moves the cutting position by 1 (PRIVATE).

        For internal use only.

        location is an integer corresponding to the location of the match for
        the enzyme pattern in the sequence.
        _modify returns the real place where the enzyme will cut.

        Example::

            EcoRI pattern : GAATTC
            EcoRI will cut after the G.
            so in the sequence:
                     ______
            GAATACACGGAATTCGA
                     |
                     10
            dna.finditer(GAATTC, 6) will return 10 as G is the 10th base
            EcoRI cut after the G so:
            EcoRI._modify(10) -> 11.

        if the enzyme cut twice _modify will returns two integer corresponding
        to each cutting site.
        """
        yield location + cls.fst5

    @classmethod
    def _rev_modify(cls, location):
        """Return a generator that moves the cutting position by 1 (PRIVATE).

        For internal use only.

        As _modify for site situated on the antiparallel strand when the
        enzyme is not palindromic
        """
        yield location - cls.fst3

    @classmethod
    def characteristic(cls):
        """Return a list of the enzyme's characteristics as tuple.

        The tuple contains the attributes:

        - fst5 -> first 5' cut ((current strand) or None
        - fst3 -> first 3' cut (complementary strand) or None
        - scd5 -> second 5' cut (current strand) or None
        - scd5 -> second 3' cut (complementary strand) or None
        - site -> recognition site.

        """
        return cls.fst5, cls.fst3, None, None, cls.site


class TwoCuts(AbstractCut):
    """Implement the methods for enzymes that cut the DNA twice.

    Correspond to ncuts values of 4 in emboss_e.###

    Internal use only. Not meant to be instantiated.
    """

    @classmethod
    def cut_once(cls):
        """Return if the cutting pattern has one cut.

        True if the enzyme cut the sequence one time on each strand.
        """
        return False

    @classmethod
    def cut_twice(cls):
        """Return if the cutting pattern has two cuts.

        True if the enzyme cut the sequence twice on each strand.
        """
        return True

    @classmethod
    def _modify(cls, location):
        """Return a generator that moves the cutting position by 1 (PRIVATE).

        For internal use only.

        location is an integer corresponding to the location of the match for
        the enzyme pattern in the sequence.
        _modify returns the real place where the enzyme will cut.

        example::

            EcoRI pattern : GAATTC
            EcoRI will cut after the G.
            so in the sequence:
                     ______
            GAATACACGGAATTCGA
                     |
                     10
            dna.finditer(GAATTC, 6) will return 10 as G is the 10th base
            EcoRI cut after the G so:
            EcoRI._modify(10) -> 11.

        if the enzyme cut twice _modify will returns two integer corresponding
        to each cutting site.
        """
        yield location + cls.fst5
        yield location + cls.scd5

    @classmethod
    def _rev_modify(cls, location):
        """Return a generator that moves the cutting position by 1 (PRIVATE).

        for internal use only.

        as _modify for site situated on the antiparallel strand when the
        enzyme is not palindromic
        """
        yield location - cls.fst3
        yield location - cls.scd3

    @classmethod
    def characteristic(cls):
        """Return a list of the enzyme's characteristics as tuple.

        the tuple contains the attributes:

        - fst5 -> first 5' cut ((current strand) or None
        - fst3 -> first 3' cut (complementary strand) or None
        - scd5 -> second 5' cut (current strand) or None
        - scd5 -> second 3' cut (complementary strand) or None
        - site -> recognition site.

        """
        return cls.fst5, cls.fst3, cls.scd5, cls.scd3, cls.site


class Meth_Dep(AbstractCut):
    """Implement the information about methylation.

    Enzymes of this class possess a site which is methylable.
    """

    @classmethod
    def is_methylable(cls):
        """Return if recognition site can be methylated.

        True if the recognition site is a methylable.
        """
        return True


class Meth_Undep(AbstractCut):
    """Implement information about methylation sensitibility.

    Enzymes of this class are not sensible to methylation.
    """

    @classmethod
    def is_methylable(cls):
        """Return if recognition site can be methylated.

        True if the recognition site is a methylable.
        """
        return False


class Palindromic(AbstractCut):
    """Implement methods for enzymes with palindromic recognition sites.

    palindromic means : the recognition site and its reverse complement are
                        identical.
    Remarks     : an enzyme with a site CGNNCG is palindromic even if some
                  of the sites that it will recognise are not.
                  for example here : CGAACG

    Internal use only. Not meant to be instantiated.
    """

    @classmethod
    def _search(cls):
        """Return a list of cutting sites of the enzyme in the sequence (PRIVATE).

        For internal use only.

        Implement the search method for palindromic enzymes.
        """
        siteloc = cls.dna.finditer(cls.compsite, cls.size)
        cls.results = [r for s, g in siteloc for r in cls._modify(s)]
        if cls.results:
            cls._drop()
        return cls.results

    @classmethod
    def is_palindromic(cls):
        """Return if the enzyme has a palindromic recoginition site."""
        return True


class NonPalindromic(AbstractCut):
    """Implement methods for enzymes with non-palindromic recognition sites.

    Palindromic means : the recognition site and its reverse complement are
                        identical.

    Internal use only. Not meant to be instantiated.
    """

    @classmethod
    def _search(cls):
        """Return a list of cutting sites of the enzyme in the sequence (PRIVATE).

        For internal use only.

        Implement the search method for non palindromic enzymes.
        """
        iterator = cls.dna.finditer(cls.compsite, cls.size)
        cls.results = []
        modif = cls._modify
        revmodif = cls._rev_modify
        s = str(cls)
        cls.on_minus = []

        for start, group in iterator:
            if group(s):
                cls.results += list(modif(start))
            else:
                cls.on_minus += list(revmodif(start))
        cls.results += cls.on_minus

        if cls.results:
            cls.results.sort()
            cls._drop()
        return cls.results

    @classmethod
    def is_palindromic(cls):
        """Return if the enzyme has a palindromic recoginition site."""
        return False


class Unknown(AbstractCut):
    """Implement methods for enzymes that produce unknown overhangs.

    These enzymes are also NotDefined and NoCut.

    Internal use only. Not meant to be instantiated.
    """

    @classmethod
    def catalyse(cls, dna, linear=True):
        """List the sequence fragments after cutting dna with enzyme.

        Return a tuple of dna as will be produced by using RE to restrict the
        dna.

        dna must be a Bio.Seq.Seq instance or a Bio.Seq.MutableSeq instance.

        If linear is False, the sequence is considered to be circular and the
        output will be modified accordingly.
        """
        raise NotImplementedError(f"{cls.__name__} restriction is unknown.")

    catalyze = catalyse

    @classmethod
    def is_blunt(cls):
        """Return if the enzyme produces blunt ends.

        True if the enzyme produces blunt end.

        Related methods:

        - RE.is_3overhang()
        - RE.is_5overhang()
        - RE.is_unknown()

        """
        return False

    @classmethod
    def is_5overhang(cls):
        """Return if the enzymes produces 5' overhanging ends.

        True if the enzyme produces 5' overhang sticky end.

        Related methods:

        - RE.is_3overhang()
        - RE.is_blunt()
        - RE.is_unknown()

        """
        return False

    @classmethod
    def is_3overhang(cls):
        """Return if the enzyme produces 3' overhanging ends.

        True if the enzyme produces 3' overhang sticky end.

        Related methods:

        - RE.is_5overhang()
        - RE.is_blunt()
        - RE.is_unknown()

        """
        return False

    @classmethod
    def overhang(cls):
        """Return the type of the enzyme's overhang as string.

        Can be "3' overhang", "5' overhang", "blunt", "unknown".
        """
        return "unknown"

    @classmethod
    def compatible_end(cls):
        """List all enzymes that produce compatible ends for the enzyme."""
        return []

    @classmethod
    def _mod1(cls, other):
        """Test if other enzyme produces compatible ends for enzyme (PRIVATE).

        For internal use only.

        Test for the compatibility of restriction ending of RE and other.
        """
        return False


class Blunt(AbstractCut):
    """Implement methods for enzymes that produce blunt ends.

    The enzyme cuts the + strand and the - strand of the DNA at the same
    place.

    Internal use only. Not meant to be instantiated.
    """

    @classmethod
    def catalyse(cls, dna, linear=True):
        """List the sequence fragments after cutting dna with enzyme.

        Return a tuple of dna as will be produced by using RE to restrict the
        dna.

        dna must be a Bio.Seq.Seq instance or a Bio.Seq.MutableSeq instance.

        If linear is False, the sequence is considered to be circular and the
        output will be modified accordingly.
        """
        r = cls.search(dna, linear)
        d = cls.dna
        if not r:
            return (d[1:],)
        fragments = []
        length = len(r) - 1
        if d.is_linear():
            #
            #   START of the sequence to FIRST site.
            #
            fragments.append(d[1 : r[0]])
            if length:
                #
                #   if more than one site add them.
                #
                fragments += [d[r[x] : r[x + 1]] for x in range(length)]
            #
            #   LAST site to END of the sequence.
            #
            fragments.append(d[r[-1] :])
        else:
            #
            #   circular : bridge LAST site to FIRST site.
            #
            fragments.append(d[r[-1] :] + d[1 : r[0]])
            if not length:
                #
                #   one site we finish here.
                #
                return tuple(fragments)
            #
            #   add the others.
            #
            fragments += [d[r[x] : r[x + 1]] for x in range(length)]
        return tuple(fragments)

    catalyze = catalyse

    @classmethod
    def is_blunt(cls):
        """Return if the enzyme produces blunt ends.

        True if the enzyme produces blunt end.

        Related methods:

        - RE.is_3overhang()
        - RE.is_5overhang()
        - RE.is_unknown()

        """
        return True

    @classmethod
    def is_5overhang(cls):
        """Return if the enzymes produces 5' overhanging ends.

        True if the enzyme produces 5' overhang sticky end.

        Related methods:

        - RE.is_3overhang()
        - RE.is_blunt()
        - RE.is_unknown()

        """
        return False

    @classmethod
    def is_3overhang(cls):
        """Return if the enzyme produces 3' overhanging ends.

        True if the enzyme produces 3' overhang sticky end.

        Related methods:

        - RE.is_5overhang()
        - RE.is_blunt()
        - RE.is_unknown()

        """
        return False

    @classmethod
    def overhang(cls):
        """Return the type of the enzyme's overhang as string.

        Can be "3' overhang", "5' overhang", "blunt", "unknown".
        """
        return "blunt"

    @classmethod
    def compatible_end(cls, batch=None):
        """List all enzymes that produce compatible ends for the enzyme."""
        if not batch:
            batch = AllEnzymes
        r = sorted(x for x in iter(AllEnzymes) if x.is_blunt())
        return r

    @staticmethod
    def _mod1(other):
        """Test if other enzyme produces compatible ends for enzyme (PRIVATE).

        For internal use only

        Test for the compatibility of restriction ending of RE and other.
        """
        return issubclass(other, Blunt)


class Ov5(AbstractCut):
    """Implement methods for enzymes that produce 5' overhanging ends.

    The enzyme cuts the + strand after the - strand of the DNA.

    Internal use only. Not meant to be instantiated.
    """

    @classmethod
    def catalyse(cls, dna, linear=True):
        """List the sequence fragments after cutting dna with enzyme.

        Return a tuple of dna as will be produced by using RE to restrict the
        dna.

        dna must be a Bio.Seq.Seq instance or a Bio.Seq.MutableSeq instance.

        If linear is False, the sequence is considered to be circular and the
        output will be modified accordingly.
        """
        r = cls.search(dna, linear)
        d = cls.dna
        if not r:
            return (d[1:],)
        length = len(r) - 1
        fragments = []
        if d.is_linear():
            #
            #   START of the sequence to FIRST site.
            #
            fragments.append(d[1 : r[0]])
            if length:
                #
                #   if more than one site add them.
                #
                fragments += [d[r[x] : r[x + 1]] for x in range(length)]
            #
            #   LAST site to END of the sequence.
            #
            fragments.append(d[r[-1] :])
        else:
            #
            #   circular : bridge LAST site to FIRST site.
            #
            fragments.append(d[r[-1] :] + d[1 : r[0]])
            if not length:
                #
                #   one site we finish here.
                #
                return tuple(fragments)
            #
            #   add the others.
            #
            fragments += [d[r[x] : r[x + 1]] for x in range(length)]
        return tuple(fragments)

    catalyze = catalyse

    @classmethod
    def is_blunt(cls):
        """Return if the enzyme produces blunt ends.

        True if the enzyme produces blunt end.

        Related methods:

        - RE.is_3overhang()
        - RE.is_5overhang()
        - RE.is_unknown()

        """
        return False

    @classmethod
    def is_5overhang(cls):
        """Return if the enzymes produces 5' overhanging ends.

        True if the enzyme produces 5' overhang sticky end.

        Related methods:

        - RE.is_3overhang()
        - RE.is_blunt()
        - RE.is_unknown()

        """
        return True

    @classmethod
    def is_3overhang(cls):
        """Return if the enzyme produces 3' overhanging ends.

        True if the enzyme produces 3' overhang sticky end.

        Related methods:

        - RE.is_5overhang()
        - RE.is_blunt()
        - RE.is_unknown()

        """
        return False

    @classmethod
    def overhang(cls):
        """Return the type of the enzyme's overhang as string.

        Can be "3' overhang", "5' overhang", "blunt", "unknown".
        """
        return "5' overhang"

    @classmethod
    def compatible_end(cls, batch=None):
        """List all enzymes that produce compatible ends for the enzyme."""
        if not batch:
            batch = AllEnzymes
        r = sorted(x for x in iter(AllEnzymes) if x.is_5overhang() and x % cls)
        return r

    @classmethod
    def _mod1(cls, other):
        """Test if other enzyme produces compatible ends for enzyme (PRIVATE).

        For internal use only.

        Test for the compatibility of restriction ending of RE and other.
        """
        if issubclass(other, Ov5):
            return cls._mod2(other)
        else:
            return False


class Ov3(AbstractCut):
    """Implement methods for enzymes that produce 3' overhanging ends.

    The enzyme cuts the - strand after the + strand of the DNA.

    Internal use only. Not meant to be instantiated.
    """

    @classmethod
    def catalyse(cls, dna, linear=True):
        """List the sequence fragments after cutting dna with enzyme.

        Return a tuple of dna as will be produced by using RE to restrict the
        dna.

        dna must be a Bio.Seq.Seq instance or a Bio.Seq.MutableSeq instance.

        If linear is False, the sequence is considered to be circular and the
        output will be modified accordingly.
        """
        r = cls.search(dna, linear)
        d = cls.dna
        if not r:
            return (d[1:],)
        fragments = []
        length = len(r) - 1
        if d.is_linear():
            #
            #   START of the sequence to FIRST site.
            #
            fragments.append(d[1 : r[0]])
            if length:
                #
                #   if more than one site add them.
                #
                fragments += [d[r[x] : r[x + 1]] for x in range(length)]
            #
            #   LAST site to END of the sequence.
            #
            fragments.append(d[r[-1] :])
        else:
            #
            #   circular : bridge LAST site to FIRST site.
            #
            fragments.append(d[r[-1] :] + d[1 : r[0]])
            if not length:
                #
                #   one site we finish here.
                #
                return tuple(fragments)
            #
            #   add the others.
            #
            fragments += [d[r[x] : r[x + 1]] for x in range(length)]
        return tuple(fragments)

    catalyze = catalyse

    @classmethod
    def is_blunt(cls):
        """Return if the enzyme produces blunt ends.

        True if the enzyme produces blunt end.

        Related methods:

        - RE.is_3overhang()
        - RE.is_5overhang()
        - RE.is_unknown()

        """
        return False

    @classmethod
    def is_5overhang(cls):
        """Return if the enzymes produces 5' overhanging ends.

        True if the enzyme produces 5' overhang sticky end.

        Related methods:

        - RE.is_3overhang()
        - RE.is_blunt()
        - RE.is_unknown()

        """
        return False

    @classmethod
    def is_3overhang(cls):
        """Return if the enzyme produces 3' overhanging ends.

        True if the enzyme produces 3' overhang sticky end.

        Related methods:

        - RE.is_5overhang()
        - RE.is_blunt()
        - RE.is_unknown()

        """
        return True

    @classmethod
    def overhang(cls):
        """Return the type of the enzyme's overhang as string.

        Can be "3' overhang", "5' overhang", "blunt", "unknown".
        """
        return "3' overhang"

    @classmethod
    def compatible_end(cls, batch=None):
        """List all enzymes that produce compatible ends for the enzyme."""
        if not batch:
            batch = AllEnzymes
        r = sorted(x for x in iter(AllEnzymes) if x.is_3overhang() and x % cls)
        return r

    @classmethod
    def _mod1(cls, other):
        """Test if other enzyme produces compatible ends for enzyme (PRIVATE).

        For internal use only.

        Test for the compatibility of restriction ending of RE and other.
        """
        #
        #   called by RE._mod1(other) when the one of the enzyme is ambiguous
        #
        if issubclass(other, Ov3):
            return cls._mod2(other)
        else:
            return False


class Defined(AbstractCut):
    """Implement methods for enzymes with defined recognition site and cut.

    Typical example : EcoRI -> G^AATT_C
                      The overhang will always be AATT
    Notes:
        Blunt enzymes are always defined. Even if their site is GGATCCNNN^_N
        Their overhang is always the same : blunt!

    Internal use only. Not meant to be instantiated.
    """

    @classmethod
    def _drop(cls):
        """Remove cuts that are outsite of the sequence (PRIVATE).

        For internal use only.

        Drop the site that are situated outside the sequence in linear
        sequence. Modify the index for site in circular sequences.
        """
        #
        #   remove or modify the results that are outside the sequence.
        #   This is necessary since after finding the site we add the distance
        #   from the site to the cut with the _modify and _rev_modify methods.
        #   For linear we will remove these sites altogether.
        #   For circular sequence, we modify the result rather than _drop it
        #   since the site is in the sequence.
        #
        length = len(cls.dna)
        drop = itertools.dropwhile
        take = itertools.takewhile
        if cls.dna.is_linear():
            cls.results = list(drop(lambda x: x <= 1, cls.results))
            cls.results = list(take(lambda x: x <= length, cls.results))
        else:
            for index, location in enumerate(cls.results):
                if location < 1:
                    cls.results[index] += length
                else:
                    break
            for index, location in enumerate(cls.results[::-1]):
                if location > length:
                    cls.results[-(index + 1)] -= length
                else:
                    break

    @classmethod
    def is_defined(cls):
        """Return if recognition sequence and cut are defined.

        True if the sequence recognised and cut is constant,
        i.e. the recognition site is not degenerated AND the enzyme cut inside
        the site.

        Related methods:

        - RE.is_ambiguous()
        - RE.is_unknown()

        """
        return True

    @classmethod
    def is_ambiguous(cls):
        """Return if recognition sequence and cut may be ambiguous.

        True if the sequence recognised and cut is ambiguous,
        i.e. the recognition site is degenerated AND/OR the enzyme cut outside
        the site.

        Related methods:

        - RE.is_defined()
        - RE.is_unknown()

        """
        return False

    @classmethod
    def is_unknown(cls):
        """Return if recognition sequence is unknown.

        True if the sequence is unknown,
        i.e. the recognition site has not been characterised yet.

        Related methods:

        - RE.is_defined()
        - RE.is_ambiguous()

        """
        return False

    @classmethod
    def elucidate(cls):
        """Return a string representing the recognition site and cuttings.

        Return a representation of the site with the cut on the (+) strand
        represented as '^' and the cut on the (-) strand as '_'.
        ie:

        >>> from Bio.Restriction import EcoRI, KpnI, EcoRV, SnaI
        >>> EcoRI.elucidate()   # 5' overhang
        'G^AATT_C'
        >>> KpnI.elucidate()    # 3' overhang
        'G_GTAC^C'
        >>> EcoRV.elucidate()   # blunt
        'GAT^_ATC'
        >>> SnaI.elucidate()    # NotDefined, cut profile unknown.
        '? GTATAC ?'
        >>>

        """
        f5 = cls.fst5
        f3 = cls.fst3
        site = cls.site
        if cls.cut_twice():
            re = "cut twice, not yet implemented sorry."
        elif cls.is_5overhang():
            if f5 == f3 == 0:
                re = "N^" + cls.site + "_N"
            elif f3 == 0:
                re = site[:f5] + "^" + site[f5:] + "_N"
            else:
                re = site[:f5] + "^" + site[f5:f3] + "_" + site[f3:]
        elif cls.is_blunt():
            re = site[:f5] + "^_" + site[f5:]
        else:
            if f5 == f3 == 0:
                re = "N_" + site + "^N"
            else:
                re = site[:f3] + "_" + site[f3:f5] + "^" + site[f5:]
        return re

    @classmethod
    def _mod2(cls, other):
        """Test if other enzyme produces compatible ends for enzyme (PRIVATE).

        For internal use only.

        Test for the compatibility of restriction ending of RE and other.
        """
        #
        #   called by RE._mod1(other) when the one of the enzyme is ambiguous
        #
        if other.ovhgseq == cls.ovhgseq:
            return True
        elif issubclass(other, Ambiguous):
            return other._mod2(cls)
        else:
            return False


class Ambiguous(AbstractCut):
    """Implement methods for enzymes that produce variable overhangs.

    Typical example : BstXI -> CCAN_NNNN^NTGG
                      The overhang can be any sequence of 4 bases.

    Notes:
        Blunt enzymes are always defined. Even if their site is GGATCCNNN^_N
        Their overhang is always the same : blunt!

    Internal use only. Not meant to be instantiated.

    """

    @classmethod
    def _drop(cls):
        """Remove cuts that are outsite of the sequence (PRIVATE).

        For internal use only.

        Drop the site that are situated outside the sequence in linear
        sequence. Modify the index for site in circular sequences.
        """
        length = len(cls.dna)
        drop = itertools.dropwhile
        take = itertools.takewhile
        if cls.dna.is_linear():
            cls.results = list(drop(lambda x: x <= 1, cls.results))
            cls.results = list(take(lambda x: x <= length, cls.results))
        else:
            for index, location in enumerate(cls.results):
                if location < 1:
                    cls.results[index] += length
                else:
                    break
            for index, location in enumerate(cls.results[::-1]):
                if location > length:
                    cls.results[-(index + 1)] -= length
                else:
                    break

    @classmethod
    def is_defined(cls):
        """Return if recognition sequence and cut are defined.

        True if the sequence recognised and cut is constant,
        i.e. the recognition site is not degenerated AND the enzyme cut inside
        the site.

        Related methods:

        - RE.is_ambiguous()
        - RE.is_unknown()

        """
        return False

    @classmethod
    def is_ambiguous(cls):
        """Return if recognition sequence and cut may be ambiguous.

        True if the sequence recognised and cut is ambiguous,
        i.e. the recognition site is degenerated AND/OR the enzyme cut outside
        the site.

        Related methods:

        - RE.is_defined()
        - RE.is_unknown()

        """
        return True

    @classmethod
    def is_unknown(cls):
        """Return if recognition sequence is unknown.

        True if the sequence is unknown,
        i.e. the recognition site has not been characterised yet.

        Related methods:

        - RE.is_defined()
        - RE.is_ambiguous()

        """
        return False

    @classmethod
    def _mod2(cls, other):
        """Test if other enzyme produces compatible ends for enzyme (PRIVATE).

        For internal use only.

        Test for the compatibility of restriction ending of RE and other.
        """
        #
        #   called by RE._mod1(other) when the one of the enzyme is ambiguous
        #
        if len(cls.ovhgseq) != len(other.ovhgseq):
            return False
        else:
            se = cls.ovhgseq
            for base in se:
                if base in "ATCG":
                    pass
                if base in "N":
                    se = ".".join(se.split("N"))
                if base in "RYWMSKHDBV":
                    expand = "[" + matching[base] + "]"
                    se = expand.join(se.split(base))
            if re.match(se, other.ovhgseq):
                return True
            else:
                return False

    @classmethod
    def elucidate(cls):
        """Return a string representing the recognition site and cuttings.

        Return a representation of the site with the cut on the (+) strand
        represented as '^' and the cut on the (-) strand as '_'.
        ie:

        >>> from Bio.Restriction import EcoRI, KpnI, EcoRV, SnaI
        >>> EcoRI.elucidate()   # 5' overhang
        'G^AATT_C'
        >>> KpnI.elucidate()    # 3' overhang
        'G_GTAC^C'
        >>> EcoRV.elucidate()   # blunt
        'GAT^_ATC'
        >>> SnaI.elucidate()     # NotDefined, cut profile unknown.
        '? GTATAC ?'
        >>>

        """
        f5 = cls.fst5
        f3 = cls.fst3
        length = len(cls)
        site = cls.site
        if cls.cut_twice():
            re = "cut twice, not yet implemented sorry."
        elif cls.is_5overhang():
            if f3 == f5 == 0:
                re = "N^" + site + "_N"
            elif 0 <= f5 <= length and 0 <= f3 + length <= length:
                re = site[:f5] + "^" + site[f5:f3] + "_" + site[f3:]
            elif 0 <= f5 <= length:
                re = site[:f5] + "^" + site[f5:] + f3 * "N" + "_N"
            elif 0 <= f3 + length <= length:
                re = "N^" + abs(f5) * "N" + site[:f3] + "_" + site[f3:]
            elif f3 + length < 0:
                re = "N^" + abs(f5) * "N" + "_" + abs(length + f3) * "N" + site
            elif f5 > length:
                re = site + (f5 - length) * "N" + "^" + (length + f3 - f5) * "N" + "_N"
            else:
                re = "N^" + abs(f5) * "N" + site + f3 * "N" + "_N"
        elif cls.is_blunt():
            if f5 < 0:
                re = "N^_" + abs(f5) * "N" + site
            elif f5 > length:
                re = site + (f5 - length) * "N" + "^_N"
            else:
                raise ValueError("%s.easyrepr() : error f5=%i" % (cls.name, f5))
        else:
            if f3 == 0:
                if f5 == 0:
                    re = "N_" + site + "^N"
                else:
                    re = site + "_" + (f5 - length) * "N" + "^N"
            elif 0 < f3 + length <= length and 0 <= f5 <= length:
                re = site[:f3] + "_" + site[f3:f5] + "^" + site[f5:]
            elif 0 < f3 + length <= length:
                re = site[:f3] + "_" + site[f3:] + (f5 - length) * "N" + "^N"
            elif 0 <= f5 <= length:
                re = "N_" + "N" * (f3 + length) + site[:f5] + "^" + site[f5:]
            elif f3 > 0:
                re = site + f3 * "N" + "_" + (f5 - f3 - length) * "N" + "^N"
            elif f5 < 0:
                re = "N_" + abs(f3 - f5 + length) * "N" + "^" + abs(f5) * "N" + site
            else:
                re = "N_" + abs(f3 + length) * "N" + site + (f5 - length) * "N" + "^N"
        return re


class NotDefined(AbstractCut):
    """Implement methods for enzymes with non-characterized overhangs.

    Correspond to NoCut and Unknown.

    Internal use only. Not meant to be instantiated.
    """

    @classmethod
    def _drop(cls):
        """Remove cuts that are outsite of the sequence (PRIVATE).

        For internal use only.

        Drop the site that are situated outside the sequence in linear
        sequence. Modify the index for site in circular sequences.
        """
        if cls.dna.is_linear():
            return
        else:
            length = len(cls.dna)
            for index, location in enumerate(cls.results):
                if location < 1:
                    cls.results[index] += length
                else:
                    break
            for index, location in enumerate(cls.results[:-1]):
                if location > length:
                    cls.results[-(index + 1)] -= length
                else:
                    break

    @classmethod
    def is_defined(cls):
        """Return if recognition sequence and cut are defined.

        True if the sequence recognised and cut is constant,
        i.e. the recognition site is not degenerated AND the enzyme cut inside
        the site.

        Related methods:

        - RE.is_ambiguous()
        - RE.is_unknown()

        """
        return False

    @classmethod
    def is_ambiguous(cls):
        """Return if recognition sequence and cut may be ambiguous.

        True if the sequence recognised and cut is ambiguous,
        i.e. the recognition site is degenerated AND/OR the enzyme cut outside
        the site.

        Related methods:

        - RE.is_defined()
        - RE.is_unknown()

        """
        return False

    @classmethod
    def is_unknown(cls):
        """Return if recognition sequence is unknown.

        True if the sequence is unknown,
        i.e. the recognition site has not been characterised yet.

        Related methods:

        - RE.is_defined()
        - RE.is_ambiguous()

        """
        return True

    @classmethod
    def _mod2(cls, other):
        """Test if other enzyme produces compatible ends for enzyme (PRIVATE).

        For internal use only.

        Test for the compatibility of restriction ending of RE and other.
        """
        #
        #   Normally we should not arrive here. But well better safe than
        #   sorry.
        #   the overhang is not defined we are compatible with nobody.
        #   could raise an Error may be rather than return quietly.
        #
        # return False
        raise ValueError(
            "%s.mod2(%s), %s : NotDefined. pas glop pas glop!"
            % (str(cls), str(other), str(cls))
        )

    @classmethod
    def elucidate(cls):
        """Return a string representing the recognition site and cuttings.

        Return a representation of the site with the cut on the (+) strand
        represented as '^' and the cut on the (-) strand as '_'.
        ie:

        >>> from Bio.Restriction import EcoRI, KpnI, EcoRV, SnaI
        >>> EcoRI.elucidate()   # 5' overhang
        'G^AATT_C'
        >>> KpnI.elucidate()    # 3' overhang
        'G_GTAC^C'
        >>> EcoRV.elucidate()   # blunt
        'GAT^_ATC'
        >>> SnaI.elucidate()     # NotDefined, cut profile unknown.
        '? GTATAC ?'
        >>>

        """
        return f"? {cls.site} ?"


class Commercially_available(AbstractCut):
    """Implement methods for enzymes which are commercially available.

    Internal use only. Not meant to be instantiated.
    """

    #
    #   Recent addition to Rebase make this naming convention uncertain.
    #   May be better to says enzymes which have a supplier.
    #

    @classmethod
    def suppliers(cls):
        """Print a list of suppliers of the enzyme."""
        for s in cls.suppl:
            print(suppliers_dict[s][0] + ",")

    @classmethod
    def supplier_list(cls):
        """Return a list of suppliers of the enzyme."""
        return [v[0] for k, v in suppliers_dict.items() if k in cls.suppl]

    @classmethod
    def buffers(cls, supplier):
        """Return the recommended buffer of the supplier for this enzyme.

        Not implemented yet.
        """

    @classmethod
    def is_comm(cls):
        """Return if enzyme is commercially available.

        True if RE has suppliers.
        """
        return True


class Not_available(AbstractCut):
    """Implement methods for enzymes which are not commercially available.

    Internal use only. Not meant to be instantiated.
    """

    @staticmethod
    def suppliers():
        """Print a list of suppliers of the enzyme."""
        return None

    @classmethod
    def supplier_list(cls):
        """Return a list of suppliers of the enzyme."""
        return []

    @classmethod
    def buffers(cls, supplier):
        """Return the recommended buffer of the supplier for this enzyme.

        Not implemented yet.
        """
        raise TypeError("Enzyme not commercially available.")

    @classmethod
    def is_comm(cls):
        """Return if enzyme is commercially available.

        True if RE has suppliers.
        """
        return False


###############################################################################
#                                                                             #
#                       Restriction Batch                                     #
#                                                                             #
###############################################################################


class RestrictionBatch(set):
    """Class for operations on more than one enzyme."""

    def __init__(self, first=(), suppliers=()):
        """Initialize empty RB or pre-fill with enzymes (from supplier)."""
        first = [self.format(x) for x in first]
        first += [eval(x) for n in suppliers for x in suppliers_dict[n][1]]
        set.__init__(self, first)
        self.mapping = dict.fromkeys(self)
        self.already_mapped = None
        self.suppliers = [x for x in suppliers if x in suppliers_dict]

    def __str__(self):
        """Return a readable representation of the ``RestrictionBatch``."""
        if len(self) < 5:
            return "+".join(self.elements())
        else:
            return "...".join(
                ("+".join(self.elements()[:2]), "+".join(self.elements()[-2:]))
            )

    def __repr__(self):
        """Represent ``RestrictionBatch`` class as a string for debugging."""
        return f"RestrictionBatch({self.elements()})"

    def __contains__(self, other):
        """Implement ``in`` for ``RestrictionBatch``."""
        try:
            other = self.format(other)
        except ValueError:  # other is not a restriction enzyme
            return False
        return set.__contains__(self, other)

    def __div__(self, other):
        """Override '/' operator to use as search method."""
        return self.search(other)

    def __rdiv__(self, other):
        """Override division with reversed operands to use as search method."""
        return self.search(other)

    def __truediv__(self, other):
        """Override Python 3 division operator to use as search method.

        Like __div__.
        """
        return self.search(other)

    def __rtruediv__(self, other):
        """As __truediv___, with reversed operands.

        Like __rdiv__.
        """
        return self.search(other)

    def get(self, enzyme, add=False):
        """Check if enzyme is in batch and return it.

        If add is True and enzyme is not in batch add enzyme to batch.
        If add is False (which is the default) only return enzyme.
        If enzyme is not a RestrictionType or can not be evaluated to
        a RestrictionType, raise a ValueError.
        """
        e = self.format(enzyme)
        if e in self:
            return e
        elif add:
            self.add(e)
            return e
        else:
            raise ValueError(f"enzyme {e.__name__} is not in RestrictionBatch")

    def lambdasplit(self, func):
        """Filter enzymes in batch with supplied function.

        The new batch will contain only the enzymes for which
        func return True.
        """
        d = list(filter(func, self))
        new = RestrictionBatch()
        new._data = dict(zip(d, [True] * len(d)))
        return new

    def add_supplier(self, letter):
        """Add all enzymes from a given supplier to batch.

        letter represents the suppliers as defined in the dictionary
        RestrictionDictionary.suppliers
        Returns None.
        Raise a KeyError if letter is not a supplier code.
        """
        supplier = suppliers_dict[letter]
        self.suppliers.append(letter)
        for x in supplier[1]:
            self.add_nocheck(eval(x))

    def current_suppliers(self):
        """List the current suppliers for the restriction batch.

        Return a sorted list of the suppliers which have been used to
        create the batch.
        """
        suppl_list = sorted(suppliers_dict[x][0] for x in self.suppliers)
        return suppl_list

    def __iadd__(self, other):
        """Override '+=' for use with sets.

        b += other -> add other to b, check the type of other.
        """
        self.add(other)
        return self

    def __add__(self, other):
        """Override '+' for use with sets.

        b + other -> new RestrictionBatch.
        """
        new = self.__class__(self)
        new.add(other)
        return new

    def remove(self, other):
        """Remove enzyme from restriction batch.

        Safe set.remove method. Verify that other is a RestrictionType or can
        be evaluated to a RestrictionType.
        Raise a ValueError if other can not be evaluated to a RestrictionType.
        Raise a KeyError if other is not in B.
        """
        return set.remove(self, self.format(other))

    def add(self, other):
        """Add a restriction enzyme to the restriction batch.

        Safe set.add method. Verify that other is a RestrictionType or can be
        evaluated to a RestrictionType.
        Raise a ValueError if other can not be evaluated to a RestrictionType.
        """
        return set.add(self, self.format(other))

    def add_nocheck(self, other):
        """Add restriction enzyme to batch without checking its type."""
        return set.add(self, other)

    def format(self, y):
        """Evaluate enzyme (name) and return it (as RestrictionType).

        If y is a RestrictionType return y.
        If y can be evaluated to a RestrictionType return eval(y).
        Raise a ValueError in all other case.
        """
        try:
            if isinstance(y, RestrictionType):
                return y
            elif isinstance(eval(str(y)), RestrictionType):
                return eval(y)
        except (NameError, SyntaxError):
            pass
        raise ValueError(f"{y.__class__} is not a RestrictionType")

    def is_restriction(self, y):
        """Return if enzyme (name) is a known enzyme.

        True if y or eval(y) is a RestrictionType.
        """
        return isinstance(y, RestrictionType) or isinstance(
            eval(str(y)), RestrictionType
        )

    def split(self, *classes, **bool):
        """Extract enzymes of a certain class and put in new RestrictionBatch.

        It works but it is slow, so it has really an interest when splitting
        over multiple conditions.
        """

        def splittest(element):
            for klass in classes:
                b = bool.get(klass.__name__, True)
                if issubclass(element, klass):
                    if b:
                        continue
                    else:
                        return False
                elif b:
                    return False
                else:
                    continue
            return True

        d = list(filter(splittest, self))
        new = RestrictionBatch()
        new._data = dict(zip(d, [True] * len(d)))
        return new

    def elements(self):
        """List the enzymes of the RestrictionBatch as list of strings.

        Give all the names of the enzymes in B sorted alphabetically.
        """
        return sorted(str(e) for e in self)

    def as_string(self):
        """List the names of the enzymes of the RestrictionBatch.

        Return a list of the name of the elements of the batch.
        """
        return [str(e) for e in self]

    @classmethod
    def suppl_codes(cls):
        """Return a dictionary with supplier codes.

        Letter code for the suppliers.
        """
        supply = {k: v[0] for k, v in suppliers_dict.items()}
        return supply

    @classmethod
    def show_codes(cls):
        """Print a list of supplier codes."""
        supply = [" = ".join(i) for i in cls.suppl_codes().items()]
        print("\n".join(supply))

    def search(self, dna, linear=True):
        """Return a dic of cutting sites in the seq for the batch enzymes."""
        #
        #   here we replace the search method of the individual enzymes
        #   with one unique testing method.
        #
        if not hasattr(self, "already_mapped"):
            # TODO - Why does this happen!
            # Try the "doctest" at the start of PrintFormat.py
            self.already_mapped = None
        if isinstance(dna, DNA):
            # For the searching, we just care about the sequence as a string,
            # if that is the same we can use the cached search results.
            # At the time of writing, Seq == method isn't implemented,
            # and therefore does object identity which is stricter.
            if (str(dna), linear) == self.already_mapped:
                return self.mapping
            else:
                self.already_mapped = str(dna), linear
                fseq = FormattedSeq(dna, linear)
                self.mapping = {x: x.search(fseq) for x in self}
                return self.mapping
        elif isinstance(dna, FormattedSeq):
            if (str(dna), dna.linear) == self.already_mapped:
                return self.mapping
            else:
                self.already_mapped = str(dna), dna.linear
                self.mapping = {x: x.search(dna) for x in self}
                return self.mapping
        raise TypeError(f"Expected Seq or MutableSeq instance, got {type(dna)} instead")


###############################################################################
#                                                                             #
#                       Restriction Analysis                                  #
#                                                                             #
###############################################################################

_empty_DNA = DNA("")
_restrictionbatch = RestrictionBatch()


class Analysis(RestrictionBatch, PrintFormat):
    """Provide methods for enhanced analysis and pretty printing."""

    def __init__(
        self, restrictionbatch=_restrictionbatch, sequence=_empty_DNA, linear=True
    ):
        """Initialize an Analysis with RestrictionBatch and sequence.

        For most of the methods of this class if a dictionary is given it will
        be used as the base to calculate the results.
        If no dictionary is given a new analysis using the RestrictionBatch
        which has been given when the Analysis class has been instantiated,
        will be carried out and used.
        """
        RestrictionBatch.__init__(self, restrictionbatch)
        self.rb = restrictionbatch
        self.sequence = sequence
        self.linear = linear
        if self.sequence:
            self.search(self.sequence, self.linear)

    def __repr__(self):
        """Represent ``Analysis`` class as a string."""
        return f"Analysis({self.rb!r},{self.sequence!r},{self.linear})"

    def _sub_set(self, wanted):
        """Filter result for keys which are in wanted (PRIVATE).

        Internal use only. Returns a dict.

        Screen the results through wanted set.
        Keep only the results for which the enzymes is in wanted set.
        """
        # It seems that this method is not used in the whole class!
        return {k: v for k, v in self.mapping.items() if k in wanted}

    def _boundaries(self, start, end):
        """Set boundaries to correct values (PRIVATE).

        Format the boundaries for use with the methods that limit the
        search to only part of the sequence given to analyse.
        """
        if not isinstance(start, int):
            raise TypeError(f"expected int, got {type(start)} instead")
        if not isinstance(end, int):
            raise TypeError(f"expected int, got {type(end)} instead")
        if start < 1:  # Looks like this tries to do python list like indexing
            start += len(self.sequence)
        if end < 1:
            end += len(self.sequence)
        if start < end:
            pass
        else:
            start, end = end, start
        if start < end:
            return start, end, self._test_normal

    def _test_normal(self, start, end, site):
        """Test if site is between start and end (PRIVATE).

        Internal use only
        """
        return start <= site < end

    def _test_reverse(self, start, end, site):
        """Test if site is between end and start, for circular sequences (PRIVATE).

        Internal use only.
        """
        return start <= site <= len(self.sequence) or 1 <= site < end

    def format_output(self, dct=None, title="", s1=""):
        """Collect data and pass to PrintFormat.

        If dct is not given the full dictionary is used.
        """
        if not dct:
            dct = self.mapping
        return PrintFormat.format_output(self, dct, title, s1)

    def print_that(self, dct=None, title="", s1=""):
        """Print the output of the analysis.

        If dct is not given the full dictionary is used.
        s1: Title for non-cutting enzymes
        This method prints the output of A.format_output() and it is here
        for backwards compatibility.
        """
        print(self.format_output(dct, title, s1))

    def change(self, **what):
        """Change parameters of print output.

        It is possible to change the width of the shell by setting
        self.ConsoleWidth to what you want.
        self.NameWidth refer to the maximal length of the enzyme name.

        Changing one of these parameters here might not give the results
        you expect. In which case, you can settle back to a 80 columns shell
        or try to change self.Cmodulo and self.PrefWidth in PrintFormat until
        you get it right.
        """
        for k, v in what.items():
            if k in ("NameWidth", "ConsoleWidth"):
                setattr(self, k, v)
                self.Cmodulo = self.ConsoleWidth % self.NameWidth
                self.PrefWidth = self.ConsoleWidth - self.Cmodulo
            elif k == "sequence":
                setattr(self, "sequence", v)
                self.search(self.sequence, self.linear)
            elif k == "rb":
                self = Analysis.__init__(self, v, self.sequence, self.linear)
            elif k == "linear":
                setattr(self, "linear", v)
                self.search(self.sequence, v)
            elif k in ("Indent", "Maxsize"):
                setattr(self, k, v)
            elif k in ("Cmodulo", "PrefWidth"):
                raise AttributeError(
                    f"To change {k}, change NameWidth and/or ConsoleWidth"
                )
            else:
                raise AttributeError(f"Analysis has no attribute {k}")

    def full(self, linear=True):
        """Perform analysis with all enzymes of batch and return all results.

        Full Restriction Map of the sequence, as a dictionary.
        """
        return self.mapping

    def blunt(self, dct=None):
        """Return only cuts that have blunt ends."""
        if not dct:
            dct = self.mapping
        return {k: v for k, v in dct.items() if k.is_blunt()}

    def overhang5(self, dct=None):
        """Return only cuts that have 5' overhangs."""
        if not dct:
            dct = self.mapping
        return {k: v for k, v in dct.items() if k.is_5overhang()}

    def overhang3(self, dct=None):
        """Return only cuts that have 3' overhangs."""
        if not dct:
            dct = self.mapping
        return {k: v for k, v in dct.items() if k.is_3overhang()}

    def defined(self, dct=None):
        """Return only results from enzymes that produce defined overhangs."""
        if not dct:
            dct = self.mapping
        return {k: v for k, v in dct.items() if k.is_defined()}

    def with_sites(self, dct=None):
        """Return only results from enzyme with at least one cut."""
        if not dct:
            dct = self.mapping
        return {k: v for k, v in dct.items() if v}

    def without_site(self, dct=None):
        """Return only results from enzymes that don't cut the sequence."""
        if not dct:
            dct = self.mapping
        return {k: v for k, v in dct.items() if not v}

    def with_N_sites(self, N, dct=None):
        """Return only results from enzymes that cut the sequence N times."""
        if not dct:
            dct = self.mapping
        return {k: v for k, v in dct.items() if len(v) == N}

    def with_number_list(self, list, dct=None):
        """Return only results from enzymes that cut (x,y,z,...) times."""
        if not dct:
            dct = self.mapping
        return {k: v for k, v in dct.items() if len(v) in list}

    def with_name(self, names, dct=None):
        """Return only results from enzymes which names are listed."""
        for i, enzyme in enumerate(names):
            if enzyme not in AllEnzymes:
                warnings.warn(f"no data for the enzyme: {enzyme}", BiopythonWarning)
                del names[i]
        if not dct:
            return RestrictionBatch(names).search(self.sequence, self.linear)
        return {n: dct[n] for n in names if n in dct}

    def with_site_size(self, site_size, dct=None):
        """Return only results form enzymes with a given site size."""
        sites = [name for name in self if name.size == site_size]
        if not dct:
            return RestrictionBatch(sites).search(self.sequence)
        return {k: v for k, v in dct.items() if k in site_size}

    def only_between(self, start, end, dct=None):
        """Return only results from enzymes that only cut within start, end."""
        start, end, test = self._boundaries(start, end)
        if not dct:
            dct = self.mapping
        d = dict(dct)
        for key, sites in dct.items():
            if not sites:
                del d[key]
                continue
            for site in sites:
                if test(start, end, site):
                    continue
                else:
                    del d[key]
                    break
        return d

    def between(self, start, end, dct=None):
        """Return only results from enzymes that cut at least within borders.

        Enzymes that cut the sequence at least in between start and end.
        They may cut outside as well.
        """
        start, end, test = self._boundaries(start, end)
        d = {}
        if not dct:
            dct = self.mapping
        for key, sites in dct.items():
            for site in sites:
                if test(start, end, site):
                    d[key] = sites
                    break
                continue
        return d

    def show_only_between(self, start, end, dct=None):
        """Return only results from within start, end.

        Enzymes must cut inside start/end and may also cut outside. However,
        only the cutting positions within start/end will be returned.
        """
        d = []
        if start <= end:
            d = [
                (k, [vv for vv in v if start <= vv <= end])
                for k, v in self.between(start, end, dct).items()
            ]
        else:
            d = [
                (k, [vv for vv in v if start <= vv or vv <= end])
                for k, v in self.between(start, end, dct).items()
            ]
        return dict(d)

    def only_outside(self, start, end, dct=None):
        """Return only results from enzymes that only cut outside start, end.

        Enzymes that cut the sequence outside of the region
        in between start and end but do not cut inside.
        """
        start, end, test = self._boundaries(start, end)
        if not dct:
            dct = self.mapping
        d = dict(dct)
        for key, sites in dct.items():
            if not sites:
                del d[key]
                continue
            for site in sites:
                if test(start, end, site):
                    del d[key]
                    break
                else:
                    continue
        return d

    def outside(self, start, end, dct=None):
        """Return only results from enzymes that at least cut outside borders.

        Enzymes that cut outside the region in between start and end.
        They may cut inside as well.
        """
        start, end, test = self._boundaries(start, end)
        if not dct:
            dct = self.mapping
        d = {}
        for key, sites in dct.items():
            for site in sites:
                if test(start, end, site):
                    continue
                else:
                    d[key] = sites
                    break
        return d

    def do_not_cut(self, start, end, dct=None):
        """Return only results from enzymes that don't cut between borders."""
        if not dct:
            dct = self.mapping
        d = self.without_site()
        d.update(self.only_outside(start, end, dct))
        return d


#
#   The restriction enzyme classes are created dynamically when the module is
#   imported. Here is the magic which allow the creation of the
#   restriction-enzyme classes.
#
#   The reason for the two dictionaries in Restriction_Dictionary
#   one for the types (which will be called pseudo-type as they really
#   correspond to the values that instances of RestrictionType can take)
#   and one for the enzymes is efficiency as the bases are evaluated
#   once per pseudo-type.
#
#   However Restriction is still a very inefficient module at import. But
#   remember that around 660 classes (which is more or less the size of Rebase)
#   have to be created dynamically. However, this processing take place only
#   once.
#   This inefficiency is however largely compensated by the use of metaclass
#   which provide a very efficient layout for the class themselves mostly
#   alleviating the need of if/else loops in the class methods.
#
#   It is essential to run Restriction with doc string optimisation (-OO
#   switch) as the doc string of 660 classes take a lot of processing.
#
CommOnly = RestrictionBatch()  # commercial enzymes
NonComm = RestrictionBatch()  # not available commercially
for TYPE, (bases, enzymes) in typedict.items():
    #
    #   The keys are the pseudo-types TYPE (stored as type1, type2...)
    #   The names are not important and are only present to differentiate
    #   the keys in the dict. All the pseudo-types are in fact RestrictionType.
    #   These names will not be used after and the pseudo-types are not
    #   kept in the locals() dictionary. It is therefore impossible to
    #   import them.
    #   Now, if you have look at the dictionary, you will see that not all the
    #   types are present as those without corresponding enzymes have been
    #   removed by Dictionary_Builder().
    #
    #   The values are tuples which contain
    #   as first element a tuple of bases (as string) and
    #   as second element the names of the enzymes.
    #
    #   First eval the bases.
    #
    bases = tuple(eval(x) for x in bases)
    #
    #   now create the particular value of RestrictionType for the classes
    #   in enzymes.
    #
    T = type.__new__(RestrictionType, "RestrictionType", bases, {})
    for k in enzymes:
        #
        #   Now, we go through all the enzymes and assign them their type.
        #   enzymedict[k] contains the values of the attributes for this
        #   particular class (self.site, self.ovhg,....).
        #
        newenz = T(k, bases, enzymedict[k])
        #
        #   we add the enzymes to the corresponding batch.
        #
        #   No need to verify the enzyme is a RestrictionType -> add_nocheck
        #
        if newenz.is_comm():
            CommOnly.add_nocheck(newenz)
        else:
            NonComm.add_nocheck(newenz)
#
#   AllEnzymes is a RestrictionBatch with all the enzymes from Rebase.
#
AllEnzymes = RestrictionBatch(CommOnly)
AllEnzymes.update(NonComm)
#
#   Now, place the enzymes in locals so they can be imported.
#
names = [str(x) for x in AllEnzymes]
locals().update(dict(zip(names, AllEnzymes)))
__all__ = (
    "FormattedSeq",
    "Analysis",
    "RestrictionBatch",
    "AllEnzymes",
    "CommOnly",
    "NonComm",
) + tuple(names)
del k, enzymes, TYPE, bases, names
