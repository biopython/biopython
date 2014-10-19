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
        ----------------------------------------------------------------------------
            AbstractCut implements some methods that are common to all enzymes.
        ----------------------------------------------------------------------------
            NoCut, OneCut,TwoCuts   represent the number of double strand cuts
                                    produced by the enzyme.
                                    they correspond to the 4th field of the rebase
                                    record emboss_e.NNN.
                    0->NoCut    : the enzyme is not characterised.
                    2->OneCut   : the enzyme produce one double strand cut.
                    4->TwoCuts  : two double strand cuts.
        ----------------------------------------------------------------------------
            Meth_Dep, Meth_Undep    represent the methylation susceptibility to
                                    the enzyme.
                                    Not implemented yet.
        ----------------------------------------------------------------------------
            Palindromic,            if the site is palindromic or not.
            NotPalindromic          allow some optimisations of the code.
                                    No need to check the reverse strand
                                    with palindromic sites.
        ----------------------------------------------------------------------------
            Unknown, Blunt,         represent the overhang.
            Ov5, Ov3                Unknown is here for symetry reasons and
                                    correspond to enzymes that are not characterised
                                    in rebase.
        ----------------------------------------------------------------------------
            Defined, Ambiguous,     represent the sequence of the overhang.
            NotDefined
                                    NotDefined is for enzymes not characterised in
                                    rebase.

                                    Defined correspond to enzymes that display a
                                    constant overhang whatever the sequence.
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
                   So the enzyme ApoI (RAATTY) is defined even if its restriction
                   site is ambiguous.

                        ApoI R^AATTY -> overhang : AATT -> Defined
                             YTTAA^R
                   Accordingly, blunt enzymes are always Defined even
                   when they cut outside their restriction site.
        ----------------------------------------------------------------------------
            Not_available,          as found in rebase file emboss_r.NNN files.
            Commercially_available
                                    allow the selection of the enzymes according to
                                    their suppliers to reduce the quantity
                                    of results.
                                    Also will allow the implementation of buffer
                                    compatibility tables. Not implemented yet.

                                    the list of suppliers is extracted from
                                    emboss_s.NNN
        ----------------------------------------------------------------------------
"""

from __future__ import print_function
from Bio._py3k import zip
from Bio._py3k import filter
from Bio._py3k import range

import re
import itertools

from Bio.Seq import Seq, MutableSeq
from Bio.Alphabet import IUPAC

from Bio.Restriction.Restriction_Dictionary import rest_dict as enzymedict
from Bio.Restriction.Restriction_Dictionary import typedict
from Bio.Restriction.Restriction_Dictionary import suppliers as suppliers_dict
# TODO: Consider removing this wildcard import.
from Bio.Restriction.RanaConfig import *
from Bio.Restriction.PrintFormat import PrintFormat

__docformat__ = "restructuredtext en"

# Used to use Bio.Restriction.DNAUtils.check_bases (and expose it under this
# namespace), but have deprecated that module.
def _check_bases(seq_string):
    """Check characters in a string (PRIVATE).

    Remove digits and white space present in string. Allows any valid ambiguous
    IUPAC DNA single letters codes (ABCDGHKMNRSTVWY, lower case are converted).

    Other characters (e.g. symbols) trigger a TypeError.

    Returns the string WITH A LEADING SPACE (!). This is for backwards
    compatibility, and may in part be explained by the fact that
    Bio.Restriction doesn't use zero based counting.
    """
    # Remove white space and make upper case:
    seq_string = "".join(seq_string.split()).upper()
    # Remove digits
    for c in "0123456789":
        seq_string = seq_string.replace(c, "")
    # Check only allowed IUPAC letters
    if not set(seq_string).issubset(set("ABCDGHKMNRSTVWY")):
        raise TypeError("Invalid character found in %s" % repr(seq_string))
    return " " + seq_string


matching = {'A': 'ARWMHVDN', 'C': 'CYSMHBVN', 'G': 'GRSKBVDN',
            'T': 'TYWKHBDN', 'R': 'ABDGHKMNSRWV', 'Y': 'CBDHKMNSTWVY',
            'W': 'ABDHKMNRTWVY', 'S': 'CBDGHKMNSRVY', 'M': 'ACBDHMNSRWVY',
            'K': 'BDGHKNSRTWVY', 'H': 'ACBDHKMNSRTWVY',
            'B': 'CBDGHKMNSRTWVY', 'V': 'ACBDGHKMNSRWVY',
            'D': 'ABDGHKMNSRTWVY', 'N': 'ACBDGHKMNSRTWVY'}

DNA = Seq


class FormattedSeq(object):
    """FormattedSeq(seq, [linear=True])-> new FormattedSeq.

    Translate a Bio.Seq into a formatted sequence to be used with Restriction.

    Roughly:
        remove anything which is not IUPAC alphabet and then add a space
        in front of the sequence to get a biological index instead of a
        python index (i.e. index of the first base is 1 not 0).

        Retains information about the shape of the molecule linear (default)
        or circular. Restriction sites are search over the edges of circular
        sequence."""

    def __init__(self, seq, linear=True):
        """FormattedSeq(seq, [linear=True])-> new FormattedSeq.

        seq is either a Bio.Seq, Bio.MutableSeq or a FormattedSeq.
        if seq is a FormattedSeq, linear will have no effect on the
        shape of the sequence."""
        if isinstance(seq, Seq) or isinstance(seq, MutableSeq):
            stringy = str(seq)
            self.lower = stringy.islower()
            # Note this adds a leading space to the sequence (!)
            self.data = _check_bases(stringy)
            self.linear = linear
            self.klass = seq.__class__
            self.alphabet = seq.alphabet
        elif isinstance(seq, FormattedSeq):
            self.lower = seq.lower
            self.data = seq.data
            self.linear = seq.linear
            self.alphabet = seq.alphabet
            self.klass = seq.klass
        else:
            raise TypeError('expected Seq or MutableSeq, got %s' % type(seq))

    def __len__(self):
        return len(self.data) - 1

    def __repr__(self):
        return 'FormattedSeq(%s, linear=%s)' %(repr(self[1:]), repr(self.linear))

    def __eq__(self, other):
        if isinstance(other, FormattedSeq):
            if repr(self) == repr(other):
                return True
            else:
                return False
        return False

    def circularise(self):
        """FS.circularise() -> circularise FS"""
        self.linear = False
        return

    def linearise(self):
        """FS.linearise() -> linearise FS"""
        self.linear = True
        return

    def to_linear(self):
        """FS.to_linear() -> new linear FS instance"""
        new = self.__class__(self)
        new.linear = True
        return new

    def to_circular(self):
        """FS.to_circular() -> new circular FS instance"""
        new = self.__class__(self)
        new.linear = False
        return new

    def is_linear(self):
        """FS.is_linear() -> bool.

        True if the sequence will analysed as a linear sequence."""
        return self.linear

    def finditer(self, pattern, size):
        """FS.finditer(pattern, size) -> list.

        return a list of pattern into the sequence.
        the list is made of tuple (location, pattern.group).
        the latter is used with non palindromic sites.
        pattern is the regular expression pattern corresponding to the
        enzyme restriction site.
        size is the size of the restriction enzyme recognition-site size."""
        if self.is_linear():
            data = self.data
        else:
            data = self.data + self.data[1:size]
        return [(i.start(), i.group) for i in re.finditer(pattern, data)]

    def __getitem__(self, i):
        if self.lower:
            return self.klass((self.data[i]).lower(), self.alphabet)
        return self.klass(self.data[i], self.alphabet)


class RestrictionType(type):
    """RestrictionType. Type from which derives all enzyme classes.

    Implement the operator methods."""

    def __init__(cls, name='', bases=(), dct={}):
        """RE(name, bases, dct) -> RestrictionType instance.

        Not intended to be used in normal operation. The enzymes are
        instantiated when importing the module.

        see below."""
        if "-" in name:
            raise ValueError("Problem with hyphen in %s as enzyme name"
                             % repr(name))
        # 2011/11/26 - Nobody knows what this call was supposed to accomplish,
        # but all unit tests seem to pass without it.
        # super(RestrictionType, cls).__init__(cls, name, bases, dct)
        try:
            cls.compsite = re.compile(cls.compsite)
        except Exception as err:
            raise ValueError("Problem with regular expression, re.compiled(%s)"
                             % repr(cls.compsite))

    def __add__(cls, other):
        """RE.__add__(other) -> RestrictionBatch().

        if other is an enzyme returns a batch of the two enzymes.
        if other is already a RestrictionBatch add enzyme to it."""
        if isinstance(other, RestrictionType):
            return RestrictionBatch([cls, other])
        elif isinstance(other, RestrictionBatch):
            return other.add_nocheck(cls)
        else:
            raise TypeError

    def __div__(cls, other):
        """RE.__div__(other) -> list.

        RE/other
        returns RE.search(other)."""
        return cls.search(other)

    def __rdiv__(cls, other):
        """RE.__rdiv__(other) -> list.

        other/RE
        returns RE.search(other)."""
        return cls.search(other)

    def __truediv__(cls, other):
        """RE.__truediv__(other) -> list.

        RE/other
        returns RE.search(other)."""
        return cls.search(other)

    def __rtruediv__(cls, other):
        """RE.__rtruediv__(other) -> list.

        other/RE
        returns RE.search(other)."""
        return cls.search(other)

    def __floordiv__(cls, other):
        """RE.__floordiv__(other) -> list.

        RE//other
        returns RE.catalyse(other)."""
        return cls.catalyse(other)

    def __rfloordiv__(cls, other):
        """RE.__rfloordiv__(other) -> list.

        other//RE
        returns RE.catalyse(other)."""
        return cls.catalyse(other)

    def __str__(cls):
        """RE.__str__() -> str.

        return the name of the enzyme."""
        return cls.__name__

    def __repr__(cls):
        """RE.__repr__() -> str.

        used with eval or exec will instantiate the enzyme."""
        return "%s" % cls.__name__

    def __len__(cls):
        """RE.__len__() -> int.

        length of the recognition site."""
        return cls.size

    def __hash__(cls):
        # Python default is to use id(...)
        # This is consistent with the __eq__ implementation
        return id(cls)

    def __eq__(cls, other):
        """RE == other -> bool

        True if RE and other are the same enzyme.

        Specifically this checks they are the same Python object.
        """
        # assert (id(cls)==id(other)) == (other is cls) == (cls is other)
        return id(cls)==id(other)

    def __ne__(cls, other):
        """RE != other -> bool.
        isoschizomer strict, same recognition site, same restriction -> False
        all the other-> True

        WARNING - This is not the inverse of the __eq__ method.
        """
        if not isinstance(other, RestrictionType):
            return True
        elif cls.charac == other.charac:
            return False
        else:
            return True

    def __rshift__(cls, other):
        """RE >> other -> bool.

        neoschizomer : same recognition site, different restriction. -> True
        all the others :                                             -> False"""
        if not isinstance(other, RestrictionType):
            return False
        elif cls.site == other.site and cls.charac != other.charac:
            return True
        else:
            return False

    def __mod__(cls, other):
        """a % b -> bool.

        Test compatibility of the overhang of a and b.
        True if a and b have compatible overhang."""
        if not isinstance(other, RestrictionType):
            raise TypeError(
                  'expected RestrictionType, got %s instead' % type(other))
        return cls._mod1(other)

    def __ge__(cls, other):
        """a >= b -> bool.

        a is greater or equal than b if the a site is longer than b site.
        if their site have the same length sort by alphabetical order of their
        names."""
        if not isinstance(other, RestrictionType):
            raise NotImplementedError
        if len(cls) > len(other):
            return True
        elif cls.size == len(other) and cls.__name__ >= other.__name__:
            return True
        else:
            return False

    def __gt__(cls, other):
        """a > b -> bool.

        sorting order:
                    1. size of the recognition site.
                    2. if equal size, alphabetical order of the names."""
        if not isinstance(other, RestrictionType):
            raise NotImplementedError
        if len(cls) > len(other):
            return True
        elif cls.size == len(other) and cls.__name__ > other.__name__:
            return True
        else:
            return False

    def __le__(cls, other):
        """a <= b -> bool.

        sorting order:
                    1. size of the recognition site.
                    2. if equal size, alphabetical order of the names."""
        if not isinstance(other, RestrictionType):
            raise NotImplementedError
        elif len(cls) < len(other):
            return True
        elif len(cls) == len(other) and cls.__name__ <= other.__name__:
            return True
        else:
            return False

    def __lt__(cls, other):
        """a < b -> bool.

        sorting order:
                    1. size of the recognition site.
                    2. if equal size, alphabetical order of the names."""
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

    For internal use only. Not meant to be instantiate."""

    @classmethod
    def search(cls, dna, linear=True):
        """RE.search(dna, linear=True) -> list.

        return a list of all the site of RE in dna. Compensate for circular
        sequences and so on.

        dna must be a Bio.Seq.Seq instance or a Bio.Seq.MutableSeq instance.

        if linear is False, the restriction sites than span over the boundaries
        will be included.

        The positions are the first base of the 3' fragment,
        i.e. the first base after the position the enzyme will cut. """
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
    def all_suppliers(self):
        """RE.all_suppliers -> print all the suppliers of R"""
        supply = sorted(x[0] for x in suppliers_dict.values())
        print(",\n".join(supply))
        return

    @classmethod
    def is_equischizomer(self, other):
        """RE.is_equischizomers(other) -> bool.

        True if other is an isoschizomer of RE.
        False else.

        equischizomer <=> same site, same position of restriction."""
        return not self != other

    @classmethod
    def is_neoschizomer(self, other):
        """RE.is_neoschizomers(other) -> bool.

        True if other is an isoschizomer of RE.
        False else.

        neoschizomer <=> same site, different position of restriction."""
        return self >> other

    @classmethod
    def is_isoschizomer(self, other):
        """RE.is_isoschizomers(other) -> bool.

        True if other is an isoschizomer of RE.
        False else.

        isoschizomer <=> same site."""
        return (not self != other) or self >> other

    @classmethod
    def equischizomers(self, batch=None):
        """RE.equischizomers([batch]) -> list.

        return a tuple of all the isoschizomers of RE.
        if batch is supplied it is used instead of the default AllEnzymes.

        equischizomer <=> same site, same position of restriction."""
        if not batch:
            batch = AllEnzymes
        r = [x for x in batch if not self != x]
        i = r.index(self)
        del r[i]
        r.sort()
        return r

    @classmethod
    def neoschizomers(self, batch=None):
        """RE.neoschizomers([batch]) -> list.

        return a tuple of all the neoschizomers of RE.
        if batch is supplied it is used instead of the default AllEnzymes.

        neoschizomer <=> same site, different position of restriction."""
        if not batch:
            batch = AllEnzymes
        r = sorted(x for x in batch if self >> x)
        return r

    @classmethod
    def isoschizomers(self, batch=None):
        """RE.isoschizomers([batch]) -> list.

        return a tuple of all the equischizomers and neoschizomers of RE.
        if batch is supplied it is used instead of the default AllEnzymes."""
        if not batch:
            batch = AllEnzymes
        r = [x for x in batch if (self >> x) or (not self != x)]
        i = r.index(self)
        del r[i]
        r.sort()
        return r

    @classmethod
    def frequency(self):
        """RE.frequency() -> int.

        frequency of the site."""
        return self.freq


class NoCut(AbstractCut):
    """Implement the methods specific to the enzymes that do not cut.

    These enzymes are generally enzymes that have been only partially
    characterised and the way they cut the DNA is unknow or enzymes for
    which the pattern of cut is to complex to be recorded in Rebase
    (ncuts values of 0 in emboss_e.###).

    When using search() with these enzymes the values returned are at the start of
    the restriction site.

    Their catalyse() method returns a TypeError.

    Unknown and NotDefined are also part of the base classes of these enzymes.

    Internal use only. Not meant to be instantiated."""

    @classmethod
    def cut_once(self):
        """RE.cut_once() -> bool.

        True if the enzyme cut the sequence one time on each strand."""
        return False

    @classmethod
    def cut_twice(self):
        """RE.cut_twice() -> bool.

        True if the enzyme cut the sequence twice on each strand."""
        return False

    @classmethod
    def _modify(self, location):
        """RE._modify(location) -> int.

        for internal use only.

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
        yield location

    @classmethod
    def _rev_modify(self, location):
        """RE._rev_modify(location) -> generator of int.

        for internal use only.

        as _modify for site situated on the antiparallel strand when the
        enzyme is not palindromic
        """
        yield location

    @classmethod
    def characteristic(self):
        """RE.characteristic() -> tuple.

        the tuple contains the attributes:
            fst5 -> first 5' cut ((current strand) or None
            fst3 -> first 3' cut (complementary strand) or None
            scd5 -> second 5' cut (current strand) or None
            scd5 -> second 3' cut (complementary strand) or None
            site -> recognition site."""
        return None, None, None, None, self.site


class OneCut(AbstractCut):
    """Implement the methods specific to the enzymes that cut the DNA only once

    Correspond to ncuts values of 2 in emboss_e.###

    Internal use only. Not meant to be instantiated."""

    @classmethod
    def cut_once(self):
        """RE.cut_once() -> bool.

        True if the enzyme cut the sequence one time on each strand."""
        return True

    @classmethod
    def cut_twice(self):
        """RE.cut_twice() -> bool.

        True if the enzyme cut the sequence twice on each strand."""
        return False

    @classmethod
    def _modify(self, location):
        """RE._modify(location) -> int.

        for internal use only.

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
        yield location + self.fst5

    @classmethod
    def _rev_modify(self, location):
        """RE._rev_modify(location) -> generator of int.

        for internal use only.

        as _modify for site situated on the antiparallel strand when the
        enzyme is not palindromic
        """
        yield location - self.fst3

    @classmethod
    def characteristic(self):
        """RE.characteristic() -> tuple.

        the tuple contains the attributes:
            fst5 -> first 5' cut ((current strand) or None
            fst3 -> first 3' cut (complementary strand) or None
            scd5 -> second 5' cut (current strand) or None
            scd5 -> second 3' cut (complementary strand) or None
            site -> recognition site."""
        return self.fst5, self.fst3, None, None, self.site


class TwoCuts(AbstractCut):
    """Implement the methods specific to the enzymes that cut the DNA twice

    Correspond to ncuts values of 4 in emboss_e.###

    Internal use only. Not meant to be instantiated."""

    @classmethod
    def cut_once(self):
        """RE.cut_once() -> bool.

        True if the enzyme cut the sequence one time on each strand."""
        return False

    @classmethod
    def cut_twice(self):
        """RE.cut_twice() -> bool.

        True if the enzyme cut the sequence twice on each strand."""
        return True

    @classmethod
    def _modify(self, location):
        """RE._modify(location) -> int.

        for internal use only.

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
        yield location + self.fst5
        yield location + self.scd5

    @classmethod
    def _rev_modify(self, location):
        """RE._rev_modify(location) -> generator of int.

        for internal use only.

        as _modify for site situated on the antiparallel strand when the
        enzyme is not palindromic
        """
        yield location - self.fst3
        yield location - self.scd3

    @classmethod
    def characteristic(self):
        """RE.characteristic() -> tuple.

        the tuple contains the attributes:
            fst5 -> first 5' cut ((current strand) or None
            fst3 -> first 3' cut (complementary strand) or None
            scd5 -> second 5' cut (current strand) or None
            scd5 -> second 3' cut (complementary strand) or None
            site -> recognition site."""
        return self.fst5, self.fst3, self.scd5, self.scd3, self.site


class Meth_Dep(AbstractCut):
    """Implement the information about methylation.

    Enzymes of this class possess a site which is methylable."""

    @classmethod
    def is_methylable(self):
        """RE.is_methylable() -> bool.

        True if the recognition site is a methylable."""
        return True


class Meth_Undep(AbstractCut):
    """Implement information about methylation sensitibility.

    Enzymes of this class are not sensible to methylation."""

    @classmethod
    def is_methylable(self):
        """RE.is_methylable() -> bool.

        True if the recognition site is a methylable."""
        return False


class Palindromic(AbstractCut):
    """Implement the methods specific to the enzymes which are palindromic

    palindromic means : the recognition site and its reverse complement are
                        identical.
    Remarks     : an enzyme with a site CGNNCG is palindromic even if some
                  of the sites that it will recognise are not.
                  for example here : CGAACG

    Internal use only. Not meant to be instantiated."""

    @classmethod
    def _search(self):
        """RE._search() -> list.

        for internal use only.

        implement the search method for palindromic and non palindromic enzyme.
        """
        siteloc = self.dna.finditer(self.compsite, self.size)
        self.results = [r for s, g in siteloc for r in self._modify(s)]
        if self.results:
            self._drop()
        return self.results

    @classmethod
    def is_palindromic(self):
        """RE.is_palindromic() -> bool.

        True if the recognition site is a palindrom."""
        return True


class NonPalindromic(AbstractCut):
    """Implement the methods specific to the enzymes which are not palindromic

    palindromic means : the recognition site and its reverse complement are
                        identical.

    Internal use only. Not meant to be instantiated."""

    @classmethod
    def _search(self):
        """RE._search() -> list.

        for internal use only.

        implement the search method for palindromic and non palindromic enzyme.
        """
        iterator = self.dna.finditer(self.compsite, self.size)
        self.results = []
        modif = self._modify
        revmodif = self._rev_modify
        s = str(self)
        self.on_minus = []
        for start, group in iterator:
            if group(s):
                self.results += [r for r in modif(start)]
            else:
                self.on_minus += [r for r in revmodif(start)]
        self.results += self.on_minus
        if self.results:
            self.results.sort()
            self._drop()
        return self.results

    @classmethod
    def is_palindromic(self):
        """RE.is_palindromic() -> bool.

        True if the recognition site is a palindrom."""
        return False


class Unknown(AbstractCut):
    """Implement the methods specific to the enzymes for which the overhang
    is unknown.

    These enzymes are also NotDefined and NoCut.

    Internal use only. Not meant to be instantiated."""

    @classmethod
    def catalyse(self, dna, linear=True):
        """RE.catalyse(dna, linear=True) -> tuple of DNA.
        RE.catalyze(dna, linear=True) -> tuple of DNA.

        return a tuple of dna as will be produced by using RE to restrict the
        dna.

        dna must be a Bio.Seq.Seq instance or a Bio.Seq.MutableSeq instance.

        if linear is False, the sequence is considered to be circular and the
        output will be modified accordingly."""
        raise NotImplementedError('%s restriction is unknown.'
                                  % self.__name__)
    catalyze = catalyse

    @classmethod
    def is_blunt(self):
        """RE.is_blunt() -> bool.

        True if the enzyme produces blunt end.

        see also:
            RE.is_3overhang()
            RE.is_5overhang()
            RE.is_unknown()"""
        return False

    @classmethod
    def is_5overhang(self):
        """RE.is_5overhang() -> bool.

        True if the enzyme produces 5' overhang sticky end.

        see also:
            RE.is_3overhang()
            RE.is_blunt()
            RE.is_unknown()"""
        return False

    @classmethod
    def is_3overhang(self):
        """RE.is_3overhang() -> bool.

        True if the enzyme produces 3' overhang sticky end.

        see also:
            RE.is_5overhang()
            RE.is_blunt()
            RE.is_unknown()"""
        return False

    @classmethod
    def overhang(self):
        """RE.overhang() -> str. type of overhang of the enzyme.,

        can be "3' overhang", "5' overhang", "blunt", "unknown"   """
        return 'unknown'

    @classmethod
    def compatible_end(self):
        """RE.compatible_end() -> list.

        list of all the enzymes that share compatible end with RE."""
        return []

    @classmethod
    def _mod1(self, other):
        """RE._mod1(other) -> bool.

        for internal use only

        test for the compatibility of restriction ending of RE and other."""
        return False


class Blunt(AbstractCut):
    """Implement the methods specific to the enzymes for which the overhang
    is blunt.

    The enzyme cuts the + strand and the - strand of the DNA at the same
    place.

    Internal use only. Not meant to be instantiated."""

    @classmethod
    def catalyse(self, dna, linear=True):
        """RE.catalyse(dna, linear=True) -> tuple of DNA.
        RE.catalyze(dna, linear=True) -> tuple of DNA.

        return a tuple of dna as will be produced by using RE to restrict the
        dna.

        dna must be a Bio.Seq.Seq instance or a Bio.Seq.MutableSeq instance.

        if linear is False, the sequence is considered to be circular and the
        output will be modified accordingly."""
        r = self.search(dna, linear)
        d = self.dna
        if not r:
            return d[1:],
        fragments = []
        length = len(r)-1
        if d.is_linear():
            #
            #   START of the sequence to FIRST site.
            #
            fragments.append(d[1:r[0]])
            if length:
                #
                #   if more than one site add them.
                #
                fragments += [d[r[x]:r[x+1]] for x in range(length)]
            #
            #   LAST site to END of the sequence.
            #
            fragments.append(d[r[-1]:])
        else:
            #
            #   circular : bridge LAST site to FIRST site.
            #
            fragments.append(d[r[-1]:]+d[1:r[0]])
            if not length:
                #
                #   one site we finish here.
                #
                return tuple(fragments)
            #
            #   add the others.
            #
            fragments += [d[r[x]:r[x+1]] for x in range(length)]
        return tuple(fragments)
    catalyze = catalyse

    @classmethod
    def is_blunt(self):
        """RE.is_blunt() -> bool.

        True if the enzyme produces blunt end.

        see also:
            RE.is_3overhang()
            RE.is_5overhang()
            RE.is_unknown()"""
        return True

    @classmethod
    def is_5overhang(self):
        """RE.is_5overhang() -> bool.

        True if the enzyme produces 5' overhang sticky end.

        see also:
            RE.is_3overhang()
            RE.is_blunt()
            RE.is_unknown()"""
        return False

    @classmethod
    def is_3overhang(self):
        """RE.is_3overhang() -> bool.

        True if the enzyme produces 3' overhang sticky end.

        see also:
            RE.is_5overhang()
            RE.is_blunt()
            RE.is_unknown()"""
        return False

    @classmethod
    def overhang(self):
        """RE.overhang() -> str. type of overhang of the enzyme.,

        can be "3' overhang", "5' overhang", "blunt", "unknown"   """
        return 'blunt'

    @classmethod
    def compatible_end(self, batch=None):
        """RE.compatible_end() -> list.

        list of all the enzymes that share compatible end with RE."""
        if not batch:
            batch = AllEnzymes
        r = sorted(x for x in iter(AllEnzymes) if x.is_blunt())
        return r

    @staticmethod
    def _mod1(other):
        """RE._mod1(other) -> bool.

        for internal use only

        test for the compatibility of restriction ending of RE and other."""
        return issubclass(other, Blunt)


class Ov5(AbstractCut):
    """Implement the methods specific to the enzymes for which the overhang
    is recessed in 3'.

    The enzyme cuts the + strand after the - strand of the DNA.

    Internal use only. Not meant to be instantiated."""

    @classmethod
    def catalyse(self, dna, linear=True):
        """RE.catalyse(dna, linear=True) -> tuple of DNA.
        RE.catalyze(dna, linear=True) -> tuple of DNA.

        return a tuple of dna as will be produced by using RE to restrict the
        dna.

        dna must be a Bio.Seq.Seq instance or a Bio.Seq.MutableSeq instance.

        if linear is False, the sequence is considered to be circular and the
        output will be modified accordingly."""
        r = self.search(dna, linear)
        d = self.dna
        if not r:
            return d[1:],
        length = len(r)-1
        fragments = []
        if d.is_linear():
            #
            #   START of the sequence to FIRST site.
            #
            fragments.append(d[1:r[0]])
            if length:
                #
                #   if more than one site add them.
                #
                fragments += [d[r[x]:r[x+1]] for x in range(length)]
            #
            #   LAST site to END of the sequence.
            #
            fragments.append(d[r[-1]:])
        else:
            #
            #   circular : bridge LAST site to FIRST site.
            #
            fragments.append(d[r[-1]:]+d[1:r[0]])
            if not length:
                #
                #   one site we finish here.
                #
                return tuple(fragments)
            #
            #   add the others.
            #
            fragments += [d[r[x]:r[x+1]] for x in range(length)]
        return tuple(fragments)
    catalyze = catalyse

    @classmethod
    def is_blunt(self):
        """RE.is_blunt() -> bool.

        True if the enzyme produces blunt end.

        see also:
            RE.is_3overhang()
            RE.is_5overhang()
            RE.is_unknown()"""
        return False

    @classmethod
    def is_5overhang(self):
        """RE.is_5overhang() -> bool.

        True if the enzyme produces 5' overhang sticky end.

        see also:
            RE.is_3overhang()
            RE.is_blunt()
            RE.is_unknown()"""
        return True

    @classmethod
    def is_3overhang(self):
        """RE.is_3overhang() -> bool.

        True if the enzyme produces 3' overhang sticky end.

        see also:
            RE.is_5overhang()
            RE.is_blunt()
            RE.is_unknown()"""
        return False

    @classmethod
    def overhang(self):
        """RE.overhang() -> str. type of overhang of the enzyme.,

        can be "3' overhang", "5' overhang", "blunt", "unknown"   """
        return "5' overhang"

    @classmethod
    def compatible_end(self, batch=None):
        """RE.compatible_end() -> list.

        list of all the enzymes that share compatible end with RE."""
        if not batch:
            batch = AllEnzymes
        r = sorted(x for x in iter(AllEnzymes) if x.is_5overhang() and x % self)
        return r

    @classmethod
    def _mod1(self, other):
        """RE._mod1(other) -> bool.

        for internal use only

        test for the compatibility of restriction ending of RE and other."""
        if issubclass(other, Ov5):
            return self._mod2(other)
        else:
            return False


class Ov3(AbstractCut):
    """Implement the methods specific to the enzymes for which the overhang
    is recessed in 5'.

    The enzyme cuts the - strand after the + strand of the DNA.

    Internal use only. Not meant to be instantiated."""

    @classmethod
    def catalyse(self, dna, linear=True):
        """RE.catalyse(dna, linear=True) -> tuple of DNA.
        RE.catalyze(dna, linear=True) -> tuple of DNA.

        return a tuple of dna as will be produced by using RE to restrict the
        dna.

        dna must be a Bio.Seq.Seq instance or a Bio.Seq.MutableSeq instance.

        if linear is False, the sequence is considered to be circular and the
        output will be modified accordingly."""
        r = self.search(dna, linear)
        d = self.dna
        if not r:
            return d[1:],
        fragments = []
        length = len(r)-1
        if d.is_linear():
            #
            #   START of the sequence to FIRST site.
            #
            fragments.append(d[1:r[0]])
            if length:
                #
                #   if more than one site add them.
                #
                fragments += [d[r[x]:r[x+1]] for x in range(length)]
            #
            #   LAST site to END of the sequence.
            #
            fragments.append(d[r[-1]:])
        else:
            #
            #   circular : bridge LAST site to FIRST site.
            #
            fragments.append(d[r[-1]:]+d[1:r[0]])
            if not length:
                #
                #   one site we finish here.
                #
                return tuple(fragments)
            #
            #   add the others.
            #
            fragments += [d[r[x]:r[x+1]] for x in range(length)]
        return tuple(fragments)
    catalyze = catalyse

    @classmethod
    def is_blunt(self):
        """RE.is_blunt() -> bool.

        True if the enzyme produces blunt end.

        see also:
            RE.is_3overhang()
            RE.is_5overhang()
            RE.is_unknown()"""
        return False

    @classmethod
    def is_5overhang(self):
        """RE.is_5overhang() -> bool.

        True if the enzyme produces 5' overhang sticky end.

        see also:
            RE.is_3overhang()
            RE.is_blunt()
            RE.is_unknown()"""
        return False

    @classmethod
    def is_3overhang(self):
        """RE.is_3overhang() -> bool.

        True if the enzyme produces 3' overhang sticky end.

        see also:
            RE.is_5overhang()
            RE.is_blunt()
            RE.is_unknown()"""
        return True

    @classmethod
    def overhang(self):
        """RE.overhang() -> str. type of overhang of the enzyme.,

        can be "3' overhang", "5' overhang", "blunt", "unknown"   """
        return "3' overhang"

    @classmethod
    def compatible_end(self, batch=None):
        """RE.compatible_end() -> list.

        list of all the enzymes that share compatible end with RE."""
        if not batch:
            batch = AllEnzymes
        r = sorted(x for x in iter(AllEnzymes) if x.is_3overhang() and x % self)
        return r

    @classmethod
    def _mod1(self, other):
        """RE._mod1(other) -> bool.

        for internal use only

        test for the compatibility of restriction ending of RE and other."""
        #
        #   called by RE._mod1(other) when the one of the enzyme is ambiguous
        #
        if issubclass(other, Ov3):
            return self._mod2(other)
        else:
            return False


class Defined(AbstractCut):
    """Implement the methods specific to the enzymes for which the overhang
    and the cut are not variable.

    Typical example : EcoRI -> G^AATT_C
                      The overhang will always be AATT
    Notes:
        Blunt enzymes are always defined. even if there site is GGATCCNNN^_N
        There overhang is always the same : blunt!

    Internal use only. Not meant to be instantiated."""

    @classmethod
    def _drop(self):
        """RE._drop() -> list.

        for internal use only.

        drop the site that are situated outside the sequence in linear sequence.
        modify the index for site in circular sequences."""
        #
        #   remove or modify the results that are outside the sequence.
        #   This is necessary since after finding the site we add the distance
        #   from the site to the cut with the _modify and _rev_modify methods.
        #   For linear we will remove these sites altogether.
        #   For circular sequence, we modify the result rather than _drop it
        #   since the site is in the sequence.
        #
        length = len(self.dna)
        drop = itertools.dropwhile
        take = itertools.takewhile
        if self.dna.is_linear():
            self.results = [x for x in drop(lambda x:x<1, self.results)]
            self.results = [x for x in take(lambda x:x<length, self.results)]
        else:
            for index, location in enumerate(self.results):
                if location < 1:
                    self.results[index] += length
                else:
                    break
            for index, location in enumerate(self.results[::-1]):
                if location > length:
                    self.results[-(index+1)] -= length
                else:
                    break
        return

    @classmethod
    def is_defined(self):
        """RE.is_defined() -> bool.

        True if the sequence recognised and cut is constant,
        i.e. the recognition site is not degenerated AND the enzyme cut inside
        the site.

        see also:
            RE.is_ambiguous()
            RE.is_unknown()"""
        return True

    @classmethod
    def is_ambiguous(self):
        """RE.is_ambiguous() -> bool.

        True if the sequence recognised and cut is ambiguous,
        i.e. the recognition site is degenerated AND/OR the enzyme cut outside
        the site.

        see also:
            RE.is_defined()
            RE.is_unknown()"""
        return False

    @classmethod
    def is_unknown(self):
        """RE.is_unknown() -> bool.

        True if the sequence is unknown,
        i.e. the recognition site has not been characterised yet.

        see also:
            RE.is_defined()
            RE.is_ambiguous()"""
        return False

    @classmethod
    def elucidate(self):
        """RE.elucidate() -> str

        return a representation of the site with the cut on the (+) strand
        represented as '^' and the cut on the (-) strand as '_'.
        ie:
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
        f5 = self.fst5
        f3 = self.fst3
        site = self.site
        if self.cut_twice():
            re = 'cut twice, not yet implemented sorry.'
        elif self.is_5overhang():
            if f5 == f3 == 0:
                re = 'N^'+ self.site + '_N'
            elif f3 == 0:
                re = site[:f5] + '^' + site[f5:] + '_N'
            else:
                re = site[:f5] + '^' + site[f5:f3] + '_' + site[f3:]
        elif self.is_blunt():
            re = site[:f5] + '^_' + site[f5:]
        else:
            if f5 == f3 == 0:
                re = 'N_'+ site + '^N'
            else:
                re = site[:f3] + '_' + site[f3:f5] +'^'+ site[f5:]
        return re

    @classmethod
    def _mod2(self, other):
        """RE._mod2(other) -> bool.

        for internal use only

        test for the compatibility of restriction ending of RE and other."""
        #
        #   called by RE._mod1(other) when the one of the enzyme is ambiguous
        #
        if other.ovhgseq == self.ovhgseq:
            return True
        elif issubclass(other, Ambiguous):
            return other._mod2(self)
        else:
            return False


class Ambiguous(AbstractCut):
    """Implement the methods specific to the enzymes for which the overhang
    is variable.

    Typical example : BstXI -> CCAN_NNNN^NTGG
                      The overhang can be any sequence of 4 bases.
    Notes:
        Blunt enzymes are always defined. even if there site is GGATCCNNN^_N
        There overhang is always the same : blunt!

    Internal use only. Not meant to be instantiated."""

    @classmethod
    def _drop(self):
        """RE._drop() -> list.

        for internal use only.

        drop the site that are situated outside the sequence in linear sequence.
        modify the index for site in circular sequences."""
        length = len(self.dna)
        drop = itertools.dropwhile
        take = itertools.takewhile
        if self.dna.is_linear():
            self.results = [x for x in drop(lambda x: x < 1, self.results)]
            self.results = [x for x in take(lambda x: x <length, self.results)]
        else:
            for index, location in enumerate(self.results):
                if location < 1:
                    self.results[index] += length
                else:
                    break
            for index, location in enumerate(self.results[::-1]):
                if location > length:
                    self.results[-(index+1)] -= length
                else:
                    break
        return

    @classmethod
    def is_defined(self):
        """RE.is_defined() -> bool.

        True if the sequence recognised and cut is constant,
        i.e. the recognition site is not degenerated AND the enzyme cut inside
        the site.

        see also:
            RE.is_ambiguous()
            RE.is_unknown()"""
        return False

    @classmethod
    def is_ambiguous(self):
        """RE.is_ambiguous() -> bool.

        True if the sequence recognised and cut is ambiguous,
        i.e. the recognition site is degenerated AND/OR the enzyme cut outside
        the site.

        see also:
            RE.is_defined()
            RE.is_unknown()"""
        return True

    @classmethod
    def is_unknown(self):
        """RE.is_unknown() -> bool.

        True if the sequence is unknown,
        i.e. the recognition site has not been characterised yet.

        see also:
            RE.is_defined()
            RE.is_ambiguous()"""
        return False

    @classmethod
    def _mod2(self, other):
        """RE._mod2(other) -> bool.

        for internal use only

        test for the compatibility of restriction ending of RE and other."""
        #
        #   called by RE._mod1(other) when the one of the enzyme is ambiguous
        #
        if len(self.ovhgseq) != len(other.ovhgseq):
            return False
        else:
            se = self.ovhgseq
            for base in se:
                if base in 'ATCG':
                    pass
                if base in 'N':
                    se = '.'.join(se.split('N'))
                if base in 'RYWMSKHDBV':
                    expand = '['+ matching[base] + ']'
                    se = expand.join(se.split(base))
            if re.match(se, other.ovhgseq):
                return True
            else:
                return False

    @classmethod
    def elucidate(self):
        """RE.elucidate() -> str

        return a representation of the site with the cut on the (+) strand
        represented as '^' and the cut on the (-) strand as '_'.
        ie:
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
        f5 = self.fst5
        f3 = self.fst3
        length = len(self)
        site = self.site
        if self.cut_twice():
            re = 'cut twice, not yet implemented sorry.'
        elif self.is_5overhang():
            if f3 == f5 == 0:
                re = 'N^' + site +'_N'
            elif 0 <= f5 <= length and 0 <= f3+length <= length:
                re = site[:f5] + '^' + site[f5:f3] + '_' + site[f3:]
            elif 0 <= f5 <= length:
                re = site[:f5] + '^' + site[f5:] + f3*'N' + '_N'
            elif 0 <= f3+length <= length:
                re = 'N^' + abs(f5) * 'N' + site[:f3] + '_' + site[f3:]
            elif f3+length < 0:
                re = 'N^'*abs(f5)*'N' + '_' + abs(length+f3)*'N' + site
            elif f5 > length:
                re = site + (f5-length)*'N'+'^'+(length+f3-f5)*'N'+'_N'
            else:
                re = 'N^' + abs(f5) * 'N' + site + f3*'N' + '_N'
        elif self.is_blunt():
            if f5 < 0:
                re = 'N^_' + abs(f5)*'N' + site
            elif f5 > length:
                re = site + (f5-length)*'N' + '^_N'
            else:
                raise ValueError('%s.easyrepr() : error f5=%i'
                                 % (self.name, f5))
        else:
            if f3 == 0:
                if f5 == 0:
                    re = 'N_' + site + '^N'
                else:
                    re = site + '_' + (f5-length)*'N' + '^N'
            elif 0 < f3+length <= length and 0 <= f5 <= length:
                re = site[:f3] + '_' + site[f3:f5] + '^' + site[f5:]
            elif 0 < f3+length <= length:
                re = site[:f3] + '_' + site[f3:] + (f5-length)*'N' + '^N'
            elif 0 <= f5 <= length:
                re = 'N_' +'N'*(f3+length) + site[:f5] + '^' + site[f5:]
            elif f3 > 0:
                re = site + f3*'N' + '_' + (f5-f3-length)*'N' + '^N'
            elif f5 < 0:
                re = 'N_' + abs(f3-f5+length)*'N' + '^' + abs(f5)*'N' + site
            else:
                re = 'N_' + abs(f3+length)*'N' + site + (f5-length)*'N' + '^N'
        return re


class NotDefined(AbstractCut):
    """Implement the methods specific to the enzymes for which the overhang
    is not characterised.

    Correspond to NoCut and Unknown.

    Internal use only. Not meant to be instantiated."""

    @classmethod
    def _drop(self):
        """RE._drop() -> list.

        for internal use only.

        drop the site that are situated outside the sequence in linear sequence.
        modify the index for site in circular sequences."""
        if self.dna.is_linear():
            return
        else:
            length = len(self.dna)
            for index, location in enumerate(self.results):
                if location < 1:
                    self.results[index] += length
                else:
                    break
            for index, location in enumerate(self.results[:-1]):
                if location > length:
                    self.results[-(index+1)] -= length
                else:
                    break
        return

    @classmethod
    def is_defined(self):
        """RE.is_defined() -> bool.

        True if the sequence recognised and cut is constant,
        i.e. the recognition site is not degenerated AND the enzyme cut inside
        the site.

        see also:
            RE.is_ambiguous()
            RE.is_unknown()"""
        return False

    @classmethod
    def is_ambiguous(self):
        """RE.is_ambiguous() -> bool.

        True if the sequence recognised and cut is ambiguous,
        i.e. the recognition site is degenerated AND/OR the enzyme cut outside
        the site.

        see also:
            RE.is_defined()
            RE.is_unknown()"""
        return False

    @classmethod
    def is_unknown(self):
        """RE.is_unknown() -> bool.

        True if the sequence is unknown,
        i.e. the recognition site has not been characterised yet.

        see also:
            RE.is_defined()
            RE.is_ambiguous()"""
        return True

    @classmethod
    def _mod2(self, other):
        """RE._mod2(other) -> bool.

        for internal use only

        test for the compatibility of restriction ending of RE and other."""
        #
        #   Normally we should not arrive here. But well better safe than sorry.
        #   the overhang is not defined we are compatible with nobody.
        #   could raise an Error may be rather than return quietly.
        #
        # return False
        raise ValueError("%s.mod2(%s), %s : NotDefined. pas glop pas glop!"
                         % (str(self), str(other), str(self)))

    @classmethod
    def elucidate(self):
        """RE.elucidate() -> str

        return a representation of the site with the cut on the (+) strand
        represented as '^' and the cut on the (-) strand as '_'.
        ie:
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
        return '? %s ?' % self.site


class Commercially_available(AbstractCut):
    #
    #   Recent addition to Rebase make this naming convention uncertain.
    #   May be better to says enzymes which have a supplier.
    #
    """Implement the methods specific to the enzymes which are commercially
    available.

    Internal use only. Not meant to be instantiated."""

    @classmethod
    def suppliers(self):
        """RE.suppliers() -> print the suppliers of RE."""
        for s in self.suppliers_dict():
            print(s + ',')
        return

    @classmethod
    def supplier_list(self):
        """RE.supplier_list() -> list.

        list of the supplier names for RE."""
        return [v[0] for k, v in suppliers_dict.items() if k in self.suppl]

    @classmethod
    def buffers(self, supplier):
        """RE.buffers(supplier) -> string.

        not implemented yet."""
        return

    @classmethod
    def is_comm(self):
        """RE.iscomm() -> bool.

        True if RE has suppliers."""
        return True


class Not_available(AbstractCut):
    """Implement the methods specific to the enzymes which are not commercially
    available.

    Internal use only. Not meant to be instantiated."""

    @staticmethod
    def suppliers():
        """RE.suppliers() -> print the suppliers of RE."""
        return None

    @classmethod
    def supplier_list(self):
        """RE.supplier_list() -> list.

        list of the supplier names for RE."""
        return []

    @classmethod
    def buffers(self, supplier):
        """RE.buffers(supplier) -> string.

        not implemented yet."""
        raise TypeError("Enzyme not commercially available.")

    @classmethod
    def is_comm(self):
        """RE.iscomm() -> bool.

        True if RE has suppliers."""
        return False


###############################################################################
#                                                                             #
#                       Restriction Batch                                     #
#                                                                             #
###############################################################################


class RestrictionBatch(set):

    def __init__(self, first=[], suppliers=[]):
        """RestrictionBatch([sequence]) -> new RestrictionBatch."""
        first = [self.format(x) for x in first]
        first += [eval(x) for n in suppliers for x in suppliers_dict[n][1]]
        set.__init__(self, first)
        self.mapping = dict.fromkeys(self)
        self.already_mapped = None

    def __str__(self):
        if len(self) < 5:
            return '+'.join(self.elements())
        else:
            return '...'.join(('+'.join(self.elements()[:2]),
                               '+'.join(self.elements()[-2:])))

    def __repr__(self):
        return 'RestrictionBatch(%s)' % self.elements()

    def __contains__(self, other):
        try:
            other = self.format(other)
        except ValueError: # other is not a restriction enzyme
            return False
        return set.__contains__(self, other)

    def __div__(self, other):
        return self.search(other)

    def __rdiv__(self, other):
        return self.search(other)

    def get(self, enzyme, add=False):
        """B.get(enzyme[, add]) -> enzyme class.

        if add is True and enzyme is not in B add enzyme to B.
        if add is False (which is the default) only return enzyme.
        if enzyme is not a RestrictionType or can not be evaluated to
        a RestrictionType, raise a ValueError."""
        e = self.format(enzyme)
        if e in self:
            return e
        elif add:
            self.add(e)
            return e
        else:
            raise ValueError('enzyme %s is not in RestrictionBatch'
                             % e.__name__)

    def lambdasplit(self, func):
        """B.lambdasplit(func) -> RestrictionBatch .

        the new batch will contains only the enzymes for which
        func return True."""
        d = [x for x in filter(func, self)]
        new = RestrictionBatch()
        new._data = dict(zip(d, [True]*len(d)))
        return new

    def add_supplier(self, letter):
        """B.add_supplier(letter) -> add a new set of enzyme to B.

        letter represents the suppliers as defined in the dictionary
        RestrictionDictionary.suppliers
        return None.
        raise a KeyError if letter is not a supplier code."""
        supplier = suppliers_dict[letter]
        self.suppliers.append(letter)
        for x in supplier[1]:
            self.add_nocheck(eval(x))
        return

    def current_suppliers(self):
        """B.current_suppliers() -> add a new set of enzyme to B.

        return a sorted list of the suppliers which have been used to
        create the batch."""
        suppl_list = sorted(suppliers_dict[x][0] for x in self.suppliers)
        return suppl_list

    def __iadd__(self, other):
        """ b += other -> add other to b, check the type of other."""
        self.add(other)
        return self

    def __add__(self, other):
        """ b + other -> new RestrictionBatch."""
        new = self.__class__(self)
        new.add(other)
        return new

    def remove(self, other):
        """B.remove(other) -> remove other from B if other is a RestrictionType.

        Safe set.remove method. Verify that other is a RestrictionType or can be
        evaluated to a RestrictionType.
        raise a ValueError if other can not be evaluated to a RestrictionType.
        raise a KeyError if other is not in B."""
        return set.remove(self, self.format(other))

    def add(self, other):
        """B.add(other) -> add other to B if other is a RestrictionType.

        Safe set.add method. Verify that other is a RestrictionType or can be
        evaluated to a RestrictionType.
        raise a ValueError if other can not be evaluated to a RestrictionType.
        """
        return set.add(self, self.format(other))

    def add_nocheck(self, other):
        """B.add_nocheck(other) -> add other to B. don't check type of other.
        """
        return set.add(self, other)

    def format(self, y):
        """B.format(y) -> RestrictionType or raise ValueError.

        if y is a RestrictionType return y
        if y can be evaluated to a RestrictionType return eval(y)
        raise a Value Error in all other case."""
        try:
            if isinstance(y, RestrictionType):
                return y
            elif isinstance(eval(str(y)), RestrictionType):
                return eval(y)
            else:
                pass
        except (NameError, SyntaxError):
            pass
        raise ValueError('%s is not a RestrictionType' % y.__class__)

    def is_restriction(self, y):
        """B.is_restriction(y) -> bool.

        True is y or eval(y) is a RestrictionType."""
        return isinstance(y, RestrictionType) or \
               isinstance(eval(str(y)), RestrictionType)

    def split(self, *classes, **bool):
        """B.split(class, [class.__name__ = True]) -> new RestrictionBatch.

        it works but it is slow, so it has really an interest when splitting
        over multiple conditions."""
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
        d = [k for k in filter(splittest, self)]
        new = RestrictionBatch()
        new._data = dict(zip(d, [True]*len(d)))
        return new

    def elements(self):
        """B.elements() -> tuple.

        give all the names of the enzymes in B sorted alphabetically."""
        l = sorted(str(e) for e in self)
        return l

    def as_string(self):
        """B.as_string() -> list.

        return a list of the name of the elements of B."""
        return [str(e) for e in self]

    @classmethod
    def suppl_codes(self):
        """B.suppl_codes() -> dict

        letter code for the suppliers"""
        supply = dict((k, v[0]) for k, v in suppliers_dict.items())
        return supply

    @classmethod
    def show_codes(self):
        """B.show_codes() -> letter codes for the suppliers"""
        supply = [' = '.join(i) for i in self.suppl_codes().items()]
        print('\n'.join(supply))
        return

    def search(self, dna, linear=True):
        """B.search(dna) -> dict."""
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
                self.mapping = dict((x, x.search(fseq)) for x in self)
                return self.mapping
        elif isinstance(dna, FormattedSeq):
            if (str(dna), dna.linear) == self.already_mapped:
                return self.mapping
            else:
                self.already_mapped = str(dna), dna.linear
                self.mapping = dict((x, x.search(dna)) for x in self)
                return self.mapping
        raise TypeError("Expected Seq or MutableSeq instance, got %s instead"
                        %type(dna))

###############################################################################
#                                                                             #
#                       Restriction Analysis                                  #
#                                                                             #
###############################################################################


class Analysis(RestrictionBatch, PrintFormat):

    def __init__(self, restrictionbatch=RestrictionBatch(), sequence=DNA(''),
                 linear=True):
        """Analysis([restrictionbatch [, sequence] linear=True]) -> New Analysis class.

        For most of the method of this class if a dictionary is given it will
        be used as the base to calculate the results.
        If no dictionary is given a new analysis using the Restriction Batch
        which has been given when the Analysis class has been instantiated."""
        RestrictionBatch.__init__(self, restrictionbatch)
        self.rb = restrictionbatch
        self.sequence = sequence
        self.linear = linear
        if self.sequence:
            self.search(self.sequence, self.linear)

    def __repr__(self):
        return 'Analysis(%s,%s,%s)'%\
               (repr(self.rb), repr(self.sequence), self.linear)

    def _sub_set(self, wanted):
        """A._sub_set(other_set) -> dict.

        Internal use only.

        screen the results through wanted set.
        Keep only the results for which the enzymes is in wanted set.
        """
        return dict((k, v) for k, v in self.mapping.items() if k in wanted)

    def _boundaries(self, start, end):
        """A._boundaries(start, end) -> tuple.

        Format the boundaries for use with the methods that limit the
        search to only part of the sequence given to analyse.
        """
        if not isinstance(start, int):
            raise TypeError('expected int, got %s instead' % type(start))
        if not isinstance(end, int):
            raise TypeError('expected int, got %s instead' % type(end))
        if start < 1:
            start += len(self.sequence)
        if end < 1:
            end += len(self.sequence)
        if start < end:
            pass
        else:
            start, end == end, start
        if start < 1:
            start == 1
        if start < end:
            return start, end, self._test_normal
        else:
            return start, end, self._test_reverse

    def _test_normal(self, start, end, site):
        """A._test_normal(start, end, site) -> bool.

        Internal use only
        Test if site is in between start and end.
        """
        return start <= site < end

    def _test_reverse(self, start, end, site):
        """A._test_reverse(start, end, site) -> bool.

        Internal use only
        Test if site is in between end and start (for circular sequences).
        """
        return start <= site <= len(self.sequence) or 1 <= site < end

    def print_that(self, dct=None, title='', s1=''):
        """A.print_that([dct[, title[, s1]]]) -> print the results from dct.

        If dct is not given the full dictionary is used.
        """
        if not dct:
            dct = self.mapping
        print("")
        return PrintFormat.print_that(self, dct, title, s1)

    def change(self, **what):
        """A.change(**attribute_name) -> Change attribute of Analysis.

        It is possible to change the width of the shell by setting
        self.ConsoleWidth to what you want.
        self.NameWidth refer to the maximal length of the enzyme name.

        Changing one of these parameters here might not give the results
        you expect. In which case, you can settle back to a 80 columns shell
        or try to change self.Cmodulo and self.PrefWidth in PrintFormat until
        you get it right."""
        for k, v in what.items():
            if k in ('NameWidth', 'ConsoleWidth'):
                setattr(self, k, v)
                self.Cmodulo = self.ConsoleWidth % self.NameWidth
                self.PrefWidth = self.ConsoleWidth - self.Cmodulo
            elif k is 'sequence':
                setattr(self, 'sequence', v)
                self.search(self.sequence, self.linear)
            elif k is 'rb':
                self = Analysis.__init__(self, v, self.sequence, self.linear)
            elif k is 'linear':
                setattr(self, 'linear', v)
                self.search(self.sequence, v)
            elif k in ('Indent', 'Maxsize'):
                setattr(self, k, v)
            elif k in ('Cmodulo', 'PrefWidth'):
                raise AttributeError(
                    'To change %s, change NameWidth and/or ConsoleWidth'
                    % name)
            else:
                raise AttributeError(
                    'Analysis has no attribute %s' % name)
        return

    def full(self, linear=True):
        """A.full() -> dict.

        Full Restriction Map of the sequence."""
        return self.mapping

    def blunt(self, dct=None):
        """A.blunt([dct]) -> dict.

        Only the enzymes which have a 3'overhang restriction site."""
        if not dct:
            dct = self.mapping
        return dict((k, v) for k, v in dct.items() if k.is_blunt())

    def overhang5(self, dct=None):
        """A.overhang5([dct]) -> dict.

        Only the enzymes which have a 5' overhang restriction site."""
        if not dct:
            dct = self.mapping
        return dict((k, v) for k, v in dct.items() if k.is_5overhang())

    def overhang3(self, dct=None):
        """A.Overhang3([dct]) -> dict.

        Only the enzymes which have a 3'overhang restriction site."""
        if not dct:
            dct = self.mapping
        return dict((k, v) for k, v in dct.items() if k.is_3overhang())

    def defined(self, dct=None):
        """A.defined([dct]) -> dict.

        Only the enzymes that have a defined restriction site in Rebase."""
        if not dct:
            dct = self.mapping
        return dict((k, v) for k, v in dct.items() if k.is_defined())

    def with_sites(self, dct=None):
        """A.with_sites([dct]) -> dict.

        Enzymes which have at least one site in the sequence."""
        if not dct:
            dct = self.mapping
        return dict((k, v) for k, v in dct.items() if v)

    def without_site(self, dct=None):
        """A.without_site([dct]) -> dict.

        Enzymes which have no site in the sequence."""
        if not dct:
            dct = self.mapping
        return dict((k, v) for k, v in dct.items() if not v)

    def with_N_sites(self, N, dct=None):
        """A.With_N_Sites(N [, dct]) -> dict.

        Enzymes which cut N times the sequence."""
        if not dct:
            dct = self.mapping
        return dict((k, v) for k, v in dct.items()if len(v) == N)

    def with_number_list(self, list, dct=None):
        if not dct:
            dct = self.mapping
        return dict((k, v) for k, v in dct.items() if len(v) in list)

    def with_name(self, names, dct=None):
        """A.with_name(list_of_names [, dct]) ->

         Limit the search to the enzymes named in list_of_names."""
        for i, enzyme in enumerate(names):
            if enzyme not in AllEnzymes:
                print("no data for the enzyme: %s" % name)
                del names[i]
        if not dct:
            return RestrictionBatch(names).search(self.sequence)
        return dict((n, dct[n]) for n in names if n in dct)

    def with_site_size(self, site_size, dct=None):
        """A.with_site_size(site_size [, dct]) ->

         Limit the search to the enzymes whose site is of size <site_size>."""
        sites = [name for name in self if name.size == site_size]
        if not dct:
            return RestrictionBatch(sites).search(self.sequence)
        return dict((k, v) for k, v in dct.items() if k in site_size)

    def only_between(self, start, end, dct=None):
        """A.only_between(start, end[, dct]) -> dict.

        Enzymes that cut the sequence only in between start and end."""
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
        """A.between(start, end [, dct]) -> dict.

        Enzymes that cut the sequence at least in between start and end.
        They may cut outside as well."""
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
        """A.show_only_between(start, end [, dct]) -> dict.

        Enzymes that cut the sequence outside of the region
        in between start and end but do not cut inside."""
        d = []
        if start <= end:
            d = [(k, [vv for vv in v if start<=vv<=end])
                 for v in self.between(start, end, dct)]
        else:
            d = [(k, [vv for vv in v if start<=vv or vv <= end])
                 for v in self.between(start, end, dct)]
        return dict(d)

    def only_outside(self, start, end, dct=None):
        """A.only_outside(start, end [, dct]) -> dict.

        Enzymes that cut the sequence outside of the region
        in between start and end but do not cut inside."""
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
        """A.outside((start, end [, dct]) -> dict.

        Enzymes that cut outside the region in between start and end.
        No test is made to know if they cut or not inside this region."""
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
        """A.do_not_cut(start, end [, dct]) -> dict.

        Enzymes that do not cut the region in between start and end."""
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
#   It is essential to run Restriction with doc string optimisation (-OO switch)
#   as the doc string of 660 classes take a lot of processing.
#
CommOnly = RestrictionBatch()    # commercial enzymes
NonComm = RestrictionBatch()     # not available commercially
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
    T = type.__new__(RestrictionType, 'RestrictionType', bases, {})
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
AllEnzymes = CommOnly | NonComm
#
#   Now, place the enzymes in locals so they can be imported.
#
names = [str(x) for x in AllEnzymes]
try:
    del x
except NameError:
    # Scoping changed in Python 3, the variable isn't leaked
    pass
locals().update(dict(zip(names, AllEnzymes)))
__all__=['FormattedSeq', 'Analysis', 'RestrictionBatch', 'AllEnzymes', 'CommOnly', 'NonComm']+names
del k, enzymes, TYPE, bases, names
