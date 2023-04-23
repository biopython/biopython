# Copyright 2000 Andrew Dalke.
# Copyright 2000-2002 Brad Chapman.
# Copyright 2004-2005, 2010 by M de Hoon.
# Copyright 2007-2020 by Peter Cock.
# All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Provide objects to represent biological sequences.

See also the Seq_ wiki and the chapter in our tutorial:
 - `HTML Tutorial`_
 - `PDF Tutorial`_

.. _Seq: http://biopython.org/wiki/Seq
.. _`HTML Tutorial`: http://biopython.org/DIST/docs/tutorial/Tutorial.html
.. _`PDF Tutorial`: http://biopython.org/DIST/docs/tutorial/Tutorial.pdf

"""
import array
import numbers
import warnings

from abc import ABC
from abc import abstractmethod

from Bio import BiopythonDeprecationWarning
from Bio import BiopythonWarning
from Bio.Data import CodonTable
from Bio.Data import IUPACData


def _maketrans(complement_mapping):
    """Make a python string translation table (PRIVATE).

    Arguments:
     - complement_mapping - a dictionary such as ambiguous_dna_complement
       and ambiguous_rna_complement from Data.IUPACData.

    Returns a translation table (a bytes object of length 256) for use with
    the python string's translate method to use in a (reverse) complement.

    Compatible with lower case and upper case sequences.

    For internal use only.
    """
    keys = "".join(complement_mapping.keys()).encode("ASCII")
    values = "".join(complement_mapping.values()).encode("ASCII")
    return bytes.maketrans(keys + keys.lower(), values + values.lower())


ambiguous_dna_complement = dict(IUPACData.ambiguous_dna_complement)
ambiguous_dna_complement["U"] = ambiguous_dna_complement["T"]
_dna_complement_table = _maketrans(ambiguous_dna_complement)
del ambiguous_dna_complement
ambiguous_rna_complement = dict(IUPACData.ambiguous_rna_complement)
ambiguous_rna_complement["T"] = ambiguous_rna_complement["U"]
_rna_complement_table = _maketrans(ambiguous_rna_complement)
del ambiguous_rna_complement


class SequenceDataAbstractBaseClass(ABC):
    """Abstract base class for sequence content providers.

    Most users will not need to use this class. It is used internally as a base
    class for sequence content provider classes such as _UndefinedSequenceData
    defined in this module, and _TwoBitSequenceData in Bio.SeqIO.TwoBitIO.
    Instances of these classes can be used instead of a ``bytes`` object as the
    data argument when creating a Seq object, and provide the sequence content
    only when requested via ``__getitem__``. This allows lazy parsers to load
    and parse sequence data from a file only for the requested sequence regions,
    and _UndefinedSequenceData instances to raise an exception when undefined
    sequence data are requested.

    Future implementations of lazy parsers that similarly provide on-demand
    parsing of sequence data should use a subclass of this abstract class and
    implement the abstract methods ``__len__`` and ``__getitem__``:

    * ``__len__`` must return the sequence length;
    * ``__getitem__`` must return

      * a ``bytes`` object for the requested region; or
      * a new instance of the subclass for the requested region; or
      * raise an ``UndefinedSequenceError``.

      Calling ``__getitem__`` for a sequence region of size zero should always
      return an empty ``bytes`` object.
      Calling ``__getitem__`` for the full sequence (as in data[:]) should
      either return a ``bytes`` object with the full sequence, or raise an
      ``UndefinedSequenceError``.

    Subclasses of SequenceDataAbstractBaseClass must call ``super().__init__()``
    as part of their ``__init__`` method.
    """

    __slots__ = ()

    def __init__(self):
        """Check if ``__getitem__`` returns a bytes-like object."""
        assert self[:0] == b""

    @abstractmethod
    def __len__(self):
        pass

    @abstractmethod
    def __getitem__(self, key):
        pass

    def __bytes__(self):
        return self[:]

    def __hash__(self):
        return hash(bytes(self))

    def __eq__(self, other):
        return bytes(self) == other

    def __lt__(self, other):
        return bytes(self) < other

    def __le__(self, other):
        return bytes(self) <= other

    def __gt__(self, other):
        return bytes(self) > other

    def __ge__(self, other):
        return bytes(self) >= other

    def __add__(self, other):
        try:
            return bytes(self) + bytes(other)
        except UndefinedSequenceError:
            return NotImplemented
            # will be handled by _UndefinedSequenceData.__radd__ or
            # by _PartiallyDefinedSequenceData.__radd__

    def __radd__(self, other):
        return other + bytes(self)

    def __mul__(self, other):
        return other * bytes(self)

    def __contains__(self, item):
        return bytes(self).__contains__(item)

    def decode(self, encoding="utf-8"):
        """Decode the data as bytes using the codec registered for encoding.

        encoding
          The encoding with which to decode the bytes.
        """
        return bytes(self).decode(encoding)

    def count(self, sub, start=None, end=None):
        """Return the number of non-overlapping occurrences of sub in data[start:end].

        Optional arguments start and end are interpreted as in slice notation.
        This method behaves as the count method of Python strings.
        """
        return bytes(self).count(sub, start, end)

    def find(self, sub, start=None, end=None):
        """Return the lowest index in data where subsection sub is found.

        Return the lowest index in data where subsection sub is found,
        such that sub is contained within data[start,end].  Optional
        arguments start and end are interpreted as in slice notation.

        Return -1 on failure.
        """
        return bytes(self).find(sub, start, end)

    def rfind(self, sub, start=None, end=None):
        """Return the highest index in data where subsection sub is found.

        Return the highest index in data where subsection sub is found,
        such that sub is contained within data[start,end].  Optional
        arguments start and end are interpreted as in slice notation.

        Return -1 on failure.
        """
        return bytes(self).rfind(sub, start, end)

    def index(self, sub, start=None, end=None):
        """Return the lowest index in data where subsection sub is found.

        Return the lowest index in data where subsection sub is found,
        such that sub is contained within data[start,end].  Optional
        arguments start and end are interpreted as in slice notation.

        Raises ValueError when the subsection is not found.
        """
        return bytes(self).index(sub, start, end)

    def rindex(self, sub, start=None, end=None):
        """Return the highest index in data where subsection sub is found.

        Return the highest index in data where subsection sub is found,
        such that sub is contained within data[start,end].  Optional
        arguments start and end are interpreted as in slice notation.

        Raise ValueError when the subsection is not found.
        """
        return bytes(self).rindex(sub, start, end)

    def startswith(self, prefix, start=None, end=None):
        """Return True if data starts with the specified prefix, False otherwise.

        With optional start, test data beginning at that position.
        With optional end, stop comparing data at that position.
        prefix can also be a tuple of bytes to try.
        """
        return bytes(self).startswith(prefix, start, end)

    def endswith(self, suffix, start=None, end=None):
        """Return True if data ends with the specified suffix, False otherwise.

        With optional start, test data beginning at that position.
        With optional end, stop comparing data at that position.
        suffix can also be a tuple of bytes to try.
        """
        return bytes(self).endswith(suffix, start, end)

    def split(self, sep=None, maxsplit=-1):
        """Return a list of the sections in the data, using sep as the delimiter.

        sep
          The delimiter according which to split the data.
          None (the default value) means split on ASCII whitespace characters
          (space, tab, return, newline, formfeed, vertical tab).
        maxsplit
          Maximum number of splits to do.
          -1 (the default value) means no limit.
        """
        return bytes(self).split(sep, maxsplit)

    def rsplit(self, sep=None, maxsplit=-1):
        """Return a list of the sections in the data, using sep as the delimiter.

        sep
          The delimiter according which to split the data.
          None (the default value) means split on ASCII whitespace characters
          (space, tab, return, newline, formfeed, vertical tab).
        maxsplit
          Maximum number of splits to do.
          -1 (the default value) means no limit.

        Splitting is done starting at the end of the data and working to the front.
        """
        return bytes(self).rsplit(sep, maxsplit)

    def strip(self, chars=None):
        """Strip leading and trailing characters contained in the argument.

        If the argument is omitted or None, strip leading and trailing ASCII whitespace.
        """
        return bytes(self).strip(chars)

    def lstrip(self, chars=None):
        """Strip leading characters contained in the argument.

        If the argument is omitted or None, strip leading ASCII whitespace.
        """
        return bytes(self).lstrip(chars)

    def rstrip(self, chars=None):
        """Strip trailing characters contained in the argument.

        If the argument is omitted or None, strip trailing ASCII whitespace.
        """
        return bytes(self).rstrip(chars)

    def upper(self):
        """Return a copy of data with all ASCII characters converted to uppercase."""
        return bytes(self).upper()

    def lower(self):
        """Return a copy of data with all ASCII characters converted to lowercase."""
        return bytes(self).lower()

    def isupper(self):
        """Return True if all ASCII characters in data are uppercase.

        If there are no cased characters, the method returns False.
        """
        return bytes(self).isupper()

    def islower(self):
        """Return True if all ASCII characters in data are lowercase.

        If there are no cased characters, the method returns False.
        """
        return bytes(self).islower()

    def replace(self, old, new):
        """Return a copy with all occurrences of substring old replaced by new."""
        return bytes(self).replace(old, new)

    def translate(self, table, delete=b""):
        """Return a copy with each character mapped by the given translation table.

          table
            Translation table, which must be a bytes object of length 256.

        All characters occurring in the optional argument delete are removed.
        The remaining characters are mapped through the given translation table.
        """
        return bytes(self).translate(table, delete)

    @property
    def defined(self):
        """Return True if the sequence is defined, False if undefined or partially defined.

        Zero-length sequences are always considered to be defined.
        """
        return True

    @property
    def defined_ranges(self):
        """Return a tuple of the ranges where the sequence contents is defined.

        The return value has the format ((start1, end1), (start2, end2), ...).
        """
        length = len(self)
        if length > 0:
            return ((0, length),)
        else:
            return ()


class _SeqAbstractBaseClass(ABC):
    """Abstract base class for the Seq and MutableSeq classes (PRIVATE).

    Most users will not need to use this class. It is used internally as an
    abstract base class for Seq and MutableSeq, as most of their methods are
    identical.
    """

    __slots__ = ("_data",)
    __array_ufunc__ = None  # turn off numpy Ufuncs

    @abstractmethod
    def __init__(self):
        pass

    def __bytes__(self):
        return bytes(self._data)

    def __repr__(self):
        """Return (truncated) representation of the sequence."""
        data = self._data
        if isinstance(data, _UndefinedSequenceData):
            return f"Seq(None, length={len(self)})"
        if isinstance(data, _PartiallyDefinedSequenceData):
            d = {}
            for position, seq in data._data.items():
                if len(seq) > 60:
                    start = seq[:54].decode("ASCII")
                    end = seq[-3:].decode("ASCII")
                    seq = f"{start}...{end}"
                else:
                    seq = seq.decode("ASCII")
                d[position] = seq
            return "Seq(%r, length=%d)" % (d, len(self))
        if len(data) > 60:
            # Shows the last three letters as it is often useful to see if
            # there is a stop codon at the end of a sequence.
            # Note total length is 54+3+3=60
            start = data[:54].decode("ASCII")
            end = data[-3:].decode("ASCII")
            return f"{self.__class__.__name__}('{start}...{end}')"
        else:
            data = data.decode("ASCII")
            return f"{self.__class__.__name__}('{data}')"

    def __str__(self):
        """Return the full sequence as a python string."""
        return self._data.decode("ASCII")

    def __eq__(self, other):
        """Compare the sequence to another sequence or a string.

        Sequences are equal to each other if their sequence contents is
        identical:

        >>> from Bio.Seq import Seq, MutableSeq
        >>> seq1 = Seq("ACGT")
        >>> seq2 = Seq("ACGT")
        >>> mutable_seq = MutableSeq("ACGT")
        >>> seq1 == seq2
        True
        >>> seq1 == mutable_seq
        True
        >>> seq1 == "ACGT"
        True

        Note that the sequence objects themselves are not identical to each
        other:

        >>> id(seq1) == id(seq2)
        False
        >>> seq1 is seq2
        False

        Sequences can also be compared to strings, ``bytes``, and ``bytearray``
        objects:

        >>> seq1 == "ACGT"
        True
        >>> seq1 == b"ACGT"
        True
        >>> seq1 == bytearray(b"ACGT")
        True
        """
        if isinstance(other, _SeqAbstractBaseClass):
            return self._data == other._data
        elif isinstance(other, str):
            return self._data == other.encode("ASCII")
        else:
            return self._data == other

    def __lt__(self, other):
        """Implement the less-than operand."""
        if isinstance(other, _SeqAbstractBaseClass):
            return self._data < other._data
        elif isinstance(other, str):
            return self._data < other.encode("ASCII")
        else:
            return self._data < other

    def __le__(self, other):
        """Implement the less-than or equal operand."""
        if isinstance(other, _SeqAbstractBaseClass):
            return self._data <= other._data
        elif isinstance(other, str):
            return self._data <= other.encode("ASCII")
        else:
            return self._data <= other

    def __gt__(self, other):
        """Implement the greater-than operand."""
        if isinstance(other, _SeqAbstractBaseClass):
            return self._data > other._data
        elif isinstance(other, str):
            return self._data > other.encode("ASCII")
        else:
            return self._data > other

    def __ge__(self, other):
        """Implement the greater-than or equal operand."""
        if isinstance(other, _SeqAbstractBaseClass):
            return self._data >= other._data
        elif isinstance(other, str):
            return self._data >= other.encode("ASCII")
        else:
            return self._data >= other

    def __len__(self):
        """Return the length of the sequence."""
        return len(self._data)

    def __iter__(self):
        """Return an iterable of the sequence."""
        return self._data.decode("ASCII").__iter__()

    def __getitem__(self, index):
        """Return a subsequence as a single letter or as a sequence object.

        If the index is an integer, a single letter is returned as a Python
        string:

        >>> seq = Seq('ACTCGACGTCG')
        >>> seq[5]
        'A'

        Otherwise, a new sequence object of the same class is returned:

        >>> seq[5:8]
        Seq('ACG')
        >>> mutable_seq = MutableSeq('ACTCGACGTCG')
        >>> mutable_seq[5:8]
        MutableSeq('ACG')
        """
        if isinstance(index, numbers.Integral):
            # Return a single letter as a string
            return chr(self._data[index])
        else:
            # Return the (sub)sequence as another Seq/MutableSeq object
            return self.__class__(self._data[index])

    def __add__(self, other):
        """Add a sequence or string to this sequence.

        >>> from Bio.Seq import Seq, MutableSeq
        >>> Seq("MELKI") + "LV"
        Seq('MELKILV')
        >>> MutableSeq("MELKI") + "LV"
        MutableSeq('MELKILV')
        """
        if isinstance(other, _SeqAbstractBaseClass):
            return self.__class__(self._data + other._data)
        elif isinstance(other, str):
            return self.__class__(self._data + other.encode("ASCII"))
        else:
            # If other is a SeqRecord, then SeqRecord's __radd__ will handle
            # this. If not, returning NotImplemented will trigger a TypeError.
            return NotImplemented

    def __radd__(self, other):
        """Add a sequence string on the left.

        >>> from Bio.Seq import Seq, MutableSeq
        >>> "LV" + Seq("MELKI")
        Seq('LVMELKI')
        >>> "LV" + MutableSeq("MELKI")
        MutableSeq('LVMELKI')

        Adding two sequence objects is handled via the __add__ method.
        """
        if isinstance(other, str):
            return self.__class__(other.encode("ASCII") + self._data)
        else:
            return NotImplemented

    def __mul__(self, other):
        """Multiply sequence by integer.

        >>> from Bio.Seq import Seq, MutableSeq
        >>> Seq('ATG') * 2
        Seq('ATGATG')
        >>> MutableSeq('ATG') * 2
        MutableSeq('ATGATG')
        """
        if not isinstance(other, numbers.Integral):
            raise TypeError(f"can't multiply {self.__class__.__name__} by non-int type")
        # we would like to simply write
        # data = self._data * other
        # here, but currently that causes a bug on PyPy if self._data is a
        # bytearray and other is a numpy integer. Using this workaround:
        data = self._data.__mul__(other)
        return self.__class__(data)

    def __rmul__(self, other):
        """Multiply integer by sequence.

        >>> from Bio.Seq import Seq
        >>> 2 * Seq('ATG')
        Seq('ATGATG')
        """
        if not isinstance(other, numbers.Integral):
            raise TypeError(f"can't multiply {self.__class__.__name__} by non-int type")
        # we would like to simply write
        # data = self._data * other
        # here, but currently that causes a bug on PyPy if self._data is a
        # bytearray and other is a numpy integer. Using this workaround:
        data = self._data.__mul__(other)
        return self.__class__(data)

    def __imul__(self, other):
        """Multiply the sequence object by other and assign.

        >>> from Bio.Seq import Seq
        >>> seq = Seq('ATG')
        >>> seq *= 2
        >>> seq
        Seq('ATGATG')

        Note that this is different from in-place multiplication. The ``seq``
        variable is reassigned to the multiplication result, but any variable
        pointing to ``seq`` will remain unchanged:

        >>> seq = Seq('ATG')
        >>> seq2 = seq
        >>> id(seq) == id(seq2)
        True
        >>> seq *= 2
        >>> seq
        Seq('ATGATG')
        >>> seq2
        Seq('ATG')
        >>> id(seq) == id(seq2)
        False
        """
        if not isinstance(other, numbers.Integral):
            raise TypeError(f"can't multiply {self.__class__.__name__} by non-int type")
        # we would like to simply write
        # data = self._data * other
        # here, but currently that causes a bug on PyPy if self._data is a
        # bytearray and other is a numpy integer. Using this workaround:
        data = self._data.__mul__(other)
        return self.__class__(data)

    def count(self, sub, start=None, end=None):
        """Return a non-overlapping count, like that of a python string.

        The number of occurrences of substring argument sub in the
        (sub)sequence given by [start:end] is returned as an integer.
        Optional arguments start and end are interpreted as in slice
        notation.

        Arguments:
         - sub - a string or another Seq object to look for
         - start - optional integer, slice start
         - end - optional integer, slice end

        e.g.

        >>> from Bio.Seq import Seq
        >>> my_seq = Seq("AAAATGA")
        >>> print(my_seq.count("A"))
        5
        >>> print(my_seq.count("ATG"))
        1
        >>> print(my_seq.count(Seq("AT")))
        1
        >>> print(my_seq.count("AT", 2, -1))
        1

        HOWEVER, please note because the ``count`` method of Seq and MutableSeq
        objects, like that of Python strings, do a non-overlapping search, this
        may not give the answer you expect:

        >>> "AAAA".count("AA")
        2
        >>> print(Seq("AAAA").count("AA"))
        2

        For an overlapping search, use the ``count_overlap`` method:

        >>> print(Seq("AAAA").count_overlap("AA"))
        3
        """
        if isinstance(sub, MutableSeq):
            sub = sub._data
        elif isinstance(sub, Seq):
            sub = bytes(sub)
        elif isinstance(sub, str):
            sub = sub.encode("ASCII")
        elif not isinstance(sub, (bytes, bytearray)):
            raise TypeError(
                "a Seq, MutableSeq, str, bytes, or bytearray object is required, not '%s'"
                % type(sub)
            )
        return self._data.count(sub, start, end)

    def count_overlap(self, sub, start=None, end=None):
        """Return an overlapping count.

        Returns an integer, the number of occurrences of substring
        argument sub in the (sub)sequence given by [start:end].
        Optional arguments start and end are interpreted as in slice
        notation.

        Arguments:
         - sub - a string or another Seq object to look for
         - start - optional integer, slice start
         - end - optional integer, slice end

        e.g.

        >>> from Bio.Seq import Seq
        >>> print(Seq("AAAA").count_overlap("AA"))
        3
        >>> print(Seq("ATATATATA").count_overlap("ATA"))
        4
        >>> print(Seq("ATATATATA").count_overlap("ATA", 3, -1))
        1

        For a non-overlapping search, use the ``count`` method:

        >>> print(Seq("AAAA").count("AA"))
        2

        Where substrings do not overlap, ``count_overlap`` behaves the same as
        the ``count`` method:

        >>> from Bio.Seq import Seq
        >>> my_seq = Seq("AAAATGA")
        >>> print(my_seq.count_overlap("A"))
        5
        >>> my_seq.count_overlap("A") == my_seq.count("A")
        True
        >>> print(my_seq.count_overlap("ATG"))
        1
        >>> my_seq.count_overlap("ATG") == my_seq.count("ATG")
        True
        >>> print(my_seq.count_overlap(Seq("AT")))
        1
        >>> my_seq.count_overlap(Seq("AT")) == my_seq.count(Seq("AT"))
        True
        >>> print(my_seq.count_overlap("AT", 2, -1))
        1
        >>> my_seq.count_overlap("AT", 2, -1) == my_seq.count("AT", 2, -1)
        True

        HOWEVER, do not use this method for such cases because the
        count() method is much for efficient.
        """
        if isinstance(sub, MutableSeq):
            sub = sub._data
        elif isinstance(sub, Seq):
            sub = bytes(sub)
        elif isinstance(sub, str):
            sub = sub.encode("ASCII")
        elif not isinstance(sub, (bytes, bytearray)):
            raise TypeError(
                "a Seq, MutableSeq, str, bytes, or bytearray object is required, not '%s'"
                % type(sub)
            )
        data = self._data
        overlap_count = 0
        while True:
            start = data.find(sub, start, end) + 1
            if start != 0:
                overlap_count += 1
            else:
                return overlap_count

    def __contains__(self, item):
        """Return True if item is a subsequence of the sequence, and False otherwise.

        e.g.

        >>> from Bio.Seq import Seq, MutableSeq
        >>> my_dna = Seq("ATATGAAATTTGAAAA")
        >>> "AAA" in my_dna
        True
        >>> Seq("AAA") in my_dna
        True
        >>> MutableSeq("AAA") in my_dna
        True
        """
        if isinstance(item, _SeqAbstractBaseClass):
            item = bytes(item)
        elif isinstance(item, str):
            item = item.encode("ASCII")
        return item in self._data

    def find(self, sub, start=None, end=None):
        """Return the lowest index in the sequence where subsequence sub is found.

        With optional arguments start and end, return the lowest index in the
        sequence such that the subsequence sub is contained within the sequence
        region [start:end].

        Arguments:
         - sub - a string or another Seq or MutableSeq object to search for
         - start - optional integer, slice start
         - end - optional integer, slice end

        Returns -1 if the subsequence is NOT found.

        e.g. Locating the first typical start codon, AUG, in an RNA sequence:

        >>> from Bio.Seq import Seq
        >>> my_rna = Seq("GUCAUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAGUUG")
        >>> my_rna.find("AUG")
        3

        The next typical start codon can then be found by starting the search
        at position 4:

        >>> my_rna.find("AUG", 4)
        15
        """
        if isinstance(sub, _SeqAbstractBaseClass):
            sub = bytes(sub)
        elif isinstance(sub, str):
            sub = sub.encode("ASCII")
        elif not isinstance(sub, (bytes, bytearray)):
            raise TypeError(
                "a Seq, MutableSeq, str, bytes, or bytearray object is required, not '%s'"
                % type(sub)
            )
        return self._data.find(sub, start, end)

    def rfind(self, sub, start=None, end=None):
        """Return the highest index in the sequence where subsequence sub is found.

        With optional arguments start and end, return the highest index in the
        sequence such that the subsequence sub is contained within the sequence
        region [start:end].

        Arguments:
         - sub - a string or another Seq or MutableSeq object to search for
         - start - optional integer, slice start
         - end - optional integer, slice end

        Returns -1 if the subsequence is NOT found.

        e.g. Locating the last typical start codon, AUG, in an RNA sequence:

        >>> from Bio.Seq import Seq
        >>> my_rna = Seq("GUCAUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAGUUG")
        >>> my_rna.rfind("AUG")
        15

        The location of the typical start codon before that can be found by
        ending the search at position 15:

        >>> my_rna.rfind("AUG", end=15)
        3
        """
        if isinstance(sub, _SeqAbstractBaseClass):
            sub = bytes(sub)
        elif isinstance(sub, str):
            sub = sub.encode("ASCII")
        elif not isinstance(sub, (bytes, bytearray)):
            raise TypeError(
                "a Seq, MutableSeq, str, bytes, or bytearray object is required, not '%s'"
                % type(sub)
            )
        return self._data.rfind(sub, start, end)

    def index(self, sub, start=None, end=None):
        """Return the lowest index in the sequence where subsequence sub is found.

        With optional arguments start and end, return the lowest index in the
        sequence such that the subsequence sub is contained within the sequence
        region [start:end].

        Arguments:
         - sub - a string or another Seq or MutableSeq object to search for
         - start - optional integer, slice start
         - end - optional integer, slice end

        Raises a ValueError if the subsequence is NOT found.

        e.g. Locating the first typical start codon, AUG, in an RNA sequence:

        >>> from Bio.Seq import Seq
        >>> my_rna = Seq("GUCAUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAGUUG")
        >>> my_rna.index("AUG")
        3

        The next typical start codon can then be found by starting the search
        at position 4:

        >>> my_rna.index("AUG", 4)
        15

        This method performs the same search as the ``find`` method.  However,
        if the subsequence is not found, ``find`` returns -1 which ``index``
        raises a ValueError:

        >>> my_rna.index("T")
        Traceback (most recent call last):
                   ...
        ValueError: ...
        >>> my_rna.find("T")
        -1
        """
        if isinstance(sub, MutableSeq):
            sub = sub._data
        elif isinstance(sub, Seq):
            sub = bytes(sub)
        elif isinstance(sub, str):
            sub = sub.encode("ASCII")
        elif not isinstance(sub, (bytes, bytearray)):
            raise TypeError(
                "a Seq, MutableSeq, str, bytes, or bytearray object is required, not '%s'"
                % type(sub)
            )
        return self._data.index(sub, start, end)

    def rindex(self, sub, start=None, end=None):
        """Return the highest index in the sequence where subsequence sub is found.

        With optional arguments start and end, return the highest index in the
        sequence such that the subsequence sub is contained within the sequence
        region [start:end].

        Arguments:
         - sub - a string or another Seq or MutableSeq object to search for
         - start - optional integer, slice start
         - end - optional integer, slice end

        Returns -1 if the subsequence is NOT found.

        e.g. Locating the last typical start codon, AUG, in an RNA sequence:

        >>> from Bio.Seq import Seq
        >>> my_rna = Seq("GUCAUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAGUUG")
        >>> my_rna.rindex("AUG")
        15

        The location of the typical start codon before that can be found by
        ending the search at position 15:

        >>> my_rna.rindex("AUG", end=15)
        3

        This method performs the same search as the ``rfind`` method.  However,
        if the subsequence is not found, ``rfind`` returns -1 which ``rindex``
        raises a ValueError:

        >>> my_rna.rindex("T")
        Traceback (most recent call last):
                   ...
        ValueError: ...
        >>> my_rna.rfind("T")
        -1
        """
        if isinstance(sub, MutableSeq):
            sub = sub._data
        elif isinstance(sub, Seq):
            sub = bytes(sub)
        elif isinstance(sub, str):
            sub = sub.encode("ASCII")
        elif not isinstance(sub, (bytes, bytearray)):
            raise TypeError(
                "a Seq, MutableSeq, str, bytes, or bytearray object is required, not '%s'"
                % type(sub)
            )
        return self._data.rindex(sub, start, end)

    def startswith(self, prefix, start=None, end=None):
        """Return True if the sequence starts with the given prefix, False otherwise.

        Return True if the sequence starts with the specified prefix
        (a string or another Seq object), False otherwise.
        With optional start, test sequence beginning at that position.
        With optional end, stop comparing sequence at that position.
        prefix can also be a tuple of strings to try.  e.g.

        >>> from Bio.Seq import Seq
        >>> my_rna = Seq("GUCAUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAGUUG")
        >>> my_rna.startswith("GUC")
        True
        >>> my_rna.startswith("AUG")
        False
        >>> my_rna.startswith("AUG", 3)
        True
        >>> my_rna.startswith(("UCC", "UCA", "UCG"), 1)
        True
        """
        if isinstance(prefix, tuple):
            prefix = tuple(
                bytes(p) if isinstance(p, _SeqAbstractBaseClass) else p.encode("ASCII")
                for p in prefix
            )
        elif isinstance(prefix, _SeqAbstractBaseClass):
            prefix = bytes(prefix)
        elif isinstance(prefix, str):
            prefix = prefix.encode("ASCII")
        return self._data.startswith(prefix, start, end)

    def endswith(self, suffix, start=None, end=None):
        """Return True if the sequence ends with the given suffix, False otherwise.

        Return True if the sequence ends with the specified suffix
        (a string or another Seq object), False otherwise.
        With optional start, test sequence beginning at that position.
        With optional end, stop comparing sequence at that position.
        suffix can also be a tuple of strings to try.  e.g.

        >>> from Bio.Seq import Seq
        >>> my_rna = Seq("GUCAUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAGUUG")
        >>> my_rna.endswith("UUG")
        True
        >>> my_rna.endswith("AUG")
        False
        >>> my_rna.endswith("AUG", 0, 18)
        True
        >>> my_rna.endswith(("UCC", "UCA", "UUG"))
        True
        """
        if isinstance(suffix, tuple):
            suffix = tuple(
                bytes(p) if isinstance(p, _SeqAbstractBaseClass) else p.encode("ASCII")
                for p in suffix
            )
        elif isinstance(suffix, _SeqAbstractBaseClass):
            suffix = bytes(suffix)
        elif isinstance(suffix, str):
            suffix = suffix.encode("ASCII")
        return self._data.endswith(suffix, start, end)

    def split(self, sep=None, maxsplit=-1):
        """Return a list of subsequences when splitting the sequence by separator sep.

        Return a list of the subsequences in the sequence (as Seq objects),
        using sep as the delimiter string.  If maxsplit is given, at
        most maxsplit splits are done.  If maxsplit is omitted, all
        splits are made.

        For consistency with the ``split`` method of Python strings, any
        whitespace (tabs, spaces, newlines) is a separator if sep is None, the
        default value

        e.g.

        >>> from Bio.Seq import Seq
        >>> my_rna = Seq("GUCAUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAGUUG")
        >>> my_aa = my_rna.translate()
        >>> my_aa
        Seq('VMAIVMGR*KGAR*L')
        >>> for pep in my_aa.split("*"):
        ...     pep
        Seq('VMAIVMGR')
        Seq('KGAR')
        Seq('L')
        >>> for pep in my_aa.split("*", 1):
        ...     pep
        Seq('VMAIVMGR')
        Seq('KGAR*L')

        See also the rsplit method, which splits the sequence starting from the
        end:

        >>> for pep in my_aa.rsplit("*", 1):
        ...     pep
        Seq('VMAIVMGR*KGAR')
        Seq('L')
        """
        if isinstance(sep, _SeqAbstractBaseClass):
            sep = bytes(sep)
        elif isinstance(sep, str):
            sep = sep.encode("ASCII")
        return [Seq(part) for part in self._data.split(sep, maxsplit)]

    def rsplit(self, sep=None, maxsplit=-1):
        """Return a list of subsequences by splitting the sequence from the right.

        Return a list of the subsequences in the sequence (as Seq objects),
        using sep as the delimiter string.  If maxsplit is given, at
        most maxsplit splits are done.  If maxsplit is omitted, all
        splits are made.

        For consistency with the ``rsplit`` method of Python strings, any
        whitespace (tabs, spaces, newlines) is a separator if sep is None, the
        default value

        e.g.

        >>> from Bio.Seq import Seq
        >>> my_rna = Seq("GUCAUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAGUUG")
        >>> my_aa = my_rna.translate()
        >>> my_aa
        Seq('VMAIVMGR*KGAR*L')
        >>> for pep in my_aa.rsplit("*"):
        ...     pep
        Seq('VMAIVMGR')
        Seq('KGAR')
        Seq('L')
        >>> for pep in my_aa.rsplit("*", 1):
        ...     pep
        Seq('VMAIVMGR*KGAR')
        Seq('L')

        See also the split method, which splits the sequence starting from the
        beginning:

        >>> for pep in my_aa.split("*", 1):
        ...     pep
        Seq('VMAIVMGR')
        Seq('KGAR*L')
        """
        if isinstance(sep, _SeqAbstractBaseClass):
            sep = bytes(sep)
        elif isinstance(sep, str):
            sep = sep.encode("ASCII")
        return [Seq(part) for part in self._data.rsplit(sep, maxsplit)]

    def strip(self, chars=None, inplace=False):
        """Return a sequence object with leading and trailing ends stripped.

        With default arguments, leading and trailing whitespace is removed:

        >>> seq = Seq(" ACGT ")
        >>> seq.strip()
        Seq('ACGT')
        >>> seq
        Seq(' ACGT ')

        If ``chars`` is given and not ``None``, remove characters in ``chars``
        instead.  The order of the characters to be removed is not important:

        >>> Seq("ACGTACGT").strip("TGCA")
        Seq('')

        A copy of the sequence is returned if ``inplace`` is ``False`` (the
        default value).  If ``inplace`` is ``True``, the sequence is stripped
        in-place and returned.

        >>> seq = MutableSeq(" ACGT ")
        >>> seq.strip(inplace=False)
        MutableSeq('ACGT')
        >>> seq
        MutableSeq(' ACGT ')
        >>> seq.strip(inplace=True)
        MutableSeq('ACGT')
        >>> seq
        MutableSeq('ACGT')

        As ``Seq`` objects are immutable, a ``TypeError`` is raised if ``strip``
        is called on a ``Seq`` object with ``inplace=True``.

        See also the lstrip and rstrip methods.
        """
        if isinstance(chars, _SeqAbstractBaseClass):
            chars = bytes(chars)
        elif isinstance(chars, str):
            chars = chars.encode("ASCII")
        try:
            data = self._data.strip(chars)
        except TypeError:
            raise TypeError(
                "argument must be None or a string, Seq, MutableSeq, or bytes-like object"
            ) from None
        if inplace:
            if not isinstance(self._data, bytearray):
                raise TypeError("Sequence is immutable")
            self._data[:] = data
            return self
        else:
            return self.__class__(data)

    def lstrip(self, chars=None, inplace=False):
        """Return a sequence object with leading and trailing ends stripped.

        With default arguments, leading whitespace is removed:

        >>> seq = Seq(" ACGT ")
        >>> seq.lstrip()
        Seq('ACGT ')
        >>> seq
        Seq(' ACGT ')

        If ``chars`` is given and not ``None``, remove characters in ``chars``
        from the leading end instead.  The order of the characters to be removed
        is not important:

        >>> Seq("ACGACGTTACG").lstrip("GCA")
        Seq('TTACG')

        A copy of the sequence is returned if ``inplace`` is ``False`` (the
        default value).  If ``inplace`` is ``True``, the sequence is stripped
        in-place and returned.

        >>> seq = MutableSeq(" ACGT ")
        >>> seq.lstrip(inplace=False)
        MutableSeq('ACGT ')
        >>> seq
        MutableSeq(' ACGT ')
        >>> seq.lstrip(inplace=True)
        MutableSeq('ACGT ')
        >>> seq
        MutableSeq('ACGT ')

        As ``Seq`` objects are immutable, a ``TypeError`` is raised if
        ``lstrip`` is called on a ``Seq`` object with ``inplace=True``.

        See also the strip and rstrip methods.
        """
        if isinstance(chars, _SeqAbstractBaseClass):
            chars = bytes(chars)
        elif isinstance(chars, str):
            chars = chars.encode("ASCII")
        try:
            data = self._data.lstrip(chars)
        except TypeError:
            raise TypeError(
                "argument must be None or a string, Seq, MutableSeq, or bytes-like object"
            ) from None
        if inplace:
            if not isinstance(self._data, bytearray):
                raise TypeError("Sequence is immutable")
            self._data[:] = data
            return self
        else:
            return self.__class__(data)

    def rstrip(self, chars=None, inplace=False):
        """Return a sequence object with trailing ends stripped.

        With default arguments, trailing whitespace is removed:

        >>> seq = Seq(" ACGT ")
        >>> seq.rstrip()
        Seq(' ACGT')
        >>> seq
        Seq(' ACGT ')

        If ``chars`` is given and not ``None``, remove characters in ``chars``
        from the trailing end instead.  The order of the characters to be
        removed is not important:

        >>> Seq("ACGACGTTACG").rstrip("GCA")
        Seq('ACGACGTT')

        A copy of the sequence is returned if ``inplace`` is ``False`` (the
        default value).  If ``inplace`` is ``True``, the sequence is stripped
        in-place and returned.

        >>> seq = MutableSeq(" ACGT ")
        >>> seq.rstrip(inplace=False)
        MutableSeq(' ACGT')
        >>> seq
        MutableSeq(' ACGT ')
        >>> seq.rstrip(inplace=True)
        MutableSeq(' ACGT')
        >>> seq
        MutableSeq(' ACGT')

        As ``Seq`` objects are immutable, a ``TypeError`` is raised if
        ``rstrip`` is called on a ``Seq`` object with ``inplace=True``.

        See also the strip and lstrip methods.
        """
        if isinstance(chars, _SeqAbstractBaseClass):
            chars = bytes(chars)
        elif isinstance(chars, str):
            chars = chars.encode("ASCII")
        try:
            data = self._data.rstrip(chars)
        except TypeError:
            raise TypeError(
                "argument must be None or a string, Seq, MutableSeq, or bytes-like object"
            ) from None
        if inplace:
            if not isinstance(self._data, bytearray):
                raise TypeError("Sequence is immutable")
            self._data[:] = data
            return self
        else:
            return self.__class__(data)

    def upper(self, inplace=False):
        """Return the sequence in upper case.

        An upper-case copy of the sequence is returned if inplace is False,
        the default value:

        >>> from Bio.Seq import Seq, MutableSeq
        >>> my_seq = Seq("VHLTPeeK*")
        >>> my_seq
        Seq('VHLTPeeK*')
        >>> my_seq.lower()
        Seq('vhltpeek*')
        >>> my_seq.upper()
        Seq('VHLTPEEK*')
        >>> my_seq
        Seq('VHLTPeeK*')

        The sequence is modified in-place and returned if inplace is True:

        >>> my_seq = MutableSeq("VHLTPeeK*")
        >>> my_seq
        MutableSeq('VHLTPeeK*')
        >>> my_seq.lower()
        MutableSeq('vhltpeek*')
        >>> my_seq.upper()
        MutableSeq('VHLTPEEK*')
        >>> my_seq
        MutableSeq('VHLTPeeK*')

        >>> my_seq.lower(inplace=True)
        MutableSeq('vhltpeek*')
        >>> my_seq
        MutableSeq('vhltpeek*')
        >>> my_seq.upper(inplace=True)
        MutableSeq('VHLTPEEK*')
        >>> my_seq
        MutableSeq('VHLTPEEK*')

        As ``Seq`` objects are immutable, a ``TypeError`` is raised if
        ``upper`` is called on a ``Seq`` object with ``inplace=True``.

        See also the ``lower`` method.
        """
        data = self._data.upper()
        if inplace:
            if not isinstance(self._data, bytearray):
                raise TypeError("Sequence is immutable")
            self._data[:] = data
            return self
        else:
            return self.__class__(data)

    def lower(self, inplace=False):
        """Return the sequence in lower case.

        An lower-case copy of the sequence is returned if inplace is False,
        the default value:

        >>> from Bio.Seq import Seq, MutableSeq
        >>> my_seq = Seq("VHLTPeeK*")
        >>> my_seq
        Seq('VHLTPeeK*')
        >>> my_seq.lower()
        Seq('vhltpeek*')
        >>> my_seq.upper()
        Seq('VHLTPEEK*')
        >>> my_seq
        Seq('VHLTPeeK*')

        The sequence is modified in-place and returned if inplace is True:

        >>> my_seq = MutableSeq("VHLTPeeK*")
        >>> my_seq
        MutableSeq('VHLTPeeK*')
        >>> my_seq.lower()
        MutableSeq('vhltpeek*')
        >>> my_seq.upper()
        MutableSeq('VHLTPEEK*')
        >>> my_seq
        MutableSeq('VHLTPeeK*')

        >>> my_seq.lower(inplace=True)
        MutableSeq('vhltpeek*')
        >>> my_seq
        MutableSeq('vhltpeek*')
        >>> my_seq.upper(inplace=True)
        MutableSeq('VHLTPEEK*')
        >>> my_seq
        MutableSeq('VHLTPEEK*')

        As ``Seq`` objects are immutable, a ``TypeError`` is raised if
        ``lower`` is called on a ``Seq`` object with ``inplace=True``.

        See also the ``upper`` method.
        """
        data = self._data.lower()
        if inplace:
            if not isinstance(self._data, bytearray):
                raise TypeError("Sequence is immutable")
            self._data[:] = data
            return self
        else:
            return self.__class__(data)

    def isupper(self):
        """Return True if all ASCII characters in data are uppercase.

        If there are no cased characters, the method returns False.
        """
        return self._data.isupper()

    def islower(self):
        """Return True if all ASCII characters in data are lowercase.

        If there are no cased characters, the method returns False.
        """
        return self._data.islower()

    def translate(
        self, table="Standard", stop_symbol="*", to_stop=False, cds=False, gap="-"
    ):
        """Turn a nucleotide sequence into a protein sequence by creating a new sequence object.

        This method will translate DNA or RNA sequences. It should not
        be used on protein sequences as any result will be biologically
        meaningless.

        Arguments:
         - table - Which codon table to use?  This can be either a name
           (string), an NCBI identifier (integer), or a CodonTable
           object (useful for non-standard genetic codes).  This
           defaults to the "Standard" table.
         - stop_symbol - Single character string, what to use for
           terminators.  This defaults to the asterisk, "*".
         - to_stop - Boolean, defaults to False meaning do a full
           translation continuing on past any stop codons (translated as the
           specified stop_symbol).  If True, translation is terminated at
           the first in frame stop codon (and the stop_symbol is not
           appended to the returned protein sequence).
         - cds - Boolean, indicates this is a complete CDS.  If True,
           this checks the sequence starts with a valid alternative start
           codon (which will be translated as methionine, M), that the
           sequence length is a multiple of three, and that there is a
           single in frame stop codon at the end (this will be excluded
           from the protein sequence, regardless of the to_stop option).
           If these tests fail, an exception is raised.
         - gap - Single character string to denote symbol used for gaps.
           Defaults to the minus sign.

        A ``Seq`` object is returned if ``translate`` is called on a ``Seq``
        object; a ``MutableSeq`` object is returned if ``translate`` is called
        pn a ``MutableSeq`` object.

        e.g. Using the standard table:

        >>> coding_dna = Seq("GTGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
        >>> coding_dna.translate()
        Seq('VAIVMGR*KGAR*')
        >>> coding_dna.translate(stop_symbol="@")
        Seq('VAIVMGR@KGAR@')
        >>> coding_dna.translate(to_stop=True)
        Seq('VAIVMGR')

        Now using NCBI table 2, where TGA is not a stop codon:

        >>> coding_dna.translate(table=2)
        Seq('VAIVMGRWKGAR*')
        >>> coding_dna.translate(table=2, to_stop=True)
        Seq('VAIVMGRWKGAR')

        In fact, GTG is an alternative start codon under NCBI table 2, meaning
        this sequence could be a complete CDS:

        >>> coding_dna.translate(table=2, cds=True)
        Seq('MAIVMGRWKGAR')

        It isn't a valid CDS under NCBI table 1, due to both the start codon
        and also the in frame stop codons:

        >>> coding_dna.translate(table=1, cds=True)
        Traceback (most recent call last):
            ...
        Bio.Data.CodonTable.TranslationError: First codon 'GTG' is not a start codon

        If the sequence has no in-frame stop codon, then the to_stop argument
        has no effect:

        >>> coding_dna2 = Seq("TTGGCCATTGTAATGGGCCGC")
        >>> coding_dna2.translate()
        Seq('LAIVMGR')
        >>> coding_dna2.translate(to_stop=True)
        Seq('LAIVMGR')

        NOTE - Ambiguous codons like "TAN" or "NNN" could be an amino acid
        or a stop codon.  These are translated as "X".  Any invalid codon
        (e.g. "TA?" or "T-A") will throw a TranslationError.

        NOTE - This does NOT behave like the python string's translate
        method.  For that use str(my_seq).translate(...) instead
        """
        try:
            data = str(self)
        except UndefinedSequenceError:
            # translating an undefined sequence yields an undefined
            # sequence with the length divided by 3
            n = len(self)
            if n % 3 != 0:
                warnings.warn(
                    "Partial codon, len(sequence) not a multiple of three. "
                    "This may become an error in future.",
                    BiopythonWarning,
                )
            return Seq(None, n // 3)

        return self.__class__(
            _translate_str(str(self), table, stop_symbol, to_stop, cds, gap=gap)
        )

    def complement(self, inplace=None):
        """Return the complement as a DNA sequence.

        >>> Seq("CGA").complement()
        Seq('GCT')

        Any U in the sequence is treated as a T:

        >>> Seq("CGAUT").complement(inplace=False)
        Seq('GCTAA')

        In contrast, ``complement_rna`` returns an RNA sequence:

        >>> Seq("CGAUT").complement_rna()
        Seq('GCUAA')

        The sequence is modified in-place and returned if inplace is True:

        >>> my_seq = MutableSeq("CGA")
        >>> my_seq
        MutableSeq('CGA')
        >>> my_seq.complement(inplace=False)
        MutableSeq('GCT')
        >>> my_seq
        MutableSeq('CGA')

        >>> my_seq.complement(inplace=True)
        MutableSeq('GCT')
        >>> my_seq
        MutableSeq('GCT')

        As ``Seq`` objects are immutable, a ``TypeError`` is raised if
        ``complement_rna`` is called on a ``Seq`` object with ``inplace=True``.
        """
        ttable = _dna_complement_table
        try:
            if inplace is None:
                # deprecated
                if isinstance(self._data, bytearray):  # MutableSeq
                    warnings.warn(
                        "mutable_seq.complement() will change in the near "
                        "future and will no longer change the sequence in-"
                        "place by default. Please use\n"
                        "\n"
                        "mutable_seq.complement(inplace=True)\n"
                        "\n"
                        "if you want to continue to use this method to change "
                        "a mutable sequence in-place.",
                        BiopythonDeprecationWarning,
                    )
                    inplace = True
                if isinstance(self._data, _PartiallyDefinedSequenceData):
                    for seq in self._data._data.values():
                        if b"U" in seq or b"u" in seq:
                            warnings.warn(
                                "seq.complement() will change in the near "
                                "future to always return DNA nucleotides only. "
                                "Please use\n"
                                "\n"
                                "seq.complement_rna()\n"
                                "\n"
                                "if you want to receive an RNA sequence instead.",
                                BiopythonDeprecationWarning,
                            )
                            for seq in self._data._data.values():
                                if b"t" in seq or b"T" in seq:
                                    raise ValueError("Mixed RNA/DNA found")
                            ttable = _rna_complement_table
                            break

                elif b"U" in self._data or b"u" in self._data:
                    warnings.warn(
                        "seq.complement() will change in the near future to "
                        "always return DNA nucleotides only. "
                        "Please use\n"
                        "\n"
                        "seq.complement_rna()\n"
                        "\n"
                        "if you want to receive an RNA sequence instead.",
                        BiopythonDeprecationWarning,
                    )
                    if b"t" in self._data or b"T" in self._data:
                        raise ValueError("Mixed RNA/DNA found")
                    ttable = _rna_complement_table
            data = self._data.translate(ttable)
        except UndefinedSequenceError:
            # complement of an undefined sequence is an undefined sequence
            # of the same length
            return self
        if inplace:
            if not isinstance(self._data, bytearray):
                raise TypeError("Sequence is immutable")
            self._data[:] = data
            return self
        return self.__class__(data)

    def complement_rna(self, inplace=False):
        """Return the complement as an RNA sequence.

        >>> Seq("CGA").complement_rna()
        Seq('GCU')

        Any T in the sequence is treated as a U:

        >>> Seq("CGAUT").complement_rna()
        Seq('GCUAA')

        In contrast, ``complement`` returns a DNA sequence by default:

        >>> Seq("CGA").complement()
        Seq('GCT')

        The sequence is modified in-place and returned if inplace is True:

        >>> my_seq = MutableSeq("CGA")
        >>> my_seq
        MutableSeq('CGA')
        >>> my_seq.complement_rna()
        MutableSeq('GCU')
        >>> my_seq
        MutableSeq('CGA')

        >>> my_seq.complement_rna(inplace=True)
        MutableSeq('GCU')
        >>> my_seq
        MutableSeq('GCU')

        As ``Seq`` objects are immutable, a ``TypeError`` is raised if
        ``complement_rna`` is called on a ``Seq`` object with ``inplace=True``.
        """
        try:
            data = self._data.translate(_rna_complement_table)
        except UndefinedSequenceError:
            # complement of an undefined sequence is an undefined sequence
            # of the same length
            return self
        if inplace:
            if not isinstance(self._data, bytearray):
                raise TypeError("Sequence is immutable")
            self._data[:] = data
            return self
        return self.__class__(data)

    def reverse_complement(self, inplace=None):
        """Return the reverse complement as a DNA sequence.

        >>> Seq("CGA").reverse_complement(inplace=False)
        Seq('TCG')

        Any U in the sequence is treated as a T:

        >>> Seq("CGAUT").reverse_complement(inplace=False)
        Seq('AATCG')

        In contrast, ``reverse_complement_rna`` returns an RNA sequence:

        >>> Seq("CGA").reverse_complement_rna()
        Seq('UCG')

        The sequence is modified in-place and returned if inplace is True:

        >>> my_seq = MutableSeq("CGA")
        >>> my_seq
        MutableSeq('CGA')
        >>> my_seq.reverse_complement(inplace=False)
        MutableSeq('TCG')
        >>> my_seq
        MutableSeq('CGA')

        >>> my_seq.reverse_complement(inplace=True)
        MutableSeq('TCG')
        >>> my_seq
        MutableSeq('TCG')

        As ``Seq`` objects are immutable, a ``TypeError`` is raised if
        ``reverse_complement`` is called on a ``Seq`` object with
        ``inplace=True``.
        """
        try:
            if inplace is None:
                # deprecated
                if isinstance(self._data, bytearray):  # MutableSeq
                    warnings.warn(
                        "mutable_seq.reverse_complement() will change in the "
                        "near future and will no longer change the sequence in-"
                        "place by default. Please use\n"
                        "\n"
                        "mutable_seq.reverse_complement(inplace=True)\n"
                        "\n"
                        "if you want to continue to use this method to change "
                        "a mutable sequence in-place.",
                        BiopythonDeprecationWarning,
                    )
                    inplace = True
                else:
                    inplace = False
                if isinstance(self._data, _PartiallyDefinedSequenceData):
                    for seq in self._data._data.values():
                        if b"U" in seq or b"u" in seq:
                            warnings.warn(
                                "seq.reverse_complement() will change in the near "
                                "future to always return DNA nucleotides only. "
                                "Please use\n"
                                "\n"
                                "seq.reverse_complement_rna()\n"
                                "\n"
                                "if you want to receive an RNA sequence instead.",
                                BiopythonDeprecationWarning,
                            )
                            for seq in self._data._data.values():
                                if b"t" in seq or b"T" in seq:
                                    raise ValueError("Mixed RNA/DNA found")
                            return self.reverse_complement_rna(inplace=inplace)
                elif b"U" in self._data or b"u" in self._data:
                    warnings.warn(
                        "seq.reverse_complement() will change in the near "
                        "future to always return DNA nucleotides only. "
                        "Please use\n"
                        "\n"
                        "seq.reverse_complement_rna()\n"
                        "\n"
                        "if you want to receive an RNA sequence instead.",
                        BiopythonDeprecationWarning,
                    )
                    if b"t" in self._data or b"T" in self._data:
                        raise ValueError("Mixed RNA/DNA found")
                    return self.reverse_complement_rna(inplace=inplace)
            data = self._data.translate(_dna_complement_table)
        except UndefinedSequenceError:
            # reverse complement of an undefined sequence is an undefined sequence
            # of the same length
            return self
        if inplace:
            if not isinstance(self._data, bytearray):
                raise TypeError("Sequence is immutable")
            self._data[::-1] = data
            return self
        return self.__class__(data[::-1])

    def reverse_complement_rna(self, inplace=False):
        """Return the reverse complement as an RNA sequence.

        >>> Seq("CGA").reverse_complement_rna()
        Seq('UCG')

        Any T in the sequence is treated as a U:

        >>> Seq("CGAUT").reverse_complement_rna()
        Seq('AAUCG')

        In contrast, ``reverse_complement`` returns a DNA sequence:

        >>> Seq("CGA").reverse_complement(inplace=False)
        Seq('TCG')

        The sequence is modified in-place and returned if inplace is True:

        >>> my_seq = MutableSeq("CGA")
        >>> my_seq
        MutableSeq('CGA')
        >>> my_seq.reverse_complement_rna()
        MutableSeq('UCG')
        >>> my_seq
        MutableSeq('CGA')

        >>> my_seq.reverse_complement_rna(inplace=True)
        MutableSeq('UCG')
        >>> my_seq
        MutableSeq('UCG')

        As ``Seq`` objects are immutable, a ``TypeError`` is raised if
        ``reverse_complement_rna`` is called on a ``Seq`` object with
        ``inplace=True``.
        """
        try:
            data = self._data.translate(_rna_complement_table)
        except UndefinedSequenceError:
            # reverse complement of an undefined sequence is an undefined sequence
            # of the same length
            return self
        if inplace:
            if not isinstance(self._data, bytearray):
                raise TypeError("Sequence is immutable")
            self._data[::-1] = data
            return self
        return self.__class__(data[::-1])

    def transcribe(self, inplace=False):
        """Transcribe a DNA sequence into RNA and return the RNA sequence as a new Seq object.

        >>> from Bio.Seq import Seq
        >>> coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
        >>> coding_dna
        Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')
        >>> coding_dna.transcribe()
        Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')

        The sequence is modified in-place and returned if inplace is True:

        >>> sequence = MutableSeq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
        >>> sequence
        MutableSeq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')
        >>> sequence.transcribe()
        MutableSeq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')
        >>> sequence
        MutableSeq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')

        >>> sequence.transcribe(inplace=True)
        MutableSeq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')
        >>> sequence
        MutableSeq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')

        As ``Seq`` objects are immutable, a ``TypeError`` is raised if
        ``transcribe`` is called on a ``Seq`` object with ``inplace=True``.

        Trying to transcribe an RNA sequence has no effect.
        If you have a nucleotide sequence which might be DNA or RNA
        (or even a mixture), calling the transcribe method will ensure
        any T becomes U.

        Trying to transcribe a protein sequence will replace any
        T for Threonine with U for Selenocysteine, which has no
        biologically plausible rational.

        >>> from Bio.Seq import Seq
        >>> my_protein = Seq("MAIVMGRT")
        >>> my_protein.transcribe()
        Seq('MAIVMGRU')
        """
        data = self._data.replace(b"T", b"U").replace(b"t", b"u")
        if inplace:
            if not isinstance(self._data, bytearray):
                raise TypeError("Sequence is immutable")
            self._data[:] = data
            return self
        return self.__class__(data)

    def back_transcribe(self, inplace=False):
        """Return the DNA sequence from an RNA sequence by creating a new Seq object.

        >>> from Bio.Seq import Seq
        >>> messenger_rna = Seq("AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG")
        >>> messenger_rna
        Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')
        >>> messenger_rna.back_transcribe()
        Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')

        The sequence is modified in-place and returned if inplace is True:

        >>> sequence = MutableSeq("AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG")
        >>> sequence
        MutableSeq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')
        >>> sequence.back_transcribe()
        MutableSeq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')
        >>> sequence
        MutableSeq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')

        >>> sequence.back_transcribe(inplace=True)
        MutableSeq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')
        >>> sequence
        MutableSeq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')

        As ``Seq`` objects are immutable, a ``TypeError`` is raised if
        ``transcribe`` is called on a ``Seq`` object with ``inplace=True``.

        Trying to back-transcribe DNA has no effect, If you have a nucleotide
        sequence which might be DNA or RNA (or even a mixture), calling the
        back-transcribe method will ensure any U becomes T.

        Trying to back-transcribe a protein sequence will replace any U for
        Selenocysteine with T for Threonine, which is biologically meaningless.

        >>> from Bio.Seq import Seq
        >>> my_protein = Seq("MAIVMGRU")
        >>> my_protein.back_transcribe()
        Seq('MAIVMGRT')
        """
        data = self._data.replace(b"U", b"T").replace(b"u", b"t")
        if inplace:
            if not isinstance(self._data, bytearray):
                raise TypeError("Sequence is immutable")
            self._data[:] = data
            return self
        return self.__class__(data)

    def join(self, other):
        """Return a merge of the sequences in other, spaced by the sequence from self.

        Accepts a Seq object, MutableSeq object, or string (and iterates over
        the letters), or an iterable containing Seq, MutableSeq, or string
        objects. These arguments will be concatenated with the calling sequence
        as the spacer:

        >>> concatenated = Seq('NNNNN').join([Seq("AAA"), Seq("TTT"), Seq("PPP")])
        >>> concatenated
        Seq('AAANNNNNTTTNNNNNPPP')

        Joining the letters of a single sequence:

        >>> Seq('NNNNN').join(Seq("ACGT"))
        Seq('ANNNNNCNNNNNGNNNNNT')
        >>> Seq('NNNNN').join("ACGT")
        Seq('ANNNNNCNNNNNGNNNNNT')
        """
        if isinstance(other, _SeqAbstractBaseClass):
            return self.__class__(str(self).join(str(other)))
        elif isinstance(other, str):
            return self.__class__(str(self).join(other))

        from Bio.SeqRecord import SeqRecord  # Lazy to avoid circular imports

        if isinstance(other, SeqRecord):
            raise TypeError("Iterable cannot be a SeqRecord")

        for c in other:
            if isinstance(c, SeqRecord):
                raise TypeError("Iterable cannot contain SeqRecords")
            elif not isinstance(c, (str, _SeqAbstractBaseClass)):
                raise TypeError(
                    "Input must be an iterable of Seq objects, MutableSeq objects, or strings"
                )
        return self.__class__(str(self).join([str(_) for _ in other]))

    def replace(self, old, new, inplace=False):
        """Return a copy with all occurrences of subsequence old replaced by new.

        >>> s = Seq("ACGTAACCGGTT")
        >>> t = s.replace("AC", "XYZ")
        >>> s
        Seq('ACGTAACCGGTT')
        >>> t
        Seq('XYZGTAXYZCGGTT')

        For mutable sequences, passing inplace=True will modify the sequence in place:

        >>> m = MutableSeq("ACGTAACCGGTT")
        >>> t = m.replace("AC", "XYZ")
        >>> m
        MutableSeq('ACGTAACCGGTT')
        >>> t
        MutableSeq('XYZGTAXYZCGGTT')

        >>> m = MutableSeq("ACGTAACCGGTT")
        >>> t = m.replace("AC", "XYZ", inplace=True)
        >>> m
        MutableSeq('XYZGTAXYZCGGTT')
        >>> t
        MutableSeq('XYZGTAXYZCGGTT')

        As ``Seq`` objects are immutable, a ``TypeError`` is raised if
        ``replace`` is called on a ``Seq`` object with ``inplace=True``.
        """
        if isinstance(old, _SeqAbstractBaseClass):
            old = bytes(old)
        elif isinstance(old, str):
            old = old.encode("ASCII")
        if isinstance(new, _SeqAbstractBaseClass):
            new = bytes(new)
        elif isinstance(new, str):
            new = new.encode("ASCII")
        data = self._data.replace(old, new)
        if inplace:
            if not isinstance(self._data, bytearray):
                raise TypeError("Sequence is immutable")
            self._data[:] = data
            return self
        return self.__class__(data)

    @property
    def defined(self):
        """Return True if the sequence is defined, False if undefined or partially defined.

        Zero-length sequences are always considered to be defined.
        """
        if isinstance(self._data, (bytes, bytearray)):
            return True
        else:
            return self._data.defined

    @property
    def defined_ranges(self):
        """Return a tuple of the ranges where the sequence contents is defined.

        The return value has the format ((start1, end1), (start2, end2), ...).
        """
        if isinstance(self._data, (bytes, bytearray)):
            length = len(self)
            if length > 0:
                return ((0, length),)
            else:
                return ()
        else:
            return self._data.defined_ranges


class Seq(_SeqAbstractBaseClass):
    """Read-only sequence object (essentially a string with biological methods).

    Like normal python strings, our basic sequence object is immutable.
    This prevents you from doing my_seq[5] = "A" for example, but does allow
    Seq objects to be used as dictionary keys.

    The Seq object provides a number of string like methods (such as count,
    find, split and strip).

    The Seq object also provides some biological methods, such as complement,
    reverse_complement, transcribe, back_transcribe and translate (which are
    not applicable to protein sequences).
    """

    def __init__(self, data, length=None):
        """Create a Seq object.

        Arguments:
         - data - Sequence, required (string)
         - length - Sequence length, used only if data is None or a dictionary (integer)

        You will typically use Bio.SeqIO to read in sequences from files as
        SeqRecord objects, whose sequence will be exposed as a Seq object via
        the seq property.

        However, you can also create a Seq object directly:

        >>> from Bio.Seq import Seq
        >>> my_seq = Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF")
        >>> my_seq
        Seq('MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF')
        >>> print(my_seq)
        MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF

        To create a Seq object with for a sequence of known length but
        unknown sequence contents, use None for the data argument and pass
        the sequence length for the length argument. Trying to access the
        sequence contents of a Seq object created in this way will raise
        an UndefinedSequenceError:

        >>> my_undefined_sequence = Seq(None, 20)
        >>> my_undefined_sequence
        Seq(None, length=20)
        >>> len(my_undefined_sequence)
        20
        >>> print(my_undefined_sequence)
        Traceback (most recent call last):
        ...
        Bio.Seq.UndefinedSequenceError: Sequence content is undefined

        If the sequence contents is known for parts of the sequence only, use
        a dictionary for the data argument to pass the known sequence segments:

        >>> my_partially_defined_sequence = Seq({3: "ACGT"}, 10)
        >>> my_partially_defined_sequence
        Seq({3: 'ACGT'}, length=10)
        >>> len(my_partially_defined_sequence)
        10
        >>> print(my_partially_defined_sequence)
        Traceback (most recent call last):
        ...
        Bio.Seq.UndefinedSequenceError: Sequence content is only partially defined
        >>> my_partially_defined_sequence[3:7]
        Seq('ACGT')
        >>> print(my_partially_defined_sequence[3:7])
        ACGT
        """
        if data is None:
            if length is None:
                raise ValueError("length must not be None if data is None")
            elif length == 0:
                self._data = b""
            elif length < 0:
                raise ValueError("length must not be negative.")
            else:
                self._data = _UndefinedSequenceData(length)
        elif isinstance(data, (bytes, SequenceDataAbstractBaseClass)):
            self._data = data
        elif isinstance(data, (bytearray, _SeqAbstractBaseClass)):
            self._data = bytes(data)
        elif isinstance(data, str):
            self._data = bytes(data, encoding="ASCII")
        elif isinstance(data, dict):
            if length is None:
                raise ValueError("length must not be None if data is a dictionary")
            elif length == 0:
                self._data = b""
            elif length < 0:
                raise ValueError("length must not be negative.")
            else:
                end = -1
                starts = sorted(data.keys())
                _data = {}
                for start in starts:
                    seq = data[start]
                    if isinstance(seq, str):
                        seq = bytes(seq, encoding="ASCII")
                    else:
                        try:
                            seq = bytes(seq)
                        except Exception:
                            raise ValueError("Expected bytes-like objects or strings")
                    if start < end:
                        raise ValueError("Sequence data are overlapping.")
                    elif start == end:
                        _data[current] += seq  # noqa: F821
                    else:
                        _data[start] = seq
                        current = start
                    end = start + len(seq)
                if end > length:
                    raise ValueError(
                        "Provided sequence data extend beyond sequence length."
                    )
                elif end == length and current == 0:
                    # sequence is fully defined
                    self._data = _data[current]
                else:
                    self._data = _PartiallyDefinedSequenceData(length, _data)
        else:
            raise TypeError(
                "data should be a string, bytes, bytearray, Seq, or MutableSeq object"
            )

    def __hash__(self):
        """Hash of the sequence as a string for comparison.

        See Seq object comparison documentation (method ``__eq__`` in
        particular) as this has changed in Biopython 1.65. Older versions
        would hash on object identity.
        """
        return hash(self._data)

    def ungap(self, gap="-"):
        """Return a copy of the sequence without the gap character(s) (DEPRECATED).

        The gap character now defaults to the minus sign, and can only
        be specified via the method argument. This is no longer possible
        via the sequence's alphabet (as was possible up to Biopython 1.77):

        >>> from Bio.Seq import Seq
        >>> my_dna = Seq("-ATA--TGAAAT-TTGAAAA")
        >>> my_dna
        Seq('-ATA--TGAAAT-TTGAAAA')
        >>> my_dna.ungap("-")
        Seq('ATATGAAATTTGAAAA')

        This method is DEPRECATED; please use my_dna.replace(gap, "") instead.
        """
        warnings.warn(
            """\
myseq.ungap(gap) is deprecated; please use myseq.replace(gap, "") instead.""",
            BiopythonDeprecationWarning,
        )
        if not gap:
            raise ValueError("Gap character required.")
        elif len(gap) != 1 or not isinstance(gap, str):
            raise ValueError(f"Unexpected gap character, {gap!r}")
        return self.replace(gap, b"")


class MutableSeq(_SeqAbstractBaseClass):
    """An editable sequence object.

    Unlike normal python strings and our basic sequence object (the Seq class)
    which are immutable, the MutableSeq lets you edit the sequence in place.
    However, this means you cannot use a MutableSeq object as a dictionary key.

    >>> from Bio.Seq import MutableSeq
    >>> my_seq = MutableSeq("ACTCGTCGTCG")
    >>> my_seq
    MutableSeq('ACTCGTCGTCG')
    >>> my_seq[5]
    'T'
    >>> my_seq[5] = "A"
    >>> my_seq
    MutableSeq('ACTCGACGTCG')
    >>> my_seq[5]
    'A'
    >>> my_seq[5:8] = "NNN"
    >>> my_seq
    MutableSeq('ACTCGNNNTCG')
    >>> len(my_seq)
    11

    Note that the MutableSeq object does not support as many string-like
    or biological methods as the Seq object.
    """

    def __init__(self, data):
        """Create a MutableSeq object."""
        if isinstance(data, bytearray):
            self._data = data
        elif isinstance(data, bytes):
            self._data = bytearray(data)
        elif isinstance(data, str):
            self._data = bytearray(data, "ASCII")
        elif isinstance(data, MutableSeq):
            self._data = data._data[:]  # Take a copy
        elif isinstance(data, Seq):
            # Make no assumptions about the Seq subclass internal storage
            self._data = bytearray(bytes(data))
        else:
            raise TypeError(
                "data should be a string, bytearray object, Seq object, or a "
                "MutableSeq object"
            )

    def __setitem__(self, index, value):
        """Set a subsequence of single letter via value parameter.

        >>> my_seq = MutableSeq('ACTCGACGTCG')
        >>> my_seq[0] = 'T'
        >>> my_seq
        MutableSeq('TCTCGACGTCG')
        """
        if isinstance(index, numbers.Integral):
            # Replacing a single letter with a new string
            self._data[index] = ord(value)
        else:
            # Replacing a sub-sequence
            if isinstance(value, MutableSeq):
                self._data[index] = value._data
            elif isinstance(value, Seq):
                self._data[index] = bytes(value)
            elif isinstance(value, str):
                self._data[index] = value.encode("ASCII")
            else:
                raise TypeError(f"received unexpected type '{type(value).__name__}'")

    def __delitem__(self, index):
        """Delete a subsequence of single letter.

        >>> my_seq = MutableSeq('ACTCGACGTCG')
        >>> del my_seq[0]
        >>> my_seq
        MutableSeq('CTCGACGTCG')
        """
        # Could be deleting a single letter, or a slice
        del self._data[index]

    def append(self, c):
        """Add a subsequence to the mutable sequence object.

        >>> my_seq = MutableSeq('ACTCGACGTCG')
        >>> my_seq.append('A')
        >>> my_seq
        MutableSeq('ACTCGACGTCGA')

        No return value.
        """
        self._data.append(ord(c.encode("ASCII")))

    def insert(self, i, c):
        """Add a subsequence to the mutable sequence object at a given index.

        >>> my_seq = MutableSeq('ACTCGACGTCG')
        >>> my_seq.insert(0,'A')
        >>> my_seq
        MutableSeq('AACTCGACGTCG')
        >>> my_seq.insert(8,'G')
        >>> my_seq
        MutableSeq('AACTCGACGGTCG')

        No return value.
        """
        self._data.insert(i, ord(c.encode("ASCII")))

    def pop(self, i=(-1)):
        """Remove a subsequence of a single letter at given index.

        >>> my_seq = MutableSeq('ACTCGACGTCG')
        >>> my_seq.pop()
        'G'
        >>> my_seq
        MutableSeq('ACTCGACGTC')
        >>> my_seq.pop()
        'C'
        >>> my_seq
        MutableSeq('ACTCGACGT')

        Returns the last character of the sequence.
        """
        c = self._data[i]
        del self._data[i]
        return chr(c)

    def remove(self, item):
        """Remove a subsequence of a single letter from mutable sequence.

        >>> my_seq = MutableSeq('ACTCGACGTCG')
        >>> my_seq.remove('C')
        >>> my_seq
        MutableSeq('ATCGACGTCG')
        >>> my_seq.remove('A')
        >>> my_seq
        MutableSeq('TCGACGTCG')

        No return value.
        """
        codepoint = ord(item)
        try:
            self._data.remove(codepoint)
        except ValueError:
            raise ValueError("value not found in MutableSeq") from None

    def reverse(self):
        """Modify the mutable sequence to reverse itself.

        No return value.
        """
        self._data.reverse()

    def extend(self, other):
        """Add a sequence to the original mutable sequence object.

        >>> my_seq = MutableSeq('ACTCGACGTCG')
        >>> my_seq.extend('A')
        >>> my_seq
        MutableSeq('ACTCGACGTCGA')
        >>> my_seq.extend('TTT')
        >>> my_seq
        MutableSeq('ACTCGACGTCGATTT')

        No return value.
        """
        if isinstance(other, MutableSeq):
            self._data.extend(other._data)
        elif isinstance(other, Seq):
            self._data.extend(bytes(other))
        elif isinstance(other, str):
            self._data.extend(other.encode("ASCII"))
        else:
            raise TypeError("expected a string, Seq or MutableSeq")


class UndefinedSequenceError(ValueError):
    """Sequence contents is undefined."""


class _UndefinedSequenceData(SequenceDataAbstractBaseClass):
    """Stores the length of a sequence with an undefined sequence contents (PRIVATE).

    Objects of this class can be used to create a Seq object to represent
    sequences with a known length, but an unknown sequence contents.
    Calling __len__ returns the sequence length, calling __getitem__ raises an
    UndefinedSequenceError except for requests of zero size, for which it
    returns an empty bytes object.
    """

    __slots__ = ("_length",)

    def __init__(self, length):
        """Initialize the object with the sequence length.

        The calling function is responsible for ensuring that the length is
        greater than zero.
        """
        self._length = length
        super().__init__()

    def __getitem__(self, key):
        if isinstance(key, slice):
            start, end, step = key.indices(self._length)
            size = len(range(start, end, step))
            if size == 0:
                return b""
            return _UndefinedSequenceData(size)
        else:
            raise UndefinedSequenceError("Sequence content is undefined")

    def __len__(self):
        return self._length

    def __bytes__(self):
        raise UndefinedSequenceError("Sequence content is undefined")

    def __add__(self, other):
        length = len(self) + len(other)
        try:
            other = bytes(other)
        except UndefinedSequenceError:
            if isinstance(other, _UndefinedSequenceData):
                return _UndefinedSequenceData(length)
            else:
                return NotImplemented
                # _PartiallyDefinedSequenceData.__radd__ will handle this
        else:
            data = {len(self): other}
            return _PartiallyDefinedSequenceData(length, data)

    def __radd__(self, other):
        data = {0: bytes(other)}
        length = len(other) + len(self)
        return _PartiallyDefinedSequenceData(length, data)

    def upper(self):
        """Return an upper case copy of the sequence."""
        # An upper case copy of an undefined sequence is an undefined
        # sequence of the same length
        return _UndefinedSequenceData(self._length)

    def lower(self):
        """Return a lower case copy of the sequence."""
        # A lower case copy of an undefined sequence is an undefined
        # sequence of the same length
        return _UndefinedSequenceData(self._length)

    def isupper(self):
        """Return True if all ASCII characters in data are uppercase.

        If there are no cased characters, the method returns False.
        """
        # Character case is irrelevant for an undefined sequence
        raise UndefinedSequenceError("Sequence content is undefined")

    def islower(self):
        """Return True if all ASCII characters in data are lowercase.

        If there are no cased characters, the method returns False.
        """
        # Character case is irrelevant for an undefined sequence
        raise UndefinedSequenceError("Sequence content is undefined")

    def replace(self, old, new):
        """Return a copy with all occurrences of substring old replaced by new."""
        # Replacing substring old by new in an undefined sequence will result
        # in an undefined sequence of the same length, if old and new have the
        # number of characters.
        if len(old) != len(new):
            raise UndefinedSequenceError("Sequence content is undefined")
        return _UndefinedSequenceData(self._length)

    @property
    def defined(self):
        """Return False, as the sequence is not defined and has a non-zero length."""
        return False

    @property
    def defined_ranges(self):
        """Return a tuple of the ranges where the sequence contents is defined.

        As the sequence contents of an _UndefinedSequenceData object is fully
        undefined, the return value is always an empty tuple.
        """
        return ()


class _PartiallyDefinedSequenceData(SequenceDataAbstractBaseClass):
    """Stores the length of a sequence with an undefined sequence contents (PRIVATE).

    Objects of this class can be used to create a Seq object to represent
    sequences with a known length, but with a sequence contents that is only
    partially known.
    Calling __len__ returns the sequence length, calling __getitem__ returns
    the sequence contents if known, otherwise an UndefinedSequenceError is
    raised.
    """

    __slots__ = ("_length", "_data")

    def __init__(self, length, data):
        """Initialize with the sequence length and defined sequence segments.

        The calling function is responsible for ensuring that the length is
        greater than zero.
        """
        self._length = length
        self._data = data
        super().__init__()

    def __getitem__(self, key):
        if isinstance(key, slice):
            start, end, step = key.indices(self._length)
            size = len(range(start, end, step))
            if size == 0:
                return b""
            data = {}
            for s, d in self._data.items():
                indices = range(-s, -s + self._length)[key]
                e = indices.stop
                if step > 0:
                    if e <= 0:
                        continue
                    if indices.start < 0:
                        s = indices.start % step
                    else:
                        s = indices.start
                else:  # step < 0
                    if e < 0:
                        e = None
                    end = len(d) - 1
                    if indices.start > end:
                        s = end + (indices.start - end) % step
                    else:
                        s = indices.start
                    if s < 0:
                        continue
                start = (s - indices.start) // step
                d = d[s:e:step]
                if d:
                    data[start] = d
            if len(data) == 0:  # Fully undefined sequence
                return _UndefinedSequenceData(size)
            # merge adjacent sequence segments
            end = -1
            previous = None  # not needed here, but it keeps flake happy
            items = data.items()
            data = {}
            for start, seq in items:
                if end == start:
                    data[previous] += seq
                else:
                    data[start] = seq
                    previous = start
                end = start + len(seq)
            if len(data) == 1:
                seq = data.get(0)
                if seq is not None and len(seq) == size:
                    return seq  # Fully defined sequence; return bytes
            if step < 0:
                # use this after we drop Python 3.7:
                # data = {start: data[start] for start in reversed(data)}
                # use this as long as we support Python 3.7:
                data = {start: data[start] for start in reversed(list(data.keys()))}
            return _PartiallyDefinedSequenceData(size, data)
        elif self._length <= key:
            raise IndexError("sequence index out of range")
        else:
            for start, seq in self._data.items():
                if start <= key and key < start + len(seq):
                    return seq[key - start]
            raise UndefinedSequenceError("Sequence at position %d is undefined" % key)

    def __len__(self):
        return self._length

    def __bytes__(self):
        raise UndefinedSequenceError("Sequence content is only partially defined")

    def __add__(self, other):
        length = len(self) + len(other)
        data = dict(self._data)
        items = list(self._data.items())
        start, seq = items[-1]
        end = start + len(seq)
        try:
            other = bytes(other)
        except UndefinedSequenceError:
            if isinstance(other, _UndefinedSequenceData):
                pass
            elif isinstance(other, _PartiallyDefinedSequenceData):
                other_items = list(other._data.items())
                if end == len(self):
                    other_start, other_seq = other_items.pop(0)
                    if other_start == 0:
                        data[start] += other_seq
                    else:
                        data[len(self) + other_start] = other_seq
                for other_start, other_seq in other_items:
                    data[len(self) + other_start] = other_seq
        else:
            if end == len(self):
                data[start] += other
            else:
                data[len(self)] = other
        return _PartiallyDefinedSequenceData(length, data)

    def __radd__(self, other):
        length = len(other) + len(self)
        try:
            other = bytes(other)
        except UndefinedSequenceError:
            data = {len(other) + start: seq for start, seq in self._data.items()}
        else:
            data = {0: other}
            items = list(self._data.items())
            start, seq = items.pop(0)
            if start == 0:
                data[0] += seq
            else:
                data[len(other) + start] = seq
            for start, seq in items:
                data[len(other) + start] = seq
        return _PartiallyDefinedSequenceData(length, data)

    def __mul__(self, other):
        length = self._length
        items = self._data.items()
        data = {}
        end = -1
        previous = None  # not needed here, but it keeps flake happy
        for i in range(other):
            for start, seq in items:
                start += i * length
                if end == start:
                    data[previous] += seq
                else:
                    data[start] = seq
                    previous = start
            end = start + len(seq)
        return _PartiallyDefinedSequenceData(length * other, data)

    def upper(self):
        """Return an upper case copy of the sequence."""
        data = {start: seq.upper() for start, seq in self._data.items()}
        return _PartiallyDefinedSequenceData(self._length, data)

    def lower(self):
        """Return a lower case copy of the sequence."""
        data = {start: seq.lower() for start, seq in self._data.items()}
        return _PartiallyDefinedSequenceData(self._length, data)

    def isupper(self):
        """Return True if all ASCII characters in data are uppercase.

        If there are no cased characters, the method returns False.
        """
        # Character case is irrelevant for an undefined sequence
        raise UndefinedSequenceError("Sequence content is only partially defined")

    def islower(self):
        """Return True if all ASCII characters in data are lowercase.

        If there are no cased characters, the method returns False.
        """
        # Character case is irrelevant for an undefined sequence
        raise UndefinedSequenceError("Sequence content is only partially defined")

    def translate(self, table, delete=b""):
        """Return a copy with each character mapped by the given translation table.

          table
            Translation table, which must be a bytes object of length 256.

        All characters occurring in the optional argument delete are removed.
        The remaining characters are mapped through the given translation table.
        """
        items = self._data.items()
        data = {start: seq.translate(table, delete) for start, seq in items}
        return _PartiallyDefinedSequenceData(self._length, data)

    def replace(self, old, new):
        """Return a copy with all occurrences of substring old replaced by new."""
        # Replacing substring old by new in the undefined sequence segments
        # will result in an undefined sequence segment of the same length, if
        # old and new have the number of characters. If not, an error is raised,
        # as the correct start positions cannot be calculated reliably.
        if len(old) != len(new):
            raise UndefinedSequenceError(
                "Sequence content is only partially defined; substring \n"
                "replacement cannot be performed reliably"
            )
        items = self._data.items()
        data = {start: seq.replace(old, new) for start, seq in items}
        return _PartiallyDefinedSequenceData(self._length, data)

    @property
    def defined(self):
        """Return False, as the sequence is not fully defined and has a non-zero length."""
        return False

    @property
    def defined_ranges(self):
        """Return a tuple of the ranges where the sequence contents is defined.

        The return value has the format ((start1, end1), (start2, end2), ...).
        """
        return tuple((start, start + len(seq)) for start, seq in self._data.items())


# The transcribe, backward_transcribe, and translate functions are
# user-friendly versions of the corresponding Seq/MutableSeq methods.
# The functions work both on Seq objects, and on strings.


def transcribe(dna):
    """Transcribe a DNA sequence into RNA.

    If given a string, returns a new string object.

    Given a Seq or MutableSeq, returns a new Seq object.

    e.g.

    >>> transcribe("ACTGN")
    'ACUGN'
    """
    if isinstance(dna, Seq):
        return dna.transcribe()
    elif isinstance(dna, MutableSeq):
        return Seq(dna).transcribe()
    else:
        return dna.replace("T", "U").replace("t", "u")


def back_transcribe(rna):
    """Return the RNA sequence back-transcribed into DNA.

    If given a string, returns a new string object.

    Given a Seq or MutableSeq, returns a new Seq object.

    e.g.

    >>> back_transcribe("ACUGN")
    'ACTGN'
    """
    if isinstance(rna, Seq):
        return rna.back_transcribe()
    elif isinstance(rna, MutableSeq):
        return Seq(rna).back_transcribe()
    else:
        return rna.replace("U", "T").replace("u", "t")


def _translate_str(
    sequence, table, stop_symbol="*", to_stop=False, cds=False, pos_stop="X", gap=None
):
    """Translate nucleotide string into a protein string (PRIVATE).

    Arguments:
     - sequence - a string
     - table - Which codon table to use?  This can be either a name (string),
       an NCBI identifier (integer), or a CodonTable object (useful for
       non-standard genetic codes).  This defaults to the "Standard" table.
     - stop_symbol - a single character string, what to use for terminators.
     - to_stop - boolean, should translation terminate at the first
       in frame stop codon?  If there is no in-frame stop codon
       then translation continues to the end.
     - pos_stop - a single character string for a possible stop codon
       (e.g. TAN or NNN)
     - cds - Boolean, indicates this is a complete CDS.  If True, this
       checks the sequence starts with a valid alternative start
       codon (which will be translated as methionine, M), that the
       sequence length is a multiple of three, and that there is a
       single in frame stop codon at the end (this will be excluded
       from the protein sequence, regardless of the to_stop option).
       If these tests fail, an exception is raised.
     - gap - Single character string to denote symbol used for gaps.
       Defaults to None.

    Returns a string.

    e.g.

    >>> from Bio.Data import CodonTable
    >>> table = CodonTable.ambiguous_dna_by_id[1]
    >>> _translate_str("AAA", table)
    'K'
    >>> _translate_str("TAR", table)
    '*'
    >>> _translate_str("TAN", table)
    'X'
    >>> _translate_str("TAN", table, pos_stop="@")
    '@'
    >>> _translate_str("TA?", table)
    Traceback (most recent call last):
       ...
    Bio.Data.CodonTable.TranslationError: Codon 'TA?' is invalid

    In a change to older versions of Biopython, partial codons are now
    always regarded as an error (previously only checked if cds=True)
    and will trigger a warning (likely to become an exception in a
    future release).

    If **cds=True**, the start and stop codons are checked, and the start
    codon will be translated at methionine. The sequence must be an
    while number of codons.

    >>> _translate_str("ATGCCCTAG", table, cds=True)
    'MP'
    >>> _translate_str("AAACCCTAG", table, cds=True)
    Traceback (most recent call last):
       ...
    Bio.Data.CodonTable.TranslationError: First codon 'AAA' is not a start codon
    >>> _translate_str("ATGCCCTAGCCCTAG", table, cds=True)
    Traceback (most recent call last):
       ...
    Bio.Data.CodonTable.TranslationError: Extra in frame stop codon 'TAG' found.
    """
    try:
        table_id = int(table)
    except ValueError:
        # Assume it's a table name
        # The same table can be used for RNA or DNA
        try:
            codon_table = CodonTable.ambiguous_generic_by_name[table]
        except KeyError:
            if isinstance(table, str):
                raise ValueError(
                    "The Bio.Seq translate methods and function DO NOT "
                    "take a character string mapping table like the python "
                    "string object's translate method. "
                    "Use str(my_seq).translate(...) instead."
                ) from None
            else:
                raise TypeError("table argument must be integer or string") from None
    except (AttributeError, TypeError):
        # Assume it's a CodonTable object
        if isinstance(table, CodonTable.CodonTable):
            codon_table = table
        else:
            raise ValueError("Bad table argument") from None
    else:
        # Assume it's a table ID
        # The same table can be used for RNA or DNA
        codon_table = CodonTable.ambiguous_generic_by_id[table_id]
    sequence = sequence.upper()
    amino_acids = []
    forward_table = codon_table.forward_table
    stop_codons = codon_table.stop_codons
    if codon_table.nucleotide_alphabet is not None:
        valid_letters = set(codon_table.nucleotide_alphabet.upper())
    else:
        # Assume the worst case, ambiguous DNA or RNA:
        valid_letters = set(
            IUPACData.ambiguous_dna_letters.upper()
            + IUPACData.ambiguous_rna_letters.upper()
        )
    n = len(sequence)

    # Check for tables with 'ambiguous' (dual-coding) stop codons:
    dual_coding = [c for c in stop_codons if c in forward_table]
    if dual_coding:
        c = dual_coding[0]
        if to_stop:
            raise ValueError(
                "You cannot use 'to_stop=True' with this table as it contains"
                f" {len(dual_coding)} codon(s) which can be both STOP and an"
                f" amino acid (e.g. '{c}' -> '{forward_table[c]}' or STOP)."
            )
        warnings.warn(
            f"This table contains {len(dual_coding)} codon(s) which code(s) for"
            f" both STOP and an amino acid (e.g. '{c}' -> '{forward_table[c]}'"
            " or STOP). Such codons will be translated as amino acid.",
            BiopythonWarning,
        )

    if cds:
        if str(sequence[:3]).upper() not in codon_table.start_codons:
            raise CodonTable.TranslationError(
                f"First codon '{sequence[:3]}' is not a start codon"
            )
        if n % 3 != 0:
            raise CodonTable.TranslationError(
                f"Sequence length {n} is not a multiple of three"
            )
        if str(sequence[-3:]).upper() not in stop_codons:
            raise CodonTable.TranslationError(
                f"Final codon '{sequence[-3:]}' is not a stop codon"
            )
        # Don't translate the stop symbol, and manually translate the M
        sequence = sequence[3:-3]
        n -= 6
        amino_acids = ["M"]
    elif n % 3 != 0:
        warnings.warn(
            "Partial codon, len(sequence) not a multiple of three. "
            "Explicitly trim the sequence or add trailing N before "
            "translation. This may become an error in future.",
            BiopythonWarning,
        )
    if gap is not None:
        if not isinstance(gap, str):
            raise TypeError("Gap character should be a single character string.")
        elif len(gap) > 1:
            raise ValueError("Gap character should be a single character string.")

    for i in range(0, n - n % 3, 3):
        codon = sequence[i : i + 3]
        try:
            amino_acids.append(forward_table[codon])
        except (KeyError, CodonTable.TranslationError):
            if codon in codon_table.stop_codons:
                if cds:
                    raise CodonTable.TranslationError(
                        f"Extra in frame stop codon '{codon}' found."
                    ) from None
                if to_stop:
                    break
                amino_acids.append(stop_symbol)
            elif valid_letters.issuperset(set(codon)):
                # Possible stop codon (e.g. NNN or TAN)
                amino_acids.append(pos_stop)
            elif gap is not None and codon == gap * 3:
                # Gapped translation
                amino_acids.append(gap)
            else:
                raise CodonTable.TranslationError(
                    f"Codon '{codon}' is invalid"
                ) from None
    return "".join(amino_acids)


def translate(
    sequence, table="Standard", stop_symbol="*", to_stop=False, cds=False, gap=None
):
    """Translate a nucleotide sequence into amino acids.

    If given a string, returns a new string object. Given a Seq or
    MutableSeq, returns a Seq object.

    Arguments:
     - table - Which codon table to use?  This can be either a name
       (string), an NCBI identifier (integer), or a CodonTable object
       (useful for non-standard genetic codes).  Defaults to the "Standard"
       table.
     - stop_symbol - Single character string, what to use for any
       terminators, defaults to the asterisk, "*".
     - to_stop - Boolean, defaults to False meaning do a full
       translation continuing on past any stop codons
       (translated as the specified stop_symbol).  If
       True, translation is terminated at the first in
       frame stop codon (and the stop_symbol is not
       appended to the returned protein sequence).
     - cds - Boolean, indicates this is a complete CDS.  If True, this
       checks the sequence starts with a valid alternative start
       codon (which will be translated as methionine, M), that the
       sequence length is a multiple of three, and that there is a
       single in frame stop codon at the end (this will be excluded
       from the protein sequence, regardless of the to_stop option).
       If these tests fail, an exception is raised.
     - gap - Single character string to denote symbol used for gaps.
       Defaults to None.

    A simple string example using the default (standard) genetic code:

    >>> coding_dna = "GTGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
    >>> translate(coding_dna)
    'VAIVMGR*KGAR*'
    >>> translate(coding_dna, stop_symbol="@")
    'VAIVMGR@KGAR@'
    >>> translate(coding_dna, to_stop=True)
    'VAIVMGR'

    Now using NCBI table 2, where TGA is not a stop codon:

    >>> translate(coding_dna, table=2)
    'VAIVMGRWKGAR*'
    >>> translate(coding_dna, table=2, to_stop=True)
    'VAIVMGRWKGAR'

    In fact this example uses an alternative start codon valid under NCBI
    table 2, GTG, which means this example is a complete valid CDS which
    when translated should really start with methionine (not valine):

    >>> translate(coding_dna, table=2, cds=True)
    'MAIVMGRWKGAR'

    Note that if the sequence has no in-frame stop codon, then the to_stop
    argument has no effect:

    >>> coding_dna2 = "GTGGCCATTGTAATGGGCCGC"
    >>> translate(coding_dna2)
    'VAIVMGR'
    >>> translate(coding_dna2, to_stop=True)
    'VAIVMGR'

    NOTE - Ambiguous codons like "TAN" or "NNN" could be an amino acid
    or a stop codon.  These are translated as "X".  Any invalid codon
    (e.g. "TA?" or "T-A") will throw a TranslationError.

    It will however translate either DNA or RNA.

    NOTE - Since version 1.71 Biopython contains codon tables with 'ambiguous
    stop codons'. These are stop codons with unambiguous sequence but which
    have a context dependent coding as STOP or as amino acid. With these tables
    'to_stop' must be False (otherwise a ValueError is raised). The dual
    coding codons will always be translated as amino acid, except for
    'cds=True', where the last codon will be translated as STOP.

    >>> coding_dna3 = "ATGGCACGGAAGTGA"
    >>> translate(coding_dna3)
    'MARK*'

    >>> translate(coding_dna3, table=27)  # Table 27: TGA -> STOP or W
    'MARKW'

    It will however raise a BiopythonWarning (not shown).

    >>> translate(coding_dna3, table=27, cds=True)
    'MARK'

    >>> translate(coding_dna3, table=27, to_stop=True)
    Traceback (most recent call last):
       ...
    ValueError: You cannot use 'to_stop=True' with this table ...
    """
    if isinstance(sequence, Seq):
        return sequence.translate(table, stop_symbol, to_stop, cds)
    elif isinstance(sequence, MutableSeq):
        # Return a Seq object
        return Seq(sequence).translate(table, stop_symbol, to_stop, cds)
    else:
        # Assume it's a string, return a string
        return _translate_str(sequence, table, stop_symbol, to_stop, cds, gap=gap)


def reverse_complement(sequence, inplace=None):
    """Return the reverse complement as a DNA sequence.

    If given a string, returns a new string object.
    Given a Seq object, returns a new Seq object.
    Given a MutableSeq, returns a new MutableSeq object.
    Given a SeqRecord object, returns a new SeqRecord object.

    >>> my_seq = "CGA"
    >>> reverse_complement(my_seq, inplace=False)
    'TCG'
    >>> my_seq = Seq("CGA")
    >>> reverse_complement(my_seq, inplace=False)
    Seq('TCG')
    >>> my_seq = MutableSeq("CGA")
    >>> reverse_complement(my_seq, inplace=False)
    MutableSeq('TCG')
    >>> my_seq
    MutableSeq('CGA')

    Any U in the sequence is treated as a T:

    >>> reverse_complement(Seq("CGAUT"), inplace=False)
    Seq('AATCG')

    In contrast, ``reverse_complement_rna`` returns an RNA sequence:

    >>> reverse_complement_rna(Seq("CGAUT"))
    Seq('AAUCG')

    Supports and lower- and upper-case characters, and unambiguous and
    ambiguous nucleotides. All other characters are not converted:

    >>> reverse_complement("ACGTUacgtuXYZxyz", inplace=False)
    'zrxZRXaacgtAACGT'

    The sequence is modified in-place and returned if inplace is True:

    >>> my_seq = MutableSeq("CGA")
    >>> reverse_complement(my_seq, inplace=True)
    MutableSeq('TCG')
    >>> my_seq
    MutableSeq('TCG')

    As strings and ``Seq`` objects are immutable, a ``TypeError`` is
    raised if ``reverse_complement`` is called on a ``Seq`` object with
    ``inplace=True``.
    """
    from Bio.SeqRecord import SeqRecord  # Lazy to avoid circular imports

    if inplace is None:
        # deprecated
        if isinstance(sequence, Seq):
            if b"U" in sequence._data or b"u" in sequence._data:
                warnings.warn(
                    "reverse_complement(sequence) will change in the "
                    "near future to always return DNA nucleotides only. "
                    "Please use\n"
                    "\n"
                    "reverse_complement_rna(sequence)\n"
                    "\n"
                    "if you want to receive an RNA sequence instead.",
                    BiopythonDeprecationWarning,
                )
                if b"T" in sequence._data or b"t" in sequence._data:
                    raise ValueError("Mixed RNA/DNA found")
                return sequence.reverse_complement_rna()
        elif isinstance(sequence, MutableSeq):
            # Return a Seq
            # Don't use the MutableSeq reverse_complement method as it is
            # 'in place'.
            warnings.warn(
                "reverse_complement(mutable_seq) will change in the near "
                "future to return a MutableSeq object instead of a Seq object.",
                BiopythonDeprecationWarning,
            )
            return Seq(sequence).reverse_complement()
        else:  # str
            if "U" in sequence or "u" in sequence:
                warnings.warn(
                    "reverse_complement(sequence) will change in the "
                    "near future to always return DNA nucleotides only. "
                    "Please use\n"
                    "\n"
                    "reverse_complement_rna(sequence)\n"
                    "\n"
                    "if you want to receive an RNA sequence instead.",
                    BiopythonDeprecationWarning,
                )
                if "T" in sequence or "t" in sequence:
                    raise ValueError("Mixed RNA/DNA found")
                sequence = sequence.encode("ASCII")
                sequence = sequence.translate(_rna_complement_table)
                return sequence.decode("ASCII")[::-1]
    if isinstance(sequence, (Seq, MutableSeq)):
        return sequence.reverse_complement(inplace)
    if isinstance(sequence, SeqRecord):
        if inplace:
            raise TypeError("SeqRecords are immutable")
        return sequence.reverse_complement()
    # Assume it's a string.
    if inplace:
        raise TypeError("strings are immutable")
    sequence = sequence.encode("ASCII")
    sequence = sequence.translate(_dna_complement_table)
    sequence = sequence.decode("ASCII")
    return sequence[::-1]


def reverse_complement_rna(sequence, inplace=False):
    """Return the reverse complement as an RNA sequence.

    If given a string, returns a new string object.
    Given a Seq object, returns a new Seq object.
    Given a MutableSeq, returns a new MutableSeq object.
    Given a SeqRecord object, returns a new SeqRecord object.

    >>> my_seq = "CGA"
    >>> reverse_complement_rna(my_seq)
    'UCG'
    >>> my_seq = Seq("CGA")
    >>> reverse_complement_rna(my_seq)
    Seq('UCG')
    >>> my_seq = MutableSeq("CGA")
    >>> reverse_complement_rna(my_seq)
    MutableSeq('UCG')
    >>> my_seq
    MutableSeq('CGA')

    Any T in the sequence is treated as a U:

    >>> reverse_complement_rna(Seq("CGAUT"))
    Seq('AAUCG')

    In contrast, ``reverse_complement`` returns a DNA sequence:

    >>> reverse_complement(Seq("CGAUT"), inplace=False)
    Seq('AATCG')

    Supports and lower- and upper-case characters, and unambiguous and
    ambiguous nucleotides. All other characters are not converted:

    >>> reverse_complement_rna("ACGTUacgtuXYZxyz")
    'zrxZRXaacguAACGU'

    The sequence is modified in-place and returned if inplace is True:

    >>> my_seq = MutableSeq("CGA")
    >>> reverse_complement_rna(my_seq, inplace=True)
    MutableSeq('UCG')
    >>> my_seq
    MutableSeq('UCG')

    As strings and ``Seq`` objects are immutable, a ``TypeError`` is
    raised if ``reverse_complement`` is called on a ``Seq`` object with
    ``inplace=True``.
    """
    from Bio.SeqRecord import SeqRecord  # Lazy to avoid circular imports

    if isinstance(sequence, (Seq, MutableSeq)):
        return sequence.reverse_complement_rna(inplace)
    if isinstance(sequence, SeqRecord):
        if inplace:
            raise TypeError("SeqRecords are immutable")
        return sequence.reverse_complement_rna()
    # Assume it's a string.
    if inplace:
        raise TypeError("strings are immutable")
    sequence = sequence.encode("ASCII")
    sequence = sequence.translate(_rna_complement_table)
    sequence = sequence.decode("ASCII")
    return sequence[::-1]


def complement(sequence, inplace=None):
    """Return the complement as a DNA sequence.

    If given a string, returns a new string object.
    Given a Seq object, returns a new Seq object.
    Given a MutableSeq, returns a new MutableSeq object.
    Given a SeqRecord object, returns a new SeqRecord object.

    >>> my_seq = "CGA"
    >>> complement(my_seq, inplace=False)
    'GCT'
    >>> my_seq = Seq("CGA")
    >>> complement(my_seq, inplace=False)
    Seq('GCT')
    >>> my_seq = MutableSeq("CGA")
    >>> complement(my_seq, inplace=False)
    MutableSeq('GCT')
    >>> my_seq
    MutableSeq('CGA')

    Any U in the sequence is treated as a T:

    >>> complement(Seq("CGAUT"), inplace=False)
    Seq('GCTAA')

    In contrast, ``complement_rna`` returns an RNA sequence:

    >>> complement_rna(Seq("CGAUT"))
    Seq('GCUAA')

    Supports and lower- and upper-case characters, and unambiguous and
    ambiguous nucleotides. All other characters are not converted:

    >>> complement("ACGTUacgtuXYZxyz", inplace=False)
    'TGCAAtgcaaXRZxrz'

    The sequence is modified in-place and returned if inplace is True:

    >>> my_seq = MutableSeq("CGA")
    >>> complement(my_seq, inplace=True)
    MutableSeq('GCT')
    >>> my_seq
    MutableSeq('GCT')

    As strings and ``Seq`` objects are immutable, a ``TypeError`` is
    raised if ``reverse_complement`` is called on a ``Seq`` object with
    ``inplace=True``.
    """
    from Bio.SeqRecord import SeqRecord  # Lazy to avoid circular imports

    if inplace is None:
        # deprecated
        if isinstance(sequence, Seq):
            # Return a Seq
            if b"U" in sequence._data or b"u" in sequence._data:
                warnings.warn(
                    "complement(sequence) will change in the near "
                    "future to always return DNA nucleotides only. "
                    "Please use\n"
                    "\n"
                    "complement_rna(sequence)\n"
                    "\n"
                    "if you want to receive an RNA sequence instead.",
                    BiopythonDeprecationWarning,
                )
                if b"T" in sequence._data or b"t" in sequence._data:
                    raise ValueError("Mixed RNA/DNA found")
                return sequence.complement_rna()
        elif isinstance(sequence, MutableSeq):
            # Return a Seq
            # Don't use the MutableSeq reverse_complement method as it is
            # 'in place'.
            warnings.warn(
                "complement(mutable_seq) will change in the near future"
                "to return a MutableSeq object instead of a Seq object.",
                BiopythonDeprecationWarning,
            )
            return Seq(sequence).complement()
        else:
            if "U" in sequence or "u" in sequence:
                warnings.warn(
                    "complement(sequence) will change in the near "
                    "future to always return DNA nucleotides only. "
                    "Please use\n"
                    "\n"
                    "complement_rna(sequence)\n"
                    "\n"
                    "if you want to receive an RNA sequence instead.",
                    BiopythonDeprecationWarning,
                )
                if "T" in sequence or "t" in sequence:
                    raise ValueError("Mixed RNA/DNA found")
                ttable = _rna_complement_table
                sequence = sequence.encode("ASCII")
                sequence = sequence.translate(ttable)
                return sequence.decode("ASCII")
    if isinstance(sequence, (Seq, MutableSeq)):
        return sequence.complement(inplace)
    if isinstance(sequence, SeqRecord):
        if inplace:
            raise TypeError("SeqRecords are immutable")
        return sequence.complement()
    # Assume it's a string.
    if inplace:
        raise TypeError("strings are immutable")
    sequence = sequence.encode("ASCII")
    sequence = sequence.translate(_dna_complement_table)
    return sequence.decode("ASCII")


def complement_rna(sequence, inplace=False):
    """Return the complement as an RNA sequence.

    If given a string, returns a new string object.
    Given a Seq object, returns a new Seq object.
    Given a MutableSeq, returns a new MutableSeq object.
    Given a SeqRecord object, returns a new SeqRecord object.

    >>> my_seq = "CGA"
    >>> complement_rna(my_seq)
    'GCU'
    >>> my_seq = Seq("CGA")
    >>> complement_rna(my_seq)
    Seq('GCU')
    >>> my_seq = MutableSeq("CGA")
    >>> complement_rna(my_seq)
    MutableSeq('GCU')
    >>> my_seq
    MutableSeq('CGA')

    Any T in the sequence is treated as a U:

    >>> complement_rna(Seq("CGAUT"))
    Seq('GCUAA')

    In contrast, ``complement`` returns a DNA sequence:

    >>> complement(Seq("CGAUT"),inplace=False)
    Seq('GCTAA')

    Supports and lower- and upper-case characters, and unambiguous and
    ambiguous nucleotides. All other characters are not converted:

    >>> complement_rna("ACGTUacgtuXYZxyz")
    'UGCAAugcaaXRZxrz'

    The sequence is modified in-place and returned if inplace is True:

    >>> my_seq = MutableSeq("CGA")
    >>> complement(my_seq, inplace=True)
    MutableSeq('GCT')
    >>> my_seq
    MutableSeq('GCT')

    As strings and ``Seq`` objects are immutable, a ``TypeError`` is
    raised if ``reverse_complement`` is called on a ``Seq`` object with
    ``inplace=True``.
    """
    from Bio.SeqRecord import SeqRecord  # Lazy to avoid circular imports

    if isinstance(sequence, (Seq, MutableSeq)):
        return sequence.complement_rna(inplace)
    if isinstance(sequence, SeqRecord):
        if inplace:
            raise TypeError("SeqRecords are immutable")
        return sequence.complement_rna()
    # Assume it's a string.
    if inplace:
        raise TypeError("strings are immutable")
    sequence = sequence.encode("ASCII")
    sequence = sequence.translate(_rna_complement_table)
    return sequence.decode("ASCII")


def _test():
    """Run the Bio.Seq module's doctests (PRIVATE)."""
    print("Running doctests...")
    import doctest

    doctest.testmod(optionflags=doctest.IGNORE_EXCEPTION_DETAIL)
    print("Done")


if __name__ == "__main__":
    _test()
