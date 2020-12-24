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
import warnings
from abc import ABC, abstractmethod

from Bio import BiopythonWarning, BiopythonDeprecationWarning
from Bio.Data import IUPACData, CodonTable


def _maketrans(complement_mapping):
    """Make a python string translation table (PRIVATE).

    Arguments:
     - complement_mapping - a dictionary such as ambiguous_dna_complement
       and ambiguous_rna_complement from Data.IUPACData.

    Returns a translation table (a string of length 256) for use with the
    python string's translate method to use in a (reverse) complement.

    Compatible with lower case and upper case sequences.

    For internal use only.
    """
    keys = "".join(complement_mapping.keys()).encode("ASCII")
    values = "".join(complement_mapping.values()).encode("ASCII")
    return bytes.maketrans(keys + keys.lower(), values + values.lower())


_dna_complement_table = _maketrans(IUPACData.ambiguous_dna_complement)
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
        return bytes(self) + other

    def __radd__(self, other):
        return other + bytes(self)

    def __mul__(self, other):
        return bytes(self) * other

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
        return bytes(self).translate(table)


class Seq:
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
         - length - Sequence length, used only if data is None (integer)

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

        >>> my_undefined_seq = Seq(None, 20)
        >>> my_undefined_seq
        Seq(None, length=20)
        >>> len(my_undefined_seq)
        20
        >>> print(my_undefined_seq)
        Traceback (most recent call last):
        ...
        Bio.Seq.UndefinedSequenceError: Sequence content is undefined
        """
        if length is None:
            if isinstance(data, (bytes, SequenceDataAbstractBaseClass)):
                self._data = data
            elif isinstance(data, (bytearray, Seq, MutableSeq)):
                self._data = bytes(data)
            elif isinstance(data, str):
                self._data = bytes(data, encoding="ASCII")
            else:
                raise TypeError(
                    "data should be a string, bytes, bytearray, Seq, or MutableSeq object"
                )
        else:
            if data is not None:
                raise ValueError("length should be None if data is None")
            self._data = _UndefinedSequenceData(length)

    def __bytes__(self):
        return bytes(self._data)

    def __repr__(self):
        """Return (truncated) representation of the sequence for debugging."""
        data = self._data
        if isinstance(data, _UndefinedSequenceData):
            return f"Seq(None, length={len(self)})"
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
        """Return the full sequence as a python string, use str(my_seq).

        Note that Biopython 1.44 and earlier would give a truncated
        version of repr(my_seq) for str(my_seq).  If you are writing code
        which need to be backwards compatible with really old Biopython,
        you should continue to use my_seq.tostring() as follows::

            try:
                # The old way, removed in Biopython 1.73
                as_string = seq_obj.tostring()
            except AttributeError:
                # The new way, needs Biopython 1.45 or later.
                # Don't use this on Biopython 1.44 or older as truncates
                as_string = str(seq_obj)

        """
        return self._data.decode("ASCII")

    def __hash__(self):
        """Hash of the sequence as a string for comparison.

        See Seq object comparison documentation (method ``__eq__`` in
        particular) as this has changed in Biopython 1.65. Older versions
        would hash on object identity.
        """
        return hash(self._data)

    def __eq__(self, other):
        """Compare the sequence to another sequence or a string (README).

        Historically comparing Seq objects has done Python object comparison.
        After considerable discussion (keeping in mind constraints of the
        Python language, hashes and dictionary support), Biopython now uses
        simple string comparison (with a warning about the change).

        If you still need to support releases prior to Biopython 1.65, please
        just do explicit comparisons:

        >>> from Bio.Seq import Seq
        >>> seq1 = Seq("ACGT")
        >>> seq2 = Seq("ACGT")
        >>> id(seq1) == id(seq2)
        False
        >>> str(seq1) == str(seq2)
        True

        The new behaviour is to use string-like equality:

        >>> from Bio.Seq import Seq
        >>> seq1 == seq2
        True
        >>> seq1 == "ACGT"
        True
        """
        if isinstance(other, (Seq, MutableSeq)):
            return self._data == other._data
        elif isinstance(other, str):
            return self._data == other.encode("ASCII")
        else:
            return self._data == other

    def __lt__(self, other):
        """Implement the less-than operand."""
        if isinstance(other, (Seq, MutableSeq)):
            return self._data < other._data
        elif isinstance(other, str):
            return self._data < other.encode("ASCII")
        else:
            return self._data < other

    def __le__(self, other):
        """Implement the less-than or equal operand."""
        if isinstance(other, (Seq, MutableSeq)):
            return self._data <= other._data
        elif isinstance(other, str):
            return self._data <= other.encode("ASCII")
        else:
            return self._data <= other

    def __gt__(self, other):
        """Implement the greater-than operand."""
        if isinstance(other, (Seq, MutableSeq)):
            return self._data > other._data
        elif isinstance(other, str):
            return self._data > other.encode("ASCII")
        else:
            return self._data > other

    def __ge__(self, other):
        """Implement the greater-than or equal operand."""
        if isinstance(other, (Seq, MutableSeq)):
            return self._data >= other._data
        elif isinstance(other, str):
            return self._data >= other.encode("ASCII")
        else:
            return self._data >= other

    def __len__(self):
        """Return the length of the sequence, use len(my_seq)."""
        return len(self._data)

    def __getitem__(self, index):
        """Return a subsequence of single letter, use my_seq[index].

        >>> my_seq = Seq('ACTCGACGTCG')
        >>> my_seq[5]
        'A'
        """
        if isinstance(index, int):
            # Return a single letter as a string
            return chr(self._data[index])
        else:
            # Return the (sub)sequence as another Seq object
            return self.__class__(self._data[index])

    def __add__(self, other):
        """Add another sequence or string to this sequence.

        >>> from Bio.Seq import Seq
        >>> Seq("MELKI") + "LV"
        Seq('MELKILV')
        """
        if isinstance(other, (Seq, MutableSeq)):
            return self.__class__(self._data + other._data)
        elif isinstance(other, str):
            return self.__class__(self._data + other.encode("ASCII"))

        from Bio.SeqRecord import SeqRecord  # Lazy to avoid circular imports

        if isinstance(other, SeqRecord):
            # Get the SeqRecord's __radd__ to handle this
            return NotImplemented
        else:
            raise TypeError

    def __radd__(self, other):
        """Add a sequence on the left.

        >>> from Bio.Seq import Seq
        >>> "LV" + Seq("MELKI")
        Seq('LVMELKI')

        Adding two Seq (like) objects is handled via the __add__ method.
        """
        if isinstance(other, (Seq, MutableSeq)):
            return self.__class__(other._data + self._data)
        elif isinstance(other, str):
            return self.__class__(other.encode("ASCII") + self._data)
        else:
            raise TypeError

    def __mul__(self, other):
        """Multiply Seq by integer.

        >>> from Bio.Seq import Seq
        >>> Seq('ATG') * 2
        Seq('ATGATG')
        """
        if not isinstance(other, int):
            raise TypeError(f"can't multiply {self.__class__.__name__} by non-int type")
        return self.__class__(self._data * other)

    def __rmul__(self, other):
        """Multiply integer by Seq.

        >>> from Bio.Seq import Seq
        >>> 2 * Seq('ATG')
        Seq('ATGATG')
        """
        if not isinstance(other, int):
            raise TypeError(f"can't multiply {self.__class__.__name__} by non-int type")
        return self.__class__(self._data * other)

    def __imul__(self, other):
        """Multiply Seq in-place.

        >>> from Bio.Seq import Seq
        >>> seq = Seq('ATG')
        >>> seq *= 2
        >>> seq
        Seq('ATGATG')
        """
        if not isinstance(other, int):
            raise TypeError(f"can't multiply {self.__class__.__name__} by non-int type")
        return self.__class__(self._data * other)

    def tomutable(self):
        """Return the full sequence as a MutableSeq object.

        >>> from Bio.Seq import Seq
        >>> my_seq = Seq("MKQHKAMIVALIVICITAVVAAL")
        >>> my_seq
        Seq('MKQHKAMIVALIVICITAVVAAL')
        >>> my_seq.tomutable()
        MutableSeq('MKQHKAMIVALIVICITAVVAAL')
        """
        warnings.warn(
            "myseq.tomutable() is deprecated; please use MutableSeq(myseq) instead.",
            BiopythonDeprecationWarning,
        )
        return MutableSeq(self)

    def count(self, sub, start=None, end=None):
        """Return a non-overlapping count, like that of a python string.

        This behaves like the python string method of the same name,
        which does a non-overlapping count!

        For an overlapping search use the newer count_overlap() method.

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
        >>> my_seq = Seq("AAAATGA")
        >>> print(my_seq.count("A"))
        5
        >>> print(my_seq.count("ATG"))
        1
        >>> print(my_seq.count(Seq("AT")))
        1
        >>> print(my_seq.count("AT", 2, -1))
        1

        HOWEVER, please note because python strings and Seq objects (and
        MutableSeq objects) do a non-overlapping search, this may not give
        the answer you expect:

        >>> "AAAA".count("AA")
        2
        >>> print(Seq("AAAA").count("AA"))
        2

        An overlapping search, as implemented in .count_overlap(),
        would give the answer as three!
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

        For a non-overlapping search use the count() method.

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

        Where substrings do not overlap, should behave the same as
        the count() method:

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
        """Implement the 'in' keyword, like a python string.

        e.g.

        >>> from Bio.Seq import Seq
        >>> my_dna = Seq("ATATGAAATTTGAAAA")
        >>> "AAA" in my_dna
        True
        >>> Seq("AAA") in my_dna
        True
        """
        if isinstance(item, (Seq, MutableSeq)):
            item = bytes(item)
        elif isinstance(item, str):
            item = item.encode("ASCII")
        return item in self._data

    def find(self, sub, start=None, end=None):
        """Find method, like that of a python string.

        This behaves like the python string method of the same name.

        Returns an integer, the index of the first occurrence of substring
        argument sub in the (sub)sequence given by [start:end].

        Arguments:
         - sub - a string or another Seq object to look for
         - start - optional integer, slice start
         - end - optional integer, slice end

        Returns -1 if the subsequence is NOT found.

        e.g. Locating the first typical start codon, AUG, in an RNA sequence:

        >>> from Bio.Seq import Seq
        >>> my_rna = Seq("GUCAUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAGUUG")
        >>> my_rna.find("AUG")
        3
        """
        if isinstance(sub, (Seq, MutableSeq)):
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
        """Find from right method, like that of a python string.

        This behaves like the python string method of the same name.

        Returns an integer, the index of the last (right most) occurrence of
        substring argument sub in the (sub)sequence given by [start:end].

        Arguments:
         - sub - a string or another Seq object to look for
         - start - optional integer, slice start
         - end - optional integer, slice end

        Returns -1 if the subsequence is NOT found.

        e.g. Locating the last typical start codon, AUG, in an RNA sequence:

        >>> from Bio.Seq import Seq
        >>> my_rna = Seq("GUCAUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAGUUG")
        >>> my_rna.rfind("AUG")
        15
        """
        if isinstance(sub, (Seq, MutableSeq)):
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
        """Like find() but raise ValueError when the substring is not found.

        >>> from Bio.Seq import Seq
        >>> my_rna = Seq("GUCAUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAGUUG")
        >>> my_rna.find("T")
        -1
        >>> my_rna.index("T")
        Traceback (most recent call last):
                   ...
        ValueError: ...
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
        """Like rfind() but raise ValueError when the substring is not found."""
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
        """Return True if the Seq starts with the given prefix, False otherwise.

        This behaves like the python string method of the same name.

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
                bytes(p) if isinstance(p, (Seq, MutableSeq)) else p.encode("ASCII")
                for p in prefix
            )
        elif isinstance(prefix, (Seq, MutableSeq)):
            prefix = bytes(prefix)
        elif isinstance(prefix, str):
            prefix = prefix.encode("ASCII")
        return self._data.startswith(prefix, start, end)

    def endswith(self, suffix, start=None, end=None):
        """Return True if the Seq ends with the given suffix, False otherwise.

        This behaves like the python string method of the same name.

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
                bytes(p) if isinstance(p, (Seq, MutableSeq)) else p.encode("ASCII")
                for p in suffix
            )
        elif isinstance(suffix, (Seq, MutableSeq)):
            suffix = bytes(suffix)
        elif isinstance(suffix, str):
            suffix = suffix.encode("ASCII")
        return self._data.endswith(suffix, start, end)

    def split(self, sep=None, maxsplit=-1):
        """Split method, like that of a python string.

        This behaves like the python string method of the same name.

        Return a list of the 'words' in the string (as Seq objects),
        using sep as the delimiter string.  If maxsplit is given, at
        most maxsplit splits are done.  If maxsplit is omitted, all
        splits are made.

        Following the python string method, sep will by default be any
        white space (tabs, spaces, newlines) but this is unlikely to
        apply to biological sequences.

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

        See also the rsplit method:

        >>> for pep in my_aa.rsplit("*", 1):
        ...     pep
        Seq('VMAIVMGR*KGAR')
        Seq('L')
        """
        if isinstance(sep, (Seq, MutableSeq)):
            sep = bytes(sep)
        elif isinstance(sep, str):
            sep = sep.encode("ASCII")
        return [Seq(part) for part in self._data.split(sep, maxsplit)]

    def rsplit(self, sep=None, maxsplit=-1):
        """Do a right split method, like that of a python string.

        This behaves like the python string method of the same name.

        Return a list of the 'words' in the string (as Seq objects),
        using sep as the delimiter string.  If maxsplit is given, at
        most maxsplit splits are done COUNTING FROM THE RIGHT.
        If maxsplit is omitted, all splits are made.

        Following the python string method, sep will by default be any
        white space (tabs, spaces, newlines) but this is unlikely to
        apply to biological sequences.

        e.g. print(my_seq.rsplit("*",1))

        See also the split method.
        """
        if isinstance(sep, (Seq, MutableSeq)):
            sep = bytes(sep)
        elif isinstance(sep, str):
            sep = sep.encode("ASCII")
        return [Seq(part) for part in self._data.rsplit(sep, maxsplit)]

    def strip(self, chars=None):
        """Return a new Seq object with leading and trailing ends stripped.

        This behaves like the python string method of the same name.

        Optional argument chars defines which characters to remove.  If
        omitted or None (default) then as for the python string method,
        this defaults to removing any white space.

        e.g.

        >>> Seq("ACGT ").strip()
        Seq('ACGT')
        >>> Seq("ACGT ").strip(" ")
        Seq('ACGT')

        Just like the Python string, the order of the characters to be
        removed is not important:

        >>> Seq("ACGTACGT").strip("TGCA")
        Seq('')

        As with the Python string, an inappropriate argument
        will give a TypeError:

        >>> Seq("ACGT ").strip(7)
        Traceback (most recent call last):
           ...
        TypeError: argument must be None or a string, Seq, MutableSeq, or bytes-like object

        See also the lstrip and rstrip methods.
        """
        if isinstance(chars, (Seq, MutableSeq)):
            chars = bytes(chars)
        elif isinstance(chars, str):
            chars = chars.encode("ASCII")
        try:
            data = self._data.strip(chars)
        except TypeError:
            raise TypeError(
                "argument must be None or a string, Seq, MutableSeq, or bytes-like object"
            ) from None
        return Seq(data)

    def lstrip(self, chars=None):
        """Return a new Seq object with leading (left) end stripped.

        This behaves like the python string method of the same name.

        Optional argument chars defines which characters to remove.  If
        omitted or None (default) then as for the python string method,
        this defaults to removing any white space.

        >>> Seq("AAACGTA").lstrip("A")
        Seq('CGTA')

        See also the strip and rstrip methods.
        """
        if isinstance(chars, (Seq, MutableSeq)):
            chars = bytes(chars)
        elif isinstance(chars, str):
            chars = chars.encode("ASCII")
        try:
            data = self._data.lstrip(chars)
        except TypeError:
            raise TypeError(
                "argument must be None or a string, Seq, MutableSeq, or bytes-like object"
            ) from None
        return Seq(data)

    def rstrip(self, chars=None):
        """Return a new Seq object with trailing (right) end stripped.

        This behaves like the python string method of the same name.

        Optional argument chars defines which characters to remove.  If
        omitted or None (default) then as for the python string method,
        this defaults to removing any white space.

        e.g. Removing a nucleotide sequence's polyadenylation (poly-A tail):

        >>> from Bio.Seq import Seq
        >>> my_seq = Seq("CGGTACGCTTATGTCACGTAGAAAAAA")
        >>> my_seq
        Seq('CGGTACGCTTATGTCACGTAGAAAAAA')
        >>> my_seq.rstrip("A")
        Seq('CGGTACGCTTATGTCACGTAG')

        See also the strip and lstrip methods.
        """
        if isinstance(chars, (Seq, MutableSeq)):
            chars = bytes(chars)
        elif isinstance(chars, str):
            chars = chars.encode("ASCII")
        try:
            data = self._data.rstrip(chars)
        except TypeError:
            raise TypeError(
                "argument must be None or a string, Seq, MutableSeq, or bytes-like object"
            ) from None
        return Seq(data)

    def upper(self):
        """Return an upper case copy of the sequence.

        >>> from Bio.Seq import Seq
        >>> my_seq = Seq("VHLTPeeK*")
        >>> my_seq
        Seq('VHLTPeeK*')
        >>> my_seq.lower()
        Seq('vhltpeek*')
        >>> my_seq.upper()
        Seq('VHLTPEEK*')
        """
        return Seq(self._data.upper())

    def lower(self):
        """Return a lower case copy of the sequence.

        >>> from Bio.Seq import Seq
        >>> my_seq = Seq("CGGTACGCTTATGTCACGTAGAAAAAA")
        >>> my_seq
        Seq('CGGTACGCTTATGTCACGTAGAAAAAA')
        >>> my_seq.lower()
        Seq('cggtacgcttatgtcacgtagaaaaaa')

        See also the upper method.
        """
        return Seq(self._data.lower())

    def encode(self, encoding="utf-8", errors="strict"):
        """Return an encoded version of the sequence as a bytes object.

        The Seq object aims to match the interface of a Python string.

        This is essentially to save you doing str(my_seq).encode() when
        you need a bytes string, for example for computing a hash:

        >>> from Bio.Seq import Seq
        >>> Seq("ACGT").encode("ascii")
        b'ACGT'
        """
        warnings.warn(
            "myseq.encode has been deprecated; please use bytes(myseq) instead.",
            BiopythonDeprecationWarning,
        )
        return str(self).encode(encoding, errors)

    def complement(self):
        """Return the complement sequence by creating a new Seq object.

        This method is intended for use with DNA sequences:

        >>> from Bio.Seq import Seq
        >>> my_dna = Seq("CCCCCGATAG")
        >>> my_dna
        Seq('CCCCCGATAG')
        >>> my_dna.complement()
        Seq('GGGGGCTATC')

        You can of course used mixed case sequences,

        >>> from Bio.Seq import Seq
        >>> my_dna = Seq("CCCCCgatA-GD")
        >>> my_dna
        Seq('CCCCCgatA-GD')
        >>> my_dna.complement()
        Seq('GGGGGctaT-CH')

        Note in the above example, ambiguous character D denotes
        G, A or T so its complement is H (for C, T or A).

        Note that if the sequence contains neither T nor U, we
        assume it is DNA and map any A character to T:

        >>> Seq("CGA").complement()
        Seq('GCT')
        >>> Seq("CGAT").complement()
        Seq('GCTA')

        If you actually have RNA, this currently works but we
        may deprecate this later. We recommend using the new
        complement_rna method instead:

        >>> Seq("CGAU").complement()
        Seq('GCUA')
        >>> Seq("CGAU").complement_rna()
        Seq('GCUA')

        If the sequence contains both T and U, an exception is
        raised:

        >>> Seq("CGAUT").complement()
        Traceback (most recent call last):
           ...
        ValueError: Mixed RNA/DNA found

        Trying to complement a protein sequence gives a meaningless
        sequence:

        >>> my_protein = Seq("MAIVMGR")
        >>> my_protein.complement()
        Seq('KTIBKCY')

        Here "M" was interpreted as the IUPAC ambiguity code for
        "A" or "C", with complement "K" for "T" or "G". Likewise
        "A" has complement "T". The letter "I" has no defined
        meaning under the IUPAC convention, and is unchanged.
        """
        if isinstance(self._data, _UndefinedSequenceData):
            # complement of an undefined sequence is an undefined sequence
            # of the same length
            return self
        if (b"U" in self._data or b"u" in self._data) and (
            b"T" in self._data or b"t" in self._data
        ):
            # TODO - Handle this cleanly?
            raise ValueError("Mixed RNA/DNA found")
        elif b"U" in self._data or b"u" in self._data:
            ttable = _rna_complement_table
        else:
            ttable = _dna_complement_table
        # Much faster on really long sequences than the previous loop based
        # one. Thanks to Michael Palmer, University of Waterloo.
        return Seq(self._data.translate(ttable))

    def reverse_complement(self):
        """Return the reverse complement sequence by creating a new Seq object.

        This method is intended for use with DNA sequences:

        >>> from Bio.Seq import Seq
        >>> my_dna = Seq("CCCCCGATAGNR")
        >>> my_dna
        Seq('CCCCCGATAGNR')
        >>> my_dna.reverse_complement()
        Seq('YNCTATCGGGGG')

        Note in the above example, since R = G or A, its complement
        is Y (which denotes C or T).

        You can of course used mixed case sequences,

        >>> from Bio.Seq import Seq
        >>> my_dna = Seq("CCCCCgatA-G")
        >>> my_dna
        Seq('CCCCCgatA-G')
        >>> my_dna.reverse_complement()
        Seq('C-TatcGGGGG')

        As discussed for the complement method, if the sequence
        contains neither T nor U, is is assumed to be DNA and
        will map any letter A to T.

        If you are dealing with RNA you should use the new
        reverse_complement_rna method instead

        >>> Seq("CGA").reverse_complement()  # defaults to DNA
        Seq('TCG')
        >>> Seq("CGA").reverse_complement_rna()
        Seq('UCG')

        If the sequence contains both T and U, an exception is raised:

        >>> Seq("CGAUT").reverse_complement()
        Traceback (most recent call last):
           ...
        ValueError: Mixed RNA/DNA found

        Trying to reverse complement a protein sequence will give
        a meaningless sequence:

        >>> from Bio.Seq import Seq
        >>> my_protein = Seq("MAIVMGR")
        >>> my_protein.reverse_complement()
        Seq('YCKBITK')

        Here "M" was interpretted as the IUPAC ambiguity code for
        "A" or "C", with complement "K" for "T" or "G" - and so on.
        """
        # Use -1 stride/step to reverse the complement
        return self.complement()[::-1]

    def complement_rna(self):
        """Complement of an RNA sequence.

        >>> Seq("CGA").complement()  # defaults to DNA
        Seq('GCT')
        >>> Seq("CGA").complement_rna()
        Seq('GCU')

        Any T in the sequence is treated as a U:

        >>> Seq("CGAUT").complement_rna()
        Seq('GCUAA')
        """
        if isinstance(self._data, _UndefinedSequenceData):
            # complement of an undefined sequence is an undefined sequence
            # of the same length
            return self
        return Seq(self._data.translate(_rna_complement_table))

    def reverse_complement_rna(self):
        """Reverse complement of an RNA sequence.

        >>> from Bio.Seq import Seq
        >>> Seq("ACG").reverse_complement_rna()
        Seq('CGU')
        """
        # Use -1 stride/step to reverse the complement
        return self.complement_rna()[::-1]

    def transcribe(self):
        """Return the RNA sequence from a DNA sequence by creating a new Seq object.

        >>> from Bio.Seq import Seq
        >>> coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
        >>> coding_dna
        Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')
        >>> coding_dna.transcribe()
        Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')

        Trying to transcribe an RNA sequence should have no effect.
        If you have a nucleotide sequence which might be DNA or RNA
        (or even a mixture), calling the transcribe method will ensure
        any T becomes U.

        Trying to transcribe a protein sequence will replace any
        T for Threonine with U for Selenocysteine, which has no
        biologically plausible rational. Older versions of Biopython
        would throw an exception.

        >>> from Bio.Seq import Seq
        >>> my_protein = Seq("MAIVMGRT")
        >>> my_protein.transcribe()
        Seq('MAIVMGRU')
        """
        try:
            data = self._data.replace(b"T", b"U").replace(b"t", b"u")
        except UndefinedSequenceError:
            # transcribing an undefined sequence yields an undefined sequence
            # of the same length
            return self
        return Seq(data)

    def back_transcribe(self):
        """Return the DNA sequence from an RNA sequence by creating a new Seq object.

        >>> from Bio.Seq import Seq
        >>> messenger_rna = Seq("AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG")
        >>> messenger_rna
        Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')
        >>> messenger_rna.back_transcribe()
        Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')

        Trying to back-transcribe DNA has no effect, If you have a nucleotide
        sequence which might be DNA or RNA (or even a mixture), calling the
        back-transcribe method will ensure any T becomes U.

        Trying to back-transcribe a protein sequence will replace any U for
        Selenocysteine with T for Threonine, which is biologically meaningless.
        Older versions of Biopython would raise an exception here:

        >>> from Bio.Seq import Seq
        >>> my_protein = Seq("MAIVMGRU")
        >>> my_protein.back_transcribe()
        Seq('MAIVMGRT')
        """
        try:
            data = self._data.replace(b"U", b"T").replace(b"u", b"t")
        except UndefinedSequenceError:
            # back-transcribing an undefined sequence yields an undefined
            # sequence of the same length
            return self
        return Seq(data)

    def translate(
        self, table="Standard", stop_symbol="*", to_stop=False, cds=False, gap="-"
    ):
        """Turn a nucleotide sequence into a protein sequence by creating a new Seq object.

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
        if isinstance(table, str) and len(table) == 256:
            raise ValueError(
                "The Seq object translate method DOES NOT take "
                "a 256 character string mapping table like "
                "the python string object's translate method. "
                "Use str(my_seq).translate(...) instead."
            )
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
            _translate_str(data, table, stop_symbol, to_stop, cds, gap=gap)
        )

    def ungap(self, gap="-"):
        """Return a copy of the sequence without the gap character(s).

        The gap character now defaults to the minus sign, and can only
        be specified via the method argument. This is no longer possible
        via the sequence's alphabet (as was possible up to Biopython 1.77):

        >>> from Bio.Seq import Seq
        >>> my_dna = Seq("-ATA--TGAAAT-TTGAAAA")
        >>> my_dna
        Seq('-ATA--TGAAAT-TTGAAAA')
        >>> my_dna.ungap("-")
        Seq('ATATGAAATTTGAAAA')
        """
        if not gap:
            raise ValueError("Gap character required.")
        elif len(gap) != 1 or not isinstance(gap, str):
            raise ValueError(f"Unexpected gap character, {gap!r}")
        gap = gap.encode("ASCII")
        return Seq(self._data.replace(gap, b""))

    def join(self, other):
        """Return a merge of the sequences in other, spaced by the sequence from self.

        Accepts either a Seq or string (and iterates over the letters), or an
        iterable containing Seq or string objects. These arguments will be
        concatenated with the calling sequence as the spacer:

        >>> concatenated = Seq('NNNNN').join([Seq("AAA"), Seq("TTT"), Seq("PPP")])
        >>> concatenated
        Seq('AAANNNNNTTTNNNNNPPP')

        Joining the letters of a single sequence:

        >>> Seq('NNNNN').join(Seq("ACGT"))
        Seq('ANNNNNCNNNNNGNNNNNT')
        >>> Seq('NNNNN').join("ACGT")
        Seq('ANNNNNCNNNNNGNNNNNT')
        """
        if isinstance(other, (Seq, MutableSeq)):
            return self.__class__(str(self).join(str(other)))
        elif isinstance(other, str):
            return self.__class__(str(self).join(other))

        from Bio.SeqRecord import SeqRecord  # Lazy to avoid circular imports

        if isinstance(other, SeqRecord):
            raise TypeError("Iterable cannot be a SeqRecord")

        for c in other:
            if isinstance(c, SeqRecord):
                raise TypeError("Iterable cannot contain SeqRecords")
            elif not isinstance(c, (str, Seq, MutableSeq)):
                raise TypeError("Input must be an iterable of Seqs or Strings")
        return self.__class__(str(self).join([str(_) for _ in other]))


class UnknownSeq(Seq):
    """Read-only sequence object of known length but unknown contents (DEPRECATED).

    If you have an unknown sequence, you can represent this with a normal
    Seq object, for example:

    >>> my_seq = Seq("N"*5)
    >>> my_seq
    Seq('NNNNN')
    >>> len(my_seq)
    5
    >>> print(my_seq)
    NNNNN

    However, this is rather wasteful of memory (especially for large
    sequences), which is where this class is most useful:

    >>> unk_five = UnknownSeq(5)
    >>> unk_five
    UnknownSeq(5, character='?')
    >>> len(unk_five)
    5
    >>> print(unk_five)
    ?????

    You can add unknown sequence together. Provided the characters are the
    same, you get another memory saving UnknownSeq:

    >>> unk_four = UnknownSeq(4)
    >>> unk_four
    UnknownSeq(4, character='?')
    >>> unk_four + unk_five
    UnknownSeq(9, character='?')

    If the characters are different, addition gives an ordinary Seq object:

    >>> unk_nnnn = UnknownSeq(4, character="N")
    >>> unk_nnnn
    UnknownSeq(4, character='N')
    >>> unk_nnnn + unk_four
    Seq('NNNN????')

    Combining with a real Seq gives a new Seq object:

    >>> known_seq = Seq("ACGT")
    >>> unk_four + known_seq
    Seq('????ACGT')
    >>> known_seq + unk_four
    Seq('ACGT????')

    Although originally intended for unknown sequences (thus the class name),
    this can be used for homopolymer sequences like AAAAAA, and the biological
    methods will respect this:

    >>> homopolymer = UnknownSeq(6, character="A")
    >>> homopolymer.complement()
    UnknownSeq(6, character='T')
    >>> homopolymer.complement_rna()
    UnknownSeq(6, character='U')
    >>> homopolymer.translate()
    UnknownSeq(2, character='K')
    """

    def __init__(self, length, alphabet=None, character="?"):
        """Create a new UnknownSeq object.

        Arguments:
         - length - Integer, required.
         - alphabet - no longer used, must be None.
         - character - single letter string, default "?". Typically "N"
           for nucleotides, "X" for proteins, and "?" otherwise.
        """
        warnings.warn(
            "UnknownSeq(length) is deprecated; please use Seq(None, length) instead.",
            BiopythonDeprecationWarning,
        )
        if alphabet is not None:
            raise ValueError("The alphabet argument is no longer supported")
        self._length = int(length)
        if self._length < 0:
            # TODO - Block zero length UnknownSeq?  You can just use a Seq!
            raise ValueError("Length must not be negative.")
        if not character or len(character) != 1:
            raise ValueError("character argument should be a single letter string.")
        self._character = character

    def __len__(self):
        """Return the stated length of the unknown sequence."""
        return self._length

    def __bytes__(self):
        """Return the unknown sequence as full string of the given length."""
        return self._character.encode("ASCII") * self._length

    @property
    def _data(self):
        return self._character.encode("ASCII") * self._length

    def __str__(self):
        """Return the unknown sequence as full string of the given length."""
        return self._character * self._length

    def __repr__(self):
        """Return (truncated) representation of the sequence for debugging."""
        return f"UnknownSeq({self._length}, character={self._character!r})"

    def __add__(self, other):
        """Add another sequence or string to this sequence.

        Adding two UnknownSeq objects returns another UnknownSeq object
        provided the character is the same.

        >>> from Bio.Seq import UnknownSeq
        >>> UnknownSeq(10, character='X') + UnknownSeq(5, character='X')
        UnknownSeq(15, character='X')

        If the characters differ, an UnknownSeq object cannot be used, so a
        Seq object is returned:

        >>> from Bio.Seq import UnknownSeq
        >>> UnknownSeq(10, character='X') + UnknownSeq(5, character="x")
        Seq('XXXXXXXXXXxxxxx')

        If adding a string to an UnknownSeq, a new Seq is returned:

        >>> from Bio.Seq import UnknownSeq
        >>> UnknownSeq(5, character='X') + "LV"
        Seq('XXXXXLV')
        """
        if isinstance(other, UnknownSeq) and other._character == self._character:
            return UnknownSeq(len(self) + len(other), character=self._character)
        # Offload to the base class...
        return Seq(bytes(self)) + other

    def __radd__(self, other):
        """Add a sequence on the left."""
        # If other is an UnknownSeq, then __add__ would be called.
        # Offload to the base class...
        return other + Seq(bytes(self))

    def __mul__(self, other):
        """Multiply UnknownSeq by integer.

        >>> from Bio.Seq import UnknownSeq
        >>> UnknownSeq(3) * 2
        UnknownSeq(6, character='?')
        >>> UnknownSeq(3, character="N") * 2
        UnknownSeq(6, character='N')
        """
        if not isinstance(other, int):
            raise TypeError(f"can't multiply {self.__class__.__name__} by non-int type")
        return self.__class__(len(self) * other, character=self._character)

    def __rmul__(self, other):
        """Multiply integer by UnknownSeq.

        >>> from Bio.Seq import UnknownSeq
        >>> 2 * UnknownSeq(3)
        UnknownSeq(6, character='?')
        >>> 2 * UnknownSeq(3, character="N")
        UnknownSeq(6, character='N')
        """
        if not isinstance(other, int):
            raise TypeError(f"can't multiply {self.__class__.__name__} by non-int type")
        return self.__class__(len(self) * other, character=self._character)

    def __imul__(self, other):
        """Multiply UnknownSeq in-place.

        >>> from Bio.Seq import UnknownSeq
        >>> seq = UnknownSeq(3, character="N")
        >>> seq *= 2
        >>> seq
        UnknownSeq(6, character='N')
        """
        if not isinstance(other, int):
            raise TypeError(f"can't multiply {self.__class__.__name__} by non-int type")
        return self.__class__(len(self) * other, character=self._character)

    def __getitem__(self, index):
        """Get a subsequence from the UnknownSeq object.

        >>> unk = UnknownSeq(8, character="N")
        >>> print(unk[:])
        NNNNNNNN
        >>> print(unk[5:3])
        <BLANKLINE>
        >>> print(unk[1:-1])
        NNNNNN
        >>> print(unk[1:-1:2])
        NNN
        """
        if isinstance(index, int):
            if index >= -self._length and index < self._length:
                return self._character
            raise IndexError("sequence index out of range")
        start, stop, stride = index.indices(self._length)
        length = len(range(start, stop, stride))
        return UnknownSeq(length, character=self._character)

    def count(self, sub, start=None, end=None):
        """Return a non-overlapping count, like that of a python string.

        This behaves like the python string (and Seq object) method of the
        same name, which does a non-overlapping count!

        For an overlapping search use the newer count_overlap() method.

        Returns an integer, the number of occurrences of substring
        argument sub in the (sub)sequence given by [start:end].
        Optional arguments start and end are interpreted as in slice
        notation.

        Arguments:
         - sub - a string or another Seq object to look for
         - start - optional integer, slice start
         - end - optional integer, slice end

        >>> "NNNN".count("N")
        4
        >>> Seq("NNNN").count("N")
        4
        >>> UnknownSeq(4, character="N").count("N")
        4
        >>> UnknownSeq(4, character="N").count("A")
        0
        >>> UnknownSeq(4, character="N").count("AA")
        0

        HOWEVER, please note because that python strings and Seq objects (and
        MutableSeq objects) do a non-overlapping search, this may not give
        the answer you expect:

        >>> UnknownSeq(4, character="N").count("NN")
        2
        >>> UnknownSeq(4, character="N").count("NNN")
        1
        """
        if isinstance(sub, (Seq, MutableSeq)):
            sub = str(sub)
        elif not isinstance(sub, str):
            raise TypeError(
                "a Seq, MutableSeq, or string object is required, not '%s'" % type(sub)
            )
        # Handling case where subsequence not in self
        if set(sub) != set(self._character):
            return 0
        start, stop, stride = slice(start, end, len(sub)).indices(self._length)
        return len(range(start, stop - len(sub) + 1, stride))

    def count_overlap(self, sub, start=None, end=None):
        """Return an overlapping count.

        For a non-overlapping search use the count() method.

        Returns an integer, the number of occurrences of substring
        argument sub in the (sub)sequence given by [start:end].
        Optional arguments start and end are interpreted as in slice
        notation.

        Arguments:
         - sub - a string or another Seq object to look for
         - start - optional integer, slice start
         - end - optional integer, slice end

        e.g.

        >>> from Bio.Seq import UnknownSeq
        >>> UnknownSeq(4, character="N").count_overlap("NN")
        3
        >>> UnknownSeq(4, character="N").count_overlap("NNN")
        2

        Where substrings do not overlap, should behave the same as
        the count() method:

        >>> UnknownSeq(4, character="N").count_overlap("N")
        4
        >>> UnknownSeq(4, character="N").count_overlap("N") == UnknownSeq(4, character="N").count("N")
        True
        >>> UnknownSeq(4, character="N").count_overlap("A")
        0
        >>> UnknownSeq(4, character="N").count_overlap("A") == UnknownSeq(4, character="N").count("A")
        True
        >>> UnknownSeq(4, character="N").count_overlap("AA")
        0
        >>> UnknownSeq(4, character="N").count_overlap("AA") == UnknownSeq(4, character="N").count("AA")
        True
        """
        if isinstance(sub, (Seq, MutableSeq)):
            sub = str(sub)
        elif not isinstance(sub, str):
            raise TypeError(
                "a Seq, MutableSeq, or string object is required, not '%s'" % type(sub)
            )
        # Handling case where subsequence not in self
        if set(sub) != set(self._character):
            return 0
        start, stop, stride = slice(start, end).indices(self._length)
        return len(range(start, stop - len(sub) + 1, stride))

    def complement(self):
        """Return the complement assuming it is DNA.

        In typical usage this will return the same unknown sequence:

        >>> my_nuc = UnknownSeq(8, character='N')
        >>> my_nuc
        UnknownSeq(8, character='N')
        >>> print(my_nuc)
        NNNNNNNN
        >>> my_nuc.complement()
        UnknownSeq(8, character='N')
        >>> print(my_nuc.complement())
        NNNNNNNN

        If your sequence isn't actually unknown, and has a nucleotide letter
        other than N, the appropriate DNA complement base is used:

        >>> UnknownSeq(8, character="A").complement()
        UnknownSeq(8, character='T')
        """
        s = complement(self._character)
        return UnknownSeq(self._length, character=s)

    def complement_rna(self):
        """Return the complement assuming it is RNA.

        In typical usage this will return the same unknown sequence. If your
        sequence isn't actually unknown, the appropriate RNA complement base
        is used:

        >>> UnknownSeq(8, character="A").complement_rna()
        UnknownSeq(8, character='U')
        """
        s = complement_rna(self._character)
        return UnknownSeq(self._length, character=s)

    def reverse_complement(self):
        """Return the reverse complement assuming it is DNA.

        In typical usage this will return the same unknown sequence:

        >>> from Bio.Seq import UnknownSeq
        >>> example = UnknownSeq(6, character="N")
        >>> print(example)
        NNNNNN
        >>> print(example.reverse_complement())
        NNNNNN

        If your sequence isn't actually unknown, the appropriate DNA
        complement base is used:

        >>> UnknownSeq(8, character="A").reverse_complement()
        UnknownSeq(8, character='T')
        """
        return self.complement()

    def reverse_complement_rna(self):
        """Return the reverse complement assuming it is RNA.

        In typical usage this will return the same unknown sequence. If your
        sequence isn't actually unknown, the appropriate RNA complement base
        is used:

        >>> UnknownSeq(8, character="A").reverse_complement_rna()
        UnknownSeq(8, character='U')
        """
        return self.complement_rna()

    def transcribe(self):
        """Return an unknown RNA sequence from an unknown DNA sequence.

        >>> my_dna = UnknownSeq(10, character="N")
        >>> my_dna
        UnknownSeq(10, character='N')
        >>> print(my_dna)
        NNNNNNNNNN
        >>> my_rna = my_dna.transcribe()
        >>> my_rna
        UnknownSeq(10, character='N')
        >>> print(my_rna)
        NNNNNNNNNN

        In typical usage this will return the same unknown sequence. If your
        sequence isn't actually unknown, but a homopolymer of T, the standard
        DNA to RNA transcription is done, replacing T with U:

        >>> UnknownSeq(9, character="t").transcribe()
        UnknownSeq(9, character='u')
        """
        s = transcribe(self._character)
        return UnknownSeq(self._length, character=s)

    def back_transcribe(self):
        """Return an unknown DNA sequence from an unknown RNA sequence.

        >>> my_rna = UnknownSeq(20, character="N")
        >>> my_rna
        UnknownSeq(20, character='N')
        >>> print(my_rna)
        NNNNNNNNNNNNNNNNNNNN
        >>> my_dna = my_rna.back_transcribe()
        >>> my_dna
        UnknownSeq(20, character='N')
        >>> print(my_dna)
        NNNNNNNNNNNNNNNNNNNN

        In typical usage this will return the same unknown sequence. If your
        sequence is actually a U homopolymer, the standard RNA to DNA back
        translation applies, replacing U with T:

        >>> UnknownSeq(9, character="U").back_transcribe()
        UnknownSeq(9, character='T')
        """
        s = back_transcribe(self._character)
        return UnknownSeq(self._length, character=s)

    def upper(self):
        """Return an upper case copy of the sequence.

        >>> from Bio.Seq import UnknownSeq
        >>> my_seq = UnknownSeq(20, character="n")
        >>> my_seq
        UnknownSeq(20, character='n')
        >>> print(my_seq)
        nnnnnnnnnnnnnnnnnnnn
        >>> my_seq.upper()
        UnknownSeq(20, character='N')
        >>> print(my_seq.upper())
        NNNNNNNNNNNNNNNNNNNN

        See also the lower method.
        """
        return UnknownSeq(self._length, character=self._character.upper())

    def lower(self):
        """Return a lower case copy of the sequence.

        >>> from Bio.Seq import UnknownSeq
        >>> my_seq = UnknownSeq(20, character="X")
        >>> my_seq
        UnknownSeq(20, character='X')
        >>> print(my_seq)
        XXXXXXXXXXXXXXXXXXXX
        >>> my_seq.lower()
        UnknownSeq(20, character='x')
        >>> print(my_seq.lower())
        xxxxxxxxxxxxxxxxxxxx

        See also the upper method.
        """
        return UnknownSeq(self._length, character=self._character.lower())

    def translate(
        self, table="Standard", stop_symbol="*", to_stop=False, cds=False, gap="-"
    ):
        """Translate an unknown nucleotide sequence into an unknown protein.

        If your sequence makes sense as codons (e.g. a poly-A tail AAAAAA),
        it will be translated accordingly:

        >>> UnknownSeq(7, character='A').translate()
        UnknownSeq(2, character='K')

        Otherwise, it will be translated as X for unknown amino acid:

        >>> UnknownSeq(7).translate()
        UnknownSeq(2, character='X')
        """
        try:
            s = translate(
                self._character * 3,
                table=table,
                stop_symbol=stop_symbol,
                to_stop=to_stop,
                cds=cds,
                gap=gap,
            )
        except CodonTable.TranslationError:
            # Preserve historic behaviour, ??? (default character) and XXX -> X
            s = "X"
        # Don't worry about to_stop - no known stop codon is three bases the same,
        return UnknownSeq(self._length // 3, character=s)

    def ungap(self, gap="-"):
        """Return a copy of the sequence without the gap character(s).

        The gap character now defaults to the minus sign, and can only
        be specified via the method argument. This is no longer possible
        via the sequence's alphabet (as was possible up to Biopython 1.77):

        >>> from Bio.Seq import UnknownSeq
        >>> my_dna = UnknownSeq(20, character='N')
        >>> my_dna
        UnknownSeq(20, character='N')
        >>> my_dna.ungap()  # using default
        UnknownSeq(20, character='N')
        >>> my_dna.ungap("-")
        UnknownSeq(20, character='N')

        If the UnknownSeq is using the gap character, then an empty Seq is
        returned:

        >>> my_gap = UnknownSeq(20, character="-")
        >>> my_gap
        UnknownSeq(20, character='-')
        >>> my_gap.ungap()  # using default
        Seq('')
        >>> my_gap.ungap("-")
        Seq('')
        """
        # Offload the argument validation
        s = Seq(self._character).ungap(gap)
        if s:
            return UnknownSeq(self._length, character=self._character)
        else:
            return Seq("")

    def join(self, other):
        """Return a merge of the sequences in other, spaced by the sequence from self.

        Accepts either a Seq or string (and iterates over the letters), or an
        iterable containing Seq or string objects. These arguments will be
        concatenated with the calling sequence as the spacer:

        >>> concatenated = UnknownSeq(5).join([Seq("AAA"), Seq("TTT"), Seq("PPP")])
        >>> concatenated
        Seq('AAA?????TTT?????PPP')

        If all the inputs are also UnknownSeq using the same character, then it
        returns a new UnknownSeq:

        >>> UnknownSeq(5).join([UnknownSeq(3), UnknownSeq(3), UnknownSeq(3)])
        UnknownSeq(19, character='?')

        Examples taking a single sequence and joining the letters:

        >>> UnknownSeq(3).join("ACGT")
        Seq('A???C???G???T')
        >>> UnknownSeq(3).join(UnknownSeq(4))
        UnknownSeq(13, character='?')

        Will only return an UnknownSeq object if all of the objects to be joined are
        also UnknownSeqs with the same character as the spacer, similar to how the
        addition of an UnknownSeq and another UnknownSeq would work.
        """
        from Bio.SeqRecord import SeqRecord  # Lazy to avoid circular imports

        if isinstance(other, (str, Seq, MutableSeq)):
            if isinstance(other, UnknownSeq) and self._character == other._character:
                # Special case, can return an UnknownSeq
                return self.__class__(
                    len(other) + len(self) * (len(other) - 1), character=self._character
                )
            return Seq(str(self).join(str(other)))
        if isinstance(other, SeqRecord):
            raise TypeError("Iterable cannot be a SeqRecord")

        for c in other:
            if isinstance(c, SeqRecord):
                raise TypeError("Iterable cannot contain SeqRecords")
            elif not isinstance(c, (str, Seq, MutableSeq)):
                raise TypeError("Input must be an iterable of Seqs or Strings")
        temp_data = str(self).join([str(_) for _ in other])
        if temp_data.count(self._character) == len(temp_data):
            # Can return an UnknownSeq
            return self.__class__(len(temp_data), character=self._character)
        return Seq(temp_data)


class MutableSeq:
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
        if isinstance(data, array.array):
            if data.typecode != "u":
                raise ValueError(
                    "data should be a string, array of characters, Seq object, "
                    "or MutableSeq object"
                )
            warnings.warn(
                "Initializing a MutableSeq by an array has been deprecated; please "
                "use a bytearray object instead.",
                BiopythonDeprecationWarning,
            )
            data = data.tounicode()
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

    @property
    def data(self):
        """Get the data."""
        warnings.warn(
            "Accessing MutableSeq.data has been deprecated, as it is now a private "
            "attribute. Please use indexing to access the sequence contents of "
            "a MutableSeq object.",
            BiopythonDeprecationWarning,
        )
        return array.array("u", self._data.decode("ASCII"))

    @data.setter
    def data(self, value):
        """Set the data."""
        warnings.warn(
            "Accessing MutableSeq.data has been deprecated, as it is now a private "
            "attribute. Please use indexing to access the sequence contents of "
            "a MutableSeq object.",
            BiopythonDeprecationWarning,
        )
        self.__init__(value)

    def __bytes__(self):
        return bytes(self._data)

    def __repr__(self):
        """Return (truncated) representation of the sequence for debugging."""
        data = self._data
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
        """Return the full sequence as a python string.

        Note that Biopython 1.44 and earlier would give a truncated
        version of repr(my_seq) for str(my_seq).  If you are writing code
        which needs to be backwards compatible with old Biopython, you
        should continue to use my_seq.tostring() rather than str(my_seq).
        """
        return self._data.decode("ASCII")

    def __eq__(self, other):
        """Compare the sequence to another sequence or a string.

        Historically comparing DNA to RNA, or Nucleotide to Protein would
        raise an exception. This was later downgraded to a warning, but since
        Biopython 1.78 the alphabet is ignored for comparisons.

        If you need to support older Biopython versions, please just do
        explicit comparisons:

        >>> seq1 = MutableSeq("ACGT")
        >>> seq2 = MutableSeq("ACGT")
        >>> id(seq1) == id(seq2)
        False
        >>> str(seq1) == str(seq2)
        True

        Biopython now does:

        >>> seq1 == seq2
        True
        >>> seq1 == Seq("ACGT")
        True
        >>> seq1 == "ACGT"
        True

        """
        if isinstance(other, (Seq, MutableSeq)):
            return self._data == other._data
        elif isinstance(other, str):
            return self._data == other.encode("ASCII")
        else:
            return self._data == other

    def __lt__(self, other):
        """Implement the less-than operand."""
        if isinstance(other, (Seq, MutableSeq)):
            return self._data < other._data
        elif isinstance(other, str):
            return self._data < other.encode("ASCII")
        else:
            return self._data < other

    def __le__(self, other):
        """Implement the less-than or equal operand."""
        if isinstance(other, (Seq, MutableSeq)):
            return self._data <= other._data
        elif isinstance(other, str):
            return self._data <= other.encode("ASCII")
        else:
            return self._data <= other

    def __gt__(self, other):
        """Implement the greater-than operand."""
        if isinstance(other, (Seq, MutableSeq)):
            return self._data > other._data
        elif isinstance(other, str):
            return self._data > other.encode("ASCII")
        else:
            return self._data > other

    def __ge__(self, other):
        """Implement the greater-than or equal operand."""
        if isinstance(other, (Seq, MutableSeq)):
            return self._data >= other._data
        elif isinstance(other, str):
            return self._data >= other.encode("ASCII")
        else:
            return self._data >= other

    def __len__(self):
        """Return the length of the sequence, use len(my_seq)."""
        return len(self._data)

    def __getitem__(self, index):
        """Return a subsequence of single letter, use my_seq[index].

        >>> my_seq = MutableSeq('ACTCGACGTCG')
        >>> my_seq[5]
        'A'
        """
        if isinstance(index, int):
            # Return a single letter as a string
            return chr(self._data[index])
        else:
            # Return the (sub)sequence as another Seq object
            return self.__class__(self._data[index])

    def __setitem__(self, index, value):
        """Set a subsequence of single letter via value parameter.

        >>> my_seq = MutableSeq('ACTCGACGTCG')
        >>> my_seq[0] = 'T'
        >>> my_seq
        MutableSeq('TCTCGACGTCG')
        """
        if isinstance(index, int):
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
                raise TypeError("received unexpected type %s" % type(value))

    def __delitem__(self, index):
        """Delete a subsequence of single letter.

        >>> my_seq = MutableSeq('ACTCGACGTCG')
        >>> del my_seq[0]
        >>> my_seq
        MutableSeq('CTCGACGTCG')
        """
        # Could be deleting a single letter, or a slice
        del self._data[index]

    def __add__(self, other):
        """Add another sequence or string to this sequence.

        Returns a new MutableSeq object.
        """
        if isinstance(other, (Seq, MutableSeq)):
            return self.__class__(self._data + other._data)
        elif isinstance(other, str):
            return self.__class__(self._data + other.encode("ASCII"))

        from Bio.SeqRecord import SeqRecord  # Lazy to avoid circular imports

        if isinstance(other, SeqRecord):
            # Get the SeqRecord's __radd__ to handle this
            return NotImplemented
        else:
            raise TypeError

    def __radd__(self, other):
        """Add a sequence on the left.

        >>> from Bio.Seq import MutableSeq
        >>> "LV" + MutableSeq("MELKI")
        MutableSeq('LVMELKI')
        """
        if isinstance(other, (Seq, MutableSeq)):
            return self.__class__(other._data + self._data)
        elif isinstance(other, str):
            return self.__class__(other.encode("ASCII") + self._data)
        else:
            raise TypeError

    def __mul__(self, other):
        """Multiply MutableSeq by integer.

        Note this is not in-place and returns a new object,
        matching native Python list multiplication.

        >>> from Bio.Seq import MutableSeq
        >>> MutableSeq('ATG') * 2
        MutableSeq('ATGATG')
        """
        if not isinstance(other, int):
            raise TypeError(f"can't multiply {self.__class__.__name__} by non-int type")
        return self.__class__(self._data * other)

    def __rmul__(self, other):
        """Multiply integer by MutableSeq.

        Note this is not in-place and returns a new object,
        matching native Python list multiplication.

        >>> from Bio.Seq import MutableSeq
        >>> 2 * MutableSeq('ATG')
        MutableSeq('ATGATG')
        """
        if not isinstance(other, int):
            raise TypeError(f"can't multiply {self.__class__.__name__} by non-int type")
        return self.__class__(self._data * other)

    def __imul__(self, other):
        """Multiply MutableSeq in-place.

        >>> from Bio.Seq import MutableSeq
        >>> seq = MutableSeq('ATG')
        >>> seq *= 2
        >>> seq
        MutableSeq('ATGATG')
        """
        if not isinstance(other, int):
            raise TypeError(f"can't multiply {self.__class__.__name__} by non-int type")
        return self.__class__(self._data * other)

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

    def count(self, sub, start=None, end=None):
        """Return a non-overlapping count, like that of a python string.

        This behaves like the python string method of the same name,
        which does a non-overlapping count!

        For an overlapping search use the newer count_overlap() method.

        Returns an integer, the number of occurrences of substring
        argument sub in the (sub)sequence given by [start:end].
        Optional arguments start and end are interpreted as in slice
        notation.

        Arguments:
         - sub - a string or another Seq object to look for
         - start - optional integer, slice start
         - end - optional integer, slice end

        e.g.

        >>> from Bio.Seq import MutableSeq
        >>> my_mseq = MutableSeq("AAAATGA")
        >>> print(my_mseq.count("A"))
        5
        >>> print(my_mseq.count("ATG"))
        1
        >>> print(my_mseq.count(Seq("AT")))
        1
        >>> print(my_mseq.count("AT", 2, -1))
        1

        HOWEVER, please note because that python strings, Seq objects and
        MutableSeq objects do a non-overlapping search, this may not give
        the answer you expect:

        >>> "AAAA".count("AA")
        2
        >>> print(MutableSeq("AAAA").count("AA"))
        2

        An overlapping search would give the answer as three!
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

        For a non-overlapping search use the count() method.

        Returns an integer, the number of occurrences of substring
        argument sub in the (sub)sequence given by [start:end].
        Optional arguments start and end are interpreted as in slice
        notation.

        Arguments:
         - sub - a string or another Seq object to look for
         - start - optional integer, slice start
         - end - optional integer, slice end

        e.g.

        >>> from Bio.Seq import MutableSeq
        >>> print(MutableSeq("AAAA").count_overlap("AA"))
        3
        >>> print(MutableSeq("ATATATATA").count_overlap("ATA"))
        4
        >>> print(MutableSeq("ATATATATA").count_overlap("ATA", 3, -1))
        1

        Where substrings do not overlap, should behave the same as
        the count() method:

        >>> from Bio.Seq import MutableSeq
        >>> my_mseq = MutableSeq("AAAATGA")
        >>> print(my_mseq.count_overlap("A"))
        5
        >>> my_mseq.count_overlap("A") == my_mseq.count("A")
        True
        >>> print(my_mseq.count_overlap("ATG"))
        1
        >>> my_mseq.count_overlap("ATG") == my_mseq.count("ATG")
        True
        >>> print(my_mseq.count_overlap(Seq("AT")))
        1
        >>> my_mseq.count_overlap(Seq("AT")) == my_mseq.count(Seq("AT"))
        True
        >>> print(my_mseq.count_overlap("AT", 2, -1))
        1
        >>> my_mseq.count_overlap("AT", 2, -1) == my_mseq.count("AT", 2, -1)
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
        """Implement the 'in' keyword, like a python string.

        e.g.

        >>> from Bio.Seq import MutableSeq
        >>> my_dna = MutableSeq("ATATGAAATTTGAAAA")
        >>> "AAA" in my_dna
        True
        >>> MutableSeq("AAA") in my_dna
        True
        """
        if isinstance(item, (Seq, MutableSeq)):
            item = bytes(item)
        elif isinstance(item, str):
            item = item.encode("ASCII")
        return item in self._data

    def find(self, sub, start=None, end=None):
        """Find method, like that of a python string.

        This behaves like the python string method of the same name.

        Returns an integer, the index of the first occurrence of substring
        argument sub in the (sub)sequence given by [start:end].

        Arguments:
         - sub - a string or another Seq or MutableSeq object to look for
         - start - optional integer, slice start
         - end - optional integer, slice end

        Returns -1 if the subsequence is NOT found.

        e.g. Locating the first typical start codon, AUG, in an RNA sequence:

        >>> from Bio.Seq import MutableSeq
        >>> my_rna = MutableSeq("GUCAUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAGUUG")
        >>> my_rna.find("AUG")
        3
        """
        if isinstance(sub, (Seq, MutableSeq)):
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
        """Find from right method, like that of a python string.

        This behaves like the python string method of the same name.

        Returns an integer, the index of the last (right most) occurrence of
        substring argument sub in the (sub)sequence given by [start:end].

        Arguments:
         - sub - a string or another Seq or MutablSeq object to look for
         - start - optional integer, slice start
         - end - optional integer, slice end

        Returns -1 if the subsequence is NOT found.

        e.g. Locating the last typical start codon, AUG, in an RNA sequence:

        >>> from Bio.Seq import MutableSeq
        >>> my_rna = MutableSeq("GUCAUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAGUUG")
        >>> my_rna.rfind("AUG")
        15
        """
        if isinstance(sub, (Seq, MutableSeq)):
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
        """Return first occurrence position of a single entry (i.e. letter).

        >>> my_seq = MutableSeq("ACTCGACGTCG")
        >>> my_seq.index("A")
        0
        >>> my_seq.index("T")
        2
        >>> my_seq.index(Seq("T"))
        2
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
        """Like rfind() but raise ValueError when the substring is not found."""
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
        """Return True if the MutableSeq starts with the given prefix, False otherwise.

        This behaves like the python string method of the same name.

        Return True if the sequence starts with the specified prefix
        (a string or another Seq object), False otherwise.
        With optional start, test sequence beginning at that position.
        With optional end, stop comparing sequence at that position.
        prefix can also be a tuple of strings to try.  e.g.

        >>> from Bio.Seq import MutableSeq
        >>> my_rna = MutableSeq("GUCAUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAGUUG")
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
                bytes(p) if isinstance(p, (Seq, MutableSeq)) else p.encode("ASCII")
                for p in prefix
            )
        elif isinstance(prefix, (Seq, MutableSeq)):
            prefix = bytes(prefix)
        elif isinstance(prefix, str):
            prefix = prefix.encode("ASCII")
        return self._data.startswith(prefix, start, end)

    def endswith(self, suffix, start=None, end=None):
        """Return True if the MutableSeq ends with the given suffix, False otherwise.

        This behaves like the python string method of the same name.

        Return True if the sequence ends with the specified suffix
        (a string or another Seq object), False otherwise.
        With optional start, test sequence beginning at that position.
        With optional end, stop comparing sequence at that position.
        suffix can also be a tuple of strings to try.  e.g.

        >>> from Bio.Seq import MutableSeq
        >>> my_rna = MutableSeq("GUCAUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAGUUG")
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
                bytes(p) if isinstance(p, (Seq, MutableSeq)) else p.encode("ASCII")
                for p in suffix
            )
        elif isinstance(suffix, (Seq, MutableSeq)):
            suffix = bytes(suffix)
        elif isinstance(suffix, str):
            suffix = suffix.encode("ASCII")
        return self._data.endswith(suffix, start, end)

    def split(self, sep=None, maxsplit=-1):
        """Split method, like that of a python string.

        This behaves like the python string method of the same name.

        Return a list of the 'words' in the string (as MutableSeq objects),
        using sep as the delimiter string.  If maxsplit is given, at
        most maxsplit splits are done.  If maxsplit is omitted, all
        splits are made.

        Following the python string method, sep will by default be any
        white space (tabs, spaces, newlines) but this is unlikely to
        apply to biological sequences.

        e.g.

        >>> from Bio.Seq import MutableSeq
        >>> my_rna = MutableSeq("GUCAUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAGUUG")
        >>> my_aa = my_rna.translate()
        >>> my_aa
        MutableSeq('VMAIVMGR*KGAR*L')
        >>> for pep in my_aa.split("*"):
        ...     pep
        MutableSeq('VMAIVMGR')
        MutableSeq('KGAR')
        MutableSeq('L')
        >>> for pep in my_aa.split("*", 1):
        ...     pep
        MutableSeq('VMAIVMGR')
        MutableSeq('KGAR*L')

        See also the rsplit method:

        >>> for pep in my_aa.rsplit("*", 1):
        ...     pep
        MutableSeq('VMAIVMGR*KGAR')
        MutableSeq('L')
        """
        if isinstance(sep, (Seq, MutableSeq)):
            sep = bytes(sep)
        elif isinstance(sep, str):
            sep = sep.encode("ASCII")
        return [self.__class__(part) for part in self._data.split(sep, maxsplit)]

    def rsplit(self, sep=None, maxsplit=-1):
        """Do a right split method, like that of a python string.

        This behaves like the python string method of the same name.

        Return a list of the 'words' in the string (as MutableSeq objects),
        using sep as the delimiter string.  If maxsplit is given, at
        most maxsplit splits are done COUNTING FROM THE RIGHT.
        If maxsplit is omitted, all splits are made.

        Following the python string method, sep will by default be any
        white space (tabs, spaces, newlines) but this is unlikely to
        apply to biological sequences.

        e.g. print(my_seq.rsplit("*",1))

        See also the split method.
        """
        if isinstance(sep, (Seq, MutableSeq)):
            sep = bytes(sep)
        elif isinstance(sep, str):
            sep = sep.encode("ASCII")
        return [self.__class__(part) for part in self._data.rsplit(sep, maxsplit)]

    def strip(self, chars=None, inplace=False):
        """Return a MutableSeq object with leading and trailing ends stripped.

        This behaves like the python string method of the same name.

        A copy of the sequence is returned if ``inplace`` is `False` (the
        default value). If ``inplace`` is `True`, the sequence is stripped
        in-place and returned:

        >>> seq = MutableSeq("ACGT ")
        >>> seq.strip()
        MutableSeq('ACGT')
        >>> seq
        MutableSeq('ACGT ')
        >>> seq.strip(inplace=True)
        MutableSeq('ACGT')
        >>> seq
        MutableSeq('ACGT')

        Optional argument chars defines which characters to remove.  If
        omitted or None (default) then as for the python string method,
        this defaults to removing any white space.

        e.g.

        >>> MutableSeq("ACGT ").strip()
        MutableSeq('ACGT')
        >>> MutableSeq("ACGT ").strip(" ")
        MutableSeq('ACGT')

        Just like the Python string, the order of the characters to be
        removed is not important:

        >>> MutableSeq("ACGTACGT").strip("TGCA")
        MutableSeq('')

        As with the Python string, an inappropriate argument
        will give a TypeError:

        >>> MutableSeq("ACGT ").strip(7)
        Traceback (most recent call last):
           ...
        TypeError: argument must be None or a string, Seq, MutableSeq, or bytes-like object

        See also the lstrip and rstrip methods.
        """
        if isinstance(chars, (Seq, MutableSeq)):
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
            self._data[:] = data
            return self
        else:
            return MutableSeq(data)

    def lstrip(self, chars=None, inplace=False):
        """Return a MutableSeq object with leading (left) end stripped.

        This behaves like the python string method of the same name.

        A copy of the sequence is returned if ``inplace`` is `False` (the
        default value). If ``inplace`` is `True`, the sequence is stripped
        in-place and returned:

        >>> seq = MutableSeq(" ACGT ")
        >>> seq.lstrip()
        MutableSeq('ACGT ')
        >>> seq
        MutableSeq(' ACGT ')
        >>> seq.lstrip(inplace=True)
        MutableSeq('ACGT ')
        >>> seq
        MutableSeq('ACGT ')

        Optional argument chars defines which characters to remove.  If
        omitted or None (default) then as for the python string method,
        this defaults to removing any white space.

        >>> MutableSeq("AAACGTA").lstrip("A")
        MutableSeq('CGTA')

        See also the strip and rstrip methods.
        """
        if isinstance(chars, (Seq, MutableSeq)):
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
            self._data[:] = data
            return self
        else:
            return MutableSeq(data)

    def rstrip(self, chars=None, inplace=False):
        """Return a MutableSeq object with trailing (right) end stripped.

        This behaves like the python string method of the same name.

        A copy of the sequence is returned if ``inplace`` is `False` (the
        default value). If ``inplace`` is `True`, the sequence is stripped
        in-place and returned:

        >>> seq = MutableSeq(" ACGT ")
        >>> seq.rstrip()
        MutableSeq(' ACGT')
        >>> seq
        MutableSeq(' ACGT ')
        >>> seq.rstrip(inplace=True)
        MutableSeq(' ACGT')
        >>> seq
        MutableSeq(' ACGT')

        Optional argument chars defines which characters to remove.  If
        omitted or None (default) then as for the python string method,
        this defaults to removing any white space.

        e.g. Removing a nucleotide sequence's polyadenylation (poly-A tail):

        >>> from Bio.Seq import MutableSeq
        >>> my_seq = MutableSeq("CGGTACGCTTATGTCACGTAGAAAAAA")
        >>> my_seq
        MutableSeq('CGGTACGCTTATGTCACGTAGAAAAAA')
        >>> my_seq.rstrip("A")
        MutableSeq('CGGTACGCTTATGTCACGTAG')

        See also the strip and lstrip methods.
        """
        if isinstance(chars, (Seq, MutableSeq)):
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
            self._data[:] = data
            return self
        else:
            return MutableSeq(data)

    def upper(self, inplace=False):
        """Return the sequence in upper case.

        An upper-case copy of the sequence is returned if inplace is False,
        the default value:

        >>> from Bio.Seq import MutableSeq
        >>> my_seq = MutableSeq("VHLTPeeK*")
        >>> my_seq
        MutableSeq('VHLTPeeK*')
        >>> my_seq.lower()
        MutableSeq('vhltpeek*')
        >>> my_seq.upper()
        MutableSeq('VHLTPEEK*')
        >>> my_seq
        MutableSeq('VHLTPeeK*')

        The sequence is modified in-place and returned if inplace is True:

        >>> my_seq.lower(inplace=True)
        MutableSeq('vhltpeek*')
        >>> my_seq
        MutableSeq('vhltpeek*')
        >>> my_seq.upper(inplace=True)
        MutableSeq('VHLTPEEK*')
        >>> my_seq
        MutableSeq('VHLTPEEK*')
        """
        data = self._data.upper()
        if inplace:
            self._data[:] = data
            return self
        else:
            return MutableSeq(data)

    def lower(self, inplace=False):
        """Return the sequence in lower case.

        A lower-case copy of the sequence is returned if inplace is False,
        the default value:

        >>> from Bio.Seq import MutableSeq
        >>> my_seq = MutableSeq("VHLTPeeK*")
        >>> my_seq
        MutableSeq('VHLTPeeK*')
        >>> my_seq.lower()
        MutableSeq('vhltpeek*')
        >>> my_seq.upper()
        MutableSeq('VHLTPEEK*')
        >>> my_seq
        MutableSeq('VHLTPeeK*')

        The sequence is modified in-place and returned if inplace is True:

        >>> my_seq.lower(inplace=True)
        MutableSeq('vhltpeek*')
        >>> my_seq
        MutableSeq('vhltpeek*')
        >>> my_seq.upper(inplace=True)
        MutableSeq('VHLTPEEK*')
        >>> my_seq
        MutableSeq('VHLTPEEK*')
        """
        data = self._data.lower()
        if inplace:
            self._data[:] = data
            return self
        else:
            return MutableSeq(data)

    def translate(
        self, table="Standard", stop_symbol="*", to_stop=False, cds=False, gap="-"
    ):
        """Turn a nucleotide sequence into a protein sequence by creating a new MutableSeq object.

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

        e.g. Using the standard table:

        >>> coding_dna = MutableSeq("GTGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
        >>> coding_dna.translate()
        MutableSeq('VAIVMGR*KGAR*')
        >>> coding_dna.translate(stop_symbol="@")
        MutableSeq('VAIVMGR@KGAR@')
        >>> coding_dna.translate(to_stop=True)
        MutableSeq('VAIVMGR')

        Now using NCBI table 2, where TGA is not a stop codon:

        >>> coding_dna.translate(table=2)
        MutableSeq('VAIVMGRWKGAR*')
        >>> coding_dna.translate(table=2, to_stop=True)
        MutableSeq('VAIVMGRWKGAR')

        In fact, GTG is an alternative start codon under NCBI table 2, meaning
        this sequence could be a complete CDS:

        >>> coding_dna.translate(table=2, cds=True)
        MutableSeq('MAIVMGRWKGAR')

        It isn't a valid CDS under NCBI table 1, due to both the start codon
        and also the in frame stop codons:

        >>> coding_dna.translate(table=1, cds=True)
        Traceback (most recent call last):
            ...
        Bio.Data.CodonTable.TranslationError: First codon 'GTG' is not a start codon

        If the sequence has no in-frame stop codon, then the to_stop argument
        has no effect:

        >>> coding_dna2 = MutableSeq("TTGGCCATTGTAATGGGCCGC")
        >>> coding_dna2.translate()
        MutableSeq('LAIVMGR')
        >>> coding_dna2.translate(to_stop=True)
        MutableSeq('LAIVMGR')

        NOTE - Ambiguous codons like "TAN" or "NNN" could be an amino acid
        or a stop codon.  These are translated as "X".  Any invalid codon
        (e.g. "TA?" or "T-A") will throw a TranslationError.

        NOTE - This does NOT behave like the python string's translate
        method.  For that use str(my_seq).translate(...) instead
        """
        if isinstance(table, str) and len(table) == 256:
            raise ValueError(
                "The MutableSeq object translate method DOES NOT "
                "take a 256 character string mapping table like "
                "the python string object's translate method. "
                "Use str(my_seq).translate(...) instead."
            )

        return self.__class__(
            _translate_str(str(self), table, stop_symbol, to_stop, cds, gap=gap)
        )

    def reverse(self):
        """Modify the mutable sequence to reverse itself.

        No return value.
        """
        self._data.reverse()

    def complement(self):
        """Modify the mutable sequence to take on its complement.

        No return value.

        If the sequence contains neither T nor U, DNA is assumed
        and any A will be mapped to T.

        If the sequence contains both T and U, an exception is raised.
        """
        if ord("U") in self._data and ord("T") in self._data:
            raise ValueError("Mixed RNA/DNA found")
        elif ord("U") in self._data:
            table = _rna_complement_table
        else:
            table = _dna_complement_table
        self._data = self._data.translate(table)

    def reverse_complement(self):
        """Modify the mutable sequence to take on its reverse complement.

        No return value.
        """
        self.complement()
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

    def toseq(self):
        """Return the full sequence as a new immutable Seq object.

        >>> from Bio.Seq import MutableSeq
        >>> my_mseq = MutableSeq("MKQHKAMIVALIVICITAVVAAL")
        >>> my_mseq
        MutableSeq('MKQHKAMIVALIVICITAVVAAL')
        >>> my_mseq.toseq()
        Seq('MKQHKAMIVALIVICITAVVAAL')
        """
        warnings.warn(
            "myseq.toseq() is deprecated; please use Seq(myseq) instead.",
            BiopythonDeprecationWarning,
        )
        return Seq(self)

    def join(self, other):
        """Return a merge of the sequences in other, spaced by the sequence from self.

        Accepts all Seq objects and Strings as objects to be concatenated with the spacer

        >>> concatenated = MutableSeq('NNNNN').join([Seq("AAA"), Seq("TTT"), Seq("PPP")])
        >>> concatenated
        Seq('AAANNNNNTTTNNNNNPPP')

        Throws error if other is not an iterable and if objects inside of the iterable
        are not Seq or String objects
        """
        # returns Seq object instead of MutableSeq
        return Seq(self).join(other)


class UndefinedSequenceError(ValueError):
    """Sequence contents is undefined."""


class _UndefinedSequenceData(SequenceDataAbstractBaseClass):
    """Stores the length of a sequence with an undefined sequence contents (PRIVATE).

    Objects of this class can be used to create a Seq object to represent
    sequences with a known length, but an unknown sequence contents.
    Calling __len__ returns the sequence length, calling __getitem__ raises a
    ValueError except for requests of zero size, for which it returns an empty
    bytes object.
    """

    def __init__(self, length):
        """Initialize the object with the sequence length."""
        if length < 0:
            raise ValueError("Length must not be negative.")
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
        if self._length == 0:
            return b""
        raise UndefinedSequenceError("Sequence content is undefined")

    def __add__(self, other):
        if isinstance(other, _UndefinedSequenceData):
            return _UndefinedSequenceData(self._length + other._length)
        raise TypeError


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
    Bio.Data.CodonTable.TranslationError: Extra in frame stop codon found.
    """
    try:
        table_id = int(table)
    except ValueError:
        # Assume it's a table name
        # The same table can be used for RNA or DNA
        codon_table = CodonTable.ambiguous_generic_by_name[table]

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
                        "Extra in frame stop codon found."
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


def reverse_complement(sequence):
    """Return the reverse complement sequence of a nucleotide string.

    If given a string, returns a new string object.
    Given a Seq or a MutableSeq, returns a new Seq object.

    Supports unambiguous and ambiguous nucleotide sequences.

    e.g.

    >>> reverse_complement("ACTG-NH")
    'DN-CAGT'

    If neither T nor U is present, DNA is assumed and A is mapped to T:

    >>> reverse_complement("A")
    'T'
    """
    return complement(sequence)[::-1]


def complement(sequence):
    """Return the complement sequence of a DNA string.

    If given a string, returns a new string object.

    Given a Seq or a MutableSeq, returns a new Seq object.

    Supports unambiguous and ambiguous nucleotide sequences.

    e.g.

    >>> complement("ACTG-NH")
    'TGAC-ND'

    If neither T nor U is present, DNA is assumed and A is mapped to T:

    >>> complement("A")
    'T'

    However, this may not be supported in future. Please use the
    complement_rna function if you have RNA.
    """
    if isinstance(sequence, Seq):
        # Return a Seq
        return sequence.complement()
    elif isinstance(sequence, MutableSeq):
        # Return a Seq
        # Don't use the MutableSeq reverse_complement method as it is
        # 'in place'.
        return Seq(sequence).complement()

    # Assume it's a string.
    # In order to avoid some code duplication, the old code would turn the
    # string into a Seq, use the reverse_complement method, and convert back
    # to a string.
    # This worked, but is over five times slower on short sequences!
    sequence = sequence.encode("ASCII")
    if (b"U" in sequence or b"u" in sequence) and (
        b"T" in sequence or b"t" in sequence
    ):  # ugly but this is what black wants
        raise ValueError("Mixed RNA/DNA found")
    elif b"U" in sequence or b"u" in sequence:
        # TODO - warning or exception in future?
        ttable = _rna_complement_table
    else:
        ttable = _dna_complement_table
    sequence = sequence.translate(ttable)
    return sequence.decode("ASCII")


def complement_rna(sequence):
    """Return the complement sequence of an RNA string.

    >>> complement("ACG")  # assumed DNA
    'TGC'
    >>> complement_rna("ACG")
    'UGC'

    Any T in the sequence is treated as a U.
    """
    if isinstance(sequence, Seq):
        # Return a Seq
        return sequence.complement_rna()
    elif isinstance(sequence, MutableSeq):
        # Return a Seq
        return Seq(sequence).complement_rna()
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
