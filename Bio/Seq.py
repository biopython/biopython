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

import sys
import copy
from array import array
import warnings
from io import StringIO

from Bio import StreamModeError
from Bio import BiopythonWarning
from Bio import _seqobject
from Bio.Data.IUPACData import ambiguous_dna_complement, ambiguous_rna_complement
from Bio.Data.IUPACData import ambiguous_dna_letters as _ambiguous_dna_letters
from Bio.Data.IUPACData import ambiguous_rna_letters as _ambiguous_rna_letters
from Bio.Data import CodonTable


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
    before = "".join(complement_mapping.keys())
    after = "".join(complement_mapping.values())
    before += before.lower()
    after += after.lower()
    return str.maketrans(before, after)


_dna_complement_table = _maketrans(ambiguous_dna_complement)
_dna_complement_table[ord(b'U')] = _dna_complement_table[ord(b'T')]
_rna_complement_table = _maketrans(ambiguous_rna_complement)
_rna_complement_table[ord(b'T')] = _rna_complement_table[ord(b'U')]


class Seq(_seqobject.Seq):
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

    __slots__ = ()

    def __new__(cls, data=None, id="", name=None, description=None,
                annotations=None, features=None, dbxrefs=None,
                letter_annotations=None):
        """Create a Seq object.

        Arguments:
         - data - Sequence, required (string)

        You will typically use Bio.SeqIO to read in sequences from files as
        Seq objects.

        However, will often want to create your own Seq objects directly:

        >>> from Bio.Seq import Seq
        >>> my_seq = Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF")
        >>> my_seq
        Seq('MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF')
        >>> print(my_seq)
        MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF
        """
        if isinstance(data, int):
            raise TypeError
        if isinstance(data, str):
            data = bytes(data, "ASCII")
        elif data is not None:
            data = bytes(data)
        return super().__new__(cls, data)

    def __init__(self, data=None, id="", name=None, description=None,
                annotations=None, features=None, dbxrefs=None,
                letter_annotations=None):
        self.id = id
        if name is not None:
            self.name = name
        if description is not None:
            self.description = description
        if annotations is not None:
            self.annotations = annotations
        if features is not None:
            self.features = features
        if dbxrefs is not None:
            self.dbxrefs = dbxrefs
        if letter_annotations is not None:
            self.letter_annotations.update(letter_annotations)

    def __repr__(self):
        """Return (truncated) representation of the sequence."""
        if len(self) > 60:
            # Shows the last three letters as it is often useful to see if
            # there is a stop codon at the end of a sequence.
            # Note total length is 54+3+3=60
            term = f"'{self[:54]!s}...{self[-3:]!s}'"
        else:
            term = f"'{self!s}'"
        terms = [term]
        if self.id:
            terms.append(f"id={self.id!r}")
        if self.name:
            terms.append(f"name={self.name!r}")
        if self.description:
            terms.append(f"description={self.description!r}")
        return f"{self.__class__.__name__}({', '.join(terms)})"

    def __format__(self, format_spec):
        r"""Return the sequence as a string in the specified file format.

        This method supports the Python format() function and f-strings.
        The format_spec should be a lower case string supported by
        Bio.SeqIO as a text output file format.

        >>> from Bio.Seq import Seq
        >>> record = Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF",
        ...              id="YP_025292.1", name="HokC",
        ...              description="toxic membrane protein")
        ...
        >>> format(record, "fasta")
        '>YP_025292.1 toxic membrane protein\nMKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF\n'
        >>> print(f"Here is {record.id} in FASTA format:\n{record:fasta}")
        Here is YP_025292.1 in FASTA format:
        >YP_025292.1 toxic membrane protein
        MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF
        <BLANKLINE>
        """
        if not format_spec:
            # Follow python convention and default to using __str__
            return str(self)

        if format_spec == "debug":
            lines = []
            # Don't want to include the entire sequence
            if len(self) > 60:
                line = f"{str(self[:54])}...{str(self[-3:])}"
            else:
                line = str(self)
            lines.append(line)
            if self.id:
                lines.append(f"ID: {self.id}")
            if self.name:
                lines.append(f"Name: {self.name}")
            if self.description:
                lines.append(f"Description: {self.description}")
            if self.dbxrefs:
                lines.append("Database cross-references: " + ", ".join(self.dbxrefs))
            lines.append(f"Number of features: {len(self.features)}")
            for a in self.annotations:
                lines.append(f"/{a}={str(self.annotations[a])}")
            if self.letter_annotations:
                lines.append(
                    "Per letter annotation for: " + ", ".join(self.letter_annotations)
                )
            return "\n".join(lines)

        from Bio import SeqIO
        # Easy case, can call string-building function directly
        if format_spec in SeqIO._FormatToString:
            return SeqIO._FormatToString[format_spec](self)

        # Harder case, make a temp handle instead
        handle = StringIO()
        try:
            SeqIO.write(self, handle, format_spec)
        except StreamModeError:
            raise ValueError(
                "Binary format %s cannot be used with Seq format method"
                % format_spec
            ) from None
        return handle.getvalue()

    def __hash__(self):
        """Hash of the sequence as a string for comparison.

        See Seq object comparison documentation (method ``__eq__`` in
        particular).
        """
        return hash(str(self))

    def __getattr__(self, name):
        if name == "letter_annotations":
            self.letter_annotations = TrackDict(len(self))
            return self.letter_annotations
        elif name == "defined":
            return True
        raise AttributeError("'%s' object has no attribute '%s'" % (self.__class__.__name__, name))

    def __setattr__(self, name, value):
        if name == "letter_annotations":
            if not isinstance(value, TrackDict):
                raise ValueError("attribute 'letter_annotations' should be a TrackDict instance")
        super().__setattr__(name, value)

    def __getitem__(self, index):  # Seq API requirement
        """Return a subsequence of single letter, use my_seq[index].

        >>> my_seq = Seq('ACTCGACGTCG')
        >>> my_seq[5]
        'A'
        """
        letters = super().__getitem__(index)
        if isinstance(index, int):
            # Return a single letter as a string
            return letters.decode()
        else:
            # Return the (sub)sequence as another Seq object
            sequence = self.__class__(
                    letters,
                    id=self.id,
                    name=self.name,
                    description=self.description,
                )
            start, stop, step = index.indices(len(self))
            if step == 1:
                # Select relevant features, add them with shifted locations
                for f in self.features:
                    if f.ref or f.ref_db:
                        # TODO - Implement this (with lots of tests)?
                        import warnings

                        warnings.warn(
                            "When slicing SeqRecord objects, any "
                            "SeqFeature referencing other sequences (e.g. "
                            "from segmented GenBank records) are ignored."
                        )
                        continue
                    try:
                        feature_start = int(f.location.start)
                    except TypeError:  # e.g. UnknownPosition
                        continue
                    try:
                        feature_end = int(f.location.end)
                    except TypeError:  # e.g. UnknownPosition
                        continue
                    if (
                        start <= feature_start
                        and feature_end <= stop
                    ):
                        sequence.features.append(f._shift(-start))

            # Slice all the values to match the sliced sequence
            # (this should also work with strides, even negative strides):
            for key, value in self.letter_annotations.items():
                sequence.letter_annotations[key] = value[index]

            return sequence

    def __add__(self, other):
        """Add another sequence or string to this sequence.

        >>> from Bio.Seq import Seq
        >>> Seq("MELKI") + "LV"
        Seq('MELKILV')
        """
        if isinstance(other, str):
            sequence = self.__class__(bytes(self) + bytes(other, "ASCII"))
            sequence.id = self.id
            sequence.name = self.name
            sequence.description = self.description
            sequence.features = self.features[:]
            sequence.annotations = self.annotations.copy()
            sequence.dbxrefs = self.dbxrefs[:]
            return sequence

        if isinstance(other, (Seq, MutableSeq)):
            sequence = self.__class__(bytes(self) + bytes(other))
            sequence.features = self.features[:]
            sequence.dbxrefs = self.dbxrefs[:]
            # Will take all the features and all the db cross refs,
            length = len(self)
            for f in other.features:
                sequence.features.append(f._shift(length))
            for ref in other.dbxrefs:
                if ref not in sequence.dbxrefs:
                    sequence.dbxrefs.append(ref)
            # Take common id/name/description/annotation
            if self.id == other.id:
                sequence.id = self.id
            if self.name == other.name:
                sequence.name = self.name
            if self.description == other.description:
                sequence.description = self.description
            for k, v in self.annotations.items():
                if k in other.annotations and other.annotations[k] == v:
                    sequence.annotations[k] = v
            # Can append matching per-letter-annotation
            for k, v in self.letter_annotations.items():
                if k in other.letter_annotations:
                    sequence.letter_annotations[k] = v + other.letter_annotations[k]
            return sequence

        raise TypeError

    def __radd__(self, other):
        """Add a sequence on the left.

        >>> from Bio.Seq import Seq
        >>> "LV" + Seq("MELKI")
        Seq('LVMELKI')

        Adding two Seq (like) objects is handled via the __add__ method.
        """
        # Assume it is a string or a Seq.
        # Note can't transfer any per-letter-annotations
        offset = len(other)
        return self.__class__(
            other + str(self),
            id=self.id,
            name=self.name,
            description=self.description,
            features=[f._shift(offset) for f in self.features],
            annotations=self.annotations.copy(),
            dbxrefs=self.dbxrefs[:],
        )

    def __mul__(self, other):
        """Multiply Seq by integer.

        >>> from Bio.Seq import Seq
        >>> Seq('ATG') * 2
        Seq('ATGATG')
        """
        if not isinstance(other, int):
            raise TypeError(f"can't multiply {self.__class__.__name__} by non-int type")
        return self.__class__(bytes(self) * other)

    def __rmul__(self, other):
        """Multiply integer by Seq.

        >>> from Bio.Seq import Seq
        >>> 2 * Seq('ATG')
        Seq('ATGATG')
        """
        if not isinstance(other, int):
            raise TypeError(f"can't multiply {self.__class__.__name__} by non-int type")
        return self.__class__(bytes(self) * other)

    def count_overlap(self, sub, start=0, end=sys.maxsize):
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
        sub_str = str(sub)
        self_str = str(self)
        overlap_count = 0
        while True:
            start = self_str.find(sub_str, start, end) + 1
            if start != 0:
                overlap_count += 1
            else:
                return overlap_count

    def __contains__(self, char):
        """Implement the 'in' keyword, like a python string.

        e.g.

        >>> from Bio.Seq import Seq
        >>> my_dna = Seq("ATATGAAATTTGAAAA")
        >>> "AAA" in my_dna
        True
        >>> Seq("AAA") in my_dna
        True
        """
        if isinstance(char, str):
            char = char.encode("ASCII")
        return char in bytes(self)

    def find(self, sub, start=0, end=sys.maxsize):
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
        if isinstance(sub, str):
            sub = sub.encode("ASCII")
        return bytes(self).find(sub, start, end)

    def rfind(self, sub, start=0, end=sys.maxsize):
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
        if isinstance(sub, str):
            sub = sub.encode("ASCII")
        return bytes(self).rfind(sub, start, end)

    def index(self, sub, start=0, end=sys.maxsize):
        """Like find() but raise ValueError when the substring is not found.

        >>> from Bio.Seq import Seq
        >>> my_rna = Seq("GUCAUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAGUUG")
        >>> my_rna.find("T")
        -1
        >>> my_rna.index("T")
        Traceback (most recent call last):
                   ...
        ValueError: substring not found...
        """
        if isinstance(sub, str):
            sub = sub.encode("ASCII")
        try:
            return bytes(self).index(sub, start, end)
        except ValueError:
            raise ValueError("substring not found")

    def rindex(self, sub, start=0, end=sys.maxsize):
        """Like rfind() but raise ValueError when the substring is not found."""
        if isinstance(sub, str):
            sub = sub.encode("ASCII")
        try:
            return bytes(self).rindex(sub, start, end)
        except ValueError:
            raise ValueError("Seq.rindex(item): item not in list")

    def startswith(self, prefix, start=0, end=sys.maxsize):
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
            prefix_strs = tuple(str(p) for p in prefix)
            return str(self).startswith(prefix_strs, start, end)
        else:
            return str(self).startswith(str(prefix), start, end)

    def endswith(self, suffix, start=0, end=sys.maxsize):
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
            suffix_strs = tuple(str(p) for p in suffix)
            return str(self).endswith(suffix_strs, start, end)
        else:
            return str(self).endswith(str(suffix), start, end)

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
        return [Seq(part) for part in str(self).split(str(sep), maxsplit)]

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
        return [Seq(part) for part in str(self).rsplit(str(sep), maxsplit)]

    def strip(self, chars=None):
        """Return a new Seq object with leading and trailing ends stripped.

        This behaves like the python string method of the same name.

        Optional argument chars defines which characters to remove.  If
        omitted or None (default) then as for the python string method,
        this defaults to removing any white space.

        e.g. print(my_seq.strip("-"))

        See also the lstrip and rstrip methods.
        """
        return Seq(str(self).strip(str(chars)))

    def lstrip(self, chars=None):
        """Return a new Seq object with leading (left) end stripped.

        This behaves like the python string method of the same name.

        Optional argument chars defines which characters to remove.  If
        omitted or None (default) then as for the python string method,
        this defaults to removing any white space.

        e.g. print(my_seq.lstrip("-"))

        See also the strip and rstrip methods.
        """
        return Seq(str(self).lstrip(str(chars)))

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
        return Seq(str(self).rstrip(str(chars)))

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
        return Seq(str(self).upper())

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
        return Seq(str(self).lower())

    def encode(self, encoding="utf-8", errors="strict"):
        """Return an encoded version of the sequence as a bytes object.

        The Seq object aims to match the interface of a Python string.

        This is essentially to save you doing str(my_seq).encode() when
        you need a bytes string, for example for computing a hash:

        >>> from Bio.Seq import Seq
        >>> Seq("ACGT").encode("ascii")
        b'ACGT'
        """
        return str(self).encode(encoding, errors)

    def complement(self, id=None, name=None, description=None, dbxrefs=None, annotations=None, features=None, letter_annotations=None):
        """Return the complement sequence by creating a new Seq object.

        This method will return a DNA nucleotide sequence; use
        ``rna_complement`` to generate an RNA nucleotide sequence instead.

        >>> from Bio.Seq import Seq
        >>> my_dna = Seq("CCCCCGATAG")
        >>> my_dna
        Seq('CCCCCGATAG')
        >>> my_dna.complement()
        Seq('GGGGGCTATC')
        >>> my_dna.rna_complement()
        Seq('GGGGGCUAUC')

        >>> from Bio.Seq import Seq
        >>> my_rna = Seq("CCCCCGAUAG")
        >>> my_rna
        Seq('CCCCCGAUAG')
        >>> my_rna.complement()
        Seq('GGGGGCTATC')
        >>> my_rna.rna_complement()
        Seq('GGGGGCUAUC')

        You can of course use mixed case sequences,

        >>> from Bio.Seq import Seq
        >>> my_dna = Seq("CCCCCgatANGD")
        >>> my_dna
        Seq('CCCCCgatANGD')
        >>> my_dna.complement()
        Seq('GGGGGctaTNCH')

        Note in the above example, ambiguous character D denotes
        G, A or T so its complement is H (for C, T or A).

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
        sequence = Seq(super().complement())
        if id is True:
            sequence.id = self.id
        elif id is not None and id is not False:
            sequence.id = id
        if name is True:
            sequence.name = self.name
        elif name is not None and name is not False:
            sequence.name = name
        if description is True:
            sequence.description = self.description
        elif description is not None and description is not False:
            sequence.description = description
        if annotations is True:
            sequence.annotations = copy.deepcopy(self.annotations)
        elif annotations is not None and annotations is not False:
            sequence.annotations = annotations
        if dbxrefs is True:
            sequence.dbxrefs = copy.deepcopy(self.dbxrefs)
        elif dbxrefs is not None and dbxrefs is not False:
            sequence.dbxrefs = dbxrefs
        if features is True:
            sequence.features = copy.deepcopy(self.features)
        elif features is not None and features is not False:
            sequence.features = features
        if letter_annotations is True:
            sequence.letter_annotations = copy.deepcopy(self.letter_annotations)
        elif letter_annotations is not None and letter_annotations is not False:
            sequence.letter_annotations = letter_annotations
        return sequence

    def reverse_complement(self, id=None, name=None, description=None, dbxrefs=None, annotations=None, features=True, letter_annotations=None):
        """Return the reverse complement sequence by creating a new Seq object.

        This method will return a DNA nucleotide sequence; use
        ``rna_complement`` to generate an RNA nucleotide sequence instead.

        >>> from Bio.Seq import Seq
        >>> my_dna = Seq("CCCCCGATAG")
        >>> my_dna
        Seq('CCCCCGATAG')
        >>> my_dna.reverse_complement()
        Seq('CTATCGGGGG')
        >>> my_dna.rna_reverse_complement()
        Seq('CUAUCGGGGG')

        >>> from Bio.Seq import Seq
        >>> my_rna = Seq("CCCCCGAUAG")
        >>> my_rna
        Seq('CCCCCGAUAG')
        >>> my_rna.reverse_complement()
        Seq('CTATCGGGGG')
        >>> my_rna.rna_reverse_complement()
        Seq('CUAUCGGGGG')

        You can use mixed case and ambiguous sequences:

        >>> from Bio.Seq import Seq
        >>> my_dna = Seq("CCCCCGATAGNR")
        >>> my_dna
        Seq('CCCCCGATAGNR')
        >>> my_dna.reverse_complement()
        Seq('YNCTATCGGGGG')

        Note in the above example, since R = G or A, its complement
        is Y (which denotes C or T).

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
        ttable = _dna_complement_table
        # Much faster on really long sequences than the previous loop based
        # one. Thanks to Michael Palmer, University of Waterloo.
        sequence = Seq(str(self).translate(ttable))[::-1]
        if id is True:
            sequence.id = self.id
        elif id is not None and id is not False:
            sequence.id = id
        if name is True:
            sequence.name = self.name
        elif name is not None and name is not False:
            sequence.name = name
        if description is True:
            sequence.description = self.description
        elif description is not None and description is not False:
            sequence.description = description
        if annotations is True:
            sequence.annotations = copy.deepcopy(self.annotations)
        elif annotations is not None and annotations is not False:
            sequence.annotations = annotations
        if dbxrefs is True:
            sequence.dbxrefs = copy.deepcopy(self.dbxrefs)
        elif dbxrefs is not None and dbxrefs is not False:
            sequence.dbxrefs = dbxrefs
        if features is True:
            length = len(sequence)
            sequence.features = [feature._flip(length) for feature in self.features]
        elif features is not None and features is not False:
            sequence.features = features
        if letter_annotations is True:
            for key, value in self.letter_annotations.items():
                sequence.letter_annotations[key] = copy.copy(value[::-1])
        elif letter_annotations is not None and letter_annotations is not False:
            sequence.letter_annotations.update(letter_annotations)
        return sequence


    def rna_complement(self):
        """RNA complement of a nucleotide sequence.

        >>> Seq("CGA").complement()  # returns a DNA sequence
        Seq('GCT')
        >>> Seq("CGA").rna_complement()  # returns an RNA sequence
        Seq('GCU')

        """
        return Seq(str(self).translate(_rna_complement_table))

    def rna_reverse_complement(self):
        """RNA reverse complement of a nucleotide sequence.

        >>> from Bio.Seq import Seq
        >>> Seq("ACG").rna_reverse_complement()
        Seq('CGU')
        """
        # Use -1 stride/step to reverse the complement
        return self.rna_complement()[::-1]

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
        return Seq(str(self).replace("T", "U").replace("t", "u"))

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
        return Seq(str(self).replace("U", "T").replace("u", "t"))

    def translate(
        self, table="Standard", stop_symbol="*", to_stop=False, cds=False, gap="-", id=None, name=None, description=None, dbxrefs=None, annotations=None, features=None, letter_annotations=None
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
        if "protein" == self.annotations.get("molecule_type", ""):
            raise ValueError("Proteins cannot be translated!")
        if isinstance(table, str) and len(table) == 256:
            raise ValueError(
                "The Seq object translate method DOES NOT take "
                "a 256 character string mapping table like "
                "the python string object's translate method. "
                "Use str(my_seq).translate(...) instead."
            )
        try:
            table_id = int(table)
        except ValueError:
            # Assume its a table name
            # The same table can be used for RNA or DNA
            codon_table = CodonTable.ambiguous_generic_by_name[table]

        except (AttributeError, TypeError):
            # Assume its a CodonTable object
            if isinstance(table, CodonTable.CodonTable):
                codon_table = table
            else:
                raise ValueError("Bad table argument") from None
        else:
            # Assume its a table ID
            # The same table can be used for RNA or DNA
            codon_table = CodonTable.ambiguous_generic_by_id[table_id]

        sequence = Seq(
            _translate_str(str(self), codon_table, stop_symbol, to_stop, cds, gap=gap)
        )
        if id is True:
            sequence.id = self.id
        elif id is not None and id is not False:
            sequence.id = id
        if name is True:
            sequence.name = self.name
        elif name is not None and name is not False:
            sequence.name = name
        if description is True:
            sequence.description = self.description
        elif description is not None and description is not False:
            sequence.description = description
        if annotations is True:
            sequence.annotations = copy.deepcopy(self.annotations)
        elif annotations is not None and annotations is not False:
            sequence.annotations = annotations
        if dbxrefs is True:
            sequence.dbxrefs = copy.deepcopy(self.dbxrefs)
        elif dbxrefs is not None and dbxrefs is not False:
            sequence.dbxrefs = dbxrefs
        if features is True:
            raise TypeError("features cannot be preserved during translation")
        elif features is not None and features is not False:
            sequence.features = features
        if letter_annotations is not None:
            sequence.letter_annotations.update(letter_annotations)
        sequence.annotations["molecule_type"] = "protein"

        return sequence

    def ungap(self, gap="-"):
        """Return a copy of the sequence without the gap character(s).

        The gap character now defaults to the minus sign, and can only
        be specified via the method argument. This is no longer possible
        via sequence's alphabet (as was possible up to Biopython 1.77):

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
        return Seq(str(self).replace(gap, ""))

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
        if isinstance(other, (str, Seq, MutableSeq)):
            return self.__class__(str(self).join(str(other)))

        for c in other:
            if not isinstance(c, (str, Seq, MutableSeq)):
                raise TypeError("Input must be an iterable of Seqs or Strings")
        return self.__class__(str(self).join([str(_) for _ in other]))


class UnknownSeq(Seq):
    """Read-only sequence object of known length but unknown contents.

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
    """

    def __new__(cls, length, character="?",
                     id=None, name=None, description=None, annotations=None,
                     features=None, dbxrefs=None, letter_annotations=None):
        """Create a new UnknownSeq object.

        Arguments:
         - length - Integer, required.
         - character - single letter string, default "?". Typically "N"
           for nucleotides, "X" for proteins, and "?" otherwise.
        """
        self = super(Seq, cls).__new__(cls, data=length, character=character)
        self._character = character
        return self

    def __init__(self, length, character="?",
                     id=None, name=None, description=None, annotations=None,
                     features=None, dbxrefs=None, letter_annotations=None):
        super().__init__(None, id=id, name=name, description=description,
                         annotations=annotations, features=features,
                         dbxrefs=dbxrefs, letter_annotations=letter_annotations)

    def __repr__(self):
        """Return (truncated) representation of the sequence for debugging."""
        return f"UnknownSeq({len(self)}, character={self._character!r})"

    def __getattr__(self, name):
        if name == "defined":
            return False
        return super().__getattr__(name)

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
        return Seq(str(self)) + other

    def __radd__(self, other):
        """Add a sequence on the left."""
        # If other is an UnknownSeq, then __add__ would be called.
        # Offload to the base class...
        return other + Seq(str(self))

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
            # TODO - Check the bounds without wasting memory
            return str(self)[index]
        old_length = len(self)
        step = index.step
        if step is None or step == 1:
            # This calculates the length you'd get from ("N"*old_length)[index]
            start = index.start
            end = index.stop
            if start is None:
                start = 0
            elif start < 0:
                start = max(0, old_length + start)
            elif start > old_length:
                start = old_length
            if end is None:
                end = old_length
            elif end < 0:
                end = max(0, old_length + end)
            elif end > old_length:
                end = old_length
            new_length = max(0, end - start)
        elif step == 0:
            raise ValueError("slice step cannot be zero")
        else:
            # TODO - handle step efficiently
            new_length = len(("X" * old_length)[index])
        # assert new_length == len(("X"*old_length)[index]), \
        #       (index, start, end, step, old_length,
        #        new_length, len(("X"*old_length)[index]))
        sequence = UnknownSeq(new_length, character=self._character,
                              id=self.id, name=self.name,
                              description=self.description)
        start, stop, step = index.indices(len(self))
        if step == 1:
            # Select relevant features, add them with shifted locations
            for f in self.features:
                if f.ref or f.ref_db:
                    # TODO - Implement this (with lots of tests)?
                    import warnings

                    warnings.warn(
                        "When slicing SeqRecord objects, any "
                        "SeqFeature referencing other sequences (e.g. "
                        "from segmented GenBank records) are ignored."
                    )
                    continue
                try:
                    feature_start = int(f.location.start)
                except TypeError:  # e.g. UnknownPosition
                    continue
                try:
                    feature_end = int(f.location.end)
                except TypeError:  # e.g. UnknownPosition
                    continue
                if (
                    start <= feature_start
                    and feature_end <= stop
                ):
                    sequence.features.append(f._shift(-start))

        # Slice all the values to match the sliced sequence
        # (this should also work with strides, even negative strides):
        for key, value in self.letter_annotations.items():
            sequence.letter_annotations[key] = value[index]

        return sequence

    def complement(self):
        """Return the complement of an unknown nucleotide equals itself.

        >>> my_nuc = UnknownSeq(8)
        >>> my_nuc
        UnknownSeq(8, character='?')
        >>> print(my_nuc)
        ????????
        >>> my_nuc.complement()
        UnknownSeq(8, character='?')
        >>> print(my_nuc.complement())
        ????????
        """
        return self

    def rna_complement(self):
        """Return the complement assuming it is RNA."""
        return self

    def reverse_complement(self):
        """Return the reverse complement of an unknown sequence.

        The reverse complement of an unknown nucleotide equals itself:

        >>> from Bio.Seq import UnknownSeq
        >>> example = UnknownSeq(6, character="N")
        >>> print(example)
        NNNNNN
        >>> print(example.reverse_complement())
        NNNNNN
        """
        return self

    def rna_reverse_complement(self):
        """Return the reverse complement assuming it is RNA."""
        return self

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
        """
        s = Seq(self._character).transcribe()
        return UnknownSeq(len(self), character=str(s))

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
        """
        s = Seq(self._character).back_transcribe()
        return UnknownSeq(len(self), character=str(s))

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
        return UnknownSeq(len(self), character=self._character.upper())

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
        return UnknownSeq(len(self), character=self._character.lower())

    def translate(
        self, table="Standard", stop_symbol="*", to_stop=False, cds=False, gap="-"
    ):
        """Translate an unknown nucleotide sequence into an unknown protein.

        e.g.

        >>> my_seq = UnknownSeq(9, character="N")
        >>> print(my_seq)
        NNNNNNNNN
        >>> my_protein = my_seq.translate()
        >>> my_protein
        UnknownSeq(3, character='X')
        >>> print(my_protein)
        XXX

        In comparison, using a normal Seq object:

        >>> my_seq = Seq("NNNNNNNNN")
        >>> print(my_seq)
        NNNNNNNNN
        >>> my_protein = my_seq.translate()
        >>> my_protein
        Seq('XXX')
        >>> print(my_protein)
        XXX

        """
        return UnknownSeq(len(self) // 3, character="X")

    def ungap(self, gap="-"):
        """Return a copy of the sequence without the gap character(s).

        The gap character now defaults to the minus sign, and can only
        be specified via the method argument. This is no longer possible
        via sequence's alphabet (as was possible up to Biopython 1.77):

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
            return UnknownSeq(len(self), character=self._character)
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
        if isinstance(other, (str, Seq, MutableSeq)):
            if isinstance(other, UnknownSeq) and self._character == other._character:
                # Special case, can return an UnknownSeq
                return self.__class__(
                    len(other) + len(self) * (len(other) - 1), character=self._character
                )
            return Seq(str(self).join(str(other)))

        for c in other:
            if not isinstance(c, (str, Seq, MutableSeq)):
                raise TypeError("Input must be an iterable of Seqs or Strings")
        temp_data = str(self).join([str(_) for _ in other])
        if temp_data.count(self._character) == len(temp_data):
            # Can return an UnknownSeq
            return self.__class__(len(temp_data), character=self._character)
        return Seq(temp_data)


class MutableSeq(Seq):
    """An editable sequence object.

    Unlike normal python strings and our basic sequence object (the Seq class),
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

    def __new__(cls, data, id=None, name=None, description=None,
                     annotations=None, features=None, dbxrefs=None,
                     letter_annotations=None):
        """Initialize the class."""
        if isinstance(data, int):
            raise TypeError
        if isinstance(data, str):
            data = bytearray(data, "ASCII")
        else:
            data = bytearray(data)
        self = super(Seq, cls).__new__(cls, data)
        if id is not None:
            self.id = id
        if name is not None:
            self.name = name
        if description is not None:
            self.description = description
        if annotations is not None:
            self.annotations = annotations
        if features is not None:
            self.features = features
        if letter_annotations is not None:
            self.letter_annotations.update(letter_annotations)
        return self

    def __repr__(self):
        """Return (truncated) representation of the sequence for debugging."""
        if len(self) > 60:
            # Shows the last three letters as it is often useful to see if
            # there is a stop codon at the end of a sequence.
            # Note total length is 54+3+3=60
            return f"{self.__class__.__name__}('{str(self[:54])}...{str(self[-3:])}')"
        else:
            return f"{self.__class__.__name__}('{str(self)}')"

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
        try:
            index = self.index(item)
        except ValueError:
            raise ValueError("MutableSeq.remove(x): x not in list")
        del self[index]

    def complement(self):
        """Modify the mutable sequence to take on its DNA complement.

        No return value.
        """
        super(Seq, self).complement()

    def rna_complement(self):
        """Modify the mutable sequence into its RNA complement.

        No return value.
        """
        super(Seq, self).rna_complement()

    def reverse_complement(self, inplace=True):
        """Modify the mutable sequence to take on its reverse complement.

        No return value.
        """
        if inplace:
            self.complement()
            self.reverse()
        else:
            sequence = MutableSeq(self)
            sequence.complement()
            sequence.reverse()
            return sequence

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


# The transcribe, backward_transcribe, and translate functions are
# user-friendly versions of the corresponding functions in Bio.Transcribe
# and Bio.Translate. The functions work both on Seq objects, and on strings.


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
     - table - a CodonTable object (NOT a table name or id number)
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
    sequence = sequence.upper()
    amino_acids = []
    forward_table = table.forward_table
    stop_codons = table.stop_codons
    if table.nucleotide_alphabet is not None:
        valid_letters = set(table.nucleotide_alphabet.upper())
    else:
        # Assume the worst case, ambiguous DNA or RNA:
        valid_letters = set(
            _ambiguous_dna_letters.upper() + _ambiguous_rna_letters.upper()
        )
    n = len(sequence)

    # Check for tables with 'ambiguous' (dual-coding) stop codons:
    dual_coding = [c for c in stop_codons if c in forward_table]
    if dual_coding:
        c = dual_coding[0]
        if to_stop:
            raise ValueError(
                "You cannot use 'to_stop=True' with this table as it contains"
                f" {len(dual_coding)} codon(s) which can be both  STOP and an"
                f" amino acid (e.g. '{c}' -> '{forward_table[c]}' or STOP)."
            )
        warnings.warn(
            f"This table contains {len(dual_coding)} codon(s) which code(s) for"
            f" both STOP and an amino acid (e.g. '{c}' -> '{forward_table[c]}'"
            " or STOP). Such codons will be translated as amino acid.",
            BiopythonWarning,
        )

    if cds:
        if str(sequence[:3]).upper() not in table.start_codons:
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
            if codon in table.stop_codons:
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
        # Assume its a string, return a string
        try:
            codon_table = CodonTable.ambiguous_generic_by_id[int(table)]
        except ValueError:
            codon_table = CodonTable.ambiguous_generic_by_name[table]
        except (AttributeError, TypeError):
            if isinstance(table, CodonTable.CodonTable):
                codon_table = table
            else:
                raise ValueError("Bad table argument") from None
        return _translate_str(sequence, codon_table, stop_symbol, to_stop, cds, gap=gap)


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
    if isinstance(sequence, MutableSeq):
        # Return a Seq
        # Don't use the MutableSeq reverse_complement method as it is
        # 'in place'.
        return Seq(sequence).complement()

    elif isinstance(sequence, Seq):
        # Return a Seq
        return sequence.complement()
    # Assume its a string.
    # In order to avoid some code duplication, the old code would turn the
    # string into a Seq, use the reverse_complement method, and convert back
    # to a string.
    # This worked, but is over five times slower on short sequences!
    ttable = _dna_complement_table
    return sequence.translate(ttable)


def rna_complement(sequence):
    """Return the complement sequence of an RNA string.

    >>> complement("ACG")
    'TGC'
    >>> rna_complement("ACG")
    'UGC'
    """
    if isinstance(sequence, Seq):
        # Return a Seq
        return sequence.rna_complement()
    elif isinstance(sequence, MutableSeq):
        # Return a Seq
        return Seq(sequence).rna_complement()
    return sequence.translate(_rna_complement_table)


def rna_reverse_complement(sequence):
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
    return rna_complement(sequence)[::-1]

class TrackDict(dict):
    """Dict which only allows objects of given length as values (PRIVATE).

    This simple subclass of the Python dictionary is used in Seq objects
    for holding per-letter-annotations.  Only python objects with a fixed
    length (e.g. strings, tuples, bytes, and Python arrays) matching the
    length of the sequence stored in the Seq object.

    >>> x = TrackDict(5)
    >>> x["test"] = "hello"
    >>> x
    {'test': 'hello'}

    Attempting to store an entry which doesn't have the expected length
    raises a ValueError:

    >>> x["test"] = "hello world"
    Traceback (most recent call last):
    ...
    ValueError: value has length 11 (expected 5)

    Attempting to store an entry with a mutable length raise a TypeError:

    >>> x["test"] = [8, 6, 4, 1, 9]
    Traceback (most recent call last):
    ...
    TypeError: only python objects with a fixed length (strings, tuples, bytes, or Python arrays) can be stored

    The expected length is stored as an attribute:

    >>> x.length
    5

    In order that the SeqRecord (and other objects using this class) can be
    pickled, for example for use in the multiprocessing library, we need to
    be able to pickle the restricted dictionary objects.

    Using the default protocol, which is 3 on Python 3,

    >>> import pickle
    >>> y = pickle.loads(pickle.dumps(x))
    >>> y
    {'test': 'hello'}
    >>> y.length
    5

    Using the highest protocol, which is 4 on Python 3,

    >>> import pickle
    >>> z = pickle.loads(pickle.dumps(x, pickle.HIGHEST_PROTOCOL))
    >>> z
    {'test': 'hello'}
    >>> z.length
    5
    """

    def __new__(cls, length, items=None):
        """Create an empty TrackDict for the specified sequence length."""
        self = super().__new__(cls)
        self.length = int(length)
        return self

    def __init__(self, length, items=None):
        if items is not None:
            self.update(items)

    def __reduce__(self):
        arguments = (self.length, dict(self))
        return (TrackDict, arguments)

    def __setitem__(self, key, value):
        if not isinstance(value, (str, tuple, bytes, array)):
            raise TypeError(
                "only python objects with a fixed length (strings, tuples, "
                "bytes, or Python arrays) can be stored"
            )
        if len(value) != self.length:
            raise ValueError("value has length %d (expected %d)" % (len(value), self.length))

        dict.__setitem__(self, key, value)

    def update(self, new_dict):
        # Force this to go via our strict __setitem__ method
        if not isinstance(new_dict, dict):
            raise TypeError("expected a dictionary object; got %s" % new_dict.__class__.__name__)
        for (key, value) in new_dict.items():
            self[key] = value


def _test():
    """Run the Bio.Seq module's doctests (PRIVATE)."""
    print("Running doctests...")
    import doctest

    doctest.testmod(optionflags=doctest.IGNORE_EXCEPTION_DETAIL)
    print("Done")


if __name__ == "__main__":
    _test()
