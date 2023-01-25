# Copyright 2000, 2004 by Brad Chapman.
# Revisions copyright 2010-2013, 2015-2018 by Peter Cock.
# All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Code for dealing with sequence alignments.

One of the most important things in this module is the MultipleSeqAlignment
class, used in the Bio.AlignIO module.

"""

import sys
import collections
import copy
import importlib
import warnings
import numbers
from itertools import zip_longest

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Please install numpy if you want to use Bio.Align. "
        "See http://www.numpy.org/"
    ) from None

from Bio import BiopythonDeprecationWarning
from Bio.Align import _aligners
from Bio.Align import substitution_matrices
from Bio.Seq import Seq, MutableSeq, reverse_complement, UndefinedSequenceError
from Bio.SeqRecord import SeqRecord, _RestrictedDict

# Import errors may occur here if a compiled aligners.c file
# (_aligners.pyd or _aligners.so) is missing or if the user is
# importing from within the Biopython source tree, see PR #2007:
# https://github.com/biopython/biopython/pull/2007


AlignmentCounts = collections.namedtuple(
    "AlignmentCounts", ["gaps", "identities", "mismatches"]
)


class MultipleSeqAlignment:
    """Represents a classical multiple sequence alignment (MSA).

    By this we mean a collection of sequences (usually shown as rows) which
    are all the same length (usually with gap characters for insertions or
    padding). The data can then be regarded as a matrix of letters, with well
    defined columns.

    You would typically create an MSA by loading an alignment file with the
    AlignIO module:

    >>> from Bio import AlignIO
    >>> align = AlignIO.read("Clustalw/opuntia.aln", "clustal")
    >>> print(align)
    Alignment with 7 rows and 156 columns
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273285|gb|AF191659.1|AF191
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273284|gb|AF191658.1|AF191
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273287|gb|AF191661.1|AF191
    TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273286|gb|AF191660.1|AF191
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273290|gb|AF191664.1|AF191
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273289|gb|AF191663.1|AF191
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273291|gb|AF191665.1|AF191

    In some respects you can treat these objects as lists of SeqRecord objects,
    each representing a row of the alignment. Iterating over an alignment gives
    the SeqRecord object for each row:

    >>> len(align)
    7
    >>> for record in align:
    ...     print("%s %i" % (record.id, len(record)))
    ...
    gi|6273285|gb|AF191659.1|AF191 156
    gi|6273284|gb|AF191658.1|AF191 156
    gi|6273287|gb|AF191661.1|AF191 156
    gi|6273286|gb|AF191660.1|AF191 156
    gi|6273290|gb|AF191664.1|AF191 156
    gi|6273289|gb|AF191663.1|AF191 156
    gi|6273291|gb|AF191665.1|AF191 156

    You can also access individual rows as SeqRecord objects via their index:

    >>> print(align[0].id)
    gi|6273285|gb|AF191659.1|AF191
    >>> print(align[-1].id)
    gi|6273291|gb|AF191665.1|AF191

    And extract columns as strings:

    >>> print(align[:, 1])
    AAAAAAA

    Or, take just the first ten columns as a sub-alignment:

    >>> print(align[:, :10])
    Alignment with 7 rows and 10 columns
    TATACATTAA gi|6273285|gb|AF191659.1|AF191
    TATACATTAA gi|6273284|gb|AF191658.1|AF191
    TATACATTAA gi|6273287|gb|AF191661.1|AF191
    TATACATAAA gi|6273286|gb|AF191660.1|AF191
    TATACATTAA gi|6273290|gb|AF191664.1|AF191
    TATACATTAA gi|6273289|gb|AF191663.1|AF191
    TATACATTAA gi|6273291|gb|AF191665.1|AF191

    Combining this alignment slicing with alignment addition allows you to
    remove a section of the alignment. For example, taking just the first
    and last ten columns:

    >>> print(align[:, :10] + align[:, -10:])
    Alignment with 7 rows and 20 columns
    TATACATTAAGTGTACCAGA gi|6273285|gb|AF191659.1|AF191
    TATACATTAAGTGTACCAGA gi|6273284|gb|AF191658.1|AF191
    TATACATTAAGTGTACCAGA gi|6273287|gb|AF191661.1|AF191
    TATACATAAAGTGTACCAGA gi|6273286|gb|AF191660.1|AF191
    TATACATTAAGTGTACCAGA gi|6273290|gb|AF191664.1|AF191
    TATACATTAAGTATACCAGA gi|6273289|gb|AF191663.1|AF191
    TATACATTAAGTGTACCAGA gi|6273291|gb|AF191665.1|AF191

    Note - This object does NOT attempt to model the kind of alignments used
    in next generation sequencing with multiple sequencing reads which are
    much shorter than the alignment, and where there is usually a consensus or
    reference sequence with special status.
    """

    def __init__(
        self, records, alphabet=None, annotations=None, column_annotations=None
    ):
        """Initialize a new MultipleSeqAlignment object.

        Arguments:
         - records - A list (or iterator) of SeqRecord objects, whose
                     sequences are all the same length.  This may be an be an
                     empty list.
         - alphabet - For backward compatibility only; its value should always
                      be None.
         - annotations - Information about the whole alignment (dictionary).
         - column_annotations - Per column annotation (restricted dictionary).
                      This holds Python sequences (lists, strings, tuples)
                      whose length matches the number of columns. A typical
                      use would be a secondary structure consensus string.

        You would normally load a MSA from a file using Bio.AlignIO, but you
        can do this from a list of SeqRecord objects too:

        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> from Bio.Align import MultipleSeqAlignment
        >>> a = SeqRecord(Seq("AAAACGT"), id="Alpha")
        >>> b = SeqRecord(Seq("AAA-CGT"), id="Beta")
        >>> c = SeqRecord(Seq("AAAAGGT"), id="Gamma")
        >>> align = MultipleSeqAlignment([a, b, c],
        ...                              annotations={"tool": "demo"},
        ...                              column_annotations={"stats": "CCCXCCC"})
        >>> print(align)
        Alignment with 3 rows and 7 columns
        AAAACGT Alpha
        AAA-CGT Beta
        AAAAGGT Gamma
        >>> align.annotations
        {'tool': 'demo'}
        >>> align.column_annotations
        {'stats': 'CCCXCCC'}
        """
        if alphabet is not None:
            raise ValueError("The alphabet argument is no longer supported")

        self._records = []
        if records:
            self.extend(records)

        # Annotations about the whole alignment
        if annotations is None:
            annotations = {}
        elif not isinstance(annotations, dict):
            raise TypeError("annotations argument should be a dict")
        self.annotations = annotations

        # Annotations about each column of the alignment
        if column_annotations is None:
            column_annotations = {}
        # Handle this via the property set function which will validate it
        self.column_annotations = column_annotations

    def _set_per_column_annotations(self, value):
        if not isinstance(value, dict):
            raise TypeError(
                "The per-column-annotations should be a (restricted) dictionary."
            )
        # Turn this into a restricted-dictionary (and check the entries)
        if len(self):
            # Use the standard method to get the length
            expected_length = self.get_alignment_length()
            self._per_col_annotations = _RestrictedDict(length=expected_length)
            self._per_col_annotations.update(value)
        else:
            # Bit of a problem case... number of columns is undefined
            self._per_col_annotations = None
            if value:
                raise ValueError(
                    "Can't set per-column-annotations without an alignment"
                )

    def _get_per_column_annotations(self):
        if self._per_col_annotations is None:
            # This happens if empty at initialisation
            if len(self):
                # Use the standard method to get the length
                expected_length = self.get_alignment_length()
            else:
                # Should this raise an exception? Compare SeqRecord behaviour...
                expected_length = 0
            self._per_col_annotations = _RestrictedDict(length=expected_length)
        return self._per_col_annotations

    column_annotations = property(
        fget=_get_per_column_annotations,
        fset=_set_per_column_annotations,
        doc="""Dictionary of per-letter-annotation for the sequence.""",
    )

    def _str_line(self, record, length=50):
        """Return a truncated string representation of a SeqRecord (PRIVATE).

        This is a PRIVATE function used by the __str__ method.
        """
        if record.seq.__class__.__name__ == "CodonSeq":
            if len(record.seq) <= length:
                return f"{record.seq} {record.id}"
            else:
                return "%s...%s %s" % (
                    record.seq[: length - 3],
                    record.seq[-3:],
                    record.id,
                )
        else:
            if len(record.seq) <= length:
                return f"{record.seq} {record.id}"
            else:
                return "%s...%s %s" % (
                    record.seq[: length - 6],
                    record.seq[-3:],
                    record.id,
                )

    def __str__(self):
        """Return a multi-line string summary of the alignment.

        This output is intended to be readable, but large alignments are
        shown truncated.  A maximum of 20 rows (sequences) and 50 columns
        are shown, with the record identifiers.  This should fit nicely on a
        single screen. e.g.

        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> from Bio.Align import MultipleSeqAlignment
        >>> a = SeqRecord(Seq("ACTGCTAGCTAG"), id="Alpha")
        >>> b = SeqRecord(Seq("ACT-CTAGCTAG"), id="Beta")
        >>> c = SeqRecord(Seq("ACTGCTAGATAG"), id="Gamma")
        >>> align = MultipleSeqAlignment([a, b, c])
        >>> print(align)
        Alignment with 3 rows and 12 columns
        ACTGCTAGCTAG Alpha
        ACT-CTAGCTAG Beta
        ACTGCTAGATAG Gamma

        See also the alignment's format method.
        """
        rows = len(self._records)
        lines = [
            "Alignment with %i rows and %i columns"
            % (rows, self.get_alignment_length())
        ]
        if rows <= 20:
            lines.extend(self._str_line(rec) for rec in self._records)
        else:
            lines.extend(self._str_line(rec) for rec in self._records[:18])
            lines.append("...")
            lines.append(self._str_line(self._records[-1]))
        return "\n".join(lines)

    def __repr__(self):
        """Return a representation of the object for debugging.

        The representation cannot be used with eval() to recreate the object,
        which is usually possible with simple python objects.  For example:

        <Bio.Align.MultipleSeqAlignment instance (2 records of length 14)
        at a3c184c>

        The hex string is the memory address of the object, see help(id).
        This provides a simple way to visually distinguish alignments of
        the same size.
        """
        # A doctest for __repr__ would be nice, but __class__ comes out differently
        # if run via the __main__ trick.
        return "<%s instance (%i records of length %i) at %x>" % (
            self.__class__,
            len(self._records),
            self.get_alignment_length(),
            id(self),
        )
        # This version is useful for doing eval(repr(alignment)),
        # but it can be VERY long:
        # return "%s(%r)" \
        #       % (self.__class__, self._records)

    def __format__(self, format_spec):
        """Return the alignment as a string in the specified file format.

        The format should be a lower case string supported as an output
        format by Bio.AlignIO (such as "fasta", "clustal", "phylip",
        "stockholm", etc), which is used to turn the alignment into a
        string.

        e.g.

        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> from Bio.Align import MultipleSeqAlignment
        >>> a = SeqRecord(Seq("ACTGCTAGCTAG"), id="Alpha", description="")
        >>> b = SeqRecord(Seq("ACT-CTAGCTAG"), id="Beta", description="")
        >>> c = SeqRecord(Seq("ACTGCTAGATAG"), id="Gamma", description="")
        >>> align = MultipleSeqAlignment([a, b, c])
        >>> print(format(align, "fasta"))
        >Alpha
        ACTGCTAGCTAG
        >Beta
        ACT-CTAGCTAG
        >Gamma
        ACTGCTAGATAG
        <BLANKLINE>
        >>> print(format(align, "phylip"))
         3 12
        Alpha      ACTGCTAGCT AG
        Beta       ACT-CTAGCT AG
        Gamma      ACTGCTAGAT AG
        <BLANKLINE>
        """
        if format_spec:
            from io import StringIO
            from Bio import AlignIO

            handle = StringIO()
            AlignIO.write([self], handle, format_spec)
            return handle.getvalue()
        else:
            # Follow python convention and default to using __str__
            return str(self)

    def __iter__(self):
        """Iterate over alignment rows as SeqRecord objects.

        e.g.

        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> from Bio.Align import MultipleSeqAlignment
        >>> a = SeqRecord(Seq("ACTGCTAGCTAG"), id="Alpha")
        >>> b = SeqRecord(Seq("ACT-CTAGCTAG"), id="Beta")
        >>> c = SeqRecord(Seq("ACTGCTAGATAG"), id="Gamma")
        >>> align = MultipleSeqAlignment([a, b, c])
        >>> for record in align:
        ...    print(record.id)
        ...    print(record.seq)
        ...
        Alpha
        ACTGCTAGCTAG
        Beta
        ACT-CTAGCTAG
        Gamma
        ACTGCTAGATAG
        """
        return iter(self._records)

    def __len__(self):
        """Return the number of sequences in the alignment.

        Use len(alignment) to get the number of sequences (i.e. the number of
        rows), and alignment.get_alignment_length() to get the length of the
        longest sequence (i.e. the number of columns).

        This is easy to remember if you think of the alignment as being like a
        list of SeqRecord objects.
        """
        return len(self._records)

    def get_alignment_length(self):
        """Return the maximum length of the alignment.

        All objects in the alignment should (hopefully) have the same
        length. This function will go through and find this length
        by finding the maximum length of sequences in the alignment.

        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> from Bio.Align import MultipleSeqAlignment
        >>> a = SeqRecord(Seq("ACTGCTAGCTAG"), id="Alpha")
        >>> b = SeqRecord(Seq("ACT-CTAGCTAG"), id="Beta")
        >>> c = SeqRecord(Seq("ACTGCTAGATAG"), id="Gamma")
        >>> align = MultipleSeqAlignment([a, b, c])
        >>> align.get_alignment_length()
        12

        If you want to know the number of sequences in the alignment,
        use len(align) instead:

        >>> len(align)
        3

        """
        max_length = 0

        for record in self._records:
            if len(record.seq) > max_length:
                max_length = len(record.seq)

        return max_length

    def extend(self, records):
        """Add more SeqRecord objects to the alignment as rows.

        They must all have the same length as the original alignment. For
        example,

        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> from Bio.Align import MultipleSeqAlignment
        >>> a = SeqRecord(Seq("AAAACGT"), id="Alpha")
        >>> b = SeqRecord(Seq("AAA-CGT"), id="Beta")
        >>> c = SeqRecord(Seq("AAAAGGT"), id="Gamma")
        >>> d = SeqRecord(Seq("AAAACGT"), id="Delta")
        >>> e = SeqRecord(Seq("AAA-GGT"), id="Epsilon")

        First we create a small alignment (three rows):

        >>> align = MultipleSeqAlignment([a, b, c])
        >>> print(align)
        Alignment with 3 rows and 7 columns
        AAAACGT Alpha
        AAA-CGT Beta
        AAAAGGT Gamma

        Now we can extend this alignment with another two rows:

        >>> align.extend([d, e])
        >>> print(align)
        Alignment with 5 rows and 7 columns
        AAAACGT Alpha
        AAA-CGT Beta
        AAAAGGT Gamma
        AAAACGT Delta
        AAA-GGT Epsilon

        Because the alignment object allows iteration over the rows as
        SeqRecords, you can use the extend method with a second alignment
        (provided its sequences have the same length as the original alignment).
        """
        if len(self):
            # Use the standard method to get the length
            expected_length = self.get_alignment_length()
        else:
            # Take the first record's length
            records = iter(records)  # records arg could be list or iterator
            try:
                rec = next(records)
            except StopIteration:
                # Special case, no records
                return
            expected_length = len(rec)
            self._append(rec, expected_length)
            # Can now setup the per-column-annotations as well, set to None
            # while missing the length:
            self.column_annotations = {}
            # Now continue to the rest of the records as usual

        for rec in records:
            self._append(rec, expected_length)

    def append(self, record):
        """Add one more SeqRecord object to the alignment as a new row.

        This must have the same length as the original alignment (unless this is
        the first record).

        >>> from Bio import AlignIO
        >>> align = AlignIO.read("Clustalw/opuntia.aln", "clustal")
        >>> print(align)
        Alignment with 7 rows and 156 columns
        TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273285|gb|AF191659.1|AF191
        TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273284|gb|AF191658.1|AF191
        TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273287|gb|AF191661.1|AF191
        TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273286|gb|AF191660.1|AF191
        TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273290|gb|AF191664.1|AF191
        TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273289|gb|AF191663.1|AF191
        TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273291|gb|AF191665.1|AF191
        >>> len(align)
        7

        We'll now construct a dummy record to append as an example:

        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> dummy = SeqRecord(Seq("N"*156), id="dummy")

        Now append this to the alignment,

        >>> align.append(dummy)
        >>> print(align)
        Alignment with 8 rows and 156 columns
        TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273285|gb|AF191659.1|AF191
        TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273284|gb|AF191658.1|AF191
        TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273287|gb|AF191661.1|AF191
        TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273286|gb|AF191660.1|AF191
        TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273290|gb|AF191664.1|AF191
        TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273289|gb|AF191663.1|AF191
        TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273291|gb|AF191665.1|AF191
        NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN...NNN dummy
        >>> len(align)
        8

        """
        if self._records:
            self._append(record, self.get_alignment_length())
        else:
            self._append(record)

    def _append(self, record, expected_length=None):
        """Validate and append a record (PRIVATE)."""
        if not isinstance(record, SeqRecord):
            raise TypeError("New sequence is not a SeqRecord object")

        # Currently the get_alignment_length() call is expensive, so we need
        # to avoid calling it repeatedly for __init__ and extend, hence this
        # private _append method
        if expected_length is not None and len(record) != expected_length:
            # TODO - Use the following more helpful error, but update unit tests
            # raise ValueError("New sequence is not of length %i"
            #                  % self.get_alignment_length())
            raise ValueError("Sequences must all be the same length")

        self._records.append(record)

    def __add__(self, other):
        """Combine two alignments with the same number of rows by adding them.

        If you have two multiple sequence alignments (MSAs), there are two ways to think
        about adding them - by row or by column. Using the extend method adds by row.
        Using the addition operator adds by column. For example,

        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> from Bio.Align import MultipleSeqAlignment
        >>> a1 = SeqRecord(Seq("AAAAC"), id="Alpha")
        >>> b1 = SeqRecord(Seq("AAA-C"), id="Beta")
        >>> c1 = SeqRecord(Seq("AAAAG"), id="Gamma")
        >>> a2 = SeqRecord(Seq("GT"), id="Alpha")
        >>> b2 = SeqRecord(Seq("GT"), id="Beta")
        >>> c2 = SeqRecord(Seq("GT"), id="Gamma")
        >>> left = MultipleSeqAlignment([a1, b1, c1],
        ...                             annotations={"tool": "demo", "name": "start"},
        ...                             column_annotations={"stats": "CCCXC"})
        >>> right = MultipleSeqAlignment([a2, b2, c2],
        ...                             annotations={"tool": "demo", "name": "end"},
        ...                             column_annotations={"stats": "CC"})

        Now, let's look at these two alignments:

        >>> print(left)
        Alignment with 3 rows and 5 columns
        AAAAC Alpha
        AAA-C Beta
        AAAAG Gamma
        >>> print(right)
        Alignment with 3 rows and 2 columns
        GT Alpha
        GT Beta
        GT Gamma

        And add them:

        >>> combined = left + right
        >>> print(combined)
        Alignment with 3 rows and 7 columns
        AAAACGT Alpha
        AAA-CGT Beta
        AAAAGGT Gamma

        For this to work, both alignments must have the same number of records (here
        they both have 3 rows):

        >>> len(left)
        3
        >>> len(right)
        3
        >>> len(combined)
        3

        The individual rows are SeqRecord objects, and these can be added together. Refer
        to the SeqRecord documentation for details of how the annotation is handled. This
        example is a special case in that both original alignments shared the same names,
        meaning when the rows are added they also get the same name.

        Any common annotations are preserved, but differing annotation is lost. This is
        the same behaviour used in the SeqRecord annotations and is designed to prevent
        accidental propagation of inappropriate values:

        >>> combined.annotations
        {'tool': 'demo'}

        Similarly any common per-column-annotations are combined:

        >>> combined.column_annotations
        {'stats': 'CCCXCCC'}

        """
        if not isinstance(other, MultipleSeqAlignment):
            raise NotImplementedError
        if len(self) != len(other):
            raise ValueError(
                "When adding two alignments they must have the same length"
                " (i.e. same number or rows)"
            )
        merged = (left + right for left, right in zip(self, other))
        # Take any common annotation:
        annotations = {}
        for k, v in self.annotations.items():
            if k in other.annotations and other.annotations[k] == v:
                annotations[k] = v
        column_annotations = {}
        for k, v in self.column_annotations.items():
            if k in other.column_annotations:
                column_annotations[k] = v + other.column_annotations[k]
        return MultipleSeqAlignment(
            merged, annotations=annotations, column_annotations=column_annotations
        )

    def __getitem__(self, index):
        """Access part of the alignment.

        Depending on the indices, you can get a SeqRecord object
        (representing a single row), a Seq object (for a single columns),
        a string (for a single characters) or another alignment
        (representing some part or all of the alignment).

        align[r,c] gives a single character as a string
        align[r] gives a row as a SeqRecord
        align[r,:] gives a row as a SeqRecord
        align[:,c] gives a column as a Seq

        align[:] and align[:,:] give a copy of the alignment

        Anything else gives a sub alignment, e.g.
        align[0:2] or align[0:2,:] uses only row 0 and 1
        align[:,1:3] uses only columns 1 and 2
        align[0:2,1:3] uses only rows 0 & 1 and only cols 1 & 2

        We'll use the following example alignment here for illustration:

        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> from Bio.Align import MultipleSeqAlignment
        >>> a = SeqRecord(Seq("AAAACGT"), id="Alpha")
        >>> b = SeqRecord(Seq("AAA-CGT"), id="Beta")
        >>> c = SeqRecord(Seq("AAAAGGT"), id="Gamma")
        >>> d = SeqRecord(Seq("AAAACGT"), id="Delta")
        >>> e = SeqRecord(Seq("AAA-GGT"), id="Epsilon")
        >>> align = MultipleSeqAlignment([a, b, c, d, e])

        You can access a row of the alignment as a SeqRecord using an integer
        index (think of the alignment as a list of SeqRecord objects here):

        >>> first_record = align[0]
        >>> print("%s %s" % (first_record.id, first_record.seq))
        Alpha AAAACGT
        >>> last_record = align[-1]
        >>> print("%s %s" % (last_record.id, last_record.seq))
        Epsilon AAA-GGT

        You can also access use python's slice notation to create a sub-alignment
        containing only some of the SeqRecord objects:

        >>> sub_alignment = align[2:5]
        >>> print(sub_alignment)
        Alignment with 3 rows and 7 columns
        AAAAGGT Gamma
        AAAACGT Delta
        AAA-GGT Epsilon

        This includes support for a step, i.e. align[start:end:step], which
        can be used to select every second sequence:

        >>> sub_alignment = align[::2]
        >>> print(sub_alignment)
        Alignment with 3 rows and 7 columns
        AAAACGT Alpha
        AAAAGGT Gamma
        AAA-GGT Epsilon

        Or to get a copy of the alignment with the rows in reverse order:

        >>> rev_alignment = align[::-1]
        >>> print(rev_alignment)
        Alignment with 5 rows and 7 columns
        AAA-GGT Epsilon
        AAAACGT Delta
        AAAAGGT Gamma
        AAA-CGT Beta
        AAAACGT Alpha

        You can also use two indices to specify both rows and columns. Using simple
        integers gives you the entry as a single character string. e.g.

        >>> align[3, 4]
        'C'

        This is equivalent to:

        >>> align[3][4]
        'C'

        or:

        >>> align[3].seq[4]
        'C'

        To get a single column (as a string) use this syntax:

        >>> align[:, 4]
        'CCGCG'

        Or, to get part of a column,

        >>> align[1:3, 4]
        'CG'

        However, in general you get a sub-alignment,

        >>> print(align[1:5, 3:6])
        Alignment with 4 rows and 3 columns
        -CG Beta
        AGG Gamma
        ACG Delta
        -GG Epsilon

        This should all seem familiar to anyone who has used the NumPy
        array or matrix objects.
        """
        if isinstance(index, int):
            # e.g. result = align[x]
            # Return a SeqRecord
            return self._records[index]
        elif isinstance(index, slice):
            # e.g. sub_align = align[i:j:k]
            new = MultipleSeqAlignment(self._records[index])
            if self.column_annotations and len(new) == len(self):
                # All rows kept (although could have been reversed)
                # Preserve the column annotations too,
                for k, v in self.column_annotations.items():
                    new.column_annotations[k] = v
            return new
        elif len(index) != 2:
            raise TypeError("Invalid index type.")

        # Handle double indexing
        row_index, col_index = index
        if isinstance(row_index, int):
            # e.g. row_or_part_row = align[6, 1:4], gives a SeqRecord
            return self._records[row_index][col_index]
        elif isinstance(col_index, int):
            # e.g. col_or_part_col = align[1:5, 6], gives a string
            return "".join(rec[col_index] for rec in self._records[row_index])
        else:
            # e.g. sub_align = align[1:4, 5:7], gives another alignment
            new = MultipleSeqAlignment(
                rec[col_index] for rec in self._records[row_index]
            )
            if self.column_annotations and len(new) == len(self):
                # All rows kept (although could have been reversed)
                # Preserve the column annotations too,
                for k, v in self.column_annotations.items():
                    new.column_annotations[k] = v[col_index]
            return new

    def sort(self, key=None, reverse=False):
        """Sort the rows (SeqRecord objects) of the alignment in place.

        This sorts the rows alphabetically using the SeqRecord object id by
        default. The sorting can be controlled by supplying a key function
        which must map each SeqRecord to a sort value.

        This is useful if you want to add two alignments which use the same
        record identifiers, but in a different order. For example,

        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> from Bio.Align import MultipleSeqAlignment
        >>> align1 = MultipleSeqAlignment([
        ...              SeqRecord(Seq("ACGT"), id="Human"),
        ...              SeqRecord(Seq("ACGG"), id="Mouse"),
        ...              SeqRecord(Seq("ACGC"), id="Chicken"),
        ...          ])
        >>> align2 = MultipleSeqAlignment([
        ...              SeqRecord(Seq("CGGT"), id="Mouse"),
        ...              SeqRecord(Seq("CGTT"), id="Human"),
        ...              SeqRecord(Seq("CGCT"), id="Chicken"),
        ...          ])

        If you simple try and add these without sorting, you get this:

        >>> print(align1 + align2)
        Alignment with 3 rows and 8 columns
        ACGTCGGT <unknown id>
        ACGGCGTT <unknown id>
        ACGCCGCT Chicken

        Consult the SeqRecord documentation which explains why you get a
        default value when annotation like the identifier doesn't match up.
        However, if we sort the alignments first, then add them we get the
        desired result:

        >>> align1.sort()
        >>> align2.sort()
        >>> print(align1 + align2)
        Alignment with 3 rows and 8 columns
        ACGCCGCT Chicken
        ACGTCGTT Human
        ACGGCGGT Mouse

        As an example using a different sort order, you could sort on the
        GC content of each sequence.

        >>> from Bio.SeqUtils import gc_fraction
        >>> print(align1)
        Alignment with 3 rows and 4 columns
        ACGC Chicken
        ACGT Human
        ACGG Mouse
        >>> align1.sort(key = lambda record: gc_fraction(record.seq))
        >>> print(align1)
        Alignment with 3 rows and 4 columns
        ACGT Human
        ACGC Chicken
        ACGG Mouse

        There is also a reverse argument, so if you wanted to sort by ID
        but backwards:

        >>> align1.sort(reverse=True)
        >>> print(align1)
        Alignment with 3 rows and 4 columns
        ACGG Mouse
        ACGT Human
        ACGC Chicken

        """
        if key is None:
            self._records.sort(key=lambda r: r.id, reverse=reverse)
        else:
            self._records.sort(key=key, reverse=reverse)

    @property
    def substitutions(self):
        """Return an Array with the number of substitutions of letters in the alignment.

        As an example, consider a multiple sequence alignment of three DNA sequences:

        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> from Bio.Align import MultipleSeqAlignment
        >>> seq1 = SeqRecord(Seq("ACGT"), id="seq1")
        >>> seq2 = SeqRecord(Seq("A--A"), id="seq2")
        >>> seq3 = SeqRecord(Seq("ACGT"), id="seq3")
        >>> seq4 = SeqRecord(Seq("TTTC"), id="seq4")
        >>> alignment = MultipleSeqAlignment([seq1, seq2, seq3, seq4])
        >>> print(alignment)
        Alignment with 4 rows and 4 columns
        ACGT seq1
        A--A seq2
        ACGT seq3
        TTTC seq4

        >>> m = alignment.substitutions
        >>> print(m)
            A   C   G   T
        A 3.0 0.5 0.0 2.5
        C 0.5 1.0 0.0 2.0
        G 0.0 0.0 1.0 1.0
        T 2.5 2.0 1.0 1.0
        <BLANKLINE>

        Note that the matrix is symmetric, with counts divided equally on both
        sides of the diagonal. For example, the total number of substitutions
        between A and T in the alignment is 3.5 + 3.5 = 7.

        Any weights associated with the sequences are taken into account when
        calculating the substitution matrix.  For example, given the following
        multiple sequence alignment::

            GTATC  0.5
            AT--C  0.8
            CTGTC  1.0

        For the first column we have::

            ('A', 'G') : 0.5 * 0.8 = 0.4
            ('C', 'G') : 0.5 * 1.0 = 0.5
            ('A', 'C') : 0.8 * 1.0 = 0.8

        """
        letters = set.union(*(set(record.seq) for record in self))
        try:
            letters.remove("-")
        except KeyError:
            pass
        letters = "".join(sorted(letters))
        m = substitution_matrices.Array(letters, dims=2)
        for rec_num1, alignment1 in enumerate(self):
            seq1 = alignment1.seq
            weight1 = alignment1.annotations.get("weight", 1.0)
            for rec_num2, alignment2 in enumerate(self):
                if rec_num1 == rec_num2:
                    break
                seq2 = alignment2.seq
                weight2 = alignment2.annotations.get("weight", 1.0)
                for residue1, residue2 in zip(seq1, seq2):
                    if residue1 == "-":
                        continue
                    if residue2 == "-":
                        continue
                    m[(residue1, residue2)] += weight1 * weight2

        m += m.transpose()
        m /= 2.0

        return m


class Alignment:
    """Represents a sequence alignment.

    An Alignment object has a `.sequences` attribute storing the sequences
    (Seq, MutableSeq, SeqRecord, or string objects) that were aligned, as well
    as a `.coordinates` attribute storing the sequence coordinates defining the
    alignment as a numpy array.

    Other commonly used attributes (which may or may not be present) are:
         - annotations        - A dictionary with annotations describing the
                                alignment;
         - column_annotations - A dictionary with annotations describing each
                                column in the alignment;
         - score              - The alignment score.
    """

    @classmethod
    def infer_coordinates(cls, lines, skipped_columns=None):
        """Infer the coordinates from a printed alignment.

        This method is primarily employed in Biopython's alignment parsers,
        though it may be useful for other purposes.

        For an alignment consisting of N sequences, printed as N lines with
        the same number of columns, where gaps are represented by dashes,
        this method will calculate the sequence coordinates that define the
        alignment. The coordinates are returned as a numpy array of integers,
        and can be used to create an Alignment object.

        The argument skipped columns should be None (the default) or an empty
        list. If skipped_columns is a list, then the indices of any columns in
        the alignment with a gap in all lines are appended to skipped_columns.

        This is an example for the alignment of three sequences TAGGCATACGTG,
        AACGTACGT, and ACGCATACTTG, with gaps in the second and third sequence:

        >>> from Bio.Align import Alignment
        >>> lines = ["TAGGCATACGTG",
        ...          "AACG--TACGT-",
        ...          "-ACGCATACTTG",
        ...         ]
        >>> sequences = [line.replace("-", "") for line in lines]
        >>> sequences
        ['TAGGCATACGTG', 'AACGTACGT', 'ACGCATACTTG']
        >>> coordinates = Alignment.infer_coordinates(lines)
        >>> coordinates
        array([[ 0,  1,  4,  6, 11, 12],
               [ 0,  1,  4,  4,  9,  9],
               [ 0,  0,  3,  5, 10, 11]])
        >>> alignment = Alignment(sequences, coordinates)
        """
        n = len(lines)
        m = len(lines[0])
        for line in lines:
            assert m == len(line)
        path = []
        if m > 0:
            indices = [0] * n
            current_state = [None] * n
            for i in range(m):
                next_state = [line[i] != "-" for line in lines]
                if not any(next_state):
                    # skip columns in which all rows have a gap
                    if skipped_columns is not None:
                        skipped_columns.append(i)
                elif next_state == current_state:
                    step += 1  # noqa: F821
                else:
                    indices = [
                        index + step if state else index
                        for index, state in zip(indices, current_state)
                    ]
                    path.append(indices)
                    step = 1
                    current_state = next_state
            indices = [
                index + step if state else index
                for index, state in zip(indices, current_state)
            ]
            path.append(indices)
        coordinates = numpy.array(path).transpose()
        return coordinates

    def __init__(self, sequences, coordinates=None):
        """Initialize a new Alignment object.

        Arguments:
         - sequences   - A list of the sequences (Seq, MutableSeq, SeqRecord,
                         or string objects) that were aligned.
         - coordinates - The sequence coordinates that define the alignment.
                         If None (the default value), assume that the sequences
                         align to each other without any gaps.
        """
        self.sequences = sequences
        if coordinates is None:
            try:
                lengths = {len(sequence) for sequence in sequences}
            except TypeError:
                # this may happen if sequences contain a SeqRecord where
                # the seq attribute is None, as neither the sequence nor
                # its length are known.
                pass
            else:
                if len(lengths) == 0:
                    coordinates = numpy.empty((0, 0), dtype=int)
                elif len(lengths) == 1:
                    length = lengths.pop()
                    coordinates = numpy.array([[0, length]] * len(sequences))
                else:
                    raise ValueError(
                        "sequences must have the same length if coordinates is None"
                    )
        self.coordinates = coordinates

    def __array__(self, dtype=None):
        coordinates = self.coordinates.copy()
        sequences = list(self.sequences)
        steps = numpy.diff(self.coordinates, 1)
        aligned = sum(steps != 0, 0) > 1
        # True for steps in which at least two sequences align, False if a gap
        for i, sequence in enumerate(sequences):
            row = steps[i, aligned]
            if (row >= 0).all():
                pass
            elif (row <= 0).all():
                sequences[i] = reverse_complement(sequence, inplace=False)
                coordinates[i, :] = len(sequence) - coordinates[i, :]
                steps[i, :] = -steps[i, :]
            else:
                raise ValueError(f"Inconsistent steps in row {i}")
        gaps = steps.max(0)
        if not ((steps == gaps) | (steps <= 0)).all():
            raise ValueError("Unequal step sizes in alignment")
        n = len(steps)
        m = sum(gaps)
        data = numpy.empty((n, m), "S1")
        for i in range(n):
            sequence = sequences[i]
            k = coordinates[i, 0]
            m = 0
            for step, gap in zip(steps[i], gaps):
                if step > 0:
                    j = k + step
                    n = m + step
                    try:
                        subsequence = bytes(sequence[k:j])
                    except TypeError:  # str
                        subsequence = bytes(sequence[k:j], "UTF8")
                    data[i, :].data.cast("B")[m:n] = subsequence
                    k = j
                    m = n
                elif step < 0:
                    k += step
                else:  # step == 0
                    n = m + gap
                    data[i, m:n] = b"-"
                    m = n
        if dtype is not None:
            data = numpy.array(data, dtype)
        return data

    @property
    def target(self):
        """Return self.sequences[0] for a pairwise alignment."""
        n = len(self.sequences)
        if n != 2:
            raise ValueError(
                "self.target is defined for pairwise alignments only (found alignment of %d sequences)"
                % n
            )
        return self.sequences[0]

    @target.setter
    def target(self, value):
        """For a pairwise alignment, set self.sequences[0]."""
        n = len(self.sequences)
        if n != 2:
            raise ValueError(
                "self.target is defined for pairwise alignments only (found alignment of %d sequences)"
                % n
            )
        self.sequences[0] = value

    @property
    def query(self):
        """Return self.sequences[1] for a pairwise alignment."""
        n = len(self.sequences)
        if n != 2:
            raise ValueError(
                "self.query is defined for pairwise alignments only (found alignment of %d sequences)"
                % n
            )
        return self.sequences[1]

    @query.setter
    def query(self, value):
        """For a pairwise alignment, set self.sequences[1]."""
        n = len(self.sequences)
        if n != 2:
            raise ValueError(
                "self.query is defined for pairwise alignments only (found alignment of %d sequences)"
                % n
            )
        self.sequences[1] = value

    def __eq__(self, other):
        """Check if two Alignment objects specify the same alignment."""
        for left, right in zip_longest(self.sequences, other.sequences):
            try:
                left = left.seq
            except AttributeError:
                pass
            try:
                right = right.seq
            except AttributeError:
                pass
            if left != right:
                return False
        return numpy.array_equal(self.coordinates, other.coordinates)

    def __ne__(self, other):
        """Check if two Alignment objects have different alignments."""
        for left, right in zip_longest(self.sequences, other.sequences):
            try:
                left = left.seq
            except AttributeError:
                pass
            try:
                right = right.seq
            except AttributeError:
                pass
            if left != right:
                return True

        return not numpy.array_equal(self.coordinates, other.coordinates)

    def __lt__(self, other):
        """Check if self should come before other."""
        for left, right in zip_longest(self.sequences, other.sequences):
            try:
                left = left.seq
            except AttributeError:
                pass
            try:
                right = right.seq
            except AttributeError:
                pass
            if left < right:
                return True
            if left > right:
                return False
        for left, right in zip(
            self.coordinates.transpose(), other.coordinates.transpose()
        ):
            left, right = tuple(left), tuple(right)
            if left < right:
                return True
            if left > right:
                return False
        return False

    def __le__(self, other):
        """Check if self should come before or is equal to other."""
        for left, right in zip_longest(self.sequences, other.sequences):
            try:
                left = left.seq
            except AttributeError:
                pass
            try:
                right = right.seq
            except AttributeError:
                pass
            if left < right:
                return True
            if left > right:
                return False
        for left, right in zip(
            self.coordinates.transpose(), other.coordinates.transpose()
        ):
            left, right = tuple(left), tuple(right)
            if left < right:
                return True
            if left > right:
                return False
        return True

    def __gt__(self, other):
        """Check if self should come after other."""
        for left, right in zip_longest(self.sequences, other.sequences):
            try:
                left = left.seq
            except AttributeError:
                pass
            try:
                right = right.seq
            except AttributeError:
                pass
            if left < right:
                return False
            if left > right:
                return True
        for left, right in zip(
            self.coordinates.transpose(), other.coordinates.transpose()
        ):
            left, right = tuple(left), tuple(right)
            if left > right:
                return True
            if left < right:
                return False
        return False

    def __ge__(self, other):
        """Check if self should come after or is equal to other."""
        for left, right in zip_longest(self.sequences, other.sequences):
            try:
                left = left.seq
            except AttributeError:
                pass
            try:
                right = right.seq
            except AttributeError:
                pass
            if left < right:
                return False
            if left > right:
                return True
        for left, right in zip(
            self.coordinates.transpose(), other.coordinates.transpose()
        ):
            left, right = tuple(left), tuple(right)
            if left > right:
                return True
            if left < right:
                return False
        return True

    @property
    def path(self):
        """Return the path through the trace matrix."""
        warnings.warn(
            "The path attribute is deprecated; please use the coordinates "
            "attribute instead. The coordinates attribute is a numpy array "
            "containing the same values as the path attributes, after "
            "transposition.",
            BiopythonDeprecationWarning,
        )
        return tuple(tuple(row) for row in self.coordinates.transpose())

    @path.setter
    def path(self, value):
        warnings.warn(
            "The path attribute is deprecated; please use the coordinates "
            "attribute instead. The coordinates attribute is a numpy array "
            "containing the same values as the path attributes, after "
            "transposition.",
            BiopythonDeprecationWarning,
        )
        self.coordinates = numpy.array(value).transpose()

    def _get_row(self, index):
        """Return self[index], where index is an integer (PRIVATE).

        This method is called by __getitem__ for invocations of the form

        self[row]

        where row is an integer.
        Return value is a string if the aligned sequences are string, Seq,
        or SeqRecord objects, otherwise the return value is a list.
        """
        steps = numpy.diff(self.coordinates, 1)
        n = len(steps)
        if index < 0:
            index += n
            if index < 0:
                raise IndexError("row index out of range")
        elif index >= n:
            raise IndexError("row index out of range")
        aligned = sum(steps != 0, 0) > 1
        # True for steps in which at least two sequences align, False if a gap
        coordinates = self.coordinates[index, :]
        sequence = self.sequences[index]
        for i in range(n):
            row = steps[i, aligned]
            if (row >= 0).all():
                pass
            elif (row <= 0).all():
                steps[i, :] = -steps[i, :]
                if i == index:
                    sequence = reverse_complement(sequence, inplace=False)
                    coordinates = len(sequence) - coordinates
            else:
                raise ValueError(f"Inconsistent steps in row {index}")
        gaps = steps.max(0)
        if not ((steps == gaps) | (steps <= 0)).all():
            raise ValueError("Unequal step sizes in alignment")
        try:
            sequence = sequence.seq  # SeqRecord confusion
        except AttributeError:
            pass
        steps = steps[index]
        k = coordinates[0]
        if isinstance(sequence, (str, Seq)):
            line = ""
            for step, gap in zip(steps, gaps):
                if step > 0:
                    j = k + step
                    line += str(sequence[k:j])
                    k = j
                elif step < 0:
                    k += step
                else:  # step == 0
                    line += "-" * gap
        else:
            line = []
            for step, gap in zip(steps, gaps):
                if step > 0:
                    j = k + step
                    line.extend(sequence[k:j])
                    k = j
                else:
                    line.extend([None] * gap)
        return line

    def _get_rows(self, key):
        """Return self[key], where key is a slice object (PRIVATE).

        This method is called by __getitem__ for invocations of the form

        self[rows]

        where rows is a slice object. Return value is an Alignment object.
        """
        sequences = self.sequences[key]
        coordinates = self.coordinates[key].copy()
        alignment = Alignment(sequences, coordinates)
        if numpy.array_equal(self.coordinates, coordinates):
            try:
                alignment.score = self.score
            except AttributeError:
                pass
            try:
                alignment.column_annotations = self.column_annotations
            except AttributeError:
                pass
        return alignment

    def _get_row_col(self, j, col, steps, gaps, sequence):
        """Return the sequence contents at alignment column j (PRIVATE).

        This method is called by __getitem__ for invocations of the form

        self[row, col]

        where both row and col are integers.
        Return value is a string of length 1.
        """
        indices = gaps.cumsum()
        index = indices.searchsorted(col, side="right")
        if steps[index]:
            offset = col - indices[index]
            j += sum(steps[: index + 1]) + offset
            return sequence[j]
        else:
            return "-"

    def _get_row_cols_slice(
        self, coordinate, start_index, stop_index, steps, gaps, sequence
    ):
        """Return the alignment contents of one row and consecutive columns (PRIVATE).

        This method is called by __getitem__ for invocations of the form

        self[row, cols]

        where row is an integer and cols is a slice object with step 1.
        Return value is a string if the aligned sequences are string, Seq,
        or SeqRecord objects, otherwise the return value is a list.
        """
        indices = gaps.cumsum()
        i = indices.searchsorted(start_index, side="right")
        j = i + indices[i:].searchsorted(stop_index, side="right")
        try:
            sequence = sequence.seq  # stupid SeqRecord
        except AttributeError:
            pass
        if isinstance(sequence, (str, Seq)):
            if i == j:
                length = stop_index - start_index
                if steps[i] == 0:
                    line = "-" * length
                else:
                    start = coordinate[i] + start_index - indices[i - 1]
                    stop = start + length
                    line = str(sequence[start:stop])
            else:
                length = indices[i] - start_index
                if steps[i] == 0:
                    line = "-" * length
                else:
                    stop = coordinate[i + 1]
                    start = stop - length
                    line = str(sequence[start:stop])
                i += 1
                while i < j:
                    step = gaps[i]
                    if steps[i] == 0:
                        line += "-" * step
                    else:
                        start = coordinate[i]
                        stop = coordinate[i + 1]
                        line += str(sequence[start:stop])
                    i += 1
                length = stop_index - indices[i - 1]
                if length > 0:
                    if steps[i] == 0:
                        line += "-" * length
                    else:
                        start = coordinate[i]
                        stop = start + length
                        line += str(sequence[start:stop])
        else:
            if i == j:
                length = stop_index - start_index
                if steps[i] == 0:
                    line = [None] * length
                else:
                    start = coordinate[i] + start_index - indices[i - 1]
                    stop = start + length
                    line = sequence[start:stop]
            else:
                length = indices[i] - start_index
                if steps[i] == 0:
                    line = [None] * length
                else:
                    stop = coordinate[i + 1]
                    start = stop - length
                    line = sequence[start:stop]
                i += 1
                while i < j:
                    step = gaps[i]
                    if steps[i] == 0:
                        line.extend([None] * step)
                    else:
                        start = coordinate[i]
                        stop = coordinate[i + 1]
                        line.extend(sequence[start:stop])
                    i += 1
                length = stop_index - indices[i - 1]
                if length > 0:
                    if steps[j] == 0:
                        line.extend([None] * length)
                    else:
                        start = coordinate[i]
                        stop = start + length
                        line.extend(sequence[start:stop])
        return line

    def _get_row_cols_iterable(self, coordinate, cols, gaps, sequence):
        """Return the alignment contents of one row and multiple columns (PRIVATE).

        This method is called by __getitem__ for invocations of the form

        self[row, cols]

        where row is an integer and cols is an iterable of integers.
        Return value is a string if the aligned sequences are string, Seq,
        or SeqRecord objects, otherwise the return value is a list.
        """
        try:
            sequence = sequence.seq  # stupid SeqRecord
        except AttributeError:
            pass
        if isinstance(sequence, (str, Seq)):
            line = ""
            start = coordinate[0]
            for end, gap in zip(coordinate[1:], gaps):
                if start < end:
                    line += str(sequence[start:end])
                else:
                    line += "-" * gap
                start = end
            try:
                line = "".join(line[col] for col in cols)
            except IndexError:
                raise
            except Exception:
                raise TypeError(
                    "second index must be an integer, slice, or iterable of integers"
                ) from None
        else:
            line = []
            start = coordinate[0]
            for end, gap in zip(coordinate[1:], gaps):
                if start < end:
                    line.extend(sequence[start:end])
                else:
                    line.extend([None] * gap)
                start = end
            try:
                line = [line[col] for col in cols]
            except IndexError:
                raise
            except Exception:
                raise TypeError(
                    "second index must be an integer, slice, or iterable of integers"
                ) from None
        return line

    def _get_rows_col(self, coordinates, col, steps, gaps, sequences):
        """Return the alignment contents of multiple rows and one column (PRIVATE).

        This method is called by __getitem__ for invocations of the form

        self[rows, col]

        where rows is a slice object, and col is an integer.
        Return value is a string.
        """
        indices = gaps.cumsum()
        j = indices.searchsorted(col, side="right")
        offset = indices[j] - col
        line = ""
        for sequence, coordinate, step in zip(sequences, coordinates, steps):
            if step[j] == 0:
                line += "-"
            else:
                index = coordinate[j] + step[j] - offset
                line += sequence[index]
        return line

    def _get_rows_cols_slice(
        self, coordinates, row, start_index, stop_index, steps, gaps
    ):
        """Return a subalignment of multiple rows and consecutive columns (PRIVATE).

        This method is called by __getitem__ for invocations of the form

        self[rows, cols]

        where rows is an arbitrary slice object, and cols is a slice object
        with step 1, allowing the alignment sequences to be reused in the
        subalignment. Return value is an Alignment object.
        """
        rcs = numpy.any(coordinates != self.coordinates[row], axis=1)
        indices = gaps.cumsum()
        i = indices.searchsorted(start_index, side="right")
        j = i + indices[i:].searchsorted(stop_index, side="left") + 1
        offset = steps[:, i] - indices[i] + start_index
        coordinates[:, i] += offset * (steps[:, i] > 0)
        offset = indices[j - 1] - stop_index
        coordinates[:, j] -= offset * (steps[:, j - 1] > 0)
        coordinates = coordinates[:, i : j + 1]
        sequences = self.sequences[row]
        for coordinate, rc, sequence in zip(coordinates, rcs, sequences):
            if rc:
                # mapped to reverse strand
                coordinate[:] = len(sequence) - coordinate[:]
        alignment = Alignment(sequences, coordinates)
        if numpy.array_equal(self.coordinates, coordinates):
            try:
                alignment.score = self.score
            except AttributeError:
                pass
        try:
            column_annotations = self.column_annotations
        except AttributeError:
            pass
        else:
            alignment.column_annotations = {}
            for key, value in column_annotations.items():
                value = value[start_index:stop_index]
                try:
                    value = value.copy()
                except AttributeError:
                    # immutable tuples like str, tuple
                    pass
                alignment.column_annotations[key] = value
        return alignment

    def _get_rows_cols_iterable(self, coordinates, col, steps, gaps, sequences):
        """Return a subalignment of multiple rows and columns (PRIVATE).

        This method is called by __getitem__ for invocations of the form

        self[rows, cols]

        where rows is a slice object and cols is an iterable of integers.
        This method will create new sequences for use by the subalignment
        object. Return value is an Alignment object.
        """
        indices = tuple(col)
        lines = []
        for i, sequence in enumerate(sequences):
            try:
                s = sequence.seq  # stupid SeqRecord
            except AttributeError:
                s = sequence
            line = ""
            k = coordinates[i, 0]
            for step, gap in zip(steps[i], gaps):
                if step:
                    j = k + step
                    line += str(s[k:j])
                    k = j
                else:
                    line += "-" * gap
            try:
                line = "".join(line[index] for index in indices)
            except IndexError:
                raise
            except Exception:
                raise TypeError(
                    "second index must be an integer, slice, or iterable of integers"
                ) from None
            lines.append(line)
            line = line.replace("-", "")
            s = s.__class__(line)
            try:
                sequence.seq  # stupid SeqRecord
            except AttributeError:
                sequence = s
            else:
                sequence = copy.deepcopy(sequence)
                sequence.seq = s
            sequences[i] = sequence
        coordinates = self.infer_coordinates(lines)
        alignment = Alignment(sequences, coordinates)
        try:
            column_annotations = self.column_annotations
        except AttributeError:
            pass
        else:
            alignment.column_annotations = {}
            for key, value in column_annotations.items():
                values = (value[index] for index in indices)
                if isinstance(value, str):
                    value = "".join(values)
                else:
                    value = value.__class__(values)
                alignment.column_annotations[key] = value
        return alignment

    def __getitem__(self, key):
        """Return self[key].

        Indices of the form

        self[:, :]

        return a copy of the Alignment object;

        self[:, i:]
        self[:, :j]
        self[:, i:j]
        self[:, iterable] (where iterable returns integers)

        return a new Alignment object spanning the selected columns;

        self[k, i]
        self[k, i:]
        self[k, :j]
        self[k, i:j]
        self[k, iterable] (where iterable returns integers)
        self[k] (equivalent to self[k, :])

        return a string with the aligned sequence (including gaps) for the
        selected columns, where k = 0 represents the target and k = 1
        represents the query sequence; and

        self[:, i]

        returns a string with the selected column in the alignment.

        >>> from Bio.Align import PairwiseAligner
        >>> aligner = PairwiseAligner()
        >>> alignments = aligner.align("ACCGGTTT", "ACGGGTT")
        >>> alignment = alignments[0]
        >>> print(alignment)
        target            0 ACCGG-TTT 8
                          0 ||-||-||- 9
        query             0 AC-GGGTT- 7
        <BLANKLINE>
        >>> alignment[0, :]
        'ACCGG-TTT'
        >>> alignment[1, :]
        'AC-GGGTT-'
        >>> alignment[0]
        'ACCGG-TTT'
        >>> alignment[1]
        'AC-GGGTT-'
        >>> alignment[0, 1:-2]
        'CCGG-T'
        >>> alignment[1, 1:-2]
        'C-GGGT'
        >>> alignment[0, (1, 5, 2)]
        'C-C'
        >>> alignment[1, ::2]
        'A-GT-'
        >>> alignment[1, range(0, 9, 2)]
        'A-GT-'
        >>> alignment[:, 0]
        'AA'
        >>> alignment[:, 5]
        '-G'
        >>> alignment[:, 1:]  # doctest:+ELLIPSIS
        <Alignment object (2 rows x 8 columns) at 0x...>
        >>> print(alignment[:, 1:])
        target            1 CCGG-TTT 8
                          0 |-||-||- 8
        query             1 C-GGGTT- 7
        <BLANKLINE>
        >>> print(alignment[:, 2:])
        target            2 CGG-TTT 8
                          0 -||-||- 7
        query             2 -GGGTT- 7
        <BLANKLINE>
        >>> print(alignment[:, 3:])
        target            3 GG-TTT 8
                          0 ||-||- 6
        query             2 GGGTT- 7
        <BLANKLINE>
        >>> print(alignment[:, 3:-1])
        target            3 GG-TT 7
                          0 ||-|| 5
        query             2 GGGTT 7
        <BLANKLINE>
        >>> print(alignment[:, ::2])
        target            0 ACGTT 5
                          0 |-||- 5
        query             0 A-GT- 3
        <BLANKLINE>
        >>> print(alignment[:, range(1, 9, 2)])
        target            0 CG-T 3
                          0 ||-| 4
        query             0 CGGT 4
        <BLANKLINE>
        >>> print(alignment[:, (2, 7, 3)])
        target            0 CTG 3
                          0 -|| 3
        query             0 -TG 2
        <BLANKLINE>
        """
        if isinstance(key, numbers.Integral):
            return self._get_row(key)
        if isinstance(key, slice):
            return self._get_rows(key)
        sequences = list(self.sequences)
        coordinates = self.coordinates.copy()
        steps = numpy.diff(coordinates, 1)
        aligned = sum(steps != 0, 0) > 1
        # True for steps in which at least two sequences align, False if a gap
        for i, sequence in enumerate(sequences):
            row = steps[i, aligned]
            if (row >= 0).all():
                pass
            elif (row <= 0).all():
                steps[i, :] = -steps[i, :]
                coordinates[i, :] = len(sequence) - coordinates[i, :]
                sequences[i] = reverse_complement(sequence, inplace=False)
                try:
                    sequences[i].id = sequence.id
                except AttributeError:
                    pass
            else:
                raise ValueError(f"Inconsistent steps in row {i}")
        gaps = steps.max(0)
        if not ((steps == gaps) | (steps <= 0)).all():
            raise ValueError("Unequal step sizes in alignment")
        m = sum(gaps)
        if isinstance(key, tuple):
            try:
                row, col = key
            except ValueError:
                raise ValueError("only tuples of length 2 can be alignment indices")
        else:
            raise TypeError("alignment indices must be integers, slices, or tuples")
        if isinstance(col, numbers.Integral):
            if col < 0:
                col += m
            if col < 0 or col >= m:
                raise IndexError(
                    "column index %d is out of bounds (%d columns)" % (col, m)
                )
        steps = steps[row]
        if isinstance(row, numbers.Integral):
            sequence = sequences[row]
            if isinstance(col, numbers.Integral):
                return self._get_row_col(
                    coordinates[row, 0], col, steps, gaps, sequence
                )
            coordinate = coordinates[row, :]
            if isinstance(col, slice):
                start_index, stop_index, step = col.indices(m)
                if start_index < stop_index and step == 1:
                    return self._get_row_cols_slice(
                        coordinate, start_index, stop_index, steps, gaps, sequence
                    )
                # make an iterable if step != 1
                col = range(start_index, stop_index, step)
            return self._get_row_cols_iterable(coordinate, col, gaps, sequence)
        if isinstance(row, slice):
            sequences = sequences[row]
            coordinates = coordinates[row]
            if isinstance(col, numbers.Integral):
                return self._get_rows_col(coordinates, col, steps, gaps, sequences)
            if isinstance(col, slice):
                start_index, stop_index, step = col.indices(m)
                if start_index < stop_index and step == 1:
                    return self._get_rows_cols_slice(
                        coordinates,
                        row,
                        start_index,
                        stop_index,
                        steps,
                        gaps,
                    )
                # make an iterable if step != 1
                col = range(start_index, stop_index, step)
            # try if we can use col as an iterable
            return self._get_rows_cols_iterable(
                coordinates, col, steps, gaps, sequences
            )
        raise TypeError("first index must be an integer or slice")

    def _convert_sequence_string(self, sequence):
        """Convert given sequence to string using the appropriate method (PRIVATE)."""
        if isinstance(sequence, (bytes, bytearray)):
            return sequence.decode()
        if isinstance(sequence, str):
            return sequence
        if isinstance(sequence, Seq):
            return str(sequence)
        try:  # check if target is a SeqRecord
            sequence = sequence.seq
        except AttributeError:
            pass
        else:
            return str(sequence)
        try:
            view = memoryview(sequence)
        except TypeError:
            pass
        else:
            if view.format == "c":
                return str(sequence)
        return None

    def __format__(self, format_spec):
        """Return the alignment as a string in the specified file format.

        Wrapper for self.format().
        """
        return self.format(format_spec)

    def format(self, fmt="", *args, **kwargs):
        """Return the alignment as a string in the specified file format.

        Arguments:
         - fmt       - File format. Acceptable values are an empty string to
                       create a human-readable representation of the alignment,
                       or any of the alignment file formats supported by
                       `Bio.Align` (some have not yet been implemented).

        All other arguments are passed to the format-specific writer functions:
         - mask      - PSL format only. Specify if repeat regions in the target
                       sequence are masked and should be reported in the
                       `repMatches` field of the PSL file instead of in the
                       `matches` field. Acceptable values are
                       None   : no masking (default);
                       "lower": masking by lower-case characters;
                       "upper": masking by upper-case characters.
         - wildcard  - PSL format only. Report alignments to the wildcard
                       character in the target or query sequence in the
                       `nCount` field of the PSL file instead of in the
                       `matches`, `misMatches`, or `repMatches` fields.
                       Default value is 'N'.
         - md        - SAM format only. If True, calculate the MD tag from
                       the alignment and include it in the output. If False
                       (default), do not include the MD tag in the output.
        """
        if fmt == "":
            return self._format_pretty()
        module = _load(fmt)
        try:
            writer = module.AlignmentWriter(None, *args, **kwargs)
        except AttributeError:
            if module.AlignmentIterator.mode == "b":
                raise ValueError(f"{fmt} is a binary file format")
            raise ValueError(
                f"Formatting alignments has not yet been implemented for the {fmt} format"
            ) from None
        return writer.format_alignment(self)

    def _format_pretty(self):
        """Return default string representation (PRIVATE).

        Helper for self.format().
        """
        n = len(self.sequences)
        if n == 2:
            write_pattern = True
        else:
            write_pattern = False
        steps = numpy.diff(self.coordinates, 1)
        aligned = sum(steps != 0, 0) > 1
        # True for steps in which at least two sequences align, False if a gap
        name_width = 10
        names = []
        seqs = []
        indices = numpy.zeros(self.coordinates.shape, int)
        for i, (seq, positions, row) in enumerate(
            zip(self.sequences, self.coordinates, indices)
        ):
            try:
                name = seq.id
                if name is None:
                    raise AttributeError
            except AttributeError:
                if n == 2:
                    if i == 0:
                        name = "target"
                    else:
                        name = "query"
                else:
                    name = ""
            else:
                name = name[: name_width - 1]
            name = name.ljust(name_width)
            names.append(name)
            try:
                seq = seq.seq  # SeqRecord confusion
            except AttributeError:
                pass
            start = min(positions)
            end = max(positions)
            seq = seq[start:end]
            aligned_steps = steps[i, aligned]
            if len(aligned_steps) == 0:
                aligned_steps = steps[i]
            if (aligned_steps >= 0).all():
                start = min(positions)
                row[:] = positions - start
            elif (aligned_steps <= 0).all():
                steps[i, :] = -steps[i, :]
                seq = reverse_complement(seq, inplace=False)
                end = max(positions)
                row[:] = end - positions
            else:
                raise ValueError(f"Inconsistent steps in row {i}")
            if isinstance(seq, str):
                if not seq.isascii():
                    return self._format_unicode()
            elif isinstance(seq, (Seq, MutableSeq)):
                try:
                    seq = bytes(seq)
                except UndefinedSequenceError:
                    s = bytearray(b"?" * (end - start))
                    for start, end in seq.defined_ranges:
                        s[start:end] = bytes(seq[start:end])
                    seq = s
                seq = seq.decode()
            else:
                return self._format_generalized()
            seqs.append(seq)
        minstep = steps.min(0)
        maxstep = steps.max(0)
        steps = numpy.where(-minstep > maxstep, minstep, maxstep)
        for i, row in enumerate(indices):
            row_steps = numpy.diff(row)
            row_aligned = (row_steps > 0) & aligned
            row_steps = row_steps[row_aligned]
            aligned_steps = steps[row_aligned]
            if (row_steps == aligned_steps).all():
                pass
            elif (3 * row_steps == aligned_steps).all():
                row[:] *= 3
                seqs[i] = "  ".join(seqs[i]) + "  "
                write_pattern = False
            else:
                raise ValueError("Inconsistent coordinates")
        prefix_width = 10
        position_width = 10
        line_width = 80
        lines = []
        steps = indices[:, 1:] - indices[:, :-1]
        minstep = steps.min(0)
        maxstep = steps.max(0)
        steps = numpy.where(-minstep > maxstep, minstep, maxstep)
        for name, seq, positions, row in zip(names, seqs, self.coordinates, indices):
            start = positions[0]
            column = line_width
            start_index = row[0]
            for step, end, end_index in zip(steps, positions[1:], row[1:]):
                if step < 0:
                    if prefix_width + position_width < column:
                        position_text = str(start)
                        offset = position_width - len(position_text) - 1
                        if offset < 0:
                            lines[-1] += " .." + position_text[-offset + 3 :]
                        else:
                            lines[-1] += " " + position_text
                    column = line_width
                    start = end
                    start_index = end_index
                    continue
                elif end_index == start_index:
                    s = "-" * step
                else:
                    s = seq[start_index:end_index]
                while column + len(s) >= line_width:
                    rest = line_width - column
                    if rest > 0:
                        lines[-1] += s[:rest]
                        s = s[rest:]
                        if start != end:
                            if (end_index - start_index) == abs(end - start):
                                step = rest
                            else:
                                # protein to dna alignment;
                                # integer division, but round up:
                                step = -(rest // -3)
                            if start < end:
                                start += step
                            else:
                                start -= step
                        start_index += rest
                    line = name
                    position_text = str(start)
                    offset = position_width - len(position_text) - 1
                    if offset < 0:
                        line += " .." + position_text[-offset + 3 :]
                    else:
                        line += " " * offset + position_text
                    line += " "
                    lines.append(line)
                    column = name_width + position_width
                lines[-1] += s
                if start_index != end_index:
                    start_index = end_index
                    start = end
                column += len(s)
        if write_pattern is True:
            dash = "-"
            position = 0
            m = len(lines) // 2
            lines1 = lines[:m]
            lines2 = lines[m:]
            pattern_lines = []
            for line1, line2 in zip(lines1, lines2):
                aligned_seq1 = line1[name_width + position_width :]
                aligned_seq2 = line2[name_width + position_width :]
                pattern = ""
                for c1, c2 in zip(aligned_seq1, aligned_seq2):
                    if c1 == c2:
                        if c1 == " ":
                            break
                        c = "|"
                    elif c1 == dash or c2 == dash:
                        c = "-"
                    else:
                        c = "."
                    pattern += c
                pattern_line = "          %9d %s" % (position, pattern)
                pattern_lines.append(pattern_line)
                position += len(pattern)
            final_position_width = len(str(max(max(self.coordinates[:, -1]), position)))
            if column + final_position_width <= line_width:
                if prefix_width + position_width < column:
                    fmt = f" %{final_position_width}d"
                    lines1[-1] += fmt % self.coordinates[0, -1]
                    lines2[-1] += fmt % self.coordinates[1, -1]
                    pattern_lines[-1] += fmt % position
            else:
                name1, name2 = names
                fmt = "%s%9d"
                line = name1 + format(self.coordinates[0, -1], "9d")
                lines1.append(line)
                line = fmt % ("          ", position)
                pattern_lines.append(line)
                line = fmt % (name2, self.coordinates[1, -1])
                lines2.append(line)
                lines.append("")
            return "\n".join(
                f"{line1}\n{pattern_line}\n{line2}\n"
                for (line1, line2, pattern_line) in zip(lines1, lines2, pattern_lines)
            )
        else:
            m = len(lines) // n
            final_position_width = len(str(max(self.coordinates[:, -1])))
            if column + final_position_width < line_width:
                if prefix_width + position_width < column:
                    fmt = f" %{final_position_width}d"
                    for i in range(n):
                        lines[m - 1 + i * m] += fmt % self.coordinates[i, -1]
                blocks = ["\n".join(lines[j::m]) + "\n" for j in range(m)]
            else:
                blocks = ["\n".join(lines[j::m]) + "\n" for j in range(m)]
                lines = []
                fmt = "%s%9d"
                for i in range(n):
                    line = names[i] + format(self.coordinates[i, -1], "9d")
                    lines.append(line)
                block = "\n".join(lines) + "\n"
                blocks.append(block)
            return "\n".join(blocks)

    def _format_unicode(self):
        """Return default string representation (PRIVATE).

        Helper for self.format().
        """
        seqs = []
        names = []
        coordinates = self.coordinates.copy()
        for seq, row in zip(self.sequences, coordinates):
            seq = self._convert_sequence_string(seq)
            if seq is None:
                return self._format_generalized()
            if row[0] > row[-1]:  # mapped to reverse strand
                row[:] = len(seq) - row[:]
                seq = reverse_complement(seq, inplace=False)
            seqs.append(seq)
            try:
                name = seq.id
            except AttributeError:
                if len(self.sequences) == 2:
                    if len(names) == 0:
                        name = "target"
                    else:
                        name = "query"
                else:
                    name = ""
            else:
                name = name[:9]
            name = name.ljust(10)
            names.append(name)
        steps = numpy.diff(coordinates, 1).max(0)
        aligned_seqs = []
        for row, seq in zip(coordinates, seqs):
            aligned_seq = ""
            start = row[0]
            for step, end in zip(steps, row[1:]):
                if end == start:
                    aligned_seq += "-" * step
                else:
                    aligned_seq += seq[start:end]
                start = end
            aligned_seqs.append(aligned_seq)
        if len(seqs) > 2:
            return "\n".join(aligned_seqs) + "\n"
        aligned_seq1, aligned_seq2 = aligned_seqs
        pattern = ""
        for c1, c2 in zip(aligned_seq1, aligned_seq2):
            if c1 == c2:
                c = "|"
            elif c1 == "-" or c2 == "-":
                c = "-"
            else:
                c = "."
            pattern += c
        return f"{aligned_seq1}\n{pattern}\n{aligned_seq2}\n"

    def _format_generalized(self):
        """Return generalized string representation (PRIVATE).

        Helper for self._format_pretty().
        """
        seq1, seq2 = self.sequences
        aligned_seq1 = []
        aligned_seq2 = []
        pattern = []
        end1, end2 = self.coordinates[:, 0]
        if end1 > 0 or end2 > 0:
            if end1 <= end2:
                for c2 in seq2[: end2 - end1]:
                    s2 = str(c2)
                    s1 = " " * len(s2)
                    aligned_seq1.append(s1)
                    aligned_seq2.append(s2)
                    pattern.append(s1)
            else:  # end1 > end2
                for c1 in seq1[: end1 - end2]:
                    s1 = str(c1)
                    s2 = " " * len(s1)
                    aligned_seq1.append(s1)
                    aligned_seq2.append(s2)
                    pattern.append(s2)
        start1 = end1
        start2 = end2
        for end1, end2 in self.coordinates[:, 1:].transpose():
            if end1 == start1:
                for c2 in seq2[start2:end2]:
                    s2 = str(c2)
                    s1 = "-" * len(s2)
                    aligned_seq1.append(s1)
                    aligned_seq2.append(s2)
                    pattern.append(s1)
                start2 = end2
            elif end2 == start2:
                for c1 in seq1[start1:end1]:
                    s1 = str(c1)
                    s2 = "-" * len(s1)
                    aligned_seq1.append(s1)
                    aligned_seq2.append(s2)
                    pattern.append(s2)
                start1 = end1
            else:
                t1 = seq1[start1:end1]
                t2 = seq2[start2:end2]
                if len(t1) != len(t2):
                    raise ValueError("Unequal step sizes in alignment")
                for c1, c2 in zip(t1, t2):
                    s1 = str(c1)
                    s2 = str(c2)
                    m1 = len(s1)
                    m2 = len(s2)
                    if c1 == c2:
                        p = "|"
                    else:
                        p = "."
                    if m1 < m2:
                        space = (m2 - m1) * " "
                        s1 += space
                        pattern.append(p * m1 + space)
                    elif m1 > m2:
                        space = (m1 - m2) * " "
                        s2 += space
                        pattern.append(p * m2 + space)
                    else:
                        pattern.append(p * m1)
                    aligned_seq1.append(s1)
                    aligned_seq2.append(s2)
                start1 = end1
                start2 = end2
        aligned_seq1 = " ".join(aligned_seq1)
        aligned_seq2 = " ".join(aligned_seq2)
        pattern = " ".join(pattern)
        return f"{aligned_seq1}\n{pattern}\n{aligned_seq2}\n"

    def __str__(self):
        """Return a human-readable string representation of the alignment.

        For sequence alignments, each line has at most 80 columns.
        The first 10 columns show the (possibly truncated) sequence name,
        which may be the id attribute of a SeqRecord, or otherwise 'target'
        or 'query' for pairwise alignments.
        The next 10 columns show the sequence coordinate, using zero-based
        counting as usual in Python.
        The remaining 60 columns shown the sequence, using dashes to represent
        gaps.
        At the end of the alignment, the end coordinates are shown on the right
        of the sequence, again in zero-based coordinates.

        Pairwise alignments have an additional line between the two sequences
        showing whether the sequences match ('|') or mismatch ('.'), or if
        there is a gap ('-').
        The coordinates shown for this line are the column indices, which can
        be useful when extracting a subalignment.

        For example,

        >>> from Bio.Align import PairwiseAligner
        >>> aligner = PairwiseAligner()

        >>> seqA = "TTAACCCCATTTG"
        >>> seqB = "AAGCCCCTTT"
        >>> seqC = "AAAGGGGCTT"

        >>> alignments = aligner.align(seqA, seqB)
        >>> len(alignments)
        1
        >>> alignment = alignments[0]
        >>> print(alignment)
        target            0 TTAA-CCCCATTTG 13
                          0 --||-||||-|||- 14
        query             0 --AAGCCCC-TTT- 10
        <BLANKLINE>

        Note that seqC is the reverse complement of seqB. Aligning it to the
        reverse strand gives the same alignment, but the query coordinates are
        switched:

        >>> alignments = aligner.align(seqA, seqC, strand="-")
        >>> len(alignments)
        1
        >>> alignment = alignments[0]
        >>> print(alignment)
        target            0 TTAA-CCCCATTTG 13
                          0 --||-||||-|||- 14
        query            10 --AAGCCCC-TTT-  0
        <BLANKLINE>

        """
        return self.format()

    def __repr__(self):
        """Return a representation of the alignment, including its shape.

        The representation cannot be used with eval() to recreate the object,
        which is usually possible with simple python objects.  For example:

        <Alignment object (2 rows x 14 columns) at 0x10403d850>

        The hex string is the memory address of the object and can be used to
        distinguish different Alignment objects.  See help(id) for more
        information.

        >>> import numpy
        >>> from Bio.Align import Alignment
        >>> alignment = Alignment(("ACCGT", "ACGT"),
        ...                       coordinates = numpy.array([[0, 2, 3, 5],
        ...                                                  [0, 2, 2, 4],
        ...                                                 ]))
        >>> print(alignment)
        target            0 ACCGT 5
                          0 ||-|| 5
        query             0 AC-GT 4
        <BLANKLINE>
        >>> alignment  # doctest:+ELLIPSIS
        <Alignment object (2 rows x 5 columns) at 0x...>
        """
        if self.coordinates is None:
            return "<%s object at 0x%x>" % (
                self.__class__.__name__,
                id(self),
            )
        n, m = self.shape
        return "<%s object (%i rows x %i columns) at 0x%x>" % (
            self.__class__.__name__,
            n,
            m,
            id(self),
        )

    def __len__(self):
        """Return the number of sequences in the alignment."""
        return len(self.sequences)

    @property
    def shape(self):
        """Return the shape of the alignment as a tuple of two integer values.

        The first integer value is the number of sequences in the alignment as
        returned by len(alignment), which is always 2 for pairwise alignments.

        The second integer value is the number of columns in the alignment when
        it is printed, and is equal to the sum of the number of matches, number
        of mismatches, and the total length of gaps in the target and query.
        Sequence sections beyond the aligned segment are not included in the
        number of columns.

        For example,

        >>> from Bio import Align
        >>> aligner = Align.PairwiseAligner()
        >>> aligner.mode = "global"
        >>> alignments = aligner.align("GACCTG", "CGATCG")
        >>> alignment = alignments[0]
        >>> print(alignment)
        target            0 -GACCT-G 6
                          0 -||--|-| 8
        query             0 CGA--TCG 6
        <BLANKLINE>
        >>> len(alignment)
        2
        >>> alignment.shape
        (2, 8)
        >>> aligner.mode = "local"
        >>> alignments = aligner.align("GACCTG", "CGATCG")
        >>> alignment = alignments[0]
        >>> print(alignment)
        target            0 GACCT-G 6
                          0 ||--|-| 7
        query             1 GA--TCG 6
        <BLANKLINE>
        >>> len(alignment)
        2
        >>> alignment.shape
        (2, 7)
        """
        n = len(self.coordinates)
        if n == 0:  # no sequences
            return (0, 0)
        steps = numpy.diff(self.coordinates, 1)
        aligned = sum(steps != 0, 0) > 1
        # True for steps in which at least two sequences align, False if a gap
        for i in range(n):
            row = steps[i, aligned]
            if (row >= 0).all():
                pass
            elif (row <= 0).all():
                steps[i, :] = -steps[i, :]
            else:
                raise ValueError(f"Inconsistent steps in row {i}")
        gaps = steps.max(0)
        if not ((steps == gaps) | (steps <= 0)).all():
            raise ValueError("Unequal step sizes in alignment")
        m = sum(gaps)
        return (n, m)

    @property
    def aligned(self):
        """Return the indices of subsequences aligned to each other.

        This property returns the start and end indices of subsequences
        in the target and query sequence that were aligned to each other.
        If the alignment between target (t) and query (q) consists of N
        chunks, you get two tuples of length N:

            (((t_start1, t_end1), (t_start2, t_end2), ..., (t_startN, t_endN)),
             ((q_start1, q_end1), (q_start2, q_end2), ..., (q_startN, q_endN)))

        For example,

        >>> from Bio import Align
        >>> aligner = Align.PairwiseAligner()
        >>> alignments = aligner.align("GAACT", "GAT")
        >>> alignment = alignments[0]
        >>> print(alignment)
        target            0 GAACT 5
                          0 ||--| 5
        query             0 GA--T 3
        <BLANKLINE>
        >>> alignment.aligned
        array([[[0, 2],
                [4, 5]],
        <BLANKLINE>
               [[0, 2],
                [2, 3]]])
        >>> alignment = alignments[1]
        >>> print(alignment)
        target            0 GAACT 5
                          0 |-|-| 5
        query             0 G-A-T 3
        <BLANKLINE>
        >>> alignment.aligned
        array([[[0, 1],
                [2, 3],
                [4, 5]],
        <BLANKLINE>
               [[0, 1],
                [1, 2],
                [2, 3]]])

        Note that different alignments may have the same subsequences
        aligned to each other. In particular, this may occur if alignments
        differ from each other in terms of their gap placement only:

        >>> aligner.mismatch_score = -10
        >>> alignments = aligner.align("AAACAAA", "AAAGAAA")
        >>> len(alignments)
        2
        >>> print(alignments[0])
        target            0 AAAC-AAA 7
                          0 |||--||| 8
        query             0 AAA-GAAA 7
        <BLANKLINE>
        >>> alignments[0].aligned
        array([[[0, 3],
                [4, 7]],
        <BLANKLINE>
               [[0, 3],
                [4, 7]]])
        >>> print(alignments[1])
        target            0 AAA-CAAA 7
                          0 |||--||| 8
        query             0 AAAG-AAA 7
        <BLANKLINE>
        >>> alignments[1].aligned
        array([[[0, 3],
                [4, 7]],
        <BLANKLINE>
               [[0, 3],
                [4, 7]]])

        The property can be used to identify alignments that are identical
        to each other in terms of their aligned sequences.
        """
        if len(self.sequences) > 2:
            raise NotImplementedError(
                "aligned is currently implemented for pairwise alignments only"
            )
        coordinates = self.coordinates.copy()
        steps = numpy.diff(coordinates, 1)
        aligned = sum(steps != 0, 0) > 1
        # True for steps in which at least two sequences align, False if a gap
        for i, sequence in enumerate(self.sequences):
            row = steps[i, aligned]
            if (row >= 0).all():
                pass
            elif (row <= 0).all():
                steps[i, :] = -steps[i, :]
                coordinates[i, :] = len(sequence) - coordinates[i, :]
            else:
                raise ValueError(f"Inconsistent steps in row {i}")
        coordinates = coordinates.transpose()
        steps = numpy.diff(coordinates, axis=0)
        steps = abs(steps).min(1)
        indices = numpy.flatnonzero(steps)
        starts = coordinates[indices, :]
        ends = coordinates[indices + 1, :]
        segments = numpy.stack([starts, ends], axis=0).transpose()
        steps = numpy.diff(self.coordinates, 1)
        for i, sequence in enumerate(self.sequences):
            row = steps[i, aligned]
            if (row >= 0).all():
                pass
            elif (row <= 0).all():  # mapped to reverse strand
                segments[i, :] = len(sequence) - segments[i, :]
            else:
                raise ValueError(f"Inconsistent steps in row {i}")
        return segments

    @property
    def indices(self):
        """Return the sequence index of each lettter in the alignment.

        This property returns a 2D numpy array with the sequence index of each
        letter in the alignment. Gaps are indicated by -1.  The array has the
        same number of rows and columns as the alignment, as given by
        `self.shape`.

        For example,

        >>> from Bio import Align
        >>> aligner = Align.PairwiseAligner()
        >>> aligner.mode = "local"

        >>> alignments = aligner.align("GAACTGG", "AATG")
        >>> alignment = alignments[0]
        >>> print(alignment)
        target            1 AACTG 6
                          0 ||-|| 5
        query             0 AA-TG 4
        <BLANKLINE>
        >>> alignment.indices
        array([[ 1,  2,  3,  4,  5],
               [ 0,  1, -1,  2,  3]])
        >>> alignment = alignments[1]
        >>> print(alignment)
        target            1 AACTGG 7
                          0 ||-|-| 6
        query             0 AA-T-G 4
        <BLANKLINE>
        >>> alignment.indices
        array([[ 1,  2,  3,  4,  5,  6],
               [ 0,  1, -1,  2, -1,  3]])

        >>> alignments = aligner.align("GAACTGG", "CATT", strand="-")
        >>> alignment = alignments[0]
        >>> print(alignment)
        target            1 AACTG 6
                          0 ||-|| 5
        query             4 AA-TG 0
        <BLANKLINE>
        >>> alignment.indices
        array([[ 1,  2,  3,  4,  5],
               [ 3,  2, -1,  1,  0]])
        >>> alignment = alignments[1]
        >>> print(alignment)
        target            1 AACTGG 7
                          0 ||-|-| 6
        query             4 AA-T-G 0
        <BLANKLINE>
        >>> alignment.indices
        array([[ 1,  2,  3,  4,  5,  6],
               [ 3,  2, -1,  1, -1,  0]])

        """
        a = -numpy.ones(self.shape, int)
        n, m = self.coordinates.shape
        steps = numpy.diff(self.coordinates, 1)
        aligned = sum(steps != 0, 0) > 1
        # True for steps in which at least two sequences align, False if a gap
        steps = steps[:, aligned]
        rcs = numpy.zeros(n, bool)
        for i, row in enumerate(steps):
            if (row >= 0).all():
                rcs[i] = False
            elif (row <= 0).all():
                rcs[i] = True
            else:
                raise ValueError(f"Inconsistent steps in row {i}")
        i = 0
        j = 0
        ends = self.coordinates[:, 0]
        for k in range(1, m):
            starts = ends
            ends = self.coordinates[:, k]
            for row, start, end, rc in zip(a, starts, ends, rcs):
                if rc == False and start < end:  # noqa: 712
                    j = i + end - start
                    row[i:j] = range(start, end)
                elif rc == True and start > end:  # noqa: 712
                    j = i + start - end
                    row[i:j] = range(start - 1, end - 1, -1)
            i = j
        return a

    @property
    def inverse_indices(self):
        """Return the alignment column index for each letter in each sequence.

        This property returns a list of 1D numpy arrays; the number of arrays
        is equal to the number of aligned sequences, and the length of each
        array is equal to the length of the corresponding sequence. For each
        letter in each sequence, the array contains the corresponding column
        index in the alignment. Letters not included in the alignment are
        indicated by -1.

        For example,

        >>> from Bio import Align
        >>> aligner = Align.PairwiseAligner()
        >>> aligner.mode = "local"

        >>> alignments = aligner.align("GAACTGG", "AATG")
        >>> alignment = alignments[0]
        >>> print(alignment)
        target            1 AACTG 6
                          0 ||-|| 5
        query             0 AA-TG 4
        <BLANKLINE>
        >>> alignment.inverse_indices
        [array([-1,  0,  1,  2,  3,  4, -1]), array([0, 1, 3, 4])]
        >>> alignment = alignments[1]
        >>> print(alignment)
        target            1 AACTGG 7
                          0 ||-|-| 6
        query             0 AA-T-G 4
        <BLANKLINE>
        >>> alignment.inverse_indices
        [array([-1,  0,  1,  2,  3,  4,  5]), array([0, 1, 3, 5])]
        >>> alignments = aligner.align("GAACTGG", "CATT", strand="-")
        >>> alignment = alignments[0]
        >>> print(alignment)
        target            1 AACTG 6
                          0 ||-|| 5
        query             4 AA-TG 0
        <BLANKLINE>
        >>> alignment.inverse_indices
        [array([-1,  0,  1,  2,  3,  4, -1]), array([4, 3, 1, 0])]
        >>> alignment = alignments[1]
        >>> print(alignment)
        target            1 AACTGG 7
                          0 ||-|-| 6
        query             4 AA-T-G 0
        <BLANKLINE>
        >>> alignment.inverse_indices
        [array([-1,  0,  1,  2,  3,  4,  5]), array([5, 3, 1, 0])]

        """
        a = [-numpy.ones(len(sequence), int) for sequence in self.sequences]
        n, m = self.coordinates.shape
        steps = numpy.diff(self.coordinates, 1)
        aligned = sum(steps != 0, 0) > 1
        # True for steps in which at least two sequences align, False if a gap
        steps = steps[:, aligned]
        rcs = numpy.zeros(n, bool)
        for i, row in enumerate(steps):
            if (row >= 0).all():
                rcs[i] = False
            elif (row <= 0).all():
                rcs[i] = True
            else:
                raise ValueError(f"Inconsistent steps in row {i}")
        i = 0
        j = 0
        for k in range(m - 1):
            starts = self.coordinates[:, k]
            ends = self.coordinates[:, k + 1]
            for row, start, end, rc in zip(a, starts, ends, rcs):
                if rc == False and start < end:  # noqa: 712
                    j = i + end - start
                    row[start:end] = range(i, j)
                elif rc == True and start > end:  # noqa: 712
                    j = i + start - end
                    if end > 0:
                        row[start - 1 : end - 1 : -1] = range(i, j)
                    elif start > 0:
                        row[start - 1 :: -1] = range(i, j)
            i = j
        return a

    def sort(self, key=None, reverse=False):
        """Sort the sequences of the alignment in place.

        By default, this sorts the sequences alphabetically using their id
        attribute if available, or by their sequence contents otherwise.
        For example,

        >>> from Bio.Align import PairwiseAligner
        >>> aligner = PairwiseAligner()
        >>> aligner.gap_score = -1
        >>> alignments = aligner.align("AATAA", "AAGAA")
        >>> len(alignments)
        1
        >>> alignment = alignments[0]
        >>> print(alignment)
        target            0 AATAA 5
                          0 ||.|| 5
        query             0 AAGAA 5
        <BLANKLINE>
        >>> alignment.sort()
        >>> print(alignment)
        target            0 AAGAA 5
                          0 ||.|| 5
        query             0 AATAA 5
        <BLANKLINE>

        Alternatively, a key function can be supplied that maps each sequence
        to a sort value.  For example, you could sort on the GC content of each
        sequence.

        >>> from Bio.SeqUtils import gc_fraction
        >>> alignment.sort(key=gc_fraction)
        >>> print(alignment)
        target            0 AATAA 5
                          0 ||.|| 5
        query             0 AAGAA 5
        <BLANKLINE>

        You can reverse the sort order by passing `reverse=True`:

        >>> alignment.sort(key=gc_fraction, reverse=True)
        >>> print(alignment)
        target            0 AAGAA 5
                          0 ||.|| 5
        query             0 AATAA 5
        <BLANKLINE>

        The sequences are now sorted by decreasing GC content value.
        """
        sequences = self.sequences
        if key is None:
            try:
                values = [sequence.id for sequence in sequences]
            except AttributeError:
                values = sequences
        else:
            values = [key(sequence) for sequence in sequences]
        indices = sorted(range(len(sequences)), key=values.__getitem__, reverse=reverse)
        self.sequences = [sequences[index] for index in indices]
        self.coordinates = self.coordinates.take(indices, 0)

    def map(self, alignment):
        r"""Map the alignment to self.target and return the resulting alignment.

        Here, self.query and alignment.target are the same sequence.

        A typical example is where self is the pairwise alignment between a
        chromosome and a transcript, the argument is the pairwise alignment
        between the transcript and a sequence (e.g., as obtained by RNA-seq),
        and we want to find the alignment of the sequence to the chromosome:

        >>> from Bio import Align
        >>> aligner = Align.PairwiseAligner()
        >>> aligner.mode = 'local'
        >>> aligner.open_gap_score = -1
        >>> aligner.extend_gap_score = 0
        >>> chromosome = "AAAAAAAACCCCCCCAAAAAAAAAAAGGGGGGAAAAAAAA"
        >>> transcript = "CCCCCCCGGGGGG"
        >>> alignments1 = aligner.align(chromosome, transcript)
        >>> len(alignments1)
        1
        >>> alignment1 = alignments1[0]
        >>> print(alignment1)
        target            8 CCCCCCCAAAAAAAAAAAGGGGGG 32
                          0 |||||||-----------|||||| 24
        query             0 CCCCCCC-----------GGGGGG 13
        <BLANKLINE>
        >>> sequence = "CCCCGGGG"
        >>> alignments2 = aligner.align(transcript, sequence)
        >>> len(alignments2)
        1
        >>> alignment2 = alignments2[0]
        >>> print(alignment2)
        target            3 CCCCGGGG 11
                          0 ||||||||  8
        query             0 CCCCGGGG  8
        <BLANKLINE>
        >>> alignment = alignment1.map(alignment2)
        >>> print(alignment)
        target           11 CCCCAAAAAAAAAAAGGGG 30
                          0 ||||-----------|||| 19
        query             0 CCCC-----------GGGG  8
        <BLANKLINE>
        >>> format(alignment, "psl")
        '8\t0\t0\t0\t0\t0\t1\t11\t+\tquery\t8\t0\t8\ttarget\t40\t11\t30\t2\t4,4,\t0,4,\t11,26,\n'

        Mapping the alignment does not depend on the sequence contents. If we
        delete the sequence contents, the same alignment is found in PSL format
        (though we obviously lose the ability to print the sequence alignment):

        >>> alignment1.target = Seq(None, len(alignment1.target))
        >>> alignment1.query = Seq(None, len(alignment1.query))
        >>> alignment2.target = Seq(None, len(alignment2.target))
        >>> alignment2.query = Seq(None, len(alignment2.query))
        >>> alignment = alignment1.map(alignment2)
        >>> format(alignment, "psl")
        '8\t0\t0\t0\t0\t0\t1\t11\t+\tquery\t8\t0\t8\ttarget\t40\t11\t30\t2\t4,4,\t0,4,\t11,26,\n'
        """
        alignment1, alignment2 = self, alignment
        if len(alignment1.query) != len(alignment2.target):
            raise ValueError(
                "length of alignment1 query sequence (%d) != length of alignment2 target sequence (%d)"
                % (len(alignment1.query), len(alignment2.target))
            )
        target = alignment1.target
        query = alignment2.query
        coordinates1 = alignment1.coordinates
        coordinates2 = alignment2.coordinates
        n1 = len(alignment1.query)
        n2 = len(alignment2.query)
        steps1 = numpy.diff(coordinates1, 1)
        row = numpy.prod(numpy.sign(steps1), 0)
        if (row >= 0).all():
            strand1 = "+"
        elif (row <= 0).all():
            strand1 = "-"
        else:
            raise ValueError("Inconsistent steps in the first alignment")
        steps2 = numpy.diff(coordinates2, 1)
        row = numpy.prod(numpy.sign(steps2), 0)
        if (row >= 0).all():
            strand2 = "+"
        elif (row <= 0).all():
            strand2 = "-"
        else:
            raise ValueError("Inconsistent steps in the second alignment")
        if strand1 == "+":
            if strand2 == "-":  # mapped to reverse strand
                coordinates2 = coordinates2.copy()
                coordinates2[1, :] = n2 - coordinates2[1, :]
        else:  # mapped to reverse strand
            coordinates1 = coordinates1.copy()
            coordinates1[1, :] = n1 - coordinates1[1, :]
            coordinates2 = coordinates2.copy()
            coordinates2[0, :] = n1 - coordinates2[0, ::-1]
            if strand2 == "+":
                coordinates2[1, :] = n2 - coordinates2[1, ::-1]
            else:  # mapped to reverse strand
                coordinates2[1, :] = coordinates2[1, ::-1]
        steps1 = numpy.diff(coordinates1, 1)
        gaps1 = steps1.max(0)
        if not ((steps1 == gaps1) | (steps1 <= 0)).all():
            raise ValueError("Unequal step sizes in first alignment")
        steps2 = numpy.diff(coordinates2, 1)
        gaps2 = steps2.max(0)
        if not ((steps2 == gaps2) | (steps2 <= 0)).all():
            raise ValueError("Unequal step sizes in second alignment")
        path = []
        tEnd, qEnd = sys.maxsize, sys.maxsize
        coordinates1 = iter(coordinates1.transpose())
        tStart1, qStart1 = sys.maxsize, sys.maxsize
        for tEnd1, qEnd1 in coordinates1:
            if tStart1 < tEnd1 and qStart1 < qEnd1:
                break
            tStart1, qStart1 = tEnd1, qEnd1
        tStart2, qStart2 = sys.maxsize, sys.maxsize
        for tEnd2, qEnd2 in coordinates2.transpose():
            while qStart2 < qEnd2 and tStart2 < tEnd2:
                while True:
                    if tStart2 < qStart1:
                        if tEnd2 < qStart1:
                            size = tEnd2 - tStart2
                        else:
                            size = qStart1 - tStart2
                        break
                    elif tStart2 < qEnd1:
                        offset = tStart2 - qStart1
                        if tEnd2 > qEnd1:
                            size = qEnd1 - tStart2
                        else:
                            size = tEnd2 - tStart2
                        qStart = qStart2
                        tStart = tStart1 + offset
                        if tStart > tEnd and qStart > qEnd:
                            # adding a gap both in target and in query;
                            # add gap to target first:
                            path.append([tStart, qEnd])
                        qEnd = qStart2 + size
                        tEnd = tStart + size
                        path.append([tStart, qStart])
                        path.append([tEnd, qEnd])
                        break
                    tStart1, qStart1 = sys.maxsize, sys.maxsize
                    for tEnd1, qEnd1 in coordinates1:
                        if tStart1 < tEnd1 and qStart1 < qEnd1:
                            break
                        tStart1, qStart1 = tEnd1, qEnd1
                    else:
                        size = qEnd2 - qStart2
                        break
                qStart2 += size
                tStart2 += size
            tStart2, qStart2 = tEnd2, qEnd2
        coordinates = numpy.array(path).transpose()
        if strand1 != strand2:
            coordinates[1, :] = n2 - coordinates[1, :]
        sequences = [target, query]
        alignment = Alignment(sequences, coordinates)
        return alignment

    @property
    def substitutions(self):
        """Return an Array with the number of substitutions of letters in the alignment.

        As an example, consider a sequence alignment of two RNA sequences:

        >>> from Bio.Align import PairwiseAligner
        >>> target = "ATACTTACCTGGCAGGGGAGATACCATGATCACGAAGGTGGTTTTCCCAGGGCGAGGCTTATCCATTGCACTCCGGATGTGCTGACCCCTGCGATTTCCCCAAATGTGGGAAACTCGACTGCATAATTTGTGGTAGTGGGGGACTGCGTTCGCGCTTTCCCCTG"  # human spliceosomal small nuclear RNA U1
        >>> query = "ATACTTACCTGACAGGGGAGGCACCATGATCACACAGGTGGTCCTCCCAGGGCGAGGCTCTTCCATTGCACTGCGGGAGGGTTGACCCCTGCGATTTCCCCAAATGTGGGAAACTCGACTGTATAATTTGTGGTAGTGGGGGACTGCGTTCGCGCTATCCCCCG"  # sea lamprey spliceosomal small RNA U1
        >>> aligner = PairwiseAligner()
        >>> aligner.gap_score = -10
        >>> alignments = aligner.align(target, query)
        >>> len(alignments)
        1
        >>> alignment = alignments[0]
        >>> print(alignment)
        target            0 ATACTTACCTGGCAGGGGAGATACCATGATCACGAAGGTGGTTTTCCCAGGGCGAGGCTT
                          0 |||||||||||.||||||||..|||||||||||..|||||||..|||||||||||||||.
        query             0 ATACTTACCTGACAGGGGAGGCACCATGATCACACAGGTGGTCCTCCCAGGGCGAGGCTC
        <BLANKLINE>
        target           60 ATCCATTGCACTCCGGATGTGCTGACCCCTGCGATTTCCCCAAATGTGGGAAACTCGACT
                         60 .|||||||||||.|||..|.|.||||||||||||||||||||||||||||||||||||||
        query            60 TTCCATTGCACTGCGGGAGGGTTGACCCCTGCGATTTCCCCAAATGTGGGAAACTCGACT
        <BLANKLINE>
        target          120 GCATAATTTGTGGTAGTGGGGGACTGCGTTCGCGCTTTCCCCTG 164
                        120 |.||||||||||||||||||||||||||||||||||.|||||.| 164
        query           120 GTATAATTTGTGGTAGTGGGGGACTGCGTTCGCGCTATCCCCCG 164
        <BLANKLINE>
        >>> m = alignment.substitutions
        >>> print(m)
             A    C    G    T
        A 28.0  1.0  2.0  1.0
        C  0.0 39.0  1.0  2.0
        G  2.0  0.0 45.0  0.0
        T  2.0  5.0  1.0 35.0
        <BLANKLINE>

        Note that the matrix is not symmetric: rows correspond to the target
        sequence, and columns to the query sequence.  For example, the number
        of T's in the target sequence that are aligned to a C in the query
        sequence is

        >>> m['T', 'C']
        5.0

        and the number of C's in the query sequence tat are aligned to a T in
        the query sequence is

        >>> m['C', 'T']
        2.0

        For some applications (for example, to define a scoring matrix from
        the substitution matrix), a symmetric matrix may be preferred, which
        can be calculated as follows:

        >>> m += m.transpose()
        >>> m /= 2.0
        >>> print(m)
             A    C    G    T
        A 28.0  0.5  2.0  1.5
        C  0.5 39.0  0.5  3.5
        G  2.0  0.5 45.0  0.5
        T  1.5  3.5  0.5 35.0
        <BLANKLINE>

        The matrix is now symmetric, with counts divided equally on both sides
        of the diagonal:

        >>> m['C', 'T']
        3.5
        >>> m['T', 'C']
        3.5

        The total number of substitutions between T's and C's in the alignment
        is 3.5 + 3.5 = 7.
        """
        coordinates = self.coordinates.copy()
        sequences = list(self.sequences)
        steps = numpy.diff(self.coordinates, 1)
        aligned = sum(steps != 0, 0) > 1
        # True for steps in which at least two sequences align, False if a gap
        for i, sequence in enumerate(sequences):
            row = steps[i, aligned]
            if (row >= 0).all():
                pass
            elif (row <= 0).all():
                sequences[i] = reverse_complement(sequence, inplace=False)
                coordinates[i, :] = len(sequence) - coordinates[i, :]
            else:
                raise ValueError(f"Inconsistent steps in row {i}")
        letters = set()
        for sequence in sequences:
            try:
                s = set(sequence)
            except UndefinedSequenceError:
                try:
                    sequence = sequence.seq  # SeqRecord confusion
                except AttributeError:
                    pass
                for start, end in sequence.defined_ranges:
                    s = set(sequence[start:end])
                    letters.update(s)
            else:
                letters.update(s)
        letters = "".join(sorted(letters))
        m = substitution_matrices.Array(letters, dims=2)
        n = len(sequences)
        for i1 in range(n):
            sequence1 = sequences[i1]
            coordinates1 = coordinates[i1, :]
            for i2 in range(i1 + 1, n):
                sequence2 = sequences[i2]
                coordinates2 = coordinates[i2, :]
                start1, start2 = sys.maxsize, sys.maxsize
                for end1, end2 in zip(coordinates1, coordinates2):
                    if start1 < end1 and start2 < end2:  # aligned
                        segment1 = sequence1[start1:end1]
                        segment2 = sequence2[start2:end2]
                        if len(segment1) != len(segment2):
                            raise ValueError("Unequal step sizes in alignment")
                        for c1, c2 in zip(segment1, segment2):
                            m[c1, c2] += 1.0
                    start1, start2 = end1, end2
        return m

    def counts(self):
        """Return number of identities, mismatches, and gaps, of a pairwise alignment.

        >>> aligner = PairwiseAligner(mode='global', match_score=2, mismatch_score=-1)
        >>> for alignment in aligner.align("TACCG", "ACG"):
        ...     print("Score = %.1f:" % alignment.score)
        ...     c = alignment.counts()  # namedtuple
        ...     print(f"{c.gaps} gaps, {c.identities} identities, {c.mismatches} mismatches")
        ...     print(alignment)
        ...
        Score = 6.0:
        2 gaps, 3 identities, 0 mismatches
        target            0 TACCG 5
                          0 -||-| 5
        query             0 -AC-G 3
        <BLANKLINE>
        Score = 6.0:
        2 gaps, 3 identities, 0 mismatches
        target            0 TACCG 5
                          0 -|-|| 5
        query             0 -A-CG 3
        <BLANKLINE>

        This classifies each pair of letters in a pairwise alignment into gaps,
        perfect matches, or mismatches. It has been defined as a method (not a
        property) so that it may in future take optional argument(s) allowing
        the behaviour to be customised. These three values are returned as a
        namedtuple. This is calculated for all the pairs of sequences in the
        alignment.
        """
        gaps = identities = mismatches = 0
        for i, seq1 in enumerate(self):
            for j, seq2 in enumerate(self):
                if i == j:
                    # Don't count seq1 vs seq2 and seq2 vs seq1
                    break
                for a, b in zip(seq1, seq2):
                    if a == "-" or b == "-":
                        gaps += 1
                    elif a == b:
                        identities += 1
                    else:
                        mismatches += 1
        return AlignmentCounts(gaps, identities, mismatches)


class PairwiseAlignments:
    """Implements an iterator over pairwise alignments returned by the aligner.

    This class also supports indexing, which is fast for increasing indices,
    but may be slow for random access of a large number of alignments.

    Note that pairwise aligners can return an astronomical number of alignments,
    even for relatively short sequences, if they align poorly to each other. We
    therefore recommend to first check the number of alignments, accessible as
    len(alignments), which can be calculated quickly even if the number of
    alignments is very large.
    """

    def __init__(self, seqA, seqB, score, paths):
        """Initialize a new PairwiseAlignments object.

        Arguments:
         - seqA  - The first sequence, as a plain string, without gaps.
         - seqB  - The second sequence, as a plain string, without gaps.
         - score - The alignment score.
         - paths - An iterator over the paths in the traceback matrix;
                   each path defines one alignment.

        You would normally obtain a PairwiseAlignments object by calling
        aligner.align(seqA, seqB), where aligner is a PairwiseAligner object.
        """
        self.sequences = [seqA, seqB]
        self.score = score
        self._paths = paths
        self._index = -1

    def __len__(self):
        """Return the number of alignments."""
        return len(self._paths)

    def __getitem__(self, index):
        if not isinstance(index, int):
            raise TypeError(f"index must be an integer, not {index.__class__.__name__}")
        if index < 0:
            index += len(self._paths)
        if index == self._index:
            return self._alignment
        if index < self._index:
            self._paths.reset()
            self._index = -1
        while True:
            try:
                alignment = next(self)
            except StopIteration:
                raise IndexError("index out of range") from None
            if self._index == index:
                break
        return alignment

    def __iter__(self):
        self._paths.reset()
        self._index = -1
        return self

    def __next__(self):
        path = next(self._paths)
        self._index += 1
        coordinates = numpy.array(path)
        alignment = Alignment(self.sequences, coordinates)
        alignment.score = self.score
        self._alignment = alignment
        return alignment


class PairwiseAligner(_aligners.PairwiseAligner):
    """Performs pairwise sequence alignment using dynamic programming.

    This provides functions to get global and local alignments between two
    sequences.  A global alignment finds the best concordance between all
    characters in two sequences.  A local alignment finds just the
    subsequences that align the best.

    To perform a pairwise sequence alignment, first create a PairwiseAligner
    object.  This object stores the match and mismatch scores, as well as the
    gap scores.  Typically, match scores are positive, while mismatch scores
    and gap scores are negative or zero.  By default, the match score is 1,
    and the mismatch and gap scores are zero.  Based on the values of the gap
    scores, a PairwiseAligner object automatically chooses the appropriate
    alignment algorithm (the Needleman-Wunsch, Smith-Waterman, Gotoh, or
    Waterman-Smith-Beyer global or local alignment algorithm).

    Calling the "score" method on the aligner with two sequences as arguments
    will calculate the alignment score between the two sequences.
    Calling the "align" method on the aligner with two sequences as arguments
    will return a generator yielding the alignments between the two
    sequences.

    Some examples:

    >>> from Bio import Align
    >>> aligner = Align.PairwiseAligner()
    >>> alignments = aligner.align("TACCG", "ACG")
    >>> for alignment in sorted(alignments):
    ...     print("Score = %.1f:" % alignment.score)
    ...     print(alignment)
    ...
    Score = 3.0:
    target            0 TACCG 5
                      0 -|-|| 5
    query             0 -A-CG 3
    <BLANKLINE>
    Score = 3.0:
    target            0 TACCG 5
                      0 -||-| 5
    query             0 -AC-G 3
    <BLANKLINE>

    Specify the aligner mode as local to generate local alignments:

    >>> aligner.mode = 'local'
    >>> alignments = aligner.align("TACCG", "ACG")
    >>> for alignment in sorted(alignments):
    ...     print("Score = %.1f:" % alignment.score)
    ...     print(alignment)
    ...
    Score = 3.0:
    target            1 ACCG 5
                      0 |-|| 4
    query             0 A-CG 3
    <BLANKLINE>
    Score = 3.0:
    target            1 ACCG 5
                      0 ||-| 4
    query             0 AC-G 3
    <BLANKLINE>

    Do a global alignment.  Identical characters are given 2 points,
    1 point is deducted for each non-identical character.

    >>> aligner.mode = 'global'
    >>> aligner.match_score = 2
    >>> aligner.mismatch_score = -1
    >>> for alignment in aligner.align("TACCG", "ACG"):
    ...     print("Score = %.1f:" % alignment.score)
    ...     print(alignment)
    ...
    Score = 6.0:
    target            0 TACCG 5
                      0 -||-| 5
    query             0 -AC-G 3
    <BLANKLINE>
    Score = 6.0:
    target            0 TACCG 5
                      0 -|-|| 5
    query             0 -A-CG 3
    <BLANKLINE>

    Same as above, except now 0.5 points are deducted when opening a
    gap, and 0.1 points are deducted when extending it.

    >>> aligner.open_gap_score = -0.5
    >>> aligner.extend_gap_score = -0.1
    >>> aligner.target_end_gap_score = 0.0
    >>> aligner.query_end_gap_score = 0.0
    >>> for alignment in aligner.align("TACCG", "ACG"):
    ...     print("Score = %.1f:" % alignment.score)
    ...     print(alignment)
    ...
    Score = 5.5:
    target            0 TACCG 5
                      0 -|-|| 5
    query             0 -A-CG 3
    <BLANKLINE>
    Score = 5.5:
    target            0 TACCG 5
                      0 -||-| 5
    query             0 -AC-G 3
    <BLANKLINE>

    The alignment function can also use known matrices already included in
    Biopython:

    >>> from Bio.Align import substitution_matrices
    >>> aligner = Align.PairwiseAligner()
    >>> aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    >>> alignments = aligner.align("KEVLA", "EVL")
    >>> alignments = list(alignments)
    >>> print("Number of alignments: %d" % len(alignments))
    Number of alignments: 1
    >>> alignment = alignments[0]
    >>> print("Score = %.1f" % alignment.score)
    Score = 13.0
    >>> print(alignment)
    target            0 KEVLA 5
                      0 -|||- 5
    query             0 -EVL- 3
    <BLANKLINE>

    You can also set the value of attributes directly during construction
    of the PairwiseAligner object by providing them as keyword arguments:

    >>> aligner = Align.PairwiseAligner(mode='global', match_score=2, mismatch_score=-1)
    >>> for alignment in aligner.align("TACCG", "ACG"):
    ...     print("Score = %.1f:" % alignment.score)
    ...     print(alignment)
    ...
    Score = 6.0:
    target            0 TACCG 5
                      0 -||-| 5
    query             0 -AC-G 3
    <BLANKLINE>
    Score = 6.0:
    target            0 TACCG 5
                      0 -|-|| 5
    query             0 -A-CG 3
    <BLANKLINE>

    """

    def __init__(self, scoring=None, **kwargs):
        """Initialize a new PairwiseAligner as specified by the keyword arguments.

        If scoring is None, use the default scoring scheme match = 1.0,
        mismatch = 0.0, gap score = 0.0
        If scoring is "blastn", "megablast", or "blastp", use the default
        substitution matrix and gap scores for BLASTN, MEGABLAST, or BLASTP,
        respectively.

        Loops over the remaining keyword arguments and sets them as attributes
        on the object.
        """
        super().__init__()
        if scoring is None:
            # use default values:
            # match = 1.0
            # mismatch = 0.0
            # gap_score = 0.0
            pass
        elif scoring == "blastn":
            self.substitution_matrix = substitution_matrices.load("BLASTN")
            self.open_gap_score = -7.0
            self.extend_gap_score = -2.0
        elif scoring == "megablast":
            self.substitution_matrix = substitution_matrices.load("MEGABLAST")
            self.open_gap_score = -2.5
            self.extend_gap_score = -2.5
        elif scoring == "blastp":
            self.substitution_matrix = substitution_matrices.load("BLASTP")
            self.open_gap_score = -12.0
            self.extend_gap_score = -1.0
        else:
            raise ValueError("Unknown scoring scheme '%s'" % scoring)
        for name, value in kwargs.items():
            setattr(self, name, value)

    def __setattr__(self, key, value):
        if key not in dir(_aligners.PairwiseAligner):
            # To prevent confusion, don't allow users to create new attributes.
            # On CPython, __slots__ can be used for this, but currently
            # __slots__ does not behave the same way on PyPy at least.
            raise AttributeError("'PairwiseAligner' object has no attribute '%s'" % key)
        _aligners.PairwiseAligner.__setattr__(self, key, value)

    def align(self, seqA, seqB, strand="+"):
        """Return the alignments of two sequences using PairwiseAligner."""
        if isinstance(seqA, (Seq, MutableSeq, SeqRecord)):
            sA = bytes(seqA)
        else:
            sA = seqA
        if strand == "+":
            sB = seqB
        else:  # strand == "-":
            sB = reverse_complement(seqB, inplace=False)
        if isinstance(seqB, (Seq, MutableSeq, SeqRecord)):
            sB = bytes(sB)
        score, paths = _aligners.PairwiseAligner.align(self, sA, sB, strand)
        alignments = PairwiseAlignments(seqA, seqB, score, paths)
        return alignments

    def score(self, seqA, seqB, strand="+"):
        """Return the alignments score of two sequences using PairwiseAligner."""
        if isinstance(seqA, (Seq, MutableSeq, SeqRecord)):
            seqA = bytes(seqA)
        if strand == "-":
            seqB = reverse_complement(seqB, inplace=False)
        if isinstance(seqB, (Seq, MutableSeq, SeqRecord)):
            seqB = bytes(seqB)
        return _aligners.PairwiseAligner.score(self, seqA, seqB, strand)

    def __getstate__(self):
        state = {
            "wildcard": self.wildcard,
            "target_internal_open_gap_score": self.target_internal_open_gap_score,
            "target_internal_extend_gap_score": self.target_internal_extend_gap_score,
            "target_left_open_gap_score": self.target_left_open_gap_score,
            "target_left_extend_gap_score": self.target_left_extend_gap_score,
            "target_right_open_gap_score": self.target_right_open_gap_score,
            "target_right_extend_gap_score": self.target_right_extend_gap_score,
            "query_internal_open_gap_score": self.query_internal_open_gap_score,
            "query_internal_extend_gap_score": self.query_internal_extend_gap_score,
            "query_left_open_gap_score": self.query_left_open_gap_score,
            "query_left_extend_gap_score": self.query_left_extend_gap_score,
            "query_right_open_gap_score": self.query_right_open_gap_score,
            "query_right_extend_gap_score": self.query_right_extend_gap_score,
            "mode": self.mode,
        }
        if self.substitution_matrix is None:
            state["match_score"] = self.match_score
            state["mismatch_score"] = self.mismatch_score
        else:
            state["substitution_matrix"] = self.substitution_matrix
        return state

    def __setstate__(self, state):
        self.wildcard = state["wildcard"]
        self.target_internal_open_gap_score = state["target_internal_open_gap_score"]
        self.target_internal_extend_gap_score = state[
            "target_internal_extend_gap_score"
        ]
        self.target_left_open_gap_score = state["target_left_open_gap_score"]
        self.target_left_extend_gap_score = state["target_left_extend_gap_score"]
        self.target_right_open_gap_score = state["target_right_open_gap_score"]
        self.target_right_extend_gap_score = state["target_right_extend_gap_score"]
        self.query_internal_open_gap_score = state["query_internal_open_gap_score"]
        self.query_internal_extend_gap_score = state["query_internal_extend_gap_score"]
        self.query_left_open_gap_score = state["query_left_open_gap_score"]
        self.query_left_extend_gap_score = state["query_left_extend_gap_score"]
        self.query_right_open_gap_score = state["query_right_open_gap_score"]
        self.query_right_extend_gap_score = state["query_right_extend_gap_score"]
        self.mode = state["mode"]
        substitution_matrix = state.get("substitution_matrix")
        if substitution_matrix is None:
            self.match_score = state["match_score"]
            self.mismatch_score = state["mismatch_score"]
        else:
            self.substitution_matrix = substitution_matrix


class PairwiseAlignment(Alignment):
    """Represents a pairwise sequence alignment.

    Internally, the pairwise alignment is stored as the path through
    the traceback matrix, i.e. a tuple of pairs of indices corresponding
    to the vertices of the path in the traceback matrix.
    """

    def __init__(self, target, query, path, score):
        """Initialize a new PairwiseAlignment object.

        Arguments:
         - target  - The first sequence, as a plain string, without gaps.
         - query   - The second sequence, as a plain string, without gaps.
         - path    - The path through the traceback matrix, defining an
                     alignment.
         - score   - The alignment score.

        You would normally obtain a PairwiseAlignment object by iterating
        over a PairwiseAlignments object.
        """
        warnings.warn(
            "The PairwiseAlignment class is deprecated; please use the "
            "Alignment class instead.  Note that the coordinates attribute of "
            "an Alignment object is a numpy array and the transpose of the "
            "path attribute of a PairwiseAlignment object.",
            BiopythonDeprecationWarning,
        )
        sequences = [target, query]
        coordinates = numpy.array(path).transpose()
        super().__init__(sequences, coordinates)
        self.score = score


# fmt: off
formats = (
    "a2m",        # A2M files created by align2model or hmmscore
    "bed",        # BED (Browser Extensible Data) files
    "bigbed",     # bigBed format
    "bigmaf",     # MAF file saved as a bigBed file
    "bigpsl",     # PSL file saved as a bigBed file
    "clustal",    # clustal output from CLUSTAL W and other tools.
    "emboss",     # emboss output from EMBOSS tools such as needle, water
    "exonerate",  # Exonerate pairwise alignment output
    "fasta",      # FASTA format with gaps represented by dashes
    "hhr",        # hhr files generated by HHsearch, HHblits in HH-suite
    "maf",        # MAF (Multiple Alignment Format) format.
    "mauve",      # xmfa output from Mauve/ProgressiveMauve
    "msf",        # MSF format produced by GCG PileUp and LocalPileUp
    "nexus",      # Nexus file format
    "phylip",     # Alignment format for input files for PHYLIP tools
    "psl",        # Pattern Space Layout (PSL) format generated by Blat
    "sam",        # Sequence Alignment/Map (SAM) format
    "stockholm",  # Stockholm file format used by PFAM and others
    "tabular",    # Tabular output from BLAST or FASTA
)
# fmt: on

_modules = {}


def _load(fmt):
    fmt = fmt.lower()
    try:
        return _modules[fmt]
    except KeyError:
        pass
    if fmt not in formats:
        raise ValueError("Unknown file format %s" % fmt)
    module = importlib.import_module(f"Bio.Align.{fmt}")
    _modules[fmt] = module
    return module


def write(alignments, target, fmt, *args, **kwargs):
    """Write alignments to a file.

    Arguments:
     - alignments - List (or iterator) of Alignment objects, or a single
       Alignment.
     - target     - File or file-like object to write to, or filename as string.
     - fmt        - String describing the file format (case-insensitive).

    Note if providing a file or file-like object, your code should close the
    target after calling this function, or call .flush(), to ensure the data
    gets flushed to disk.

    Returns the number of alignments written (as an integer).
    """
    if isinstance(alignments, Alignment):
        alignments = [alignments]

    module = _load(fmt)
    try:
        writer = module.AlignmentWriter
    except AttributeError:
        raise ValueError(
            f"File writing has not yet been implemented for the {fmt} format"
        )
    return writer(target, *args, **kwargs).write_file(alignments)


def parse(source, fmt):
    """Parse an alignment file and return an iterator over alignments.

    Arguments:
     - source - File or file-like object to read from, or filename as string.
     - fmt    - String describing the file format (case-insensitive).

    Typical usage, opening a file to read in, and looping over the aligments:

    >>> from Bio import Align
    >>> filename = "Exonerate/exn_22_m_ner_cigar.exn"
    >>> for alignment in Align.parse(filename, "exonerate"):
    ...    print("Number of sequences in alignment", len(alignment))
    ...    print("Alignment score:", alignment.score)
    Number of sequences in alignment 2
    Alignment score: 6150.0
    Number of sequences in alignment 2
    Alignment score: 502.0
    Number of sequences in alignment 2
    Alignment score: 440.0

    For lazy-loading file formats such as bigMaf, for which the file contents
    is read on demand only, ensure that the file remains open while extracting
    alignment data.

    You can use the Bio.Align.read(...) function when the file contains only
    one alignment.
    """
    module = _load(fmt)
    alignments = module.AlignmentIterator(source)
    return alignments


def read(handle, fmt):
    """Parse a file containing one alignment, and return it.

    Arguments:
     - source - File or file-like object to read from, or filename as string.
     - fmt    - String describing the file format (case-insensitive).

    This function is for use parsing alignment files containing exactly one
    alignment.  For example, reading a Clustal file:

    >>> from Bio import Align
    >>> alignment = Align.read("Clustalw/opuntia.aln", "clustal")
    >>> print("Alignment shape:", alignment.shape)
    Alignment shape: (7, 156)
    >>> for sequence in alignment.sequences:
    ...     print(sequence.id, len(sequence))
    gi|6273285|gb|AF191659.1|AF191 146
    gi|6273284|gb|AF191658.1|AF191 148
    gi|6273287|gb|AF191661.1|AF191 146
    gi|6273286|gb|AF191660.1|AF191 146
    gi|6273290|gb|AF191664.1|AF191 150
    gi|6273289|gb|AF191663.1|AF191 150
    gi|6273291|gb|AF191665.1|AF191 156

    If the file contains no records, or more than one record, an exception is
    raised.  For example:

    >>> from Bio import Align
    >>> filename = "Exonerate/exn_22_m_ner_cigar.exn"
    >>> alignment = Align.read(filename, "exonerate")
    Traceback (most recent call last):
        ...
    ValueError: More than one alignment found in file

    Use the Bio.Align.parse function if you want to read a file containing
    more than one alignment.
    """
    alignments = parse(handle, fmt)
    try:
        alignment = next(alignments)
    except StopIteration:
        raise ValueError("No alignments found in file") from None
    try:
        next(alignments)
        raise ValueError("More than one alignment found in file")
    except StopIteration:
        pass
    return alignment


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
