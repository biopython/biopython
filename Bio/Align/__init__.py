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
import warnings
import numbers

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

    Note - This object replaced the older Alignment object defined in module
    Bio.Align.Generic but is not fully backwards compatible with it.

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

        >>> from Bio.SeqUtils import GC
        >>> print(align1)
        Alignment with 3 rows and 4 columns
        ACGC Chicken
        ACGT Human
        ACGG Mouse
        >>> align1.sort(key = lambda record: GC(record.seq))
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

    Internally, the alignment is stored as a numpy array containing the
    sequence coordinates defining the alignment.

    Commonly used attributes (which may or may not be present):
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
         - sequences   - A list of the sequences that were aligned.
         - coordinates - The sequence coordinates that define the alignment.
                         If None (the default value), assume that the sequences
                         align to each other without any gaps.
        """
        self.sequences = sequences
        if coordinates is None:
            lengths = {len(sequence) for sequence in sequences}
            if len(lengths) != 1:
                raise ValueError(
                    "sequences must have the same length if coordinates is None"
                )
            length = lengths.pop()
            coordinates = numpy.array([[0, length]] * len(sequences))
        self.coordinates = coordinates

    @property
    def target(self):
        """Return self.sequences[0] for a pairwise alignment."""
        n = len(self.sequences)
        if n != 2:
            raise ValueError(
                "self.target is defined for pairwise alignments only (found alignment of % sequences)"
                % n
            )
        return self.sequences[0]

    @target.setter
    def target(self, value):
        """For a pairwise alignment, set self.sequences[0]."""
        n = len(self.sequences)
        if n != 2:
            raise ValueError(
                "self.target is defined for pairwise alignments only (found alignment of % sequences)"
                % n
            )
        self.sequences[0] = value

    @property
    def query(self):
        """Return self.sequences[1] for a pairwise alignment."""
        n = len(self.sequences)
        if n != 2:
            raise ValueError(
                "self.query is defined for pairwise alignments only (found alignment of % sequences)"
                % n
            )
        return self.sequences[1]

    @query.setter
    def query(self, value):
        """For a pairwise alignment, set self.sequences[1]."""
        n = len(self.sequences)
        if n != 2:
            raise ValueError(
                "self.query is defined for pairwise alignments only (found alignment of % sequences)"
                % n
            )
        self.sequences[1] = value

    def __eq__(self, other):
        """Check if two Alignment objects specify the same alignment."""
        from itertools import zip_longest

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
        from itertools import zip_longest

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
        from itertools import zip_longest

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
        from itertools import zip_longest

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
        from itertools import zip_longest

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
        from itertools import zip_longest

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
        """
        coordinates = self.coordinates.copy()
        for sequence, row in zip(self.sequences, coordinates):
            if row[0] > row[-1]:  # reverse strand
                row[:] = len(sequence) - row[:]
        steps = numpy.diff(coordinates, 1)
        gaps = steps.max(0)
        if not ((steps == gaps) | (steps == 0)).all():
            raise ValueError("Unequal step sizes in alignment")
        sequence = self.sequences[index]
        try:
            sequence = sequence.seq  # SeqRecord confusion
        except AttributeError:
            pass
        coordinates = coordinates[index]
        if self.coordinates[index, 0] > self.coordinates[index, -1]:
            # reverse strand
            sequence = reverse_complement(sequence, inplace=False)
        line = ""
        steps = steps[index]
        i = coordinates[0]
        for step, gap in zip(steps, gaps):
            if step:
                j = i + step
                line += str(sequence[i:j])
                i = j
            else:
                line += "-" * gap
        return line

    def _get_rows(self, key):
        """Return self[key], where key is a slice object (PRIVATE).

        This method is called by __getitem__ for invocations of the form

        self[rows]

        where rows is a slice object.
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
        """
        sequence_indices = coordinate + steps.cumsum()
        indices = gaps.cumsum()
        i = indices.searchsorted(start_index, side="right")
        j = i + indices[i:].searchsorted(stop_index, side="right")
        if i == j:
            length = stop_index - start_index
            if steps[i] == 0:
                line = "-" * length
            else:
                offset = start_index - indices[i]
                start = sequence_indices[i] + offset
                stop = start + length
                line = str(sequence[start:stop])
        else:
            length = indices[i] - start_index
            stop = sequence_indices[i]
            if steps[i] == 0:
                line = "-" * length
            else:
                start = stop - length
                line = str(sequence[start:stop])
            i += 1
            while i < j:
                step = gaps[i]
                if steps[i] == 0:
                    line += "-" * step
                else:
                    start = stop
                    stop = start + step
                    line += str(sequence[start:stop])
                i += 1
            length = stop_index - indices[j - 1]
            if length > 0:
                if steps[j] == 0:
                    line += "-" * length
                else:
                    start = stop
                    stop = start + length
                    line += str(sequence[start:stop])
        return line

    def _get_row_cols_iterable(self, i, cols, steps, gaps, sequence):
        """Return the alignment contents of one row and multiple columns (PRIVATE).

        This method is called by __getitem__ for invocations of the form

        self[row, cols]

        where row is an integer and cols is an iterable of integers.
        """
        line = ""
        for step, gap in zip(steps, gaps):
            if step:
                j = i + step
                line += str(sequence[i:j])
                i = j
            else:
                line += "-" * gap
        try:
            line = "".join(line[col] for col in cols)
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
        self, coordinates, row, start_index, stop_index, steps, gaps, sequences
    ):
        """Return a subalignment of multiple rows and consecutive columns (PRIVATE).

        This method is called by __getitem__ for invocations of the form

        self[rows, cols]

        where rows is an arbitrary slice object, and cols is a slice object
        with step 1, allowing the alignment sequences to be reused in the
        subalignment.
        """
        indices = gaps.cumsum()
        i = indices.searchsorted(start_index, side="right")
        j = i + indices[i:].searchsorted(stop_index, side="left") + 1
        offset = steps[:, i] - indices[i] + start_index
        coordinates[:, i] += offset * (steps[:, i] > 0)
        offset = indices[j - 1] - stop_index
        coordinates[:, j] -= offset * (steps[:, j - 1] > 0)
        coordinates = coordinates[:, i : j + 1]
        reverse_complemented = self.coordinates[row, 0] > self.coordinates[row, -1]
        for i, sequence in enumerate(sequences):
            if reverse_complemented[i]:
                # mapped to reverse strand
                coordinates[i, :] = len(sequence) - coordinates[i, :]
        alignment = Alignment(self.sequences[row], coordinates)
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
        object.
        """
        indices = tuple(col)
        lines = []
        for i, sequence in enumerate(sequences):
            line = ""
            k = coordinates[i, 0]
            for step, gap in zip(steps[i], gaps):
                if step:
                    j = k + step
                    line += str(sequence[k:j])
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
        sequences = [line.replace("-", "") for line in lines]
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
        ACCGG-TTT
        ||-||-||-
        AC-GGGTT-
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
        <Bio.Align.Alignment object (2 rows x 8 columns) at 0x...>
        >>> print(alignment[:, 1:])
        ACCGG-TTT
         |-||-||-
        AC-GGGTT-
        <BLANKLINE>
        >>> print(alignment[:, 2:])
        ACCGG-TTT
          -||-||-
        AC-GGGTT-
        <BLANKLINE>
        >>> print(alignment[:, 3:])
        ACCGG-TTT
           ||-||-
         ACGGGTT-
        <BLANKLINE>
        >>> print(alignment[:, 3:-1])
        ACCGG-TTT
           ||-||
         ACGGGTT
        <BLANKLINE>
        >>> print(alignment[:, ::2])
        ACGTT
        |-||-
        A-GT-
        <BLANKLINE>
        >>> print(alignment[:, range(1, 9, 2)])
        CG-T
        ||-|
        CGGT
        <BLANKLINE>
        >>> print(alignment[:, (2, 7, 3)])
        CTG
        -||
        -TG
        <BLANKLINE>
        """
        if isinstance(key, numbers.Integral):
            return self._get_row(key)
        if isinstance(key, slice):
            return self._get_rows(key)
        sequences = list(self.sequences)
        coordinates = self.coordinates.copy()
        for i, sequence in enumerate(sequences):
            try:
                sequence = sequence.seq  # SeqRecord confusion
            except AttributeError:
                pass
            if coordinates[i, 0] > coordinates[i, -1]:  # reverse strand
                coordinates[i, :] = len(sequence) - coordinates[i, :]
                sequence = reverse_complement(sequence, inplace=False)
            sequences[i] = sequence
        steps = numpy.diff(coordinates, 1)
        gaps = steps.max(0)
        if not ((steps == gaps) | (steps == 0)).all():
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
            coordinate = coordinates[row, 0]
            if isinstance(col, numbers.Integral):
                return self._get_row_col(coordinate, col, steps, gaps, sequence)
            if isinstance(col, slice):
                start_index, stop_index, step = col.indices(m)
                if start_index < stop_index and step == 1:
                    return self._get_row_cols_slice(
                        coordinate, start_index, stop_index, steps, gaps, sequence
                    )
                # make an iterable if step != 1
                col = range(start_index, stop_index, step)
            return self._get_row_cols_iterable(coordinate, col, steps, gaps, sequence)
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
                        sequences,
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

        Wrapper for self.format() .
        """
        return self.format(format_spec)

    def format(self, fmt="", **kwargs):
        """Return the alignment as a string in the specified file format.

        Arguments:
         - fmt       - File format. Acceptable values are
                       ""   : create a human-readable representation of the
                              alignment (default);
                       "bed": create a line representing the alignment in
                              the Browser Extensible Data (BED) file format;
                       "psl": create a line representing the alignment in
                              the Pattern Space Layout (PSL) file format as
                              generated by BLAT;
                       "sam": create a line representing the alignment in
                              the Sequence Alignment/Map (SAM) format.
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
        """
        if len(self.sequences) > 2:
            raise NotImplementedError(
                "format is currently implemented for pairwise alignments only"
            )
        if fmt == "":
            return self._format_pretty(**kwargs)
        elif fmt == "psl":
            return self._format_psl(**kwargs)
        elif fmt == "bed":
            return self._format_bed(**kwargs)
        elif fmt == "sam":
            return self._format_sam(**kwargs)
        else:
            raise ValueError("Unknown format %s" % fmt)

    def _format_pretty(self):
        """Return default string representation (PRIVATE).

        Helper for self.format() .
        """
        target, query = self.sequences
        seq1 = self._convert_sequence_string(target)
        if seq1 is None:
            return self._format_generalized()
        seq2 = self._convert_sequence_string(query)
        if seq2 is None:
            return self._format_generalized()
        n1 = len(seq1)
        n2 = len(seq2)
        aligned_seq1 = ""
        aligned_seq2 = ""
        pattern = ""
        coordinates = self.coordinates
        if coordinates[0, 0] > coordinates[0, -1]:  # mapped to reverse strand
            coordinates = coordinates.copy()
            coordinates[0, :] = n1 - coordinates[0, :]
            seq2 = reverse_complement(seq2, inplace=False)
        if coordinates[1, 0] > coordinates[1, -1]:  # mapped to reverse strand
            coordinates = coordinates.copy()
            coordinates[1, :] = n2 - coordinates[1, :]
            seq2 = reverse_complement(seq2, inplace=False)
        coordinates = coordinates.transpose()
        end1, end2 = coordinates[0, :]
        if end1 > 0 or end2 > 0:
            end = max(end1, end2)
            aligned_seq1 += " " * (end - end1) + seq1[:end1]
            aligned_seq2 += " " * (end - end2) + seq2[:end2]
            pattern += " " * end
        start1 = end1
        start2 = end2
        for end1, end2 in coordinates[1:]:
            if end1 == start1:
                gap = end2 - start2
                aligned_seq1 += "-" * gap
                aligned_seq2 += seq2[start2:end2]
                pattern += "-" * gap
            elif end2 == start2:
                gap = end1 - start1
                aligned_seq1 += seq1[start1:end1]
                aligned_seq2 += "-" * gap
                pattern += "-" * gap
            else:
                s1 = seq1[start1:end1]
                s2 = seq2[start2:end2]
                if len(s1) != len(s2):
                    raise ValueError("Unequal step sizes in alignment")
                aligned_seq1 += s1
                aligned_seq2 += s2
                for c1, c2 in zip(s1, s2):
                    if c1 == c2:
                        pattern += "|"
                    else:
                        pattern += "."
            start1 = end1
            start2 = end2
        aligned_seq1 += seq1[end1:]
        aligned_seq2 += seq2[end2:]
        return f"{aligned_seq1}\n{pattern}\n{aligned_seq2}\n"

    def _format_generalized(self):
        """Return generalized string representation (PRIVATE).

        Helper for self._format_pretty() .
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

    def _format_bed(self, bedN=12):
        """Return BED file format string representation (PRIVATE).

        Helper for self.format().
        """
        from . import bed

        writer = bed.AlignmentWriter(None, bedN=bedN)
        return writer.format_alignment(self)

    def _format_psl(self, mask=False, wildcard="N"):
        """Return PSL file format string representation (PRIVATE).

        Helper for self.format() .
        """
        from . import psl

        writer = psl.AlignmentWriter(None, header=False, mask=mask, wildcard=wildcard)
        return writer.format_alignment(self)

    def _format_sam(self):
        """Return SAM file format string representation (PRIVATE).

        Helper for self.format() .
        """
        target, query = self.sequences
        try:
            qName = query.id
        except AttributeError:
            qName = "query"
        else:
            query = query.seq
        try:
            rName = target.id
        except AttributeError:
            rName = "target"
        else:
            target = target.seq
        n1 = len(target)
        n2 = len(query)
        pos = None
        tSize = n1
        cigar = []
        coordinates = self.coordinates
        if coordinates[1, 0] < coordinates[1, -1]:  # mapped to forward strand
            flag = 0
            seq = query
        else:  # mapped to reverse strand
            flag = 16
            seq = reverse_complement(query, inplace=False)
            coordinates = coordinates.copy()
            coordinates[1, :] = n2 - coordinates[1, :]
        try:
            seq = bytes(seq)
        except TypeError:  # string
            pass
        else:
            seq = str(seq, "ASCII")
        tStart, qStart = coordinates[:, 0]
        for tEnd, qEnd in coordinates[:, 1:].transpose():
            tCount = tEnd - tStart
            qCount = qEnd - qStart
            if tCount == 0:
                length = qCount
                if pos is None or tEnd == tSize:
                    operation = "S"
                else:
                    operation = "I"
                qStart = qEnd
            elif qCount == 0:
                if tStart > 0 and tEnd < tSize:
                    length = tCount
                    operation = "D"
                else:
                    operation = None
                tStart = tEnd
            else:
                if tCount != qCount:
                    raise ValueError("Unequal step sizes in alignment")
                if pos is None:
                    pos = tStart
                tStart = tEnd
                qStart = qEnd
                operation = "M"
                length = tCount
            if operation is not None:
                cigar.append(str(length) + operation)
        mapQ = 255  # not available
        rNext = "*"
        pNext = 0
        tLen = 0
        qual = "*"
        cigar = "".join(cigar)
        tag = "AS:i:%d" % int(round(self.score))
        words = [
            qName,
            str(flag),
            rName,
            str(pos + 1),  # 1-based coordinates
            str(mapQ),
            cigar,
            rNext,
            str(pNext),
            str(tLen),
            seq,
            qual,
            tag,
        ]
        line = "\t".join(words) + "\n"
        return line

    def __str__(self):
        """Return a string representation of the Alignment object.

        Wrapper for self.format().
        """
        if len(self.sequences) > 2:
            raise NotImplementedError(
                "__str__ is currently implemented for pairwise alignments only"
            )
        return self.format()

    def __repr__(self):
        """Return a representation of the alignment, including its shape.

        The representation cannot be used with eval() to recreate the object,
        which is usually possible with simple python objects.  For example:

        <Bio.Align.Alignment object (2 rows x 14 columns) at 0x10403d850>

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
        ACCGT
        ||-||
        AC-GT
        <BLANKLINE>
        >>> alignment  # doctest:+ELLIPSIS
        <Bio.Align.Alignment object (2 rows x 5 columns) at 0x...>
        """
        n, m = self.shape
        return "<%s.%s object (%i rows x %i columns) at 0x%x>" % (
            self.__module__,
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
        -GACCT-G
        -||--|-|
        CGA--TCG
        <BLANKLINE>
        >>> len(alignment)
        2
        >>> alignment.shape
        (2, 8)
        >>> aligner.mode = "local"
        >>> alignments = aligner.align("GACCTG", "CGATCG")
        >>> alignment = alignments[0]
        >>> print(alignment)
         GACCT-G
         ||--|-|
        CGA--TCG
        <BLANKLINE>
        >>> len(alignment)
        2
        >>> alignment.shape
        (2, 7)
        """
        coordinates = self.coordinates.copy()
        n = len(coordinates)
        for i in range(n):
            if coordinates[i, 0] > coordinates[i, -1]:  # mapped to reverse strand
                k = len(self.sequences[i])
                coordinates[i, :] = k - coordinates[i, :]
        steps = numpy.diff(coordinates, 1)
        gaps = steps.max(0)
        if not ((steps == gaps) | (steps == 0)).all():
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
        GAACT
        ||--|
        GA--T
        <BLANKLINE>
        >>> alignment.aligned
        array([[[0, 2],
                [4, 5]],
        <BLANKLINE>
               [[0, 2],
                [2, 3]]])
        >>> alignment = alignments[1]
        >>> print(alignment)
        GAACT
        |-|-|
        G-A-T
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
        AAAC-AAA
        |||--|||
        AAA-GAAA
        <BLANKLINE>
        >>> alignments[0].aligned
        array([[[0, 3],
                [4, 7]],
        <BLANKLINE>
               [[0, 3],
                [4, 7]]])
        >>> print(alignments[1])
        AAA-CAAA
        |||--|||
        AAAG-AAA
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
        for i, sequence in enumerate(self.sequences):
            if coordinates[i, 0] > coordinates[i, -1]:  # mapped to reverse strand
                n = len(sequence)
                coordinates[i, :] = n - coordinates[i, :]
        coordinates = coordinates.transpose()
        steps = numpy.diff(coordinates, axis=0).min(1)
        indices = numpy.flatnonzero(steps)
        starts = coordinates[indices, :]
        ends = coordinates[indices + 1, :]
        segments = numpy.stack([starts, ends], axis=0).transpose()
        for i, sequence in enumerate(self.sequences):
            if self.coordinates[i, 0] > self.coordinates[i, -1]:
                # mapped to reverse strand
                n = len(sequence)
                segments[i, :] = n - segments[i, :]
        return segments

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
        AATAA
        ||.||
        AAGAA
        <BLANKLINE>
        >>> alignment.sort()
        >>> print(alignment)
        AAGAA
        ||.||
        AATAA
        <BLANKLINE>

        Alternatively, a key function can be supplied that maps each sequence
        to a sort value.  For example, you could sort on the GC content of each
        sequence.

        >>> from Bio.SeqUtils import GC
        >>> alignment.sort(key=GC)
        >>> print(alignment)
        AATAA
        ||.||
        AAGAA
        <BLANKLINE>

        You can reverse the sort order by passing `reverse=True`:

        >>> alignment.sort(key=GC, reverse=True)
        >>> print(alignment)
        AAGAA
        ||.||
        AATAA
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
        AAAAAAAACCCCCCCAAAAAAAAAAAGGGGGGAAAAAAAA
                |||||||-----------||||||
                CCCCCCC-----------GGGGGG
        <BLANKLINE>
        >>> sequence = "CCCCGGGG"
        >>> alignments2 = aligner.align(transcript, sequence)
        >>> len(alignments2)
        1
        >>> alignment2 = alignments2[0]
        >>> print(alignment2)
        CCCCCCCGGGGGG
           ||||||||
           CCCCGGGG
        <BLANKLINE>
        >>> alignment = alignment1.map(alignment2)
        >>> print(alignment)
        AAAAAAAACCCCCCCAAAAAAAAAAAGGGGGGAAAAAAAA
                   ||||-----------||||
                   CCCC-----------GGGG
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
        if coordinates1[1, 0] < coordinates1[1, -1]:  # mapped to forward strand
            strand1 = "+"
        else:  # mapped to reverse strand
            strand1 = "-"
        if coordinates2[1, 0] < coordinates2[1, -1]:  # mapped to forward strand
            strand2 = "+"
        else:  # mapped to reverse strand
            strand2 = "-"
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
        if not ((steps1 == gaps1) | (steps1 == 0)).all():
            raise ValueError("Unequal step sizes in first alignment")
        steps2 = numpy.diff(coordinates2, 1)
        gaps2 = steps2.max(0)
        if not ((steps2 == gaps2) | (steps2 == 0)).all():
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
        ATACTTACCTGGCAGGGGAGATACCATGATCACGAAGGTGGTTTTCCCAGGGCGAGGCTTATCCATTGCACTCCGGATGTGCTGACCCCTGCGATTTCCCCAAATGTGGGAAACTCGACTGCATAATTTGTGGTAGTGGGGGACTGCGTTCGCGCTTTCCCCTG
        |||||||||||.||||||||..|||||||||||..|||||||..|||||||||||||||..|||||||||||.|||..|.|.|||||||||||||||||||||||||||||||||||||||.||||||||||||||||||||||||||||||||||.|||||.|
        ATACTTACCTGACAGGGGAGGCACCATGATCACACAGGTGGTCCTCCCAGGGCGAGGCTCTTCCATTGCACTGCGGGAGGGTTGACCCCTGCGATTTCCCCAAATGTGGGAAACTCGACTGTATAATTTGTGGTAGTGGGGGACTGCGTTCGCGCTATCCCCCG
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
        sequences = self.sequences
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
            coordinates1 = self.coordinates[i1, :]
            if coordinates1[0] > coordinates1[-1]:
                sequence1 = reverse_complement(sequence1)
                coordinates1 = len(sequence1) - coordinates1
            for i2 in range(i1 + 1, n):
                sequence2 = sequences[i2]
                coordinates2 = self.coordinates[i2, :]
                if coordinates2[0] > coordinates2[-1]:
                    sequence2 = reverse_complement(sequence2)
                    coordinates2 = len(sequence2) - coordinates2
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
        self.paths = paths
        self.index = -1

    def __len__(self):
        """Return the number of alignments."""
        return len(self.paths)

    def __getitem__(self, index):
        if index == self.index:
            return self.alignment
        if index < self.index:
            self.paths.reset()
            self.index = -1
        while self.index < index:
            try:
                alignment = next(self)
            except StopIteration:
                raise IndexError("index out of range") from None
        return alignment

    def __iter__(self):
        self.paths.reset()
        self.index = -1
        return self

    def __next__(self):
        path = next(self.paths)
        self.index += 1
        coordinates = numpy.array(path)
        alignment = Alignment(self.sequences, coordinates)
        alignment.score = self.score
        self.alignment = alignment
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
    TACCG
    -|-||
    -A-CG
    <BLANKLINE>
    Score = 3.0:
    TACCG
    -||-|
    -AC-G
    <BLANKLINE>

    Specify the aligner mode as local to generate local alignments:

    >>> aligner.mode = 'local'
    >>> alignments = aligner.align("TACCG", "ACG")
    >>> for alignment in sorted(alignments):
    ...     print("Score = %.1f:" % alignment.score)
    ...     print(alignment)
    ...
    Score = 3.0:
    TACCG
     |-||
     A-CG
    <BLANKLINE>
    Score = 3.0:
    TACCG
     ||-|
     AC-G
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
    TACCG
    -||-|
    -AC-G
    <BLANKLINE>
    Score = 6.0:
    TACCG
    -|-||
    -A-CG
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
    TACCG
    -|-||
    -A-CG
    <BLANKLINE>
    Score = 5.5:
    TACCG
    -||-|
    -AC-G
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
    KEVLA
    -|||-
    -EVL-
    <BLANKLINE>

    You can also set the value of attributes directly during construction
    of the PairwiseAligner object by providing them as keyword arguments:

    >>> aligner = Align.PairwiseAligner(mode='global', match_score=2, mismatch_score=-1)
    >>> for alignment in aligner.align("TACCG", "ACG"):
    ...     print("Score = %.1f:" % alignment.score)
    ...     print(alignment)
    ...
    Score = 6.0:
    TACCG
    -||-|
    -AC-G
    <BLANKLINE>
    Score = 6.0:
    TACCG
    -|-||
    -A-CG
    <BLANKLINE>

    """

    def __init__(self, **kwargs):
        """Initialize a new PairwiseAligner with the keyword arguments as attributes.

        Loops over the keyword arguments and sets them as attributes on the object.
        """
        super().__init__()
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
        if isinstance(seqA, (Seq, MutableSeq)):
            sA = bytes(seqA)
        else:
            sA = seqA
        if strand == "+":
            sB = seqB
        else:  # strand == "-":
            sB = reverse_complement(seqB, inplace=False)
        if isinstance(sB, (Seq, MutableSeq)):
            sB = bytes(sB)
        score, paths = _aligners.PairwiseAligner.align(self, sA, sB, strand)
        alignments = PairwiseAlignments(seqA, seqB, score, paths)
        return alignments

    def score(self, seqA, seqB, strand="+"):
        """Return the alignments score of two sequences using PairwiseAligner."""
        if isinstance(seqA, (Seq, MutableSeq)):
            seqA = bytes(seqA)
        if strand == "-":
            seqB = reverse_complement(seqB, inplace=False)
        if isinstance(seqB, (Seq, MutableSeq)):
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


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
