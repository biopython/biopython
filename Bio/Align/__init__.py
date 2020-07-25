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


import warnings

from Bio import BiopythonDeprecationWarning
from Bio.Align import _aligners
from Bio.Align import substitution_matrices
from Bio.Seq import Seq
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
                return "%s %s" % (record.seq, record.id)
            else:
                return "%s...%s %s" % (
                    record.seq[: length - 3],
                    record.seq[-3:],
                    record.id,
                )
        else:
            if len(record.seq) <= length:
                return "%s %s" % (record.seq, record.id)
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

        >>> from Bio.Align import MultipleSeqAlignment
        >>> align = MultipleSeqAlignment([])
        >>> align.add_sequence("Alpha", "ACTGCTAGCTAG")
        >>> align.add_sequence("Beta",  "ACT-CTAGCTAG")
        >>> align.add_sequence("Gamma", "ACTGCTAGATAG")
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

    def format(self, format_spec):
        """Return the alignment as a string in the specified file format [DEPRECATED].

        This method is deprecated; instead of alignment.format(format_spec),
        please use format(alignment, format_spec) or an f-string.
        """
        warnings.warn(
            "alignment.format has been deprecated, and we intend to remove it in a "
            "future release of Biopython. Instead of alignment.format(format_spec), "
            "please use format(alignment, format_spec) or an f-string.",
            BiopythonDeprecationWarning,
        )
        return self.__format__(format_spec)

    def __format__(self, format_spec):
        """Return the alignment as a string in the specified file format.

        The format should be a lower case string supported as an output
        format by Bio.AlignIO (such as "fasta", "clustal", "phylip",
        "stockholm", etc), which is used to turn the alignment into a
        string.

        e.g.

        >>> from Bio.Align import MultipleSeqAlignment
        >>> align = MultipleSeqAlignment([])
        >>> align.add_sequence("Alpha", "ACTGCTAGCTAG")
        >>> align.add_sequence("Beta",  "ACT-CTAGCTAG")
        >>> align.add_sequence("Gamma", "ACTGCTAGATAG")
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
        >>> align = MultipleSeqAlignment([])
        >>> align.add_sequence("Alpha", "ACTGCTAGCTAG")
        >>> align.add_sequence("Beta",  "ACT-CTAGCTAG")
        >>> align.add_sequence("Gamma", "ACTGCTAGATAG")
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
        >>> align = MultipleSeqAlignment([])
        >>> align.add_sequence("Alpha", "ACTGCTAGCTAG")
        >>> align.add_sequence("Beta",  "ACT-CTAGCTAG")
        >>> align.add_sequence("Gamma", "ACTGCTAGATAG")
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

    def add_sequence(self, descriptor, sequence, start=None, end=None, weight=1.0):
        """Add a sequence to the alignment.

        This doesn't do any kind of alignment, it just adds in the sequence
        object, which is assumed to be prealigned with the existing
        sequences.

        Arguments:
            - descriptor - The descriptive id of the sequence being added.
              This will be used as the resulting SeqRecord's
              .id property (and, for historical compatibility,
              also the .description property)
            - sequence - A string with sequence info.
            - start - You can explicitly set the start point of the sequence.
              This is useful (at least) for BLAST alignments, which can
              just be partial alignments of sequences.
            - end - Specify the end of the sequence, which is important
              for the same reason as the start.
            - weight - The weight to place on the sequence in the alignment.
              By default, all sequences have the same weight. (0.0 =>
              no weight, 1.0 => highest weight)

        In general providing a SeqRecord and calling .append is preferred.
        """
        new_seq = Seq(sequence)

        # We are now effectively using the SeqRecord's .id as
        # the primary identifier (e.g. in Bio.SeqIO) so we should
        # populate it with the descriptor.
        # For backwards compatibility, also store this in the
        # SeqRecord's description property.
        new_record = SeqRecord(new_seq, id=descriptor, description=descriptor)

        # hack! We really need to work out how to deal with annotations
        # and features in biopython. Right now, I'll just use the
        # generic annotations dictionary we've got to store the start
        # and end, but we should think up something better. I don't know
        # if I'm really a big fan of the LocatableSeq thing they've got
        # in BioPerl, but I'm not positive what the best thing to do on
        # this is...
        if start:
            new_record.annotations["start"] = start
        if end:
            new_record.annotations["end"] = end

        # another hack to add weight information to the sequence
        new_record.annotations["weight"] = weight

        self._records.append(new_record)

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
                (rec[col_index] for rec in self._records[row_index])
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
        letters = set.union(*[set(record.seq) for record in self])
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


class PairwiseAlignment:
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
        self.target = target
        self.query = query
        self.score = score
        self.path = path

    def __eq__(self, other):
        return self.path == other.path

    def __ne__(self, other):
        return self.path != other.path

    def __lt__(self, other):
        return self.path < other.path

    def __le__(self, other):
        return self.path <= other.path

    def __gt__(self, other):
        return self.path > other.path

    def __ge__(self, other):
        return self.path >= other.path

    def _convert_sequence_string(self, sequence):
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
        """Create a human-readable representation of the alignment."""
        if format_spec == "psl":
            return self._format_psl()
        seq1 = self._convert_sequence_string(self.target)
        if seq1 is None:
            return self._format_generalized()
        seq2 = self._convert_sequence_string(self.query)
        if seq2 is None:
            return self._format_generalized()
        n1 = len(seq1)
        n2 = len(seq2)
        aligned_seq1 = ""
        aligned_seq2 = ""
        pattern = ""
        path = self.path
        end1, end2 = path[0]
        if end1 > 0 or end2 > 0:
            end = max(end1, end2)
            aligned_seq1 += " " * (end - end1) + seq1[:end1]
            aligned_seq2 += " " * (end - end2) + seq2[:end2]
            pattern += " " * end
        start1 = end1
        start2 = end2
        for end1, end2 in path[1:]:
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
                aligned_seq1 += s1
                aligned_seq2 += s2
                for c1, c2 in zip(s1, s2):
                    if c1 == c2:
                        pattern += "|"
                    else:
                        pattern += "."
            start1 = end1
            start2 = end2
        n1 -= end1
        n2 -= end2
        n = max(n1, n2)
        aligned_seq1 += seq1[end1:] + " " * (n - n1)
        aligned_seq2 += seq2[end2:] + " " * (n - n2)
        pattern += " " * n
        return "%s\n%s\n%s\n" % (aligned_seq1, pattern, aligned_seq2)

    def _format_generalized(self):
        seq1 = self.target
        seq2 = self.query
        aligned_seq1 = []
        aligned_seq2 = []
        pattern = []
        path = self.path
        end1, end2 = path[0]
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
        for end1, end2 in path[1:]:
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
                for c1, c2 in zip(seq1[start1:end1], seq2[start2:end2]):
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
        return "%s\n%s\n%s\n" % (aligned_seq1, pattern, aligned_seq2)

    def _format_psl(self):
        query = self.query
        target = self.target
        try:
            qName = query.id
        except AttributeError:
            qName = "query"
        else:
            query = query.seq
        try:
            tName = target.id
        except AttributeError:
            tName = "target"
        else:
            target = target.seq
        seq1 = str(target)
        seq2 = str(query)
        n1 = len(seq1)
        n2 = len(seq2)
        # variable names follow those in the PSL file format specification
        matches = 0
        misMatches = 0
        repMatches = 0
        nCount = 0
        qNumInsert = 0
        qBaseInsert = 0
        tNumInsert = 0
        tBaseInsert = 0
        qSize = n2
        tSize = n1
        blockSizes = []
        qStarts = []
        tStarts = []
        strand = "+"
        start1, start2 = self.path[0]
        tStart = start1
        qStart = start2
        for end1, end2 in self.path[1:]:
            count1 = end1 - start1
            count2 = end2 - start2
            if count1 == 0:
                qNumInsert += 1
                qBaseInsert += count2
                start2 = end2
            elif count2 == 0:
                tNumInsert += 1
                tBaseInsert += count1
                start1 = end1
            else:
                assert count1 == count2
                tStarts.append(start1)
                qStarts.append(start2)
                blockSizes.append(count1)
                for c1, c2 in zip(seq1[start1:end1], seq2[start2:end2]):
                    if c1 == "N" or c2 == "N":
                        nCount += 1
                    elif c1 == c2:
                        matches += 1
                    else:
                        misMatches += 1
                start1 = end1
                start2 = end2
        tEnd = end1
        qEnd = end2
        blockcount = len(blockSizes)
        blockSizes = ",".join(map(str, blockSizes)) + ","
        qStarts = ",".join(map(str, qStarts)) + ","
        tStarts = ",".join(map(str, tStarts)) + ","
        words = [
            str(matches),
            str(misMatches),
            str(repMatches),
            str(nCount),
            str(qNumInsert),
            str(qBaseInsert),
            str(tNumInsert),
            str(tBaseInsert),
            strand,
            qName,
            str(qSize),
            str(qStart),
            str(qEnd),
            tName,
            str(tSize),
            str(tStart),
            str(tEnd),
            str(blockcount),
            blockSizes,
            qStarts,
            tStarts,
        ]
        line = "\t".join(words) + "\n"
        return line

    def format(self):
        """Create a human-readable representation of the alignment (DEPRECATED).

        This method is deprecated; instead of alignment.format(), please use
        format(alignment) or an f-string.
        """
        warnings.warn(
            "alignment.format has been deprecated, and we intend to remove it in a "
            "future release of Biopython. Instead of alignment.format(), please use "
            "format(alignment) or an f-string.",
            BiopythonDeprecationWarning,
        )
        return self.__format__(None)

    def __str__(self):
        return self.__format__(None)

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
        (((0, 2), (4, 5)), ((0, 2), (2, 3)))
        >>> alignment = alignments[1]
        >>> print(alignment)
        GAACT
        |-|-|
        G-A-T
        <BLANKLINE>
        >>> alignment.aligned
        (((0, 1), (2, 3), (4, 5)), ((0, 1), (1, 2), (2, 3)))

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
        (((0, 3), (4, 7)), ((0, 3), (4, 7)))
        >>> print(alignments[1])
        AAA-CAAA
        |||--|||
        AAAG-AAA
        <BLANKLINE>
        >>> alignments[1].aligned
        (((0, 3), (4, 7)), ((0, 3), (4, 7)))

        The property can be used to identify alignments that are identical
        to each other in terms of their aligned sequences.
        """
        segments1 = []
        segments2 = []
        i1, i2 = self.path[0]
        for node in self.path[1:]:
            j1, j2 = node
            if j1 > i1 and j2 > i2:
                segment1 = (i1, j1)
                segment2 = (i2, j2)
                segments1.append(segment1)
                segments2.append(segment2)
            i1, i2 = j1, j2
        return tuple(segments1), tuple(segments2)


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

        You would normally obtain an PairwiseAlignments object by calling
        aligner.align(seqA, seqB), where aligner is a PairwiseAligner object.
        """
        self.seqA = seqA
        self.seqB = seqB
        self.score = score
        self.paths = paths
        self.index = -1

    def __len__(self):
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
        alignment = PairwiseAlignment(self.seqA, self.seqB, path, self.score)
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
    of the PairwiseAligner object by providing them as keyword arguemnts:

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

        This function subclasses `_aligners.PairwiseAligner` and loops over all
        the keyword arguments that are given in the constructor to set them
        as attributes on the object. This will call the `__setattr__` method to
        do that.
        """
        super().__init__()
        for name, value in kwargs.items():
            setattr(self, name, value)

    def __setattr__(self, key, value):
        if key not in dir(_aligners.PairwiseAligner):
            # To prevent confusion, don't allow users to create new attributes
            raise AttributeError("PairwiseAligner object has no attribute '%s'" % key)
        _aligners.PairwiseAligner.__setattr__(self, key, value)

    def align(self, seqA, seqB):
        """Return the alignments of two sequences using PairwiseAligner."""
        if isinstance(seqA, Seq):
            seqA = str(seqA)
        if isinstance(seqB, Seq):
            seqB = str(seqB)
        score, paths = _aligners.PairwiseAligner.align(self, seqA, seqB)
        alignments = PairwiseAlignments(seqA, seqB, score, paths)
        return alignments

    def score(self, seqA, seqB):
        """Return the alignments score of two sequences using PairwiseAligner."""
        if isinstance(seqA, Seq):
            seqA = str(seqA)
        if isinstance(seqB, Seq):
            seqB = str(seqB)
        return _aligners.PairwiseAligner.score(self, seqA, seqB)


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
