# Copyright 2008-2011 by Peter Cock.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Code for dealing with sequence alignments.

One of the most important things in this module is the MultipleSeqAlignment
class, used in the Bio.AlignIO module.

"""
from __future__ import print_function

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord, _RestrictedDict
from Bio import Alphabet


class MultipleSeqAlignment(object):
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
    SingleLetterAlphabet() alignment with 7 rows and 156 columns
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
    SingleLetterAlphabet() alignment with 7 rows and 10 columns
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
    SingleLetterAlphabet() alignment with 7 rows and 20 columns
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

    def __init__(self, records, alphabet=None,
                 annotations=None, column_annotations=None):
        """Initialize a new MultipleSeqAlignment object.

        Arguments:
         - records - A list (or iterator) of SeqRecord objects, whose
                     sequences are all the same length.  This may be an be an
                     empty list.
         - alphabet - The alphabet for the whole alignment, typically a gapped
                      alphabet, which should be a super-set of the individual
                      record alphabets.  If omitted, a consensus alphabet is
                      used.
         - annotations - Information about the whole alignment (dictionary).
         - column_annotations - Per column annotation (restricted dictionary).
                      This holds Python sequences (lists, strings, tuples)
                      whose length matches the number of columns. A typical
                      use would be a secondary structure consensus string.

        You would normally load a MSA from a file using Bio.AlignIO, but you
        can do this from a list of SeqRecord objects too:

        >>> from Bio.Alphabet import generic_dna
        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> from Bio.Align import MultipleSeqAlignment
        >>> a = SeqRecord(Seq("AAAACGT", generic_dna), id="Alpha")
        >>> b = SeqRecord(Seq("AAA-CGT", generic_dna), id="Beta")
        >>> c = SeqRecord(Seq("AAAAGGT", generic_dna), id="Gamma")
        >>> align = MultipleSeqAlignment([a, b, c],
        ...                              annotations={"tool": "demo"},
        ...                              column_annotations={"stats": "CCCXCCC"})
        >>> print(align)
        DNAAlphabet() alignment with 3 rows and 7 columns
        AAAACGT Alpha
        AAA-CGT Beta
        AAAAGGT Gamma
        >>> align.annotations
        {'tool': 'demo'}
        >>> align.column_annotations
        {'stats': 'CCCXCCC'}

        NOTE - The older Bio.Align.Generic.Alignment class only accepted a
        single argument, an alphabet.  This is still supported via a backwards
        compatible "hack" so as not to disrupt existing scripts and users, but
        is deprecated and will be removed in a future release.
        """
        if isinstance(records, (Alphabet.Alphabet, Alphabet.AlphabetEncoder)):
            if alphabet is None:
                # TODO - Remove this backwards compatible mode!
                alphabet = records
                records = []
                import warnings
                from Bio import BiopythonDeprecationWarning
                warnings.warn("Invalid records argument: While the old "
                              "Bio.Align.Generic.Alignment class only "
                              "accepted a single argument (the alphabet), the "
                              "newer Bio.Align.MultipleSeqAlignment class "
                              "expects a list/iterator of SeqRecord objects "
                              "(which can be an empty list) and an optional "
                              "alphabet argument", BiopythonDeprecationWarning)
            else:
                raise ValueError("Invalid records argument")
        if alphabet is not None:
            if not isinstance(alphabet, (Alphabet.Alphabet, Alphabet.AlphabetEncoder)):
                raise ValueError("Invalid alphabet argument")
            self._alphabet = alphabet
        else:
            # Default while we add sequences, will take a consensus later
            self._alphabet = Alphabet.single_letter_alphabet

        self._records = []
        if records:
            self.extend(records)
            if alphabet is None:
                # No alphabet was given, take a consensus alphabet
                self._alphabet = Alphabet._consensus_alphabet(rec.seq.alphabet for
                                                              rec in self._records
                                                              if rec.seq is not None)

        # Annotations about the whole alignment
        if annotations is None:
            annotations = {}
        elif not isinstance(annotations, dict):
            raise TypeError("annotations argument should be a dict")
        self.annotations = annotations

        # Annotations about each colum of the alignment
        if column_annotations is None:
            column_annotations = {}
        # Handle this via the property set function which will validate it
        self.column_annotations = column_annotations

    def _set_per_column_annotations(self, value):
        if not isinstance(value, dict):
            raise TypeError("The per-column-annotations should be a "
                            "(restricted) dictionary.")
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
                raise ValueError("Can't set per-column-annotations without an alignment")

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
        doc="""Dictionary of per-letter-annotation for the sequence.""")

    def _str_line(self, record, length=50):
        """Return a truncated string representation of a SeqRecord (PRIVATE).

        This is a PRIVATE function used by the __str__ method.
        """
        if record.seq.__class__.__name__ == "CodonSeq":
            if len(record.seq) <= length:
                return "%s %s" % (record.seq, record.id)
            else:
                return "%s...%s %s" \
                    % (record.seq[:length - 3], record.seq[-3:], record.id)
        else:
            if len(record.seq) <= length:
                return "%s %s" % (record.seq, record.id)
            else:
                return "%s...%s %s" \
                    % (record.seq[:length - 6], record.seq[-3:], record.id)

    def __str__(self):
        """Return a multi-line string summary of the alignment.

        This output is intended to be readable, but large alignments are
        shown truncated.  A maximum of 20 rows (sequences) and 50 columns
        are shown, with the record identifiers.  This should fit nicely on a
        single screen. e.g.

        >>> from Bio.Alphabet import IUPAC, Gapped
        >>> from Bio.Align import MultipleSeqAlignment
        >>> align = MultipleSeqAlignment([], Gapped(IUPAC.unambiguous_dna, "-"))
        >>> align.add_sequence("Alpha", "ACTGCTAGCTAG")
        >>> align.add_sequence("Beta",  "ACT-CTAGCTAG")
        >>> align.add_sequence("Gamma", "ACTGCTAGATAG")
        >>> print(align)
        Gapped(IUPACUnambiguousDNA(), '-') alignment with 3 rows and 12 columns
        ACTGCTAGCTAG Alpha
        ACT-CTAGCTAG Beta
        ACTGCTAGATAG Gamma

        See also the alignment's format method.
        """
        rows = len(self._records)
        lines = ["%s alignment with %i rows and %i columns"
                 % (str(self._alphabet), rows, self.get_alignment_length())]
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
        which is usually possible with simple python ojects.  For example:

        <Bio.Align.MultipleSeqAlignment instance (2 records of length 14,
        SingleLetterAlphabet()) at a3c184c>

        The hex string is the memory address of the object, see help(id).
        This provides a simple way to visually distinguish alignments of
        the same size.
        """
        # A doctest for __repr__ would be nice, but __class__ comes out differently
        # if run via the __main__ trick.
        return "<%s instance (%i records of length %i, %s) at %x>" % \
            (self.__class__, len(self._records),
             self.get_alignment_length(), repr(self._alphabet), id(self))
        # This version is useful for doing eval(repr(alignment)),
        # but it can be VERY long:
        # return "%s(%s, %s)" \
        #       % (self.__class__, repr(self._records), repr(self._alphabet))

    def format(self, format):
        """Return the alignment as a string in the specified file format.

        The format should be a lower case string supported as an output
        format by Bio.AlignIO (such as "fasta", "clustal", "phylip",
        "stockholm", etc), which is used to turn the alignment into a
        string.

        e.g.

        >>> from Bio.Alphabet import IUPAC, Gapped
        >>> from Bio.Align import MultipleSeqAlignment
        >>> align = MultipleSeqAlignment([], Gapped(IUPAC.unambiguous_dna, "-"))
        >>> align.add_sequence("Alpha", "ACTGCTAGCTAG")
        >>> align.add_sequence("Beta",  "ACT-CTAGCTAG")
        >>> align.add_sequence("Gamma", "ACTGCTAGATAG")
        >>> print(align.format("fasta"))
        >Alpha
        ACTGCTAGCTAG
        >Beta
        ACT-CTAGCTAG
        >Gamma
        ACTGCTAGATAG
        <BLANKLINE>
        >>> print(align.format("phylip"))
         3 12
        Alpha      ACTGCTAGCT AG
        Beta       ACT-CTAGCT AG
        Gamma      ACTGCTAGAT AG
        <BLANKLINE>

        For Python 2.6, 3.0 or later see also the built in format() function.
        """
        # See also the __format__ added for Python 2.6 / 3.0, PEP 3101
        # See also the SeqRecord class and its format() method using Bio.SeqIO
        return self.__format__(format)

    def __format__(self, format_spec):
        """Return the alignment as a string in the specified file format.

        This method supports the python format() function added in
        Python 2.6/3.0.  The format_spec should be a lower case
        string supported by Bio.AlignIO as an output file format.
        See also the alignment's format() method.
        """
        if format_spec:
            from Bio._py3k import StringIO
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

        >>> from Bio.Alphabet import IUPAC, Gapped
        >>> from Bio.Align import MultipleSeqAlignment
        >>> align = MultipleSeqAlignment([], Gapped(IUPAC.unambiguous_dna, "-"))
        >>> align.add_sequence("Alpha", "ACTGCTAGCTAG")
        >>> align.add_sequence("Beta",  "ACT-CTAGCTAG")
        >>> align.add_sequence("Gamma", "ACTGCTAGATAG")
        >>> for record in align:
        ...    print(record.id)
        ...    print(record.seq)
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

        >>> from Bio.Alphabet import IUPAC, Gapped
        >>> from Bio.Align import MultipleSeqAlignment
        >>> align = MultipleSeqAlignment([], Gapped(IUPAC.unambiguous_dna, "-"))
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

    def add_sequence(self, descriptor, sequence, start=None, end=None,
                     weight=1.0):
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
        new_seq = Seq(sequence, self._alphabet)

        # We are now effectively using the SeqRecord's .id as
        # the primary identifier (e.g. in Bio.SeqIO) so we should
        # populate it with the descriptor.
        # For backwards compatibility, also store this in the
        # SeqRecord's description property.
        new_record = SeqRecord(new_seq,
                               id=descriptor,
                               description=descriptor)

        # hack! We really need to work out how to deal with annotations
        # and features in biopython. Right now, I'll just use the
        # generic annotations dictionary we've got to store the start
        # and end, but we should think up something better. I don't know
        # if I'm really a big fan of the LocatableSeq thing they've got
        # in BioPerl, but I'm not positive what the best thing to do on
        # this is...
        if start:
            new_record.annotations['start'] = start
        if end:
            new_record.annotations['end'] = end

        # another hack to add weight information to the sequence
        new_record.annotations['weight'] = weight

        self._records.append(new_record)

    def extend(self, records):
        """Add more SeqRecord objects to the alignment as rows.

        They must all have the same length as the original alignment, and have
        alphabets compatible with the alignment's alphabet. For example,

        >>> from Bio.Alphabet import generic_dna
        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> from Bio.Align import MultipleSeqAlignment
        >>> a = SeqRecord(Seq("AAAACGT", generic_dna), id="Alpha")
        >>> b = SeqRecord(Seq("AAA-CGT", generic_dna), id="Beta")
        >>> c = SeqRecord(Seq("AAAAGGT", generic_dna), id="Gamma")
        >>> d = SeqRecord(Seq("AAAACGT", generic_dna), id="Delta")
        >>> e = SeqRecord(Seq("AAA-GGT", generic_dna), id="Epsilon")

        First we create a small alignment (three rows):

        >>> align = MultipleSeqAlignment([a, b, c])
        >>> print(align)
        DNAAlphabet() alignment with 3 rows and 7 columns
        AAAACGT Alpha
        AAA-CGT Beta
        AAAAGGT Gamma

        Now we can extend this alignment with another two rows:

        >>> align.extend([d, e])
        >>> print(align)
        DNAAlphabet() alignment with 5 rows and 7 columns
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
        the first record), and have an alphabet compatible with the alignment's
        alphabet.

        >>> from Bio import AlignIO
        >>> align = AlignIO.read("Clustalw/opuntia.aln", "clustal")
        >>> print(align)
        SingleLetterAlphabet() alignment with 7 rows and 156 columns
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
        SingleLetterAlphabet() alignment with 8 rows and 156 columns
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
            # raise ValueError("New sequence is not of length %i" \
            #                 % self.get_alignment_length())
            raise ValueError("Sequences must all be the same length")

        # Using not self.alphabet.contains(record.seq.alphabet) needs fixing
        # for AlphabetEncoders (e.g. gapped versus ungapped).
        if not Alphabet._check_type_compatible([self._alphabet, record.seq.alphabet]):
            raise ValueError("New sequence's alphabet is incompatible")
        self._records.append(record)

    def __add__(self, other):
        """Combine two alignments with the same number of rows by adding them.

        If you have two multiple sequence alignments (MSAs), there are two ways to think
        about adding them - by row or by column. Using the extend method adds by row.
        Using the addition operator adds by column. For example,

        >>> from Bio.Alphabet import generic_dna
        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> from Bio.Align import MultipleSeqAlignment
        >>> a1 = SeqRecord(Seq("AAAAC", generic_dna), id="Alpha")
        >>> b1 = SeqRecord(Seq("AAA-C", generic_dna), id="Beta")
        >>> c1 = SeqRecord(Seq("AAAAG", generic_dna), id="Gamma")
        >>> a2 = SeqRecord(Seq("GT", generic_dna), id="Alpha")
        >>> b2 = SeqRecord(Seq("GT", generic_dna), id="Beta")
        >>> c2 = SeqRecord(Seq("GT", generic_dna), id="Gamma")
        >>> left = MultipleSeqAlignment([a1, b1, c1],
        ...                             annotations={"tool": "demo", "name": "start"},
        ...                             column_annotations={"stats": "CCCXC"})
        >>> right = MultipleSeqAlignment([a2, b2, c2],
        ...                             annotations={"tool": "demo", "name": "end"},
        ...                             column_annotations={"stats": "CC"})

        Now, let's look at these two alignments:

        >>> print(left)
        DNAAlphabet() alignment with 3 rows and 5 columns
        AAAAC Alpha
        AAA-C Beta
        AAAAG Gamma
        >>> print(right)
        DNAAlphabet() alignment with 3 rows and 2 columns
        GT Alpha
        GT Beta
        GT Gamma

        And add them:

        >>> combined = left + right
        >>> print(combined)
        DNAAlphabet() alignment with 3 rows and 7 columns
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
            raise ValueError("When adding two alignments they must have the same length"
                             " (i.e. same number or rows)")
        alpha = Alphabet._consensus_alphabet([self._alphabet, other._alphabet])
        merged = (left + right for left, right in zip(self, other))
        # Take any common annotation:
        annotations = dict()
        for k, v in self.annotations.items():
            if k in other.annotations and other.annotations[k] == v:
                annotations[k] = v
        column_annotations = dict()
        for k, v in self.column_annotations.items():
            if k in other.column_annotations:
                column_annotations[k] = v + other.column_annotations[k]
        return MultipleSeqAlignment(merged, alpha, annotations, column_annotations)

    def __getitem__(self, index):
        """Access part of the alignment.

        Depending on the indices, you can get a SeqRecord object
        (representing a single row), a Seq object (for a single columns),
        a string (for a single characters) or another alignment
        (representing some part or all of the alignment).

        align[r,c] gives a single character as a string
        align[r] gives a row as a SeqRecord
        align[r,:] gives a row as a SeqRecord
        align[:,c] gives a column as a Seq (using the alignment's alphabet)

        align[:] and align[:,:] give a copy of the alignment

        Anything else gives a sub alignment, e.g.
        align[0:2] or align[0:2,:] uses only row 0 and 1
        align[:,1:3] uses only columns 1 and 2
        align[0:2,1:3] uses only rows 0 & 1 and only cols 1 & 2

        We'll use the following example alignment here for illustration:

        >>> from Bio.Alphabet import generic_dna
        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> from Bio.Align import MultipleSeqAlignment
        >>> a = SeqRecord(Seq("AAAACGT", generic_dna), id="Alpha")
        >>> b = SeqRecord(Seq("AAA-CGT", generic_dna), id="Beta")
        >>> c = SeqRecord(Seq("AAAAGGT", generic_dna), id="Gamma")
        >>> d = SeqRecord(Seq("AAAACGT", generic_dna), id="Delta")
        >>> e = SeqRecord(Seq("AAA-GGT", generic_dna), id="Epsilon")
        >>> align = MultipleSeqAlignment([a, b, c, d, e], generic_dna)

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
        DNAAlphabet() alignment with 3 rows and 7 columns
        AAAAGGT Gamma
        AAAACGT Delta
        AAA-GGT Epsilon

        This includes support for a step, i.e. align[start:end:step], which
        can be used to select every second sequence:

        >>> sub_alignment = align[::2]
        >>> print(sub_alignment)
        DNAAlphabet() alignment with 3 rows and 7 columns
        AAAACGT Alpha
        AAAAGGT Gamma
        AAA-GGT Epsilon

        Or to get a copy of the alignment with the rows in reverse order:

        >>> rev_alignment = align[::-1]
        >>> print(rev_alignment)
        DNAAlphabet() alignment with 5 rows and 7 columns
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
        DNAAlphabet() alignment with 4 rows and 3 columns
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
            new = MultipleSeqAlignment(self._records[index], self._alphabet)
            if self.column_annotations and len(new) == len(self):
                # All rows kept (although could have been reversed)
                # Perserve the column annotations too,
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
            new = MultipleSeqAlignment((rec[col_index] for rec in self._records[row_index]),
                                       self._alphabet)
            if self.column_annotations and len(new) == len(self):
                # All rows kept (although could have been reversed)
                # Perserve the column annotations too,
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

        >>> from Bio.Alphabet import generic_dna
        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> from Bio.Align import MultipleSeqAlignment
        >>> align1 = MultipleSeqAlignment([
        ...              SeqRecord(Seq("ACGT", generic_dna), id="Human"),
        ...              SeqRecord(Seq("ACGG", generic_dna), id="Mouse"),
        ...              SeqRecord(Seq("ACGC", generic_dna), id="Chicken"),
        ...          ])
        >>> align2 = MultipleSeqAlignment([
        ...              SeqRecord(Seq("CGGT", generic_dna), id="Mouse"),
        ...              SeqRecord(Seq("CGTT", generic_dna), id="Human"),
        ...              SeqRecord(Seq("CGCT", generic_dna), id="Chicken"),
        ...          ])

        If you simple try and add these without sorting, you get this:

        >>> print(align1 + align2)
        DNAAlphabet() alignment with 3 rows and 8 columns
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
        DNAAlphabet() alignment with 3 rows and 8 columns
        ACGCCGCT Chicken
        ACGTCGTT Human
        ACGGCGGT Mouse

        As an example using a different sort order, you could sort on the
        GC content of each sequence.

        >>> from Bio.SeqUtils import GC
        >>> print(align1)
        DNAAlphabet() alignment with 3 rows and 4 columns
        ACGC Chicken
        ACGT Human
        ACGG Mouse
        >>> align1.sort(key = lambda record: GC(record.seq))
        >>> print(align1)
        DNAAlphabet() alignment with 3 rows and 4 columns
        ACGT Human
        ACGC Chicken
        ACGG Mouse

        There is also a reverse argument, so if you wanted to sort by ID
        but backwards:

        >>> align1.sort(reverse=True)
        >>> print(align1)
        DNAAlphabet() alignment with 3 rows and 4 columns
        ACGG Mouse
        ACGT Human
        ACGC Chicken

        """
        if key is None:
            self._records.sort(key=lambda r: r.id, reverse=reverse)
        else:
            self._records.sort(key=key, reverse=reverse)


if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()
