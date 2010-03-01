# Copyright 2008-2010 by Peter Cock.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Code for dealing with sequence alignments.
"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Alphabet

#We only import this and subclass it for some limited backward compatibilty.
from Bio.Align.Generic import Alignment as _Alignment
class MultipleSeqAlignment(_Alignment):
    """"Represents a classical multiple sequence alignment (MSA).

    By this we mean a collection of sequences (usually shown as rows) which
    are all the same length (usually with gap characters for insertions of
    padding). The data can then be regarded as a matrix of letters, with well
    defined columns.

    You would typically create an MSA by loading an alignment file with the
    AlignIO module:

    >>> from Bio import AlignIO
    >>> align = AlignIO.read("Clustalw/opuntia.aln", "clustal")
    >>> print align
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
    ...     print record.id, len(record)
    gi|6273285|gb|AF191659.1|AF191 156
    gi|6273284|gb|AF191658.1|AF191 156
    gi|6273287|gb|AF191661.1|AF191 156
    gi|6273286|gb|AF191660.1|AF191 156
    gi|6273290|gb|AF191664.1|AF191 156
    gi|6273289|gb|AF191663.1|AF191 156
    gi|6273291|gb|AF191665.1|AF191 156

    You can also access individual rows as SeqRecord objects via their index:

    >>> print align[0].id
    gi|6273285|gb|AF191659.1|AF191
    >>> print align[-1].id
    gi|6273291|gb|AF191665.1|AF191

    Or, take just the first ten columns as a sub-alignment:

    >>> print align[:,:10]
    SingleLetterAlphabet() alignment with 7 rows and 10 columns
    TATACATTAA gi|6273285|gb|AF191659.1|AF191
    TATACATTAA gi|6273284|gb|AF191658.1|AF191
    TATACATTAA gi|6273287|gb|AF191661.1|AF191
    TATACATAAA gi|6273286|gb|AF191660.1|AF191
    TATACATTAA gi|6273290|gb|AF191664.1|AF191
    TATACATTAA gi|6273289|gb|AF191663.1|AF191
    TATACATTAA gi|6273291|gb|AF191665.1|AF191
    
    Combining this alignment slicing with alignment addtion allows you to
    remove a section of the alignment. For example, taking just the first
    and last ten columns:

    >>> print align[:,:10] + align[:,-10:]
    SingleLetterAlphabet() alignment with 7 rows and 20 columns
    TATACATTAAGTGTACCAGA gi|6273285|gb|AF191659.1|AF191
    TATACATTAAGTGTACCAGA gi|6273284|gb|AF191658.1|AF191
    TATACATTAAGTGTACCAGA gi|6273287|gb|AF191661.1|AF191
    TATACATAAAGTGTACCAGA gi|6273286|gb|AF191660.1|AF191
    TATACATTAAGTGTACCAGA gi|6273290|gb|AF191664.1|AF191
    TATACATTAAGTATACCAGA gi|6273289|gb|AF191663.1|AF191
    TATACATTAAGTGTACCAGA gi|6273291|gb|AF191665.1|AF191
    
    Note - This object is intended to replace the existing Alignment object
    defined in module Bio.Align.Generic but is not fully backwards compatible
    with it.

    Note - This object does NOT attempt to model the kind of alignments used
    in next generation sequencing with multiple sequencing reads which are
    much shorter than the alignment, and where there is usually a consensus or
    reference sequence with special status.
    """

    def __init__(self, records, alphabet=None):
        """Initialize a new MultipleSeqAlignment object.

        Arguments:
        records - A list (or iterator) of SeqRecord objects, whose sequences
                  are all the same length.  This may be an be an empty list.
        alphabet - The alphabet for the whole alignment, typically a gapped
                  alphabet, which should be a super-set of the individual
                  record alphabets.  If omitted, a consensus alphabet is used.

        You would normally load a MSA from a file using Bio.AlignIO, but you
        can do this from a list of SeqRecord objects too:

        >>> from Bio.Alphabet import generic_dna
        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> a = SeqRecord(Seq("AAAACGT", generic_dna), id="Alpha")
        >>> b = SeqRecord(Seq("AAA-CGT", generic_dna), id="Beta")
        >>> c = SeqRecord(Seq("AAAAGGT", generic_dna), id="Gamma")
        >>> align = MultipleSeqAlignment([a, b, c])
        >>> print align
        DNAAlphabet() alignment with 3 rows and 7 columns
        AAAACGT Alpha
        AAA-CGT Beta
        AAAAGGT Gamma

        NOTE - The older Bio.Align.Generic.Alignment class only accepted a
        single argument, an alphabet.  This is still supported via a backwards
        compatible "hack" so as not to disrupt existing scripts and users, but
        this will in future be deprecated.
        """
        if isinstance(records, Alphabet.Alphabet) \
        or isinstance(records, Alphabet.AlphabetEncoder):
            if alphabet is None:
                #TODO - Deprecate this backwards compatible mode!                
                alphabet = records
                records = []
            else :
                raise ValueError("Invalid records argument")
        if alphabet is not None :
            if not (isinstance(alphabet, Alphabet.Alphabet) \
            or isinstance(alphabet, Alphabet.AlphabetEncoder)):
                raise ValueError("Invalid alphabet argument")
            self._alphabet = alphabet
        else :
            #Default while we add sequences, will take a consensus later
            self._alphabet = Alphabet.single_letter_alphabet

        self._records = []
        if records:
            self.extend(records)
            if alphabet is None:
                #No alphabet was given, take a consensus alphabet
                self._alphabet = Alphabet._consensus_alphabet(rec.seq.alphabet for \
                                                              rec in self._records \
                                                              if rec.seq is not None)

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
        >>> print align
        DNAAlphabet() alignment with 3 rows and 7 columns
        AAAACGT Alpha
        AAA-CGT Beta
        AAAAGGT Gamma

        Now we can extend this alignment with another two rows:

        >>> align.extend([d, e])
        >>> print align
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
        for rec in records:
            self.append(rec)

    def append(self, record):
        """Add one more SeqRecord object to the alignment as a new row.

        This must have the same length as the original alignment (unless this is
        the first record), and have an alphabet compatible with the alignment's
        alphabet.

        >>> from Bio import AlignIO
        >>> align = AlignIO.read("Clustalw/opuntia.aln", "clustal")
        >>> print align
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
        >>> print align
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
        if not isinstance(record, SeqRecord):
            raise TypeError("New sequence is not a SeqRecord object")
        if self._records and len(record) != self.get_alignment_length():
            #TODO - Use the following more helpful error, but update unit tests
            #raise ValueError("New sequence is not of length %i" \
            #                 % self.get_alignment_length())
            raise ValueError("Sequences must all be the same length")
        #Using not self.alphabet.contains(record.seq.alphabet) needs fixing
        #for AlphabetEncoders (e.g. gapped versus ungapped).
        if not Alphabet._check_type_compatible([self._alphabet, record.seq.alphabet]):
            raise ValueError("New sequence's alphabet is incompatible")
        self._records.append(record)

    def __add__(self, other):
        """Combines to alignments with the same number of rows by adding them.

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
        >>> left = MultipleSeqAlignment([a1, b1, c1])
        >>> right = MultipleSeqAlignment([a2, b2, c2])

        Now, let's look at these two alignments:

        >>> print left
        DNAAlphabet() alignment with 3 rows and 5 columns
        AAAAC Alpha
        AAA-C Beta
        AAAAG Gamma
        >>> print right
        DNAAlphabet() alignment with 3 rows and 2 columns
        GT Alpha
        GT Beta
        GT Gamma

        And add them:

        >>> print left + right
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

        The individual rows are SeqRecord objects, and these can be added together. Refer
        to the SeqRecord documentation for details of how the annotation is handled. This
        example is a special case in that both original alignments shared the same names,
        meaning when the rows are added they also get the same name.
        """
        if not isinstance(other, MultipleSeqAlignment):
            raise NotImplementedError
        if len(self) != len(other):
            raise ValueError("When adding two alignments they must have the same length"
                             " (i.e. same number or rows)")
        alpha = Alphabet._consensus_alphabet([self._alphabet, other._alphabet])
        merged = (left+right for left,right in zip(self, other))
        return MultipleSeqAlignment(merged, alpha)

    def __getitem__(self, index):
        """Access part of the alignment.

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
        >>> print first_record.id, first_record.seq
        Alpha AAAACGT
        >>> last_record = align[-1]
        >>> print last_record.id, last_record.seq
        Epsilon AAA-GGT

        You can also access use python's slice notation to create a sub-alignment
        containing only some of the SeqRecord objects:

        >>> sub_alignment = align[2:5]
        >>> print sub_alignment
        DNAAlphabet() alignment with 3 rows and 7 columns
        AAAAGGT Gamma
        AAAACGT Delta
        AAA-GGT Epsilon

        This includes support for a step, i.e. align[start:end:step], which
        can be used to select every second sequence:

        >>> sub_alignment = align[::2]
        >>> print sub_alignment
        DNAAlphabet() alignment with 3 rows and 7 columns
        AAAACGT Alpha
        AAAAGGT Gamma
        AAA-GGT Epsilon

        Or to get a copy of the alignment with the rows in reverse order:

        >>> rev_alignment = align[::-1]
        >>> print rev_alignment
        DNAAlphabet() alignment with 5 rows and 7 columns
        AAA-GGT Epsilon
        AAAACGT Delta
        AAAAGGT Gamma
        AAA-CGT Beta
        AAAACGT Alpha
    
        You can also use two indices to specify both rows and columns. e.g.

        >>> print align[3,4]
        C

        This is equivalent to:

        >>> print align[3][4]
        C

        or:

        >>> print align[3].seq[4]
        C

        To get a single column, use this syntax:

        >>> print align[:,4]
        CCGCG

        In general you get a sub-alignment,

        >>> print align[1:5,3:6]
        DNAAlphabet() alignment with 4 rows and 3 columns
        -CG Beta
        AGG Gamma
        ACG Delta
        -GG Epsilon

        This should all seem familiar to anyone who has used the NumPy
        array or matrix operations.
        """
        if isinstance(index, int):
            #e.g. result = align[x]
            #Return a SeqRecord
            return self._records[index]
        elif isinstance(index, slice):
            #e.g. sub_align = align[i:j:k]
            return MultipleSeqAlignment(self._records[index], self._alphabet)
        elif len(index)!=2:
            raise TypeError("Invalid index type.")

        #Handle double indexing
        row_index, col_index = index
        if isinstance(row_index, int):
            #e.g. row_or_part_row = align[6, 1:4], gives a SeqRecord
            return self._records[row_index][col_index]
        elif isinstance(col_index, int):
            #e.g. col_or_part_col = align[1:5, 6], gives a Seq
            return Seq("".join(rec[col_index] for rec in self._records[row_index]),
                       self._alphabet)
        else:
            #e.g. sub_align = align[1:4, 5:7], gives another alignment
            return MultipleSeqAlignment((rec[col_index] for rec in self._records[row_index]),
                                        self._alphabet)

    
def _test():
    """Run the Bio.Align module's doctests.

    This will try and locate the unit tests directory, and run the doctests
    from there in order that the relative paths used in the examples work.
    """
    import doctest
    import os
    if os.path.isdir(os.path.join("..", "..", "Tests", "Clustalw")):
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join("..", "..", "Tests"))
        doctest.testmod()
        os.chdir(cur_dir)
        del cur_dir
        print "Done"
    elif os.path.isdir(os.path.join("Tests", "Clustalw")):
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join("Tests"))
        doctest.testmod()
        os.chdir(cur_dir)
        del cur_dir
        print "Done"

if __name__ == "__main__":
    #Run the doctests
    _test()
