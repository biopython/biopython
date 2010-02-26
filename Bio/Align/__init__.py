# Copyright 2008-2010 by Peter Cock.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Code for dealing with sequence alignments.
"""

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

        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> a = SeqRecord(Seq("AAAACGT"), id="Alpha")
        >>> b = SeqRecord(Seq("AAA-CGT"), id="Beta")
        >>> c = SeqRecord(Seq("AAAAGGT"), id="Gamma")
        >>> align = MultipleSeqAlignment([a, b, c])
        >>> print align
        SingleLetterAlphabet() alignment with 3 rows and 7 columns
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
            if alphabet is None :
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
            assert len(records) == len(self)
            if alphabet is None :
                #No alphabet was given, take a consensus alphabet
                self.alphabet = Alphabet._consensus_alphabet(rec.seq.alphabet for \
                                                             rec in self._records)

    def extend(self, records) :
        """Add more SeqRecord objects to the alignment as rows.

        They must all have the same length as the original alignment, and have
        alphabets compatible with the alignment's alphabet. For example,

        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> a = SeqRecord(Seq("AAAACGT"), id="Alpha")
        >>> b = SeqRecord(Seq("AAA-CGT"), id="Beta")
        >>> c = SeqRecord(Seq("AAAAGGT"), id="Gamma")
        >>> d = SeqRecord(Seq("AAAACGT"), id="Delta")
        >>> e = SeqRecord(Seq("AAA-GGT"), id="Epsilon")

        First we create a small alignment (three rows):

        >>> align = MultipleSeqAlignment([a, b, c])
        >>> print align
        SingleLetterAlphabet() alignment with 3 rows and 7 columns
        AAAACGT Alpha
        AAA-CGT Beta
        AAAAGGT Gamma

        Now we can extend this alignment with another two rows:

        >>> align.extend([d, e])
        >>> print align
        SingleLetterAlphabet() alignment with 5 rows and 7 columns
        AAAACGT Alpha
        AAA-CGT Beta
        AAAAGGT Gamma
        AAAACGT Delta
        AAA-GGT Epsilon

        Because the alignment object allows iteration over the rows as
        SeqRecords, you can use the extend method with a second alignment
        (provided its sequences have the same length as the original alignment).
        """
        for rec in records :
            self.append(rec)

    def append(self, record) :
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
            raise ValueError("New sequence is not of length %i" \
                             % self.get_alignment_length())
        #Using not self._alphabet.contains(record.seq.alphabet) needs fixing
        #for AlphabetEncoders (e.g. gapped versus ungapped).
        if not Alphabet._check_type_compatible([self._alphabet, record.seq.alphabet]):
            raise ValueError("New sequence's alphabet is incompatible")
        self._records.append(record)
    
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
