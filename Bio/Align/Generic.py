# Copyright 2000-2004 Brad Chapman.
# Copyright 2001 Iddo Friedberg.
# Copyright 2007-2010 by Peter Cock.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
Contains classes to deal with generic sequence alignment stuff not
specific to a particular program or format.

Classes:
 - Alignment
"""
__docformat__ = "epytext en" #Don't just use plain text in epydoc API pages!

# biopython
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Alphabet

class Alignment(object):
    """Represent a set of alignments (DEPRECATED).

    This is a base class to represent alignments, which can be subclassed
    to deal with an alignment in a specific format.

    With the introduction of the MultipleSeqAlignment class in Bio.Align,
    this base class is deprecated and is likely to be removed in future
    releases of Biopython.
    """
    def __init__(self, alphabet):
        """Initialize a new Alignment object.

        Arguments:
         - alphabet - The alphabet to use for the sequence objects that are
                      created. This alphabet must be a gapped type.

        e.g.

        >>> from Bio.Alphabet import IUPAC, Gapped
        >>> align = Alignment(Gapped(IUPAC.unambiguous_dna, "-"))
        >>> align.add_sequence("Alpha", "ACTGCTAGCTAG")
        >>> align.add_sequence("Beta",  "ACT-CTAGCTAG")
        >>> align.add_sequence("Gamma", "ACTGCTAGATAG")
        >>> print align
        Gapped(IUPACUnambiguousDNA(), '-') alignment with 3 rows and 12 columns
        ACTGCTAGCTAG Alpha
        ACT-CTAGCTAG Beta
        ACTGCTAGATAG Gamma
        """
        import warnings
        import Bio
        warnings.warn("With the introduction of the MultipleSeqAlignment class in Bio.Align, this base class is deprecated and is likely to be removed in a future release of Biopython.", Bio.BiopythonDeprecationWarning)
        if not (isinstance(alphabet, Alphabet.Alphabet) \
        or isinstance(alphabet, Alphabet.AlphabetEncoder)):
            raise ValueError("Invalid alphabet argument")
        self._alphabet = alphabet
        # hold everything at a list of SeqRecord objects
        self._records = []

    def _str_line(self, record):
        """Returns a truncated string representation of a SeqRecord (PRIVATE).

        This is a PRIVATE function used by the __str__ method.
        """
        if len(record.seq) <= 50:
            return "%s %s" % (record.seq, record.id)
        else:
            return "%s...%s %s" \
                   % (record.seq[:44], record.seq[-3:], record.id)

    def __str__(self):
        """Returns a multi-line string summary of the alignment.

        This output is intended to be readable, but large alignments are
        shown truncated.  A maximum of 20 rows (sequences) and 50 columns
        are shown, with the record identifiers.  This should fit nicely on a
        single screen.  e.g.

        >>> from Bio.Alphabet import IUPAC, Gapped
        >>> align = Alignment(Gapped(IUPAC.unambiguous_dna, "-"))
        >>> align.add_sequence("Alpha", "ACTGCTAGCTAG")
        >>> align.add_sequence("Beta",  "ACT-CTAGCTAG")
        >>> align.add_sequence("Gamma", "ACTGCTAGATAG")
        >>> print align
        Gapped(IUPACUnambiguousDNA(), '-') alignment with 3 rows and 12 columns
        ACTGCTAGCTAG Alpha
        ACT-CTAGCTAG Beta
        ACTGCTAGATAG Gamma

        See also the alignment's format method.
        """
        rows = len(self._records)
        lines = ["%s alignment with %i rows and %i columns" \
                 % (str(self._alphabet), rows, self.get_alignment_length())]
        if rows <= 20:
            lines.extend([self._str_line(rec) for rec in self._records])
        else:
            lines.extend([self._str_line(rec) for rec in self._records[:18]])
            lines.append("...")
            lines.append(self._str_line(self._records[-1]))
        return "\n".join(lines)

    def __repr__(self):
        """Returns a representation of the object for debugging.

        The representation cannot be used with eval() to recreate the object,
        which is usually possible with simple python ojects.  For example:

        <Bio.Align.Generic.Alignment instance (2 records of length 14,
        SingleLetterAlphabet()) at a3c184c>

        The hex string is the memory address of the object, see help(id).
        This provides a simple way to visually distinguish alignments of
        the same size.
        """
        #A doctest for __repr__ would be nice, but __class__ comes out differently
        #if run via the __main__ trick.
        return "<%s instance (%i records of length %i, %s) at %x>" % \
               (self.__class__, len(self._records),
                self.get_alignment_length(), repr(self._alphabet), id(self))
        #This version is useful for doing eval(repr(alignment)),
        #but it can be VERY long:
        #return "%s(%s, %s)" \
        #       % (self.__class__, repr(self._records), repr(self._alphabet))

    def format(self, format):
        """Returns the alignment as a string in the specified file format.

        The format should be a lower case string supported as an output
        format by Bio.AlignIO (such as "fasta", "clustal", "phylip",
        "stockholm", etc), which is used to turn the alignment into a
        string.

        e.g.

        >>> from Bio.Alphabet import IUPAC, Gapped
        >>> align = Alignment(Gapped(IUPAC.unambiguous_dna, "-"))
        >>> align.add_sequence("Alpha", "ACTGCTAGCTAG")
        >>> align.add_sequence("Beta",  "ACT-CTAGCTAG")
        >>> align.add_sequence("Gamma", "ACTGCTAGATAG")
        >>> print align.format("fasta")
        >Alpha
        ACTGCTAGCTAG
        >Beta
        ACT-CTAGCTAG
        >Gamma
        ACTGCTAGATAG
        <BLANKLINE>
        >>> print align.format("phylip")
         3 12
        Alpha      ACTGCTAGCT AG
        Beta       ACT-CTAGCT AG
        Gamma      ACTGCTAGAT AG
        <BLANKLINE>

        For Python 2.6, 3.0 or later see also the built in format() function.
        """
        #See also the __format__ added for Python 2.6 / 3.0, PEP 3101
        #See also the SeqRecord class and its format() method using Bio.SeqIO
        return self.__format__(format)


    def __format__(self, format_spec):
        """Returns the alignment as a string in the specified file format.

        This method supports the python format() function added in
        Python 2.6/3.0.  The format_spec should be a lower case
        string supported by Bio.AlignIO as an output file format.
        See also the alignment's format() method."""
        if format_spec:
            from StringIO import StringIO
            from Bio import AlignIO
            handle = StringIO()
            AlignIO.write([self], handle, format_spec)
            return handle.getvalue()
        else:
            #Follow python convention and default to using __str__
            return str(self)    

    def get_all_seqs(self):
        """Return all of the sequences involved in the alignment (DEPRECATED).

        The return value is a list of SeqRecord objects.

        This method is deprecated, as the Alignment object itself now offers
        much of the functionality of a list of SeqRecord objects (e.g.
        iteration or slicing to create a sub-alignment). Instead use the
        Python builtin function list, i.e. my_list = list(my_align)
        """
        import warnings
        import Bio
        warnings.warn("This method is deprecated, since the alignment object"
                      "now acts more like a list. Instead of calling "
                      "align.get_all_seqs() you can use list(align)",
                      Bio.BiopythonDeprecationWarning)
        return self._records

    def __iter__(self):
        """Iterate over alignment rows as SeqRecord objects.

        e.g.

        >>> from Bio.Alphabet import IUPAC, Gapped
        >>> align = Alignment(Gapped(IUPAC.unambiguous_dna, "-"))
        >>> align.add_sequence("Alpha", "ACTGCTAGCTAG")
        >>> align.add_sequence("Beta",  "ACT-CTAGCTAG")
        >>> align.add_sequence("Gamma", "ACTGCTAGATAG")
        >>> for record in align:
        ...    print record.id
        ...    print record.seq
        Alpha
        ACTGCTAGCTAG
        Beta
        ACT-CTAGCTAG
        Gamma
        ACTGCTAGATAG
        """
        return iter(self._records) 

    def get_seq_by_num(self, number):
        """Retrieve a sequence by row number (DEPRECATED).

        Returns:
         - A Seq object for the requested sequence.

        Raises:
         - IndexError - If the specified number is out of range.

        NOTE: This is a legacy method.  In new code where you need to access
        the rows of the alignment (i.e. the sequences) consider iterating
        over them or accessing them as SeqRecord objects.
        """
        import warnings
        import Bio
        warnings.warn("This is a legacy method and is likely to be removed in a future release of Biopython. In new code where you need to access the rows of the alignment (i.e. the sequences) consider iterating over them or accessing them as SeqRecord objects.", Bio.BiopythonDeprecationWarning)
        return self._records[number].seq

    def __len__(self):
        """Returns the number of sequences in the alignment.

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
        >>> align = Alignment(Gapped(IUPAC.unambiguous_dna, "-"))
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

    def add_sequence(self, descriptor, sequence, start = None, end = None,
                     weight = 1.0):
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
        """
        new_seq = Seq(sequence, self._alphabet)

        #We are now effectively using the SeqRecord's .id as
        #the primary identifier (e.g. in Bio.SeqIO) so we should
        #populate it with the descriptor.
        #For backwards compatibility, also store this in the
        #SeqRecord's description property.
        new_record = SeqRecord(new_seq,
                               id = descriptor,
                               description = descriptor)

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
        
    def get_column(self,col):
        """Returns a string containing a given column.

        e.g.

        >>> from Bio.Alphabet import IUPAC, Gapped
        >>> align = Alignment(Gapped(IUPAC.unambiguous_dna, "-"))
        >>> align.add_sequence("Alpha", "ACTGCTAGCTAG")
        >>> align.add_sequence("Beta",  "ACT-CTAGCTAG")
        >>> align.add_sequence("Gamma", "ACTGCTAGATAG")
        >>> align.get_column(0)
        'AAA'
        >>> align.get_column(3)
        'G-G'
        """
        #TODO - Support negative indices?
        col_str = ''
        assert col >= 0 and col <= self.get_alignment_length()
        for rec in self._records:
            col_str += rec.seq[col]
        return col_str

    def __getitem__(self, index):
        """Access part of the alignment.

        We'll use the following example alignment here for illustration:

        >>> from Bio.Alphabet import IUPAC, Gapped
        >>> align = Alignment(Gapped(IUPAC.unambiguous_dna, "-"))
        >>> align.add_sequence("Alpha",  "ACTGCTAGCTAG")
        >>> align.add_sequence("Beta",   "ACT-CTAGCTAG")
        >>> align.add_sequence("Gamma",  "ACTGCTAGATAG")
        >>> align.add_sequence("Delta",  "ACTGCTTGCTAG")
        >>> align.add_sequence("Epsilon","ACTGCTTGATAG")
        
        You can access a row of the alignment as a SeqRecord using an integer
        index (think of the alignment as a list of SeqRecord objects here):

        >>> first_record = align[0]
        >>> print first_record.id, first_record.seq
        Alpha ACTGCTAGCTAG
        >>> last_record = align[-1]
        >>> print last_record.id, last_record.seq
        Epsilon ACTGCTTGATAG

        You can also access use python's slice notation to create a sub-alignment
        containing only some of the SeqRecord objects:

        >>> sub_alignment = align[2:5]
        >>> print sub_alignment
        Gapped(IUPACUnambiguousDNA(), '-') alignment with 3 rows and 12 columns
        ACTGCTAGATAG Gamma
        ACTGCTTGCTAG Delta
        ACTGCTTGATAG Epsilon

        This includes support for a step, i.e. align[start:end:step], which
        can be used to select every second sequence:

        >>> sub_alignment = align[::2]
        >>> print sub_alignment
        Gapped(IUPACUnambiguousDNA(), '-') alignment with 3 rows and 12 columns
        ACTGCTAGCTAG Alpha
        ACTGCTAGATAG Gamma
        ACTGCTTGATAG Epsilon

        Or to get a copy of the alignment with the rows in reverse order:

        >>> rev_alignment = align[::-1]
        >>> print rev_alignment
        Gapped(IUPACUnambiguousDNA(), '-') alignment with 5 rows and 12 columns
        ACTGCTTGATAG Epsilon
        ACTGCTTGCTAG Delta
        ACTGCTAGATAG Gamma
        ACT-CTAGCTAG Beta
        ACTGCTAGCTAG Alpha

        Right now, these are the ONLY indexing operations supported.  The use of
        a second column based index is under discussion for a future update.
        """
        if isinstance(index, int):
            #e.g. result = align[x]
            #Return a SeqRecord
            return self._records[index]
        elif isinstance(index, slice):
            #e.g. sub_aling = align[i:j:k]
            #Return a new Alignment using only the specified records.
            #TODO - See Bug 2554 for changing the __init__ method
            #to allow us to do this more cleanly.
            sub_align = Alignment(self._alphabet)
            sub_align._records = self._records[index]
            return sub_align
        elif len(index)==2:
            raise TypeError("Row and Column indexing is not currently supported,"\
                            +"but may be in future.")
        else:
            raise TypeError("Invalid index type.")

def _test():
    """Run the Bio.Align.Generic module's doctests."""
    print "Running doctests..."
    import doctest
    doctest.testmod()
    print "Done"

if __name__ == "__main__":
    _test()
