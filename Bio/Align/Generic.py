"""
Contains classes to deal with generic sequence alignment stuff not
specific to a particular program or format.

classes:
o Alignment
"""

# standard library
import string

# biopython
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Alphabet
from Bio.Alphabet import IUPAC

class Alignment:
    """Represent a set of alignments.

    This is a base class to represent alignments, which can be subclassed
    to deal with an alignment in a specific format.
    """
    def __init__(self, alphabet):
        """Initialize a new Alignment object.

        Arguments:
        o alphabet - The alphabet to use for the sequence objects that are
        created. This alphabet must be a gapped type.
        """
        if not (isinstance(alphabet, Alphabet.Alphabet) \
        or isinstance(alphabet, Alphabet.AlphabetEncoder)):
            raise ValueError("Invalid alphabet argument")
        self._alphabet = alphabet
        # hold everything at a list of seq record objects
        self._records = []

    def _str_line(self, record) :
        """Returns a truncated string representation of a SeqRecord (PRIVATE).

        This is a PRIVATE function used by the __str__ method.
        """
        if len(record.seq) <= 50 :
            return "%s %s" % (record.seq, record.id)
        else :
            return "%s...%s %s" \
                   % (record.seq[:44], record.seq[-3:], record.id)

    def __str__(self) :
        """Returns a multi-line string summary of the alignment.

        This output is intended to be readable, but large alignments are
        shown truncated.  A maximum of 20 rows (sequences) and 50 columns
        are shown, with the record identifiers.  This should fit nicely on a
        single screen.  e.g.

        DNAAlphabet() alignment with 3 rows and 14 columns
        ACGATCAGCTAGCT Alpha
        CCGATCAGCTAGCT Beta
        ACGATGAGCTAGCT Gamma

        See also the alignment's format method.
        """
        rows = len(self._records)
        lines = ["%s alignment with %i rows and %i columns" \
                 % (str(self._alphabet), rows, self.get_alignment_length())]
        if rows <= 20 :
            lines.extend([self._str_line(rec) for rec in self._records])
        else :
            lines.extend([self._str_line(rec) for rec in self._records[:18]])
            lines.append("...")
            lines.append(self._str_line(self._records[-1]))
        return "\n".join(lines)

    def __repr__(self) :
        """Returns a representation of the object for debugging.

        The representation cannot be used with eval() to recreate the object,
        which is usually possible with simple python ojects.  For example:

        <Bio.Align.Generic.Alignment instance (2 records of length 14,
        SingleLetterAlphabet()) at a3c184c>

        The hex string is the memory address of the object, see help(id).
        This provides a simple way to visually distinguish alignments of
        the same size.
        """
        return "<%s instance (%i records of length %i, %s) at %x>" % \
               (self.__class__, len(self._records),
                self.get_alignment_length(), repr(self._alphabet), id(self))
        #This version is useful for doing eval(repr(alignment)),
        #but it can be VERY long:
        #return "%s(%s, %s)" \
        #       % (self.__class__, repr(self._records), repr(self._alphabet))

    def format(self, format) :
        """Returns the alignment as a string in the specified file format.

        The format should be a lower case string supported as an output
        format by Bio.AlignIO, which is used to turn the alignment into a
        string.

        e.g.
        print my_alignment.format("clustal")
        print my_alignment.format("fasta")
        """
        #See also the __format__ added for Python 2.6 / 3.0, PEP 3101
        #See also the SeqRecord class and its format() method using Bio.SeqIO
        return self.__format__(format)


    def __format__(self, format_spec) :
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
            handle.seek(0)
            return handle.read()
        else :
            #Follow python convention and default to using __str__
            return str(self)    

    def get_all_seqs(self):
        """Return all of the sequences involved in the alignment.

        The return value is a list of SeqRecord objects.

        This method is semi-obsolete, as the Alignment object itself offers
        much of the functionality of a list of SeqRecord objects (e.g. iteration
        or slicing to create a sub-alignment).
        """
        return self._records

    def __iter__(self) :
        """Iterate over alignment rows as SeqRecord objects.

        e.g.

        for record in align :
            print record.id
            print record.seq
        """
        return iter(self._records) 

    def get_seq_by_num(self, number):
        """Retrieve a sequence by row number (OBSOLETE).

        Returns:
        o A Seq object for the requested sequence.

        Raises:
        o IndexError - If the specified number is out of range.

        NOTE: This is a legacy method.  In new code where you need to access
        the rows of the alignment (i.e. the sequences) consider iterating
        over them or accessing them as SeqRecord objects.  e.g.

        for record in alignment :
            print record.id
            print record.seq
        first_record = alignment[0]
        last_record = alignment[-1]
        """
        return self._records[number].seq

    def get_alignment_length(self):
        """Return the maximum length of the alignment.

        All objects in the alignment should (hopefully) have the same
        length. This function will go through and find this length
        by finding the maximum length of sequences in the alignment.
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
        o descriptor - The descriptive id of the sequence being added.
                       This will be used as the resulting SeqRecord's
                       .id property (and, for historical compatibility,
                       also the .description property)
        o sequence - A string with sequence info.
        o start - You can explicitly set the start point of the sequence.
        This is useful (at least) for BLAST alignments, which can just
        be partial alignments of sequences.
        o end - Specify the end of the sequence, which is important
        for the same reason as the start.
        o weight - The weight to place on the sequence in the alignment.
        By default, all sequences have the same weight. (0.0 => no weight,
        1.0 => highest weight)
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
        """Returns a string containing a given column."""
        #TODO - Support negative indices?
        col_str = ''
        assert col >= 0 and col <= self.get_alignment_length()
        for rec in self._records:
            col_str += rec.seq[col]
        return col_str

    def __getitem__(self, index) :
        """Access part of the alignment.

        You can access a row of the alignment as a SeqRecord using an integer
        index (think of the alignment as a list of SeqRecord objects here):

            first_record = my_alignment[0]
            last_record = my_alignment[-1]

        You can also access use python's slice notation to create a sub-alignment
        containing only some of the SeqRecord objects:

            sub_alignment = my_alignment[2:20]

        This includes support for a step,

            sub_alignment = my_alignment[start:end:step]

        For example to select every second sequence:

            sub_alignment = my_alignment[::2]

        Or to reverse the row order:

            rev_alignment = my_alignment[::-1]

        Right now, these are the ONLY indexing operations supported.  The use of
        a second column based index is under discussion for a future update.
        """
        if isinstance(index, int) :
            #e.g. result = align[x]
            #Return a SeqRecord
            return self._records[index]
        elif isinstance(index, slice) :
            #e.g. sub_aling = align[i:j:k]
            #Return a new Alignment using only the specified records.
            #TODO - See Bug 2554 for changing the __init__ method
            #to allow us to do this more cleanly.
            sub_align = Alignment(self._alphabet)
            sub_align._records = self._records[index]
            return sub_align
        elif len(index)==2 :
            raise TypeError("Row and Column indexing is not currently supported,"\
                            +"but may be in future.")
        else :
            raise TypeError("Invalid index type.")

if __name__ == "__main__" :
    print "Mini self test..."

    raw_data = ["ACGATCAGCTAGCT", "CCGATCAGCTAGCT", "ACGATGAGCTAGCT"]
    a = Alignment(Alphabet.generic_dna)
    a.add_sequence("Alpha", raw_data[0], weight=2)
    a.add_sequence("Beta",  raw_data[1])
    a.add_sequence("Gamma", raw_data[2])

    print
    print "str(a):"
    print str(a)
    print
    print "repr(a):"
    print repr(a)
    print

    #Iterating over the rows...
    for rec in a :
        assert isinstance(rec, SeqRecord)
    for r,rec in enumerate(a) :
        assert isinstance(rec, SeqRecord)
        assert raw_data[r] == rec.seq.tostring()
        if r==0 : assert rec.annotations['weight']==2
    print "Alignment iteration as SeqRecord OK"

    print
    print "SeqRecord access by row:"
    print a[0].id, "...", a[-1].id

    print
    for format in ["fasta","phylip","clustal"] :
        print "="*60
        print "Using .format('%s')," % format
        print "="*60
        print a.format(format)

    print
    print "Row slicing the alignment:"
    print a[1:3]
    print
    print "Reversing the row order:"
    print a[::-1]
