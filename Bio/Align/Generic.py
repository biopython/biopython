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

    This is a base class to represent alignments, which should be subclassed
    to deal with an alignment in a specific format.
    """
    def __init__(self, alphabet):
        """Initialize a new Alignment object.

        Arguments:
        o alphabet - The alphabet to use for the sequence objects that are
        created. This alphabet must be a gapped type.
        """
        self._alphabet = alphabet
        # hold everything at a list of seq record objects
        self._records = []

    def get_all_seqs(self):
        """Return all of the sequences involved in the alignment.

        The return value is a list of SeqRecord objects.
        """
        return self._records

    def get_seq_by_num(self, number):
        """Retrieve a sequence by the number of the sequence in the consensus.

        Returns:
        o A Seq object for the requested sequence.

        Raises:
        o IndexError - If the specified number is out of range.
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
        new_record = SeqRecord(new_seq, description = descriptor)

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
        """Returns a string containing a given column"""
        col_str = ''
        assert col >= 0 and col <= self.get_alignment_length()
        for rec in self._records:
            col_str += rec.seq[col]
        return col_str
                
        
        

