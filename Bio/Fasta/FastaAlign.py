"""
Code to deal with alignments written in Fasta format (DEPRECATED).

This module is considered obsolete and has been deprecated. It will be
removed in a future release of Biopython. Please use Bio.AlignIO instead
for reading and writing alignments in FASTA format.

This mostly just uses the regular Fasta parsing stuff written by Jeff
to deal with all of the input and output formats.

functions:
o parse_file()

classes:
FastaAlignment"""
# standard library
import os

# biopython
from Bio.Align.Generic import Alignment
from Bio import Alphabet
from Bio.Alphabet import IUPAC
from Bio import Fasta

def parse_file(file_name, type = 'DNA'):
    """Parse the given file into a FastaAlignment object.

    Arguments:
    o file_name - The location of the file to parse.
    o type - The type of information contained in the file.
    """
    if type.upper() == 'DNA':
        alphabet = IUPAC.ambiguous_dna
    elif type.upper() == 'RNA':
        alphabet = IUPAC.ambiguous_rna
    elif type.upper() == 'PROTEIN':
        alphabet = IUPAC.protein
    else:
        raise ValueError("Invalid type %s passed. Need DNA, RNA or PROTEIN"
                         % type)

    # create a new alignment object
    fasta_align = FastaAlignment(Alphabet.Gapped(alphabet))

    # now parse the file and fill up the alignment object
    align_file = open(file_name, 'r')

    parser = Fasta.RecordParser()
    iterator = Fasta.Iterator(align_file, parser)

    cur_align = iterator.next()
    while cur_align:
        fasta_align.add_sequence(cur_align.title, cur_align.sequence)

        cur_align = iterator.next()

    return fasta_align

class FastaAlignment(Alignment):
    """Work with the Fasta Alignment format.

    The fasta alignment format is basically the same as the regular ol'
    Fasta format we know and love, except the sequences have gaps
    (represented by -'s).
    """
    def __init__(self, alphabet = Alphabet.Gapped(IUPAC.ambiguous_dna)):
        Alignment.__init__(self, alphabet)

    def __str__(self):
        """Print out a fasta version of the alignment info."""
        return_string = ''
        for item in self._records:
            new_f_record = Fasta.Record()
            new_f_record.title = item.description
            new_f_record.sequence = item.seq.data

            return_string = return_string + str(new_f_record) + os.linesep + os.linesep

        # have a extra newline, so strip two off and add one before returning
        return return_string.rstrip() + os.linesep

            
            
        

