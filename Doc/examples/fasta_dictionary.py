# The New Way
# ===========
# This next bit of code use Bio.SeqIO to parse a FASTA file, and turn
# it into an in-memory python dictionary.
# Note Bio.SeqIO will use a generic alphabet for the SeqRecord objects.

from Bio import SeqIO

def get_accession_num(seq_record):
    accession_atoms = seq_record.id.split('|')
    gb_name = accession_atoms[3]
    # strip the version info before returning
    return gb_name[:-2]

orchid_dict = SeqIO.to_dict(SeqIO.parse(open("ls_orchid.fasta"),"fasta"), \
                            get_accession_num)

for id_num in orchid_dict.keys():
    print 'id number:', id_num
    print 'description:', orchid_dict[id_num].description
    print 'sequence:', orchid_dict[id_num].seq

# The Old Way
# ===========
# This next bit of code still works fine. It uses Bio.Fasta instead,
# and builds an index as a set of files on disc in the sub-directory
# my_orchid_dict.idx
# Note that the alphabet is explicitly defined for the sequences.

import os
from Bio import Fasta
from Bio.Alphabet import IUPAC

def get_accession_num(fasta_record):
    title_atoms = fasta_record.title.split()
    accession_atoms = title_atoms[0].split('|')
    gb_name = accession_atoms[3]
    # strip the version info before returning
    return gb_name[:-2]

if not os.path.isdir("my_orchid_dict.idx") :
    #Build a new index
    Fasta.index_file("ls_orchid.fasta", "my_orchid_dict.idx",
                     get_accession_num)
else :
    print "Reusing existing index"

dna_parser = Fasta.SequenceParser(IUPAC.ambiguous_dna)

orchid_dict = Fasta.Dictionary("my_orchid_dict.idx", dna_parser)

for id_num in orchid_dict.keys():
    print 'id number:', id_num
    print 'description:', orchid_dict[id_num].description
    print 'sequence:', orchid_dict[id_num].seq
