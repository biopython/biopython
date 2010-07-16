# In Memory
# =========
# This next bit of code uses Bio.SeqIO.parse() to load a FASTA file,
# and then turns it into an in-memory python dictionary.
# This is *not* suitable for FASTA files with millions of entries.

from Bio.Alphabet import generic_dna
from Bio import SeqIO

def get_accession_num(seq_record):
    accession_atoms = seq_record.id.split('|')
    gb_name = accession_atoms[3]
    # strip the version info before returning
    return gb_name[:-2]

rec_iterator = SeqIO.parse("ls_orchid.fasta","fasta", generic_dna)
orchid_dict = SeqIO.to_dict(rec_iterator, get_accession_num)

for id_num in orchid_dict:
    print 'id number:', id_num
    print 'description:', orchid_dict[id_num].description
    print 'sequence:', orchid_dict[id_num].seq


# Indexed
# =======
# This next version uses the Bio.SeqIO.index() function which will index
# the FASTA file without loading all the records into memory at once.
# This is suitable for FASTA files with millions of entries.

from Bio.Alphabet import generic_dna
from Bio import SeqIO

def get_accession_num(record_id):
    accession_atoms = record_id.split('|')
    gb_name = accession_atoms[3]
    # strip the version info before returning
    return gb_name[:-2]

orchid_dict = SeqIO.index("ls_orchid.fasta","fasta", generic_dna)

for id_num in orchid_dict:
    print 'id number:', id_num
    print 'description:', orchid_dict[id_num].description
    print 'sequence:', orchid_dict[id_num].seq
