import string
from Bio import Fasta
from Bio.Alphabet import IUPAC

def get_accession_num(fasta_record):
    title_atoms = string.split(fasta_record.title)

    accession_atoms = string.split(title_atoms[0], '|')

    gb_name = accession_atoms[3]

    # strip the version info before returning
    return gb_name[:-2]


Fasta.index_file("ls_orchid.fasta", "my_orchid_dict.idx",
                 get_accession_num)

dna_parser = Fasta.SequenceParser(IUPAC.ambiguous_dna)

orchid_dict = Fasta.Dictionary("my_orchid_dict.idx", dna_parser)

for id_num in orchid_dict.keys():
    print 'id number:', id_num
    print 'description:', orchid_dict[id_num].description
    print 'sequence:', orchid_dict[id_num].seq
