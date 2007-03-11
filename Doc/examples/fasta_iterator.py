# The New Way
# ===========
# This next bit of code use Bio.SeqIO to parse a FASTA file

from Bio import SeqIO

def extract_organisms(file_to_parse, format):
    all_species = []
    for cur_record in SeqIO.parse(open(file_to_parse), format) :
        # extract the info from the description
        new_species = cur_record.description.split()[1]

        # append the new species to the list if it isn't there
        if new_species not in all_species:
            all_species.append(new_species)

    return all_species

if __name__ == "__main__":
    print "Using Bio.SeqIO on a FASTA file"
    all_species = extract_organisms("ls_orchid.fasta", "fasta")
    print "number of species:", len(all_species)
    print 'species names:', all_species


# The Old Way
# ===========
# This next bit of code still works fine, it uses Bio.Fasta instead

from Bio import Fasta

def extract_organisms(file_to_parse):
    # set up the parser and iterator
    parser = Fasta.RecordParser()
    file = open(file_to_parse, 'r')
    iterator = Fasta.Iterator(file, parser)

    all_species = []

    while 1:
        cur_record = iterator.next()

        if cur_record is None:
            break
        
        # extract the info from the title
        new_species = cur_record.title.split()[1]

        # append the new species to the list if it isn't there
        if new_species not in all_species:
            all_species.append(new_species)

    return all_species

if __name__ == "__main__":
    print "Using Bio.Fasta"
    all_species = extract_organisms("ls_orchid.fasta")
    print "number of species:", len(all_species)
    print 'species names:', all_species
    
