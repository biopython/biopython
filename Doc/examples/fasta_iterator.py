import string
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
        title_atoms = string.split(cur_record.title)
        new_species = title_atoms[1]

        # append the new species to the list if it isn't there
        if new_species not in all_species:
            all_species.append(new_species)

    return all_species

if __name__ == "__main__":
    all_species = extract_organisms("ls_orchid.fasta")
    print "number of species:", len(all_species)
    print 'species names:', all_species
    
