from Bio.ParserSupport import AbstractConsumer
from Bio import Fasta
from Bio.File import UndoHandle

import string

class SpeciesExtractor(AbstractConsumer):

    def __init__(self):
        self.species_list = []

    def title(self, title_info):
        title_atoms = string.split(title_info)
        new_species = title_atoms[1]

        if new_species not in self.species_list:
            self.species_list.append(new_species)


def extract_organisms(file, num_records):
    scanner = Fasta._Scanner()
    consumer = SpeciesExtractor()

    file_to_parse = UndoHandle(open(file, 'r'))

    for fasta_record in range(num_records):
        scanner.feed(file_to_parse, consumer)

    file_to_parse.close()

    return consumer.species_list


if __name__ == "__main__":
    all_species = extract_organisms("ls_orchid.fasta", 94)
    print "number of species:", len(all_species)
    print 'species names:', all_species


