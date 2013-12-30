# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#

"""Example using Bio.SeqIO to parse a FASTA file."""

from __future__ import print_function

from Bio import SeqIO


def extract_organisms(file_to_parse, format):
    all_species = []
    for cur_record in SeqIO.parse(open(file_to_parse), format):
        # extract the info from the description
        new_species = cur_record.description.split()[1]

        # append the new species to the list if it isn't there
        if new_species not in all_species:
            all_species.append(new_species)

    return all_species

if __name__ == "__main__":
    print("Using Bio.SeqIO on a FASTA file")
    all_species = extract_organisms("ls_orchid.fasta", "fasta")
    print("number of species: %i" % len(all_species))
    print("species names: %s" % all_species)
