# Copyright 2000 Brad Chapman.  All rights reserved.
# Revisions copyright 2007 Peter Cock.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Example using Bio.SeqIO to load a FASTA file as a dictionary.

An example function (get_accession_num) is defined to demonstrate
a non-trivial naming scheme where the dictionary key is based on
the record identifier.

The first version uses Bio.SeqIO.parse() and loads the entire
FASTA file into memory as a Python dictionary of SeqRecord
objects. This is *not* suitable for large files.

The second version used Bio.SeqIO.index() which is suitable
for FASTA files with millions of records.

See also Bio.SeqIO.index_db() and the examples in the main tutorial.
"""

from __future__ import print_function

from Bio.Alphabet import generic_dna
from Bio import SeqIO


def get_accession_num(seq_record):
    accession_atoms = seq_record.id.split('|')
    gb_name = accession_atoms[3]
    # strip the version info before returning
    return gb_name[:-2]


# In Memory
# =========
# This next bit of code uses Bio.SeqIO.parse() to load a FASTA file,
# and then turns it into an in-memory python dictionary.
# This is *not* suitable for FASTA files with millions of entries.


rec_iterator = SeqIO.parse("ls_orchid.fasta", "fasta", generic_dna)
orchid_dict = SeqIO.to_dict(rec_iterator, get_accession_num)

for id_num in orchid_dict:
    print('id number: %s' % id_num)
    print('description: %s' % orchid_dict[id_num].description)
    print('sequence: %s' % orchid_dict[id_num].seq)


# Indexed
# =======
# This next version uses the Bio.SeqIO.index() function which will index
# the FASTA file without loading all the records into memory at once.
# This is suitable for FASTA files with millions of entries.


orchid_dict = SeqIO.index("ls_orchid.fasta", "fasta", generic_dna)

for id_num in orchid_dict:
    print('id number: %s' % id_num)
    print('description: %s' % orchid_dict[id_num].description)
    print('sequence: %s' % orchid_dict[id_num].seq)
