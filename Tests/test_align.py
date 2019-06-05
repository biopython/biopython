#!/usr/bin/env python
# Copyright 2000-2001 by Brad Chapman.  All rights reserved.
# Revisions copyright 2007-2003 by Peter Cock. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Test alignment stuff.

Right now we've got tests for:

- Reading and Writing clustal format
- Reading and Writing fasta format
- Converting between formats

"""

# standard library
from __future__ import print_function

import os
import unittest

# biopython
from Bio import Alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Align import AlignInfo
from Bio import AlignIO
from Bio.SubsMat import FreqTable
from Bio.Align import MultipleSeqAlignment

# Very simple tests on an empty alignment
alignment = MultipleSeqAlignment([], Alphabet.generic_alphabet)
assert alignment.get_alignment_length() == 0
assert len(alignment) == 0
del alignment

# Basic tests on simple three string alignment
alignment = MultipleSeqAlignment([], Alphabet.generic_alphabet)
letters = "AbcDefGhiJklMnoPqrStuVwxYz"
alignment.append(SeqRecord(Seq(letters), id="mixed"))
alignment.append(SeqRecord(Seq(letters.lower()), id="lower"))
alignment.append(SeqRecord(Seq(letters.upper()), id="upper"))
assert alignment.get_alignment_length() == 26
assert len(alignment) == 3
assert str(alignment[0].seq) == letters
assert str(alignment[1].seq) == letters.lower()
assert str(alignment[2].seq) == letters.upper()
assert alignment[0].id == "mixed"
assert alignment[1].id == "lower"
assert alignment[2].id == "upper"
for (col, letter) in enumerate(letters):
    assert alignment[:, col] == letter + letter.lower() + letter.upper()
# Check row extractions:
assert alignment[0].id == "mixed"
assert alignment[-1].id == "upper"
# Check sub-alignment extraction by row slicing:
assert isinstance(alignment[::-1], MultipleSeqAlignment)
assert alignment[::-1][0].id == "upper"
assert alignment[::-1][2].id == "mixed"

del alignment
del letters

print("testing reading and writing clustal format...")
test_dir = os.path.join(os.getcwd(), 'Clustalw')
test_names = ['opuntia.aln', 'cw02.aln']

test_files = []
for name in test_names:
    test_files.append(os.path.join(test_dir, name))

test_file = test_files[0]
# parse the alignment file and get an aligment object
alignment = AlignIO.read(test_file, "clustal")
# print the alignment back out
print(alignment.format("clustal"))

test_file = test_files[1]
# parse the alignment file and get an aligment object
alignment = AlignIO.read(test_file, "clustal")
# print the alignment back out
print(alignment.format("clustal"))

alignment = AlignIO.read(os.path.join(test_dir, test_names[0]), "clustal",
                         alphabet=Alphabet.Gapped(IUPAC.unambiguous_dna))

# test the base alignment stuff
print('all_seqs...')
assert len(alignment) == 7
seq_record = alignment[0]
print('description: %s' % seq_record.description)
print('seq: %r' % seq_record.seq)
seq_record = alignment[1]
print('description: %s' % seq_record.description)
print('seq: %r' % seq_record.seq)
seq_record = alignment[2]
print('description: %s' % seq_record.description)
print('seq: %r' % seq_record.seq)
seq_record = alignment[3]
print('description: %s' % seq_record.description)
print('seq: %r' % seq_record.seq)
seq_record = alignment[4]
print('description: %s' % seq_record.description)
print('seq: %r' % seq_record.seq)
seq_record = alignment[5]
print('description: %s' % seq_record.description)
print('seq: %r' % seq_record.seq)
seq_record = alignment[6]
print('description: %s' % seq_record.description)
print('seq: %r' % seq_record.seq)


print('length: %i' % alignment.get_alignment_length())

print('Calculating summary information...')
align_info = AlignInfo.SummaryInfo(alignment)
consensus = align_info.dumb_consensus()
assert isinstance(consensus, Seq)
print('consensus: %r' % consensus)


print('Replacement dictionary')
ks = sorted(align_info.replacement_dictionary(['N']))
assert len(ks) == 16
key = ('A', 'A')
print("%s : %s" % (key, align_info.replacement_dictionary(['N'])[key]))
key = ('A', 'C')
print("%s : %s" % (key, align_info.replacement_dictionary(['N'])[key]))
key = ('A', 'G')
print("%s : %s" % (key, align_info.replacement_dictionary(['N'])[key]))
key = ('A', 'T')
print("%s : %s" % (key, align_info.replacement_dictionary(['N'])[key]))
key = ('C', 'A')
print("%s : %s" % (key, align_info.replacement_dictionary(['N'])[key]))
key = ('C', 'C')
print("%s : %s" % (key, align_info.replacement_dictionary(['N'])[key]))
key = ('C', 'G')
print("%s : %s" % (key, align_info.replacement_dictionary(['N'])[key]))
key = ('C', 'T')
print("%s : %s" % (key, align_info.replacement_dictionary(['N'])[key]))
key = ('G', 'A')
print("%s : %s" % (key, align_info.replacement_dictionary(['N'])[key]))
key = ('G', 'C')
print("%s : %s" % (key, align_info.replacement_dictionary(['N'])[key]))
key = ('G', 'G')
print("%s : %s" % (key, align_info.replacement_dictionary(['N'])[key]))
key = ('G', 'T')
print("%s : %s" % (key, align_info.replacement_dictionary(['N'])[key]))
key = ('T', 'A')
print("%s : %s" % (key, align_info.replacement_dictionary(['N'])[key]))
key = ('T', 'C')
print("%s : %s" % (key, align_info.replacement_dictionary(['N'])[key]))
key = ('T', 'G')
print("%s : %s" % (key, align_info.replacement_dictionary(['N'])[key]))
key = ('T', 'T')
print("%s : %s" % (key, align_info.replacement_dictionary(['N'])[key]))

print('position specific score matrix.')
print('with a supplied consensus sequence...')
print(align_info.pos_specific_score_matrix(consensus, ['N']))

print('defaulting to a consensus sequence...')
print(align_info.pos_specific_score_matrix(chars_to_ignore=['N']))

print('with a selected sequence...')
second_seq = alignment[1].seq
print(align_info.pos_specific_score_matrix(second_seq, ['N']))

print('information content')
print('part of alignment: %0.2f'
      % align_info.information_content(5, 50, chars_to_ignore=['N']))
print('entire alignment: %0.2f'
      % align_info.information_content(chars_to_ignore=['N']))

print('relative information content')
e_freq = {'G': 0.25,
          'C': 0.25,
          'A': 0.25,
          'T': 0.25}

e_freq_table = FreqTable.FreqTable(e_freq, FreqTable.FREQ,
                                   IUPAC.unambiguous_dna)

print('relative information: %0.2f'
      % align_info.information_content(e_freq_table=e_freq_table,
                                       chars_to_ignore=['N']))

print('Column 1: %s' % align_info.get_column(1))
print('IC for column 1: %0.2f' % align_info.ic_vector[1])
print('Column 7: %s' % align_info.get_column(7))
print('IC for column 7: %0.2f' % align_info.ic_vector[7])
print('test print_info_content')
AlignInfo.print_info_content(align_info)
print("testing reading and writing fasta format...")

to_parse = os.path.join(os.curdir, 'Quality', 'example.fasta')

alignment = AlignIO.read(to_parse, "fasta",
                         alphabet=Alphabet.Gapped(IUPAC.ambiguous_dna))

# test the base alignment stuff
print('all_seqs...')
assert len(alignment) == 3
seq_record = alignment[0]
print('description: %s' % seq_record.description)
print('seq: %r' % seq_record.seq)
seq_record = alignment[1]
print('description: %s' % seq_record.description)
print('seq: %r' % seq_record.seq)
seq_record = alignment[2]
print('description: %s' % seq_record.description)
print('seq: %r' % seq_record.seq)

print('length: %i' % alignment.get_alignment_length())
align_info = AlignInfo.SummaryInfo(alignment)
consensus = align_info.dumb_consensus(ambiguous="N", threshold=0.6)
assert isinstance(consensus, Seq)
print('consensus: %r' % consensus)

print(alignment)


print("Test format conversion...")

# parse the alignment file and get an aligment object
alignment = AlignIO.read(os.path.join(os.curdir, 'Clustalw', 'opuntia.aln'),
                         'clustal')

print("As FASTA:")
print(alignment.format("fasta"))
print("As Clustal:")
print(alignment.format("clustal"))


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
