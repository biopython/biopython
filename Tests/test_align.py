#!/usr/bin/env python
"""test_align.py

A script to test alignment stuff.

Right now we've got tests for:
o Reading and Writing clustal format
o Reading and Writing fasta format
o Converting between formats"""

# standard library
import os 

# biopython
from Bio.Clustalw import Clustalw
from Bio.Align.FormatConvert import FormatConverter
from Bio.Fasta import FastaAlign

print "testing reading and writing clustal format..."
test_dir = os.path.join(os.getcwd(), 'Clustalw')
test_names = ['opuntia.aln', 'cw02.aln']

test_files = []
for name in test_names:
    test_files.append(os.path.join(test_dir, name))

for test_file in test_files:
    # parse the alignment file and get an aligment object
    alignment = Clustalw.parse_file(test_file)

    # print the alignment back out
    print alignment

alignment = Clustalw.parse_file(os.path.join(test_dir, test_names[0]))

# test the base alignment stuff
print 'all_seqs...'
all_seqs = alignment.get_all_seqs()
for seq_record in all_seqs:
    print 'description:', seq_record.description
    print 'seq:', seq_record.seq
print 'length:', alignment.get_alignment_length()
print 'consensus:', alignment.dumb_consensus()


print "testing reading and writing fasta format..."

to_parse = os.path.join(os.curdir, 'Fasta', 'fa01')

alignment = FastaAlign.parse_file(to_parse, 'PROTEIN')

# test the base alignment stuff
print 'all_seqs...'
all_seqs = alignment.get_all_seqs()
for seq_record in all_seqs:
    print 'description:', seq_record.description
    print 'seq:', seq_record.seq

print 'length:', alignment.get_alignment_length()
print 'consensus:', alignment.dumb_consensus(ambiguous = "X")

print alignment


print "Test format conversion..."

# parse the alignment file and get an aligment object
alignment = Clustalw.parse_file(os.path.join(os.curdir, 'Clustalw',
                                             'opuntia.aln'))

converter = FormatConverter(alignment)

fasta_align = converter.to_fasta()
clustal_align = converter.to_clustal()

print fasta_align
print clustal_align



