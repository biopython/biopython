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
from Bio.Alphabet import IUPAC
from Bio import Clustalw
from Bio.Align.FormatConvert import FormatConverter
from Bio.Align import AlignInfo
from Bio.Fasta import FastaAlign
from Bio.SubsMat import FreqTable

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

print 'Calculating summary information...'
align_info = AlignInfo.SummaryInfo(alignment)
consensus = align_info.dumb_consensus()
print 'consensus:', consensus


print 'Replacement dictionary'
print align_info.replacement_dictionary(['N'])

print 'position specific score matrix.'
print 'with a supplied consensus sequence...'
print align_info.pos_specific_score_matrix(consensus, ['N'])

print 'defaulting to a consensus sequence...'
print align_info.pos_specific_score_matrix(chars_to_ignore = ['N'])

print 'with a selected sequence...'
second_seq = alignment.get_seq_by_num(1)
print align_info.pos_specific_score_matrix(second_seq, ['N'])

print 'information content'
print 'part of alignment:', align_info.information_content(5, 50,
                                chars_to_ignore = ['N'])
print 'entire alignment:', align_info.information_content(
                                chars_to_ignore = ['N'])

print 'relative information content'
e_freq = {'G' : 0.25,
          'C' : 0.25,
          'A' : 0.25,
          'T' : 0.25}

e_freq_table = FreqTable.FreqTable(e_freq, FreqTable.FREQ,
                                   IUPAC.unambiguous_dna)

print 'relative information:', align_info.information_content(
                                   e_freq_table = e_freq_table,
                                   chars_to_ignore = ['N'])

print 'Column 1:', align_info.get_column(1)
print 'IC for column 1:', align_info.ic_vector[1]
print 'Column 7:', align_info.get_column(7)
print 'IC for column 7:', align_info.ic_vector[7]
print 'test print_info_content'
AlignInfo.print_info_content(align_info)
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
align_info = AlignInfo.SummaryInfo(alignment)
consensus = align_info.dumb_consensus(ambiguous = "X")
print 'consensus:', consensus

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



