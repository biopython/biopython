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
from Bio import Alphabet
from Bio import Seq
from Bio.Alphabet import IUPAC
from Bio import Clustalw
from Bio.Align.FormatConvert import FormatConverter
from Bio.Align import AlignInfo
from Bio.Fasta import FastaAlign
from Bio.SubsMat import FreqTable
from Bio.Align.Generic import Alignment

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
ks = align_info.replacement_dictionary(['N']).keys()
ks.sort()
for key in ks:
    print "%s : %s" % (key, align_info.replacement_dictionary(['N'])[key])

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

"""
# test to find a position in an original sequence given a
# column position in an alignment
print "Testing finding column positions..."
alignment_info = ["GATC--CGATC--G",
                  "GA--CCCG-TC--G",
                  "GAT--CC--TC--G"]

gapped_unambiguous = Alphabet.Gapped(IUPAC.unambiguous_dna)

alignment = Alignment(gapped_unambiguous)
for seq in alignment_info:
    alignment.add_sequence("Blah", seq)

test_seq_1 = Seq.Seq("GATCCGATCG")
orig_pos = alignment.original_sequence_pos(3, test_seq_1, 0)
assert orig_pos == 3, "Got unexpected position: %s" % orig_pos
orig_pos = alignment.original_sequence_pos(7, test_seq_1, 0)
assert orig_pos == 5, "Got unexpected position: %s" % orig_pos
orig_pos = alignment.original_sequence_pos(0, test_seq_1, 0)
assert orig_pos == 0, "Got unexpected position: %s" % orig_pos
orig_pos = alignment.original_sequence_pos(13, test_seq_1, 0)
assert orig_pos == 9, "Got unexpected position: %s" % orig_pos

try:
    orig_pos = alignment.original_sequence_pos(5, test_seq_1, 0)
    raise AssertionError("Did not fail with a junk position")
except AssertionError:
    pass
"""
