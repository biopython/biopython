#!/usr/bin/env python
"""clustal_run.py

Example code to show how to create a clustalw command line, run clustalw
and parse the results into an object that can be dealt with easily."""
# standard library
import os

# biopython
from Bio.Alphabet import IUPAC
from Bio import Clustalw
from Bio.Clustalw import MultipleAlignCL
from Bio.Align import AlignInfo
from Bio.SubsMat import FreqTable

# create the command line to run clustalw
# this assumes you've got clustalw somewhere on your path, otherwise
# you need to pass a second argument to MultipleAlignCL with the complete
# path to clustalw
cline = MultipleAlignCL(os.path.join(os.curdir, 'opuntia.fasta'))
cline.set_output('test.aln')

# actually perform the alignment and get back an alignment object
alignment = Clustalw.do_alignment(cline)

# get the records in the alignment
all_records = alignment.get_all_seqs()

print 'description:', all_records[0].description
print 'sequence:', all_records[0].seq

# get the length of the alignment
print 'length', alignment.get_alignment_length()

print alignment

# print out interesting information about the alignment
summary_align = AlignInfo.SummaryInfo(alignment)

consensus = summary_align.dumb_consensus()
print 'consensus', consensus

my_pssm = summary_align.pos_specific_score_matrix(consensus,
                                                  chars_to_ignore = ['N'])

print my_pssm

expect_freq = {
    'A' : .3,
    'G' : .2,
    'T' : .3,
    'C' : .2}

freq_table_info = FreqTable.FreqTable(expect_freq, FreqTable.FREQ,
                                      IUPAC.unambiguous_dna)

info_content = summary_align.information_content(5, 30,
                                                 chars_to_ignore = ['N'],
                                                 e_freq_table = \
                                                 freq_table_info)

print "relative info content:", info_content

