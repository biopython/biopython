#!/usr/bin/env python
"""Example of generating a substitution matrix from an alignment.

Requires SubsMat.py and FreqTable.py from Iddo."""
# standard library
import sys

# Biopython
from Bio.SubsMat import SubsMat
from Bio.Clustalw import Clustalw
from Bio.Alphabet import IUPAC
from Bio.Align import AlignInfo

# get an alignment object from a Clustalw alignment output
c_align = Clustalw.parse_file('test.aln', IUPAC.unambiguous_dna)
summary_align = AlignInfo.SummaryInfo(c_align)

# get a replacement dictionary and accepted replacement matrix
replace_info = summary_align.replacement_dictionary(['N'])
my_arm = SubsMat.SeqMat(replace_info)

print replace_info

my_lom = SubsMat.make_log_odds_matrix(my_arm)

print 'log_odds_mat:', my_lom

my_lom.print_mat()

