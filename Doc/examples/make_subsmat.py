#!/usr/bin/env python
"""Example of generating a substitution matrix from an alignment.
"""
# standard library
import sys

# Biopython
from Bio import SubsMat
from Bio import Clustalw
from Bio.Alphabet import IUPAC
from Bio.Align import AlignInfo

# get an alignment object from a Clustalw alignment output
c_align = Clustalw.parse_file('protein.aln', IUPAC.protein)
summary_align = AlignInfo.SummaryInfo(c_align)

# get a replacement dictionary and accepted replacement matrix
# exclude all amino acids that aren't charged polar
replace_info = summary_align.replacement_dictionary(["G", "A", "V", "L", "I",
                                                     "M", "P", "F", "W", "S",
                                                     "T", "N", "Q", "Y", "C"])

my_arm = SubsMat.SeqMat(replace_info)

print replace_info

my_lom = SubsMat.make_log_odds_matrix(my_arm)

print 'log_odds_mat:', my_lom

my_lom.print_mat()

