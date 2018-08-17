#!/usr/bin/env python
# Copyright 2000 Brad Chapman.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Example of generating a substitution matrix from an alignment."""

from __future__ import print_function

# Biopython
from Bio import SubsMat
from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
from Bio.Align import AlignInfo

# get an alignment object from a Clustalw alignment output
c_align = AlignIO.read('protein.aln', 'clustal',
                       alphabet=Gapped(IUPAC.protein))
summary_align = AlignInfo.SummaryInfo(c_align)

# get a replacement dictionary and accepted replacement matrix
# exclude all amino acids that aren't charged polar
replace_info = summary_align.replacement_dictionary(["G", "A", "V", "L", "I",
                                                     "M", "P", "F", "W", "S",
                                                     "T", "N", "Q", "Y", "C"])

my_arm = SubsMat.SeqMat(replace_info)

print(replace_info)

my_lom = SubsMat.make_log_odds_matrix(my_arm)

print('log_odds_mat: %s' % my_lom)

my_lom.print_mat()
