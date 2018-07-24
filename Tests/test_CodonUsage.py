# Copyright 2003 by Iddo Friedberg.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from __future__ import print_function

from Bio.SeqUtils import CodonUsage
import os
import sys

# first make a CAI object
X = CodonUsage.CodonAdaptationIndex()
# now generate an index from a file
if os.path.exists("./CodonUsage/HighlyExpressedGenes.txt"):
    X.generate_index("./CodonUsage/HighlyExpressedGenes.txt")
elif os.path.exists("./Tests/CodonUsage/HighlyExpressedGenes.txt"):
    X.generate_index("./Tests/CodonUsage/HighlyExpressedGenes.txt")
else:
    print("Cannot find the file HighlyExpressedGene.txt\nMake sure you run the tests from within the Tests folder")
    sys.exit()
# alternatively you could use any predefined dictionary like this:
# from CaiIndices import SharpIndex # you can save your dictionary in this file.
# X.SetCaiIndex(SharpIndex)

print("The current index used:")
X.print_index()

print("-" * 60)
print("codon adaptation index for test gene: %.2f" % X.cai_for_gene("ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGA"))



""" nRCA test
    Creates a normRelativeCodonAdaptationIndex object, generates an index from
    a reference set of genes and scores a set of genes (inputfile)
    provided all in FASTA format
"""
# Copyright Ivan Erill 2018 (erill@umbc.edu)
# -*- coding: utf-8 -*-
from Bio import SeqIO  # To parse a FASTA file 

# first make a nRCA object
nrca_index = CodonUsage.normRelativeCodonAdaptationIndex()
# now generate an index from a file
if os.path.exists("./CodonUsage/refset_Sharp.fas"):
    nrca_index.generate_index("./CodonUsage/refset_Sharp.fas")
elif os.path.exists("./Tests/CodonUsage/refset_Sharp.fas"):
    nrca_index.generate_index("./Tests/CodonUsage/refset_Sharp.fas")
else:
    print("Cannot find the file refset_Sharp.fas\n \
          Make sure you run the tests from within the Tests folder")
    sys.exit()
# alternatively you could use any predefined dictionary like this:
# from CaiIndices import SharpIndex # you can save your dictionary in this file.
# X.SetCaiIndex(SharpIndex)

infilename = 'sample_genes.fas'
outfilename = 'nRCA_test_output.csv'
  
#open input and output csv
outcsvfile = open(outfilename, 'w')
outcsvfile.write('SeqID, nRCA\n')

with open(infilename, 'r') as handle: 
    # compute nRCA for sequence in input file and write to output
    for cur_record in SeqIO.parse(handle, "fasta"):
        nrca_val = nrca_index.nrca_for_gene(str(cur_record.seq).upper())
        outcsvfile.write(cur_record.id + ',' + str(nrca_val) + '\n')



