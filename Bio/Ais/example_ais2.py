"""Example of Ais.
Gene to check for mutations is XM_041014, Homo sapiens tryptophanyl-tRNA 
synthetase.  The sequence is broken into 60 residue segments.
A variant is found in SNP as rs1373 and tested against the base set..
"""

import os
import sys
import string
from urllib import FancyURLopener
from urllib import urlencode

from Bio.Seq import Seq
from Bio.Align.Generic import Alignment
from Bio.Align.AlignInfo import SummaryInfo
from Bio.Alphabet import DNAAlphabet
from Bio.Alphabet import Gapped
from Bio.Ais import match_sequence
from Bio.Ais import Lymphocyte
from Bio.Ais import Immune






align = Alignment( Gapped( DNAAlphabet() ) )
alphabet = 'agct'
fh = open( 'XM041014.txt' )
while 1:
    name = fh.readline()
    if name == '': break
    name = name.strip()
    seq = fh.readline()
    if seq == '': break
    seq = seq.strip()
    align.add_sequence( name, seq )
fh.close()
immune = Immune( align, alphabet )
mystery_seqs = [ 'CTTCTTCCCAGCCCTGCAGGGCGCCCAGACCAAAATGAGTGCCAGCGACCCCAACTCCTC',\
                 'CATCTTCCTCACCGACACGGCCAAGCAGATCAAAACCAAGGTGAGCAGTGGCCCAGGCAG', \
                 'AGGGCAGCCAGCGGGATCAGAAACATGCCCTGGGGGGCACACGCACACACACACCACTGC', \
                 'GCAGTGCTTTGGGCCTCAGCCACCCAGCCACCGTAGCAATGCTGAACTCCAGGATTGCCA', \
                 'ACTTTGCTCCCTGCAGAACTCAAGGCAGCCTGTCCTTCCCACATGAATCAGGTGAATTTG', \
                 'TTTTGTGTCCTTTCATTTAGGTAATATTTTGAATTTTTTTTCTTTTTTTCCTTTTCTGCT', \
                 'TTTCTTATCAGCACTGTTCTGCAGGCAAAGTCATGGCAATTTCCTTAAAGTTATAATTAG', \
                 'TGTGGATGTTATTAAATAATGGAAAAGATTTGTGTATATTTCTCATCATCAGCTGTGTTC', \
                 'CAGCGGCTGGACCCGGTGCCTTATTTATGTTTTCTCATGTTTCATTCTCACTGGGAGGTG', \
                 'GTGGTGCTATTCCCTTTAAAGATGAGATGATGAGGCTTAGAGAAGCCCAGCCACTTGACC', \
                 'AAGTTTCCACCACCAGGAATCAGCAGCACTGCTCTCCTGGCCCTGAGCAGCCATTACCAA', \
                 'GCAGGGGTTATGAGCTGGTCTTCTCCCAGAGAGCGTTCTCACTGAGAGACATAGAACGGA'

     ]
for mystery_seq in mystery_seqs:
    for j in range( 0, 200 ):
#    for j in range( 0, 1 ):
        suspicious = immune.random_test( mystery_seq  )
        if suspicious:
            print '%s triggered immune system' % mystery_seq.lower()




