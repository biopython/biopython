# Copyright Yair Benita Y.Benita@pharm.uu.nl
# Biopython (http://biopython.org) license applies

"""Codon adaption indxes, including Sharp and Li (1987) E. coli index.

Currently this module only defines a single codon adaption index from
Sharp & Li, Nucleic Acids Res. 1987.
"""

__docformat__ = "restructuredtext en"


SharpEcoliIndex = {
'GCA': 0.586, 'GCC': 0.122, 'GCG': 0.424, 'GCT': 1, 'AGA': 0.004, 'AGG': 0.002, 'CGA': 0.004,
'CGC': 0.356, 'CGG': 0.004, 'CGT': 1, 'AAC': 1, 'AAT': 0.051, 'GAC': 1, 'GAT': 0.434, 'TGC': 1,
'TGT': 0.5, 'CAA': 0.124, 'CAG': 1, 'GAA': 1, 'GAG': 0.259, 'GGA': 0.01, 'GGC': 0.724, 'GGG': 0.019,
'GGT': 1, 'CAC': 1, 'CAT': 0.291, 'ATA': 0.003, 'ATC': 1, 'ATT': 0.185, 'CTA': 0.007, 'CTC': 0.037,
'CTG': 1, 'CTT': 0.042, 'TTA': 0.02, 'TTG': 0.02, 'AAA': 1, 'AAG': 0.253, 'ATG': 1, 'TTC': 1, 'TTT': 0.296,
'CCA': 0.135, 'CCC': 0.012, 'CCG': 1, 'CCT': 0.07, 'AGC': 0.41, 'AGT': 0.085, 'TCA': 0.077, 'TCC': 0.744,
'TCG': 0.017, 'TCT': 1, 'ACA': 0.076, 'ACC': 1, 'ACG': 0.099, 'ACT': 0.965, 'TGG': 1, 'TAC': 1, 'TAT': 0.239,
'GTA': 0.495, 'GTC': 0.066, 'GTG': 0.221, 'GTT': 1}
