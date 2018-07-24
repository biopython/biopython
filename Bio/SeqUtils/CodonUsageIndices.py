"""Codon adaption indices, including Sharp and Li (1987) E. coli index.

This module defines the single codon adaption index from Sharp & Li, 
Nucleic Acids Res. 1987 and the Normalized Relative Codon Adaptation index 
values for Escherichia coli K-12, from O'Neill, Or and Erill 
PLoS ONE. 2013 Oct 7;8(10):e76177.
"""

# Copyright Yair Benita Y.Benita@pharm.uu.nl
# Biopython (http://biopython.org) license applies
# Sharp & Li CAI for E. coli K-12
SharpEcoliIndex = {
    'GCA': 0.586, 'GCC': 0.122, 'GCG': 0.424, 'GCT': 1, 'AGA': 0.004,
    'AGG': 0.002, 'CGA': 0.004, 'CGC': 0.356, 'CGG': 0.004, 'CGT': 1, 'AAC': 1,
    'AAT': 0.051, 'GAC': 1, 'GAT': 0.434, 'TGC': 1, 'TGT': 0.5, 'CAA': 0.124,
    'CAG': 1, 'GAA': 1, 'GAG': 0.259, 'GGA': 0.01, 'GGC': 0.724, 'GGG': 0.019,
    'GGT': 1, 'CAC': 1, 'CAT': 0.291, 'ATA': 0.003, 'ATC': 1, 'ATT': 0.185,
    'CTA': 0.007, 'CTC': 0.037, 'CTG': 1, 'CTT': 0.042, 'TTA': 0.02,
    'TTG': 0.02, 'AAA': 1, 'AAG': 0.253, 'ATG': 1, 'TTC': 1, 'TTT': 0.296,
    'CCA': 0.135, 'CCC': 0.012, 'CCG': 1, 'CCT': 0.07, 'AGC': 0.41,
    'AGT': 0.085, 'TCA': 0.077, 'TCC': 0.744, 'TCG': 0.017, 'TCT': 1,
    'ACA': 0.076, 'ACC': 1, 'ACG': 0.099, 'ACT': 0.965, 'TGG': 1, 'TAC': 1,
    'TAT': 0.239, 'GTA': 0.495, 'GTC': 0.066, 'GTG': 0.221, 'GTT': 1}

# Copyright Ivan Erill (erill@umbc.edu)
# Biopython (http://biopython.org) license applies 
# ErillLab normalized Relative Codon Adaptation for E. coli K-12
ErillLabEcoliIndex = { 
        'AAA' : 1.000, 'AAC' : 1.000, 'AAG' : 0.224, 'AAT' : 0.059, 
        'ACA' : 0.123, 'ACC' : 0.985, 'ACG' : 0.121, 'ACT' : 1.000, 
        'AGA' : 0.004, 'AGC' : 0.233, 'AGG' : 0.004, 'AGT' : 0.052, 
        'ATA' : 0.005, 'ATC' : 1.000, 'ATG' : 1.000, 'ATT' : 0.204, 
        'CAA' : 0.138, 'CAC' : 1.000, 'CAG' : 1.000, 'CAT' : 0.295, 
        'CCA' : 0.172, 'CCC' : 0.008, 'CCG' : 1.000, 'CCT' : 0.056, 
        'CGA' : 0.006, 'CGC' : 0.358, 'CGG' : 0.005, 'CGT' : 1.000, 
        'CTA' : 0.008, 'CTC' : 0.028, 'CTG' : 1.000, 'CTT' : 0.034, 
        'GAA' : 1.000, 'GAC' : 1.000, 'GAG' : 0.219, 'GAT' : 0.452, 
        'GCA' : 0.904, 'GCC' : 0.124, 'GCG' : 0.543, 'GCT' : 1.000, 
        'GGA' : 0.014, 'GGC' : 0.661, 'GGG' : 0.024, 'GGT' : 1.000, 
        'GTA' : 0.767, 'GTC' : 0.072, 'GTG' : 0.295, 'GTT' : 1.000, 
        'TAA' : 1.000, 'TAC' : 1.000, 'TAG' : 0.033, 'TAT' : 0.266, 
        'TCA' : 0.112, 'TCC' : 0.702, 'TCG' : 0.021, 'TCT' : 1.000, 
        'TGA' : 0.078, 'TGC' : 1.000, 'TGG' : 1.000, 'TGT' : 0.525, 
        'TTA' : 0.048, 'TTC' : 1.000, 'TTG' : 0.037, 'TTT' : 0.304} 
