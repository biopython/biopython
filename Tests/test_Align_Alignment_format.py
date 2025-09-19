"""Tests for controling the number of terminal columns
   Issue #5035
"""
from Bio.Align import PairwiseAligner
aligner = PairwiseAligner()
seq1="VQLQESDAELVKPGASVKISCKASGYTFTDHVIHWVKQKPEQGLEWIGYISPGNGDIKYNEKFKGKATLTADKSSSTAYMQLNSLTSEDSAVYLCKRGYY"
seq2="DVQLQESGPGLVKPSQSQSLTCTVTGYSITSDYAWNWIRQFPGNKLEWMGYMSYSGSTRYNPSLRSRISITRDTSKNQFFLQLKSVTTEDTATYFCARGW"
alignments = aligner.align(seq1, seq2)
best_alignment=alignments[0]
best_alignment.terminal_columns = 80
assert (len(str(best_alignment).split("\n")) == 8), "Expecting alignment to be output as two groups."
best_alignment.terminal_columns = 200
assert (len(str(best_alignment).split("\n")) == 4), "Expecting alignment to be output as one group."
