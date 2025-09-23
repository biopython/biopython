"""Tests for controling the number of terminal columns
   Issue #5035
"""
from Bio.Align import PairwiseAligner
aligner = PairwiseAligner()
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -1.5
aligner.extend_gap_score = -0.2
seq1 = "VQLQESDAELVKPGASVKISCKASGYTFTDHVIHWVKQKPEQGLEWIGYISPGNGDIKYNEKFKGKATLTADKSSSTAYMQLNSLTSEDSAVYLCKRGYY"
seq2 = "DVQLQESGPGLVKPSQSQSLTCTVTGYSITSDYAWNWIRQFPGNKLEWMGYMSYSGSTRYNPSLRSRISITRDTSKNQFFLQLKSVTTEDTATYFCARGW"
alignments = aligner.align(seq1, seq2)
best_alignment = alignments[0]

best_alignment.terminal_columns = 80
# target            0 -VQLQESDAELVKPGASVKIS--CKASGYTFT-DHVIHWVKQKPEQG--LEWIGYISPGNGDIKYNEKFKGKATLTADKSSSTAYMQLN-SL---------TS-------------EDSAVYLCKRGYY 100
#                   0 -||||||...||||--|...|--|...||..|-|....|..|.|--|--|||.||.|-------|--------------|.||.|---|-||---------||-------------||.|.|.|.||-. 129
# query             0 DVQLQESGPGLVKP--SQSQSLTCTVTGYSITSDYAWNWIRQFP--GNKLEWMGYMS-------Y--------------SGSTRY---NPSLRSRISITRDTSKNQFFLQLKSVTTEDTATYFCARG-W 100
# 
# (base) mini4:biopython yzhou$ python Tests/test_Align_Alignment_format.py
# target            0 -VQLQESDAELVKPGASVKIS--CKASGYTFT-DHVIHWVKQKPEQG--LEWIGYISPGN
#                   0 -||||||...||||--|...|--|...||..|-|....|..|.|--|--|||.||.|---
# query             0 DVQLQESGPGLVKP--SQSQSLTCTVTGYSITSDYAWNWIRQFP--GNKLEWMGYMS---
# 
# target           54 GDIKYNEKFKGKATLTADKSSSTAYMQLN-SL---------TS-------------EDSA
#                  60 ----|--------------|.||.|---|-||---------||-------------||.|
# query            53 ----Y--------------SGSTRY---NPSLRSRISITRDTSKNQFFLQLKSVTTEDTA
# 
# target           91 VYLCKRGYY 100
#                 120 .|.|.||-. 129
# query            92 TYFCARG-W 100
# 
assert (len(str(best_alignment).split("\n")) == 12), "Expecting alignment to be output as two groups."

best_alignment.terminal_columns = 200
# target            0 -VQLQESDAELVKPGASVKIS--CKASGYTFT-DHVIHWVKQKPEQG--LEWIGYISPGNGDIKYNEKFKGKATLTADKSSSTAYMQLN-SL---------TS-------------EDSAVYLCKRGYY 100
#                   0 -||||||...||||--|...|--|...||..|-|....|..|.|--|--|||.||.|-------|--------------|.||.|---|-||---------||-------------||.|.|.|.||-. 129
# query             0 DVQLQESGPGLVKP--SQSQSLTCTVTGYSITSDYAWNWIRQFP--GNKLEWMGYMS-------Y--------------SGSTRY---NPSLRSRISITRDTSKNQFFLQLKSVTTEDTATYFCARG-W 100
# 
assert (len(str(best_alignment).split("\n")) == 4), "Expecting alignment to be output as one group."
