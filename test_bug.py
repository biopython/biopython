from Bio import Align
aligner = Align.PairwiseAligner()

def my_gap_score_function(start, length):
    if start==2:
        return -1000
    else:
        return -1 * length

aligner.query_gap_score = my_gap_score_function
alignments = aligner.align("AACTT", "AATT")
for alignment in alignments:
    print(alignment)
