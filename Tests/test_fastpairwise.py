#!/usr/bin/env python

# Regression tests for Bio.Align.pairwise

from Bio.Align import fastpairwise


test_sequences = [
    ('GAACT', 'GAT'),
    ('AxBx', 'zABz'),
    ]
gap_penalties = [
    (0, 0),
    (-1, 0),
    (-0.1, -0.1)
    ]

for seq1, seq2 in test_sequences:
    print "Aligning '%s' to '%s'" % (seq1, seq2)
    
    for open, extend in gap_penalties:
        print "  Gap penalties are %g, %g" % (open, extend)
        alignments = fastpairwise.align_local(seq1, seq2, open, extend)
        alignments.sort()
        for align1, align2, score, begin, end in alignments:
            print "  Score=%g" % score
            print "    %s" % align1
            print "    %s%s" % (" "*begin, "|"*(end-begin))
            print "    %s" % align2
 
