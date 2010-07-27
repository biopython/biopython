#!/usr/bin/env python

from Bio import pairwise2


def _pretty_print_align(align1, align2, score, begin, end):
    print pairwise2.format_alignment(align1, align2, score, begin, end)

def _pretty_print_all(aligns):
    aligns.sort()
    for align in aligns:
        _pretty_print_align(*align)
    print

def _align_and_print(fn, *args, **keywds):
    print fn.__name__
    _pretty_print_all(fn(*args, **keywds))

a = pairwise2.align

print "### Test the function generation code"
fn = a.globalxx
print fn.__name__     # globalxx
print fn.__doc__      # globalxx(sequenceA, sequenceB) -> alignments
print a.localcd.__doc__  # localxx(...match_fn, openA, ...) -> alignments
try:
    a.blah
except:
    print "correctly failed"
    

print "#### Let's start with a simple global alignment."
_align_and_print(a.globalxx, "GAACT", "GAT")                 # 2 aligns, 3

print "#### Try a local alignment."
_align_and_print(a.localxs, "AxBx", "zABz", -0.1, 0)         # 2 aligns, 1.9

print "#### Test match score, open penalty."
_align_and_print(a.globalms, "AA", "A", 2.0, -1, -0.1, 0)    # 2 aligns, 1.9
_align_and_print(a.globalms, "GAA", "GA", 1.5, 0, -0.1, 0)   # 2 aligns, 2.9
_align_and_print(a.globalxs, "GAACT", "GAT", -0.1, 0)        # 1 align, 2.9
_align_and_print(a.globalms, "GCT", "GATA", 1, -2, -0.1, 0)  # 1 align, -0.1

print "#### Test the extend penalty."
_align_and_print(a.globalxs, "GACT", "GT", -0.2, -0.5) # 1 align, 1.3
_align_and_print(a.globalxs, "GACT", "GT", -0.2, -1.5) # 2 aligns, 0.6

print "#### Test penalize_extend_when_opening"
_align_and_print(a.globalxs, "GACT", "GT", -0.2, -1.5,
                 penalize_extend_when_opening=1)       # 1 align, -1.2

print "#### Test penalize_end_gaps"
_align_and_print(a.globalxs, "GACT", "GT", -0.2, -0.8,
                 penalize_end_gaps=0)                  # 3 aligns, 1

print "#### Test separate gap penalties"
_align_and_print(a.localxd, "GAT", "GTCT", -0.3, 0, -0.8, 0)  # 2 aligns, 1.7
_align_and_print(a.localxd, "GAT", "GTCT", -0.5, 0, -0.2, 0)  # 1 aligns, 1.8

print "#### Test separate gap penalties, with extension.  Test align list"
_align_and_print(a.localxd, list("GAAT"), list("GTCCT"), -0.1, 0, -0.1, -0.1,
                 gap_char=["-"])                              # 3 aligns, 1.9

print "#### Test match dictionary"
match_dict = {
    ("A", "A") : 1.5,
    ("A", "T") : 0.5,
    ("T", "T") : 1.0
    }
_align_and_print(a.localds, "ATAT", "ATT", match_dict, -.5, 0) # 2 aligns, 3
_align_and_print(a.localds, "ATAT", "ATT", match_dict, -1, 0)  # 1 align, 3
_align_and_print(a.localds, "ATT", "ATAT", match_dict, -1, 0)  # 1 align, 3



print "#### This used to cause errors, reported by Daishi"
_align_and_print(a.localxs, "abcde", "c", -0.3, -0.1)  # 1 align, 1.0, 1 char
_align_and_print(a.localxs, "abcce", "c", -0.3, -0.1)  # 2 aligns, 1.0, 1 char
_align_and_print(a.globalxs, "abcde", "c", -0.3, -0.1) # 1 aligns, 0.2, 1 char
