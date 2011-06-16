# Copyright 2011 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Simple 'print-and-compare' unit test for fasta-m10 parser.

Created to check for any regressions from my new implementation of the
parser.
"""
import os
from Bio import AlignIO

# test_files is a list of tuples containing:
# - string:  file format
# - integer: number of sequences per alignment
# - integer: number of alignments
# - string:  relative filename
#
# Most of the input files are also used by test_SeqIO.py,
# and by other additional tests as noted below.
test_files = [ \
    ("fasta-m10", 2, 4, 'Fasta/output001.m10'),
    ("fasta-m10", 2, 6, 'Fasta/output002.m10'),
    ("fasta-m10", 2, 3, 'Fasta/output003.m10'),
    ("fasta-m10", 2, 1, 'Fasta/output004.m10'),
    ("fasta-m10", 2, 1, 'Fasta/output005.m10'),
    ("fasta-m10", 2, 1, 'Fasta/output006.m10'),
    ("fasta-m10", 2, 9, 'Fasta/output007.m10'),
    ("fasta-m10", 2, 12,'Fasta/output008.m10'),
    ]

#Main tests...
for (t_format, t_per, t_count, t_filename) in test_files:
    assert t_format == "fasta-m10" and t_per == 2

    print "Testing reading %s format file %s with %i alignments" \
          % (t_format, t_filename, t_count)
    assert os.path.isfile(t_filename), t_filename

    #Try as an iterator using handle
    alignments  = list(AlignIO.parse(handle=open(t_filename,"r"), format=t_format))
    assert len(alignments)  == t_count, \
         "Found %i alignments but expected %i" % (len(alignments), t_count)
    for alignment in alignments:
        assert len(alignment) == t_per, \
            "Expected %i records per alignment, got %i" \
            % (t_per, len(alignment))


    #Print the alignment
    for i,alignment in enumerate(alignments):
        print "="*78
        print "Alignment %i, with %i sequences of length %i" \
              % (i,
                 len(alignment),
                 alignment.get_alignment_length())
        for k in sorted(alignment._annotations):
            print " - %s: %r" % (k, alignment._annotations[k])
        assert alignment[0].name == "query"
        assert alignment[1].name == "match"
        #Show each sequence row horizontally
        for record in alignment:
            print "-"*78
            print record.id
            print record.description
            print repr(record.seq)
            assert not record.features
            assert not record.letter_annotations
            for k in sorted(record.annotations):
                print " - %s: %r" % (k, record.annotations[k])
    print "="*78
print "Finished tested reading files"
