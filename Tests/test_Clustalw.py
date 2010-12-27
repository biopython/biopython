# Copyright 2008 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
from StringIO import StringIO
from Bio import AlignIO
import Bio

import warnings
warnings.filterwarnings("ignore", category=Bio.BiopythonDeprecationWarning)
from Bio import Clustalw
warnings.filters.pop()

from Bio.Align.Generic import Alignment


# test_files is a list of tuples containing:
# - integer: number of sequences per alignment
# - integer: number of alignments
# - string:  relative filename
test_files = [ \
    #clustalw 1.81
    (2, 1, 'Clustalw/cw02.aln'),
    (7, 1, 'Clustalw/opuntia.aln'),
    (20, 1, 'Clustalw/protein.aln'),
    # clustalw 2.0.8
    (2, 1, 'Clustalw/odd_consensus.aln'),
    # clustalw 2.0.9
    (5, 1, 'Clustalw/hedgehog.aln'),
    # PROMALS3D
    (20, 1, 'Clustalw/promals3d.aln'),
    ]


def compare(a1, a2):
    assert a1.get_alignment_length() == a2.get_alignment_length()
    assert len(a1._records) == len(a2._records)
    for r1, r2 in zip(a1,a2):
        assert r1.id == r2.id
        assert str(r1.seq) == str(r2.seq)

    if hasattr(a1, "_version") and a1._version \
    and hasattr(a2, "_version") and a2._version:
        assert a1._version == a2._version

    if hasattr(a1, "_star_info") and a1._star_info \
    and hasattr(a2, "_star_info") and a2._star_info:
        assert a1._star_info == a2._star_info \
    
    return True


print "Checking Bio.AlignIO and Bio.Clustalw can read example files..."
for (t_per, t_count, t_filename) in test_files:
    print
    print "Testing reading %s format file %s with %i alignments" \
          % ("clustal", t_filename, t_count)
    assert os.path.isfile(t_filename), t_filename

    print "Using Bio.AlignIO.parse(...)"
    alignments  = list(AlignIO.parse(t_filename, "clustal"))
    assert len(alignments)  == t_count, \
         "Found %i alignments but expected %i" % (len(alignments), t_count)
    for alignment in alignments:
        assert len(alignment) == t_per, \
            "Expected %i records per alignment, got %i" \
            % (t_per, len(alignment))
        print
        print alignment
    
    if t_count != 1 : continue
    print

    # Check Bio.AlignIO.read(...)
    alignment = AlignIO.read(t_filename, "clustal")
    assert isinstance(alignment, Alignment)
    assert compare(alignment, alignments[0])

    print "Using Bio.AlignIO.read(...)"
    #print "~" * 75
    #handle = StringIO()
    #AlignIO.write([alignment], handle, "clustal")
    #handle.seek(0)
    #print handle.read()
    #print "~" * 75

    print "Using Bio.Clustalw.parse_file(...)"
    c_alignment = Clustalw.parse_file(t_filename)
    assert isinstance(c_alignment, Alignment)
    assert isinstance(c_alignment, Clustalw.ClustalAlignment)

    #print "  Using Bio.Clustalw.parse_file(...)"
    #print "~" * 75
    #print c_alignment
    #print "~" * 75
    #print

    # Compare the two...
    assert compare(alignment, c_alignment)

    # Check Bio.AlignIO can read the Bio.Clustalw's string output
    n_alignment = AlignIO.read(StringIO(str(c_alignment)), "clustal")
    assert isinstance(alignment, Alignment)
    assert compare(n_alignment, c_alignment)

print "Finished tested reading files"
