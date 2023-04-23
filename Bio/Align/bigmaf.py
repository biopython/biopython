# Copyright 2022 by Michiel de Hoon.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.Align support for the "bigmaf" multiple alignment format.

The bigMaf format stores multiple alignments in a format compatible with
the MAF (Multiple Alignment Format) format. BigMaf files are binary and are
indexed as a bigBed file.

See https://genome.ucsc.edu/goldenPath/help/bigMaf.html
"""
from io import StringIO


from Bio.Align import interfaces, bigbed, maf


class AlignmentIterator(bigbed.AlignmentIterator, maf.AlignmentIterator):
    """Alignment iterator for bigMaf files.

    The file may contain multiple alignments, which are loaded and returned
    incrementally.

    Alignment annotations are stored in the ``.annotations`` attribute of the
    ``Alignment`` object, except for the alignment score, which is stored as an
    attribute.  Sequence information of empty parts in the alignment block
    (sequences that connect the previous alignment block to the next alignment
    block, but do not align to the current alignment block) is stored in the
    alignment annotations under the ``"empty"`` key.  Annotations specific to
    each line in the alignment are stored in the ``.annotations`` attribute of
    the corresponding sequence record.
    """

    fmt = "bigMaf"
    mode = "b"

    def _create_alignment(self, chunk):
        chromId, chromStart, chromEnd, rest = chunk
        data = rest.decode().replace(";", "\n")
        stream = StringIO()
        stream.write(data)
        stream.seek(0)
        alignment = maf.AlignmentIterator._create_alignment(self, stream)
        return alignment
