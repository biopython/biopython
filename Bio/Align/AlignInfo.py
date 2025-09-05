# Copyright 2000 Brad Chapman.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Extract information from alignment objects.

In order to try and avoid huge alignment objects with tons of functions,
functions which return summary type information about alignments should
be put into classes in this module.
"""
from Bio.Align import Alignment, MultipleSeqAlignment


class SummaryInfo:
    """Calculate summary info about the alignment.

    This class should be used to calculate information summarizing the
    results of an alignment. This may either be straight consensus info
    or more complicated things.
    """

    def __init__(self, alignment: Alignment | MultipleSeqAlignment):
        """Initialize with the alignment to calculate information on.

        ic_vector attribute. A list of ic content for each column number.
        """
        self.alignment = alignment
        self.ic_vector: list[float] = []

    def get_column(self, col: int) -> str:
        """Return column of alignment."""
        # TODO - Deprecate this and implement slicing?
        return self.alignment[:, col]
