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


import warnings

from Bio import BiopythonDeprecationWarning


class SummaryInfo:
    """Calculate summary info about the alignment.  (DEPRECATED)

    This class should be used to calculate information summarizing the
    results of an alignment. This may either be straight consensus info
    or more complicated things.
    """

    def __init__(self, alignment):
        """Initialize with the alignment to calculate information on.

        ic_vector attribute. A list of ic content for each column number.
        """
        warnings.warn(
            """\
The class SummaryInfo has been deprecated. Instead of

>>> align_info = AlignInfo.SummaryInfo(msa)
>>> sequence = align_info.get_column(1)

please use

>>> alignment = msa.alignment  # to get a new-style Alignment object
>>> sequence = alignment[:, 1]

Here, `msa` is a MultipleSeqAlignment object and `alignment` is an
`Alignment` object.""",
            BiopythonDeprecationWarning,
        )
        self.alignment = alignment
        self.ic_vector = []

    def get_column(self, col):
        """Return column of alignment."""
        # TODO - Deprecate this and implement slicing?
        return self.alignment[:, col]
