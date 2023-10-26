# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Legacy functionalities from other parts of Biopython used by SearchIO (DEPRECATED)."""


import warnings
from Bio import BiopythonDeprecationWarning


warnings.warn(
    "The 'Bio.SearchIO._legacy' module for parsing BLAST plain text output is "
    "deprecated and will be removed in a future release of Biopython. "
    "Consider generating your BLAST output for parsing as XML or tabular "
    "format instead.",
    BiopythonDeprecationWarning,
)
