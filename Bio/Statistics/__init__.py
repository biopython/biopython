# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Basic statistics module (DEPRECATED)."""

import warnings

from Bio import BiopythonDeprecationWarning

warnings.warn(
    "Bio.Statistics has been deprecated, and we intend to remove it "
    "in a future release of Biopython.",
    BiopythonDeprecationWarning,
)
