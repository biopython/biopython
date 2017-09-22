# Copyright 2001 by Brad Chapman.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Genetic Algorithm library (DEPRECATED)."""

import warnings
from Bio import BiopythonDeprecationWarning
warnings.warn("Bio.GA has been deprecated, and we intend to remove"
              " it in a future release of Biopython. Please consider using"
              " DEAP instead.  If you would like to"
              " continue using Bio.GA, please contact the Biopython"
              " developers via the mailing list or GitHub.",
              BiopythonDeprecationWarning)
