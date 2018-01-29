# Copyright 2001 by Brad Chapman.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Neural Network library (DEPRECATED)."""

import warnings
from Bio import BiopythonDeprecationWarning
warnings.warn("Bio.NeuralNetwork has been deprecated, and we intend to remove"
              " it in a future release of Biopython. Please consider using"
              " scikit-learn or TensorFlow instead.  If you would like to"
              " continue using Bio.NeuralNetwork, please contact the Biopython"
              " developers via the mailing list or GitHub.",
              BiopythonDeprecationWarning)
