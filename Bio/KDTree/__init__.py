# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#

"""KD tree data structure for searching N-dimensional vectors (DEPRECATED).

The KD tree data structure can be used for all kinds of searches that
involve N-dimensional vectors. For example, neighbor searches (find all points
within a radius of a given point) or finding all point pairs in a set
that are within a certain radius of each other. See "Computational Geometry:
Algorithms and Applications" (Mark de Berg, Marc van Kreveld, Mark Overmars,
Otfried Schwarzkopf).

This module is DEPRECATED; its replacement is Bio.PDB.kdtrees.
"""

from .KDTree import KDTree

import warnings
from Bio import BiopythonDeprecationWarning
warnings.warn("Bio.KDTree has been deprecated, and we intend to remove it"
              " in a future release of Biopython. Please use Bio.PDB.kdtrees"
              " instead, which is functionally very similar.",
              BiopythonDeprecationWarning)
