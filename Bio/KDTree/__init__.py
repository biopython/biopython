# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#

"""KD tree data structure for searching N-dimensional vectors.

The KD tree data structure can be used for all kinds of searches that
involve N-dimensional vectors. For example, neighbor searches (find all points
within a radius of a given point) or finding all point pairs in a set
that are within a certain radius of each other. See "Computational Geometry:
Algorithms and Applications" (Mark de Berg, Marc van Kreveld, Mark Overmars,
Otfried Schwarzkopf).
"""

from .KDTree import KDTree

__docformat__ = "restructuredtext en"
