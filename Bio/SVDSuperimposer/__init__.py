# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#

"""
SVDSuperimposer finds the best rotation and translation to put
two point sets on top of each other (minimizing the RMSD). This is
eg. useful to superimpose crystal structures. SVD stands for singular
value decomposition, which is used in the algorithm.
"""

from .SVDSuperimposer import SVDSuperimposer
