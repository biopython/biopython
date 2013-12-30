# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

#TODO - Don't use "from XXX import *"
from __future__ import print_function

try:
    from numpy import *
    from numpy import dot  # missing in old PyPy's micronumpy
    from numpy.linalg import svd, det # Missing in PyPy 2.0 numpypy
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.SVDSuperimposer.")

from Bio.SVDSuperimposer import SVDSuperimposer

# start with two coordinate sets (Nx3 arrays - Float0)

x=array([[51.65, -1.90, 50.07],
         [50.40, -1.23, 50.65],
         [50.68, -0.04, 51.54],
         [50.22, -0.02, 52.85]], 'f')

y=array([[51.30, -2.99, 46.54],
         [51.09, -1.88, 47.58],
         [52.36, -1.20, 48.03],
         [52.71, -1.18, 49.38]], 'f')

sup=SVDSuperimposer()

# set the coords
# y will be rotated and translated on x
sup.set(x, y)

# do the lsq fit
sup.run()

# get the rmsd
rms=sup.get_rms()

# get rotation (right multiplying!) and the translation
rot, tran=sup.get_rotran()

# rotate y on x manually
y_on_x1=dot(y, rot)+tran

# same thing
y_on_x2=sup.get_transformed()


def simple_matrix_print(matrix):
    """Simple string to display a floating point matrix

    This should give the same output on multiple systems.  This is
    needed because a simple "print matrix" uses scientific notation
    which varies between platforms.

    Only 4 decimal places are used to avoid false test failures due
    to slight differences in the calculation (e.g. due to different
    versions of the underlying libraries or the compilation options
    they used).
    """
    return "[%s]" % "\n ".join("[%s]" % " ".join("% 1.4f" % v for v in row)
                               for row in matrix)


# output results
print(simple_matrix_print(y_on_x1))
print("")
print(simple_matrix_print(y_on_x2))
print("")
print("%.2f" % rms)
