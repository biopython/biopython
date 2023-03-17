# Copyright (C) 2002, Thomas Hamelryck (thamelry@vub.ac.be)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Align on protein structure onto another using SVD alignment.

SVDSuperimposer finds the best rotation and translation to put
two point sets on top of each other (minimizing the RMSD). This is
eg. useful to superimpose crystal structures. SVD stands for singular
value decomposition, which is used in the algorithm.
"""


try:
    from numpy import dot, transpose, sqrt
    from numpy.linalg import svd, det
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.SVDSuperimposer."
    )


class SVDSuperimposer:
    """Class to run SVD alignment.

    SVDSuperimposer finds the best rotation and translation to put
    two point sets on top of each other (minimizing the RMSD). This is
    eg. useful to superimpose crystal structures.

    SVD stands for Singular Value Decomposition, which is used to calculate
    the superposition.

    Reference:

    Matrix computations, 2nd ed. Golub, G. & Van Loan, CF., The Johns
    Hopkins University Press, Baltimore, 1989

    start with two coordinate sets (Nx3 arrays - float)

    >>> from Bio.SVDSuperimposer import SVDSuperimposer
    >>> from numpy import array, dot, set_printoptions
    >>>
    >>> x = array([[51.65, -1.90, 50.07],
    ...      [50.40, -1.23, 50.65],
    ...      [50.68, -0.04, 51.54],
    ...      [50.22, -0.02, 52.85]], 'f')
    >>>
    >>> y = array([[51.30, -2.99, 46.54],
    ...      [51.09, -1.88, 47.58],
    ...      [52.36, -1.20, 48.03],
    ...      [52.71, -1.18, 49.38]], 'f')

    start

    >>> sup = SVDSuperimposer()

    set the coords y will be rotated and translated on x

    >>> sup.set(x, y)

    do the lsq fit

    >>> sup.run()

    get the rmsd

    >>> rms = sup.get_rms()

    get rotation (right multiplying!) and the translation

    >>> rot, tran = sup.get_rotran()

    rotate y on x

    >>> y_on_x1 = dot(y, rot) + tran

    same thing

    >>> y_on_x2 = sup.get_transformed()

    >>> set_printoptions(precision=2)
    >>> print(y_on_x1)
    [[  5.17e+01  -1.90e+00   5.01e+01]
     [  5.04e+01  -1.23e+00   5.06e+01]
     [  5.07e+01  -4.16e-02   5.15e+01]
     [  5.02e+01  -1.94e-02   5.29e+01]]
    >>> print(y_on_x2)
    [[  5.17e+01  -1.90e+00   5.01e+01]
     [  5.04e+01  -1.23e+00   5.06e+01]
     [  5.07e+01  -4.16e-02   5.15e+01]
     [  5.02e+01  -1.94e-02   5.29e+01]]
    >>> print("%.2f" % rms)
    0.00

    """

    def __init__(self):
        """Initialize the class."""
        self._clear()

    # Private methods

    def _clear(self):
        self.reference_coords = None
        self.coords = None
        self.transformed_coords = None
        self.rot = None
        self.tran = None
        self.rms = None
        self.init_rms = None

    def _rms(self, coords1, coords2):
        """Return rms deviations between coords1 and coords2 (PRIVATE)."""
        diff = coords1 - coords2
        return sqrt(sum(sum(diff * diff)) / coords1.shape[0])

    # Public methods

    def set(self, reference_coords, coords):
        """Set the coordinates to be superimposed.

        coords will be put on top of reference_coords.

        - reference_coords: an NxDIM array
        - coords: an NxDIM array

        DIM is the dimension of the points, N is the number
        of points to be superimposed.
        """
        # clear everything from previous runs
        self._clear()
        # store coordinates
        self.reference_coords = reference_coords
        self.coords = coords
        n = reference_coords.shape
        m = coords.shape
        if n != m or not (n[1] == m[1] == 3):
            raise Exception("Coordinate number/dimension mismatch.")
        self.n = n[0]

    def run(self):
        """Superimpose the coordinate sets."""
        if self.coords is None or self.reference_coords is None:
            raise Exception("No coordinates set.")
        coords = self.coords
        reference_coords = self.reference_coords
        # center on centroid
        av1 = sum(coords) / self.n
        av2 = sum(reference_coords) / self.n
        coords = coords - av1
        reference_coords = reference_coords - av2
        # correlation matrix
        a = dot(transpose(coords), reference_coords)
        u, d, vt = svd(a)
        self.rot = transpose(dot(transpose(vt), transpose(u)))
        # check if we have found a reflection
        if det(self.rot) < 0:
            vt[2] = -vt[2]
            self.rot = transpose(dot(transpose(vt), transpose(u)))
        self.tran = av2 - dot(av1, self.rot)

    def get_transformed(self):
        """Get the transformed coordinate set."""
        if self.coords is None or self.reference_coords is None:
            raise Exception("No coordinates set.")
        if self.rot is None:
            raise Exception("Nothing superimposed yet.")
        if self.transformed_coords is None:
            self.transformed_coords = dot(self.coords, self.rot) + self.tran
        return self.transformed_coords

    def get_rotran(self):
        """Right multiplying rotation matrix and translation."""
        if self.rot is None:
            raise Exception("Nothing superimposed yet.")
        return self.rot, self.tran

    def get_init_rms(self):
        """Root mean square deviation of untransformed coordinates."""
        if self.coords is None:
            raise Exception("No coordinates set yet.")
        if self.init_rms is None:
            self.init_rms = self._rms(self.coords, self.reference_coords)
        return self.init_rms

    def get_rms(self):
        """Root mean square deviation of superimposed coordinates."""
        if self.rms is None:
            transformed_coords = self.get_transformed()
            self.rms = self._rms(transformed_coords, self.reference_coords)
        return self.rms


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest(verbose=0)
