# Copyright (C) 2015, Anuj Sharma (anuj.sharma80@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
QCPSuperimposer finds the best rotation and translation to put
two point sets on top of each other (minimizing the RMSD). This is
eg. useful to superimpose crystal structures. QCP stands for
Quaternion Characteristic Polynomial, which is used in the algorithm.
"""

from __future__ import print_function

from numpy import dot, sqrt, array, matrix, inner, zeros
from .qcprotmodule import FastCalcRMSDAndRotation


class QCPSuperimposer(object):
    """Quaternion Characteristic Polynomial (QCP) Superimposer.

    QCPSuperimposer finds the best rotation and translation to put
    two point sets on top of each other (minimizing the RMSD). This is
    eg. useful to superimposing 3D structures of proteins.

    QCP stands for Quaternion Characteristic Polynomial, which is used
    in the algorithm.

    Reference:

    Douglas L Theobald (2005), "Rapid calculation of RMSDs using a
    quaternion-based characteristic polynomial.", Acta Crystallogr
    A 61(4):478-480
    """

    def __init__(self):
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
        """Return rms deviations between coords1 and coords2."""
        diff = coords1 - coords2
        l = coords1.shape[0]
        return sqrt(sum(dot(diff, diff)) / l)

    def _inner_product(self, coords1, coords2):
        G1 = inner(coords1, coords1).diagonal().sum()
        G2 = inner(coords2, coords2).diagonal().sum()
        A = dot(coords1.T, coords2)
        return ((G1 + G2) / 2, A)

    def _align(self, centered_coords1, centered_coords2):
        (E0, A) = self._inner_product(centered_coords1, centered_coords2)
        (rmsd, r0, r1, r2, r3, r4, r5, r6, r7, r8, q1, q2, q3, q4) = FastCalcRMSDAndRotation(
            A[0][0], A[0][1], A[0][2], A[1][0], A[1][1], A[1][2], A[2][0], A[2][1], A[2][2],
            E0, len(centered_coords1), -1.0)
        rot = array([r0, r1, r2, r3, r4, r5, r6, r7, r8]).reshape(3, 3)
        return (rmsd, rot.T, [q1, q2, q3, q4])

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
        # store cordinates
        self.reference_coords = reference_coords
        self.coords = coords
        n = reference_coords.shape
        m = coords.shape
        if n != m or not(n[1] == m[1] == 3):
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
        #
        (self.rms, self.rot, self.lquart) = self._align(
            coords, reference_coords)
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
            raise Exception("Nothing superimposed yet.")
        return self.rms


if __name__ == "__main__":
    from datetime import datetime

    # start with two coordinate sets (Nx3 arrays - float)

    x = array([[-2.803, -15.373, 24.556],
               [0.893, -16.062, 25.147],
               [1.368, -12.371, 25.885],
               [-1.651, -12.153, 28.177],
               [-0.440, -15.218, 30.068],
               [2.551, -13.273, 31.372],
               [0.105, -11.330, 33.567]], 'f')

    y = array([[-14.739, -18.673, 15.040],
               [-12.473, -15.810, 16.074],
               [-14.802, -13.307, 14.408],
               [-17.782, -14.852, 16.171],
               [-16.124, -14.617, 19.584],
               [-15.029, -11.037, 18.902],
               [-18.577, -10.001, 17.996]], 'f')

    t0 = datetime.now()
    for loop in range(0, 10000):
        # start!
        sup = QCPSuperimposer()

        # set the coords
        # y will be rotated and translated on x
        sup.set(x, y)

        # do the lsq fit
        sup.run()

        # get the rmsd
        rms = sup.get_rms()

        # get rotation (right multiplying!) and the translation
        rot, tran = sup.get_rotran()

        # rotate y on x
        y_on_x1 = dot(y, rot) + tran

        # same thing
        y_on_x2 = sup.get_transformed()
    t1 = datetime.now()
    dif = t1 - t0
    print("process wall time (msec): %d" % (dif.total_seconds() * 1000))

    print(y_on_x1)
    print("")
    print(y_on_x2)
    print("")
    print("%.2f" % rms)
