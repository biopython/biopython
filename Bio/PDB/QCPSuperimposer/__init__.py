# Copyright (C) 2015, Anuj Sharma (anuj.sharma80@gmail.com)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Structural alignment using Quaternion Characteristic Polynomial (QCP).

QCPSuperimposer finds the best rotation and translation to put
two point sets on top of each other (minimizing the RMSD). This is
eg. useful to superimpose crystal structures. QCP stands for
Quaternion Characteristic Polynomial, which is used in the algorithm.
"""


from numpy import dot, sqrt, array, inner
from .qcprotmodule import FastCalcRMSDAndRotation


class QCPSuperimposer:
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
        return sqrt(sum(dot(diff, diff)) / coords1.shape[0])

    def _inner_product(self, coords1, coords2):
        G1 = inner(coords1, coords1).diagonal().sum()
        G2 = inner(coords2, coords2).diagonal().sum()
        A = dot(coords1.T, coords2)
        return ((G1 + G2) / 2, A)

    def _align(self, centered_coords1, centered_coords2):
        (E0, A) = self._inner_product(centered_coords1, centered_coords2)
        (
            rmsd,
            r0,
            r1,
            r2,
            r3,
            r4,
            r5,
            r6,
            r7,
            r8,
            q1,
            q2,
            q3,
            q4,
        ) = FastCalcRMSDAndRotation(
            A[0][0],
            A[0][1],
            A[0][2],
            A[1][0],
            A[1][1],
            A[1][2],
            A[2][0],
            A[2][1],
            A[2][2],
            E0,
            len(centered_coords1),
            -1.0,
        )
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
        if n != m or n[1] != 3 or m[1] != 3:
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
        (self.rms, self.rot, self.lquart) = self._align(coords, reference_coords)
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
