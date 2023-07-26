# Copyright (C) 2022, Joao Rodrigues (j.p.g.l.m.rodrigues@gmail.com
#                     Anuj Sharma (anuj.sharma80@gmail.com)
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

Algorithm and original code described in:

Theobald DL.
Rapid calculation of RMSDs using a quaternion-based characteristic polynomial.
Acta Crystallogr A. 2005 Jul;61(Pt 4):478-80. doi: 10.1107/S0108767305015266.
Epub 2005 Jun 23. PMID: 15973002.
"""


import numpy as np

from Bio.PDB.PDBExceptions import PDBException


def qcp(coords1, coords2, natoms):
    """Implement the QCP code in Python.

    Input coordinate arrays must be centered at the origin and have
    shape Nx3.

    Variable names match (as much as possible) the C implementation.
    """
    # Original code has coords1 be the mobile. I think it makes more sense
    # for it to be the reference, so I swapped here.
    G1 = np.trace(np.dot(coords2, coords2.T))
    G2 = np.trace(np.dot(coords1, coords1.T))
    A = np.dot(coords2.T, coords1)  # referred to as M in the original paper.
    E0 = (G1 + G2) * 0.5

    Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz = A.flatten()

    Sxx2 = Sxx * Sxx
    Syy2 = Syy * Syy
    Szz2 = Szz * Szz
    Sxy2 = Sxy * Sxy
    Syz2 = Syz * Syz
    Sxz2 = Sxz * Sxz
    Syx2 = Syx * Syx
    Szy2 = Szy * Szy
    Szx2 = Szx * Szx

    SyzSzymSyySzz2 = 2.0 * (Syz * Szy - Syy * Szz)
    Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2

    C2 = -2.0 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2)
    C1 = 8.0 * (
        Sxx * Syz * Szy
        + Syy * Szx * Sxz
        + Szz * Sxy * Syx
        - Sxx * Syy * Szz
        - Syz * Szx * Sxy
        - Szy * Syx * Sxz
    )

    SxzpSzx = Sxz + Szx
    SyzpSzy = Syz + Szy
    SxypSyx = Sxy + Syx
    SyzmSzy = Syz - Szy
    SxzmSzx = Sxz - Szx
    SxymSyx = Sxy - Syx
    SxxpSyy = Sxx + Syy
    SxxmSyy = Sxx - Syy
    Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2

    negSxzpSzx = -SxzpSzx
    negSxzmSzx = -SxzmSzx
    negSxymSyx = -SxymSyx
    SxxpSyy_p_Szz = SxxpSyy + Szz

    C0 = (
        Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2
        + (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2)
        * (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2)
        + (negSxzpSzx * (SyzmSzy) + (SxymSyx) * (SxxmSyy - Szz))
        * (negSxzmSzx * (SyzpSzy) + (SxymSyx) * (SxxmSyy + Szz))
        + (negSxzpSzx * (SyzpSzy) - (SxypSyx) * (SxxpSyy - Szz))
        * (negSxzmSzx * (SyzmSzy) - (SxypSyx) * SxxpSyy_p_Szz)
        + (+(SxypSyx) * (SyzpSzy) + (SxzpSzx) * (SxxmSyy + Szz))
        * (negSxymSyx * (SyzmSzy) + (SxzpSzx) * SxxpSyy_p_Szz)
        + (+(SxypSyx) * (SyzmSzy) + (SxzmSzx) * (SxxmSyy - Szz))
        * (negSxymSyx * (SyzpSzy) + (SxzmSzx) * (SxxpSyy - Szz))
    )

    # Find largest root of the quaternion polynomial:
    #   f(x) = x ** 4 + c2 * x ** 2 + c1 * x + c0 = 0 (eq. 8 of Theobald et al.)
    #   f'(x) = 4 * x ** 3 + 2 * c2 * x + c1
    #
    # using Newton-Rhapson and E0 as initial guess. Liu et al. mentions 5
    # iterations are sufficient (on average) for convergence up to 1e-6
    # precision but original code writes 50, which we keep.
    nr_it = 50
    mxEigenV = E0  # starting guess (x in eqs above)
    evalprec = 1e-11  # convergence criterion
    for _ in range(nr_it):
        oldg = mxEigenV

        x2 = mxEigenV * mxEigenV
        b = (x2 + C2) * mxEigenV
        a = b + C1

        f = a * mxEigenV + C0
        f_prime = 2.0 * x2 * mxEigenV + b + a

        delta = f / (f_prime + evalprec)  # avoid division by zero
        mxEigenV = abs(mxEigenV - delta)
        if (mxEigenV - oldg) < (evalprec * mxEigenV):
            break  # convergence
    else:
        print(f"Newton-Rhapson did not converge after {nr_it} iterations")

    # The original code has a guard if minScore > 0 and rmsd < minScore, although
    # the default value of minScore is -1. For simplicity, we ignore that check.
    rmsd = (2.0 * abs(E0 - mxEigenV) / natoms) ** 0.5

    a11 = SxxpSyy + Szz - mxEigenV
    a12 = SyzmSzy
    a13 = negSxzmSzx
    a14 = SxymSyx
    a21 = SyzmSzy
    a22 = SxxmSyy - Szz - mxEigenV
    a23 = SxypSyx
    a24 = SxzpSzx
    a31 = a13
    a32 = a23
    a33 = Syy - Sxx - Szz - mxEigenV
    a34 = SyzpSzy
    a41 = a14
    a42 = a24
    a43 = a34
    a44 = Szz - SxxpSyy - mxEigenV
    a3344_4334 = a33 * a44 - a43 * a34
    a3244_4234 = a32 * a44 - a42 * a34
    a3243_4233 = a32 * a43 - a42 * a33
    a3143_4133 = a31 * a43 - a41 * a33
    a3144_4134 = a31 * a44 - a41 * a34
    a3142_4132 = a31 * a42 - a41 * a32
    q1 = a22 * a3344_4334 - a23 * a3244_4234 + a24 * a3243_4233
    q2 = -a21 * a3344_4334 + a23 * a3144_4134 - a24 * a3143_4133
    q3 = a21 * a3244_4234 - a22 * a3144_4134 + a24 * a3142_4132
    q4 = -a21 * a3243_4233 + a22 * a3143_4133 - a23 * a3142_4132

    qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4

    evecprec = 1e-6
    if qsqr < evecprec:
        q1 = a12 * a3344_4334 - a13 * a3244_4234 + a14 * a3243_4233
        q2 = -a11 * a3344_4334 + a13 * a3144_4134 - a14 * a3143_4133
        q3 = a11 * a3244_4234 - a12 * a3144_4134 + a14 * a3142_4132
        q4 = -a11 * a3243_4233 + a12 * a3143_4133 - a13 * a3142_4132
        qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4

        if qsqr < evecprec:
            a1324_1423 = a13 * a24 - a14 * a23
            a1224_1422 = a12 * a24 - a14 * a22
            a1223_1322 = a12 * a23 - a13 * a22
            a1124_1421 = a11 * a24 - a14 * a21
            a1123_1321 = a11 * a23 - a13 * a21
            a1122_1221 = a11 * a22 - a12 * a21

            q1 = a42 * a1324_1423 - a43 * a1224_1422 + a44 * a1223_1322
            q2 = -a41 * a1324_1423 + a43 * a1124_1421 - a44 * a1123_1321
            q3 = a41 * a1224_1422 - a42 * a1124_1421 + a44 * a1122_1221
            q4 = -a41 * a1223_1322 + a42 * a1123_1321 - a43 * a1122_1221
            qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4

            if qsqr < evecprec:
                q1 = a32 * a1324_1423 - a33 * a1224_1422 + a34 * a1223_1322
                q2 = -a31 * a1324_1423 + a33 * a1124_1421 - a34 * a1123_1321
                q3 = a31 * a1224_1422 - a32 * a1124_1421 + a34 * a1122_1221
                q4 = -a31 * a1223_1322 + a32 * a1123_1321 - a33 * a1122_1221
                qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4

                if qsqr < evecprec:
                    rot = np.eye(3)
                    return rmsd, rot, [q1, q2, q3, q4]

    normq = qsqr**0.5
    q1 /= normq
    q2 /= normq
    q3 /= normq
    q4 /= normq

    a2 = q1 * q1
    x2 = q2 * q2
    y2 = q3 * q3
    z2 = q4 * q4

    xy = q2 * q3
    az = q1 * q4
    zx = q4 * q2
    ay = q1 * q3
    yz = q3 * q4
    ax = q1 * q2

    rot = np.zeros((3, 3))

    rot[0][0] = a2 + x2 - y2 - z2
    rot[0][1] = 2 * (xy + az)
    rot[0][2] = 2 * (zx - ay)
    rot[1][0] = 2 * (xy - az)
    rot[1][1] = a2 - x2 + y2 - z2
    rot[1][2] = 2 * (yz + ax)
    rot[2][0] = 2 * (zx + ay)
    rot[2][1] = 2 * (yz - ax)
    rot[2][2] = a2 - x2 - y2 + z2

    return rmsd, rot, (q1, q2, q3, q4)


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
        self._reset_properties()

    # Private methods

    def _reset_properties(self):
        """Reset all relevant properties to None to avoid conflicts between runs."""
        self.reference_coords = None
        self.coords = None
        self.transformed_coords = None
        self.rot = None
        self.tran = None
        self.rms = None
        self.init_rms = None

    # Public methods
    def set_atoms(self, fixed, moving):
        """Prepare alignment between two atom lists.

        Put (translate/rotate) the atoms in fixed on the atoms in
        moving, in such a way that the RMSD is minimized.

        :param fixed: list of (fixed) atoms
        :param moving: list of (moving) atoms
        :type fixed,moving: [L{Atom}, L{Atom},...]
        """
        assert len(fixed) == len(moving), "Fixed and moving atom lists differ in size"

        # Grab coordinates in double precision
        fix_coord = np.array([a.get_coord() for a in fixed], dtype=np.float64)
        mov_coord = np.array([a.get_coord() for a in moving], dtype=np.float64)

        self.set(fix_coord, mov_coord)
        self.run()

        self.rms = self.get_rms()
        self.rotran = self.get_rotran()

    def apply(self, atom_list):
        """Apply the QCP rotation matrix/translation vector to a set of atoms."""
        if self.rotran is None:
            raise PDBException("No transformation has been calculated yet")

        rot, tran = self.rotran
        for atom in atom_list:
            atom.transform(rot, tran)

    # Low(er) level functions
    def set(self, reference_coords, coords):
        """Set the coordinates to be superimposed.

        coords will be put on top of reference_coords.

        - reference_coords: an NxDIM array
        - coords: an NxDIM array

        DIM is the dimension of the points, N is the number
        of points to be superimposed.
        """
        self._reset_properties()

        # store coordinates
        self.reference_coords = reference_coords
        self.coords = coords
        self._natoms, n_dim = coords.shape

        if reference_coords.shape != coords.shape:
            raise PDBException("Coordinates must have the same dimensions.")
        if n_dim != 3:
            raise PDBException("Coordinates must be Nx3 arrays.")

    def run(self):
        """Superimpose the coordinate sets."""
        if self.coords is None or self.reference_coords is None:
            raise PDBException("No coordinates set.")

        coords = self.coords.copy()
        coords_ref = self.reference_coords.copy()

        # Center Coordinates
        com_coords = np.mean(coords, axis=0)
        com_ref = np.mean(coords_ref, axis=0)

        coords -= com_coords
        coords_ref -= com_ref

        (self.rms, self.rot, _) = qcp(coords_ref, coords, self._natoms)
        self.tran = com_ref - np.dot(com_coords, self.rot)

    # Getters
    def get_transformed(self):
        """Get the transformed coordinate set."""
        if self.coords is None or self.reference_coords is None:
            raise PDBException("No coordinates set.")

        if self.rot is None:
            raise PDBException("Nothing is superimposed yet.")

        self.transformed_coords = np.dot(self.coords, self.rot) + self.tran
        return self.transformed_coords

    def get_rotran(self):
        """Return right multiplying rotation matrix and translation vector."""
        if self.rot is None:
            raise PDBException("Nothing is superimposed yet.")
        return self.rot, self.tran

    def get_init_rms(self):
        """Return the root mean square deviation of untransformed coordinates."""
        if self.coords is None:
            raise PDBException("No coordinates set yet.")

        if self.init_rms is None:
            diff = self.coords - self.reference_coords
            self.init_rms = np.sqrt(np.sum(np.sum(diff * diff, axis=1) / self._natoms))
        return self.init_rms

    def get_rms(self):
        """Root mean square deviation of superimposed coordinates."""
        if self.rms is None:
            raise PDBException("Nothing superimposed yet.")
        return self.rms
