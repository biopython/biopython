# Copyright (C) 2004, Thomas Hamelryck (thamelry@binf.ku.dk)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Vector class, including rotation-related functions."""


import numpy  # type: ignore
from typing import Tuple, Optional


def m2rotaxis(m):
    """Return angles, axis pair that corresponds to rotation matrix m.

    The case where ``m`` is the identity matrix corresponds to a singularity
    where any rotation axis is valid. In that case, ``Vector([1, 0, 0])``,
    is returned.
    """
    eps = 1e-5

    # Check for singularities a la http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToAngle/
    if (
        abs(m[0, 1] - m[1, 0]) < eps
        and abs(m[0, 2] - m[2, 0]) < eps
        and abs(m[1, 2] - m[2, 1]) < eps
    ):
        # Singularity encountered. Check if its 0 or 180 deg
        if (
            abs(m[0, 1] + m[1, 0]) < eps
            and abs(m[0, 2] + m[2, 0]) < eps
            and abs(m[1, 2] + m[2, 1]) < eps
            and abs(m[0, 0] + m[1, 1] + m[2, 2] - 3) < eps
        ):
            angle = 0
        else:
            angle = numpy.pi
    else:
        # Angle always between 0 and pi
        # Sense of rotation is defined by axis orientation
        t = 0.5 * (numpy.trace(m) - 1)
        t = max(-1, t)
        t = min(1, t)
        angle = numpy.arccos(t)

    if angle < 1e-15:
        # Angle is 0
        return 0.0, Vector(1, 0, 0)
    elif angle < numpy.pi:
        # Angle is smaller than pi
        x = m[2, 1] - m[1, 2]
        y = m[0, 2] - m[2, 0]
        z = m[1, 0] - m[0, 1]
        axis = Vector(x, y, z)
        axis.normalize()
        return angle, axis
    else:
        # Angle is pi - special case!
        m00 = m[0, 0]
        m11 = m[1, 1]
        m22 = m[2, 2]
        if m00 > m11 and m00 > m22:
            x = numpy.sqrt(m00 - m11 - m22 + 0.5)
            y = m[0, 1] / (2 * x)
            z = m[0, 2] / (2 * x)
        elif m11 > m00 and m11 > m22:
            y = numpy.sqrt(m11 - m00 - m22 + 0.5)
            x = m[0, 1] / (2 * y)
            z = m[1, 2] / (2 * y)
        else:
            z = numpy.sqrt(m22 - m00 - m11 + 0.5)
            x = m[0, 2] / (2 * z)
            y = m[1, 2] / (2 * z)
        axis = Vector(x, y, z)
        axis.normalize()
        return numpy.pi, axis


def vector_to_axis(line, point):
    """Vector to axis method.

    Return the vector between a point and
    the closest point on a line (ie. the perpendicular
    projection of the point on the line).

    :type line: L{Vector}
    :param line: vector defining a line

    :type point: L{Vector}
    :param point: vector defining the point
    """
    line = line.normalized()
    np = point.norm()
    angle = line.angle(point)
    return point - line ** (np * numpy.cos(angle))


def rotaxis2m(theta, vector):
    """Calculate left multiplying rotation matrix.

    Calculate a left multiplying rotation matrix that rotates
    theta rad around vector.

    :type theta: float
    :param theta: the rotation angle

    :type vector: L{Vector}
    :param vector: the rotation axis

    :return: The rotation matrix, a 3x3 Numeric array.

    Examples
    --------
    >>> from numpy import pi
    >>> from Bio.PDB.vectors import rotaxis2m
    >>> from Bio.PDB.vectors import Vector
    >>> m = rotaxis2m(pi, Vector(1, 0, 0))
    >>> Vector(1, 2, 3).left_multiply(m)
    <Vector 1.00, -2.00, -3.00>

    """
    vector = vector.normalized()
    c = numpy.cos(theta)
    s = numpy.sin(theta)
    t = 1 - c
    x, y, z = vector.get_array()
    rot = numpy.zeros((3, 3))
    # 1st row
    rot[0, 0] = t * x * x + c
    rot[0, 1] = t * x * y - s * z
    rot[0, 2] = t * x * z + s * y
    # 2nd row
    rot[1, 0] = t * x * y + s * z
    rot[1, 1] = t * y * y + c
    rot[1, 2] = t * y * z - s * x
    # 3rd row
    rot[2, 0] = t * x * z - s * y
    rot[2, 1] = t * y * z + s * x
    rot[2, 2] = t * z * z + c
    return rot


rotaxis = rotaxis2m


def refmat(p, q):
    """Return a (left multiplying) matrix that mirrors p onto q.

    :type p,q: L{Vector}
    :return: The mirror operation, a 3x3 Numeric array.

    Examples
    --------
    >>> from Bio.PDB.vectors import refmat
    >>> p, q = Vector(1, 2, 3), Vector(2, 3, 5)
    >>> mirror = refmat(p, q)
    >>> qq = p.left_multiply(mirror)
    >>> print(q)
    <Vector 2.00, 3.00, 5.00>
    >>> print(qq)
    <Vector 1.21, 1.82, 3.03>

    """
    p = p.normalized()
    q = q.normalized()
    if (p - q).norm() < 1e-5:
        return numpy.identity(3)
    pq = p - q
    pq.normalize()
    b = pq.get_array()
    b.shape = (3, 1)
    i = numpy.identity(3)
    ref = i - 2 * numpy.dot(b, numpy.transpose(b))
    return ref


def rotmat(p, q):
    """Return a (left multiplying) matrix that rotates p onto q.

    :param p: moving vector
    :type p: L{Vector}

    :param q: fixed vector
    :type q: L{Vector}

    :return: rotation matrix that rotates p onto q
    :rtype: 3x3 Numeric array

    Examples
    --------
    >>> from Bio.PDB.vectors import rotmat
    >>> p, q = Vector(1, 2, 3), Vector(2, 3, 5)
    >>> r = rotmat(p, q)
    >>> print(q)
    <Vector 2.00, 3.00, 5.00>
    >>> print(p)
    <Vector 1.00, 2.00, 3.00>
    >>> p.left_multiply(r)
    <Vector 1.21, 1.82, 3.03>

    """
    rot = numpy.dot(refmat(q, -p), refmat(p, -p))
    return rot


def calc_angle(v1, v2, v3):
    """Calculate angle method.

    Calculate the angle between 3 vectors
    representing 3 connected points.

    :param v1, v2, v3: the tree points that define the angle
    :type v1, v2, v3: L{Vector}

    :return: angle
    :rtype: float
    """
    v1 = v1 - v2
    v3 = v3 - v2
    return v1.angle(v3)


def calc_dihedral(v1, v2, v3, v4):
    """Calculate dihedral angle method.

    Calculate the dihedral angle between 4 vectors
    representing 4 connected points. The angle is in
    ]-pi, pi].

    :param v1, v2, v3, v4: the four points that define the dihedral angle
    :type v1, v2, v3, v4: L{Vector}
    """
    ab = v1 - v2
    cb = v3 - v2
    db = v4 - v3
    u = ab ** cb
    v = db ** cb
    w = u ** v
    angle = u.angle(v)
    # Determine sign of angle
    try:
        if cb.angle(w) > 0.001:
            angle = -angle
    except ZeroDivisionError:
        # dihedral=pi
        pass
    return angle


class Vector:
    """3D vector."""

    def __init__(self, x, y=None, z=None):
        """Initialize the class."""
        if y is None and z is None:
            # Array, list, tuple...
            if len(x) != 3:
                raise ValueError("Vector: x is not a list/tuple/array of 3 numbers")
            self._ar = numpy.array(x, "d")
        else:
            # Three numbers
            self._ar = numpy.array((x, y, z), "d")

    def __repr__(self):
        """Return vector 3D coordinates."""
        x, y, z = self._ar
        return "<Vector %.2f, %.2f, %.2f>" % (x, y, z)

    def __neg__(self):
        """Return Vector(-x, -y, -z)."""
        a = -self._ar
        return Vector(a)

    def __add__(self, other):
        """Return Vector+other Vector or scalar."""
        if isinstance(other, Vector):
            a = self._ar + other._ar
        else:
            a = self._ar + numpy.array(other)
        return Vector(a)

    def __sub__(self, other):
        """Return Vector-other Vector or scalar."""
        if isinstance(other, Vector):
            a = self._ar - other._ar
        else:
            a = self._ar - numpy.array(other)
        return Vector(a)

    def __mul__(self, other):
        """Return Vector.Vector (dot product)."""
        return sum(self._ar * other._ar)

    def __truediv__(self, x):
        """Return Vector(coords/a)."""
        a = self._ar / numpy.array(x)
        return Vector(a)

    def __pow__(self, other):
        """Return VectorxVector (cross product) or Vectorxscalar."""
        if isinstance(other, Vector):
            a, b, c = self._ar
            d, e, f = other._ar
            c1 = numpy.linalg.det(numpy.array(((b, c), (e, f))))
            c2 = -numpy.linalg.det(numpy.array(((a, c), (d, f))))
            c3 = numpy.linalg.det(numpy.array(((a, b), (d, e))))
            return Vector(c1, c2, c3)
        else:
            a = self._ar * numpy.array(other)
            return Vector(a)

    def __getitem__(self, i):
        """Return value of array index i."""
        return self._ar[i]

    def __setitem__(self, i, value):
        """Assign values to array index i."""
        self._ar[i] = value

    def __contains__(self, i):
        """Validate if i is in array."""
        return i in self._ar

    def norm(self):
        """Return vector norm."""
        return numpy.sqrt(sum(self._ar * self._ar))

    def normsq(self):
        """Return square of vector norm."""
        return abs(sum(self._ar * self._ar))

    def normalize(self):
        """Normalize the Vector object.

        Changes the state of ``self`` and doesn't return a value.
        If you need to chain function calls or create a new object
        use the ``normalized`` method.
        """
        if self.norm():
            self._ar = self._ar / self.norm()

    def normalized(self):
        """Return a normalized copy of the Vector.

        To avoid allocating new objects use the ``normalize`` method.
        """
        v = self.copy()
        v.normalize()
        return v

    def angle(self, other):
        """Return angle between two vectors."""
        n1 = self.norm()
        n2 = other.norm()
        c = (self * other) / (n1 * n2)
        # Take care of roundoff errors
        c = min(c, 1)
        c = max(-1, c)
        return numpy.arccos(c)

    def get_array(self):
        """Return (a copy of) the array of coordinates."""
        return numpy.array(self._ar)

    def left_multiply(self, matrix):
        """Return Vector=Matrix x Vector."""
        a = numpy.dot(matrix, self._ar)
        return Vector(a)

    def right_multiply(self, matrix):
        """Return Vector=Vector x Matrix."""
        a = numpy.dot(self._ar, matrix)
        return Vector(a)

    def copy(self):
        """Return a deep copy of the Vector."""
        return Vector(self._ar)


"""Homogeneous matrix geometry routines.

Rotation, translation, scale, and coordinate transformations.

Robert T. Miller 2019
"""


def homog_rot_mtx(angle_rads: float, axis: str) -> numpy.array:
    """Generate a 4x4 single-axis numpy rotation matrix.

    :param float angle_rads: the desired rotation angle in radians
    :param char axis: character specifying the rotation axis
    """
    cosang = numpy.cos(angle_rads)
    sinang = numpy.sin(angle_rads)

    if "z" == axis:
        return numpy.array(
            (
                (cosang, -sinang, 0, 0),
                (sinang, cosang, 0, 0),
                (0, 0, 1, 0),
                (0, 0, 0, 1),
            ),
            dtype=numpy.float64,
        )
    elif "y" == axis:
        return numpy.array(
            (
                (cosang, 0, sinang, 0),
                (0, 1, 0, 0),
                (-sinang, 0, cosang, 0),
                (0, 0, 0, 1),
            ),
            dtype=numpy.float64,
        )
    else:
        return numpy.array(
            (
                (1, 0, 0, 0),
                (0, cosang, -sinang, 0),
                (0, sinang, cosang, 0),
                (0, 0, 0, 1),
            ),
            dtype=numpy.float64,
        )


def set_Z_homog_rot_mtx(angle_rads: float, mtx: numpy.ndarray):
    """Update existing Z rotation matrix to new angle."""
    cosang = numpy.cos(angle_rads)
    sinang = numpy.sin(angle_rads)

    mtx[0][0] = mtx[1][1] = cosang
    mtx[1][0] = sinang
    mtx[0][1] = -sinang


def set_Y_homog_rot_mtx(angle_rads: float, mtx: numpy.ndarray):
    """Update existing Y rotation matrix to new angle."""
    cosang = numpy.cos(angle_rads)
    sinang = numpy.sin(angle_rads)

    mtx[0][0] = mtx[2][2] = cosang
    mtx[0][2] = sinang
    mtx[2][0] = -sinang


def set_X_homog_rot_mtx(angle_rads: float, mtx: numpy.ndarray):
    """Update existing X rotation matrix to new angle."""
    cosang = numpy.cos(angle_rads)
    sinang = numpy.sin(angle_rads)

    mtx[1][1] = mtx[2][2] = cosang
    mtx[2][1] = sinang
    mtx[1][2] = -sinang


def homog_trans_mtx(x: float, y: float, z: float) -> numpy.array:
    """Generate a 4x4 numpy translation matrix.

    :param x, y, z: translation in each axis
    """
    return numpy.array(
        ((1, 0, 0, x), (0, 1, 0, y), (0, 0, 1, z), (0, 0, 0, 1)), dtype=numpy.float64
    )


def set_homog_trans_mtx(x: float, y: float, z: float, mtx: numpy.ndarray):
    """Update existing translation matrix to new values."""
    mtx[0][3] = x
    mtx[1][3] = y
    mtx[2][3] = z


def homog_scale_mtx(scale: float) -> numpy.array:
    """Generate a 4x4 numpy scaling matrix.

    :param float scale: scale multiplier
    """
    return numpy.array(
        [[scale, 0, 0, 0], [0, scale, 0, 0], [0, 0, scale, 0], [0, 0, 0, 1]],
        dtype=numpy.float64,
    )


def _get_azimuth(x: float, y: float) -> float:
    sign_y = -1.0 if y < 0.0 else 1.0
    sign_x = -1.0 if x < 0.0 else 1.0
    return (
        numpy.arctan2(y, x)
        if (0 != x and 0 != y)
        else (numpy.pi / 2.0 * sign_y)  # +/-90 if X=0, Y!=0
        if 0 != y
        else numpy.pi
        if sign_x < 0.0  # 180 if Y=0, X < 0
        else 0.0  # 0 if Y=0, X >= 0
    )


def get_spherical_coordinates(xyz: numpy.array) -> Tuple[float, float, float]:
    """Compute spherical coordinates (r, azimuth, polar_angle) for X,Y,Z point.

    :param array xyz: column vector (3 row x 1 column numpy array)
    :return: tuple of r, azimuth, polar_angle for input coordinate
    """
    r = numpy.linalg.norm(xyz)
    if 0 == r:
        return (0, 0, 0)
    azimuth = _get_azimuth(xyz[0], xyz[1])
    polar_angle = numpy.arccos(xyz[2] / r)

    return (r, azimuth, polar_angle)


gtm = numpy.identity(4, dtype=numpy.float64)
gmrz = numpy.identity(4, dtype=numpy.float64)
gmry = numpy.identity(4, dtype=numpy.float64)
gmrz2 = numpy.identity(4, dtype=numpy.float64)


def coord_space(
    a0: numpy.ndarray, a1: numpy.ndarray, a2: numpy.ndarray, rev: bool = False
) -> Tuple[numpy.ndarray, Optional[numpy.ndarray]]:
    """Generate transformation matrix to coordinate space defined by 3 points.

    New coordinate space will have:
        acs[0] on XZ plane
        acs[1] origin
        acs[2] on +Z axis

    :param numpy column array x3 acs: X,Y,Z column input coordinates x3
    :param bool rev: if True, also return reverse transformation matrix
        (to return from coord_space)
    :returns: 4x4 numpy array, x2 if rev=True
    """
    # dbg = False
    # if dbg:
    #    print(a0.transpose())
    #    print(a1.transpose())
    #    print(a2.transpose())

    # a0 = acs[0]
    # a1 = acs[1]
    # a2 = acs[2]

    global gtm
    global gmry
    global gmrz, gmrz2

    tm = gtm
    mry = gmry
    mrz = gmrz
    mrz2 = gmrz2

    # tx acs[1] to origin
    # tm = homog_trans_mtx(-a1[0][0], -a1[1][0], -a1[2][0])
    set_homog_trans_mtx(-a1[0], -a1[1], -a1[2], tm)

    # directly translate a2 using a1
    p = a2 - a1
    sc = get_spherical_coordinates(p)

    # if dbg:
    #    print("p", p.transpose())
    #    print("sc", sc)

    # mrz = homog_rot_mtx(-sc[1], "z")  # rotate translated a2 -azimuth about Z
    set_Z_homog_rot_mtx(-sc[1], mrz)
    # mry = homog_rot_mtx(-sc[2], "y")  # rotate translated a2 -polar_angle about Y
    set_Y_homog_rot_mtx(-sc[2], mry)

    # mt completes a1-a2 on Z-axis, still need to align a0 with XZ plane
    # mt = mry @ mrz @ tm  # python 3.5 and later
    mt = gmry.dot(gmrz.dot(gtm))

    # if dbg:
    #    print("tm:\n", tm)
    #    print("mrz:\n", mrz)
    #    print("mry:\n", mry)
    #    # print("mt ", mt)

    p = mt.dot(a0)

    # if dbg:
    #    print("mt:\n", mt, "\na0:\n", a0, "\np:\n", p)

    # need azimuth of translated a0
    # sc2 = get_spherical_coordinates(p)
    # print(sc2)
    azimuth2 = _get_azimuth(p[0], p[1])

    # rotate a0 -azimuth2 about Z to align with X
    # mrz2 = homog_rot_mtx(-azimuth2, "z")
    set_Z_homog_rot_mtx(-azimuth2, mrz2)

    # mt = mrz2 @ mt
    mt = gmrz2.dot(mt)

    # if dbg:
    #    print("mt:", mt, "\na0:", a0, "\np:", p)
    #    # print(p, "\n", azimuth2, "\n", mrz2, "\n", mt)

    # if dbg:
    #    print("mt:\n", mt)
    #    print("<<<<<<==============================")

    if not rev:
        return mt, None

    # rev=True, so generate the reverse transformation

    # rotate a0 theta about Z, reversing alignment with X
    # mrz2 = homog_rot_mtx(azimuth2, "z")
    set_Z_homog_rot_mtx(azimuth2, mrz2)
    # rotate a2 phi about Y
    # mry = homog_rot_mtx(sc[2], "y")
    set_Y_homog_rot_mtx(sc[2], mry)
    # rotate a2 theta about Z
    # mrz = homog_rot_mtx(sc[1], "z")
    set_Z_homog_rot_mtx(sc[1], mrz)
    # translation matrix origin to a1
    # tm = homog_trans_mtx(a1[0][0], a1[1][0], a1[2][0])
    set_homog_trans_mtx(a1[0], a1[1], a1[2], tm)

    # mr = tm @ mrz @ mry @ mrz2
    mr = gtm.dot(gmrz.dot(gmry.dot(gmrz2)))
    # mr = numpy.dot(tm, numpy.dot(mrz, numpy.dot(mry, mrz2)))

    return mt, mr


def multi_rot_Z(angle_rads: numpy.ndarray) -> numpy.ndarray:
    """Create [entries] numpy Z rotation matrices for [entries] angles.

    :param entries: int number of matrices generated.
    :param angle_rads: numpy array of angles
    :returns: entries x 4 x 4 homogeneous rotation matrices
    """
    rz = numpy.empty((angle_rads.shape[0], 4, 4))
    rz[...] = numpy.identity(4)
    rz[:, 0, 0] = rz[:, 1, 1] = numpy.cos(angle_rads)
    rz[:, 1, 0] = numpy.sin(angle_rads)
    rz[:, 0, 1] = -rz[:, 1, 0]
    return rz


def multi_rot_Y(angle_rads: numpy.ndarray) -> numpy.ndarray:
    """Create [entries] numpy Y rotation matrices for [entries] angles.

    :param entries: int number of matrices generated.
    :param angle_rads: numpy array of angles
    :returns: entries x 4 x 4 homogeneous rotation matrices
    """
    ry = numpy.empty((angle_rads.shape[0], 4, 4))
    ry[...] = numpy.identity(4)
    ry[:, 0, 0] = ry[:, 2, 2] = numpy.cos(angle_rads)
    ry[:, 0, 2] = numpy.sin(angle_rads)
    ry[:, 2, 0] = -ry[:, 0, 2]

    return ry
