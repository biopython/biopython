# Copyright (C) 2004, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Vector class, including rotation-related functions."""

from __future__ import print_function

import numpy


def m2rotaxis(m):
    """Return angles, axis pair that corresponds to rotation matrix m.

    The case where `m` is the identity matrix corresponds to a singularity where any
    rotation axis is valid. In that case, `Vector([1,0,0])`, is returned.
    """
    eps = 1e-5

    # Check for singularities a la http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToAngle/
    if abs(m[0, 1] - m[1, 0]) < eps and abs(m[0, 2] - m[2, 0]) < eps and abs(m[1, 2] - m[2, 1]) < eps:
        # Singularity encountered. Check if its 0 or 180 deg
        if abs(m[0, 1] + m[1, 0]) < eps and abs(m[0, 2] + m[2, 0]) < eps and abs(m[1, 2] + m[2, 1]) < eps and abs(m[0, 0] + m[1, 1] + m[2, 2] - 3) < eps:
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
    >>> m = rotaxis(pi, Vector(1, 0, 0))
    >>> rotated_vector = any_vector.left_multiply(m)

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
    >>> mirror = refmat(p, q)
    >>> qq = p.left_multiply(mirror)
    >>> print(q)
    >>> print(qq)  # q and qq should be the same

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
    >>> r = rotmat(p, q)
    >>> print(q)
    >>> print(p.left_multiply(r))

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


class Vector(object):
    """3D vector."""

    def __init__(self, x, y=None, z=None):
        """Initialize the class."""
        if y is None and z is None:
            # Array, list, tuple...
            if len(x) != 3:
                raise ValueError("Vector: x is not a "
                                 "list/tuple/array of 3 numbers")
            self._ar = numpy.array(x, 'd')
        else:
            # Three numbers
            self._ar = numpy.array((x, y, z), 'd')

    def __repr__(self):
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

    def __div__(self, x):
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
        return self._ar[i]

    def __setitem__(self, i, value):
        self._ar[i] = value

    def __contains__(self, i):
        return i in self._ar

    def norm(self):
        """Return vector norm."""
        return numpy.sqrt(sum(self._ar * self._ar))

    def normsq(self):
        """Return square of vector norm."""
        return abs(sum(self._ar * self._ar))

    def normalize(self):
        """Normalize the Vector object.

        Changes the state of `self` and doesn't return a value. If you need to chain function
        calls or create a new object use the `normalized` method.
        """
        if self.norm():
            self._ar = self._ar / self.norm()

    def normalized(self):
        """Return a normalized copy of the Vector.

        To avoid allocating new objects use the `normalize` method.
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
