# Copyright (C) 2004, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Vector class, including rotation-related functions."""

from __future__ import print_function

import numpy


def m2rotaxis(m):
    """
    Return angles, axis pair that corresponds to rotation matrix m.
    """
    # Angle always between 0 and pi
    # Sense of rotation is defined by axis orientation
    t=0.5*(numpy.trace(m)-1)
    t=max(-1, t)
    t=min(1, t)
    angle=numpy.arccos(t)
    if angle<1e-15:
        # Angle is 0
        return 0.0, Vector(1, 0, 0)
    elif angle<numpy.pi:
        # Angle is smaller than pi
        x=m[2, 1]-m[1, 2]
        y=m[0, 2]-m[2, 0]
        z=m[1, 0]-m[0, 1]
        axis=Vector(x, y, z)
        axis.normalize()
        return angle, axis
    else:
        # Angle is pi - special case!
        m00=m[0, 0]
        m11=m[1, 1]
        m22=m[2, 2]
        if m00>m11 and m00>m22:
            x=numpy.sqrt(m00-m11-m22+0.5)
            y=m[0, 1]/(2*x)
            z=m[0, 2]/(2*x)
        elif m11>m00 and m11>m22:
            y=numpy.sqrt(m11-m00-m22+0.5)
            x=m[0, 1]/(2*y)
            z=m[1, 2]/(2*y)
        else:
            z=numpy.sqrt(m22-m00-m11+0.5)
            x=m[0, 2]/(2*z)
            y=m[1, 2]/(2*z)
        axis=Vector(x, y, z)
        axis.normalize()
        return numpy.pi, axis


def vector_to_axis(line, point):
    """
    Returns the vector between a point and
    the closest point on a line (ie. the perpendicular
    projection of the point on the line).

    @type line: L{Vector}
    @param line: vector defining a line

    @type point: L{Vector}
    @param point: vector defining the point
    """
    line=line.normalized()
    np=point.norm()
    angle=line.angle(point)
    return point-line**(np*numpy.cos(angle))


def rotaxis2m(theta, vector):
    """
    Calculate a left multiplying rotation matrix that rotates
    theta rad around vector.

    Example:

        >>> m=rotaxis(pi, Vector(1, 0, 0))
        >>> rotated_vector=any_vector.left_multiply(m)

    @type theta: float
    @param theta: the rotation angle


    @type vector: L{Vector}
    @param vector: the rotation axis

    @return: The rotation matrix, a 3x3 Numeric array.
    """
    vector=vector.copy()
    vector.normalize()
    c=numpy.cos(theta)
    s=numpy.sin(theta)
    t=1-c
    x, y, z=vector.get_array()
    rot=numpy.zeros((3, 3))
    # 1st row
    rot[0, 0]=t*x*x+c
    rot[0, 1]=t*x*y-s*z
    rot[0, 2]=t*x*z+s*y
    # 2nd row
    rot[1, 0]=t*x*y+s*z
    rot[1, 1]=t*y*y+c
    rot[1, 2]=t*y*z-s*x
    # 3rd row
    rot[2, 0]=t*x*z-s*y
    rot[2, 1]=t*y*z+s*x
    rot[2, 2]=t*z*z+c
    return rot

rotaxis=rotaxis2m


def refmat(p, q):
    """
    Return a (left multiplying) matrix that mirrors p onto q.

    Example:
        >>> mirror=refmat(p, q)
        >>> qq=p.left_multiply(mirror)
        >>> print(q)
        >>> print(qq) # q and qq should be the same

    @type p,q: L{Vector}
    @return: The mirror operation, a 3x3 Numeric array.
    """
    p.normalize()
    q.normalize()
    if (p-q).norm()<1e-5:
        return numpy.identity(3)
    pq=p-q
    pq.normalize()
    b=pq.get_array()
    b.shape=(3, 1)
    i=numpy.identity(3)
    ref=i-2*numpy.dot(b, numpy.transpose(b))
    return ref


def rotmat(p, q):
    """
    Return a (left multiplying) matrix that rotates p onto q.

    Example:
        >>> r=rotmat(p, q)
        >>> print(q)
        >>> print(p.left_multiply(r))

    @param p: moving vector
    @type p: L{Vector}

    @param q: fixed vector
    @type q: L{Vector}

    @return: rotation matrix that rotates p onto q
    @rtype: 3x3 Numeric array
    """
    rot=numpy.dot(refmat(q, -p), refmat(p, -p))
    return rot


def calc_angle(v1, v2, v3):
    """
    Calculate the angle between 3 vectors
    representing 3 connected points.

    @param v1, v2, v3: the tree points that define the angle
    @type v1, v2, v3: L{Vector}

    @return: angle
    @rtype: float
    """
    v1=v1-v2
    v3=v3-v2
    return v1.angle(v3)


def calc_dihedral(v1, v2, v3, v4):
    """
    Calculate the dihedral angle between 4 vectors
    representing 4 connected points. The angle is in
    ]-pi, pi].

    @param v1, v2, v3, v4: the four points that define the dihedral angle
    @type v1, v2, v3, v4: L{Vector}
    """
    ab=v1-v2
    cb=v3-v2
    db=v4-v3
    u=ab**cb
    v=db**cb
    w=u**v
    angle=u.angle(v)
    # Determine sign of angle
    try:
        if cb.angle(w)>0.001:
            angle=-angle
    except ZeroDivisionError:
        # dihedral=pi
        pass
    return angle


class Vector(object):
    "3D vector"

    def __init__(self, x, y=None, z=None):
        if y is None and z is None:
            # Array, list, tuple...
            if len(x)!=3:
                raise ValueError("Vector: x is not a "
                                 "list/tuple/array of 3 numbers")
            self._ar=numpy.array(x, 'd')
        else:
            # Three numbers
            self._ar=numpy.array((x, y, z), 'd')

    def __repr__(self):
        x, y, z=self._ar
        return "<Vector %.2f, %.2f, %.2f>" % (x, y, z)

    def __neg__(self):
        "Return Vector(-x, -y, -z)"
        a=-self._ar
        return Vector(a)

    def __add__(self, other):
        "Return Vector+other Vector or scalar"
        if isinstance(other, Vector):
            a=self._ar+other._ar
        else:
            a=self._ar+numpy.array(other)
        return Vector(a)

    def __sub__(self, other):
        "Return Vector-other Vector or scalar"
        if isinstance(other, Vector):
            a=self._ar-other._ar
        else:
            a=self._ar-numpy.array(other)
        return Vector(a)

    def __mul__(self, other):
        "Return Vector.Vector (dot product)"
        return sum(self._ar*other._ar)

    def __div__(self, x):
        "Return Vector(coords/a)"
        a=self._ar/numpy.array(x)
        return Vector(a)

    def __pow__(self, other):
        "Return VectorxVector (cross product) or Vectorxscalar"
        if isinstance(other, Vector):
            a, b, c=self._ar
            d, e, f=other._ar
            c1=numpy.linalg.det(numpy.array(((b, c), (e, f))))
            c2=-numpy.linalg.det(numpy.array(((a, c), (d, f))))
            c3=numpy.linalg.det(numpy.array(((a, b), (d, e))))
            return Vector(c1, c2, c3)
        else:
            a=self._ar*numpy.array(other)
            return Vector(a)

    def __getitem__(self, i):
        return self._ar[i]

    def __setitem__(self, i, value):
        self._ar[i]=value

    def __contains__(self, i):
        return (i in self._ar)

    def norm(self):
        "Return vector norm"
        return numpy.sqrt(sum(self._ar*self._ar))

    def normsq(self):
        "Return square of vector norm"
        return abs(sum(self._ar*self._ar))

    def normalize(self):
        "Normalize the Vector"
        self._ar=self._ar/self.norm()

    def normalized(self):
        "Return a normalized copy of the Vector"
        v=self.copy()
        v.normalize()
        return v

    def angle(self, other):
        "Return angle between two vectors"
        n1=self.norm()
        n2=other.norm()
        c=(self*other)/(n1*n2)
        # Take care of roundoff errors
        c=min(c, 1)
        c=max(-1, c)
        return numpy.arccos(c)

    def get_array(self):
        "Return (a copy of) the array of coordinates"
        return numpy.array(self._ar)

    def left_multiply(self, matrix):
        "Return Vector=Matrix x Vector"
        a=numpy.dot(matrix, self._ar)
        return Vector(a)

    def right_multiply(self, matrix):
        "Return Vector=Vector x Matrix"
        a=numpy.dot(self._ar, matrix)
        return Vector(a)

    def copy(self):
        "Return a deep copy of the Vector"
        return Vector(self._ar)

if __name__=="__main__":

    from numpy.random import random

    v1=Vector(0, 0, 1)
    v2=Vector(0, 0, 0)
    v3=Vector(0, 1, 0)
    v4=Vector(1, 1, 0)

    v4.normalize()

    print(v4)

    print(calc_angle(v1, v2, v3))
    dih=calc_dihedral(v1, v2, v3, v4)
    # Test dihedral sign
    assert(dih>0)
    print("DIHEDRAL %f" % dih)

    ref=refmat(v1, v3)
    rot=rotmat(v1, v3)

    print(v3)
    print(v1.left_multiply(ref))
    print(v1.left_multiply(rot))
    print(v1.right_multiply(numpy.transpose(rot)))

    # -
    print(v1-v2)
    print(v1-1)
    print(v1+(1, 2, 3))
    # +
    print(v1+v2)
    print(v1+3)
    print(v1-(1, 2, 3))
    # *
    print(v1*v2)
    # /
    print(v1/2)
    print(v1/(1, 2, 3))
    # **
    print(v1**v2)
    print(v1**2)
    print(v1**(1, 2, 3))
    # norm
    print(v1.norm())
    # norm squared
    print(v1.normsq())
    # setitem
    v1[2]=10
    print(v1)
    # getitem
    print(v1[2])

    print(numpy.array(v1))

    print("ROT")

    angle=random()*numpy.pi
    axis=Vector(random(3)-random(3))
    axis.normalize()

    m=rotaxis(angle, axis)

    cangle, caxis=m2rotaxis(m)

    print(angle-cangle)
    print(axis-caxis)
    print("")
