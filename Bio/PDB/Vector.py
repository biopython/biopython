# Copyright (C) 2004, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Numeric import array, sum, sqrt, arccos, matrixmultiply, transpose, cos, sin, zeros
from LinearAlgebra import determinant
from MLab import eye
import sys


def rotaxis(theta, vector):
    """
    Left multiplying rotation matrix that rotates
    theta along vector (which should be a unit 
    vector).

    ie. m=rotaxis(pi, Vector(1,0,0))
        rotated_vector=any_vector.left_multiply(m)
    """
    vector=vector.normalize()
    c=cos(theta)
    s=sin(theta)
    t=1-c
    x,y,z=vector.get_array()
    rot=zeros((3,3), "d")
    # 1st row
    rot[0,0]=t*x*x+c
    rot[0,1]=t*x*y-s*z
    rot[0,2]=t*x*z+s*y
    # 2nd row
    rot[1,0]=t*x*y+s*z
    rot[1,1]=t*y*y+c
    rot[1,2]=t*y*z-s*x
    # 3rd row
    rot[2,0]=t*x*z-s*y
    rot[2,1]=t*y*z+s*x
    rot[2,2]=t*z*z+c
    return rot

def refmat(p,q):
    """
    Return a (left multiplying) matrix that mirrors p onto q.
    """
    p=p.normalize()
    q=q.normalize()
    pq=p-q
    pq=pq.normalize()
    b=pq.get_array()
    b.shape=(3, 1)
    i=eye(3)
    ref=i-2*matrixmultiply(b, transpose(b))
    return ref

def rotmat(p,q):
    """
    Return a (left multiplying) matrix that rotates p onto q.
    """
    rot=matrixmultiply(refmat(q, -p), refmat(p, -p))
    return rot

def calc_angle(v1, v2, v3):
    """
    Calculate the angle between 3 vectors
    representing 3 connected points.
    """
    v1=v1-v2
    v3=v3-v2
    return v1.angle(v3)

def calc_dihedral(v1, v2, v3, v4):
    """
    Calculate the dihedral angle between 4 vectors
    representing 4 connected points. The angle is in
    ]-pi, pi].
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

class Vector:
    "3D vector"

    def __init__(self, x, y, z):
        self._ar=array((x, y, z), 'd')

    def __neg__(self):
        "Return Vector(-x, -y, -z)"
        x,y,z=-self._ar
        return Vector(x,y,z)

    def __add__(self, other):
        "Return Vector+other Vector or scalar"
        if not isinstance(other, Vector):
            x,y,z=self._ar+other
        else:
            x,y,z=self._ar+other._ar
        return Vector(x,y,z)

    def __sub__(self, other):
        "Return Vector-other Vector or scalar"
        if not isinstance(other, Vector):
            x,y,z=self._ar-other
        else:
            x,y,z=self._ar-other._ar
        return Vector(x,y,z)

    def __mul__(self, other):
        "Return Vector . Vector (dot product) or Vector*scalar"
        if not isinstance(other, Vector):
            x,y,z=self._ar*other
            return Vector(x,y,z)
        else:
            return sum(self._ar*other._ar)

    def __div__(self, a):
        "Return Vector(coords/a)"
        x,y,z=self._ar/a
        return Vector(x,y,z)

    def __pow__(self, other):
        "Return Vector x Vector (cross product)"
        a,b,c=self._ar
        d,e,f=other._ar
        c1=determinant(array(((b,c), (e,f))))
        c2=determinant(array(((a,c), (d,f))))
        c3=determinant(array(((a,b), (d,e))))
        return Vector(c1,c2,c3)

    def __str__(self):
        x,y,z=self._ar
        return "<Vector %.2f %.2f %.2f>" % (x,y,z)

    def norm(self):
        "Return vector norm"
        return sqrt(sum(self._ar*self._ar))

    def normalize(self):
        "Return normalized Vector"
        x,y,z=self._ar/self.norm()
        return Vector(x,y,z)

    def angle(self, other):
        "Return angle between two vectors"
        n1=self.norm()
        n2=other.norm()
        c=(self*other)/(n1*n2)
        # Take care of roundoff errors
        c=min(c,1)
        c=max(-1,c)
        return arccos(c)

    def get_array(self):
        "Return (a copy of) the array of coordinates"
        return array(self._ar)

    def left_multiply(self, matrix):
        "Return Vector=Matrix x Vector"
        x,y,z=matrixmultiply(matrix, self._ar)
        return Vector(x,y,z)

    def right_multiply(self, matrix):
        "Return Vector=Vector x Matrix"
        x,y,z=matrixmultiply(self._ar, matrix)
        return Vector(x,y,z)

    def copy(self):
        "Return a deep copy of the Vector"
        x,y,z=self._ar
        return Vector(x,y,z)

if __name__=="__main__":

        from math import pi
        from RandomArray import *

        v1=Vector(0.5,0.5,1.3)
        v2=Vector(0.1,0.1,0.1)
        v3=Vector(1.9,0.8,0.6)
        v4=Vector(1,-1,0)

        calc_angle(v1, v2, v3)
        calc_dihedral(v1, v2, v3, v4)

        ref=refmat(v1, v3)
        rot=rotmat(v1, v3)

        print v3
        print v1.left_multiply(ref)
        print v1.left_multiply(rot)
        print v1.right_multiply(transpose(rot))

