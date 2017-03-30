# Copyright (C) 2004, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Vector class, including rotation-related functions."""

from __future__ import print_function

import numpy


class Vector(object):
    "3D vector"

    def __init__(self, x, y=None, z=None):
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
        "Return Vector(-x, -y, -z)"
        a = -self._ar
        return Vector(a)

    def __add__(self, other):
        "Return Vector+other Vector or scalar"
        if isinstance(other, Vector):
            a = self._ar + other._ar
        else:
            a = self._ar + numpy.array(other)
        return Vector(a)

    def __sub__(self, other):
        "Return Vector-other Vector or scalar"
        if isinstance(other, Vector):
            a = self._ar - other._ar
        else:
            a = self._ar - numpy.array(other)
        return Vector(a)

    def __mul__(self, other):
        "Return Vector.Vector (dot product)"
        return sum(self._ar * other._ar)

    def __div__(self, x):
        "Return Vector(coords/a)"
        a = self._ar / numpy.array(x)
        return Vector(a)

    def __pow__(self, other):
        "Return VectorxVector (cross product) or Vectorxscalar"
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
        return (i in self._ar)

    def norm(self):
        "Return vector norm"
        return numpy.sqrt(sum(self._ar * self._ar))

    def normsq(self):
        "Return square of vector norm"
        return abs(sum(self._ar * self._ar))

    def normalize(self):
        "Normalize the Vector"
        self._ar = self._ar / self.norm()

    def normalized(self):
        "Return a normalized copy of the Vector"
        v = self.copy()
        v.normalize()
        return v

    def angle(self, other):
        "Return angle between two vectors"
        n1 = self.norm()
        n2 = other.norm()
        c = (self * other) / (n1 * n2)
        # Take care of roundoff errors
        c = min(c, 1)
        c = max(-1, c)
        return numpy.arccos(c)

    def get_array(self):
        "Return (a copy of) the array of coordinates"
        return numpy.array(self._ar)

    def left_multiply(self, matrix):
        "Return Vector=Matrix x Vector"
        a = numpy.dot(matrix, self._ar)
        return Vector(a)

    def right_multiply(self, matrix):
        "Return Vector=Vector x Matrix"
        a = numpy.dot(self._ar, matrix)
        return Vector(a)

    def copy(self):
        "Return a deep copy of the Vector"
        return Vector(self._ar)
