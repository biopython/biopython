"""Python implementation of chisqprob, to avoid SciPy dependency.

Adapted from SciPy: scipy/special/cephes/{chdtr,igam}.
"""

import math

try:
    from math import lgamma as _lgamma
except ImportError:
    # Missing in Python 2.6, using Python Python recipe from
    # http://code.activestate.com/recipes/576393/
    # with PEP8 whitespace and minor import changes
    #
    # * ====================================================
    # * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
    # *
    # * Developed at SunPro, a Sun Microsystems, Inc. business.
    # * Permission to use, copy, modify, and distribute this
    # * software is freely granted, provided that this notice
    # * is preserved.
    # * ====================================================
    two52 = 4.50359962737049600000e+15
    half = 5.00000000000000000000e-01
    one = 1.00000000000000000000e+00
    pi = 3.14159265358979311600e+00
    a0 = 7.72156649015328655494e-02
    a1 = 3.22467033424113591611e-01
    a2 = 6.73523010531292681824e-02
    a3 = 2.05808084325167332806e-02
    a4 = 7.38555086081402883957e-03
    a5 = 2.89051383673415629091e-03
    a6 = 1.19270763183362067845e-03
    a7 = 5.10069792153511336608e-04
    a8 = 2.20862790713908385557e-04
    a9 = 1.08011567247583939954e-04
    a10 = 2.52144565451257326939e-05
    a11 = 4.48640949618915160150e-05
    tc = 1.46163214496836224576e+00
    tf = -1.21486290535849611461e-01
    # /* tt = -(tail of tf) */
    tt = -3.63867699703950536541e-18
    t0 = 4.83836122723810047042e-01
    t1 = -1.47587722994593911752e-01
    t2 = 6.46249402391333854778e-02
    t3 = -3.27885410759859649565e-02
    t4 = 1.79706750811820387126e-02
    t5 = -1.03142241298341437450e-02
    t6 = 6.10053870246291332635e-03
    t7 = -3.68452016781138256760e-03
    t8 = 2.25964780900612472250e-03
    t9 = -1.40346469989232843813e-03
    t10 = 8.81081882437654011382e-04
    t11 = -5.38595305356740546715e-04
    t12 = 3.15632070903625950361e-04
    t13 = -3.12754168375120860518e-04
    t14 = 3.35529192635519073543e-04
    u0 = -7.72156649015328655494e-02
    u1 = 6.32827064025093366517e-01
    u2 = 1.45492250137234768737e+00
    u3 = 9.77717527963372745603e-01
    u4 = 2.28963728064692451092e-01
    u5 = 1.33810918536787660377e-02
    v1 = 2.45597793713041134822e+00
    v2 = 2.12848976379893395361e+00
    v3 = 7.69285150456672783825e-01
    v4 = 1.04222645593369134254e-01
    v5 = 3.21709242282423911810e-03
    s0 = -7.72156649015328655494e-02
    s1 = 2.14982415960608852501e-01
    s2 = 3.25778796408930981787e-01
    s3 = 1.46350472652464452805e-01
    s4 = 2.66422703033638609560e-02
    s5 = 1.84028451407337715652e-03
    s6 = 3.19475326584100867617e-05
    r1 = 1.39200533467621045958e+00
    r2 = 7.21935547567138069525e-01
    r3 = 1.71933865632803078993e-01
    r4 = 1.86459191715652901344e-02
    r5 = 7.77942496381893596434e-04
    r6 = 7.32668430744625636189e-06
    w0 = 4.18938533204672725052e-01
    w1 = 8.33333333333329678849e-02
    w2 = -2.77777777728775536470e-03
    w3 = 7.93650558643019558500e-04
    w4 = -5.95187557450339963135e-04
    w5 = 8.36339918996282139126e-04
    w6 = -1.63092934096575273989e-03
    zero = 0.00000000000000000000e+00

    # inf = float('inf')
    # nan = float('nan')
    inf = float(9e999)

    def _sin_pi(x):
        x = float(x)
        e, ix = math.frexp(x)
        if(abs(x) < 0.25):
            return -math.sin(pi * x)
        y = -x  # /* x is assume negative */

        # * argument reduction, make sure inexact flag not raised if input
        # * is an integer
        z = floor(y)
        if(z != y):
            y *= 0.5
            y = 2.0 * (y - floor(y))  # /* y = |x| mod 2.0 */
            n = int(y * 4.0)
        else:
            if(abs(ix) >= 53):
                y = zero
                n = 0  # /* y must be even */
            else:
                if(abs(ix) < 52):
                    z = y + two52  # /* exact */
                e, n = math.frexp(z)
                n &= 1
                y = n
                n <<= 2

        if n == 0:
            y = sin(pi * y)
        elif (n == 1 or n == 2):
            y = cos(pi * (0.5 - y))
        elif (n == 3 or n == 4):
            y = sin(pi * (one - y))
        elif (n == 5 or n == 6):
            y = -cos(pi * (y - 1.5))
        else:
            y = sin(pi * (y - 2.0))

        z = cos(pi * (z + 1.0))
        return -y * z

    def _lgamma(x):
        """Natural logarithm of gamma function of x

        raise ValueError if x is negative integer."""
        x = float(x)

        # /* purge off +-inf, NaN, +-0, and negative arguments */
        if ((x == inf) or (x == -inf)):
            return inf
        # if (x is nan):
        # return nan

        e, ix = math.frexp(x)
        nadj = 0
        signgamp = 1

        if ((e == 0.0) and (ix == 0)):
            return inf

        if (ix > 1020):
            return inf

        if ((e != 0.0) and (ix < -71)):
            if (x < 0):
                return -math.log(-x)
            else:
                return -math.log(x)

        if e < 0:
            if ix > 52:
                return inf  # one/zero
            t = sin_pi(x)
            if t == zero:
                # return inf
                raise ValueError('gamma not defined for negative integer')
            nadj = math.log(pi / fabs(t * x))
            if t < zero:
                signgamp = -1
            x = -x

        # /* purge off 1 and 2 */
        if x == 2.0 or x == 1.0:
            r = 0.0

        # /* for x < 2.0 */
        elif ix < 2:
            if x <= 0.9:  # /* lgamma(x) = lgamma(x+1)-log(x) */
                r = -math.log(x)
                if x >= 0.7316:
                    y = one - x
                    z = y * y
                    p1 = a0 + z * (a2 + z * (a4 + z * (a6 + z * (a8 + z * a10))))
                    p2 = z * (a1 + z * (a3 + z * (a5 + z * (a7 + z * (a9 + z * a11)))))
                    p = y * p1 + p2
                    r += (p - 0.5 * y)
                elif (x >= 0.23164):
                    y = x - (tc - one)
                    z = y * y
                    w = z * y
                    p1 = t0 + w * (t3 + w * (t6 + w * (t9 + w * t12)))  # /* parallel comp */
                    p2 = t1 + w * (t4 + w * (t7 + w * (t10 + w * t13)))
                    p3 = t2 + w * (t5 + w * (t8 + w * (t11 + w * t14)))
                    p = z * p1 - (tt - w * (p2 + y * p3))
                    r += (tf + p)
                else:
                    y = x
                    p1 = y * (u0 + y * (u1 + y * (u2 + y * (u3 + y * (u4 + y * u5)))))
                    p2 = one + y * (v1 + y * (v2 + y * (v3 + y * (v4 + y * v5))))
                    r += (-0.5 * y + p1 / p2)
            else:
                r = zero
                if(x >= 1.7316):
                    y = 2.0 - x  # /* [1.7316,2] */
                    z = y * y
                    p1 = a0 + z * (a2 + z * (a4 + z * (a6 + z * (a8 + z * a10))))
                    p2 = z * (a1 + z * (a3 + z * (a5 + z * (a7 + z * (a9 + z * a11)))))
                    p = y * p1 + p2
                    r += (p - 0.5 * y)
                elif(x >= 1.23164):
                    y = x - tc  # /* [1.23,1.73] */
                    z = y * y
                    w = z * y
                    p1 = t0 + w * (t3 + w * (t6 + w * (t9 + w * t12)))  # /* parallel comp */
                    p2 = t1 + w * (t4 + w * (t7 + w * (t10 + w * t13)))
                    p3 = t2 + w * (t5 + w * (t8 + w * (t11 + w * t14)))
                    p = z * p1 - (tt - w * (p2 + y * p3))
                    r += (tf + p)
                else:
                    y = x - one
                    p1 = y * (u0 + y * (u1 + y * (u2 + y * (u3 + y * (u4 + y * u5)))))
                    p2 = one + y * (v1 + y * (v2 + y * (v3 + y * (v4 + y * v5))))
                    r += (-0.5 * y + p1 / p2)

        # /* x < 8.0 */
        elif(ix < 4):
            i = int(x)
            t = zero
            y = x - i
            p = y * (s0 + y * (s1 + y * (s2 + y * (s3 + y * (s4 + y * (s5 + y * s6))))))
            q = one + y * (r1 + y * (r2 + y * (r3 + y * (r4 + y * (r5 + y * r6)))))
            r = half * y + p / q
            z = one  # /* lgamma(1+s) = log(s) + lgamma(s) */
            while (i > 2):
                i -= 1
                z *= (y + i)
                r += log(z)

        # /* 8.0 <= x < 2**58 */
        elif (ix < 58):
            t = log(x)
            z = one / x
            y = z * z
            w = w0 + z * (w1 + y * (w2 + y * (w3 + y * (w4 + y * (w5 + y * w6)))))
            r = (x - half) * (t - one) + w

        # /* 2**58 <= x <= inf */
        else:
            r = x * (log(x) - one)

        if (e < 0):
            r = nadj - r
        return signgamp * r


# Cephes Math Library Release 2.0:  April, 1987
# Copyright 1985, 1987 by Stephen L. Moshier
# Direct inquiries to 30 Frost Street, Cambridge, MA 02140
MACHEP = 0.0000001     # the machine roundoff error / tolerance
BIG = 4.503599627370496e15
BIGINV = 2.22044604925031308085e-16


def chisqprob(x, df):
    """Probability value (1-tail) for the Chi^2 probability distribution.

    Broadcasting rules apply.

    Parameters
    ----------
    x : array_like or float > 0

    df : array_like or float, probably int >= 1

    Returns
    -------
    chisqprob : ndarray
        The area from `chisq` to infinity under the Chi^2 probability
        distribution with degrees of freedom `df`.

    """
    if x <= 0:
        return 1.0
    if x == 0:
        return 0.0
    if df <= 0:
        raise ValueError("Domain error.")
    if x < 1.0 or x < df:
        return 1.0 - _igam(0.5 * df, 0.5 * x)
    return _igamc(0.5 * df, 0.5 * x)


def _igamc(a, x):
    """Complemented incomplete Gamma integral.

    SYNOPSIS:

    double a, x, y, igamc();

    y = igamc( a, x );

    DESCRIPTION:

    The function is defined by::

        igamc(a,x)   =   1 - igam(a,x)

                                inf.
                                   -
                          1       | |  -t  a-1
                    =   -----     |   e   t   dt.
                         -      | |
                        | (a)    -
                                    x

    In this implementation both arguments must be positive.
    The integral is evaluated by either a power series or
    continued fraction expansion, depending on the relative
    values of a and x.
    """
    # Compute  x**a * exp(-x) / Gamma(a)
    # TODO: Return this to math.lgamma once drop Python 2.6
    ax = math.exp(a * math.log(x) - x - _lgamma(a))

    # Continued fraction
    y = 1.0 - a
    z = x + y + 1.0
    c = 0.0
    pkm2 = 1.0
    qkm2 = x
    pkm1 = x + 1.0
    qkm1 = z * x
    ans = pkm1 / qkm1
    while True:
        c += 1.0
        y += 1.0
        z += 2.0
        yc = y * c
        pk = pkm1 * z - pkm2 * yc
        qk = qkm1 * z - qkm2 * yc
        if qk != 0:
            r = pk / qk
            t = abs((ans - r) / r)
            ans = r
        else:
            t = 1.0
        pkm2 = pkm1
        pkm1 = pk
        qkm2 = qkm1
        qkm1 = qk
        if abs(pk) > BIG:
                pkm2 *= BIGINV
                pkm1 *= BIGINV
                qkm2 *= BIGINV
                qkm1 *= BIGINV
        if t <= MACHEP:
            return ans * ax


def _igam(a, x):
    """Left tail of incomplete Gamma function.

    Computes this formula::

                 inf.      k
          a  -x   -       x
         x  e     >   ----------
                  -     -
                k=0   | (a+k+1)
    """

    # Compute  x**a * exp(-x) / Gamma(a)
    ax = math.exp(a * math.log(x) - x - math.lgamma(a))

    # Power series
    r = a
    c = 1.0
    ans = 1.0

    while True:
        r += 1.0
        c *= x / r
        ans += c
        if c / ans <= MACHEP:
            return ans * ax / a


# --- Speed ---

# try:
#    from scipy.stats import chisqprob
# except ImportError:
#    pass
