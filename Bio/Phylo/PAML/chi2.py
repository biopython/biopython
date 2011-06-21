# Copyright (C) 2011 by Brandon Invergo (b.invergo@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.
#
# This code is adapted (with permission) from the C source code of chi2.c, 
# written by Ziheng Yang and included in the PAML software package:
# http://abacus.gene.ucl.ac.uk/software/paml.html

from math import sqrt, log, exp

def cdf_chi2(df, stat):
    if df < 1:
        raise ValueError, "df must be at least 1"
    if stat < 0:
        raise ValueError, "The test statistic must be positive"
    x = 0.5 * stat
    alpha = df / 2.0
    prob = 1 - _incomplete_gamma(x, alpha)
    return prob
                      
def _ln_gamma_function(alpha):
    """Compute the log of the gamma function for a given alpha.
    
    Comments from Z. Yang:
    Returns ln(gamma(alpha)) for alpha>0, accurate to 10 decimal places.  
    Stirling's formula is used for the central polynomial part of the procedure.
    Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma function.
    Communications of the Association for Computing Machinery, 9:684
    """
    if alpha <= 0:
        raise ValueError
    x = alpha
    f = 0
    if x < 7:
        f = 1
        z = x
        while z<7:
            f *= z 
            z += 1
        x = z
        f = -log(f)
    z = 1 / (x * x)
    return  f + (x-0.5)*log(x) - x + .918938533204673             \
          + (((-.000595238095238*z+.000793650793651)*z-.002777777777778)*z \
               +.083333333333333)/x
        
def _incomplete_gamma(x, alpha):
    """Compute an incomplete gamma ratio.
    
    Comments from Z. Yang:
    Returns the incomplete gamma ratio I(x,alpha) where x is the upper 
           limit of the integration and alpha is the shape parameter.
    returns (-1) if in error
    ln_gamma_alpha = ln(Gamma(alpha)), is almost redundant.
    (1) series expansion     if alpha>x or x<=1
    (2) continued fraction   otherwise
    RATNEST FORTRAN by
    Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
    19: 285-287 (AS32)
    """
    p = alpha
    g = _ln_gamma_function(alpha)
    accurate = 1e-8
    overflow = 1e30
    gin = 0
    rn = 0
    a = 0
    b = 0
    an = 0
    dif = 0
    term = 0
    if x == 0:
       return 0
    if x < 0 or p <= 0:
        return -1
    factor = exp(p*log(x)-x-g)  
    if x > 1 and x >= p:
        a = 1 - p
        b = a + x + 1
        term = 0
        pn = [1, x, x+1, x*b, None, None]
        gin = pn[2] / pn[3]
    else:
        gin=1
        term=1
        rn=p
        while term > accurate:
            rn += 1
            term *= x / rn
            gin += term
        gin *= factor / p
        return gin
    while True:
        a += 1
        b += 2
        term += 1   
        an = a * term
        for i in range(2):
            pn[i + 4] = b * pn[i + 2] - an * pn[i]
        if pn[5] != 0:
            rn = pn[4] / pn[5]
            dif = abs(gin - rn)
            if dif > accurate:
                gin=rn
            elif dif <= accurate*rn:
                break
        for i in range(4):
            pn[i] = pn[i+2]
        if abs(pn[4]) < overflow:
            continue
        for i in range(4):
            pn[i] /= overflow
    gin = 1 - factor * gin
    return gin
