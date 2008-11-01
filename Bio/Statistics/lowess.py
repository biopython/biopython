"""
This module implements the Lowess function for nonparametric regression.

Functions:
lowess        Fit a smooth nonparametric regression curve to a scatterplot.

For more information, see

William S. Cleveland: "Robust locally weighted regression and smoothing
scatterplots", Journal of the American Statistical Association, December 1979,
volume 74, number 368, pp. 829-836.

William S. Cleveland and Susan J. Devlin: "Locally weighted regression: An
approach to regression analysis by local fitting", Journal of the American
Statistical Association, September 1988, volume 83, number 403, pp. 596-610.
"""

import numpy

try:
    from Bio.Cluster import median
    # The function median in Bio.Cluster is faster than the function median
    # in NumPy, as it does not require a full sort.
except ImportError, x:
    # Use the median function in NumPy if Bio.Cluster is not available
    from numpy import median

def lowess(x, y, f=2./3., iter=3):
  """lowess(x, y, f=2./3., iter=3) -> yest

Lowess smoother: Robust locally weighted regression.
The lowess function fits a nonparametric regression curve to a scatterplot.
The arrays x and y contain an equal number of elements; each pair
(x[i], y[i]) defines a data point in the scatterplot. The function returns
the estimated (smooth) values of y.

The smoothing span is given by f. A larger value for f will result in a
smoother curve. The number of robustifying iterations is given by iter. The
function will run faster with a smaller number of iterations."""
  n = len(x)
  r = int(numpy.ceil(f*n))
  h = [numpy.sort(abs(x-x[i]))[r] for i in range(n)]
  w = numpy.clip(abs(([x]-numpy.transpose([x]))/h),0.0,1.0)
  w = 1-w*w*w
  w = w*w*w
  yest = numpy.zeros(n)
  delta = numpy.ones(n)
  for iteration in range(iter):
    for i in range(n):
      weights = delta * w[:,i]
      b = numpy.array([sum(weights*y), sum(weights*y*x)])
      A = numpy.array([[sum(weights),   sum(weights*x)],
                       [sum(weights*x), sum(weights*x*x)]])
      beta = numpy.linalg.solve(A,b)
      yest[i] = beta[0] + beta[1]*x[i]
    residuals = y-yest
    s = numpy.median(abs(residuals))
    delta = numpy.clip(residuals/(6*s),-1,1)
    delta = 1-delta*delta
    delta = delta*delta
  return yest
