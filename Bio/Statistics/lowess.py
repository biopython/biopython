"""
This module implements the Lowess function for nonparametric regression.

Functions:
lowess        Fit a smooth nonparametric regression curve to a scatterplot.

For more information, see
William S. Cleveland: "Robust locally weighted regression and smoothing
scatterplots", Journal of the American Statistical Association, December 1979,
volume 74, number 368, pp. 829-836.
"""

try:
    from Numeric import *
    from LinearAlgebra import solve_linear_equations
except ImportError, x:
    raise ImportError, "This module requires Numeric (precursor to NumPy) with the LinearAlgebra and MLab libraries"

try:
    from Bio.Cluster import median
    # The function median in Bio.Cluster is faster than the function median
    # in Numeric's MLab, as it does not require a full sort.
except ImportError, x:
    # Use the median function in Numeric's MLab if Bio.Cluster is not available
    try:
        from MLab import median
    except ImportError, x:
        raise ImportError, "This module requires Numeric (precursor to NumPy) with the LinearAlgebra and MLab libraries"

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
  r = int(ceil(f*n))
  h = [sort(abs(x-x[i]))[r] for i in range(n)]
  w = clip(abs(([x]-transpose([x]))/h),0.0,1.0)
  w = 1-w*w*w
  w = w*w*w
  yest = zeros(n,'d')
  delta = ones(n,'d')
  for iteration in range(iter):
    for i in range(n):
      weights = delta * w[:,i]
      b = array([sum(weights*y), sum(weights*y*x)])
      A = array([[sum(weights),   sum(weights*x)],
                 [sum(weights*x), sum(weights*x*x)]])
      beta = solve_linear_equations(A,b)
      yest[i] = beta[0] + beta[1]*x[i]
    residuals = y-yest
    s = median(abs(residuals))
    delta = clip(residuals/(6*s),-1,1)
    delta = 1-delta*delta
    delta = delta*delta
  return yest
