/* The C clustering library for cDNA microarray data.
 * Copyright (C) 2002 Michiel Jan Laurens de Hoon.
 *
 * This library was written at the Laboratory of DNA Information Analysis,
 * Human Genome Center, Institute of Medical Science, University of Tokyo,
 * 4-6-1 Shirokanedai, Minato-ku, Tokyo 108-8639, Japan.
 * Contact: mdehoon@ims.u-tokyo.ac.jp
 * 
 * Permission to use, copy, modify, and distribute this software and its
 * documentation with or without modifications and for any purpose and
 * without fee is hereby granted, provided that any copyright notices
 * appear in all copies and that both those copyright notices and this
 * permission notice appear in supporting documentation, and that the
 * names of the contributors or copyright holders not be used in
 * advertising or publicity pertaining to distribution of the software
 * without specific prior permission.
 * 
 * THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
 * WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
 * CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
 * OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
 * OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
 * OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
 * OR PERFORMANCE OF THIS SOFTWARE.
 * 
 */

#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "ranlib.h"
#include "cluster.h"
#ifdef WINDOWS
#  include <windows.h>
#endif

/* ************************************************************************ */

#ifdef WINDOWS
/* Then we make a Windows DLL */
int WINAPI
clusterdll_init (HANDLE h, DWORD reason, void* foo)
{
  return 1;
}
#endif

/* ************************************************************************ */

double CALL mean(int n, double x[])
{ double result = 0.;
  int i;
  for (i = 0; i < n; i++) result += x[i];
  result /= n;
  return result;
}

/* ************************************************************************ */

double CALL median (int n, double x[])
/*
Find the median of X(1), ... , X(N), using as much of the quicksort
algorithm as is needed to isolate it.
N.B. On exit, the array X is partially ordered.
Based on Alan J. Miller's median.f90 routine.
*/

{ int i, j;
  int nr = n / 2;
  int nl = nr - 1;
  int even = 0;
  /* hi & lo are position limits encompassing the median. */
  int lo = 0;
  int hi = n-1;

  if (n==2*nr) even = 1;
  if (n<3)
  { if (n<1) return 0.;
    if (n == 1) return x[0];
    return 0.5*(x[0]+x[1]);
  }

  /* Find median of 1st, middle & last values. */
  do
  { int loop;
    int mid = (lo + hi)/2;
    double result = x[mid];
    double xlo = x[lo];
    double xhi = x[hi];
    if (xhi<xlo)
    { double temp = xlo;
      xlo = xhi;
      xhi = temp;
    }
    if (result>xhi) result = xhi;
    else if (result<xlo) result = xlo;
    /* The basic quicksort algorithm to move all values <= the sort key (XMED)
     * to the left-hand end, and all higher values to the other end.
     */
    i = lo;
    j = hi;
    do
    { while (x[i]<result) i++;
      while (x[j]>result) j--;
      loop = 0;
      if (i<j)
      { double temp = x[i];
        x[i] = x[j];
        x[j] = temp;
        i++;
        j--;
        if (i<=j) loop = 1;
      }
    } while (loop); /* Decide which half the median is in. */

    if (even)
    { if (j==nl && i==nr)
        /* Special case, n even, j = n/2 & i = j + 1, so the median is
         * between the two halves of the series.   Find max. of the first
         * half & min. of the second half, then average.
         */
        { int k;
          double xmax = x[0];
          double xmin = x[n-1];
          for (k = lo; k <= j; k++) xmax = max(xmax,x[k]);
          for (k = i; k <= hi; k++) xmin = min(xmin,x[k]);
          return 0.5*(xmin + xmax);
        }
      if (j<nl) lo = i;
      if (i>nr) hi = j;
      if (i==j)
      { if (i==nl) lo = nl;
        if (j==nr) hi = nr;
      }
    }
    else
    { if (j<nr) lo = i;
      if (i>nr) hi = j;
      /* Test whether median has been isolated. */
      if (i==j && i==nr) return result;
    }
  }
  while (lo<hi-1);

  if (even) return (0.5*(x[nl]+x[nr]));
  if (x[lo]>x[hi])
  { double temp = x[lo];
    x[lo] = x[hi];
    x[hi] = temp;
  }
  return x[nr];
}

/* *********************************************************************  */

static
int compare(const void* a, const void* b)
/* Helper function for sort. Previously, this was a nested function under sort,
 * which is not allowed under ANSI C.
 */
{ const double term1 = *(*(double**)a);
  const double term2 = *(*(double**)b);
  if (term1 < term2) return -1;
  if (term1 > term2) return +1;
  return 0;
}

void CALL sort(int n, const double data[], int index[])
/* Sets up an index table given the data, such that data[index[]] is in
 * increasing order. Sorting is done on the pointers, from which the indeces
 * are recalculated. The array data is unchanged.
 */
{ int i;
  const double** p = malloc(n*sizeof(double*));
  const double* start = data;
  for (i = 0; i < n; i++) p[i] = &(data[i]);
  qsort(p, n, sizeof(double*), compare);
  for (i = 0; i < n; i++) index[i] = (int)(p[i]-start);
  free(p);
}

static
void getrank (int n, double data[], double rank[])
/* Calculates the ranks of the elements in the array data. Two elements with the
 * same value get the same rank, equal to the average of the ranks had the
 * elements different values.
 */
{ int i;
  int* index = malloc(n*sizeof(int));
  /* Call sort to get an index table */
  sort (n, data, index);
  /* Build a rank table */
  for (i = 0; i < n; i++) rank[index[i]] = i;
  /* Fix for equal ranks */
  i = 0;
  while (i < n)
  { int m;
    double value = data[index[i]];
    int j = i + 1;
    while (j < n && data[index[j]] == value) j++;
    m = j - i; /* number of equal ranks found */
    value = rank[index[i]] + (m-1)/2.;
    for (j = i; j < i + m; j++) rank[index[j]] = value;
    i += m;
  }
  free (index);
  return;
}

/* ********************************************************************* */

void CALL svd(int m, int n, double** u, double w[], double** v, int* ierr)
/*
 *   This subroutine is a translation of the Algol procedure svd,
 *   Num. Math. 14, 403-420(1970) by Golub and Reinsch.
 *   Handbook for Auto. Comp., Vol II-Linear Algebra, 134-151(1971).
 *
 *   This subroutine determines the singular value decomposition
 *        t
 *   A=usv  of a real m by n rectangular matrix, where m is greater
 *   then or equal to n.  Householder bidiagonalization and a variant
 *   of the QR algorithm are used.
 *  
 *
 *   On input.
 *
 *      m is the number of rows of A (and u).
 *
 *      n is the number of columns of A (and u) and the order of v.
 *
 *      u contains the rectangular input matrix A to be decomposed.
 *
 *   On output.
 *
 *      w contains the n (non-negative) singular values of a (the
 *        diagonal elements of s).  they are unordered.  if an
 *        error exit is made, the singular values should be correct
 *        for indices ierr+1,ierr+2,...,n.
 *
 *      u contains the matrix u (orthogonal column vectors) of the
 *        decomposition.
 *        if an error exit is made, the columns of u corresponding
 *        to indices of correct singular values should be correct.
 *
 *      v contains the matrix v (orthogonal) of the decomposition.
 *        if an error exit is made, the columns of v corresponding
 *        to indices of correct singular values should be correct.
 *
 *      ierr is set to
 *        zero       for normal return,
 *        k          if the k-th singular value has not been
 *                   determined after 30 iterations.
 *
 *   Questions and comments should be directed to B. S. Garbow,
 *   Applied Mathematics division, Argonne National Laboratory
 *
 *   Modified to eliminate machep
 *
 *   Translated to C by Michiel de Hoon, Human Genome Center,
 *   University of Tokyo, for inclusion in the C Clustering Library.
 *   This routine is less general than the original svd routine, as
 *   it focuses on the singular value decomposition as needed for
 *   clustering. In particular,
 *     - We require m >= n
 *     - We calculate both u and v in all cases
 *     - We pass the input array A via u; this array is subsequently
 *       overwritten.
 *     - We allocate for the array rv1, used as a working space,
 *       internally in this routine, instead of passing it as an
 *       argument.
 *   2003.06.05
 */
{ int i, j, k, i1, k1, l1, its;
  double* rv1 = malloc(n*sizeof(double));
  double c,f,h,s,x,y,z;
  int l = 0;
  double g = 0.0;
  double scale = 0.0;
  double anorm = 0.0;
  *ierr = 0;
  /* Householder reduction to bidiagonal form */
  for (i = 0; i < n; i++)
  { l = i + 1;
    rv1[i] = scale * g;
    g = 0.0;
    s = 0.0;
    scale = 0.0;
    for (k = i; k < m; k++) scale += fabs(u[k][i]);
    if (scale != 0.0)
    { for (k = i; k < m; k++)
      { u[k][i] /= scale;
        s += u[k][i]*u[k][i];
      }
      f = u[i][i];
      g = (f >= 0) ? -sqrt(s) : sqrt(s);
      h = f * g - s;
      u[i][i] = f - g;
      if (i < n-1)
      { for (j = l; j < n; j++)
        { s = 0.0;
          for (k = i; k < m; k++) s += u[k][i] * u[k][j];
          f = s / h;
          for (k = i; k < m; k++) u[k][j] += f * u[k][i];
        }
      }
      for (k = i; k < m; k++) u[k][i] *= scale;
    }
    w[i] = scale * g;
    g = 0.0;
    s = 0.0;
    scale = 0.0;
    if (i<n-1)
    { for (k = l; k < n; k++) scale += fabs(u[i][k]);
      if (scale != 0.0)
      { for (k = l; k < n; k++)
        { u[i][k] /= scale;
          s += u[i][k] * u[i][k];
        }
        f = u[i][l];
        g = (f >= 0) ? -sqrt(s) : sqrt(s);
        h = f * g - s;
        u[i][l] = f - g;
        for (k = l; k < n; k++) rv1[k] = u[i][k] / h;
        for (j = l; j < m; j++)
        { s = 0.0;
          for (k = l; k < n; k++) s += u[j][k] * u[i][k];
          for (k = l; k < n; k++) u[j][k] += s * rv1[k];
        }
        for (k = l; k < n; k++)  u[i][k] *= scale;
      }
    }
    anorm = max(anorm,fabs(w[i])+fabs(rv1[i]));
  }
  /* accumulation of right-hand transformations */
  for (i = n-1; i>=0; i--)
  { if (i < n-1)
    { if (g != 0.0)
      { for (j = l; j < n; j++) v[j][i] = (u[i][j] / u[i][l]) / g;
        /* double division avoids possible underflow */
        for (j = l; j < n; j++)
        { s = 0.0;
          for (k = l; k < n; k++) s += u[i][k] * v[k][j];
          for (k = l; k < n; k++) v[k][j] += s * v[k][i];
        }
      }
    }
    for (j = l; j < n; j++)
    { v[i][j] = 0.0;
      v[j][i] = 0.0;
    }
    v[i][i] = 1.0;
    g = rv1[i];
    l = i;
  }
  /* accumulation of left-hand transformations */
  for (i = n-1; i >= 0; i--)
  { l = i + 1;
    g = w[i];
    if (i!=n-1)
      for (j = l; j < n; j++) u[i][j] = 0.0;
    if (g!=0.0)
    { if (i!=n-1)
      { for (j = l; j < n; j++)
        { s = 0.0;
          for (k = l; k < m; k++) s += u[k][i] * u[k][j];
          /* double division avoids possible underflow */
          f = (s / u[i][i]) / g;
          for (k = i; k < m; k++) u[k][j] += f * u[k][i];
        }
      }
      for (j = i; j < m; j++) u[j][i] /= g;
    }
    else
      for (j = i; j < m; j++) u[j][i] = 0.0;
    u[i][i] += 1.0;
  }
  /* diagonalization of the bidiagonal form */
  for (k = n-1; k >= 0; k--)
  { k1 = k-1;
    its = 0;
    while(1)
    /* test for splitting */
    { for (l = k; l >= 0; l--)
      { l1 = l-1;
        if (fabs(rv1[l]) + anorm == anorm) break;
        /* rv1[0] is always zero, so there is no exit
         * through the bottom of the loop */
        if (fabs(w[l1]) + anorm == anorm)
        /* cancellation of rv1[l] if l greater than 0 */
        { c = 0.0;
          s = 1.0;
          for (i = l; i <= k; i++)
          { f = s * rv1[i];
            rv1[i] *= c;
            if (fabs(f) + anorm == anorm) break;
            g = w[i];
            h = sqrt(f*f+g*g);
            w[i] = h;
            c = g / h;
            s = -f / h;
            for (j = 0; j < m; j++)
            { y = u[j][l1];
              z = u[j][i];
              u[j][l1] = y * c + z * s;
              u[j][i] = -y * s + z * c;
            }
          }
          break;
        }
      }
      /* test for convergence */
      z = w[k];
      if (l==k) /* convergence */
      { if (z < 0.0)
        /*  w[k] is made non-negative */
        { w[k] = -z;
          for (j = 0; j < n; j++) v[j][k] = -v[j][k];
        }
        break;
      }
      else if (its==30)
      { *ierr = k;
        break;
      }
      else
      /* shift from bottom 2 by 2 minor */
      { its++;
        x = w[l];
        y = w[k1];
        g = rv1[k1];
        h = rv1[k];
        f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
        g = sqrt(f*f+1.0);
        f = ((x - z) * (x + z) + h * (y / (f + (f >= 0 ? g : -g)) - h)) / x;
        /* next qr transformation */
        c = 1.0;
        s = 1.0;
        for (i1 = l; i1 <= k1; i1++)
        { i = i1 + 1;
          g = rv1[i];
          y = w[i];
          h = s * g;
          g = c * g;
          z = sqrt(f*f+h*h);
          rv1[i1] = z;
          c = f / z;
          s = h / z;
          f = x * c + g * s;
          g = -x * s + g * c;
          h = y * s;
          y = y * c;
          for (j = 0; j < n; j++)
          { x = v[j][i1];
            z = v[j][i];
            v[j][i1] = x * c + z * s;
            v[j][i] = -x * s + z * c;
          }
          z = sqrt(f*f+h*h);
          w[i1] = z;
          /* rotation can be arbitrary if z is zero */
          if (z!=0.0)
          { c = f / z;
            s = h / z;
          }
          f = c * g + s * y;
          x = -s * g + c * y;
          for (j = 0; j < m; j++)
          { y = u[j][i1];
            z = u[j][i];
            u[j][i1] = y * c + z * s;
            u[j][i] = -y * s + z * c;
          }
        }
        rv1[l] = 0.0;
        rv1[k] = f;
        w[k] = x;
      }
    }
  }
  free(rv1);
  return;
}

/* ********************************************************************* */

static
double euclid (int n, double** data1, double** data2, int** mask1, int** mask2,
  const double weight[], int index1, int index2, int transpose)
 
/*
Purpose
=======

The euclid routine calculates the weighted Euclidean distance between two
rows or columns in a matrix.

Arguments
=========

n      (input) int
The number of elements in a row or column. If transpose==0, then n is the number
of columns; otherwise, n is the number of rows.

data1  (input) double array
The data array containing the first vector.

data2  (input) double array
The data array containing the second vector.

mask1  (input) int array
This array which elements in data1 are missing. If mask1[i][j]==0, then
data1[i][j] is missing.

mask2  (input) int array
This array which elements in data2 are missing. If mask2[i][j]==0, then
data2[i][j] is missing.

weight (input) double array, dimension( n )
The weights that are used to calculate the distance.

index1     (input) int
Index of the first row or column.

index2     (input) int
Index of the second row or column.

transpose (input) int
If transpose==0, the distance between two rows in the matrix is calculated.
Otherwise, the distance between two columns in the matrix is calculated.

============================================================================
*/
{ double result = 0.;
  double tweight = 0;
  int i;
  if (transpose==0) /* Calculate the distance between two rows */
  { for (i = 0; i < n; i++)
    { if (mask1[index1][i] && mask2[index2][i])
      { double term = data1[index1][i] - data2[index2][i];
        result = result + weight[i]*term*term;
        tweight += weight[i];
      }
    }
  }
  else
  { for (i = 0; i < n; i++)
    { if (mask1[i][index1] && mask2[i][index2])
      { double term = data1[i][index1] - data2[i][index2];
        result = result + weight[i]*term*term;
        tweight += weight[i];
      }
    }
  }
  if (!tweight) return 0; /* usually due to empty clusters */
  result /= tweight;
  result *= n;
  return result;
}

/* ********************************************************************* */

static
double harmonic(int n, double** data1, double** data2, int** mask1, int** mask2,
  const double weight[], int index1, int index2, int transpose)
 
/*
Purpose
=======

The harmonic routine calculates the weighted Euclidean distance between two
rows or columns in a matrix, adding terms for the different dimensions
harmonically, i.e. summing the inverse and taking the inverse of the total.

Arguments
=========

n      (input) int
The number of elements in a row or column. If transpose==0, then n is the number
of columns; otherwise, n is the number of rows.

data1  (input) double array
The data array containing the first vector.

data2  (input) double array
The data array containing the second vector.

mask1  (input) int array
This array which elements in data1 are missing. If mask1[i][j]==0, then
data1[i][j] is missing.

mask2  (input) int array
This array which elements in data2 are missing. If mask2[i][j]==0, then
data2[i][j] is missing.

weight (input) double array, dimension( n )
The weights that are used to calculate the distance.

index1     (input) int
Index of the first row or column.

index2     (input) int
Index of the second row or column.

transpose (input) int
If transpose==0, the distance between two rows in the matrix is calculated.
Otherwise, the distance between two columns in the matrix is calculated.

============================================================================
*/
{ double result = 0.;
  double tweight = 0;
  int i;
  if (transpose==0) /* Calculate the distance between two rows */
  { for (i = 0; i < n; i++)
    { if (mask1[index1][i] && mask2[index2][i])
      { const double term = data1[index1][i] - data2[index2][i];
        if (term==0) return 0;
        result = result + weight[i]/(term*term);
        tweight += weight[i];
      }
    }
  }
  else
  { for (i = 0; i < n; i++)
    { if (mask1[i][index1] && mask2[i][index2])
      { const double term = data1[i][index1] - data2[i][index2];
        if (term==0) return 0;
        result = result + weight[i]/(term*term);
        tweight += weight[i];
      }
    }
  }
  if (!tweight) return 0; /* usually due to empty clusters */
  result /= tweight;
  result *= n;
  result = 1. / result;
  return result;
}

/* ********************************************************************* */

static
double cityblock (int n, double** data1, double** data2, int** mask1,
  int** mask2, const double weight[], int index1, int index2, int transpose)

/*
Purpose
=======

The cityblock routine calculates the weighted "City Block" distance between
two rows or columns in a matrix. City Block distance is defined as the
absolute value of X1-X2 plus the absolute value of Y1-Y2 plus..., which is
equivalent to taking an "up and over" path.

Arguments
=========

n      (input) int
The number of elements in a row or column. If transpose==0, then n is the number
of columns; otherwise, n is the number of rows.

data1  (input) double array
The data array containing the first vector.

data2  (input) double array
The data array containing the second vector.

mask1  (input) int array
This array which elements in data1 are missing. If mask1[i][j]==0, then
data1[i][j] is missing.

mask2  (input) int array
This array which elements in data2 are missing. If mask2[i][j]==0, then
data2[i][j] is missing.

weight (input) double array, dimension( n )
The weights that are used to calculate the distance.

index1     (input) int
Index of the first row or column.

index2     (input) int
Index of the second row or column.

transpose (input) int
If transpose==0, the distance between two rows in the matrix is calculated.
Otherwise, the distance between two columns in the matrix is calculated.

============================================================================ */
{ double result = 0.;
  double tweight = 0;
  int i;
  if (transpose==0) /* Calculate the distance between two rows */
  { for (i = 0; i < n; i++)
    { if (mask1[index1][i] && mask2[index2][i])
      { double term = data1[index1][i] - data2[index2][i];
        result = result + weight[i]*fabs(term);
        tweight += weight[i];
      }
    }
  }
  else
  { for (i = 0; i < n; i++)
    { if (mask1[i][index1] && mask2[i][index2])
      { double term = data1[i][index1] - data2[i][index2];
        result = result + weight[i]*fabs(term);
        tweight += weight[i];
      }
    }
  }
  if (!tweight) return 0; /* usually due to empty clusters */
  result /= tweight;
  result *= n;
  return result;
}

/* ********************************************************************* */

static
double correlation (int n, double** data1, double** data2, int** mask1,
  int** mask2, const double weight[], int index1, int index2, int transpose)
/*
Purpose
=======

The correlation routine calculates the weighted Pearson distance between two
rows or columns in a matrix. We define the Pearson distance as one minus the
Pearson correlation.
This definition yields a semi-metric: d(a,b) >= 0, and d(a,b) = 0 iff a = b.
but the triangular inequality d(a,b) + d(b,c) >= d(a,c) does not hold
(e.g., choose b = a + c).

Arguments
=========

n      (input) int
The number of elements in a row or column. If transpose==0, then n is the number
of columns; otherwise, n is the number of rows.

data1  (input) double array
The data array containing the first vector.

data2  (input) double array
The data array containing the second vector.

mask1  (input) int array
This array which elements in data1 are missing. If mask1[i][j]==0, then
data1[i][j] is missing.

mask2  (input) int array
This array which elements in data2 are missing. If mask2[i][j]==0, then
data2[i][j] is missing.

weight (input) double array, dimension( n )
The weights that are used to calculate the distance.

index1     (input) int
Index of the first row or column.

index2     (input) int
Index of the second row or column.

transpose (input) int
If transpose==0, the distance between two rows in the matrix is calculated.
Otherwise, the distance between two columns in the matrix is calculated.
============================================================================
*/
{ double result = 0.;
  double sum1 = 0.;
  double sum2 = 0.;
  double denom1 = 0.;
  double denom2 = 0.;
  double tweight = 0.;
  if (transpose==0) /* Calculate the distance between two rows */
  { int i;
    for (i = 0; i < n; i++)
    { if (mask1[index1][i] && mask2[index2][i])
      { double term1 = data1[index1][i];
        double term2 = data2[index2][i];
        double w = weight[i];
        sum1 += w*term1;
        sum2 += w*term2;
        result += w*term1*term2;
        denom1 += w*term1*term1;
        denom2 += w*term2*term2;
        tweight += w;
      }
    }
  }
  else
  { int i;
    for (i = 0; i < n; i++)
    { if (mask1[i][index1] && mask2[i][index2])
      { double term1 = data1[i][index1];
        double term2 = data2[i][index2];
        double w = weight[i];
        sum1 += w*term1;
        sum2 += w*term2;
        result += w*term1*term2;
        denom1 += w*term1*term1;
        denom2 += w*term2*term2;
        tweight += w;
      }
    }
  }
  if (!tweight) return 0; /* usually due to empty clusters */
  result -= sum1 * sum2 / tweight;
  denom1 -= sum1 * sum1 / tweight;
  denom2 -= sum2 * sum2 / tweight;
  if (denom1 <= 0) return 1; /* include '<' to deal with roundoff errors */
  if (denom2 <= 0) return 1; /* include '<' to deal with roundoff errors */
  result = result / sqrt(denom1*denom2);
  result = 1. - result;
  return result;
}

/* ********************************************************************* */

static
double acorrelation (int n, double** data1, double** data2, int** mask1,
  int** mask2, const double weight[], int index1, int index2, int transpose)
/*
Purpose
=======

The acorrelation routine calculates the weighted Pearson distance between two
rows or columns, using the absolute value of the correlation.
This definition yields a semi-metric: d(a,b) >= 0, and d(a,b) = 0 iff a = b.
but the triangular inequality d(a,b) + d(b,c) >= d(a,c) does not hold
(e.g., choose b = a + c).

Arguments
=========

n      (input) int
The number of elements in a row or column. If transpose==0, then n is the number
of columns; otherwise, n is the number of rows.

data1  (input) double array
The data array containing the first vector.

data2  (input) double array
The data array containing the second vector.

mask1  (input) int array
This array which elements in data1 are missing. If mask1[i][j]==0, then
data1[i][j] is missing.

mask2  (input) int array
This array which elements in data2 are missing. If mask2[i][j]==0, then
data2[i][j] is missing.

weight (input) double array, dimension( n )
The weights that are used to calculate the distance.

index1     (input) int
Index of the first row or column.

index2     (input) int
Index of the second row or column.

transpose (input) int
If transpose==0, the distance between two rows in the matrix is calculated.
Otherwise, the distance between two columns in the matrix is calculated.
============================================================================
*/
{ double result = 0.;
  double sum1 = 0.;
  double sum2 = 0.;
  double denom1 = 0.;
  double denom2 = 0.;
  double tweight = 0.;
  if (transpose==0) /* Calculate the distance between two rows */
  { int i;
    for (i = 0; i < n; i++)
    { if (mask1[index1][i] && mask2[index2][i])
      { double term1 = data1[index1][i];
        double term2 = data2[index2][i];
        double w = weight[i];
        sum1 += w*term1;
        sum2 += w*term2;
        result += w*term1*term2;
        denom1 += w*term1*term1;
        denom2 += w*term2*term2;
        tweight += w;
      }
    }
  }
  else
  { int i;
    for (i = 0; i < n; i++)
    { if (mask1[i][index1] && mask2[i][index2])
      { double term1 = data1[i][index1];
        double term2 = data2[i][index2];
        double w = weight[i];
        sum1 += w*term1;
        sum2 += w*term2;
        result += w*term1*term2;
        denom1 += w*term1*term1;
        denom2 += w*term2*term2;
        tweight += w;
      }
    }
  }
  if (!tweight) return 0; /* usually due to empty clusters */
  result -= sum1 * sum2 / tweight;
  denom1 -= sum1 * sum1 / tweight;
  denom2 -= sum2 * sum2 / tweight;
  if (denom1 <= 0) return 1; /* include '<' to deal with roundoff errors */
  if (denom2 <= 0) return 1; /* include '<' to deal with roundoff errors */
  result = fabs(result) / sqrt(denom1*denom2);
  result = 1. - result;
  return result;
}

/* ********************************************************************* */

static
double ucorrelation (int n, double** data1, double** data2, int** mask1,
  int** mask2, const double weight[], int index1, int index2, int transpose)
/*
Purpose
=======

The ucorrelation routine calculates the weighted Pearson distance between two
rows or columns, using the uncentered version of the Pearson correlation. In the
uncentered Pearson correlation, a zero mean is used for both vectors even if
the actual mean is nonzero.
This definition yields a semi-metric: d(a,b) >= 0, and d(a,b) = 0 iff a = b.
but the triangular inequality d(a,b) + d(b,c) >= d(a,c) does not hold
(e.g., choose b = a + c).

Arguments
=========

n      (input) int
The number of elements in a row or column. If transpose==0, then n is the number
of columns; otherwise, n is the number of rows.

data1  (input) double array
The data array containing the first vector.

data2  (input) double array
The data array containing the second vector.

mask1  (input) int array
This array which elements in data1 are missing. If mask1[i][j]==0, then
data1[i][j] is missing.

mask2  (input) int array
This array which elements in data2 are missing. If mask2[i][j]==0, then
data2[i][j] is missing.

weight (input) double array, dimension( n )
The weights that are used to calculate the distance.

index1     (input) int
Index of the first row or column.

index2     (input) int
Index of the second row or column.

transpose (input) int
If transpose==0, the distance between two rows in the matrix is calculated.
Otherwise, the distance between two columns in the matrix is calculated.
============================================================================
*/
{ double result = 0.;
  double denom1 = 0.;
  double denom2 = 0.;
  int flag = 0;
  /* flag will remain zero if no nonzero combinations of mask1 and mask2 are
   * found.
   */
  if (transpose==0) /* Calculate the distance between two rows */
  { int i;
    for (i = 0; i < n; i++)
    { if (mask1[index1][i] && mask2[index2][i])
      { double term1 = data1[index1][i];
        double term2 = data2[index2][i];
        double w = weight[i];
        result += w*term1*term2;
        denom1 += w*term1*term1;
        denom2 += w*term2*term2;
        flag = 1;
      }
    }
  }
  else
  { int i;
    for (i = 0; i < n; i++)
    { if (mask1[i][index1] && mask2[i][index2])
      { double term1 = data1[i][index1];
        double term2 = data2[i][index2];
        double w = weight[i];
        result += w*term1*term2;
        denom1 += w*term1*term1;
        denom2 += w*term2*term2;
        flag = 1;
      }
    }
  }
  if (!flag) return 0.;
  if (denom1==0.) return 1.;
  if (denom2==0.) return 1.;
  result = result / sqrt(denom1*denom2);
  result = 1. - result;
  return result;
}

/* ********************************************************************* */

static
double uacorrelation (int n, double** data1, double** data2, int** mask1,
  int** mask2, const double weight[], int index1, int index2, int transpose)
/*
Purpose
=======

The uacorrelation routine calculates the weighted Pearson distance between two
rows or columns, using the absolute value of the uncentered version of the
Pearson correlation. In the uncentered Pearson correlation, a zero mean is used
for both vectors even if the actual mean is nonzero.
This definition yields a semi-metric: d(a,b) >= 0, and d(a,b) = 0 iff a = b.
but the triangular inequality d(a,b) + d(b,c) >= d(a,c) does not hold
(e.g., choose b = a + c).

Arguments
=========

n      (input) int
The number of elements in a row or column. If transpose==0, then n is the number
of columns; otherwise, n is the number of rows.

data1  (input) double array
The data array containing the first vector.

data2  (input) double array
The data array containing the second vector.

mask1  (input) int array
This array which elements in data1 are missing. If mask1[i][j]==0, then
data1[i][j] is missing.

mask2  (input) int array
This array which elements in data2 are missing. If mask2[i][j]==0, then
data2[i][j] is missing.

weight (input) double array, dimension( n )
The weights that are used to calculate the distance.

index1     (input) int
Index of the first row or column.

index2     (input) int
Index of the second row or column.

transpose (input) int
If transpose==0, the distance between two rows in the matrix is calculated.
Otherwise, the distance between two columns in the matrix is calculated.
============================================================================
*/
{ double result = 0.;
  double denom1 = 0.;
  double denom2 = 0.;
  int flag = 0;
  /* flag will remain zero if no nonzero combinations of mask1 and mask2 are
   * found.
   */
  if (transpose==0) /* Calculate the distance between two rows */
  { int i;
    for (i = 0; i < n; i++)
    { if (mask1[index1][i] && mask2[index2][i])
      { double term1 = data1[index1][i];
        double term2 = data2[index2][i];
        double w = weight[i];
        result += w*term1*term2;
        denom1 += w*term1*term1;
        denom2 += w*term2*term2;
        flag = 1;
      }
    }
  }
  else
  { int i;
    for (i = 0; i < n; i++)
    { if (mask1[i][index1] && mask2[i][index2])
      { double term1 = data1[i][index1];
        double term2 = data2[i][index2];
        double w = weight[i];
        result += w*term1*term2;
        denom1 += w*term1*term1;
        denom2 += w*term2*term2;
        flag = 1;
      }
    }
  }
  if (!flag) return 0.;
  if (denom1==0.) return 1.;
  if (denom2==0.) return 1.;
  result = fabs(result) / sqrt(denom1*denom2);
  result = 1. - result;
  return result;
}

/* *********************************************************************  */

static
double spearman (int n, double** data1, double** data2, int** mask1,
  int** mask2, const double weight[], int index1, int index2, int transpose)
/*
Purpose
=======

The spearman routine calculates the Spearman distance between two rows or
columns. The Spearman distance is defined as one minus the Spearman rank
correlation.

Arguments
=========

n      (input) int
The number of elements in a row or column. If transpose==0, then n is the number
of columns; otherwise, n is the number of rows.

data1  (input) double array
The data array containing the first vector.

data2  (input) double array
The data array containing the second vector.

mask1  (input) int array
This array which elements in data1 are missing. If mask1[i][j]==0, then
data1[i][j] is missing.

mask2  (input) int array
This array which elements in data2 are missing. If mask2[i][j]==0, then
data2[i][j] is missing.

weight (input) double array, dimension( n )
These weights are ignored, but included for consistency with other distance
measures.

index1     (input) int
Index of the first row or column.

index2     (input) int
Index of the second row or column.

transpose (input) int
If transpose==0, the distance between two rows in the matrix is calculated.
Otherwise, the distance between two columns in the matrix is calculated.
============================================================================
*/
{ double* rank1;
  double* rank2;
  double result = 0.;
  double denom1 = 0.;
  double denom2 = 0.;
  double avgrank;
  double* tdata1 = malloc(n*sizeof(double));
  double* tdata2 = malloc(n*sizeof(double));
  int i;
  int m = 0;
  if (transpose==0)
  { for (i = 0; i < n; i++)
    { if (mask1[index1][i] && mask2[index2][i])
      { tdata1[m] = data1[index1][i];
        tdata2[m] = data2[index2][i];
        m++;
      }
    }
  }
  else
  { for (i = 0; i < n; i++)
    { if (mask1[i][index1] && mask2[i][index2])
      { tdata1[m] = data1[i][index1];
        tdata2[m] = data2[i][index2];
        m++;
      }
    }
  }
  if (m==0) return 0;
  rank1 = malloc(m*sizeof(double));
  rank2 = malloc(m*sizeof(double));
  getrank(m, tdata1, rank1);
  free(tdata1);
  getrank(m, tdata2, rank2);
  free(tdata2);
  avgrank = 0.5*(m-1); /* Average rank */
  for (i = 0; i < m; i++)
  { const double value1 = rank1[i];
    const double value2 = rank2[i];
    result += value1 * value2;
    denom1 += value1 * value1;
    denom2 += value2 * value2;
  }
  /* Note: denom1 and denom2 cannot be calculated directly from the number
   * of elements. If two elements have the same rank, the squared sum of
   * their ranks will change.
   */
  free(rank1);
  free(rank2);
  result /= m;
  denom1 /= m;
  denom2 /= m;
  result -= avgrank * avgrank;
  denom1 -= avgrank * avgrank;
  denom2 -= avgrank * avgrank;
  result = result / sqrt(denom1*denom2);
  result = 1. - result;
  return result;
}

/* *********************************************************************  */

static
double kendall (int n, double** data1, double** data2, int** mask1, int** mask2,
  const double weight[], int index1, int index2, int transpose)
/*
Purpose
=======

The kendall routine calculates the Kendall distance between two
rows or columns. The Kendall distance is defined as one minus Kendall's tau.

Arguments
=========

n      (input) int
The number of elements in a row or column. If transpose==0, then n is the number
of columns; otherwise, n is the number of rows.

data1  (input) double array
The data array containing the first vector.

data2  (input) double array
The data array containing the second vector.

mask1  (input) int array
This array which elements in data1 are missing. If mask1[i][j]==0, then
data1[i][j] is missing.

mask2  (input) int array
This array which elements in data2 are missing. If mask2[i][j]==0, then
data2[i][j] is missing.

weight (input) double array, dimension( n )
These weights are ignored, but included for consistency with other distance
measures.

index1     (input) int
Index of the first row or column.

index2     (input) int
Index of the second row or column.

transpose (input) int
If transpose==0, the distance between two rows in the matrix is calculated.
Otherwise, the distance between two columns in the matrix is calculated.
============================================================================
*/
{ int con = 0;
  int dis = 0;
  int exx = 0;
  int exy = 0;
  int flag = 0;
  /* flag will remain zero if no nonzero combinations of mask1 and mask2 are
   * found.
   */
  double denomx;
  double denomy;
  double tau;
  int i, j;
  if (transpose==0)
  { for (i = 0; i < n; i++)
    { if (mask1[index1][i] && mask2[index2][i])
      { for (j = 0; j < i; j++)
        { if (mask1[index1][j] && mask2[index2][j])
          { double x1 = data1[index1][i];
            double x2 = data1[index1][j];
            double y1 = data2[index2][i];
            double y2 = data2[index2][j];
            if (x1 < x2 && y1 < y2) con++;
            if (x1 > x2 && y1 > y2) con++;
            if (x1 < x2 && y1 > y2) dis++;
            if (x1 > x2 && y1 < y2) dis++;
            if (x1 == x2 && y1 != y2) exx++;
            if (x1 != x2 && y1 == y2) exy++;
            flag = 1;
          }
        }
      }
    }
  }
  else
  { for (i = 0; i < n; i++)
    { if (mask1[i][index1] && mask2[i][index2])
      { for (j = 0; j < i; j++)
        { if (mask1[j][index1] && mask2[j][index2])
          { double x1 = data1[i][index1];
            double x2 = data1[j][index1];
            double y1 = data2[i][index2];
            double y2 = data2[j][index2];
            if (x1 < x2 && y1 < y2) con++;
            if (x1 > x2 && y1 > y2) con++;
            if (x1 < x2 && y1 > y2) dis++;
            if (x1 > x2 && y1 < y2) dis++;
            if (x1 == x2 && y1 != y2) exx++;
            if (x1 != x2 && y1 == y2) exy++;
            flag = 1;
          }
        }
      }
    }
  }
  if (!flag) return 0.;
  denomx = con + dis + exx;
  denomy = con + dis + exy;
  if (denomx==0) return 1;
  if (denomy==0) return 1;
  tau = (con-dis)/sqrt(denomx*denomy);
  return 1.-tau;
}

/* *********************************************************************  */

static
void setmetric (char dist,
  double (**metric)
    (int,double**,double**,int**,int**, const double[],int,int,int) )
{ switch(dist)
  { case ('e'): *metric = &euclid; break;
    case ('h'): *metric = &harmonic; break;
    case ('b'): *metric = &cityblock; break;
    case ('c'): *metric = &correlation; break;
    case ('a'): *metric = &acorrelation; break;
    case ('u'): *metric = &ucorrelation; break;
    case ('x'): *metric = &uacorrelation; break;
    case ('s'): *metric = &spearman; break;
    case ('k'): *metric = &kendall; break;
    default: *metric = &euclid; break;
  }
  return;
}

/* *********************************************************************  */

void CALL initran(void)
/*
Purpose
=======

The routine initran initializes the random number generator using the current
time. The current epoch time in seconds is used as a seed for the standard C
random number generator. The first two random number generated by the standard
C random number generator are then used to initialize the ranlib random number
generator.

External Subroutines:
time.h:     time
ranlib.h:   setall
============================================================================
*/

{ int initseed = time(0);
  int iseed1, iseed2;
  srand(initseed);
  iseed1 = rand();
  iseed2 = rand();
  setall (iseed1, iseed2);
  return;
}

/* ************************************************************************ */

void CALL randomassign (int nclusters, int nelements, int clusterid[])
/*
Purpose
=======

The randomassign routine performs an initial random clustering, needed for
k-means or k-median clustering. Elements (genes or microarrays) are randomly
assigned to clusters. First, nclust elements are randomly chosen to be assigned
to the clusters 0..nclust-1 in order to guarantee that none of the clusters
are empty. The remaining elements are then randomly assigned to a cluster.

Arguments
=========

nclust  (input) int
The number of clusters.

nelements  (input) int
The number of elements to be clustered (i.e., the number of genes or microarrays
to be clustered).

clusterid  (output) int array, dimension( nelements )
The cluster number to which an element was assigned.

External Functions:
ranlib: int genprm
============================================================================
*/

{ int i;
  long* map = malloc(nelements*sizeof(long));
  /* Initialize mapping */
  for (i = 0; i < nelements; i++) map[i] = i;
  /* Create a random permutation of this mapping */
  genprm (map, nelements);

  /* Assign each of the first nclusters elements to a different cluster
   * to avoid empty clusters */
  for (i = 0; i < nclusters; i++) clusterid[map[i]] = i;

  /* Assign other elements randomly to a cluster */
  for (i = nclusters; i < nelements; i++)
    clusterid[map[i]] = ignuin (0,nclusters-1);
  free(map);
  return;
}

/* ********************************************************************* */

void getclustermean(int nclusters, int nrows, int ncolumns,
  double** data, int** mask, int clusterid[], double** cdata, int** cmask,
  int transpose)
/*
Purpose
=======

The getclustermean routine calculates the cluster centroids, given to which
cluster each element belongs. The centroid is defined as the mean over all
elements for each dimension.

Arguments
=========

nclusters  (input) int
The number of clusters.

nrows     (input) int
The number of rows in the gene expression data matrix, equal to the number of
genes.

ncolumns  (input) int
The number of columns in the gene expression data matrix, equal to the number of
microarrays.

data       (input) double array, dimension( nrows,ncolumns )
The array containing the gene expression data.

mask       (input) int array, dimension( nrows,ncolumns )
This array shows which data values are missing. If
mask[i][j] == 0, then data[i][j] is missing.

clusterid  (output) int array, dimension( nrows or ncolumns )
The cluster number to which each element belongs. If transpose==0, then the
dimension of clusterid is equal to nrows (the number of genes). Otherwise, it
is equal to ncolumns (the number of microarrays).

cdata      (output) double array, dimension( nclusters,ncolumns ) (transpose==0)
                               or dimension( nrows, nclusters) (transpose==1)
On exit of getclustermean, this array contains the cluster centroids.

cmask      (output) int array, dimension( nclusters,ncolumns ) (transpose==0)
                            or dimension( nrows, nclusters) (transpose==1)
This array shows which data values of are missing for each centroid. If
cmask[i][j] == 0, then cdata[i][j] is missing. A data value is missing for a
centroid if the corresponding data values of the cluster members are all
missing.

transpose  (input) int
If transpose==0, clusters of rows (genes) are specified. Otherwise, clusters of
columns (microarrays) are specified.

========================================================================
*/
{ int i, j, k;
  if (transpose==0)
  { int** count = malloc(nclusters*sizeof(int*));
    for (i = 0; i < nclusters; i++)
    { count[i] = calloc(ncolumns,sizeof(int));
      for (j = 0; j < ncolumns; j++) cdata[i][j] = 0.;
    }
    for (k = 0; k < nrows; k++)
    { i = clusterid[k];
      for (j = 0; j < ncolumns; j++)
        if (mask[k][j] != 0)
        { cdata[i][j] = cdata[i][j] + data[k][j];
          count[i][j] = count[i][j] + 1;
        }
    }
    for (i = 0; i < nclusters; i++)
    { for (j = 0; j < ncolumns; j++)
      { if (count[i][j]>0)
        { cdata[i][j] = cdata[i][j] / count[i][j];
          cmask[i][j] = 1;
        }
        else
          cmask[i][j] = 0;
      }
      free (count[i]);
    }
    free (count);
  }
  else
  { int** count = malloc(nrows*sizeof(int*));
    for (i = 0; i < nrows; i++)
    { count[i] = calloc(nclusters,sizeof(int));
      for (j = 0; j < nclusters; j++) cdata[i][j] = 0.;
    }
    for (k = 0; k < ncolumns; k++)
    { i = clusterid[k];
      for (j = 0; j < nrows; j++)
      { if (mask[j][k] != 0)
        { cdata[j][i] = cdata[j][i] + data[j][k];
          count[j][i] = count[j][i] + 1;
        }
      }
    }
    for (i = 0; i < nrows; i++)
    { for (j = 0; j < nclusters; j++)
      { if (count[i][j]>0)
        { cdata[i][j] = cdata[i][j] / count[i][j];
          cmask[i][j] = 1;
        }
        else
          cmask[i][j] = 0;
      }
      free (count[i]);
    }
    free (count);
  }
  return;
}

/* ********************************************************************* */

void getclustermedian(int nclusters, int nrows, int ncolumns,
  double** data, int** mask, int clusterid[], double** cdata, int** cmask,
  int transpose)
/*
Purpose
=======

The getclustermedian routine calculates the cluster centroids, given to which
cluster each element belongs. The centroid is defined as the median over all
elements for each dimension.

Arguments
=========

nclusters  (input) int
The number of clusters.

nrows     (input) int
The number of rows in the gene expression data matrix, equal to the number of
genes.

ncolumns  (input) int
The number of columns in the gene expression data matrix, equal to the number of
microarrays.

data       (input) double array, dimension( nrows,ncolumns )
The array containing the gene expression data.

mask       (input) int array, dimension( nrows,ncolumns )
This array shows which data values are missing. If
mask[i][j] == 0, then data[i][j] is missing.

clusterid  (output) int array, dimension( nrows or ncolumns )
The cluster number to which each element belongs. If transpose==0, then the
dimension of clusterid is equal to nrows (the number of genes). Otherwise, it
is equal to ncolumns (the number of microarrays).

cdata      (output) double array, dimension( nclusters,ncolumns ) (transpose==0)
                               or dimension( nrows, nclusters) (transpose==1)
On exit of getclustermedian, this array contains the cluster centroids.

cmask      (output) int array, dimension( nclusters,ncolumns ) (transpose==0)
                            or dimension( nrows, nclusters) (transpose==1)
This array shows which data values of are missing for each centroid. If
cmask[i][j] == 0, then cdata[i][j] is missing. A data value is missing for a
centroid if the corresponding data values of the cluster members are all
missing.

transpose  (input) int
If transpose==0, clusters of rows (genes) are specified. Otherwise, clusters of
columns (microarrays) are specified.

========================================================================
*/
{ int i, j, k;
  if (transpose==0)
  { double* temp = malloc(nrows*sizeof(double));
    for (i = 0; i < nclusters; i++)
    { for (j = 0; j < ncolumns; j++)
      { int count = 0;
        for (k = 0; k < nrows; k++)
          if (i==clusterid[k] && mask[k][j])
          { temp[count] = data[k][j];
            count++;
          }
        if (count>0)
        { cdata[i][j] = median (count,temp);
          cmask[i][j] = 1;
        }
        else
        { cdata[i][j] = 0.;
          cmask[i][j] = 0;
        }
      }
    }
    free (temp);
  }
  else
  { double* temp = malloc(ncolumns*sizeof(double));
    for (i = 0; i < nclusters; i++)
    { for (j = 0; j < nrows; j++)
      { int count = 0;
        for (k = 0; k < ncolumns; k++)
          if (i==clusterid[k] && mask[j][k])
          { temp[count] = data[j][k];
            count++;
          }
        if (count>0)
        { cdata[j][i] = median (count,temp);
          cmask[j][i] = 1;
        }
        else
        { cdata[j][i] = 0.;
          cmask[j][i] = 0;
        }
      }
    }
    free (temp);
  }
  return;
}

/* ********************************************************************* */

void getclustermedoid(int nclusters, int nelements, double** distance,
  int clusterid[], int centroids[], double errors[])
/*
Purpose
=======

The getclustermedoid routine calculates the cluster centroids, given to which
cluster each element belongs. The centroid is defined as the element with the
smallest sum of distances to the other elements.

Arguments
=========

nclusters  (input) int
The number of clusters.

nelements  (input) int
The total number of elements.

distmatrix (input) double array, ragged
  (number of rows is nelements, number of columns is equal to the row number)
The distance matrix. To save space, the distance matrix is given in the
form of a ragged array. The distance matrix is symmetric and has zeros
on the diagonal. See distancematrix for a description of the content.

clusterid  (output) int array, dimension( nelements )
The cluster number to which each element belongs.

centroid   (output) int array, dimension( nclusters )
The index of the element that functions as the centroid for each cluster.

errors     (output) double array, dimension( nclusters )
The within-cluster sum of distances between the items and the cluster
centroid.

========================================================================
*/
{ int i, j, k;
  for (j = 0; j < nclusters; j++) errors[j] = DBL_MAX;
  for (i = 0; i < nelements; i++)
  { double d = 0.0;
    j = clusterid[i];
    for (k = 0; k < nelements; k++)
    { if (i==k || clusterid[k]!=j) continue;
      d += (i < k ? distance[k][i] : distance[i][k]);
      if (d > errors[j]) break;
    }
    if (d < errors[j])
    { errors[j] = d;
      centroids[j] = i;
    }
  }
}

/* ********************************************************************* */

static
void emalg (int nclusters, int nrows, int ncolumns,
  double** data, int** mask, double weight[], int transpose, int init_given,
  void getclustercenter
    (int,int,int,double**,int**,int[],double**,int**,int),
  double metric (int,double**,double**,int**,int**,const double[],int,int,int),
  int clusterid[], double** cdata, int** cmask)

{ const int nobjects = (transpose==0) ? nrows : ncolumns;
  const int ndata = (transpose==0) ? ncolumns : nrows;

  int* cn = calloc(nclusters,sizeof(int));
  /* This will contain the number of elements in each cluster. This is needed
   * to check for empty clusters.
   */

  int* savedids = malloc(nobjects*sizeof(int));
  /* needed to check for periodic behavior */
  int same;

  int changed;
  int iteration = 0;
  int period = 10;
  long* order = malloc(nobjects*sizeof(long));
  int jj;
  for (jj = 0; jj < nobjects; jj++) order[jj] = jj;

  if(!init_given) randomassign (nclusters, nobjects, clusterid);

  for (jj = 0; jj < nobjects; jj++)
  { int ii = clusterid[jj];
    cn[ii]++;
  }

  /* Start the loop */
  do
  { int ii;
    if (iteration % period == 0)
    { /* save the current clustering solution */
      for (ii = 0; ii < nobjects; ii++) savedids[ii] = clusterid[ii];
      period = period * 2;
    }
    iteration += 1;

    /* Find the center */
    getclustercenter (nclusters, nrows, ncolumns, data, mask,
                      clusterid, cdata, cmask, transpose);

    /* Create a random order (except if the user specified an initial
     * clustering, in which case we run the algorithm fully
     * deterministically.  */
    if (!init_given) genprm (order, nobjects);

    changed = 0;

    for (ii = 0; ii < nobjects; ii++)
    /* Calculate the distances */
    { int i = order[ii];
      int jnow = clusterid[i];
      if (cn[jnow]>1)
      { /* No reassignment if that would lead to an empty cluster */
        /* Treat the present cluster as a special case */
        double distance =
          metric(ndata,data,cdata,mask,cmask,weight,i,jnow,transpose);
        int j;
        for (j = 0; j < jnow; j++)
        { double tdistance =
            metric(ndata,data,cdata,mask,cmask,weight,i,j,transpose);
          if (tdistance < distance)
          { distance = tdistance;
            cn[clusterid[i]]--;
            clusterid[i] = j;
            cn[j]++;
            changed = 1;
          }
        }
        for (j = jnow+1; j < nclusters; j++)
        { double tdistance =
            metric(ndata,data,cdata,mask,cmask,weight,i,j,transpose);
          if (tdistance < distance)
          { distance = tdistance;
            cn[clusterid[i]]--;
            clusterid[i] = j;
            cn[j]++;
            changed = 1;
          }
        }
      }
    }
    /* compare to the saved clustering solution */
    same = 1;
    for (ii = 0; ii < nobjects; ii++)
    { if (savedids[ii] != clusterid[ii])
      { same = 0;
        break;   /* No point in checking the other ids */
      }
    }
  } while (changed && !same);
  free (savedids);
  free (order);
  free (cn);
  return;
}

/* *********************************************************************** */

void CALL kcluster (int nclusters, int nrows, int ncolumns,
  double** data, int** mask, double weight[], int transpose,
  int npass, char method, char dist,
  int clusterid[], double** cdata, double* error, int* ifound)
/*
Purpose
=======

The kcluster routine performs k-means or k-median clustering on a given set of
elements, using the specified distance measure. The number of clusters is given
by the user. Multiple passes are being made to find the optimal clustering
solution, each time starting from a different initial clustering.


Arguments
=========

nclusters  (input) int
The number of clusters to be found.

data       (input) double array, dimension( nrows,ncolumns )
The array containing the data of the elements to be clustered (i.e., the gene
expression data).

mask       (input) int array, dimension( nrows,ncolumns )
This array shows which data values are missing. If
mask[i][j] == 0, then data[i][j] is missing.

nrows     (input) int
The number of rows in the data matrix, equal to the number of genes.

ncolumns  (input) int
The number of columns in the data matrix, equal to the number of microarrays.

weight (input) double array, dimension( n )
The weights that are used to calculate the distance.

transpose  (input) int
If transpose==0, the rows of the matrix are clustered. Otherwise, columns
of the matrix are clustered.

npass      (input) int
The number of times clustering is performed. Clustering is performed npass
times, each time starting from a different (random) initial assignment of 
genes to clusters. The clustering solution with the lowest within-cluster sum
of distances is chosen.
If npass==0, then the clustering algorithm will be run once, where the initial
assignment of elements to clusters is taken from the clusterid array.

method     (input) char
Defines whether the arithmetic mean (method=='a') or the median
(method=='m') is used to calculate the cluster center.

dist       (input) char
Defines which distance measure is used, as given by the table:
dist=='e': Euclidean distance
dist=='h': Harmonically summed Euclidean distance
dist=='b': City-block distance
dist=='c': correlation
dist=='a': absolute value of the correlation
dist=='u': uncentered correlation
dist=='x': absolute uncentered correlation
dist=='s': Spearman's rank correlation
dist=='k': Kendall's tau
For other values of dist, the default (Euclidean distance) is used.

clusterid  (output; input) int array, dimension( nrows or ncolumns )
The cluster number to which a gene or microarray was assigned. If npass==0,
then on input clusterid contains the initial clustering assignment from which
the clustering algorithm starts. On output. it contains the clustering solution
that was found.

cdata      (output) double array, dimension( nclusters,ncolumns ) (transpose==0)
                               or dimension( nrows, nclusters) (transpose==1)
This array contains the center of each cluster, as defined
as the average of the elements for each cluster (if
method=='a') or as the median (if method=='m').

error      (output) double
The sum of distances to the cluster center of each item in the optimal k-means
clustering solution that was found.

ifound     (output) int
The number of times the optimal clustering solution was
found. The value of ifound is at least 1; its maximum value is npass.

========================================================================
*/
{ const int nobjects = (transpose==0) ? nrows : ncolumns;
  const int ndata = (transpose==0) ? ncolumns : nrows;
  void (*getclustercenter)
    (int,int,int,double**,int**,int[],double**,int**,int);
  double (*metric)
    (int,double**,double**,int**,int**,const double[],int,int,int);
  const int init_given = (npass==0) ? 1 : 0;

  int i;
  int** cmask;
  int** tcmask;
  double** tcdata;
  int ipass;
  int* tclusterid;
  int* mapping;
  int* savedinitialid = NULL;

  if (nobjects < nclusters)
  { *ifound = 0;
    return;
  }
  /* More clusters asked for than objects available */

  /* First initialize the random number generator */
  initran();

  /* Set the function to find the centroid as indicated by method */
  if (method == 'm') getclustercenter = &getclustermedian;
  else getclustercenter = &getclustermean;

  /* Set the metric function as indicated by dist */
  setmetric (dist, &metric);

  /* Set the result of the first pass as the initial best clustering solution */
  *ifound = 1;

  /* Find out if the user specified an initial clustering */
  if (init_given)
  /* Save the initial clustering specified by the user */
  { savedinitialid = malloc(nobjects*sizeof(int));
    for (i = 0; i < nobjects; i++) savedinitialid[i] = clusterid[i];
  }

  if (transpose==0)
  { cmask = malloc(nclusters*sizeof(int*));
    for (i = 0; i < nclusters; i++) cmask[i] = malloc(ndata*sizeof(int));
  }
  else
  { cmask = malloc(ndata*sizeof(int*));
    for (i = 0; i < ndata; i++) cmask[i] = malloc(nclusters*sizeof(int));
  }

  *error = 0.;
  emalg(nclusters, nrows, ncolumns, data, mask, weight, transpose, init_given,
    getclustercenter, metric, clusterid, cdata, cmask);

  for (i = 0; i < nobjects; i++)
  { int j = clusterid[i];
    *error += metric(ndata, data, cdata, mask, cmask, weight, i, j, transpose);
  }
  if (transpose==0)
    for (i = 0; i < nclusters; i++) free(cmask[i]);
  else
    for (i = 0; i < ndata; i++) free(cmask[i]);
  free(cmask);

  if (npass==0) return;

  /* Create temporary space for cluster centroid information */
  if (transpose==0)
  { tcmask = malloc(nclusters*sizeof(int*));
    tcdata = malloc(nclusters*sizeof(double*));
    for (i = 0; i < nclusters; i++)
    { tcmask[i] = malloc(ndata*sizeof(int));
      tcdata[i] = malloc(ndata*sizeof(double));
    }
  }
  else
  { tcmask = malloc(ndata*sizeof(int*));
    tcdata = malloc(ndata*sizeof(double*));
    for (i = 0; i < ndata; i++)
    { tcmask[i] = malloc(nclusters*sizeof(int));
      tcdata[i] = malloc(nclusters*sizeof(double));
    }
  }

  tclusterid = malloc(nobjects*sizeof(int));
  mapping = malloc(nclusters*sizeof(int));
  for (ipass = 1; ipass < npass; ipass++)
  { double tssin = 0.;
    int same = 1;

    if (init_given)
      for (i = 0; i < nobjects; i++) tclusterid[i] = savedinitialid[i];
    emalg(nclusters, nrows, ncolumns, data, mask, weight, transpose, init_given,
      getclustercenter, metric, tclusterid, tcdata, tcmask);

    for (i = 0; i < nclusters; i++) mapping[i] = -1;
    for (i = 0; i < nobjects; i++)
    { int j = tclusterid[i];
      if (mapping[j] == -1) mapping[j] = clusterid[i];
      else if (mapping[j] != clusterid[i]) same = 0;
      tssin +=
        metric(ndata, data, tcdata, mask, tcmask, weight, i, j, transpose);
    }
    if (same) (*ifound)++;
    else if (tssin < *error)
    { int j;
      *ifound = 1;
      *error = tssin;
      for (i = 0; i < nobjects; i++) clusterid[i] = tclusterid[i];
      if (transpose==0)
      { for (i = 0; i < nclusters; i++)
          for (j = 0; j < ndata; j++)
            cdata[i][j] = tcdata[i][j];
      }
      else
      { for (i = 0; i < ndata; i++)
        { for (j = 0; j < nclusters; j++)
            cdata[i][j] = tcdata[i][j];
        }
      }
    }
  }

  /* Deallocate temporarily used space */
  free(mapping);
  free(tclusterid);
  if (savedinitialid) free(savedinitialid);

  if (transpose==0)
  { for (i = 0; i < nclusters; i++)
    { free(tcmask[i]);
      free(tcdata[i]);
    }
  }
  else
  { for (i = 0; i < ndata; i++)
    { free(tcmask[i]);
      free(tcdata[i]);
    }
  }
  free(tcmask);
  free(tcdata);

  return;
}

/* *********************************************************************** */

void CALL kmedoids (int nclusters, int nelements, double** distance,
  int npass, int clusterid[], double* error, int* ifound)
/*
Purpose
=======

The kmedoids routine performs k-medoids clustering on a given set of elements,
using the distance matrix and the number of clusters passed by the user.
Multiple passes are being made to find the optimal clustering solution, each
time starting from a different initial clustering.


Arguments
=========

nclusters  (input) int
The number of clusters to be found.

nelements  (input) int
The number of elements to be clustered.

distmatrix (input) double array, ragged
  (number of rows is nelements, number of columns is equal to the row number)
The distance matrix. To save space, the distance matrix is given in the
form of a ragged array. The distance matrix is symmetric and has zeros
on the diagonal. See distancematrix for a description of the content.

npass      (input) int
The number of times clustering is performed. Clustering is performed npass
times, each time starting from a different (random) initial assignment of genes
to clusters. The clustering solution with the lowest within-cluster sum of
distances is chosen.
If npass==0, then the clustering algorithm will be run once, where the initial
assignment of elements to clusters is taken from the clusterid array.

clusterid  (output; input) int array, dimension( nelements )
On input, if npass==0, then clusterid contains the initial clustering assignment
from which the clustering algorithm starts; all numbers in clusterid should be
between zero and nelements-1 inclusive. If npass!=0, clusterid is ignored on
input.
On output, clusterid contains the clustering solution that was found: clusterid
contains the number of the cluster to which each item was assigned. On output,
the number of a cluster is defined as the item number of the centroid of the
cluster.

error      (output) double
The sum of distances to the cluster center of each item in the optimal k-medoids
clustering solution that was found.

ifound     (output) int
The number of times the optimal clustering solution was
found. The value of ifound is at least 1; its maximum value is npass.

========================================================================
*/

{ int i, j, k, icluster, ipass;
  int* tclusterid;
  int* centroids;
  int* savedids;
  double* errors;
  int same, changed;
  int iteration = 0;
  int period = 10;
  /* needed to check for periodic behavior */

  if (nelements < nclusters)
  { *ifound = 0;
    return;
  } /* More clusters asked for than elements available */

  centroids = malloc(nclusters*sizeof(int));
  savedids = malloc(nelements*sizeof(int));
  errors = malloc(nclusters*sizeof(double));

  /* Set the result of the first pass as the initial best clustering solution */
  *ifound = 1;

  /* Find out if the user specified an initial clustering */
  if (npass)
  { initran(); /* First initialize the random number generator */
    randomassign (nclusters, nelements, clusterid);
    /* Ready for the first run */
  }

  *error = 0.;
  do /* Start the loop */
  { if (iteration % period == 0)
    { /* save the current clustering solution */
      for (i = 0; i < nelements; i++) savedids[i] = clusterid[i];
      period *= 2;
    }
    iteration++;

    /* Find the center */
    getclustermedoid (nclusters, nelements, distance, clusterid, centroids, errors);

    changed = 0;
    for (i = 0; i < nelements; i++)
    /* Find the closest cluster */
    { double d = DBL_MAX;
      for (icluster = 0; icluster < nclusters; icluster++)
      { double td;
        j = centroids[icluster];
        if (i==j)
        { d = 0.0;
          clusterid[i] = icluster;
          changed = 1;
          break;
        }
        td = (i > j) ? distance[i][j] : distance[j][i];
        if (td < d)
        { d = td;
          clusterid[i] = icluster;
          changed = 1;
        }
      }
    }
    /* compare to the saved clustering solution (periodicity check) */
    same = 1;
    for (i = 0; i < nelements; i++)
    { if (savedids[i] != clusterid[i])
      { same = 0;
        break;   /* No point in checking the other ids */
      }
    }
  } while (changed && !same);

  for (i = 0; i < nelements; i++)
  { const int j = centroids[clusterid[i]];
    /* Set the cluster number to the item number of the cluster centroid */
    clusterid[i] = j;
    if (i==j) continue;
    *error += (i > j) ? distance[i][j] : distance[j][i];
  }
  if (npass==0)
  /* Deterministic result depending on the specified initial clustering */
  { free(savedids);
    free(centroids);
    free(errors);
    return; /* Done for today */
  }

  tclusterid = malloc(nelements*sizeof(int));
  for (ipass = 1; ipass < npass; ipass++)
  { double terror = 0.0;
    same = 1;

    iteration = 0;
    period = 10;
 
    randomassign (nclusters, nelements, tclusterid);
    do /* Start the loop */
    { if (iteration % period == 0)
      { /* save the current clustering solution */
        for (i = 0; i < nelements; i++) savedids[i] = tclusterid[i];
        period = period * 2;
      }
      iteration++;

      /* Find the center */
      getclustermedoid (nclusters, nelements, distance, tclusterid, centroids, errors);

      changed = 0;
      for (i = 0; i < nelements; i++)
      /* Find the closest cluster */
      { double d = DBL_MAX;
        for (icluster = 0; icluster < nclusters; icluster++)
        { double td;
          j = centroids[icluster];
          if (i==j)
          { d = 0.0;
            tclusterid[i] = icluster;
            changed = 1;
            break;
          }
          td = (i > j) ? distance[i][j] : distance[j][i];
          if (td < d)
          { d = td;
            tclusterid[i] = icluster;
            changed = 1;
          }
        }
      }
      /* compare to the saved clustering solution */
      same = 1;
      for (i = 0; i < nelements; i++)
      { if (savedids[i] != tclusterid[i])
        { same = 0;
          break;   /* No point in checking the other ids */
        }
      }
    } while (changed && !same);

    same = 1;
    for (i = 0; i < nelements; i++)
    { k = tclusterid[i];
      j = centroids[k];
      if (j!=clusterid[i]) same = 0;
      if (i==j) continue;
      terror += (i > j) ? distance[i][j] : distance[j][i];
    }
    if (same) (*ifound)++;
    else if (terror < *error)
    { *ifound = 1;
      *error = terror;
      /* The cluster number is set to the item number of the cluster centroid */
      for (i = 0; i < nelements; i++) clusterid[i] = centroids[tclusterid[i]];
    }
  }

  /* Deallocate temporarily used space */
  free(savedids);
  free(centroids);
  free(tclusterid);
  free(errors);

  return;
}

/* ******************************************************************** */

double** CALL distancematrix (int nrows, int ncolumns, double** data,
  int** mask, double weights[], char dist, int transpose)
              
/*
Purpose
=======

The distancematrix routine calculates the distance matrix between genes or
microarrays using their measured gene expression data. Several distance measures
can be used. The routine returns a pointer to a ragged array containing the
distances between the genes. As the distance matrix is symmetric, with zeros on
the diagonal, only the lower triangular half of the distance matrix is saved.
The distancematrix routine allocates space for the distance matrix. If the
parameter transpose is set to a nonzero value, the distances between the columns
(microarrays) are calculated, otherwise distances between the rows (genes) are
calculated.
If sufficient space in memory cannot be allocated to store the distance matrix,
the routine returns a NULL pointer, and all memory allocated so far for the
distance matrix is freed.


Arguments
=========

nrows     (input) int
The number of rows in the gene expression data matrix (i.e., the number of
genes)

ncolumns   (input) int
The number of columns in the gene expression data matrix (i.e., the number of
microarrays)

data       (input) double array, dimension( nrows,ncolumns )
The array containing the gene expression data.

mask       (input) int array, dimension( nrows,ncolumns )
This array shows which data values are missing. If
mask(i,j) == 0, then data(i,j) is missing.

weight (input) double array, dimension( n )
The weights that are used to calculate the distance. The length of this vector
is equal to the number of columns if the distances between genes are calculated,
or the number of rows if the distances between microarrays are calculated.

dist       (input) char
Defines which distance measure is used, as given by the table:
dist=='e': Euclidean distance
dist=='h': Harmonically summed Euclidean distance
dist=='b': City-block distance
dist=='c': correlation
dist=='a': absolute value of the correlation
dist=='u': uncentered correlation
dist=='x': absolute uncentered correlation
dist=='s': Spearman's rank correlation
dist=='k': Kendall's tau
For other values of dist, the default (Euclidean distance) is used.

transpose  (input) int
If transpose is equal to zero, the distances between the rows is
calculated. Otherwise, the distances between the columns is calculated.
The former is needed when genes are being clustered; the latter is used
when microarrays are being clustered.

========================================================================
*/
{ /* First determine the size of the distance matrix */
  const int n = (transpose==0) ? nrows : ncolumns;
  const int ndata = (transpose==0) ? ncolumns : nrows;
  int i,j;
  double** matrix;
  double (*metric)
    (int,double**,double**,int**,int**,const double[],int,int,int);

  if (n < 2) return NULL;

  /* Set up the ragged array */
  matrix = malloc(n*sizeof(double*));
  if(matrix==NULL) return NULL; /* Not enough memory available */
  matrix[0] = NULL;
  /* The zeroth row has zero columns. We allocate it anyway for convenience.*/
  for (i = 1; i < n; i++)
  { matrix[i] = malloc(i*sizeof(double));
    if (matrix[i]==NULL) break; /* Not enough memory available */
  }
  if (i < n) /* break condition encountered */
  { j = i;
    for (i = 1; i < j; i++) free(matrix[i]);
    return NULL;
  }

  /* Set the metric function as indicated by dist */
  setmetric (dist, &metric);

  /* Calculate the distances and save them in the ragged array */
  for (i = 0; i < n; i++)
    for (j = 0; j < i; j++)
      matrix[i][j]=metric(ndata,data,data,mask,mask,weights,i,j,transpose);
  return matrix;
}


/* ******************************************************************** */

static
double getscale(int nelements, double** distmatrix, char dist)

/*
Purpose
=======

The getscale routine finds the value by which the distances should be scaled
such that all distances are between zero and two, as in case of the Pearson
distance.


Arguments
=========

nelements     (input) int
The number of elements to be clustered (i.e., the number of genes or
microarrays).

distmatrix (input) double array, ragged
  (number of rows is nelements, number of columns is equal to the row number)
The distance matrix. To save space, the distance matrix is given in the
form of a ragged array. The distance matrix is symmetric and has zeros
on the diagonal. See distancematrix for a description of the content.

dist       (input) char
Defines which distance measure is used, as given by the table:
dist=='e': Euclidean distance
dist=='h': Harmonically summed Euclidean distance
dist=='b': City-block distance
dist=='c': correlation
dist=='a': absolute value of the correlation
dist=='u': uncentered correlation
dist=='x': absolute uncentered correlation
dist=='s': Spearman's rank correlation
dist=='k': Kendall's tau
For other values of dist, no scaling is done.

========================================================================
*/
{ switch (dist)
  { case 'a':
    case 'x':
      return 0.5;
    case 'e':
    case 'h':
    case 'b':
    { int i,j;
      double maxvalue = 0.;
      for (i = 0; i < nelements; i++)
	for (j = 0; j < i; j++)
	  maxvalue = max(distmatrix[i][j], maxvalue);
      return maxvalue/2.;
    }
  }
  return 1.0;
}

/* ********************************************************************* */

void cuttree (int nelements, int tree[][2], int nclusters, int clusterid[]) 

/*
Purpose
=======

The cuttree routine takes the output of a hierarchical clustering routine, and
divides the elements in the tree structure into clusters based on the
hierarchical clustering result. The number of clusters is specified by the user.

Arguments
=========

nelements      (input) int
The number of elements that were clustered.

tree           (input) int array, dimension( nelements-1,2 )
The clustering solution. Each row in the matrix describes one linking event,
with the two columns containing the name of the nodes that were joined.
The original elements are numbered 0..nelements-1, nodes are numbered
-1..-(nelements-1). The cuttree routine checks the tree array for errors to
avoid segmentation faults. Errors in the tree array that would not cause
segmentation faults may pass undetected. If an error is found, all elements
are assigned to cluster -1, and the routine returns.

nclusters      (input) int
The number of clusters to be formed.

clusterid      (output) int array, dimensions( nelements )
The number of the cluster to which each element was assigned. Space for this
array should be allocated before calling the cuttree routine.

========================================================================
*/
{ int i, j, k;
  int icluster = 0;
  const int n = nelements-nclusters; /* number of nodes to join */
  int* nodeid;
  /* Check the tree */
  int flag = 0;
  if (nclusters > nelements || nclusters < 1) flag = 1;
  for (i = 0; i < nelements-1; i++)
  { if (tree[i][0] >= nelements || tree[i][0] < -i ||
        tree[i][1] >= nelements || tree[i][1] < -i)  
    { flag = 1;
      break;
    }
  }
  /* Assign all elements to cluster -1 and return if an error is found. */
  if (flag)
  { for (i = 0; i < nelements; i++) clusterid[i] = -1;
    return;
  }
  /* The tree array is safe to use. */
  for (i = nelements-2; i >= n; i--)
  { k = tree[i][0];
    if (k>=0)
    { clusterid[k] = icluster;
      icluster++;
    }
    k = tree[i][1];
    if (k>=0)
    { clusterid[k] = icluster;
      icluster++;
    }
  }
  nodeid = malloc(n*sizeof(int));
  for (i = 0; i < n; i++) nodeid[i] = -1;
  for (i = n-1; i >= 0; i--)
  { if(nodeid[i]<0) 
    { j = icluster;
      nodeid[i] = j;
      icluster++;
    }
    else j = nodeid[i];
    k = tree[i][0];
    if (k<0) nodeid[-k-1] = j; else clusterid[k] = j;
    k = tree[i][1];
    if (k<0) nodeid[-k-1] = j; else clusterid[k] = j;
  }
  free(nodeid);
  return;
}

/* ******************************************************************** */

static
void pclcluster (int nrows, int ncolumns, double** data, int** mask,
  double weight[], double** distmatrix, char dist, int transpose,
  int result[][2], double linkdist[])

/*

Purpose
=======

The pclcluster routine performs clustering using pairwise centroid-linking
on a given set of gene expression data, using the distance metric given by dist.

Arguments
=========

nrows     (input) int
The number of rows in the gene expression data matrix, equal to the number of
genes.

ncolumns  (input) int
The number of columns in the gene expression data matrix, equal to the number of
microarrays.

data       (input) double array, dimension( nrows,ncolumns )
The array containing the gene expression data.

mask       (input) int array, dimension( nrows,ncolumns )
This array shows which data values are missing. If
mask[i][j] == 0, then data[i][j] is missing.

weight (input) double array, dimension( n )
The weights that are used to calculate the distance. The length of this vector
is ncolumns if genes are being clustered, and nrows if microarrays are being
clustered.

transpose  (input) int
If transpose==0, the rows of the matrix are clustered. Otherwise, columns
of the matrix are clustered.

dist       (input) char
Defines which distance measure is used, as given by the table:
dist=='e': Euclidean distance
dist=='h': Harmonically summed Euclidean distance
dist=='b': City-block distance
dist=='c': correlation
dist=='a': absolute value of the correlation
dist=='u': uncentered correlation
dist=='x': absolute uncentered correlation
dist=='s': Spearman's rank correlation
dist=='k': Kendall's tau
For other values of dist, the default (Euclidean distance) is used.

distmatrix (input) double**
The distance matrix. This matrix is precalculated by the calling routine
treecluster. The pclcluster routine modifies the contents of distmatrix, but
does not deallocate it.

result  (output) int array, dimension( nelements-1,2 )
The clustering solution. Each row in the matrix describes one linking event,
with the two columns containing the name of the nodes that were joined.
The original genes are numbered 0..ngenes-1, nodes are numbered
-1..-(nelements-1), where nelements is nrows or ncolumns depending on whether
genes (rows) or microarrays (columns) are being clustered.

linkdist (output) double array, dimension(nelements-1)
For each node, the distance between the two subnodes that were joined. The
number of nodes (nnodes) is equal to the number of genes minus one if genes are
clustered, or the number of microarrays minus one if microarrays are clustered.

========================================================================
*/
{ double (*metric)
    (int,double**,double**,int**,int**,const double[],int,int,int);
  int i,j;
  const int nelements = (transpose==0) ? nrows : ncolumns;
  int* distid = malloc(nelements*sizeof(int));
  double** nodedata;
  int** nodecount;
  int inode;
  const int ndata = transpose ? nrows : ncolumns;
  const int nnodes = nelements - 1;

  /* Set the metric function as indicated by dist */
  setmetric (dist, &metric);

  for (i = 0; i < nelements; i++) distid[i] = i;
  /* To remember which row/column in the distance matrix contains what */

  /* Storage for node data */
  if (transpose)
  { nodedata = malloc(ndata*sizeof(double*));
    nodecount = malloc(ndata*sizeof(int*));
    for (i = 0; i < ndata; i++)
    { nodedata[i] = malloc(nnodes*sizeof(double));
      nodecount[i] = malloc(nnodes*sizeof(int));
    }
  }
  else
  { nodedata = malloc(nnodes*sizeof(double*));
    nodecount = malloc(nnodes*sizeof(int*));
    for (i = 0; i < nnodes; i++)
    { nodedata[i] = malloc(ndata*sizeof(double));
      nodecount[i] = malloc(ndata*sizeof(int));
    }
  }

  for (inode = 0; inode < nnodes; inode++)
  { /* Find the pair with the shortest distance */
    int isaved = 1;
    int jsaved = 0;
    double distance = distmatrix[1][0];
    for (i = 0; i < nelements-inode; i++)
      for (j = 0; j < i; j++)
      { if (distmatrix[i][j]<distance)
        { distance = distmatrix[i][j];
          isaved = i;
          jsaved = j;
        }
      }
    result[inode][0] = distid[jsaved];
    result[inode][1] = distid[isaved];
    linkdist[inode] = distance;

    /* Make node jsaved the new node */
    if (transpose)
    { for (i = 0; i < ndata; i++)
      { nodedata[i][inode] = 0.;
        nodecount[i][inode] = 0;
        if (distid[isaved]<0)
        { const int nodecolumn = -distid[isaved]-1;
          const int count = nodecount[i][nodecolumn];
          nodecount[i][inode] += count;
          nodedata[i][inode] += nodedata[i][nodecolumn] * count;
        }
        else
        { const int datacolumn = distid[isaved];
          if (mask[i][datacolumn])
          { nodecount[i][inode]++;
            nodedata[i][inode] += data[i][datacolumn];
          }
        }
        if (distid[jsaved]<0)
        { const int nodecolumn = -distid[jsaved]-1;
          const int count = nodecount[i][nodecolumn];
          nodecount[i][inode] += count;
          nodedata[i][inode] += nodedata[i][nodecolumn] * count;
        }
        else
        { const int datacolumn = distid[jsaved];
          if (mask[i][datacolumn])
          { nodecount[i][inode]++;
            nodedata[i][inode] += data[i][datacolumn];
          }
        }
        if (nodecount[i][inode] > 0) nodedata[i][inode] /= nodecount[i][inode];
      }
    }
    else
    { for (i = 0; i < ndata; i++)
      { nodedata[inode][i] = 0.;
        nodecount[inode][i] = 0;
        if (distid[isaved]<0)
        { const int noderow = -distid[isaved]-1;
          const int count = nodecount[noderow][i];
          nodecount[inode][i] += count;
          nodedata[inode][i] += nodedata[noderow][i] * count;
        }
        else
        { const int datarow = distid[isaved];
          if (mask[datarow][i])
          { nodecount[inode][i]++;
            nodedata[inode][i] += data[datarow][i];
          }
        }
        if (distid[jsaved]<0)
        { const int noderow = -distid[jsaved]-1;
          const int count = nodecount[noderow][i];
          nodecount[inode][i] += count;
          nodedata[inode][i] += nodedata[noderow][i] * count;
        }
        else
        { const int datarow = distid[jsaved];
          if (mask[datarow][i])
          { nodecount[inode][i]++;
            nodedata[inode][i] += data[datarow][i];
          }
        }
        if (nodecount[inode][i] > 0) nodedata[inode][i] /= nodecount[inode][i];
      }
    }
  
    /* Fix the distances */
    distid[isaved] = distid[nnodes-inode];
    for (i = 0; i < isaved; i++)
      distmatrix[isaved][i] = distmatrix[nnodes-inode][i];
    for (i = isaved + 1; i < nnodes-inode; i++)
      distmatrix[i][isaved] = distmatrix[nnodes-inode][i];

    distid[jsaved] = -inode-1;
    for (i = 0; i < jsaved; i++)
    { if (distid[i]<0)
      { distmatrix[jsaved][i] =
          metric(ndata,nodedata,nodedata,nodecount,nodecount,
                 weight,inode,-distid[i]-1,transpose);
      }
      else
      { distmatrix[jsaved][i] =
          metric(ndata,nodedata,data,nodecount,mask,
                 weight,inode,distid[i],transpose);
      }
    }
    for (i = jsaved + 1; i < nnodes-inode; i++)
    { if (distid[i]<0)
      { distmatrix[i][jsaved] =
          metric(ndata,nodedata,nodedata,nodecount,nodecount,
                 weight,inode,-distid[i]-1,transpose);
      }
      else
      { distmatrix[i][jsaved] =
          metric(ndata,nodedata,data,nodecount,mask,
                 weight,inode,distid[i],transpose);
      }
    }
  }

  /* Free temporarily allocated space */
  if (transpose)
  { for (i = 0; i < ndata; i++)
    { free(nodedata[i]);
      free(nodecount[i]);
    }
  }
  else
  { for (i = 0; i < nnodes; i++)
    { free(nodedata[i]);
      free(nodecount[i]);
    }
  }
  free(nodedata);
  free(nodecount);
  free(distid);
 
  return;
}

/* ********************************************************************* */

static
void pslcluster (int nelements, double** distmatrix, int result[][2],
  double linkdist[])
/*

Purpose
=======

The pslcluster routine performs clustering using pairwise single-linking
on the given distance matrix.

Arguments
=========

nelements     (input) int
The number of elements to be clustered.

distmatrix (input) double**
The distance matrix, with nelements rows, each row being filled up to the
diagonal. The elements on the diagonal are not used, as they are assumed to be
zero. The distance matrix will be modified by this routine.

result  (output) int array, dimension( nelements-1,2 )
The clustering solution. Each row in the matrix describes one linking event,
with the two columns containing the name of the nodes that were joined.
The original elements are numbered 0..nelements-1, nodes are numbered
-1..-(nelements-1).

linkdist (output) double array, dimension(nelements-1)
For each node, the distance between the two subnodes that were joined. The
number of nodes (nnodes) is equal to the number of genes minus one if genes are
clustered, or the number of microarrays minus one if microarrays are clustered.

========================================================================
*/
{ int i, j;
  int nNodes;

  /* Setup a list specifying to which cluster a gene belongs */
  int* clusterid = malloc(nelements*sizeof(int));
  for (i = 0; i < nelements; i++) clusterid[i] = i;

  for (nNodes = nelements; nNodes > 1; nNodes--)
  { int isaved = 1;
    int jsaved = 0;
    double distance = distmatrix[1][0];
    for (i = 0; i < nNodes; i++)
    { for (j = 0; j < i; j++)
      { if (distmatrix[i][j] < distance)
        { isaved = i;
          jsaved = j;
          distance = distmatrix[i][j];
        }
      }
    }
    linkdist[nelements-nNodes] = distance;

    /* Fix the distances */
    for (j = 0; j < jsaved; j++)
      distmatrix[jsaved][j] = min(distmatrix[isaved][j],distmatrix[jsaved][j]);
    for (j = jsaved+1; j < isaved; j++)
      distmatrix[j][jsaved] = min(distmatrix[isaved][j],distmatrix[j][jsaved]);
    for (j = isaved+1; j < nNodes; j++)
      distmatrix[j][jsaved] = min(distmatrix[j][isaved],distmatrix[j][jsaved]);

    for (j = 0; j < isaved; j++)
      distmatrix[isaved][j] = distmatrix[nNodes-1][j];
    for (j = isaved+1; j < nNodes-1; j++)
      distmatrix[j][isaved] = distmatrix[nNodes-1][j];

    /* Update clusterids */
    result[nelements-nNodes][0] = clusterid[isaved];
    result[nelements-nNodes][1] = clusterid[jsaved];
    clusterid[jsaved] = nNodes-nelements-1;
    clusterid[isaved] = clusterid[nNodes-1];
  }
  free(clusterid);
 
  return;
}

/* ******************************************************************** */

static
void pmlcluster (int nelements, double** distmatrix, int result[][2],
  double linkdist[])
/*

Purpose
=======

The pmlcluster routine performs clustering using pairwise maximum- (complete-)
linking on the given distance matrix.

Arguments
=========

nelements     (input) int
The number of elements to be clustered.

distmatrix (input) double**
The distance matrix, with nelements rows, each row being filled up to the
diagonal. The elements on the diagonal are not used, as they are assumed to be
zero. The distance matrix will be modified by this routine.

result  (output) int array, dimension( nelements-1,2 )
The clustering solution. Each row in the matrix describes one linking event,
with the two columns containing the name of the nodes that were joined.
The original elements are numbered 0..nelements-1, nodes are numbered
-1..-(nelements-1).

linkdist (output) double array, dimension(nelements-1)
For each node, the distance between the two subnodes that were joined. The
number of nodes (nnodes) is equal to the number of genes minus one if genes are
clustered, or the number of microarrays minus one if microarrays are clustered.

========================================================================
*/
{ int i,j;
  int nNodes;

  /* Setup a list specifying to which cluster a gene belongs */
  int* clusterid = malloc(nelements*sizeof(int));
  for (i = 0; i < nelements; i++) clusterid[i] = i;

  for (nNodes = nelements; nNodes > 1; nNodes--)
  { int isaved = 1;
    int jsaved = 0;
    double distance = distmatrix[1][0];
    for (i = 0; i < nNodes; i++)
      for (j = 0; j < i; j++)
      { if (distmatrix[i][j] < distance)
        { isaved = i;
          jsaved = j;
          distance = distmatrix[i][j];
        }
      }
    linkdist[nelements-nNodes] = distance;

    /* Fix the distances */
    for (j = 0; j < jsaved; j++)
      distmatrix[jsaved][j] = max(distmatrix[isaved][j],distmatrix[jsaved][j]);
    for (j = jsaved+1; j < isaved; j++)
      distmatrix[j][jsaved] = max(distmatrix[isaved][j],distmatrix[j][jsaved]);
    for (j = isaved+1; j < nNodes; j++)
      distmatrix[j][jsaved] = max(distmatrix[j][isaved],distmatrix[j][jsaved]);

    for (j = 0; j < isaved; j++)
      distmatrix[isaved][j] = distmatrix[nNodes-1][j];
    for (j = isaved+1; j < nNodes-1; j++)
      distmatrix[j][isaved] = distmatrix[nNodes-1][j];

    /* Update clusterids */
    result[nelements-nNodes][0] = clusterid[isaved];
    result[nelements-nNodes][1] = clusterid[jsaved];
    clusterid[jsaved] = nNodes-nelements-1;
    clusterid[isaved] = clusterid[nNodes-1];
  }
  free(clusterid);

  return;
}

/* ******************************************************************* */

static
void palcluster (int nelements, double** distmatrix, int result[][2],
  double linkdist[])
/*

Purpose
=======

The palcluster routine performs clustering using pairwise average
linking on the given distance matrix.

Arguments
=========

nelements     (input) int
The number of elements to be clustered.

distmatrix (input) double**
The distance matrix, with nelements rows, each row being filled up to the
diagonal. The elements on the diagonal are not used, as they are assumed to be
zero. The distance matrix will be modified by this routine.

result  (output) int array, dimension( nelements-1,2 )
The clustering solution. Each row in the matrix describes one linking event,
with the two columns containing the name of the nodes that were joined.
The original elements are numbered 0..nelements-1, nodes are numbered
-1..-(nelements-1).

linkdist (output) double array, dimension(nelements-1)
For each node, the distance between the two subnodes that were joined. The
number of nodes (nnodes) is equal to the number of genes minus one if genes are
clustered, or the number of microarrays minus one if microarrays are clustered.

========================================================================
*/
{ int i,j;
  int nNodes;

  /* Keep track of the number of elements in each cluster
   * (needed to calculate the average) */
  int* number = malloc(nelements*sizeof(int));
  /* Setup a list specifying to which cluster a gene belongs */
  int* clusterid = malloc(nelements*sizeof(int));
  for (i = 0; i < nelements; i++)
  { number[i] = 1;
    clusterid[i] = i;
  }

  for (nNodes = nelements; nNodes > 1; nNodes--)
  { int sum;
    int isaved = 1;
    int jsaved = 0;
    double distance = distmatrix[1][0];
    for (i = 0; i < nNodes; i++)
      for (j = 0; j < i; j++)
      { if (distmatrix[i][j] < distance)
        { isaved = i;
          jsaved = j;
          distance = distmatrix[i][j];
        }
      }

    /* Save result */
    result[nelements-nNodes][0] = clusterid[isaved];
    result[nelements-nNodes][1] = clusterid[jsaved];
    linkdist[nelements-nNodes] = distance;

    /* Fix the distances */
    sum = number[isaved] + number[jsaved];
    for (j = 0; j < jsaved; j++)
    { distmatrix[jsaved][j] = distmatrix[isaved][j]*number[isaved]
                            + distmatrix[jsaved][j]*number[jsaved];
      distmatrix[jsaved][j] /= sum;
    }
    for (j = jsaved+1; j < isaved; j++)
    { distmatrix[j][jsaved] = distmatrix[isaved][j]*number[isaved]
                            + distmatrix[j][jsaved]*number[jsaved];
      distmatrix[j][jsaved] /= sum;
    }
    for (j = isaved+1; j < nNodes; j++)
    { distmatrix[j][jsaved] = distmatrix[j][isaved]*number[isaved]
                            + distmatrix[j][jsaved]*number[jsaved];
      distmatrix[j][jsaved] /= sum;
    }

    for (j = 0; j < isaved; j++)
      distmatrix[isaved][j] = distmatrix[nNodes-1][j];
    for (j = isaved+1; j < nNodes-1; j++)
      distmatrix[j][isaved] = distmatrix[nNodes-1][j];

    /* Update number of elements in the clusters */
    number[jsaved] = sum;
    number[isaved] = number[nNodes-1];

    /* Update clusterids */
    clusterid[jsaved] = nNodes-nelements-1;
    clusterid[isaved] = clusterid[nNodes-1];
  }
  free(clusterid);
  free(number);

  return;
}

/* ******************************************************************* */

void CALL treecluster (int nrows, int ncolumns, double** data, int** mask,
  double weight[], int applyscale, int transpose, char dist, char method,
  int result[][2], double linkdist[], double** distmatrix)
/*
Purpose
=======

The treecluster routine performs hierarchical clustering using pairwise
single-, maximum-, centroid-, or average-linkage, as defined by method, on a
given set of gene expression data, using the distance metric given by dist.

Arguments
=========

nrows     (input) int
The number of rows in the data matrix, equal to the number of genes.

ncolumns  (input) int
The number of columns in the data matrix, equal to the number of microarrays.

data       (input) double array, dimension( nrows,ncolumns )
The array containing the data of the vectors to be clustered.

mask       (input) int array, dimension( nrows,ncolumns )
This array shows which data values are missing. If
mask[i][j] == 0, then data[i][j] is missing.

weight (input) double array, dimension( n )
The weights that are used to calculate the distance.

applyscale      (input) int
If applyscale is nonzero, then the distances in linkdist are scaled such
that all distances are between zero and two, as in case of the Pearson
distance. Otherwise, no scaling is applied.

transpose  (input) int
If transpose==0, the rows of the matrix are clustered. Otherwise, columns
of the matrix are clustered.

dist       (input) char
Defines which distance measure is used, as given by the table:
dist=='e': Euclidean distance
dist=='h': Harmonically summed Euclidean distance
dist=='b': City-block distance
dist=='c': correlation
dist=='a': absolute value of the correlation
dist=='u': uncentered correlation
dist=='x': absolute uncentered correlation
dist=='s': Spearman's rank correlation
dist=='k': Kendall's tau
For other values of dist, the default (Euclidean distance) is used.

method     (input) char
Defines which hierarchical clustering method is used:
method=='s': pairwise single-linkage clustering
method=='m': pairwise maximum- (or complete-) linkage clustering
method=='a': pairwise average-linkage clustering
method=='c': pairwise centroid-linkage clustering
For the first three, either the distance matrix or the gene expression data is
sufficient to perform the clustering algorithm. For pairwise centroid-linkage
clustering, however, the gene expression data are always needed, even if the
distance matrix itself is available.

result  (output) int array, dimension( nelements-1,2 )
The clustering solution. Each row in the matrix describes one linking event,
with the two columns containing the name of the nodes that were joined.
The original elements are numbered 0..nelements-1, nodes are numbered
-1..-(nelements-1), where nelements is nrows or ncolumns depending on whether
genes (rows) or microarrays (columns) are being clustered.
If the treecluster routine fails due to lack of memory, the two columns of the
first row of result are set to (0,0) before returning.

linkdist (output) double array, dimension(nelements-1)
For each node, the distance between the two subnodes that were joined. The
number of nodes (nnodes) is equal to the number of genes minus one if genes are
clustered, or the number of microarrays minus one if microarrays are clustered.

distmatrix (input) double**
The distance matrix. If the distance matrix is zero initially, the distance
matrix will be allocated and calculated from the data by treecluster, and
deallocated before treecluster returns. If the distance matrix is passed by the
calling routine, treecluster will modify the contents of the distance matrix as
part of the clustering algorithm, but will not deallocate it. The calling
routine should deallocate the distance matrix after the return from treecluster.

========================================================================
*/
{ const int nelements = (transpose==0) ? nrows : ncolumns;
  const int ldistmatrix = (distmatrix==NULL) ? 0 : 1;
  int i;

  if (nelements < 2) return;

  /* Calculate the distance matrix if the user didn't give it */
  if(!ldistmatrix)
    distmatrix =
      distancematrix (nrows, ncolumns, data, mask, weight, dist, transpose);
  if (!distmatrix) /* Insufficient memory */
  { /* Set the first clustering result to (0,0) to indicate the memory error */
    result[0][0] = 0;
    result[0][1] = 0;
    return;
  }

  switch(method)
  { case 's':
      pslcluster(nelements, distmatrix, result, linkdist);
      break;
    case 'm':
      pmlcluster(nelements, distmatrix, result, linkdist);
      break;
    case 'a':
      palcluster(nelements, distmatrix, result, linkdist);
      break;
    case 'c':
      pclcluster(nrows, ncolumns, data, mask, weight, distmatrix, dist,
		transpose, result, linkdist);
      break;
  }

  /* Scale the distances in linkdist if so requested */
  if (applyscale)
  { double scale = getscale(nelements, distmatrix, dist);
    for (i = 0; i < nelements-1; i++) linkdist[i] /= scale;
  }

  /* Deallocate space for distance matrix, if it was allocated by treecluster */
  if (!ldistmatrix)
  { for (i = 1; i < nelements; i++) free(distmatrix[i]);
    free (distmatrix);
  }
 
  return;
}

/* ******************************************************************* */

static
void somworker (int nrows, int ncolumns, double** data, int** mask,
  const double weights[], int transpose, int nxgrid, int nygrid,
  double inittau, double*** celldata, int niter, char dist)

{ const int nelements = (transpose==0) ? nrows : ncolumns;
  const int ndata = (transpose==0) ? ncolumns : nrows;
  double (*metric)
    (int,double**,double**,int**,int**,const double[],int,int,int);
  int i, j;
  double* stddata = calloc(nelements,sizeof(double));
  int** dummymask;
  int ix, iy;
  long* index;
  int iter;
  /* Maximum radius in which nodes are adjusted */
  double maxradius = sqrt(nxgrid*nxgrid+nygrid*nygrid);

  /* Initialize the random number generator */
  initran();

  /* Set the metric function as indicated by dist */
  setmetric (dist, &metric);

  /* Calculate the standard deviation for each row or column */
  if (transpose==0)
  { for (i = 0; i < nelements; i++)
    { for (j = 0; j < ndata; j++)
      { double term = data[i][j];
        term = term * term;
        stddata[i] += term;
      }
      stddata[i] = sqrt(stddata[i]);
      if (stddata[i]==0) stddata[i] = 1;
    }
  }
  else
  { for (i = 0; i < nelements; i++)
    { for (j = 0; j < ndata; j++)
      { double term = data[j][i];
        term = term * term;
        stddata[i] += term;
      }
      stddata[i] = sqrt(stddata[i]);
      if (stddata[i]==0) stddata[i] = 1;
    }
  }

  if (transpose==0)
  { dummymask = malloc(nygrid*sizeof(int*));
    for (i = 0; i < nygrid; i++)
    { dummymask[i] = malloc(ndata*sizeof(int));
      for (j = 0; j < ndata; j++) dummymask[i][j] = 1;
    }
  }
  else
  { dummymask = malloc(ndata*sizeof(int*));
    for (i = 0; i < ndata; i++)
    { dummymask[i] = malloc(sizeof(int));
      dummymask[i][0] = 1;
    }
  }

  /* Randomly initialize the nodes */
  for (ix = 0; ix < nxgrid; ix++)
  { for (iy = 0; iy < nygrid; iy++)
    { double sum = 0.;
      for (i = 0; i < ndata; i++)
      { double term = genunf(-1.,1.);
        celldata[ix][iy][i] = term;
        sum += term * term;
      }
      sum = sqrt(sum);
      for (i = 0; i < ndata; i++) celldata[ix][iy][i] /= sum;
    }
  }

  /* Randomize the order in which genes or arrays will be used */
  index = malloc(nelements*sizeof(long));
  for (i = 0; i < nelements; i++) index[i] = i;
  genprm (index, nelements);

  /* Start the iteration */
  for (iter = 0; iter < niter; iter++)
  { int ixbest = 0;
    int iybest = 0;
    long iobject = iter % nelements;
    iobject = index[iobject];
    if (transpose==0)
    { double closest = metric(ndata,data,celldata[ixbest],
        mask,dummymask,weights,iobject,iybest,transpose);
      double radius = maxradius * (1. - ((double)iter)/((double)niter));
      double tau = inittau * (1. - ((double)iter)/((double)niter));

      for (ix = 0; ix < nxgrid; ix++)
      { for (iy = 0; iy < nygrid; iy++)
        { double distance =
            metric (ndata,data,celldata[ix],
              mask,dummymask,weights,iobject,iy,transpose);
          if (distance < closest)
          { ixbest = ix;
            iybest = iy;
            closest = distance;
          }
        }
      }
      for (ix = 0; ix < nxgrid; ix++)
      { for (iy = 0; iy < nygrid; iy++)
        { if (sqrt((ix-ixbest)*(ix-ixbest)+(iy-iybest)*(iy-iybest))<radius)
          { double sum = 0.;
            for (i = 0; i < ndata; i++)
              celldata[ix][iy][i] +=
                tau * (data[iobject][i]/stddata[iobject]-celldata[ix][iy][i]);
            for (i = 0; i < ndata; i++)
            { double term = celldata[ix][iy][i];
              term = term * term;
              sum += term;
            }
            sum = sqrt(sum);
            if (sum>0)
              for (i = 0; i < ndata; i++) celldata[ix][iy][i] /= sum;
          }
        }
      }
    }
    else
    { double closest;
      double** celldatavector = malloc(ndata*sizeof(double*));
      double radius = maxradius * (1. - ((double)iter)/((double)niter));
      double tau = inittau * (1. - ((double)iter)/((double)niter));

      for (i = 0; i < ndata; i++)
        celldatavector[i] = &(celldata[ixbest][iybest][i]);
      closest = metric(ndata,data,celldatavector,
        mask,dummymask,weights,iobject,0,transpose);
      for (ix = 0; ix < nxgrid; ix++)
      { for (iy = 0; iy < nygrid; iy++)
        { double distance;
          for (i = 0; i < ndata; i++)
            celldatavector[i] = &(celldata[ixbest][iybest][i]);
          distance =
            metric (ndata,data,celldatavector,
              mask,dummymask,weights,iobject,0,transpose);
          if (distance < closest)
          { ixbest = ix;
            iybest = iy;
            closest = distance;
          }
        }
      }
      free(celldatavector);
      for (ix = 0; ix < nxgrid; ix++)
      { for (iy = 0; iy < nygrid; iy++)
        { if (sqrt((ix-ixbest)*(ix-ixbest)+(iy-iybest)*(iy-iybest))<radius)
          { double sum = 0.;
            for (i = 0; i < ndata; i++)
              celldata[ix][iy][i] +=
                tau * (data[i][iobject]/stddata[iobject]-celldata[ix][iy][i]);
            for (i = 0; i < ndata; i++)
            { double term = celldata[ix][iy][i];
              term = term * term;
              sum += term;
            }
            sum = sqrt(sum);
            if (sum>0)
              for (i = 0; i < ndata; i++) celldata[ix][iy][i] /= sum;
          }
        }
      }
    }
  }
  if (transpose==0)
    for (i = 0; i < nygrid; i++) free(dummymask[i]);
  else
    for (i = 0; i < ndata; i++) free(dummymask[i]);
  free(dummymask);
  free(stddata);
  free(index);
  return;
}

/* ******************************************************************* */

static
void somassign (int nrows, int ncolumns, double** data, int** mask,
  const double weights[], int transpose, int nxgrid, int nygrid,
  double*** celldata, char dist, int clusterid[][2])
/* Collect clusterids */
{ double (*metric)
    (int,double**,double**,int**,int**, const double[],int,int,int);
  const int ndata = (transpose==0) ? ncolumns : nrows;
  int i,j;

  setmetric (dist, &metric);

  if (transpose==0)
  { int** dummymask = malloc(nygrid*sizeof(int*));
    for (i = 0; i < nygrid; i++)
    { dummymask[i] = malloc(ncolumns*sizeof(int));
      for (j = 0; j < ncolumns; j++) dummymask[i][j] = 1;
    }
    for (i = 0; i < nrows; i++)
    { int ixbest = 0;
      int iybest = 0;
      double closest = metric(ndata,data,celldata[ixbest],
        mask,dummymask,weights,i,iybest,transpose);
      int ix, iy;
      for (ix = 0; ix < nxgrid; ix++)
      { for (iy = 0; iy < nygrid; iy++)
        { double distance =
            metric (ndata,data,celldata[ix],
              mask,dummymask,weights,i,iy,transpose);
          if (distance < closest)
          { ixbest = ix;
            iybest = iy;
            closest = distance;
          }
        }
      }
      clusterid[i][0] = ixbest;
      clusterid[i][1] = iybest;
    }
    for (i = 0; i < nygrid; i++) free(dummymask[i]);
    free(dummymask);
  }
  else
  { double** celldatavector = malloc(ndata*sizeof(double*));
    int** dummymask = malloc(nrows*sizeof(int*));
    int ixbest = 0;
    int iybest = 0;
    for (i = 0; i < nrows; i++)
    { dummymask[i] = malloc(sizeof(int));
      dummymask[i][0] = 1;
    }
    for (i = 0; i < ncolumns; i++)
    { double closest;
      int ix, iy;
      for (j = 0; j < ndata; j++)
        celldatavector[j] = &(celldata[ixbest][iybest][j]);
      closest = metric(ndata,data,celldatavector,
        mask,dummymask,weights,i,0,transpose);
      for (ix = 0; ix < nxgrid; ix++)
      { for (iy = 0; iy < nygrid; iy++)
        { double distance;
          for(j = 0; j < ndata; j++)
            celldatavector[j] = &(celldata[ix][iy][j]);
          distance = metric(ndata,data,celldatavector,
            mask,dummymask,weights,i,0,transpose);
          if (distance < closest)
          { ixbest = ix;
            iybest = iy;
            closest = distance;
          }
        }
      }
      clusterid[i][0] = ixbest;
      clusterid[i][1] = iybest;
    }
    free(celldatavector);
    for (i = 0; i < nrows; i++) free(dummymask[i]);
    free(dummymask);
  }
  return;
}

/* ******************************************************************* */

void CALL somcluster (int nrows, int ncolumns, double** data, int** mask,
  const double weight[], int transpose, int nxgrid, int nygrid,
  double inittau, int niter, char dist, double*** celldata, int clusterid[][2])
/*

Purpose
=======

The somcluster routine implements a self-organizing map (Kohonen) on a
rectangular grid, using a given set of vectors. The distance measure to be
used to find the similarity between genes and nodes is given by dist.

Arguments
=========

nrows     (input) int
The number of rows in the data matrix, equal to the number of genes.

ncolumns  (input) int
The number of columns in the data matrix, equal to the number of microarrays.

data       (input) double array, dimension( nrows,ncolumns )
The array containing the gene expression data.

mask       (input) int array, dimension( nrows,ncolumns )
This array shows which data values are missing. If
mask[i][j] == 0, then data[i][j] is missing.

weights    (input) double array, dimension( n )
The weights that are used to calculate the distance. The length of this vector
is ncolumns if genes are being clustered, or nrows if microarrays are being
clustered.

transpose  (input) int
If transpose==0, the rows of the matrix are clustered. Otherwise, columns
of the matrix are clustered.

nxgrid    (input) int
The number of grid cells horizontally in the rectangular topology of clusters.

nygrid    (input) int
The number of grid cells horizontally in the rectangular topology of clusters.

inittau    (input) double
The initial value of tau, representing the neighborhood function.

niter      (input) int
The number of iterations to be performed.

dist       (input) char
Defines which distance measure is used, as given by the table:
dist=='e': Euclidean distance
dist=='h': Harmonically summed Euclidean distance
dist=='b': City-block distance
dist=='c': correlation
dist=='a': absolute value of the correlation
dist=='u': uncentered correlation
dist=='x': absolute uncentered correlation
dist=='s': Spearman's rank correlation
dist=='k': Kendall's tau
For other values of dist, the default (Euclidean distance) is used.

celldata (output) double array,
  dimension(nxgrid, nygrid, ncolumns) if genes are being clustered
  dimension(nxgrid, nygrid, nrows) if microarrays are being clustered
The gene expression data for each node (cell) in the 2D grid. This can be
interpreted as the centroid for the cluster corresponding to that cell. If
celldata is NULL, then the centroids are not returned. If celldata is not
NULL, enough space should be allocated to store the centroid data before callingsomcluster.

clusterid (output), int[nrows][2] if genes are being clustered
                    int[ncolumns][2] if microarrays are being clustered
For each item (gene or microarray) that is clustered, the coordinates of the
cell in the 2D grid to which the item was assigned. If clusterid is NULL, the
cluster assignments are not returned. If clusterid is not NULL, enough memory
should be allocated to store the clustering information before calling
somcluster.

========================================================================
*/
{ const int nobjects = (transpose==0) ? nrows : ncolumns;
  const int ndata = (transpose==0) ? ncolumns : nrows;
  int i,j;
  const int lcelldata = (celldata==NULL) ? 0 : 1;

  if (nobjects < 2) return;

  if (lcelldata==0)
  { celldata = malloc(nxgrid*nygrid*ndata*sizeof(double**));
    for (i = 0; i < nxgrid; i++)
    { celldata[i] = malloc(nygrid*ndata*sizeof(double*));
      for (j = 0; j < nygrid; j++)
        celldata[i][j] = malloc(ndata*sizeof(double));
    }
  }

  somworker (nrows, ncolumns, data, mask, weight, transpose, nxgrid, nygrid,
    inittau, celldata, niter, dist);
  if (clusterid)
    somassign (nrows, ncolumns, data, mask, weight, transpose,
      nxgrid, nygrid, celldata, dist, clusterid);
  if(lcelldata==0)
  { for (i = 0; i < nxgrid; i++)
      for (j = 0; j < nygrid; j++)
        free(celldata[i][j]);
    for (i = 0; i < nxgrid; i++)
      free(celldata[i]);
    free(celldata);
  }
  return;
}

/* ******************************************************************** */

double CALL clusterdistance (int nrows, int ncolumns, double** data,
  int** mask, double weight[], int n1, int n2, int index1[], int index2[],
  char dist, char method, int transpose)
              
/*
Purpose
=======

The clusterdistance routine calculates the distance between two clusters
containing genes or microarrays using the measured gene expression vectors. The
distance between clusters, given the genes/microarrays in each cluster, can be
defined in several ways. Several distance measures can be used.

The routine returns the distance in double precision.
If the parameter transpose is set to a nonzero value, the clusters are
interpreted as clusters of microarrays, otherwise as clusters of gene.

Arguments
=========

nrows     (input) int
The number of rows (i.e., the number of genes) in the gene expression data
matrix.

ncolumns      (input) int
The number of columns (i.e., the number of microarrays) in the gene expression
data matrix.

data       (input) double array, dimension( nrows,ncolumns )
The array containing the data of the vectors.

mask       (input) int array, dimension( nrows,ncolumns )
This array shows which data values are missing. If
mask(i,j) == 0, then data(i,j) is missing.

weight     (input) double array, dimension( n )
The weights that are used to calculate the distance.

n1         (input) int
The number of elements in the first cluster.

n2         (input) int
The number of elements in the second cluster.

index1     (input) int array, dimension ( n1 )
Identifies which genes/microarrays belong to the first cluster.

index2     (input) int array, dimension ( n2 )
Identifies which genes/microarrays belong to the second cluster.

dist       (input) char
Defines which distance measure is used, as given by the table:
dist=='e': Euclidean distance
dist=='h': Harmonically summed Euclidean distance
dist=='b': City-block distance
dist=='c': correlation
dist=='a': absolute value of the correlation
dist=='u': uncentered correlation
dist=='x': absolute uncentered correlation
dist=='s': Spearman's rank correlation
dist=='k': Kendall's tau
For other values of dist, the default (Euclidean distance) is used.

method     (input) char
Defines how the distance between two clusters is defined, given which genes
belong to which cluster:
method=='a': the distance between the arithmetic means of the two clusters
method=='m': the distance between the medians of the two clusters
method=='s': the smallest pairwise distance between members of the two clusters
method=='x': the largest pairwise distance between members of the two clusters
method=='v': average of the pairwise distances between members of the clusters

transpose  (input) int
If transpose is equal to zero, the distances between the rows is
calculated. Otherwise, the distances between the columns is calculated.
The former is needed when genes are being clustered; the latter is used
when microarrays are being clustered.

========================================================================
*/
{ double (*metric)
    (int,double**,double**,int**,int**, const double[],int,int,int);
  /* if one or both clusters are empty, return */
  if (n1 < 1 || n2 < 1) return 0;
  /* Check the indeces */
  if (transpose==0)
  { int i;
    for (i = 0; i < n1; i++)
    { int index = index1[i];
      if (index < 0 || index >= nrows) return 0;
    }
    for (i = 0; i < n2; i++)
    { int index = index2[i];
      if (index < 0 || index >= nrows) return 0;
    }
  }
  else
  { int i;
    for (i = 0; i < n1; i++)
    { int index = index1[i];
      if (index < 0 || index >= ncolumns) return 0;
    }
    for (i = 0; i < n2; i++)
    { int index = index2[i];
      if (index < 0 || index >= ncolumns) return 0;
    }
  }
  /* Set the metric function as indicated by dist */
  setmetric (dist, &metric);
  switch (method)
  { case 'a':
    { /* Find the center */
      int i,j,k;
      if (transpose==0)
      { double distance;
        double* cdata[2];
        int* cmask[2];
        int* count[2];
        count[0] = calloc(ncolumns,sizeof(int));
        count[1] = calloc(ncolumns,sizeof(int));
        cdata[0] = calloc(ncolumns,sizeof(double));
        cdata[1] = calloc(ncolumns,sizeof(double));
        cmask[0] = malloc(ncolumns*sizeof(int));
        cmask[1] = malloc(ncolumns*sizeof(int));
        for (i = 0; i < n1; i++)
        { k = index1[i];
          for (j = 0; j < ncolumns; j++)
            if (mask[k][j] != 0)
            { cdata[0][j] = cdata[0][j] + data[k][j];
              count[0][j] = count[0][j] + 1;
            }
        }
        for (i = 0; i < n2; i++)
        { k = index2[i];
          for (j = 0; j < ncolumns; j++)
            if (mask[k][j] != 0)
            { cdata[1][j] = cdata[1][j] + data[k][j];
              count[1][j] = count[1][j] + 1;
            }
        }
        for (i = 0; i < 2; i++)
          for (j = 0; j < ncolumns; j++)
          { if (count[i][j]>0)
            { cdata[i][j] = cdata[i][j] / count[i][j];
              cmask[i][j] = 1;
            }
            else
              cmask[i][j] = 0;
          }
        distance =
          metric (ncolumns,cdata,cdata,cmask,cmask,weight,0,1,0);
        for (i = 0; i < 2; i++)
        { free (cdata[i]);
          free (cmask[i]);
          free (count[i]);
        }
        return distance;
      }
      else
      { double distance;
        int** count = malloc(nrows*sizeof(int*));
        double** cdata = malloc(nrows*sizeof(double*));
        int** cmask = malloc(nrows*sizeof(int*));
        for (i = 0; i < nrows; i++)
        { count[i] = calloc(2,sizeof(int));
          cdata[i] = calloc(2,sizeof(double));
          cmask[i] = malloc(2*sizeof(int));
        }
        for (i = 0; i < n1; i++)
        { k = index1[i];
          for (j = 0; j < nrows; j++)
          { if (mask[j][k] != 0)
            { cdata[j][0] = cdata[j][0] + data[j][k];
              count[j][0] = count[j][0] + 1;
            }
          }
        }
        for (i = 0; i < n2; i++)
        { k = index2[i];
          for (j = 0; j < nrows; j++)
          { if (mask[j][k] != 0)
            { cdata[j][1] = cdata[j][1] + data[j][k];
              count[j][1] = count[j][1] + 1;
            }
          }
        }
        for (i = 0; i < nrows; i++)
          for (j = 0; j < 2; j++)
            if (count[i][j]>0)
            { cdata[i][j] = cdata[i][j] / count[i][j];
              cmask[i][j] = 1;
            }
            else
              cmask[i][j] = 0;
        distance = metric (nrows,cdata,cdata,cmask,cmask,weight,0,1,1);
        for (i = 0; i < nrows; i++)
        { free (count[i]);
          free (cdata[i]);
          free (cmask[i]);
        }
        free (count);
        free (cdata);
        free (cmask);
        return distance;
      }
    }
    case 'm':
    { int i, j, k;
      if (transpose==0)
      { double distance;
        double* temp = malloc(nrows*sizeof(double));
        double* cdata[2];
        int* cmask[2];
        for (i = 0; i < 2; i++)
        { cdata[i] = malloc(ncolumns*sizeof(double));
          cmask[i] = malloc(ncolumns*sizeof(int));
        }
        for (j = 0; j < ncolumns; j++)
        { int count = 0;
          for (k = 0; k < n1; k++)
          { i = index1[k];
            if (mask[i][j])
            { temp[count] = data[i][j];
              count++;
            }
          }
          if (count>0)
          { cdata[0][j] = median (count,temp);
            cmask[0][j] = 1;
          }
          else
          { cdata[0][j] = 0.;
            cmask[0][j] = 0;
          }
        }
        for (j = 0; j < ncolumns; j++)
        { int count = 0;
          for (k = 0; k < n2; k++)
          { i = index2[k];
            if (mask[i][j])
            { temp[count] = data[i][j];
              count++;
            }
          }
          if (count>0)
          { cdata[1][j] = median (count,temp);
            cmask[1][j] = 1;
          }
          else
          { cdata[1][j] = 0.;
            cmask[1][j] = 0;
          }
        }
        distance = metric (ncolumns,cdata,cdata,cmask,cmask,weight,0,1,0);
        for (i = 0; i < 2; i++)
        { free (cdata[i]);
          free (cmask[i]);
        }
        free(temp);
        return distance;
      }
      else
      { double distance;
        double* temp = malloc(ncolumns*sizeof(double));
        double** cdata = malloc(nrows*sizeof(double*));
        int** cmask = malloc(nrows*sizeof(int*));
        for (i = 0; i < nrows; i++)
        { cdata[i] = malloc(2*sizeof(double));
          cmask[i] = malloc(2*sizeof(int));
        }
        for (j = 0; j < nrows; j++)
        { int count = 0;
          for (k = 0; k < n1; k++)
          { i = index1[k];
            if (mask[j][i])
            { temp[count] = data[j][i];
              count++;
            }
          }
          if (count>0)
          { cdata[j][0] = median (count,temp);
            cmask[j][0] = 1;
          }
          else
          { cdata[j][0] = 0.;
            cmask[j][0] = 0;
          }
        }
        for (j = 0; j < nrows; j++)
        { int count = 0;
          for (k = 0; k < n2; k++)
          { i = index2[k];
            if (mask[j][i])
            { temp[count] = data[j][i];
              count++;
            }
          }
          if (count>0)
          { cdata[j][1] = median (count,temp);
            cmask[j][1] = 1;
          }
          else
          { cdata[j][1] = 0.;
            cmask[j][1] = 0;
          }
        }
        distance = metric (nrows,cdata,cdata,cmask,cmask,weight,0,1,1);
        for (i = 0; i < nrows; i++)
        { free (cdata[i]);
          free (cmask[i]);
        }
        free(cdata);
        free(cmask);
        free(temp);
        return distance;
      }
    }
    case 's':
    { int i1, i2, j1, j2;
      const int n = (transpose==0) ? ncolumns : nrows;
      double mindistance = DBL_MAX;
      for (i1 = 0; i1 < n1; i1++)
        for (i2 = 0; i2 < n2; i2++)
        { double distance;
          j1 = index1[i1];
          j2 = index2[i2];
          distance = metric (n,data,data,mask,mask,weight,j1,j2,transpose);
          if (distance < mindistance) mindistance = distance;
        }
      return mindistance;
    }
    case 'x':
    { int i1, i2, j1, j2;
      const int n = (transpose==0) ? ncolumns : nrows;
      double maxdistance = 0;
      for (i1 = 0; i1 < n1; i1++)
        for (i2 = 0; i2 < n2; i2++)
        { double distance;
          j1 = index1[i1];
          j2 = index2[i2];
          distance = metric (n,data,data,mask,mask,weight,j1,j2,transpose);
          if (distance > maxdistance) maxdistance = distance;
        }
      return maxdistance;
    }
    case 'v':
    { int i1, i2, j1, j2;
      const int n = (transpose==0) ? ncolumns : nrows;
      double distance = 0;
      for (i1 = 0; i1 < n1; i1++)
        for (i2 = 0; i2 < n2; i2++)
        { j1 = index1[i1];
          j2 = index2[i2];
          distance += metric (n,data,data,mask,mask,weight,j1,j2,transpose);
        }
      distance /= (n1*n2);
      return distance;
    }
  }
  /* Never get here */
  return 0;
}
