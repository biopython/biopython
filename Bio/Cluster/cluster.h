/******************************************************************************/
/* The C Clustering Library for cDNA microarray data.
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

#ifndef C_CLUSTERING_LIB
#define C_CLUSTERING_LIB

#ifndef CALL 
# define CALL
#endif

#ifndef min
#define min(x, y)	((x) < (y) ? (x) : (y))
#endif
#ifndef max
#define	max(x, y)	((x) > (y) ? (x) : (y))
#endif

#ifdef WINDOWS
#  include <windows.h>
#endif


/* Chapter 2 */
double CALL clusterdistance (int nrows, int ncolumns, double** data, int** mask,
  double weight[], int n1, int n2, int index1[], int index2[], char dist,
  char method, int transpose);

/* Chapter 3 */
void CALL initran(void);

/* Chapter 4 */
double** CALL distancematrix (int ngenes, int ndata, double** data,
  int** mask, double* weight, char dist, int transpose);

/* Chapter 5 */
void CALL randomassign (int nclusters, int ngenes, int clusterid[]);
void getclustermean (int nclusters, int nrows, int ncolumns,
  double** data, int** mask, int clusterid[], double** cdata, int** cmask,
  int transpose);
void getclustermedian (int nclusters, int nrows, int ncolumns,
  double** data, int** mask, int clusterid[], double** cdata, int** cmask,
  int transpose);
void getclustermedoid(int nclusters, int nelements, double** distance,
  int clusterid[], int centroids[], double errors[]);
void CALL kcluster (int nclusters, int ngenes, int ndata, double** data,
  int** mask, double weight[], int transpose, int npass, char method, char dist,
  int clusterid[], double** cdata, double* error, int* ifound);
void CALL kmedoids (int nclusters, int nelements, double** distance,
  int npass, int clusterid[], double* error, int* ifound);

/* Chapter 6 */
void CALL treecluster (int nrows, int ncolumns, double** data, int** mask,
  double weight[], int applyscale, int transpose, char dist, char method,
  int result[][2], double linkdist[], double** distmatrix);
void cuttree (int nelements, int tree[][2], int nclusters, int clusterid[]);

/* Chapter 7 */
void CALL somcluster (int nrows, int ncolumns, double** data, int** mask,
  const double weight[], int transpose, int nxnodes, int nynodes,
  double inittau, int niter, char dist, double*** celldata,
  int clusterid[][2]);

/* Chapter 8 */
void CALL svd(int m, int n, double** u, double w[], double** v, int* ierr);

/* Utility routines, currently undocumented */
void CALL sort(int n, const double data[], int index[]);
double CALL mean(int n, double x[]);
double CALL median (int n, double x[]);
#endif
