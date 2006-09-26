/******************************************************************************/
/* The C Clustering Library for cDNA microarray data.
 * Copyright (C) 2002 Michiel Jan Laurens de Hoon.
 *
 * This library was written at the Laboratory of DNA Information Analysis,
 * Human Genome Center, Institute of Medical Science, University of Tokyo,
 * 4-6-1 Shirokanedai, Minato-ku, Tokyo 108-8639, Japan.
 * Contact: mdehoon@c2b2.columbia.edu
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
double** CALL distancematrix (int ngenes, int ndata, double** data,
  int** mask, double* weight, char dist, int transpose);

/* Chapter 3 */
int getclustercentroids(int nclusters, int nrows, int ncolumns,
  double** data, int** mask, int clusterid[], double** cdata, int** cmask,
  int transpose, char method);
void getclustermedoids(int nclusters, int nelements, double** distance,
  int clusterid[], int centroids[], double errors[]);
void CALL kcluster (int nclusters, int ngenes, int ndata, double** data,
  int** mask, double weight[], int transpose, int npass, char method, char dist,
  int clusterid[], double* error, int* ifound);
void CALL kmedoids (int nclusters, int nelements, double** distance,
  int npass, int clusterid[], double* error, int* ifound);

/* Chapter 4 */
typedef struct {int left; int right; double distance;} Node;
/*
 * A Node struct describes a single node in a tree created by hierarchical
 * clustering. The tree can be represented by an array of n Node structs,
 * where n is the number of elements minus one. The integers left and right
 * in each Node struct refer to the two elements or subnodes that are joined
 * in this node. The original elements are numbered 0..nelements-1, and the
 * nodes -1..-(nelements-1). For each node, distance contains the distance
 * between the two subnodes that were joined.
 */

Node* CALL treecluster (int nrows, int ncolumns, double** data, int** mask,
  double weight[], int transpose, char dist, char method, double** distmatrix);
void cuttree (int nelements, Node* tree, int nclusters, int clusterid[]);

/* Chapter 5 */
void CALL somcluster (int nrows, int ncolumns, double** data, int** mask,
  const double weight[], int transpose, int nxnodes, int nynodes,
  double inittau, int niter, char dist, double*** celldata,
  int clusterid[][2]);

/* Chapter 6 */
void CALL svd(int m, int n, double** u, double w[], double** v, int* ierr);

/* Utility routines, currently undocumented */
void CALL sort(int n, const double data[], int index[]);
double CALL mean(int n, double x[]);
double CALL median (int n, double x[]);

double* calculate_weights(int nrows, int ncolumns, double** data, int** mask,
  double weights[], int transpose, char dist, double cutoff, double exponent);
