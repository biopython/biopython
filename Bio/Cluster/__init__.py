from Numeric import *
import cluster
from cluster import mean, median
from data import readdatafile, writeclusterfiles

def kcluster(data,nclusters=2,mask=None,weight=None,transpose=0,npass=1,method='a',dist='e'):
  """returns clusterid, centroids, error, nfound.

This function implements k-means clustering.
The number of clusters is given by ncluster.
The array data is a nrows x ncolumns array containing the gene expression
data.
The array mask shows which data are missing. If mask[i][j]==0, then
data[i][j] is missing.
The array weight contains the weights to be used when calculating distances.
If transpose==0, then genes are clustered. If transpose==1, microarrays are
clustered.
The integer npass is the number of times the k-means clustering algorithm
is performed, each time with a different (random) initial condition.
The character method describes how the center of a cluster is found:
method=='a': arithmic mean
method=='m': median
The character dist defines the distance function to be used:
dist=='e': Euclidean distance
dist=='b': City Block distance
dist=='h': Harmonically summed Euclidean distance
dist=='c': correlation
dist=='a': absolute value of the correlation
dist=='u': uncentered correlation
dist=='x': absolute uncentered correlation
dist=='s': Spearman's rank correlation
dist=='k': Kendall's tau
For other values of dist, the default (Euclidean distance) is used.

Return values:
clusterid is an array containing the number of the cluster to which each
  gene/microarray was assigned;
centroids is an array containing the gene expression data for the cluster
  centroids;
error is the within-cluster sum of distances for the optimal k-means
  clustering solution;
nfound is the number of times the optimal solution was found."""

  if rank(data)!=2:
    print "Error in kcluster: data should be a two-dimensional array"
    return
  if dist not in ['e','b','h','c','a','u','x','s','k']:
    print "Error in kcluster: unknown distance function specified (dist='"+dist+"')"
    return
  (n,m) = shape(data)
  if transpose: transpose = 1
  if not mask: mask = ones((n,m))
  if not weight:
    if transpose: weight=ones(n,'d')
    else: weight=ones(m,'d')
  x = cluster.kcluster(nclusters,data,mask,weight,transpose,npass,method,dist)
  return x

def treecluster(data=None,mask=None,weight=None,applyscale=0,transpose=0,dist='e',method='m',distancematrix=None):
  """returns tree, linkdist

This function implements the pairwise single, complete, centroid, and
average linkage hierarchical clustering methods.

The nrows x ncolumns array data contains the gene expression data.
The array mask declares missing data. If mask[i][j]==0, then data[i][j]
is missing.
The array weight contains the weights to be used for the distance
calculation.
If the integer applyscale is nonzero, then the distances in linkdist are
scaled such that all distances are between zero and two (as in case of the
Pearson distance).
The integer transpose defines if rows (genes) or columns (microarrays) are
clustered. If transpose==0, then genes are clustered. If transpose==1,
microarrays are clustered.
The character dist defines the distance function to be used:
dist=='e': Euclidean distance (default)
dist=='b': City Block distance
dist=='h': Harmonically summed Euclidean distance
dist=='c': Pearson correlation
dist=='a': absolute value of the Pearson correlation
dist=='u': uncentered correlation
dist=='x': absolute uncentered correlation
dist=='s': Spearman's rank correlation
dist=='k': Kendall's tau
For other values of dist, the default (Euclidean distance) is used.
The character method specifies which linkage method is used:
method=='s': Single pairwise linkage
method=='m': Complete (maximum) pairwise linkage (default)
method=='c': Centroid linkage
method=='a': Average pairwise linkage
The 2D array distancematrix, which is square and symmetric, is the distance
matrix. Either data or distancematrix should be None. If distancematrix==None,
the hierarchical clustering solution is calculated from the gene expression
data stored in the argument data. If data==None, the hierarchical clustering
solution is calculated from the distance matrix instead. Pairwise centroid-
linkage clustering can be calculated only from the gene expression data and
not from the distance matrix. Pairwise single-, maximum-, and average-linkage
clustering can be calculated from either the gene expression data or from
the distance matrix.

Return values:
tree is an (nobject x 2) array describing the hierarchical clustering
  result. Each row in the array represents one node, with the two columns
  representing the two objects or nodes that are being joined. Objects are
  numbered 0 through (nobjects-1), while nodes are numbered -1 through
  -(nobjects-1).
linkdist is a vector with (nobjects-1) elements containing the distances
between the two subnodes that are joined at each node."""
  if data!=None and distancematrix!=None:
    print "Use either data or distancematrix, do not use both"
    return
  if data!=None:
    if not method in ['c','s','m','a']:
      print "Error in treecluster: keyword method should be 'c', 's', 'm', or 'a'"
      return
    if dist not in ['e','b','h','c','a','u','x','s','k']:
      print "Error in treecluster: unknown distance function specified (dist='"+dist+"')"
      return
    if transpose: transpose = 1
    (n,m) = shape(data)
    if not mask: mask = ones((n,m))
    if not weight:
      if transpose: weight=ones(n,'d')
      else: weight=ones(m,'d')
    return cluster.treecluster(data,mask,weight,applyscale,transpose,dist,method,0)
  if distancematrix!=None:
    if not method in ['s','m','a']:
      print "Error in treecluster: keyword method should be 's', 'm', or 'a'"
      print "(single-, maximum-, or average-linkage clustering). Centroid-"
      print "linkage clustering needs the original gene expression data"
      print "and cannot be performed using the distance matrix alone"
      return
    return cluster.treecluster(distancematrix,mask,weight,applyscale,transpose,dist,method,1)

def somcluster(data,mask=None,weight=None,transpose=0,nxgrid=2,nygrid=1,inittau=0.02,niter=1,dist='e'):
  """returns clusterid, celldata

This function implements a self-organizing map on a rectangular grid.
The nrows x ncolumns array data contains the measurement data
The array mask declares missing data. If mask[i][j]==0, then data[i][j]
is missing.
The array weights contains the weights to be used for the distance
calculation.
The integer transpose defines if rows (genes) or columns (microarrays) are
clustered. If transpose==0, then genes are clustered. If transpose==1,
microarrays are clustered.
The dimensions of the SOM map are nxgrid x nygrid.
The initial value of tau (the neighborbood function) is given by inittau.
The number of iterations is given by niter.
The character dist defines the distance function to be used:
dist=='e': Euclidean distance
dist=='b': City Block distance
dist=='h': Harmonically summed Euclidean distance
dist=='c': correlation
dist=='a': absolute value of the correlation
dist=='u': uncentered correlation
dist=='x': absolute uncentered correlation
dist=='s': Spearman's rank correlation
dist=='k': Kendall's tau
For other values of dist, the default (Euclidean distance) is used.

Return values:
clusterid is an array with two columns, while the number of rows is equal to
  the number of genes or the number of microarrays depending on whether
  genes or microarrays are being clustered. Each row in the array contains
  the x and y coordinates of the cell in the rectangular SOM grid to which
  the gene or microarray was assigned.
celldata is an array with dimensions (nxgrid, nygrid, number of microarrays)
  if genes are being clustered, or (nxgrid, nygrid, number of genes) if
  microarrays are being clustered. Each element [ix][iy] of this array is
  a 1D vector containing the gene expression data for the centroid of the
  cluster in the SOM grid cell with coordinates (ix, iy)."""

  if dist not in ['e','b','h','c','a','u','x','s','k']:
    print "Error in somcluster: unknown distance function specified (dist='"+dist+"')"
    return
  (n,m) = shape(data)
  if transpose: transpose = 1
  if not mask: mask = ones((n,m))
  if not weight:
    if transpose: weight=ones(n,'d')
    else: weight=ones(m,'d')
  x = cluster.somcluster(data,mask,weight,transpose,nxgrid,nygrid,inittau,niter,dist)
  return x

def clustercentroid(data,mask=None,clusterid=None,method='a',transpose=0):
  """The clustercentroid routine calculates the cluster centroids, given to
which cluster each element belongs. The centroid is defined as either the
mean or the median over all elements for each dimension.
The ngenes x nmicroarrays array data contains the gene expression data.
The array mask declares missing data. If mask[i][j]==0, then data[i][j] is
missing.
The integer transpose defines if rows (genes) or columns (microarrays) are
clustered. If transpose==0, then genes are clustered. If transpose==1,
microarrays are clustered.
The array clusterid contains the cluster number for each gene or microarray.
The cluster number should be non-negative.
This function returns an array cdata and an array cmask.
The array cdata contains the cluster centroids. If transpose==0, then the
dimensions of cdata are nclusters x nmicroarrays. If transpose==1, then the
dimensions of cdata are ngenes x nclusters.
The array cmask describes which elements in cdata, if any, are missing."""
  (n,m) = shape(data)
  if not mask: mask = ones((n,m))
  if transpose:
    transpose = 1
    nobjects = m
  else:
    nobjects = n
  if not clusterid: clusterid = zeros(nobjects)
  if min(clusterid)!=0:
    print 'Error in clustercentroid:\nThe cluster numbers in clusterid should be non-negative.'
    return
  nclusters = max(clusterid) + 1
  if method=='a':
    cdata, cmask = cluster.getclustermean(nclusters,data,mask,clusterid,transpose)
  else:
    cdata, cmask = cluster.getclustermedian(nclusters,data,mask,clusterid,transpose)
  return cdata, cmask

def clusterdistance(data,mask=None,weight=None,index1=[0],index2=[0],dist='e',method='a',transpose=0):
  """The distance between two clusters

The array data is a nrows x ncolumns array containing the gene expression
data.
The array mask shows which data are missing. If mask[i][j]==0, then
data[i][j] is missing.
The array weight contains the weights to be used when calculating distances.
The vector index1 identifies which genes/microarrays belong to the first
cluster.
The vector index2 identifies which genes/microarrays belong to the second
cluster.
The character dist defines the distance function to be used:
dist=='e': Euclidean distance
dist=='b': City Block distance
dist=='h': Harmonically summed Euclidean distance
dist=='c': correlation
dist=='a': absolute value of the correlation
dist=='u': uncentered correlation
dist=='x': absolute uncentered correlation
dist=='s': Spearman's rank correlation
dist=='k': Kendall's tau
For other values of dist, the default (Euclidean distance) is used.
The character method specifies how the distance between two clusters is
defined:
method=='a': the distance between the arithmic means of the two clusters
method=='m': the distance between the medians of the two clusters
method=='s': the smallest pairwise distance between members of the two
             clusters
method=='x': the largest pairwise distance between members of the two
             clusters
method=='v': average of the pairwise distances between members of the
             clusters
If transpose==0, then clusters of genes are considered. If transpose==1,
clusters of microarrays are considered."""
  if dist not in ['e','b','h','c','a','u','x','s','k']:
    print "Error in clusterdistance: unknown distance function specified (dist='"+dist+"')"
    return
  (n,m) = shape(data)
  if transpose: transpose = 1
  if not mask: mask = ones((n,m))
  if not weight:
    if transpose: weight=ones(n,'d')
    else: weight=ones(m,'d')
  x = cluster.clusterdistance(data,mask,weight,index1,index2,dist,method,transpose)
  return x


