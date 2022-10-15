# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
"""Cluster Analysis.

The Bio.Cluster provides commonly used clustering algorithms and was
designed with the application to gene expression data in mind. However,
this module can also be used for cluster analysis of other types of data.

Bio.Cluster and the underlying C Clustering Library is described in
M. de Hoon et al. (2004) https://doi.org/10.1093/bioinformatics/bth078
"""

import numbers

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Please install numpy if you want to use Bio.Cluster. "
        "See http://www.numpy.org/"
    ) from None

from . import _cluster

__all__ = (
    "Node",
    "Tree",
    "kcluster",
    "kmedoids",
    "treecluster",
    "somcluster",
    "clusterdistance",
    "clustercentroids",
    "distancematrix",
    "pca",
    "Record",
    "read",
)


__version__ = _cluster.version()


class Node(_cluster.Node):
    """Element of a hierarchical clustering tree.

    A node contains items or other Nodes(sub-nodes).
    """

    __doc__ = _cluster.Node.__doc__


class Tree(_cluster.Tree):
    """Hierarchical clustering tree.

    A Tree consists of Nodes.
    """

    def sort(self, order=None):
        """Sort the hierarchical clustering tree.

        Sort the hierarchical clustering tree by switching the left and
        right subnode of nodes such that the elements in the left-to-right
        order of the tree tend to have increasing order values.

        Return the indices of the elements in the left-to-right order in
        the hierarchical clustering tree, such that the element with index
        indices[i] occurs at position i in the dendrogram.

        """
        n = len(self) + 1
        indices = numpy.ones(n, dtype="intc")
        if order is None:
            order = numpy.ones(n, dtype="d")
        elif isinstance(order, numpy.ndarray):
            order = numpy.require(order, dtype="d", requirements="C")
        else:
            order = numpy.array(order, dtype="d")
        _cluster.Tree.sort(self, indices, order)
        return indices

    def cut(self, nclusters=None):
        """Create clusters by cutting the hierarchical clustering tree.

        Divide the elements in a hierarchical clustering result mytree
        into clusters, and return an array with the number of the cluster
        to which each element was assigned.

        Keyword arguments:
         - nclusters: The desired number of clusters.
        """
        n = len(self) + 1
        indices = numpy.ones(n, dtype="intc")
        if nclusters is None:
            nclusters = n
        _cluster.Tree.cut(self, indices, nclusters)
        return indices


def kcluster(
    data,
    nclusters=2,
    mask=None,
    weight=None,
    transpose=False,
    npass=1,
    method="a",
    dist="e",
    initialid=None,
):
    """Perform k-means clustering.

    This function performs k-means clustering on the values in data, and
    returns the cluster assignments, the within-cluster sum of distances
    of the optimal k-means clustering solution, and the number of times
    the optimal solution was found.

    Keyword arguments:
     - data: nrows x ncolumns array containing the data values.
     - nclusters: number of clusters (the 'k' in k-means).
     - mask: nrows x ncolumns array of integers, showing which data
       are missing. If mask[i,j]==0, then data[i,j] is missing.
     - weight: the weights to be used when calculating distances
     - transpose:
       - if False: rows are clustered;
       - if True: columns are clustered.
     - npass: number of times the k-means clustering algorithm is
       performed, each time with a different (random) initial
       condition.
     - method: specifies how the center of a cluster is found:
       - method == 'a': arithmetic mean;
       - method == 'm': median.
     - dist: specifies the distance function to be used:
       - dist == 'e': Euclidean distance;
       - dist == 'b': City Block distance;
       - dist == 'c': Pearson correlation;
       - dist == 'a': absolute value of the correlation;
       - dist == 'u': uncentered correlation;
       - dist == 'x': absolute uncentered correlation;
       - dist == 's': Spearman's rank correlation;
       - dist == 'k': Kendall's tau.
     - initialid: the initial clustering from which the algorithm
       should start.
       If initialid is None, the routine carries out npass
       repetitions of the EM algorithm, each time starting from a
       different random initial clustering. If initialid is given,
       the routine carries out the EM algorithm only once, starting
       from the given initial clustering and without randomizing the
       order in which items are assigned to clusters (i.e., using
       the same order as in the data matrix). In that case, the
       k-means algorithm is fully deterministic.

    Return values:
     - clusterid: array containing the index of the cluster to which each
       item was assigned in the best k-means clustering solution that was
       found in the npass runs;
     - error: the within-cluster sum of distances for the returned k-means
       clustering solution;
     - nfound: the number of times this solution was found.
    """
    data = __check_data(data)
    shape = data.shape
    if transpose:
        ndata, nitems = shape
    else:
        nitems, ndata = shape
    mask = __check_mask(mask, shape)
    weight = __check_weight(weight, ndata)
    clusterid, npass = __check_initialid(initialid, npass, nitems)
    error, nfound = _cluster.kcluster(
        data, nclusters, mask, weight, transpose, npass, method, dist, clusterid
    )
    return clusterid, error, nfound


def kmedoids(distance, nclusters=2, npass=1, initialid=None):
    """Perform k-medoids clustering.

    This function performs k-medoids clustering, and returns the cluster
    assignments, the within-cluster sum of distances of the optimal
    k-medoids clustering solution, and the number of times the optimal
    solution was found.

    Keyword arguments:
     - distance: The distance matrix between the items. There are three
       ways in which you can pass a distance matrix:
       1. a 2D numpy array (in which only the left-lower part of the array
       will be accessed);
       2. a 1D numpy array containing the distances consecutively;
       3. a list of rows containing the lower-triangular part of
       the distance matrix.

       Examples are:

           >>> from numpy import array
           >>> # option 1:
           >>> distance = array([[0.0, 1.1, 2.3],
           ...                   [1.1, 0.0, 4.5],
           ...                   [2.3, 4.5, 0.0]])
           >>> # option 2:
           >>> distance = array([1.1, 2.3, 4.5])
           >>> # option 3:
           >>> distance = [array([]),
           ...             array([1.1]),
           ...             array([2.3, 4.5])]


       These three correspond to the same distance matrix.
     - nclusters: number of clusters (the 'k' in k-medoids)
     - npass: the number of times the k-medoids clustering algorithm
       is performed, each time with a different (random) initial
       condition.
     - initialid: the initial clustering from which the algorithm should start.
       If initialid is not given, the routine carries out npass
       repetitions of the EM algorithm, each time starting from a
       different random initial clustering. If initialid is given,
       the routine carries out the EM algorithm only once, starting
       from the initial clustering specified by initialid and
       without randomizing the order in which items are assigned to
       clusters (i.e., using the same order as in the data matrix).
       In that case, the k-medoids algorithm is fully deterministic.

    Return values:
     - clusterid: array containing the index of the cluster to which each
       item was assigned in the best k-medoids clustering solution that was
       found in the npass runs; note that the index of a cluster is the index
       of the item that is the medoid of the cluster;
     - error: the within-cluster sum of distances for the returned k-medoids
       clustering solution;
     - nfound: the number of times this solution was found.
    """
    distance = __check_distancematrix(distance)
    nitems = len(distance)
    clusterid, npass = __check_initialid(initialid, npass, nitems)
    error, nfound = _cluster.kmedoids(distance, nclusters, npass, clusterid)
    return clusterid, error, nfound


def treecluster(
    data,
    mask=None,
    weight=None,
    transpose=False,
    method="m",
    dist="e",
    distancematrix=None,
):
    """Perform hierarchical clustering, and return a Tree object.

    This function implements the pairwise single, complete, centroid, and
    average linkage hierarchical clustering methods.

    Keyword arguments:
     - data: nrows x ncolumns array containing the data values.
     - mask: nrows x ncolumns array of integers, showing which data are
       missing. If mask[i][j]==0, then data[i][j] is missing.
     - weight: the weights to be used when calculating distances.
     - transpose:
       - if False, rows are clustered;
       - if True, columns are clustered.
     - dist: specifies the distance function to be used:
       - dist == 'e': Euclidean distance
       - dist == 'b': City Block distance
       - dist == 'c': Pearson correlation
       - dist == 'a': absolute value of the correlation
       - dist == 'u': uncentered correlation
       - dist == 'x': absolute uncentered correlation
       - dist == 's': Spearman's rank correlation
       - dist == 'k': Kendall's tau
     - method: specifies which linkage method is used:
       - method == 's': Single pairwise linkage
       - method == 'm': Complete (maximum) pairwise linkage (default)
       - method == 'c': Centroid linkage
       - method == 'a': Average pairwise linkage
     - distancematrix:  The distance matrix between the items. There are
       three ways in which you can pass a distance matrix:
       1. a 2D numpy array (in which only the left-lower part of the array
       will be accessed);
       2. a 1D numpy array containing the distances consecutively;
       3. a list of rows containing the lower-triangular part of
       the distance matrix.

       Examples are:

           >>> from numpy import array
           >>> # option 1:
           >>> distance = array([[0.0, 1.1, 2.3],
           ...                   [1.1, 0.0, 4.5],
           ...                   [2.3, 4.5, 0.0]])
           >>> # option 2:
           >>> distance = array([1.1, 2.3, 4.5])
           >>> # option 3:
           >>> distance = [array([]),
           ...             array([1.1]),
           ...             array([2.3, 4.5])]

       These three correspond to the same distance matrix.

       PLEASE NOTE:
       As the treecluster routine may shuffle the values in the
       distance matrix as part of the clustering algorithm, be sure
       to save this array in a different variable before calling
       treecluster if you need it later.

    Either data or distancematrix should be None. If distancematrix is None,
    the hierarchical clustering solution is calculated from the values stored
    in the argument data. If data is None, the hierarchical clustering solution
    is instead calculated from the distance matrix. Pairwise centroid-linkage
    clustering can be performed only from the data values and not from the
    distance matrix. Pairwise single-, maximum-, and average-linkage clustering
    can be calculated from the data values or from the distance matrix.

    Return value:
    treecluster returns a Tree object describing the hierarchical clustering
    result. See the description of the Tree class for more information.
    """
    if data is None and distancematrix is None:
        raise ValueError("use either data or distancematrix")
    if data is not None and distancematrix is not None:
        raise ValueError("use either data or distancematrix; do not use both")
    if data is not None:
        data = __check_data(data)
        shape = data.shape
        ndata = shape[0] if transpose else shape[1]
        mask = __check_mask(mask, shape)
        weight = __check_weight(weight, ndata)
    if distancematrix is not None:
        distancematrix = __check_distancematrix(distancematrix)
        if mask is not None:
            raise ValueError("mask is ignored if distancematrix is used")
        if weight is not None:
            raise ValueError("weight is ignored if distancematrix is used")
    tree = Tree()
    _cluster.treecluster(
        tree, data, mask, weight, transpose, method, dist, distancematrix
    )
    return tree


def somcluster(
    data,
    mask=None,
    weight=None,
    transpose=False,
    nxgrid=2,
    nygrid=1,
    inittau=0.02,
    niter=1,
    dist="e",
):
    """Calculate a Self-Organizing Map.

    This function implements a Self-Organizing Map on a rectangular grid.

    Keyword arguments:
     - data: nrows x ncolumns array containing the data values;
     - mask: nrows x ncolumns array of integers, showing which data are
       missing. If mask[i][j]==0, then data[i][j] is missing.
     - weight: the weights to be used when calculating distances
     - transpose:
       - if False: rows are clustered;
       - if True: columns are clustered.
     - nxgrid: the horizontal dimension of the rectangular SOM map
     - nygrid: the vertical dimension of the rectangular SOM map
     - inittau: the initial value of tau (the neighborbood function)
     - niter: the number of iterations
     - dist: specifies the distance function to be used:
       - dist == 'e': Euclidean distance
       - dist == 'b': City Block distance
       - dist == 'c': Pearson correlation
       - dist == 'a': absolute value of the correlation
       - dist == 'u': uncentered correlation
       - dist == 'x': absolute uncentered correlation
       - dist == 's': Spearman's rank correlation
       - dist == 'k': Kendall's tau

    Return values:

     - clusterid: array with two columns, with the number of rows equal to
       the items that are being clustered. Each row in the array contains
       the x and y coordinates of the cell in the rectangular SOM grid to
       which the item was assigned.
     - celldata:  an array with dimensions [nxgrid, nygrid, number of columns]
       if rows are being clustered, or [nxgrid, nygrid, number of rows) if
       columns are being clustered.
       Each element [ix, iy] of this array is a 1D vector containing the
       data values for the centroid of the cluster in the SOM grid cell
       with coordinates [ix, iy].
    """
    if transpose:
        ndata, nitems = data.shape
    else:
        nitems, ndata = data.shape
    data = __check_data(data)
    shape = data.shape
    mask = __check_mask(mask, shape)
    weight = __check_weight(weight, ndata)
    if nxgrid < 1:
        raise ValueError("nxgrid should be a positive integer (default is 2)")
    if nygrid < 1:
        raise ValueError("nygrid should be a positive integer (default is 1)")
    clusterids = numpy.ones((nitems, 2), dtype="intc")
    celldata = numpy.empty((nxgrid, nygrid, ndata), dtype="d")
    _cluster.somcluster(
        clusterids, celldata, data, mask, weight, transpose, inittau, niter, dist
    )
    return clusterids, celldata


def clusterdistance(
    data,
    mask=None,
    weight=None,
    index1=None,
    index2=None,
    method="a",
    dist="e",
    transpose=False,
):
    """Calculate and return the distance between two clusters.

    Keyword arguments:
     - data: nrows x ncolumns array containing the data values.
     - mask: nrows x ncolumns array of integers, showing which data are
       missing. If mask[i, j]==0, then data[i, j] is missing.
     - weight: the weights to be used when calculating distances
     - index1: 1D array identifying which items belong to the
       first cluster. If the cluster contains only one item, then
       index1 can also be written as a single integer.
     - index2: 1D array identifying which items belong to the
       second cluster. If the cluster contains only one item, then
       index2 can also be written as a single integer.
     - dist: specifies the distance function to be used:
       - dist == 'e': Euclidean distance
       - dist == 'b': City Block distance
       - dist == 'c': Pearson correlation
       - dist == 'a': absolute value of the correlation
       - dist == 'u': uncentered correlation
       - dist == 'x': absolute uncentered correlation
       - dist == 's': Spearman's rank correlation
       - dist == 'k': Kendall's tau
     - method: specifies how the distance between two clusters is defined:
       - method == 'a': the distance between the arithmetic means
       of the two clusters
       - method == 'm': the distance between the medians of the two clusters
       - method == 's': the smallest pairwise distance between members
       of the two clusters
       - method == 'x': the largest pairwise distance between members
       of the two clusters
       - method == 'v': average of the pairwise distances between members
       of the two clusters
     - transpose:
       - if False: clusters of rows are considered;
       - if True: clusters of columns are considered.
    """
    data = __check_data(data)
    shape = data.shape
    ndata = shape[0] if transpose else shape[1]
    mask = __check_mask(mask, shape)
    weight = __check_weight(weight, ndata)
    index1 = __check_index(index1)
    index2 = __check_index(index2)
    return _cluster.clusterdistance(
        data, mask, weight, index1, index2, method, dist, transpose
    )


def clustercentroids(data, mask=None, clusterid=None, method="a", transpose=False):
    """Calculate and return the centroid of each cluster.

    The clustercentroids routine calculates the cluster centroids, given to
    which cluster each item belongs. The centroid is defined as either
    the mean or the median over all items for each dimension.

    Keyword arguments:
     - data: nrows x ncolumns array containing the data values.
     - mask: nrows x ncolumns array of integers, showing which data are
       missing. If mask[i, j]==0, then data[i, j] is missing.
     - clusterid: array containing the cluster number for each item.
       The cluster number should be non-negative.
     - method: specifies whether the centroid is calculated from the
       arithmetic mean (method == 'a', default) or the median (method == 'm')
       over each dimension.
     - transpose: if False, each row contains the data for one item;
                  if True, each column contains the data for one item.

    Return values:
     - cdata: 2D array containing the cluster centroids.
       If transpose is False, then the dimensions of cdata are
       nclusters x ncolumns.
       If transpose is True, then the dimensions of cdata are
       nrows x nclusters.
     - cmask: 2D array of integers describing which items in cdata,
       if any, are missing.
    """
    data = __check_data(data)
    mask = __check_mask(mask, data.shape)
    nrows, ncolumns = data.shape
    if clusterid is None:
        n = ncolumns if transpose else nrows
        clusterid = numpy.zeros(n, dtype="intc")
        nclusters = 1
    else:
        clusterid = numpy.require(clusterid, dtype="intc", requirements="C")
        nclusters = max(clusterid + 1)
    if transpose:
        shape = (nrows, nclusters)
    else:
        shape = (nclusters, ncolumns)
    cdata = numpy.zeros(shape, dtype="d")
    cmask = numpy.zeros(shape, dtype="intc")
    _cluster.clustercentroids(data, mask, clusterid, method, transpose, cdata, cmask)
    return cdata, cmask


def distancematrix(data, mask=None, weight=None, transpose=False, dist="e"):
    """Calculate and return a distance matrix from the data.

    This function returns the distance matrix calculated from the data.

    Keyword arguments:
     - data: nrows x ncolumns array containing the data values.
     - mask: nrows x ncolumns array of integers, showing which data are
       missing. If mask[i, j]==0, then data[i, j] is missing.
     - weight: the weights to be used when calculating distances.
     - transpose: if False: the distances between rows are calculated;
                  if True:  the distances between columns are calculated.
     - dist: specifies the distance function to be used:
       - dist == 'e': Euclidean distance
       - dist == 'b': City Block distance
       - dist == 'c': Pearson correlation
       - dist == 'a': absolute value of the correlation
       - dist == 'u': uncentered correlation
       - dist == 'x': absolute uncentered correlation
       - dist == 's': Spearman's rank correlation
       - dist == 'k': Kendall's tau

    Return value:
    The distance matrix is returned as a list of 1D arrays containing the
    distance matrix calculated from the data. The number of columns in eac
    row is equal to the row number. Hence, the first row has zero length.
    For example:

    >>> from numpy import array
    >>> from Bio.Cluster import distancematrix
    >>> data = array([[0, 1,  2,  3],
    ...               [4, 5,  6,  7],
    ...               [8, 9, 10, 11],
    ...               [1, 2,  3,  4]])
    >>> distances = distancematrix(data, dist='e')
    >>> distances
    [array([], dtype=float64), array([ 16.]), array([ 64.,  16.]), array([  1.,   9.,  49.])]

    which can be rewritten as
       distances = [array([], dtype=float64),
                    array([ 16.]),
                    array([ 64.,  16.]),
                    array([  1.,   9.,  49.])]

    This corresponds to the distance matrix:

        [ 0., 16., 64.,  1.]
        [16.,  0., 16.,  9.]
        [64., 16.,  0., 49.]
        [ 1.,  9., 49.,  0.]
    """
    data = __check_data(data)
    shape = data.shape
    mask = __check_mask(mask, shape)
    if transpose:
        ndata, nitems = shape
    else:
        nitems, ndata = shape
    weight = __check_weight(weight, ndata)
    matrix = [numpy.empty(i, dtype="d") for i in range(nitems)]
    _cluster.distancematrix(data, mask, weight, transpose, dist, matrix)
    return matrix


def pca(data):
    """Perform principal component analysis.

    Keyword arguments:
     - data: nrows x ncolumns array containing the data values.

    Return value:
    This function returns an array containing the mean of each column, the
    principal components as an nmin x ncolumns array, as well as the
    coordinates (an nrows x nmin array) of the data along the principal
    components, and the associated eigenvalues. The principal components, the
    coordinates, and the eigenvalues are sorted by the magnitude of the
    eigenvalue, with the largest eigenvalues appearing first. Here, nmin is
    the smaller of nrows and ncolumns.
    Adding the column means to the dot product of the coordinates and the
    principal components recreates the data matrix:

    >>> from numpy import array, dot, amax, amin
    >>> from Bio.Cluster import pca
    >>> matrix = array([[ 0.,  0.,  0.],
    ...                 [ 1.,  0.,  0.],
    ...                 [ 7.,  3.,  0.],
    ...                 [ 4.,  2.,  6.]])
    >>> columnmean, coordinates, pc, _ = pca(matrix)
    >>> m = matrix - (columnmean + dot(coordinates, pc))
    >>> amax(m) < 1e-12 and amin(m) > -1e-12
    True

    """
    data = __check_data(data)
    nrows, ncols = data.shape
    nmin = min(nrows, ncols)
    columnmean = numpy.empty(ncols, dtype="d")
    pc = numpy.empty((nmin, ncols), dtype="d")
    coordinates = numpy.empty((nrows, nmin), dtype="d")
    eigenvalues = numpy.empty(nmin, dtype="d")
    _cluster.pca(data, columnmean, coordinates, pc, eigenvalues)
    return columnmean, coordinates, pc, eigenvalues


class Record:
    """Store gene expression data.

    A Record stores the gene expression data and related information contained
    in a data file following the file format defined for Michael Eisen's
    Cluster/TreeView program.

    Attributes:
     - data: a matrix containing the gene expression data
     - mask: a matrix containing only 1's and 0's, denoting which values
       are present (1) or missing (0). If all items of mask are
       one (no missing data), then mask is set to None.
     - geneid: a list containing a unique identifier for each gene
       (e.g., ORF name)
     - genename: a list containing an additional description for each gene
       (e.g., gene name)
     - gweight: the weight to be used for each gene when calculating the
       distance
     - gorder: an array of real numbers indicating the preferred order of the
       genes in the output file
     - expid: a list containing a unique identifier for each sample.
     - eweight: the weight to be used for each sample when calculating the
       distance
     - eorder: an array of real numbers indication the preferred order of the
       samples in the output file
     - uniqid: the string that was used instead of UNIQID in the input file.

    """

    def __init__(self, handle=None):
        """Read gene expression data from the file handle and return a Record.

        The file should be in the format defined for Michael Eisen's
        Cluster/TreeView program.
        """
        self.data = None
        self.mask = None
        self.geneid = None
        self.genename = None
        self.gweight = None
        self.gorder = None
        self.expid = None
        self.eweight = None
        self.eorder = None
        self.uniqid = None
        if not handle:
            return
        line = handle.readline().strip("\r\n").split("\t")
        n = len(line)
        self.uniqid = line[0]
        self.expid = []
        cols = {0: "GENEID"}
        for word in line[1:]:
            if word == "NAME":
                cols[line.index(word)] = word
                self.genename = []
            elif word == "GWEIGHT":
                cols[line.index(word)] = word
                self.gweight = []
            elif word == "GORDER":
                cols[line.index(word)] = word
                self.gorder = []
            else:
                self.expid.append(word)
        self.geneid = []
        self.data = []
        self.mask = []
        needmask = 0
        for line in handle:
            line = line.strip("\r\n").split("\t")
            if len(line) != n:
                raise ValueError(
                    "Line with %d columns found (expected %d)" % (len(line), n)
                )
            if line[0] == "EWEIGHT":
                i = max(cols) + 1
                self.eweight = numpy.array(line[i:], float)
                continue
            if line[0] == "EORDER":
                i = max(cols) + 1
                self.eorder = numpy.array(line[i:], float)
                continue
            rowdata = []
            rowmask = []
            n = len(line)
            for i in range(n):
                word = line[i]
                if i in cols:
                    if cols[i] == "GENEID":
                        self.geneid.append(word)
                    if cols[i] == "NAME":
                        self.genename.append(word)
                    if cols[i] == "GWEIGHT":
                        self.gweight.append(float(word))
                    if cols[i] == "GORDER":
                        self.gorder.append(float(word))
                    continue
                if not word:
                    rowdata.append(0.0)
                    rowmask.append(0)
                    needmask = 1
                else:
                    rowdata.append(float(word))
                    rowmask.append(1)
            self.data.append(rowdata)
            self.mask.append(rowmask)
        self.data = numpy.array(self.data)
        if needmask:
            self.mask = numpy.array(self.mask, int)
        else:
            self.mask = None
        if self.gweight:
            self.gweight = numpy.array(self.gweight)
        if self.gorder:
            self.gorder = numpy.array(self.gorder)

    def treecluster(self, transpose=False, method="m", dist="e"):
        """Apply hierarchical clustering and return a Tree object.

        The pairwise single, complete, centroid, and average linkage
        hierarchical clustering methods are available.

        Keyword arguments:
         - transpose: if False: rows are clustered;
                      if True: columns are clustered.
         - dist: specifies the distance function to be used:
           - dist == 'e': Euclidean distance
           - dist == 'b': City Block distance
           - dist == 'c': Pearson correlation
           - dist == 'a': absolute value of the correlation
           - dist == 'u': uncentered correlation
           - dist == 'x': absolute uncentered correlation
           - dist == 's': Spearman's rank correlation
           - dist == 'k': Kendall's tau
         - method: specifies which linkage method is used:
           - method == 's': Single pairwise linkage
           - method == 'm': Complete (maximum) pairwise linkage (default)
           - method == 'c': Centroid linkage
           - method == 'a': Average pairwise linkage

        See the description of the Tree class for more information about
        the Tree object returned by this method.
        """
        if transpose:
            weight = self.gweight
        else:
            weight = self.eweight
        return treecluster(self.data, self.mask, weight, transpose, method, dist)

    def kcluster(
        self,
        nclusters=2,
        transpose=False,
        npass=1,
        method="a",
        dist="e",
        initialid=None,
    ):
        """Apply k-means or k-median clustering.

        This method returns a tuple (clusterid, error, nfound).

        Keyword arguments:
         - nclusters: number of clusters (the 'k' in k-means)
         - transpose: if False, genes (rows) are clustered;
                      if True, samples (columns) are clustered.
         - npass: number of times the k-means clustering algorithm is
           performed, each time with a different (random) initial condition.
         - method: specifies how the center of a cluster is found:
           - method == 'a': arithmetic mean
           - method == 'm': median
         - dist: specifies the distance function to be used:
           - dist == 'e': Euclidean distance
           - dist == 'b': City Block distance
           - dist == 'c': Pearson correlation
           - dist == 'a': absolute value of the correlation
           - dist == 'u': uncentered correlation
           - dist == 'x': absolute uncentered correlation
           - dist == 's': Spearman's rank correlation
           - dist == 'k': Kendall's tau
         - initialid: the initial clustering from which the algorithm should
           start. If initialid is None, the routine carries out npass
           repetitions of the EM algorithm, each time starting from a different
           random initial clustering. If initialid is given, the routine
           carries out the EM algorithm only once, starting from the given
           initial clustering and without randomizing the order in which items
           are assigned to clusters (i.e., using the same order as in the data
           matrix). In that case, the k-means algorithm is fully deterministic.

        Return values:
         - clusterid: array containing the number of the cluster to which each
           gene/sample was assigned in the best k-means clustering
           solution that was found in the npass runs;
         - error: the within-cluster sum of distances for the returned
           k-means clustering solution;
         - nfound: the number of times this solution was found.
        """
        if transpose:
            weight = self.gweight
        else:
            weight = self.eweight
        return kcluster(
            self.data,
            nclusters,
            self.mask,
            weight,
            transpose,
            npass,
            method,
            dist,
            initialid,
        )

    def somcluster(
        self, transpose=False, nxgrid=2, nygrid=1, inittau=0.02, niter=1, dist="e"
    ):
        """Calculate a self-organizing map on a rectangular grid.

        The somcluster method returns a tuple (clusterid, celldata).

        Keyword arguments:
         - transpose: if False, genes (rows) are clustered;
                      if True,  samples (columns) are clustered.
         - nxgrid: the horizontal dimension of the rectangular SOM map
         - nygrid: the vertical dimension of the rectangular SOM map
         - inittau: the initial value of tau (the neighborbood function)
         - niter: the number of iterations
         - dist: specifies the distance function to be used:
           - dist == 'e': Euclidean distance
           - dist == 'b': City Block distance
           - dist == 'c': Pearson correlation
           - dist == 'a': absolute value of the correlation
           - dist == 'u': uncentered correlation
           - dist == 'x': absolute uncentered correlation
           - dist == 's': Spearman's rank correlation
           - dist == 'k': Kendall's tau

        Return values:
         - clusterid: array with two columns, while the number of rows is equal
           to the number of genes or the number of samples depending on
           whether genes or samples are being clustered. Each row in
           the array contains the x and y coordinates of the cell in the
           rectangular SOM grid to which the gene or samples was assigned.
         - celldata: an array with dimensions (nxgrid, nygrid, number of
           samples) if genes are being clustered, or (nxgrid, nygrid,
           number of genes) if samples are being clustered. Each item
           [ix, iy] of this array is a 1D vector containing the gene
           expression data for the centroid of the cluster in the SOM grid
           cell with coordinates [ix, iy].
        """
        if transpose:
            weight = self.gweight
        else:
            weight = self.eweight
        return somcluster(
            self.data,
            self.mask,
            weight,
            transpose,
            nxgrid,
            nygrid,
            inittau,
            niter,
            dist,
        )

    def clustercentroids(self, clusterid=None, method="a", transpose=False):
        """Calculate the cluster centroids and return a tuple (cdata, cmask).

        The centroid is defined as either the mean or the median over all
        items for each dimension.

        Keyword arguments:
         - data: nrows x ncolumns array containing the expression data
         - mask: nrows x ncolumns array of integers, showing which data
           are missing. If mask[i, j]==0, then data[i, j] is missing.
         - transpose: if False, gene (row) clusters are considered;
                      if True, sample (column) clusters are considered.
         - clusterid: array containing the cluster number for each gene or
           sample. The cluster number should be non-negative.
         - method: specifies how the centroid is calculated:
           - method == 'a': arithmetic mean over each dimension. (default)
           - method == 'm': median over each dimension.

        Return values:
         - cdata: 2D array containing the cluster centroids. If transpose
           is False, then the dimensions of cdata are nclusters x ncolumns.
           If transpose is True, then the dimensions of cdata are nrows x
           nclusters.
         - cmask: 2D array of integers describing which items in cdata,
           if any, are missing.
        """
        return clustercentroids(self.data, self.mask, clusterid, method, transpose)

    def clusterdistance(
        self, index1=0, index2=0, method="a", dist="e", transpose=False
    ):
        """Calculate the distance between two clusters.

        Keyword arguments:
         - index1: 1D array identifying which genes/samples belong to the
           first cluster. If the cluster contains only one gene, then
           index1 can also be written as a single integer.
         - index2: 1D array identifying which genes/samples belong to the
           second cluster. If the cluster contains only one gene, then
           index2 can also be written as a single integer.
         - transpose: if False, genes (rows) are clustered;
                      if True, samples (columns) are clustered.
         - dist: specifies the distance function to be used:
           - dist == 'e': Euclidean distance
           - dist == 'b': City Block distance
           - dist == 'c': Pearson correlation
           - dist == 'a': absolute value of the correlation
           - dist == 'u': uncentered correlation
           - dist == 'x': absolute uncentered correlation
           - dist == 's': Spearman's rank correlation
           - dist == 'k': Kendall's tau
         - method: specifies how the distance between two clusters is defined:
           - method == 'a': the distance between the arithmetic means
           of the two clusters
           - method == 'm': the distance between the medians of the
           two clusters
           - method == 's': the smallest pairwise distance between members
           of the two clusters
           - method == 'x': the largest pairwise distance between members
           of the two clusters
           - method == 'v': average of the pairwise distances between members
           of the two clusters
         - transpose: if False: clusters of rows are considered;
                      if True: clusters of columns are considered.
        """
        if transpose:
            weight = self.gweight
        else:
            weight = self.eweight
        return clusterdistance(
            self.data, self.mask, weight, index1, index2, method, dist, transpose
        )

    def distancematrix(self, transpose=False, dist="e"):
        """Calculate the distance matrix and return it as a list of arrays.

        Keyword arguments:
         - transpose:
             if False: calculate the distances between genes (rows);
             if True: calculate the distances between samples (columns).
         - dist: specifies the distance function to be used:
           - dist == 'e': Euclidean distance
           - dist == 'b': City Block distance
           - dist == 'c': Pearson correlation
           - dist == 'a': absolute value of the correlation
           - dist == 'u': uncentered correlation
           - dist == 'x': absolute uncentered correlation
           - dist == 's': Spearman's rank correlation
           - dist == 'k': Kendall's tau

        Return value:

        The distance matrix is returned as a list of 1D arrays containing the
        distance matrix between the gene expression data. The number of columns
        in each row is equal to the row number. Hence, the first row has zero
        length. An example of the return value is:

            matrix = [[],
                      array([1.]),
                      array([7., 3.]),
                      array([4., 2., 6.])]

        This corresponds to the distance matrix:

            [0., 1., 7., 4.]
            [1., 0., 3., 2.]
            [7., 3., 0., 6.]
            [4., 2., 6., 0.]

        """
        if transpose:
            weight = self.gweight
        else:
            weight = self.eweight
        return distancematrix(self.data, self.mask, weight, transpose, dist)

    def save(self, jobname, geneclusters=None, expclusters=None):
        """Save the clustering results.

        The saved files follow the convention for the Java TreeView program,
        which can therefore be used to view the clustering result.

        Keyword arguments:
         - jobname: The base name of the files to be saved. The filenames
           are jobname.cdt, jobname.gtr, and jobname.atr for hierarchical
           clustering, and jobname-K*.cdt, jobname-K*.kgg, jobname-K*.kag
           for k-means clustering results.
         - geneclusters: For hierarchical clustering results, geneclusters
           is a Tree object as returned by the treecluster method. For k-means
           clustering results, geneclusters is a vector containing ngenes
           integers, describing to which cluster a given gene belongs. This
           vector can be calculated by kcluster.
         - expclusters: For hierarchical clustering results, expclusters
           is a Tree object as returned by the treecluster method. For k-means
           clustering results, expclusters is a vector containing nexps
           integers, describing to which cluster a given sample belongs. This
           vector can be calculated by kcluster.
        """
        (ngenes, nexps) = numpy.shape(self.data)
        if self.gorder is None:
            gorder = numpy.arange(ngenes)
        else:
            gorder = self.gorder
        if self.eorder is None:
            eorder = numpy.arange(nexps)
        else:
            eorder = self.eorder
        if (
            geneclusters is not None
            and expclusters is not None
            and type(geneclusters) != type(expclusters)
        ):
            raise ValueError(
                "found one k-means and one hierarchical "
                "clustering solution in geneclusters and "
                "expclusters"
            )
        gid = 0
        aid = 0
        filename = jobname
        postfix = ""
        if isinstance(geneclusters, Tree):
            # This is a hierarchical clustering result.
            geneindex = self._savetree(jobname, geneclusters, gorder, False)
            gid = 1
        elif geneclusters is not None:
            # This is a k-means clustering result.
            filename = jobname + "_K"
            k = max(geneclusters) + 1
            kggfilename = "%s_K_G%d.kgg" % (jobname, k)
            geneindex = self._savekmeans(kggfilename, geneclusters, gorder, False)
            postfix = "_G%d" % k
        else:
            geneindex = numpy.argsort(gorder)
        if isinstance(expclusters, Tree):
            # This is a hierarchical clustering result.
            expindex = self._savetree(jobname, expclusters, eorder, True)
            aid = 1
        elif expclusters is not None:
            # This is a k-means clustering result.
            filename = jobname + "_K"
            k = max(expclusters) + 1
            kagfilename = "%s_K_A%d.kag" % (jobname, k)
            expindex = self._savekmeans(kagfilename, expclusters, eorder, True)
            postfix += "_A%d" % k
        else:
            expindex = numpy.argsort(eorder)
        filename = filename + postfix
        self._savedata(filename, gid, aid, geneindex, expindex)

    def _savetree(self, jobname, tree, order, transpose):
        """Save the hierarchical clustering solution (PRIVATE)."""
        if transpose:
            extension = ".atr"
            keyword = "ARRY"
        else:
            extension = ".gtr"
            keyword = "GENE"
        index = tree.sort(order)
        nnodes = len(tree)
        with open(jobname + extension, "w") as outputfile:
            nodeID = [""] * nnodes
            nodedist = numpy.array([node.distance for node in tree[:]])
            for nodeindex in range(nnodes):
                min1 = tree[nodeindex].left
                min2 = tree[nodeindex].right
                nodeID[nodeindex] = "NODE%dX" % (nodeindex + 1)
                outputfile.write(nodeID[nodeindex])
                outputfile.write("\t")
                if min1 < 0:
                    index1 = -min1 - 1
                    outputfile.write(nodeID[index1] + "\t")
                    nodedist[nodeindex] = max(nodedist[nodeindex], nodedist[index1])
                else:
                    outputfile.write("%s%dX\t" % (keyword, min1))
                if min2 < 0:
                    index2 = -min2 - 1
                    outputfile.write(nodeID[index2] + "\t")
                    nodedist[nodeindex] = max(nodedist[nodeindex], nodedist[index2])
                else:
                    outputfile.write("%s%dX\t" % (keyword, min2))
                outputfile.write(str(1.0 - nodedist[nodeindex]))
                outputfile.write("\n")
        return index

    def _savekmeans(self, filename, clusterids, order, transpose):
        """Save the k-means clustering solution (PRIVATE)."""
        if transpose:
            label = "ARRAY"
            names = self.expid
        else:
            label = self.uniqid
            names = self.geneid
        with open(filename, "w") as outputfile:
            outputfile.write(label + "\tGROUP\n")
            index = numpy.argsort(order)
            n = len(names)
            sortedindex = numpy.zeros(n, int)
            counter = 0
            cluster = 0
            while counter < n:
                for j in index:
                    if clusterids[j] == cluster:
                        outputfile.write(f"{names[j]}\t{cluster}\n")
                        sortedindex[counter] = j
                        counter += 1
                cluster += 1
        return sortedindex

    def _savedata(self, jobname, gid, aid, geneindex, expindex):
        """Save the clustered data (PRIVATE)."""
        if self.genename is None:
            genename = self.geneid
        else:
            genename = self.genename
        (ngenes, nexps) = numpy.shape(self.data)
        with open(jobname + ".cdt", "w") as outputfile:
            if self.mask is not None:
                mask = self.mask
            else:
                mask = numpy.ones((ngenes, nexps), int)
            if self.gweight is not None:
                gweight = self.gweight
            else:
                gweight = numpy.ones(ngenes)
            if self.eweight is not None:
                eweight = self.eweight
            else:
                eweight = numpy.ones(nexps)
            if gid:
                outputfile.write("GID\t")
            outputfile.write(self.uniqid)
            outputfile.write("\tNAME\tGWEIGHT")
            # Now add headers for data columns.
            for j in expindex:
                outputfile.write(f"\t{self.expid[j]}")
            outputfile.write("\n")
            if aid:
                outputfile.write("AID")
                if gid:
                    outputfile.write("\t")
                outputfile.write("\t\t")
                for j in expindex:
                    outputfile.write("\tARRY%dX" % j)
                outputfile.write("\n")
            outputfile.write("EWEIGHT")
            if gid:
                outputfile.write("\t")
            outputfile.write("\t\t")
            for j in expindex:
                outputfile.write(f"\t{eweight[j]:f}")
            outputfile.write("\n")
            for i in geneindex:
                if gid:
                    outputfile.write("GENE%dX\t" % i)
                outputfile.write(f"{self.geneid[i]}\t{genename[i]}\t{gweight[i]:f}")
                for j in expindex:
                    outputfile.write("\t")
                    if mask[i, j]:
                        outputfile.write(str(self.data[i, j]))
                outputfile.write("\n")


def read(handle):
    """Read gene expression data from the file handle and return a Record.

    The file should be in the file format defined for Michael Eisen's
    Cluster/TreeView program.
    """
    return Record(handle)


# Everything below is private
#


def __check_data(data):
    if isinstance(data, numpy.ndarray):
        data = numpy.require(data, dtype="d", requirements="C")
    else:
        data = numpy.array(data, dtype="d")
    if data.ndim != 2:
        raise ValueError("data should be 2-dimensional")
    if numpy.isnan(data).any():
        raise ValueError("data contains NaN values")
    return data


def __check_mask(mask, shape):
    if mask is None:
        return numpy.ones(shape, dtype="intc")
    elif isinstance(mask, numpy.ndarray):
        return numpy.require(mask, dtype="intc", requirements="C")
    else:
        return numpy.array(mask, dtype="intc")


def __check_weight(weight, ndata):
    if weight is None:
        return numpy.ones(ndata, dtype="d")
    if isinstance(weight, numpy.ndarray):
        weight = numpy.require(weight, dtype="d", requirements="C")
    else:
        weight = numpy.array(weight, dtype="d")
    if numpy.isnan(weight).any():
        raise ValueError("weight contains NaN values")
    return weight


def __check_initialid(initialid, npass, nitems):
    if initialid is None:
        if npass <= 0:
            raise ValueError("npass should be a positive integer")
        clusterid = numpy.empty(nitems, dtype="intc")
    else:
        npass = 0
        clusterid = numpy.array(initialid, dtype="intc")
    return clusterid, npass


def __check_index(index):
    if index is None:
        return numpy.zeros(1, dtype="intc")
    elif isinstance(index, numbers.Integral):
        return numpy.array([index], dtype="intc")
    elif isinstance(index, numpy.ndarray):
        return numpy.require(index, dtype="intc", requirements="C")
    else:
        return numpy.array(index, dtype="intc")


def __check_distancematrix(distancematrix):
    if distancematrix is None:
        return distancematrix
    if isinstance(distancematrix, numpy.ndarray):
        distancematrix = numpy.require(distancematrix, dtype="d", requirements="C")
    else:
        try:
            distancematrix = numpy.array(distancematrix, dtype="d")
        except ValueError:
            n = len(distancematrix)
            d = [None] * n
            for i, row in enumerate(distancematrix):
                if isinstance(row, numpy.ndarray):
                    row = numpy.require(row, dtype="d", requirements="C")
                else:
                    row = numpy.array(row, dtype="d")
                if row.ndim != 1:
                    raise ValueError("row %d is not one-dimensional" % i) from None
                m = len(row)
                if m != i:
                    raise ValueError(
                        "row %d has incorrect size (%d, expected %d)" % (i, m, i)
                    ) from None
                if numpy.isnan(row).any():
                    raise ValueError("distancematrix contains NaN values") from None
                d[i] = row
            return d
    if numpy.isnan(distancematrix).any():
        raise ValueError("distancematrix contains NaN values")
    return distancematrix
