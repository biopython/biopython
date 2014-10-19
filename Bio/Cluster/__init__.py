# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#

import numpy

from Bio.Cluster.cluster import *

__docformat__ = "restructuredtext en"

def _treesort(order, nodeorder, nodecounts, tree):
    # Find the order of the nodes consistent with the hierarchical clustering
    # tree, taking into account the preferred order of nodes.
    nNodes = len(tree)
    nElements = nNodes + 1
    neworder = numpy.zeros(nElements)
    clusterids = numpy.arange(nElements)
    for i in range(nNodes):
        i1 = tree[i].left
        i2 = tree[i].right
        if i1 < 0:
            order1 = nodeorder[-i1-1]
            count1 = nodecounts[-i1-1]
        else:
            order1 = order[i1]
            count1 = 1
        if i2 < 0:
            order2 = nodeorder[-i2-1]
            count2 = nodecounts[-i2-1]
        else:
            order2 = order[i2]
            count2 = 1
        # If order1 and order2 are equal, their order is determined
        # by the order in which they were clustered
        if i1 < i2:
            if order1 < order2:
                increase = count1
            else:
                increase = count2
            for j in range(nElements):
                clusterid = clusterids[j]
                if clusterid == i1 and order1 >= order2:
                    neworder[j] += increase
                if clusterid == i2 and order1 < order2:
                    neworder[j] += increase
                if clusterid == i1 or clusterid == i2:
                    clusterids[j] = -i-1
        else:
            if order1 <= order2:
                increase = count1
            else:
                increase = count2
            for j in range(nElements):
                clusterid = clusterids[j]
                if clusterid == i1 and order1 > order2:
                    neworder[j] += increase
                if clusterid == i2 and order1 <= order2:
                    neworder[j] += increase
                if clusterid == i1 or clusterid == i2:
                    clusterids[j] = -i-1
    return numpy.argsort(neworder)


def _savetree(jobname, tree, order, transpose):
    # Save the hierarchical clustering solution given by the tree, following
    # the specified order, in a file whose name is based on jobname.
    if transpose == 0:
        extension = ".gtr"
        keyword = "GENE"
    else:
        extension = ".atr"
        keyword = "ARRY"
    nnodes = len(tree)
    with open(jobname+extension, "w") as outputfile:
        nodeindex = 0
        nodeID = [''] * nnodes
        nodecounts = numpy.zeros(nnodes, int)
        nodeorder = numpy.zeros(nnodes)
        nodedist = numpy.array([node.distance for node in tree])
        for nodeindex in range(nnodes):
            min1 = tree[nodeindex].left
            min2 = tree[nodeindex].right
            nodeID[nodeindex] = "NODE%dX" % (nodeindex+1)
            outputfile.write(nodeID[nodeindex])
            outputfile.write("\t")
            if min1 < 0:
                index1 = -min1-1
                order1 = nodeorder[index1]
                counts1 = nodecounts[index1]
                outputfile.write(nodeID[index1]+"\t")
                nodedist[nodeindex] = max(nodedist[nodeindex], nodedist[index1])
            else:
                order1 = order[min1]
                counts1 = 1
                outputfile.write("%s%dX\t" % (keyword, min1))
            if min2 < 0:
                index2 = -min2-1
                order2 = nodeorder[index2]
                counts2 = nodecounts[index2]
                outputfile.write(nodeID[index2]+"\t")
                nodedist[nodeindex] = max(nodedist[nodeindex], nodedist[index2])
            else:
                order2 = order[min2]
                counts2 = 1
                outputfile.write("%s%dX\t" % (keyword, min2))
            outputfile.write(str(1.0-nodedist[nodeindex]))
            outputfile.write("\n")
            counts = counts1 + counts2
            nodecounts[nodeindex] = counts
            nodeorder[nodeindex] = (counts1*order1+counts2*order2) / counts
    # Now set up order based on the tree structure
    index = _treesort(order, nodeorder, nodecounts, tree)
    return index


class Record(object):
    """Store gene expression data.

A Record stores the gene expression data and related information contained
in a data file following the file format defined for Michael Eisen's
Cluster/TreeView program. A Record has the following members:

  - data:     a matrix containing the gene expression data
  - mask:     a matrix containing only 1's and 0's, denoting which values
    are present (1) or missing (0). If all elements of mask are
    one (no missing data), then mask is set to None.
  - geneid:   a list containing a unique identifier for each gene
    (e.g., ORF name)
  - genename: a list containing an additional description for each gene
    (e.g., gene name)
  - gweight:  the weight to be used for each gene when calculating the
    distance
  - gorder:   an array of real numbers indicating the preferred order of the
    genes in the output file
  - expid:    a list containing a unique identifier for each experimental
    condition
  - eweight:  the weight to be used for each experimental condition when
    calculating the distance
  - eorder:   an array of real numbers indication the preferred order in the
    output file of the experimental conditions
  - uniqid:   the string that was used instead of UNIQID in the input file.

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
            elif word=="GORDER":
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
                raise ValueError("Line with %d columns found (expected %d)" %
                                 (len(line), n))
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

    def treecluster(self, transpose=0, method='m', dist='e'):
        """Apply hierarchical clustering and return a Tree object.

The pairwise single, complete, centroid, and average linkage hierarchical
clustering methods are available.

- transpose: if equal to 0, genes (rows) are clustered;
  if equal to 1, microarrays (columns) are clustered.
- dist     : specifies the distance function to be used:

  - dist=='e': Euclidean distance
  - dist=='b': City Block distance
  - dist=='c': Pearson correlation
  - dist=='a': absolute value of the correlation
  - dist=='u': uncentered correlation
  - dist=='x': absolute uncentered correlation
  - dist=='s': Spearman's rank correlation
  - dist=='k': Kendall's tau

- method   : specifies which linkage method is used:

  - method=='s': Single pairwise linkage
  - method=='m': Complete (maximum) pairwise linkage (default)
  - method=='c': Centroid linkage
  - method=='a': Average pairwise linkage

See the description of the Tree class for more information about the Tree
object returned by this method.

"""
        if transpose == 0:
            weight = self.eweight
        else:
            weight = self.gweight
        return treecluster(self.data, self.mask, weight, transpose, method,
                           dist)

    def kcluster(self, nclusters=2, transpose=0, npass=1, method='a', dist='e',
                 initialid=None):
        """Apply k-means or k-median clustering.

This method returns a tuple (clusterid, error, nfound).

  - nclusters: number of clusters (the 'k' in k-means)
  - transpose: if equal to 0, genes (rows) are clustered;
    if equal to 1, microarrays (columns) are clustered.
  - npass    : number of times the k-means clustering algorithm is
    performed, each time with a different (random) initial
    condition.
  - method   : specifies how the center of a cluster is found:
    method=='a': arithmetic mean
    method=='m': median
  - dist     : specifies the distance function to be used:

    - dist=='e': Euclidean distance
    - dist=='b': City Block distance
    - dist=='c': Pearson correlation
    - dist=='a': absolute value of the correlation
    - dist=='u': uncentered correlation
    - dist=='x': absolute uncentered correlation
    - dist=='s': Spearman's rank correlation
    - dist=='k': Kendall's tau

  - initialid: the initial clustering from which the algorithm should start.
    If initialid is None, the routine carries out npass
    repetitions of the EM algorithm, each time starting from a
    different random initial clustering. If initialid is given,
    the routine carries out the EM algorithm only once, starting
    from the given initial clustering and without randomizing the
    order in which items are assigned to clusters (i.e., using
    the same order as in the data matrix). In that case, the
    k-means algorithm is fully deterministic.

Return values:
  - clusterid: array containing the number of the cluster to which each
    gene/microarray was assigned in the best k-means clustering
    solution that was found in the npass runs;
  - error:     the within-cluster sum of distances for the returned k-means
    clustering solution;
  - nfound:    the number of times this solution was found.

"""

        if transpose == 0:
            weight = self.eweight
        else:
            weight = self.gweight
        return kcluster(self.data, nclusters, self.mask, weight, transpose,
                        npass, method, dist, initialid)

    def somcluster(self, transpose=0, nxgrid=2, nygrid=1, inittau=0.02,
                   niter=1, dist='e'):
        """Calculate a self-organizing map on a rectangular grid.

The somcluster method returns a tuple (clusterid, celldata).

  - transpose: if equal to 0, genes (rows) are clustered;
    if equal to 1, microarrays (columns) are clustered.
  - nxgrid   : the horizontal dimension of the rectangular SOM map
  - nygrid   : the vertical dimension of the rectangular SOM map
  - inittau  : the initial value of tau (the neighborbood function)
  - niter    : the number of iterations
  - dist     : specifies the distance function to be used:

    - dist=='e': Euclidean distance
    - dist=='b': City Block distance
    - dist=='c': Pearson correlation
    - dist=='a': absolute value of the correlation
    - dist=='u': uncentered correlation
    - dist=='x': absolute uncentered correlation
    - dist=='s': Spearman's rank correlation
    - dist=='k': Kendall's tau

Return values:

    - clusterid: array with two columns, while the number of rows is equal to
      the number of genes or the number of microarrays depending on
      whether genes or microarrays are being clustered. Each row in
      the array contains the x and y coordinates of the cell in the
      rectangular SOM grid to which the gene or microarray was
      assigned.
    - celldata:  an array with dimensions (nxgrid, nygrid, number of
      microarrays) if genes are being clustered, or (nxgrid,
      nygrid, number of genes) if microarrays are being clustered.
      Each element [ix][iy] of this array is a 1D vector containing
      the gene expression data for the centroid of the cluster in
      the SOM grid cell with coordinates (ix, iy).

"""

        if transpose == 0:
            weight = self.eweight
        else:
            weight = self.gweight
        return somcluster(self.data, self.mask, weight, transpose,
                          nxgrid, nygrid, inittau, niter, dist)

    def clustercentroids(self, clusterid=None, method='a', transpose=0):
        """Calculate the cluster centroids and return a tuple (cdata, cmask).

The centroid is defined as either the mean or the median over all elements
for each dimension.

  - data     : nrows x ncolumns array containing the expression data
  - mask     : nrows x ncolumns array of integers, showing which data are
    missing. If mask[i][j]==0, then data[i][j] is missing.
  - transpose: if equal to 0, gene (row) clusters are considered;
    if equal to 1, microarray (column) clusters are considered.
  - clusterid: array containing the cluster number for each gene or
    microarray. The cluster number should be non-negative.
  - method   : specifies how the centroid is calculated:
    method=='a': arithmetic mean over each dimension. (default)
    method=='m': median over each dimension.

Return values:
  - cdata    : 2D array containing the cluster centroids. If transpose==0,
    then the dimensions of cdata are nclusters x ncolumns. If
    transpose==1, then the dimensions of cdata are
    nrows x nclusters.
  - cmask    : 2D array of integers describing which elements in cdata,
    if any, are missing.

"""
        return clustercentroids(self.data, self.mask, clusterid, method,
                                transpose)

    def clusterdistance(self, index1=[0], index2=[0], method='a', dist='e',
                        transpose=0):
        """Calculate the distance between two clusters.

  - index1   : 1D array identifying which genes/microarrays belong to the
    first cluster. If the cluster contains only one gene, then
    index1 can also be written as a single integer.
  - index2   : 1D array identifying which genes/microarrays belong to the
    second cluster. If the cluster contains only one gene, then
    index2 can also be written as a single integer.
  - transpose: if equal to 0, genes (rows) are clustered;
    if equal to 1, microarrays (columns) are clustered.
  - dist     : specifies the distance function to be used:

    - dist=='e': Euclidean distance
    - dist=='b': City Block distance
    - dist=='c': Pearson correlation
    - dist=='a': absolute value of the correlation
    - dist=='u': uncentered correlation
    - dist=='x': absolute uncentered correlation
    - dist=='s': Spearman's rank correlation
    - dist=='k': Kendall's tau

  - method   : specifies how the distance between two clusters is defined:

    - method=='a': the distance between the arithmetic means of the
      two clusters
    - method=='m': the distance between the medians of the two
      clusters
    - method=='s': the smallest pairwise distance between members
      of the two clusters
    - method=='x': the largest pairwise distance between members of
      the two clusters
    - method=='v': average of the pairwise distances between
      members of the clusters

  - transpose: if equal to 0: clusters of genes (rows) are considered;
    if equal to 1: clusters of microarrays (columns) are considered.

"""

        if transpose == 0:
            weight = self.eweight
        else:
            weight = self.gweight
        return clusterdistance(self.data, self.mask, weight,
                               index1, index2, method, dist, transpose)

    def distancematrix(self, transpose=0, dist='e'):
        """Calculate the distance matrix and return it as a list of arrays

  - transpose: if equal to 0: calculate the distances between genes (rows);
    if equal to 1: calculate the distances beteeen microarrays
    (columns).
  - dist     : specifies the distance function to be used:

    - dist=='e': Euclidean distance
    - dist=='b': City Block distance
    - dist=='c': Pearson correlation
    - dist=='a': absolute value of the correlation
    - dist=='u': uncentered correlation
    - dist=='x': absolute uncentered correlation
    - dist=='s': Spearman's rank correlation
    - dist=='k': Kendall's tau

Return value:
The distance matrix is returned as a list of 1D arrays containing the
distance matrix between the gene expression data. The number of columns
in each row is equal to the row number. Hence, the first row has zero
elements. An example of the return value is::

  matrix = [[],
            array([1.]),
            array([7., 3.]),
            array([4., 2., 6.])]

This corresponds to the distance matrix::

  [0., 1., 7., 4.]
  [1., 0., 3., 2.]
  [7., 3., 0., 6.]
  [4., 2., 6., 0.]

"""
        if transpose == 0:
            weight = self.eweight
        else:
            weight = self.gweight
        return distancematrix(self.data, self.mask, weight, transpose, dist)

    def save(self, jobname, geneclusters=None, expclusters=None):
        """Save the clustering results.

The saved files follow the convention for the Java TreeView program,
which can therefore be used to view the clustering result.

Arguments:

  - jobname:   The base name of the files to be saved. The filenames are
    jobname.cdt, jobname.gtr, and jobname.atr for
    hierarchical clustering, and jobname-K*.cdt,
    jobname-K*.kgg, jobname-K*.kag for k-means clustering
    results.
  - geneclusters=None:  For hierarchical clustering results, geneclusters
    is a Tree object as returned by the treecluster method.
    For k-means clustering results, geneclusters is a vector
    containing ngenes integers, describing to which cluster a
    given gene belongs. This vector can be calculated by
    kcluster.
  - expclusters=None:  For hierarchical clustering results, expclusters
    is a Tree object as returned by the treecluster method.
    For k-means clustering results, expclusters is a vector
    containing nexps integers, describing to which cluster a
    given experimental condition belongs. This vector can be
    calculated by kcluster.

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
        if geneclusters is not None and expclusters is not None and \
           type(geneclusters) != type(expclusters):
            raise ValueError("found one k-means and one hierarchical "
                           + "clustering solution in geneclusters and "
                           + "expclusters")
        gid = 0
        aid = 0
        filename = jobname
        postfix = ""
        if isinstance(geneclusters, Tree):
            # This is a hierarchical clustering result.
            geneindex = _savetree(jobname, geneclusters, gorder, 0)
            gid = 1
        elif geneclusters is not None:
            # This is a k-means clustering result.
            filename = jobname + "_K"
            k = max(geneclusters) + 1
            kggfilename = "%s_K_G%d.kgg" % (jobname, k)
            geneindex = self._savekmeans(kggfilename, geneclusters, gorder, 0)
            postfix = "_G%d" % k
        else:
            geneindex = numpy.argsort(gorder)
        if isinstance(expclusters, Tree):
            # This is a hierarchical clustering result.
            expindex = _savetree(jobname, expclusters, eorder, 1)
            aid = 1
        elif expclusters is not None:
            # This is a k-means clustering result.
            filename = jobname + "_K"
            k = max(expclusters) + 1
            kagfilename = "%s_K_A%d.kag" % (jobname, k)
            expindex = self._savekmeans(kagfilename, expclusters, eorder, 1)
            postfix += "_A%d" % k
        else:
            expindex = numpy.argsort(eorder)
        filename = filename + postfix
        self._savedata(filename, gid, aid, geneindex, expindex)

    def _savekmeans(self, filename, clusterids, order, transpose):
        # Save a k-means clustering solution
        if transpose == 0:
            label = self.uniqid
            names = self.geneid
        else:
            label = "ARRAY"
            names = self.expid
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
                        outputfile.write("%s\t%s\n" % (names[j], cluster))
                        sortedindex[counter] = j
                        counter += 1
                cluster += 1
        return sortedindex

    def _savedata(self, jobname, gid, aid, geneindex, expindex):
        # Save the clustered data.
        if self.genename is None:
            genename = self.geneid
        else:
            genename = self.genename
        (ngenes, nexps) = numpy.shape(self.data)
        with open(jobname+'.cdt', 'w') as outputfile:
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
                outputfile.write('GID\t')
            outputfile.write(self.uniqid)
            outputfile.write('\tNAME\tGWEIGHT')
            # Now add headers for data columns.
            for j in expindex:
                outputfile.write('\t%s' % self.expid[j])
            outputfile.write('\n')
            if aid:
                outputfile.write("AID")
                if gid:
                    outputfile.write('\t')
                outputfile.write("\t\t")
                for j in expindex:
                    outputfile.write('\tARRY%dX' % j)
                outputfile.write('\n')
            outputfile.write('EWEIGHT')
            if gid:
                outputfile.write('\t')
            outputfile.write('\t\t')
            for j in expindex:
                outputfile.write('\t%f' % eweight[j])
            outputfile.write('\n')
            for i in geneindex:
                if gid:
                    outputfile.write('GENE%dX\t' % i)
                outputfile.write("%s\t%s\t%f" %
                                 (self.geneid[i], genename[i], gweight[i]))
                for j in expindex:
                    outputfile.write('\t')
                    if mask[i, j]:
                        outputfile.write(str(self.data[i, j]))
                outputfile.write('\n')


def read(handle):
    """Read gene expression data from the file handle and return a Record.

The file should be in the file format defined for Michael Eisen's
Cluster/TreeView program.

"""
    return Record(handle)
