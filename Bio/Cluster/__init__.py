import numpy
from cluster import *


def _treesort(order, nodeorder, nodecounts, tree):
  nNodes = len(tree)
  nElements = nNodes + 1
  neworder = numpy.zeros(nElements)
  clusterids = range(nElements)
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
    # If order1 and order2 are equal, their order is determined by the order in which they were clustered
    if i1 < i2:
      if order1 < order2:
        increase = count1
      else:
        increase = count2
      for j in range(nElements):
        clusterid = clusterids[j]
        if clusterid==i1 and order1>=order2: neworder[j] += increase
        if clusterid==i2 and order1<order2: neworder[j] += increase
        if clusterid==i1 or clusterid==i2: clusterids[j] = -i-1
    else:
      if order1<=order2:
        increase = count1
      else:
        increase = count2
      for j in range(nElements):
        clusterid = clusterids[j]
        if clusterid==i1 and order1>order2: neworder[j] += increase
        if clusterid==i2 and order1<=order2: neworder[j] += increase
        if clusterid==i1 or clusterid==i2: clusterids[j] = -i-1
  return numpy.argsort(neworder)

def _savetree(jobname, tree, order, transpose):
  if transpose==0:
    extension = ".gtr"
    keyword = "GENE"
  else:
    extension = ".atr"
    keyword = "ARRY"
  nnodes = len(tree)
  outputfile = open(jobname+extension, "w");
  nodeindex = 0
  nodeID = [''] * (nnodes)
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
      nodedist[nodeindex] = max(nodedist[nodeindex],nodedist[index1])
    else:
      order1 = order[min1]
      counts1 = 1
      outputfile.write("%s%dX\t" % (keyword, min1))
    if min2 < 0:
      index2 = -min2-1
      order2 = nodeorder[index2]
      counts2 = nodecounts[index2]
      outputfile.write(nodeID[index2]+"\t")
      nodedist[nodeindex] = max(nodedist[nodeindex],nodedist[index2])
    else:
      order2 = order[min2];
      counts2 = 1;
      outputfile.write("%s%dX\t" % (keyword, min2))
    outputfile.write(str(1.0-nodedist[nodeindex]))
    outputfile.write("\n")
    nodecounts[nodeindex] = counts1 + counts2
    nodeorder[nodeindex] = (counts1*order1+counts2*order2) / (counts1+counts2)
  outputfile.close()
  # Now set up order based on the tree structure
  index = _treesort(order, nodeorder, nodecounts, tree)
  return index

class Record:
  """A Record stores the gene expression data and related information
     contained in a data file following the file format defined for
     Michael Eisen's Cluster/TreeView program. A Record
     has the following members:
data:     a matrix containing the gene expression data
mask:     a matrix containing only 1's and 0's, denoting which values
          are present (1) or missing (0). If all elements of mask are
          one (no missing data), then mask is set to None.
geneid:   a list containing a unique identifier for each gene
          (e.g., ORF name)
genename: a list containing an additional description for each gene
          (e.g., gene name)
gweight:  the weight to be used for each gene when calculating the
          distance
gorder:   an array of real numbers indicating the preferred order of the
          genes in the output file
expid:    a list containing a unique identifier for each experimental
          condition
eweight:  the weight to be used for each experimental condition when
          calculating the distance
eorder:   an array of real numbers indication the preferred order in the
          output file of the experimental conditions
uniqid:   the string that was used instead of UNIQID in the input file."""
  def __init__(self, handle=None):
    """Reads a data file in the format corresponding to Michael Eisen's
Cluster/TreeView program, and stores the data in a Record object"""
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
    if not handle: return
    lines = handle.readlines()
    lines = [line.strip("\r\n").split("\t") for line in lines]
    line = lines[0]
    n = len(line)
    self.uniqid = line[0]
    self.expid = []
    cols = {0: "GENEID"}
    for word in line[1:]:
      if word=="NAME":
        cols[line.index(word)] = word
        self.genename = []
      elif word=="GWEIGHT":
        cols[line.index(word)] = word
        self.gweight = []
      elif word=="GORDER":
        cols[line.index(word)] = word
        self.gorder = []
      else: self.expid.append(word)
    self.geneid = []
    self.data = []
    self.mask = []
    needmask = 0
    for line in lines[1:]:
      assert len(line)==n, "Line with %d columns found (expected %d)" % (len(line), n)
      if line[0]=="EWEIGHT":
        i = max(cols) + 1
        self.eweight = map(float, line[i:])
        continue
      if line[0]=="EORDER":
        i = max(cols) + 1
        self.eorder = map(float, line[i:])
        continue
      rowdata = []
      rowmask = []
      n = len(line)
      for i in range(n):
        word = line[i]
        if i in cols:
          if cols[i]=="GENEID": self.geneid.append(word)
          if cols[i]=="NAME": self.genename.append(word)
          if cols[i]=="GWEIGHT": self.gweight.append(float(word))
          if cols[i]=="GORDER": self.gorder.append(float(word))
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
    if needmask: self.mask = numpy.array(self.mask, int)
    else: self.mask = None
    if self.gweight: self.gweight = numpy.array(self.gweight)
    if self.gorder: self.gorder = numpy.array(self.gorder)

  def treecluster(self, transpose=0, method='m', dist='e'):
    if transpose==0: weight = self.eweight
    else: weight = self.gweight
    return treecluster(self.data, self.mask, weight, transpose, method, dist)

  def kcluster(self, nclusters=2, transpose=0, npass=1, method='a', dist='e', initialid=None):
    if transpose==0: weight = self.eweight
    else: weight = self.gweight
    clusterid, error, nfound = kcluster(self.data, nclusters, self.mask, weight, transpose, npass, method, dist, initialid)
    return clusterid, error, nfound

  def somcluster(self, transpose=0, nxgrid=2, nygrid=1, inittau=0.02, niter=1, dist='e'):
    if transpose==0: weight = self.eweight
    else: weight = self.gweight
    clusterid, celldata = somcluster(self.data, self.mask, weight, transpose, nxgrid, nygrid, inittau, niter, dist)
    return clusterid, celldata

  def clustercentroids(self, clusterid=None, method='a', transpose=0):
    cdata, cmask = clustercentroids(self.data, self.mask, clusterid, method, transpose)
    return cdata, cmask

  def clusterdistance(self, index1=[0], index2=[0], method='a', dist='e',
                      transpose=0):
    if transpose==0: weight = self.eweight
    else: weight = self.gweight
    return clusterdistance(self.data, self.mask, weight, index1, index2, method, dist, transpose)

  def distancematrix(self, transpose=0, dist='e'):
    if transpose==0: weight = self.eweight
    else: weight = self.gweight
    return distancematrix(self.data, self.mask, weight, transpose, dist)

  def save(self, jobname, geneclusters=None, expclusters=None):
    """save(jobname, geneclusters=None, expclusters=None)
saves the clustering results. The saved files follow the convention
for Java TreeView program, which can therefore be used to view the
clustering result.
Arguments:
jobname:   The base name of the files to be saved. The filenames are
           jobname.cdt, jobname.gtr, and jobname.atr for
           hierarchical clustering, and jobname-K*.cdt,
           jobname-K*.kgg, jobname-K*.kag for k-means clustering
           results.
geneclusters=None:  For hierarchical clustering results,
           geneclusters is an (ngenes-1 x 2) array that describes
           the hierarchical clustering result for genes. This array
           can be calculated by the hierarchical clustering methods
           implemented in treecluster.
           For k-means clustering results, geneclusters is a vector
           containing ngenes integers, describing to which cluster a
           given gene belongs. This vector can be calculated by
           kcluster.
expclusters=None:  For hierarchical clustering results, expclusters
           is an (nexps-1 x 2) array that describes the hierarchical
           clustering result for experimental conditions. This array
           can be calculated by the hierarchical clustering methods
           implemented in treecluster.
           For k-means clustering results, expclusters is a vector
           containing nexps integers, describing to which cluster a
           given experimental condition belongs. This vector can be
           calculated by kcluster.
"""
    (ngenes,nexps) = shape(self.data)
    if self.gorder==None: gorder = numpy.arange(ngenes)
    else: gorder = self.gorder
    if self.eorder==None: eorder = numpy.arange(nexps)
    else: eorder = self.eorder
    if geneclusters and expclusters:
      assert type(geneclusters)==type(expclusters), "found one k-means and one hierarchical clustering solution in geneclusters and expclusters"
    gid = 0
    aid = 0
    filename = jobname
    postfix = ""
    if type(geneclusters)==Tree:
      # Hierarchical clustering result
      geneindex = _savetree(jobname, geneclusters, gorder, 0)
      gid = 1
    elif geneclusters:
      # k-means clustering result
      filename = jobname + "_K"
      k = max(geneclusters+1)
      kggfilename = "%s_K_G%d.kgg" % (jobname, k)
      geneindex = self._savekmeans(kggfilename, geneclusters, gorder, 0)
      postfix = "_G%d" % k
    else:
      geneindex = numpy.argsort(gorder)
    if type(expclusters)==Tree:
      # Hierarchical clustering result
      expindex = _savetree(jobname, expclusters, eorder, 1)
      aid = 1
    elif expclusters:
      # k-means clustering result
      filename = jobname + "_K"
      k = max(expclusters+1)
      kagfilename = "%s_K_A%d.kag" % (jobname, k)
      expindex = self._savekmeans(kagfilename, expclusters, eorder, 1)
      postfix += "_A%d" % k
    else:
      expindex = numpy.argsort(eorder)
    filename = filename + postfix
    self._savedata(filename,gid,aid,geneindex,expindex)

  def _savekmeans(self, filename, clusterids, order, transpose):
    if transpose==0:
      label = self.uniqid
      names = self.geneid
    else:
      label = "ARRAY"
      names = self.expid
    outputfile = open(filename, "w");
    if not outputfile: raise "Error: Unable to open output file"
    outputfile.write(label + "\tGROUP\n")
    index = numpy.argsort(order)
    n = len(names)
    sortedindex = numpy.zeros(n, int)
    counter = 0
    cluster = 0
    while counter < n:
      for j in index:
        if clusterids[j]==cluster:
          outputfile.write("%s\t%s\n" % (names[j], cluster))
          sortedindex[counter] = j
          counter+=1
      cluster+=1
    outputfile.close();
    return sortedindex

  def _savedata(self, jobname, gid, aid, geneindex, expindex):
    if self.genename==None: genename = self.geneid
    else: genename = self.genename
    (ngenes, nexps) = numpy.shape(self.data)
    outputfile = open(jobname+'.cdt', 'w')
    if not outputfile: return "Error: Unable to open output file"
    if self.mask: mask = self.mask
    else: mask = numpy.ones((ngenes,nexps), int)
    if self.gweight: gweight = self.gweight
    else: gweight = numpy.ones(ngenes)
    if self.eweight: eweight = self.eweight
    else: eweight = numpy.ones(nexps)
    if gid: outputfile.write ('GID\t')
    outputfile.write(self.uniqid)
    outputfile.write('\tNAME\tGWEIGHT')
    # Now add headers for data columns
    for j in expindex: outputfile.write('\t%s' % self.expid[j])
    outputfile.write('\n')
    if aid:
      outputfile.write("AID")
      if gid: outputfile.write('\t')
      outputfile.write("\t\t")
      for j in expindex: outputfile.write ('\tARRY%dX' % j)
      outputfile.write('\n')
    outputfile.write('EWEIGHT')
    if gid: outputfile.write('\t')
    outputfile.write('\t\t')
    for j in expindex: outputfile.write('\t%f' % eweight[j])
    outputfile.write('\n')
    for i in geneindex:
      if gid: outputfile.write('GENE%dX\t' % i)
      outputfile.write("%s\t%s\t%f" % (self.geneid[i], genename[i], gweight[i]))
      for j in expindex:
        outputfile.write('\t')
        if mask[i][j]: outputfile.write(str(self.data[i][j]))
      outputfile.write('\n')
    outputfile.close()

def read(handle):
    return Record(handle)
