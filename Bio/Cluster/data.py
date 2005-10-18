from Numeric import *

def readdatafile(filename):
  """readdatafile reads a file containing gene expression data following
     Michael Eisen's format for Cluster/TreeView. readdatafile returns a
     tuple containing:
data:     a matrix containing the gene expression data
mask:     a matrix containing only 1's and 0's, denoting which values are
          present (1) or missing (0). If all elements of mask are one
          (no missing data), then None is returned instead of the mask
geneid:   a list containing a unique identifier for each gene
          (e.g., ORF name)
genename: a list containing an additional description for each gene
          (e.g., gene name)
gweight:  the weight to be used for each gene when calculating the distance
gorder:   an array of real numbers indicating the preferred order of the
          genes in the output file
expid:    a list containing a unique identifier for each experimental
          condition
eweight:  the weight to be used for each experimental condition when
          calculating the distance
eorder:   an array of real numbers indication the preferred order in the
          output file of the experimental conditions
uniqid:   the string that was used instead of UNIQID in the input file."""
  inputfile = open(filename)
  lines = inputfile.readlines()
  inputfile.close()
  for i in range(len(lines)):
    line = lines[i]
    line = line.replace("\r","\n") # This should work for Macintosh also
    oldline = ''
    while line != oldline:
      oldline = line
      line = line.replace("\t\t","\tmissing\t")
      line = line.replace("\t\n","\tmissing\n")
    line = line.replace("\n","")
    lines[i] = line.split("\t")
  uniqid = lines[0][0]
  column = 1
  NAME = 0
  if lines[0][column]=="NAME":
    NAME = column
    column = column + 1
  GWEIGHT = 0
  if lines[0][column]=="GWEIGHT":
    GWEIGHT = column
    column = column + 1
  GORDER = 0
  if lines[0][column]=="GORDER":
    GORDER = column
    column = column + 1
  row = 1
  EWEIGHT = 0
  if lines[row][0]=="EWEIGHT":
    EWEIGHT = row
    row = row + 1
  EORDER = 0
  if lines[row][0]=="EORDER":
    EORDER = row
    row = row + 1
  ngenes = len(lines) - row
  geneid = [''] * ngenes
  for i in range(ngenes):
    geneid[i] = lines[row+i][0]
  expid = lines[0][column:]
  nexps = len(expid)
  if EWEIGHT:
    eweight = zeros(nexps,'d')
    for i in range(nexps):
      eweight[i] = float(lines[EWEIGHT][column+i])
  else:
    eweight = None
  if EORDER:
    eorder = zeros(nexps,'d')
    for i in range(nexps):
      eorder[i] = float(lines[EORDER][column+i])
  else:
    eorder = None
  if NAME:
    genename = [''] * ngenes
    for i in range(ngenes):
      genename[i] = lines[row+i][NAME]
  else:
    genename = None
  if GWEIGHT:
    gweight = zeros(ngenes,'d')
    for i in range(ngenes):
      gweight[i] = float(lines[row+i][GWEIGHT])
  else:
    gweight = None
  if GORDER:
    gorder = zeros(ngenes,'d')
    for i in range(ngenes):
      gorder[i] = float(lines[row+i][GORDER])
  else:
    gorder = None
  data = zeros((ngenes,nexps),'d')
  mask = ones((ngenes,nexps))
  needmask = 0
  for i in range(ngenes):
    line = lines[row+i][column:]
    for j in range(nexps):
      if line[j]=="missing":
        mask[i,j] = 0
        needmask = 1
      else: data[i,j] = float(line[j])
  if not needmask: mask = None
  return data, mask, geneid, genename, gweight, gorder, expid, eweight, eorder, uniqid

def savedata(jobname, gid, aid, data, mask, geneid, expid, geneindex, expindex, genename, gweight, eweight, uniqID):
  (ngenes, nexps) = shape(data)
  outputfile = open(jobname+'.cdt', 'w')
  if not outputfile: return "Error: Unable to open output file"
  if not mask: mask = ones((ngenes,nexps))
  if not gweight: gweight = ones(ngenes)
  if not eweight: eweight = ones(nexps)
  if gid: outputfile.write ('GID\t')
  outputfile.write(uniqID)
  outputfile.write('\tNAME\tGWEIGHT')
  # Now add headers for data columns
  for column in range(nexps):
    outputfile.write('\t')
    outputfile.write(expid[expindex[column]])
  outputfile.write('\n')
  if aid:
    outputfile.write("AID")
    if gid: outputfile.write('\t')
    outputfile.write("\t\t")
    for column in range(nexps): outputfile.write ('\tARRY'+str(expindex[column]))
    outputfile.write('\n')
  outputfile.write('EWEIGHT')
  if gid: outputfile.write('\t')
  outputfile.write('\t\t')
  for column in range(nexps):
    outputfile.write('\t' + str(eweight[expindex[column]]))
  outputfile.write('\n')
  for row in range(ngenes):
    index = geneindex[row]
    if gid: outputfile.write('GENE'+str(index)+'\t')
    outputfile.write(geneid[index]+'\t'+genename[index]+'\t'+str(gweight[index]))
    for column in range(nexps):
      columnindex = expindex[column]
      outputfile.write('\t')
      if mask[index][columnindex]: outputfile.write(str(data[index][columnindex]));
    outputfile.write('\n')
  outputfile.close()

def treesort(order, nodeorder, nodecounts, NodeElement):
  nNodes = len(NodeElement)
  nElements = nNodes + 1
  neworder = zeros(nElements,'d')
  clusterids = range(nElements)
  for i in range(nNodes):
    i1 = NodeElement[i][0]
    i2 = NodeElement[i][1]
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
  return argsort(neworder)

def savetree(jobname, clusters, linkdist, order, transpose):
  if transpose==0:
    extension = ".gtr"
    keyword = "GENE"
  else:
    extension = ".atr"
    keyword = "ARRY"
  nnodes = len(clusters)
  outputfile = open(jobname+extension, "w");
  nodeindex = 0
  nodeID = [''] * (nnodes)
  nodecounts = zeros(nnodes)
  nodeorder = zeros(nnodes,'d')
  nodedist = array(linkdist)
  for nodeindex in range(nnodes):
    min1 = clusters[nodeindex][0];
    min2 = clusters[nodeindex][1];
    nodeID[nodeindex] = "NODE" + str(nodeindex+1)
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
      outputfile.write(keyword + str(min1) + "\t")
    if min2 < 0:
      index2 = -min2-1
      order2 = nodeorder[index2]
      counts2 = nodecounts[index2]
      outputfile.write(nodeID[index2]+"\t")
      nodedist[nodeindex] = max(nodedist[nodeindex],nodedist[index2])
    else:
      order2 = order[min2];
      counts2 = 1;
      outputfile.write(keyword + str(min2) + "\t")
    outputfile.write(str(1.0-nodedist[nodeindex]))
    outputfile.write("\n")
    nodecounts[nodeindex] = counts1 + counts2
    nodeorder[nodeindex] = (counts1*order1+counts2*order2) / (counts1+counts2)
  outputfile.close()
  # Now set up order based on the tree structure
  index = treesort(order, nodeorder, nodecounts, clusters)
  return index

def savekmeans(jobname, k, label, names, clusterids, order, transpose):
  if transpose==0:
    filename = jobname + "_K_G" + str(k) + ".kgg"
  else:
    filename = jobname + "_K_A" + str(k) + ".kag"
  outputfile = open(filename, "w");
  if not outputfile: raise "Error: Unable to open output file"
  outputfile.write(label + "\tGROUP\n")
  index = argsort(order)
  nElements = len(names)
  indeces = zeros(nElements)
  counter = 0
  for cluster in range(k):
    for element in range(nElements):
      j = index[element];
      if clusterids[j]==cluster:
        outputfile.write(names[j] + "\t" + str(cluster) + "\n")
        indeces[counter] = j
        counter = counter + 1
  outputfile.close();
  return indeces

def writeclusterfiles(jobname, data, geneid, expid, mask=None, geneclusters=None, genelinkdist=None, expclusters=None, explinkdist=None, gorder=None, eorder=None, genename=None, gweight=None, eweight=None, uniqid='UNIQID'):
  """writeclusterfiles saves the clustering results. The saved files
follow the convention for Java TreeView program, which can therefore be
used to view the clustering result. writeclusterfiles takes the
following arguments:
jobname:   The base name of the files to be saved. The filenames are
           jobname.cdt, jobname.gtr, and jobname.atr for hierarchical
           clustering, and jobname-K*.cdt, jobname-K*.kgg, jobname-K*.kag
           for k-means clustering results
data:      The array containing the gene expression data
geneid:    A list of strings containing an unique identifier for each
           gene (e.g., ORF name)
expid:     A list of strings containing an unique identifier for each
           experimental condition
mask=None: This array describes which elements in data are missing. If
           there are no missing data, then mask can be set to None, which
           is the default value
geneclusters=None:  For hierarchical clustering results, geneclusters
           is an (ngenes-1 x 2) array that describes the hierarchical
           clustering result for genes. This array can be calculated
           by the hierarchical clustering methods implemented in
           treecluster.
           For k-means clustering results, geneclusters is a vector
           containing ngenes integers, describing to which cluster a
           given gene belongs. This vector can be calculated by kcluster.
genelinkdist=None:  An array with (ngenes-1) elements containing the
           distance between the subnodes that were joined for each node.
           genelinkdist is calculated by the hierarchical clustering
           methods implemented in treecluster.
           genelinkdist is required only if geneclusters is given and
           contains an hierarchical clustering result. 
expclusters=None:  For hierarchical clustering results, expclusters
           is an (nexps-1 x 2) array that describes the hierarchical
           clustering result for experimental conditions. This array can
           be calculated by the hierarchical clustering methods implemented
           in treecluster.
           For k-means clustering results, expclusters is a vector
           containing nexps integers, describing to which cluster a
           given experimental condition belongs. This vector can be
           calculated by kcluster.
explinkdist=None:  An array with (nexps-1) elements containing the
           distance between the subnodes that were joined for each node.
           explinkdist is calculated by the hierarchical clustering
           methods implemented in treecluster.
           explinkdist is required only if expclusters is given and
           contains an hierarchical clustering result. 
gorder=None: an array of real numbers indicating the preferred order of
           the genes in the output file.
eorder=None: an array of real numbers indication the preferred order in
           the output file of the experimental conditions
genename:  a list containing an additional description for each gene
           (e.g., gene name)
uniqid:    the string that was used instead of UNIQID in the input file.
 """
  (ngenes,nexps) = shape(data)
  if genename==None: genename = geneid
  if gorder==None: gorder = arange(ngenes)
  if eorder==None: eorder = arange(nexps)
  if rank(geneclusters)==1 and rank(expclusters)==2:
    print "Error: k-means clustering result for genes combined with"
    print "hierarchical clustering result for experimental conditions"
    return
  if rank(geneclusters)==2 and rank(expclusters)==1:
    print "Error: Hierarchical clustering result for genes combined with"
    print "k-means clustering result for experimental conditions"
    return
  gid = 0
  aid = 0
  filename = jobname
  if rank(geneclusters)==2 or rank(expclusters)==2:
    # Hierarchical clustering result
    if geneclusters:
      if not genelinkdist:
        print "Error: For genes, hierarchical clustering result is given but"
        print "the node distances are missing (genelinkdist==None)."
        return
      geneindex = savetree(jobname, geneclusters, genelinkdist, gorder, 0)
      gid = 1
    else:
      geneindex = argsort(gorder)
    if expclusters:
      if not explinkdist:
        print "Error: For experimental conditions, hierarchical clustering"
        print "result is given but the node distances are missing"
        print "(explinkdist==None)."
        return
      expindex = savetree(jobname, expclusters, explinkdist, eorder, 1)
      aid = 1
    else:
      expindex = argsort(eorder)
  elif rank(geneclusters)==1 or rank(expclusters)==1:
    # k-means clustering result
    filename = filename + "_K"
    if geneclusters:
      k = max(geneclusters) + 1
      geneindex = savekmeans(jobname, k, uniqid, geneid, geneclusters, gorder, 0)
      filename = filename + "_G" + str(k)
    else:
      geneindex = argsort(gorder)
    if expclusters:
      k = max(expclusters) + 1
      expindex = savekmeans(jobname, k, "ARRAY", expid, expclusters, eorder, 1)
      filename = filename + "_A" + str(k)
    else:
      expindex = argsort(eorder)
  savedata(filename,gid,aid,data,mask,geneid,expid,geneindex,expindex,genename,gweight,eweight,uniqid)
