from Numeric import *

def print_row(data):
  print "[",
  for number in data[:-1]: print "%7.3f," % number,
  print "%7.3f]" % data[-1]

def print_matrix(data, mask):
  n,m = shape(data)
  for i in range(n):
    print "[",
    for j in range(m):
      if mask[i,j]: print "%7.3f" % data[i,j],
      else: print "   X   ",
    print "]"

def test_mean_median(module):
  if module=='Bio.Cluster':
    from Bio.Cluster import mean, median
  elif module=='Pycluster':
    from Pycluster import mean, median
  else:
    raise 'Unknown module name', module
  print "test_mean_median:"
  data1 = array([ 34.3, 3, 2 ])
  data2 = [ 5, 10 ,15, 20]
  data3 = [ 1, 2, 3, 5, 7, 11, 13, 17]
  data4 = [ 100, 19, 3, 1.5, 1.4, 1, 1, 1]

  for data in [data1, data2, data3, data4]:
    print "data =",
    print_row(data)
    print "mean is %7.3f; median is %7.3f" % (mean(data), median(data))
  print

def test_matrix_parse(module):
  if module=='Bio.Cluster':
    from Bio.Cluster import treecluster
  elif module=='Pycluster':
    from Pycluster import treecluster
  else:
    raise 'Unknown module name', module
  print "test_matrix_parse:"
  # Normal matrix, no errors
  data1 = array([[ 1.1, 1.2 ],
                 [ 1.4, 1.3 ],
                 [ 1.1, 1.5 ],
                 [ 2.0, 1.5 ],
                 [ 1.7, 1.9 ],
                 [ 1.7, 1.9 ],
                 [ 5.7, 5.9 ],
                 [ 5.7, 5.9 ],
                 [ 3.1, 3.3 ],
                 [ 5.4, 5.3 ],
                 [ 5.1, 5.5 ],
                 [ 5.0, 5.5 ],
                 [ 5.1, 5.2 ]])

  # Another normal matrix, no errors; written as a list
  data2 =  [[  1.1, 2.2, 3.3, 4.4, 5.5 ], 
            [  3.1, 3.2, 1.3, 2.4, 1.5 ], 
            [  4.1, 2.2, 0.3, 5.4, 0.5 ], 
            [ 12.1, 2.0, 0.0, 5.0, 0.0 ]]

  # Ragged matrix
  data3 =  [[ 91.1, 92.2, 93.3, 94.4, 95.5], 
            [ 93.1, 93.2, 91.3, 92.4 ], 
            [ 94.1, 92.2, 90.3 ], 
            [ 12.1, 92.0, 90.0, 95.0, 90.0 ]]

  # Matrix with bad cells
  data4 =  [ [ 7.1, 7.2, 7.3, 7.4, 7.5, ],
             [ 7.1, 7.2, 7.3, 7.4, 'snoopy' ], 
             [ 7.1, 7.2, 7.3, None, None]] 

  # Matrix with a bad row
  data5 =  [ [ 23.1, 23.2, 23.3, 23.4, 23.5], 
             None,
             [ 23.1, 23.0, 23.0, 23.0, 23.0]]

  # Various references that don't point to matrices at all
  data6 = "snoopy"
  data7 = {'a': [[2.3,1.2],[3.3,5.6]]}
  data8 = []
  data9 = [None]
  data10 = [[None]]

  try:
    result = treecluster(data1)
    print "Read data1 (correct)"
  except: "Error: treecluster failed to accept matrix data1"
  try:
    result = treecluster(data2)
    print "Read data2 (correct)"
  except: "Error: treecluster failed to accept matrix data2"
  try:
    result = treecluster(data3)
    print "Error: treecluster incorrectly accepted data3"
  except: print "Refused incorrect matrix data3"
  try:
    result = treecluster(data4)
    print "Error: treecluster incorrectly accepted data4"
  except: print "Refused incorrect matrix data4"
  try:
    result = treecluster(data5)
    print "Error: treecluster incorrectly accepted data5"
  except: print "Refused incorrect matrix data5"
  try:
    result = treecluster(data6)
    print "Error: treecluster incorrectly accepted data6"
  except: print "Refused incorrect matrix data6"
  try:
    result = treecluster(data7)
    print "Error: treecluster incorrectly accepted data7"
  except: print "Refused incorrect matrix data7"
  try:
    result = treecluster(data8)
    print "Error: treecluster incorrectly accepted data8"
  except: print "Refused incorrect matrix data8"
  try:
    result = treecluster(data9)
    print "Error: treecluster incorrectly accepted data9"
  except: print "Refused incorrect matrix data9"
  try:
    result = treecluster(data10)
    print "Error: treecluster incorrectly accepted data10"
  except: print "Refused incorrect matrix data10"
  print

def test_kcluster(module):
  if module=='Bio.Cluster':
    from Bio.Cluster import kcluster
  elif module=='Pycluster':
    from Pycluster import kcluster
  else:
    raise 'Unknown module name', module
  print "test_kcluster"
  nclusters = 3
  # First data set
  weight1 =  array([1,1,1,1,1])
  data1   =  array([[ 1.1, 2.2, 3.3, 4.4, 5.5],
                    [ 3.1, 3.2, 1.3, 2.4, 1.5], 
                    [ 4.1, 2.2, 0.3, 5.4, 0.5], 
                    [12.1, 2.0, 0.0, 5.0, 0.0]]) 
  mask1 =  array([[ 1, 1, 1, 1, 1], 
                  [ 1, 1, 1, 1, 1], 
                  [ 1, 1, 1, 1, 1], 
                  [ 1, 1, 1, 1, 1]]) 
  weight2 =  array([1,1])

  # Second data set
  data2 = array([[ 1.1, 1.2 ],
                 [ 1.4, 1.3 ],
                 [ 1.1, 1.5 ],
                 [ 2.0, 1.5 ],
                 [ 1.7, 1.9 ],
                 [ 1.7, 1.9 ],
                 [ 5.7, 5.9 ],
                 [ 5.7, 5.9 ],
                 [ 3.1, 3.3 ],
                 [ 5.4, 5.3 ],
                 [ 5.1, 5.5 ],
                 [ 5.0, 5.5 ],
                 [ 5.1, 5.2 ]])
  mask2 = array([[ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ]])

  # test first data set
  print "First data set"
  clusterid, centroids, error, nfound = kcluster (data1, nclusters=nclusters, mask=mask1, weight=weight1, transpose=0, npass=100, method='a', dist='e')
  print "Number of cluster ids is %d (should be %d)" % (len(clusterid), len(data1))
  print "Number of centroids is %d (should be %d)" % (len(centroids), nclusters)
  correct = [0,1,1,2]
  mapping = [clusterid[correct.index(i)] for i in range(nclusters)]
  same = 1
  for i in range(len(clusterid)):
    if clusterid[i]!=mapping[correct[i]]: same = 0
  if same: print "Correct clustering solution found."
  else: print "Wrong clustering solution found."
  print "Cluster centroids:"
  for i in range(nclusters):
    centroid = centroids[mapping[i]]
    for number in centroid:
      print "%7.3f " % number,
    print

  # test second data set
  print "Second data set"
  clusterid, centroids, error, nfound = kcluster (data2, nclusters=3, mask=mask2, weight=weight2, transpose=0, npass=100, method='a', dist='e')
  print "Number of cluster ids is %d (should be %d)" % (len(clusterid), len(data2))
  print "Number of centroids is %d (should be %d)" % (len(centroids), nclusters)
  correct = [0, 0, 0, 0, 0, 0, 1, 1, 2, 1, 1, 1, 1]
  mapping = [clusterid[correct.index(i)] for i in range(nclusters)]
  same = 1
  for i in range(len(clusterid)):
    if clusterid[i]!=mapping[correct[i]]: same = 0
  if same: print "Correct clustering solution found."
  else: print "Wrong clustering solution found."
  print "Cluster centroids:"
  for i in range(nclusters):
    centroid = centroids[mapping[i]]
    for number in centroid:
      print "%7.3f " % number,
    print
  print

def test_clusterdistance(module):
  if module=='Bio.Cluster':
    from Bio.Cluster import clusterdistance
  elif module=='Pycluster':
    from Pycluster import clusterdistance
  else:
    raise 'Unknown module name', module
  print "test_clusterdistance:"

  # First data set
  weight1 =  array([ 1,1,1,1,1 ])
  data1   =  array([[  1.1, 2.2, 3.3, 4.4, 5.5, ], 
                    [  3.1, 3.2, 1.3, 2.4, 1.5, ], 
                    [  4.1, 2.2, 0.3, 5.4, 0.5, ], 
                    [ 12.1, 2.0, 0.0, 5.0, 0.0, ]])
  mask1   = array([[ 1, 1, 1, 1, 1], 
                   [ 1, 1, 1, 1, 1], 
                   [ 1, 1, 1, 1, 1], 
                   [ 1, 1, 1, 1, 1]])

  # Second data set
  weight2 =  array([ 1,1 ])
  data2   =  array([[ 1.1, 1.2 ],
                    [ 1.4, 1.3 ],
                    [ 1.1, 1.5 ],
                    [ 2.0, 1.5 ],
                    [ 1.7, 1.9 ],
                    [ 1.7, 1.9 ],
                    [ 5.7, 5.9 ],
                    [ 5.7, 5.9 ],
                    [ 3.1, 3.3 ],
                    [ 5.4, 5.3 ],
                    [ 5.1, 5.5 ],
                    [ 5.0, 5.5 ],
                    [ 5.1, 5.2 ]])
  mask2 = array([[ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ]])
  # Cluster assignments
  c1 = [0]
  c2 = [1,2]
  c3 = [3]

  print "First data set:"
  print_matrix(data1, mask1)
  print "Clusters are", c1, c2, c3
  distance = clusterdistance(data1, mask=mask1, weight=weight1, index1=c1, index2=c2, dist='e', method='a', transpose=0);
  print "Distance between cluster", c1, "and", c2, "is %7.3f." % distance
  distance = clusterdistance(data1, mask=mask1, weight=weight1, index1=c1, index2=c3, dist='e', method='a', transpose=0);
  print "Distance between cluster", c1, "and", c3, "is %7.3f." % distance
  distance = clusterdistance(data1, mask=mask1, weight=weight1, index1=c2, index2=c3, dist='e', method='a', transpose=0);
  print "Distance between cluster", c2, "and", c3, "is %7.3f." % distance

  # Cluster assignments
  c1 = [ 0, 1, 2, 3 ]
  c2 = [ 4, 5, 6, 7 ]
  c3 = [ 8 ]
  print "Second data set:"
  print_matrix(data2, mask2)
  print "Clusters are", c1, c2, c3
  distance = clusterdistance(data2, mask=mask2, weight=weight2, index1=c1, index2=c2, dist='e', method='a', transpose=0);
  print "Distance between cluster", c1, "and", c2, "is %7.3f." % distance
  distance = clusterdistance(data2, mask=mask2, weight=weight2, index1=c1, index2=c3, dist='e', method='a', transpose=0);
  print "Distance between cluster", c1, "and", c3, "is %7.3f." % distance
  distance = clusterdistance(data2, mask=mask2, weight=weight2, index1=c2, index2=c3, dist='e', method='a', transpose=0);
  print "Distance between cluster", c2, "and", c3, "is %7.3f." % distance
  print

def test_treecluster(module):
  if module=='Bio.Cluster':
    from Bio.Cluster import treecluster
  elif module=='Pycluster':
    from Pycluster import treecluster
  else:
    raise 'Unknown module name', module
  print "test_treecluster:"
  # First data set
  weight1 =  [ 1,1,1,1,1 ]
  data1   =  array([[  1.1, 2.2, 3.3, 4.4, 5.5], 
                    [  3.1, 3.2, 1.3, 2.4, 1.5], 
                    [  4.1, 2.2, 0.3, 5.4, 0.5], 
                    [ 12.1, 2.0, 0.0, 5.0, 0.0]])
  mask1 = array([[ 1, 1, 1, 1, 1], 
                 [ 1, 1, 1, 1, 1], 
                 [ 1, 1, 1, 1, 1], 
                 [ 1, 1, 1, 1, 1]])

  # Second data set
  weight2 =  [ 1,1 ]
  data2 = array([[ 0.8223, 0.9295 ],
                 [ 1.4365, 1.3223 ],
                 [ 1.1623, 1.5364 ],
                 [ 2.1826, 1.1934 ],
                 [ 1.7763, 1.9352 ],
                 [ 1.7215, 1.9912 ],
                 [ 2.1812, 5.9935 ],
                 [ 5.3290, 5.9452 ],
                 [ 3.1491, 3.3454 ],
                 [ 5.1923, 5.3156 ],
                 [ 4.7735, 5.4012 ],
                 [ 5.1297, 5.5645 ],
                 [ 5.3934, 5.1823 ]])
  mask2 = array([[ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ]])

  # test first data set
  print "First data set:"
  print_matrix(data1, mask1)
  print "Pairwise average-linkage clustering"
  result, linkdist = treecluster(data=data1, mask=mask1, weight=weight1, applyscale=0, transpose=0, method='a', dist='e')
  print "Number of nodes is %d (should be %d)" % (len(result), len(data1)-1)
  print "Number of link distances is %d (should be %d)" % (len(linkdist), len(data1)-1)
  for i in range(len(result)):
    print "Node %3d joins node %3d with node %3d; link distance is %7.3f" % (i, result[i][0], result[i][1], linkdist[i])

  print "Pairwise single-linkage clustering"
  result, linkdist = treecluster(data=data1, mask=mask1, weight=weight1, applyscale=0, transpose=0, method='s', dist='e')
  print "Number of nodes is %d (should be %d)" % (len(result), len(data1)-1)
  print "Number of link distances is %d (should be %d)" % (len(linkdist), len(data1)-1)
  for i in range(len(result)):
    print "Node %3d joins node %3d with node %3d; link distance is %7.3f" % (i, result[i][0], result[i][1], linkdist[i])

  print "Pairwise centroid-linkage clustering"
  result, linkdist = treecluster(data=data1, mask=mask1, weight=weight1, applyscale=0, transpose=0, method='c', dist='e')
  print "Number of nodes is %d (should be %d)" % (len(result), len(data1)-1)
  print "Number of link distances is %d (should be %d)" % (len(linkdist), len(data1)-1)
  for i in range(len(result)):
    print "Node %3d joins node %3d with node %3d; link distance is %7.3f" % (i, result[i][0], result[i][1], linkdist[i])

  print "Pairwise maximum-linkage clustering"
  result, linkdist = treecluster(data=data1, mask=mask1, weight=weight1, applyscale=0, transpose=0, method='m', dist='e')
  print "Number of nodes is %d (should be %d)" % (len(result), len(data1)-1)
  print "Number of link distances is %d (should be %d)" % (len(linkdist), len(data1)-1)
  for i in range(len(result)):
    print "Node %3d joins node %3d with node %3d; link distance is %7.3f" % (i, result[i][0], result[i][1], linkdist[i])

  # Test second data set
  print "Second data set:"
  print "Pairwise average-linkage clustering"
  result, linkdist = treecluster(data=data2, mask=mask2, weight=weight2, applyscale=0, transpose=0, method='a', dist='e')
  print "Number of nodes is %d (should be %d)" % (len(result), len(data2)-1)
  print "Number of link distances is %d (should be %d)" % (len(linkdist), len(data2)-1)
  for i in range(len(result)):
    print "Node %3d joins node %3d with node %3d; link distance is %7.3f" % (i, result[i][0], result[i][1], linkdist[i])

  print "Pairwise single-linkage clustering"
  result, linkdist = treecluster(data=data2, mask=mask2, weight=weight2, applyscale=0, transpose=0, method='s', dist='e')
  print "Number of nodes is %d (should be %d)" % (len(result), len(data2)-1)
  print "Number of link distances is %d (should be %d)" % (len(linkdist), len(data2)-1)
  for i in range(len(result)):
    print "Node %3d joins node %3d with node %3d; link distance is %7.3f" % (i, result[i][0], result[i][1], linkdist[i])

  print "Pairwise centroid-linkage clustering"
  result, linkdist = treecluster(data=data2, mask=mask2, weight=weight2, applyscale=0, transpose=0, method='c', dist='e')
  print "Number of nodes is %d (should be %d)" % (len(result), len(data2)-1)
  print "Number of link distances is %d (should be %d)" % (len(linkdist), len(data2)-1)
  for i in range(len(result)):
    print "Node %3d joins node %3d with node %3d; link distance is %7.3f" % (i, result[i][0], result[i][1], linkdist[i])

  print "Pairwise maximum-linkage clustering"
  result, linkdist = treecluster(data=data2, mask=mask2, weight=weight2, applyscale=0, transpose=0, method='m', dist='e')
  print "Number of nodes is %d (should be %d)" % (len(result), len(data2)-1)
  print "Number of link distances is %d (should be %d)" % (len(linkdist), len(data2)-1)
  for i in range(len(result)):
    print "Node %3d joins node %3d with node %3d; link distance is %7.3f" % (i, result[i][0], result[i][1], linkdist[i])
  print

def test_somcluster(module):
  if module=='Bio.Cluster':
    from Bio.Cluster import somcluster
  elif module=='Pycluster':
    from Pycluster import somcluster
  else:
    raise 'Unknown module name', module
  print "test_somcluster:"

  # First data set
  weight1 = [ 1,1,1,1,1 ]
  data1 = array([[  1.1, 2.2, 3.3, 4.4, 5.5], 
                 [  3.1, 3.2, 1.3, 2.4, 1.5], 
                 [  4.1, 2.2, 0.3, 5.4, 0.5], 
                 [ 12.1, 2.0, 0.0, 5.0, 0.0]])
  mask1 = array([[ 1, 1, 1, 1, 1], 
                 [ 1, 1, 1, 1, 1], 
                 [ 1, 1, 1, 1, 1], 
                 [ 1, 1, 1, 1, 1]])

  # Second data set
  weight2 =  [ 1,1 ]
  data2 = array([[ 1.1, 1.2 ],
                 [ 1.4, 1.3 ],
                 [ 1.1, 1.5 ],
                 [ 2.0, 1.5 ],
                 [ 1.7, 1.9 ],
                 [ 1.7, 1.9 ],
                 [ 5.7, 5.9 ],
                 [ 5.7, 5.9 ],
                 [ 3.1, 3.3 ],
                 [ 5.4, 5.3 ],
                 [ 5.1, 5.5 ],
                 [ 5.0, 5.5 ],
                 [ 5.1, 5.2 ]])
  mask2 = array([[ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ],
                 [ 1, 1 ]])

  print "First data set:"
  clusterid, celldata = somcluster(data=data1, mask=mask1, weight=weight1, transpose=0, nxgrid=10, nygrid=10, inittau=0.02, niter=100, dist='e')
  print "Number of cluster ids is %d (should be %d)" % (len(clusterid), len(data1))
  print "Grid is %d-dimensional (should be 2-dimensional)" % len(clusterid[0])

  print "Second data set:"
  clusterid, celldata = somcluster(data=data2, mask=mask2, weight=weight2, transpose=0, nxgrid=10, nygrid=10, inittau=0.02, niter=100, dist='e')
  print "Number of cluster ids is %d (should be %d)" % (len(clusterid), len(data2))
  print "Grid is %d-dimensional (should be 2-dimensional)" % len(clusterid[0])
  print

def run_tests(module="Pycluster"):
  if module==[]: module = "Bio.Cluster"
  test_mean_median(module)
  test_matrix_parse(module)
  test_kcluster(module)
  test_clusterdistance(module)
  test_treecluster(module)
  test_somcluster(module)
