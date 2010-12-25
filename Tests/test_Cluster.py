# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import unittest

try:
    from Bio import Cluster
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError("If you want to use Bio.Cluster, "
                                       "install NumPy first and then "
                                       "reinstall Biopython")

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(\
        "Install NumPy if you want to use Bio.Cluster")

class TestCluster(unittest.TestCase):

    module = 'Bio.Cluster'

    def test_median_mean(self):
        if TestCluster.module=='Bio.Cluster':
            from Bio.Cluster import mean, median
        elif TestCluster.module=='Pycluster':
            from Pycluster import mean, median

        data = numpy.array([ 34.3, 3, 2 ])
        self.assertAlmostEqual(mean(data), 13.1, places=3)
        self.assertAlmostEqual(median(data), 3.0, places=3)

        data = [ 5, 10, 15, 20]
        self.assertAlmostEqual(mean(data), 12.5, places=3)
        self.assertAlmostEqual(median(data), 12.5, places=3)

        data = [ 1, 2, 3, 5, 7, 11, 13, 17]
        self.assertAlmostEqual(mean(data), 7.375, places=3)
        self.assertAlmostEqual(median(data), 6.0, places=3)

        data = [ 100, 19, 3, 1.5, 1.4, 1, 1, 1]
        self.assertAlmostEqual(mean(data), 15.988, places=3)
        self.assertAlmostEqual(median(data), 1.45, places=3)
      

    def test_matrix_parse(self):
        if TestCluster.module=='Bio.Cluster':
            from Bio.Cluster import treecluster
        elif TestCluster.module=='Pycluster':
            from Pycluster import treecluster

        # Normal matrix, no errors
        data1 = numpy.array([[ 1.1, 1.2 ],
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
      
        try:
            treecluster(data1)
        except:
            self.fail("treecluster failed to accept matrix data1")

        try:
            treecluster(data2)
        except:
            self.fail("treecluster failed to accept matrix data2")

        self.assertRaises(TypeError, lambda : treecluster(data3))
        self.assertRaises(TypeError, lambda : treecluster(data4))
        self.assertRaises(TypeError, lambda : treecluster(data5))
        self.assertRaises(TypeError, lambda : treecluster(data6))
        self.assertRaises(TypeError, lambda : treecluster(data7))
        self.assertRaises(TypeError, lambda : treecluster(data8))
        self.assertRaises(TypeError, lambda : treecluster(data9))

    def test_kcluster(self):
        if TestCluster.module=='Bio.Cluster':
            from Bio.Cluster import kcluster
        elif TestCluster.module=='Pycluster':
            from Pycluster import kcluster

        nclusters = 3
        # First data set
        weight = numpy.array([1,1,1,1,1])
        data   = numpy.array([[ 1.1, 2.2, 3.3, 4.4, 5.5],
                              [ 3.1, 3.2, 1.3, 2.4, 1.5], 
                              [ 4.1, 2.2, 0.3, 5.4, 0.5], 
                              [12.1, 2.0, 0.0, 5.0, 0.0]]) 
        mask =  numpy.array([[ 1, 1, 1, 1, 1], 
                             [ 1, 1, 1, 1, 1], 
                             [ 1, 1, 1, 1, 1], 
                             [ 1, 1, 1, 1, 1]], int) 
      
        clusterid, error, nfound = kcluster(data, nclusters=nclusters,
                                            mask=mask, weight=weight,
                                            transpose=0, npass=100,
                                            method='a', dist='e')
        self.assertEqual(len(clusterid), len(data))

        correct = [0,1,1,2]
        mapping = [clusterid[correct.index(i)] for i in range(nclusters)]
        for i in range(len(clusterid)):
            self.assertEqual(clusterid[i], mapping[correct[i]])
      
        # Second data set
        weight = numpy.array([1,1])
        data = numpy.array([[ 1.1, 1.2 ],
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
        mask = numpy.array([[ 1, 1 ],
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
                            [ 1, 1 ]], int)

        clusterid, error, nfound = kcluster(data, nclusters=3, mask=mask,
                                            weight=weight, transpose=0,
                                            npass=100, method='a', dist='e')
        self.assertEqual(len(clusterid), len(data))

        correct = [0, 0, 0, 0, 0, 0, 1, 1, 2, 1, 1, 1, 1]
        mapping = [clusterid[correct.index(i)] for i in range(nclusters)]
        for i in range(len(clusterid)):
            self.assertEqual(clusterid[i], mapping[correct[i]])

    def test_clusterdistance(self):
        if TestCluster.module=='Bio.Cluster':
            from Bio.Cluster import clusterdistance
        elif TestCluster.module=='Pycluster':
            from Pycluster import clusterdistance

        # First data set
        weight = numpy.array([ 1,1,1,1,1 ])
        data   = numpy.array([[  1.1, 2.2, 3.3, 4.4, 5.5, ], 
                              [  3.1, 3.2, 1.3, 2.4, 1.5, ], 
                              [  4.1, 2.2, 0.3, 5.4, 0.5, ], 
                              [ 12.1, 2.0, 0.0, 5.0, 0.0, ]])
        mask   = numpy.array([[ 1, 1, 1, 1, 1], 
                              [ 1, 1, 1, 1, 1], 
                              [ 1, 1, 1, 1, 1], 
                              [ 1, 1, 1, 1, 1]], int)

        # Cluster assignments
        c1 = [0]
        c2 = [1,2]
        c3 = [3]

        distance = clusterdistance(data, mask=mask, weight=weight,
                                   index1=c1, index2=c2, dist='e',
                                   method='a', transpose=0);
        self.assertAlmostEqual(distance, 6.650, places=3)
        distance = clusterdistance(data, mask=mask, weight=weight,
                                   index1=c1, index2=c3, dist='e',
                                   method='a', transpose=0);
        self.assertAlmostEqual(distance, 32.508, places=3)
        distance = clusterdistance(data, mask=mask, weight=weight,
                                   index1=c2, index2=c3, dist='e',
                                   method='a', transpose=0);
        self.assertAlmostEqual(distance, 15.118, places=3)

        # Second data set
        weight =  numpy.array([ 1,1 ])
        data   =  numpy.array([[ 1.1, 1.2 ],
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
        mask = numpy.array([[ 1, 1 ],
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
                            [ 1, 1 ]], int)

        # Cluster assignments
        c1 = [ 0, 1, 2, 3 ]
        c2 = [ 4, 5, 6, 7 ]
        c3 = [ 8 ]

        distance = clusterdistance(data, mask=mask, weight=weight,
                                   index1=c1, index2=c2, dist='e',
                                   method='a', transpose=0);
        self.assertAlmostEqual(distance, 5.833, places=3)
        distance = clusterdistance(data, mask=mask, weight=weight,
                                   index1=c1, index2=c3, dist='e',
                                   method='a', transpose=0);
        self.assertAlmostEqual(distance, 3.298, places=3)
        distance = clusterdistance(data, mask=mask, weight=weight,
                                   index1=c2, index2=c3, dist='e',
                                   method='a', transpose=0);
        self.assertAlmostEqual(distance, 0.360, places=3)


    def test_treecluster(self):
        if TestCluster.module=='Bio.Cluster':
            from Bio.Cluster import treecluster
        elif TestCluster.module=='Pycluster':
            from Pycluster import treecluster

        # First data set
        weight1 =  [ 1,1,1,1,1 ]
        data1   =  numpy.array([[  1.1, 2.2, 3.3, 4.4, 5.5], 
                                [  3.1, 3.2, 1.3, 2.4, 1.5], 
                                [  4.1, 2.2, 0.3, 5.4, 0.5], 
                                [ 12.1, 2.0, 0.0, 5.0, 0.0]])
        mask1 = numpy.array([[ 1, 1, 1, 1, 1], 
                             [ 1, 1, 1, 1, 1], 
                             [ 1, 1, 1, 1, 1], 
                             [ 1, 1, 1, 1, 1]], int)
      
        # test first data set
        # Pairwise average-linkage clustering"
        tree = treecluster(data=data1, mask=mask1, weight=weight1,
                           transpose=0, method='a', dist='e')
        self.assertEqual(len(tree), len(data1)-1)
        self.assertEqual(tree[0].left, 2)
        self.assertEqual(tree[0].right, 1)
        self.assertAlmostEqual(tree[0].distance, 2.600, places=3)
        self.assertEqual(tree[1].left, -1)
        self.assertEqual(tree[1].right, 0)
        self.assertAlmostEqual(tree[1].distance, 7.300, places=3)
        self.assertEqual(tree[2].left, 3)
        self.assertEqual(tree[2].right, -2)
        self.assertAlmostEqual(tree[2].distance, 21.348, places=3)

        # Pairwise single-linkage clustering
        tree = treecluster(data=data1, mask=mask1, weight=weight1,
                           transpose=0, method='s', dist='e')
        self.assertEqual(len(tree), len(data1)-1)
        self.assertEqual(tree[0].left, 1)
        self.assertEqual(tree[0].right, 2)
        self.assertAlmostEqual(tree[0].distance, 2.600, places=3)
        self.assertEqual(tree[1].left, 0)
        self.assertEqual(tree[1].right, -1)
        self.assertAlmostEqual(tree[1].distance, 5.800, places=3)
        self.assertEqual(tree[2].left, -2)
        self.assertEqual(tree[2].right, 3)
        self.assertAlmostEqual(tree[2].distance, 12.908, places=3)

        # Pairwise centroid-linkage clustering
        tree = treecluster(data=data1, mask=mask1, weight=weight1,
                           transpose=0, method='c', dist='e')
        self.assertEqual(len(tree), len(data1)-1)
        self.assertEqual(tree[0].left, 1)
        self.assertEqual(tree[0].right, 2)
        self.assertAlmostEqual(tree[0].distance, 2.600, places=3)
        self.assertEqual(tree[1].left, 0)
        self.assertEqual(tree[1].right, -1)
        self.assertAlmostEqual(tree[1].distance, 6.650, places=3)
        self.assertEqual(tree[2].left, -2)
        self.assertEqual(tree[2].right, 3)
        self.assertAlmostEqual(tree[2].distance, 19.437, places=3)

        # Pairwise maximum-linkage clustering
        tree = treecluster(data=data1, mask=mask1, weight=weight1,
                           transpose=0, method='m', dist='e')
        self.assertEqual(len(tree), len(data1)-1)
        self.assertEqual(tree[0].left, 2)
        self.assertEqual(tree[0].right, 1)
        self.assertAlmostEqual(tree[0].distance, 2.600, places=3)
        self.assertEqual(tree[1].left, -1)
        self.assertEqual(tree[1].right, 0)
        self.assertAlmostEqual(tree[1].distance, 8.800, places=3)
        self.assertEqual(tree[2].left, 3)
        self.assertEqual(tree[2].right, -2)
        self.assertAlmostEqual(tree[2].distance, 32.508, places=3)
      
        # Second data set
        weight2 =  [ 1,1 ]
        data2 = numpy.array([[ 0.8223, 0.9295 ],
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
        mask2 = numpy.array([[ 1, 1 ],
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
                             [ 1, 1 ]], int)
      
        # Test second data set
        # Pairwise average-linkage clustering
        tree = treecluster(data=data2, mask=mask2, weight=weight2,
                           transpose=0, method='a', dist='e')
        self.assertEqual(len(tree), len(data2)-1)
        self.assertEqual(tree[0].left, 5)
        self.assertEqual(tree[0].right, 4)
        self.assertAlmostEqual(tree[0].distance, 0.003, places=3)
        self.assertEqual(tree[1].left, 9)
        self.assertEqual(tree[1].right, 12)
        self.assertAlmostEqual(tree[1].distance, 0.029, places=3)
        self.assertEqual(tree[2].left, 2)
        self.assertEqual(tree[2].right, 1)
        self.assertAlmostEqual(tree[2].distance, 0.061, places=3)
        self.assertEqual(tree[3].left, 11)
        self.assertEqual(tree[3].right, -2)
        self.assertAlmostEqual(tree[3].distance, 0.070, places=3)
        self.assertEqual(tree[4].left, -4)
        self.assertEqual(tree[4].right, 10)
        self.assertAlmostEqual(tree[4].distance, 0.128, places=3)
        self.assertEqual(tree[5].left, 7)
        self.assertEqual(tree[5].right, -5)
        self.assertAlmostEqual(tree[5].distance, 0.224, places=3)
        self.assertEqual(tree[6].left, -3)
        self.assertEqual(tree[6].right, 0)
        self.assertAlmostEqual(tree[6].distance, 0.254, places=3)
        self.assertEqual(tree[7].left, -1)
        self.assertEqual(tree[7].right, 3)
        self.assertAlmostEqual(tree[7].distance, 0.391, places=3)
        self.assertEqual(tree[8].left, -8)
        self.assertEqual(tree[8].right, -7)
        self.assertAlmostEqual(tree[8].distance, 0.532, places=3)
        self.assertEqual(tree[9].left, 8)
        self.assertEqual(tree[9].right, -9)
        self.assertAlmostEqual(tree[9].distance, 3.234, places=3)
        self.assertEqual(tree[10].left, -6)
        self.assertEqual(tree[10].right, 6)
        self.assertAlmostEqual(tree[10].distance, 4.636, places=3)
        self.assertEqual(tree[11].left, -11)
        self.assertEqual(tree[11].right, -10)
        self.assertAlmostEqual(tree[11].distance, 12.741, places=3)
      
        # Pairwise single-linkage clustering
        tree = treecluster(data=data2, mask=mask2, weight=weight2,
                           transpose=0, method='s', dist='e')
        self.assertEqual(len(tree), len(data2)-1)
        self.assertEqual(tree[0].left, 4)
        self.assertEqual(tree[0].right, 5)
        self.assertAlmostEqual(tree[0].distance, 0.003, places=3)
        self.assertEqual(tree[1].left, 9)
        self.assertEqual(tree[1].right, 12)
        self.assertAlmostEqual(tree[1].distance, 0.029, places=3)
        self.assertEqual(tree[2].left, 11)
        self.assertEqual(tree[2].right, -2)
        self.assertAlmostEqual(tree[2].distance, 0.033, places=3)
        self.assertEqual(tree[3].left, 1)
        self.assertEqual(tree[3].right, 2)
        self.assertAlmostEqual(tree[3].distance, 0.061, places=3)
        self.assertEqual(tree[4].left, 10)
        self.assertEqual(tree[4].right, -3)
        self.assertAlmostEqual(tree[4].distance, 0.077, places=3)
        self.assertEqual(tree[5].left, 7)
        self.assertEqual(tree[5].right, -5)
        self.assertAlmostEqual(tree[5].distance, 0.092, places=3)
        self.assertEqual(tree[6].left, 0)
        self.assertEqual(tree[6].right, -4)
        self.assertAlmostEqual(tree[6].distance, 0.242, places=3)
        self.assertEqual(tree[7].left, -7)
        self.assertEqual(tree[7].right, -1)
        self.assertAlmostEqual(tree[7].distance, 0.246, places=3)
        self.assertEqual(tree[8].left, 3)
        self.assertEqual(tree[8].right, -8)
        self.assertAlmostEqual(tree[8].distance, 0.287, places=3)
        self.assertEqual(tree[9].left, -9)
        self.assertEqual(tree[9].right, 8)
        self.assertAlmostEqual(tree[9].distance, 1.936, places=3)
        self.assertEqual(tree[10].left, -10)
        self.assertEqual(tree[10].right, -6)
        self.assertAlmostEqual(tree[10].distance, 3.432, places=3)
        self.assertEqual(tree[11].left, 6)
        self.assertEqual(tree[11].right, -11)
        self.assertAlmostEqual(tree[11].distance, 3.535, places=3)
      
        # Pairwise centroid-linkage clustering
        tree = treecluster(data=data2, mask=mask2, weight=weight2,
                           transpose=0, method='c', dist='e')
        self.assertEqual(len(tree), len(data2)-1)
        self.assertEqual(tree[0].left, 4)
        self.assertEqual(tree[0].right, 5)
        self.assertAlmostEqual(tree[0].distance, 0.003, places=3)
        self.assertEqual(tree[1].left, 12)
        self.assertEqual(tree[1].right, 9)
        self.assertAlmostEqual(tree[1].distance, 0.029, places=3)
        self.assertEqual(tree[2].left, 1)
        self.assertEqual(tree[2].right, 2)
        self.assertAlmostEqual(tree[2].distance, 0.061, places=3)
        self.assertEqual(tree[3].left, -2)
        self.assertEqual(tree[3].right, 11)
        self.assertAlmostEqual(tree[3].distance, 0.063, places=3)
        self.assertEqual(tree[4].left, 10)
        self.assertEqual(tree[4].right, -4)
        self.assertAlmostEqual(tree[4].distance, 0.109, places=3)
        self.assertEqual(tree[5].left, -5)
        self.assertEqual(tree[5].right, 7)
        self.assertAlmostEqual(tree[5].distance, 0.189, places=3)
        self.assertEqual(tree[6].left, 0)
        self.assertEqual(tree[6].right, -3)
        self.assertAlmostEqual(tree[6].distance, 0.239, places=3)
        self.assertEqual(tree[7].left, 3)
        self.assertEqual(tree[7].right, -1)
        self.assertAlmostEqual(tree[7].distance, 0.390, places=3)
        self.assertEqual(tree[8].left, -7)
        self.assertEqual(tree[8].right, -8)
        self.assertAlmostEqual(tree[8].distance, 0.382, places=3)
        self.assertEqual(tree[9].left, -9)
        self.assertEqual(tree[9].right, 8)
        self.assertAlmostEqual(tree[9].distance, 3.063, places=3)
        self.assertEqual(tree[10].left, 6)
        self.assertEqual(tree[10].right, -6)
        self.assertAlmostEqual(tree[10].distance, 4.578, places=3)
        self.assertEqual(tree[11].left, -10)
        self.assertEqual(tree[11].right, -11)
        self.assertAlmostEqual(tree[11].distance, 11.536, places=3)
      
        # Pairwise maximum-linkage clustering
        tree = treecluster(data=data2, mask=mask2, weight=weight2,
                           transpose=0, method='m', dist='e')
        self.assertEqual(len(tree), len(data2)-1)
        self.assertEqual(tree[0].left, 5)
        self.assertEqual(tree[0].right, 4)
        self.assertAlmostEqual(tree[0].distance, 0.003, places=3)
        self.assertEqual(tree[1].left, 9)
        self.assertEqual(tree[1].right, 12)
        self.assertAlmostEqual(tree[1].distance, 0.029, places=3)
        self.assertEqual(tree[2].left, 2)
        self.assertEqual(tree[2].right, 1)
        self.assertAlmostEqual(tree[2].distance, 0.061, places=3)
        self.assertEqual(tree[3].left, 11)
        self.assertEqual(tree[3].right, 10)
        self.assertAlmostEqual(tree[3].distance, 0.077, places=3)
        self.assertEqual(tree[4].left, -2)
        self.assertEqual(tree[4].right, -4)
        self.assertAlmostEqual(tree[4].distance, 0.216, places=3)
        self.assertEqual(tree[5].left, -3)
        self.assertEqual(tree[5].right, 0)
        self.assertAlmostEqual(tree[5].distance, 0.266, places=3)
        self.assertEqual(tree[6].left, -5)
        self.assertEqual(tree[6].right, 7)
        self.assertAlmostEqual(tree[6].distance, 0.302, places=3)
        self.assertEqual(tree[7].left, -1)
        self.assertEqual(tree[7].right, 3)
        self.assertAlmostEqual(tree[7].distance, 0.425, places=3)
        self.assertEqual(tree[8].left, -8)
        self.assertEqual(tree[8].right, -6)
        self.assertAlmostEqual(tree[8].distance, 0.968, places=3)
        self.assertEqual(tree[9].left, 8)
        self.assertEqual(tree[9].right, 6)
        self.assertAlmostEqual(tree[9].distance, 3.975, places=3)
        self.assertEqual(tree[10].left, -10)
        self.assertEqual(tree[10].right, -7)
        self.assertAlmostEqual(tree[10].distance, 5.755, places=3)
        self.assertEqual(tree[11].left, -11)
        self.assertEqual(tree[11].right, -9)
        self.assertAlmostEqual(tree[11].distance, 22.734, places=3)

    def test_somcluster(self):
        if TestCluster.module=='Bio.Cluster':
            from Bio.Cluster import somcluster
        elif TestCluster.module=='Pycluster':
            from Pycluster import somcluster

        # First data set
        weight = [ 1,1,1,1,1 ]
        data = numpy.array([[  1.1, 2.2, 3.3, 4.4, 5.5], 
                            [  3.1, 3.2, 1.3, 2.4, 1.5], 
                            [  4.1, 2.2, 0.3, 5.4, 0.5], 
                            [ 12.1, 2.0, 0.0, 5.0, 0.0]])
        mask = numpy.array([[ 1, 1, 1, 1, 1], 
                            [ 1, 1, 1, 1, 1], 
                            [ 1, 1, 1, 1, 1], 
                            [ 1, 1, 1, 1, 1]], int)

        clusterid, celldata = somcluster(data=data, mask=mask, weight=weight,
                                         transpose=0, nxgrid=10, nygrid=10,
                                         inittau=0.02, niter=100, dist='e')
        self.assertEqual(len(clusterid), len(data))
        self.assertEqual(len(clusterid[0]), 2)

        # Second data set
        weight =  [ 1,1 ]
        data = numpy.array([[ 1.1, 1.2 ],
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
        mask = numpy.array([[ 1, 1 ],
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
                            [ 1, 1 ]], int)

        clusterid, celldata = somcluster(data=data, mask=mask, weight=weight,
                                         transpose=0, nxgrid=10, nygrid=10,
                                         inittau=0.02, niter=100, dist='e')
        self.assertEqual(len(clusterid), len(data))
        self.assertEqual(len(clusterid[0]), 2)

    def test_distancematrix_kmedoids(self):
        if TestCluster.module=='Bio.Cluster':
            from Bio.Cluster import distancematrix, kmedoids
        elif TestCluster.module=='Pycluster':
            from Pycluster import distancematrix, kmedoids

        data = numpy.array([[2.2, 3.3, 4.4],
                            [2.1, 1.4, 5.6],
                            [7.8, 9.0, 1.2],
                            [4.5, 2.3, 1.5],
                            [4.2, 2.4, 1.9],
                            [3.6, 3.1, 9.3],
                            [2.3, 1.2, 3.9],
                            [4.2, 9.6, 9.3],
                            [1.7, 8.9, 1.1]])
        mask = numpy.array([[1, 1, 1],
                            [1, 1, 1],
                            [0, 1, 1],
                            [1, 1, 1],
                            [1, 1, 1],
                            [0, 1, 0],
                            [1, 1, 1],
                            [1, 0, 1],
                            [1, 1, 1]], int)
        weight = numpy.array([2.0, 1.0, 0.5])
        matrix = distancematrix(data, mask=mask, weight=weight)

        self.assertAlmostEqual(matrix[1][0], 1.243, places=3)

        self.assertAlmostEqual(matrix[2][0], 25.073, places=3)
        self.assertAlmostEqual(matrix[2][1], 44.960, places=3)

        self.assertAlmostEqual(matrix[3][0], 4.510, places=3)
        self.assertAlmostEqual(matrix[3][1], 5.924, places=3)
        self.assertAlmostEqual(matrix[3][2], 29.957, places=3)

        self.assertAlmostEqual(matrix[4][0], 3.410, places=3)
        self.assertAlmostEqual(matrix[4][1], 4.761, places=3)
        self.assertAlmostEqual(matrix[4][2], 29.203, places=3)
        self.assertAlmostEqual(matrix[4][3], 0.077, places=3)

        self.assertAlmostEqual(matrix[5][0], 0.040, places=3)
        self.assertAlmostEqual(matrix[5][1], 2.890, places=3)
        self.assertAlmostEqual(matrix[5][2], 34.810, places=3)
        self.assertAlmostEqual(matrix[5][3], 0.640, places=3)
        self.assertAlmostEqual(matrix[5][4], 0.490, places=3)

        self.assertAlmostEqual(matrix[6][0], 1.301, places=3)
        self.assertAlmostEqual(matrix[6][1], 0.447, places=3)
        self.assertAlmostEqual(matrix[6][2], 42.990, places=3)
        self.assertAlmostEqual(matrix[6][3], 3.934, places=3)
        self.assertAlmostEqual(matrix[6][4], 3.046, places=3)
        self.assertAlmostEqual(matrix[6][5], 3.610, places=3)

        self.assertAlmostEqual(matrix[7][0], 8.002, places=3)
        self.assertAlmostEqual(matrix[7][1], 6.266, places=3)
        self.assertAlmostEqual(matrix[7][2], 65.610, places=3)
        self.assertAlmostEqual(matrix[7][3], 12.240, places=3)
        self.assertAlmostEqual(matrix[7][4], 10.952, places=3)
        self.assertAlmostEqual(matrix[7][5], 0.000, places=3)
        self.assertAlmostEqual(matrix[7][6], 8.720, places=3)

        self.assertAlmostEqual(matrix[8][0], 10.659, places=3)
        self.assertAlmostEqual(matrix[8][1], 19.056, places=3)
        self.assertAlmostEqual(matrix[8][2], 0.010, places=3)
        self.assertAlmostEqual(matrix[8][3], 16.949, places=3)
        self.assertAlmostEqual(matrix[8][4], 15.734, places=3)
        self.assertAlmostEqual(matrix[8][5], 33.640, places=3)
        self.assertAlmostEqual(matrix[8][6], 18.266, places=3)
        self.assertAlmostEqual(matrix[8][7], 18.448, places=3)
        clusterid, error, nfound = kmedoids(matrix, npass=1000)
        self.assertEqual(clusterid[0], 5)
        self.assertEqual(clusterid[1], 5)
        self.assertEqual(clusterid[2], 2)
        self.assertEqual(clusterid[3], 5)
        self.assertEqual(clusterid[4], 5)
        self.assertEqual(clusterid[5], 5)
        self.assertEqual(clusterid[6], 5)
        self.assertEqual(clusterid[7], 5)
        self.assertEqual(clusterid[8], 2)
        self.assertAlmostEqual(error, 7.680, places=3)

    def test_pca(self):
        if TestCluster.module=='Bio.Cluster':
            from Bio.Cluster import pca
        elif TestCluster.module=='Pycluster':
            from Pycluster import pca

        data = numpy.array([[ 3.1, 1.2 ],
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
                            [ 5.1, 5.2 ],
                           ])

        mean, coordinates, pc, eigenvalues =  pca(data)
        self.assertAlmostEqual(mean[0], 3.5461538461538464)
        self.assertAlmostEqual(mean[1], 3.5307692307692311)
        self.assertAlmostEqual(coordinates[0,0],  2.0323189722653883)
        self.assertAlmostEqual(coordinates[0,1],  1.2252420399694917)
        self.assertAlmostEqual(coordinates[1,0],  3.0936985166252251)
        self.assertAlmostEqual(coordinates[1,1], -0.10647619705157851)
        self.assertAlmostEqual(coordinates[2,0],  3.1453186907749426)
        self.assertAlmostEqual(coordinates[2,1], -0.46331699855941139)
        self.assertAlmostEqual(coordinates[3,0],  2.5440202962223761)
        self.assertAlmostEqual(coordinates[3,1],  0.20633980959571077)
        self.assertAlmostEqual(coordinates[4,0],  2.4468278463376221)
        self.assertAlmostEqual(coordinates[4,1], -0.28412285736824866)
        self.assertAlmostEqual(coordinates[5,0],  2.4468278463376221)
        self.assertAlmostEqual(coordinates[5,1], -0.28412285736824866)
        self.assertAlmostEqual(coordinates[6,0], -3.2018619434743254)
        self.assertAlmostEqual(coordinates[6,1],  0.019692314198662915)
        self.assertAlmostEqual(coordinates[7,0], -3.2018619434743254)
        self.assertAlmostEqual(coordinates[7,1],  0.019692314198662915)
        self.assertAlmostEqual(coordinates[8,0],  0.46978641990344067)
        self.assertAlmostEqual(coordinates[8,1], -0.17778754731982949)
        self.assertAlmostEqual(coordinates[9,0], -2.5549912731867215)
        self.assertAlmostEqual(coordinates[9,1],  0.19733897451533403)
        self.assertAlmostEqual(coordinates[10,0], -2.5033710990370044)
        self.assertAlmostEqual(coordinates[10,1], -0.15950182699250004)
        self.assertAlmostEqual(coordinates[11,0], -2.4365601663089413)
        self.assertAlmostEqual(coordinates[11,1], -0.23390813900973562)
        self.assertAlmostEqual(coordinates[12,0], -2.2801521629852974)
        self.assertAlmostEqual(coordinates[12,1],  0.0409309711916888)
        self.assertAlmostEqual(pc[0,0], -0.66810932728062988)
        self.assertAlmostEqual(pc[0,1], -0.74406312017235743)
        self.assertAlmostEqual(pc[1,0],  0.74406312017235743)
        self.assertAlmostEqual(pc[1,1], -0.66810932728062988)
        self.assertAlmostEqual(eigenvalues[0], 9.3110471246032844)
        self.assertAlmostEqual(eigenvalues[1], 1.4437456297481428)

        data = numpy.array([[ 2.3, 4.5, 1.2, 6.7, 5.3, 7.1],
                            [ 1.3, 6.5, 2.2, 5.7, 6.2, 9.1],
                            [ 3.2, 7.2, 3.2, 7.4, 7.3, 8.9],
                            [ 4.2, 5.2, 9.2, 4.4, 6.3, 7.2]])
        mean, coordinates, pc, eigenvalues =  pca(data)
        self.assertAlmostEqual(mean[0], 2.7500)
        self.assertAlmostEqual(mean[1], 5.8500)
        self.assertAlmostEqual(mean[2], 3.9500)
        self.assertAlmostEqual(mean[3], 6.0500)
        self.assertAlmostEqual(mean[4], 6.2750)
        self.assertAlmostEqual(mean[5], 8.0750)
        self.assertAlmostEqual(coordinates[0,0],  2.6460846688406905)
        self.assertAlmostEqual(coordinates[0,1], -2.1421701432732418)
        self.assertAlmostEqual(coordinates[0,2], -0.56620932754145858)
        self.assertAlmostEqual(coordinates[0,3],  0.0)
        self.assertAlmostEqual(coordinates[1,0],  2.0644120899917544)
        self.assertAlmostEqual(coordinates[1,1],  0.55542108669180323)
        self.assertAlmostEqual(coordinates[1,2],  1.4818772348457117)
        self.assertAlmostEqual(coordinates[1,3],  0.0)
        self.assertAlmostEqual(coordinates[2,0],  1.0686641862092987)
        self.assertAlmostEqual(coordinates[2,1],  1.9994412069101073)
        self.assertAlmostEqual(coordinates[2,2], -1.000720598980291)
        self.assertAlmostEqual(coordinates[2,3],  0.0)
        self.assertAlmostEqual(coordinates[3,0], -5.77916094504174)
        self.assertAlmostEqual(coordinates[3,1], -0.41269215032867046)
        self.assertAlmostEqual(coordinates[3,2],  0.085052691676038017)
        self.assertAlmostEqual(coordinates[3,3],  0.0)
        self.assertAlmostEqual(pc[0,0], -0.26379660005997291)
        self.assertAlmostEqual(pc[0,1],  0.064814972617134495)
        self.assertAlmostEqual(pc[0,2], -0.91763310094893846)
        self.assertAlmostEqual(pc[0,3],  0.26145408875373249)
        self.assertAlmostEqual(pc[1,0],  0.05073770520434398)
        self.assertAlmostEqual(pc[1,1],  0.68616983388698793)
        self.assertAlmostEqual(pc[1,2],  0.13819106187213354)
        self.assertAlmostEqual(pc[1,3],  0.19782544121828985)
        self.assertAlmostEqual(pc[2,0], -0.63000893660095947)
        self.assertAlmostEqual(pc[2,1],  0.091155993862151397)
        self.assertAlmostEqual(pc[2,2],  0.045630391256086845)
        self.assertAlmostEqual(pc[2,3], -0.67456694780914772)
        # As the last eigenvalue is zero, the corresponding eigenvector is
        # strongly affected by roundoff error, and is not being tested here.
        # For PCA, this doesn't matter since all data have a zero coefficient
        # along this eigenvector.
        self.assertAlmostEqual(eigenvalues[0], 6.7678878332578778)
        self.assertAlmostEqual(eigenvalues[1], 3.0108911400291856)
        self.assertAlmostEqual(eigenvalues[2], 1.8775592718563467)
        self.assertAlmostEqual(eigenvalues[3], 0.0)

if __name__ == "__main__":
    TestCluster.module = 'Bio.Cluster'
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
