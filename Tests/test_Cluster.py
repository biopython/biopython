# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for Cluster module."""

import unittest

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.Cluster.") from None


class TestCluster(unittest.TestCase):

    module = "Bio.Cluster"

    def test_matrix_parse(self):
        if TestCluster.module == "Bio.Cluster":
            from Bio.Cluster import treecluster
        elif TestCluster.module == "Pycluster":
            from Pycluster import treecluster

        # Normal matrix, no errors
        data1 = numpy.array([[1.1, 1.2],
                             [1.4, 1.3],
                             [1.1, 1.5],
                             [2.0, 1.5],
                             [1.7, 1.9],
                             [1.7, 1.9],
                             [5.7, 5.9],
                             [5.7, 5.9],
                             [3.1, 3.3],
                             [5.4, 5.3],
                             [5.1, 5.5],
                             [5.0, 5.5],
                             [5.1, 5.2]])

        # Another normal matrix, no errors; written as a list
        data2 = [[1.1, 2.2, 3.3, 4.4, 5.5],
                 [3.1, 3.2, 1.3, 2.4, 1.5],
                 [4.1, 2.2, 0.3, 5.4, 0.5],
                 [2.1, 2.0, 0.0, 5.0, 0.0]]

        # Rows are not contiguous
        data3 = data1[::2, :]

        # Columns are not contiguous
        data4 = numpy.array(data2)[:, ::2]

        # Matrix using float32
        data5 = numpy.array([[1.1, 2.2, 3.3, 4.4, 5.5],
                             [3.1, 3.2, 1.3, 2.4, 1.5],
                             [4.1, 2.2, 0.3, 5.4, 0.5],
                             [2.1, 2.0, 0.0, 5.0, 0.0]], numpy.float32)

        # Matrix using int
        data6 = numpy.array([[1, 2, 3, 4, 5],
                             [3, 3, 1, 2, 1],
                             [4, 2, 0, 5, 0],
                             [2, 2, 0, 5, 0]], numpy.int32)
        try:
            treecluster(data1)
        except Exception:
            self.fail("treecluster failed to accept matrix data1")

        try:
            treecluster(data2)
        except Exception:
            self.fail("treecluster failed to accept matrix data2")

        try:
            treecluster(data3)
        except Exception:
            self.fail("treecluster failed to accept matrix data3")

        try:
            treecluster(data4)
        except Exception:
            self.fail("treecluster failed to accept matrix data4")

        try:
            treecluster(data5)
        except Exception:
            self.fail("treecluster failed to accept matrix data5")

        try:
            treecluster(data6)
        except Exception:
            self.fail("treecluster failed to accept matrix data6")

        # Ragged matrix
        data7 = [[91.1, 92.2, 93.3, 94.4, 95.5],
                 [93.1, 93.2, 91.3, 92.4],
                 [94.1, 92.2, 90.3],
                 [12.1, 92.0, 90.0, 95.0, 90.0]]

        # Matrix with bad cells
        data8 = [[7.1, 7.2, 7.3, 7.4, 7.5],
                 [7.1, 7.2, 7.3, 7.4, "snoopy"],
                 [7.1, 7.2, 7.3, None, None]]

        # Matrix with a bad row
        data9 = [[23.1, 23.2, 23.3, 23.4, 23.5],
                 None,
                 [23.1, 23.0, 23.0, 23.0, 23.0]]

        # Various references that don't point to matrices at all
        data10 = "snoopy"
        data11 = {"a": [[2.3, 1.2], [3.3, 5.6]]}
        data12 = []
        data13 = [None]

        # Array of incorrect rank
        data14 = numpy.array([[[1.1, 1.2], [2.3, 1.2], [3.4, 1.6]],
                              [[1.4, 1.3], [3.2, 4.5], [9.8, 4.9]],
                              [[1.1, 1.5], [1.1, 2.3], [6.5, 0.4]]])

        # Array with non-numerical values
        data15 = numpy.array([["a", "b", "c"],
                              ["e", "f", "g"]], "c")

        # Empty array
        data16 = numpy.array([[]], "d")

        self.assertRaises(ValueError, treecluster, data7)
        self.assertRaises(ValueError, treecluster, data8)
        self.assertRaises(ValueError, treecluster, data9)
        self.assertRaises(ValueError, treecluster, data10)
        self.assertRaises(TypeError, treecluster, data11)
        self.assertRaises(ValueError, treecluster, data12)
        self.assertRaises(ValueError, treecluster, data13)
        self.assertRaises(ValueError, treecluster, data14)
        self.assertRaises(ValueError, treecluster, data15)
        self.assertRaises(ValueError, treecluster, data16)

    def test_mask_parse(self):
        if TestCluster.module == "Bio.Cluster":
            from Bio.Cluster import treecluster
        elif TestCluster.module == "Pycluster":
            from Pycluster import treecluster

        # data matrix
        data = numpy.array([[1.1, 2.2, 3.3, 4.4, 5.5],
                            [3.1, 3.2, 1.3, 2.4, 1.5],
                            [4.1, 2.2, 0.3, 5.4, 0.5],
                            [2.1, 2.0, 0.0, 5.0, 0.0]])

        # Normal mask, no errors
        mask1 = numpy.array([[1, 1, 0, 1, 0],
                             [1, 1, 1, 0, 0],
                             [1, 1, 0, 1, 1],
                             [1, 0, 1, 1, 0]])

        # Same mask, no errors; written as a list
        mask2 = [[1, 1, 0, 1, 0],
                 [1, 1, 1, 0, 0],
                 [1, 1, 0, 1, 1],
                 [1, 0, 1, 1, 0]]

        # Rows are not contiguous
        mask3 = numpy.array([[1, 1, 0, 1, 0],
                             [1, 1, 1, 0, 0],
                             [1, 1, 1, 0, 0],
                             [1, 1, 0, 1, 1],
                             [1, 1, 1, 0, 0],
                             [1, 1, 0, 1, 1],
                             [1, 1, 0, 1, 1],
                             [1, 0, 1, 1, 0]])
        mask3 = mask3[::2, :]

        # Columns are not contiguous
        mask4 = numpy.array([[1, 1, 0, 1, 0, 1, 0, 0, 1, 1],
                             [1, 1, 1, 0, 0, 1, 1, 0, 0, 1],
                             [1, 1, 0, 1, 1, 1, 0, 1, 1, 0],
                             [1, 0, 1, 1, 0, 1, 0, 0, 1, 1]])
        mask4 = mask4[:, ::2]

        # Matrix using int16
        mask5 = numpy.array([[1, 1, 0, 1, 0],
                             [1, 1, 1, 0, 0],
                             [1, 1, 1, 0, 0],
                             [1, 1, 0, 1, 1]], numpy.int16)

        # Matrix using float
        mask6 = numpy.array([[1.0, 2.2, 3.1, 4.8, 5.1],
                             [3.3, 3.3, 1.4, 2.4, 1.2],
                             [4.1, 2.2, 0.6, 5.5, 0.6],
                             [2.7, 2.5, 0.4, 5.7, 0.2]], numpy.float)
        try:
            treecluster(data, mask1)
        except Exception:
            self.fail("treecluster failed to accept matrix mask1")

        try:
            treecluster(data, mask2)
        except Exception:
            self.fail("treecluster failed to accept matrix mask2")

        try:
            treecluster(data, mask3)
        except Exception:
            self.fail("treecluster failed to accept matrix mask3")

        try:
            treecluster(data, mask4)
        except Exception:
            self.fail("treecluster failed to accept matrix mask4")

        try:
            treecluster(data, mask5)
        except Exception:
            self.fail("treecluster failed to accept matrix mask5")

        try:
            treecluster(data, mask6)
        except Exception:
            self.fail("treecluster failed to accept matrix mask6")

        # Ragged mask
        mask7 = [[1, 1, 0, 1],
                 [1, 1, 1, 0, 0],
                 [1, 1, 0, 1, 1],
                 [1, 1, 0]]

        # Mask with incorrect number of rows
        mask8 = numpy.array([[1, 1, 0, 1, 0],
                             [1, 1, 1, 0, 0],
                             [1, 1, 0, 1, 1],
                             [0, 1, 1, 0, 1],
                             [1, 0, 1, 1, 0]])

        # Mask with incorrect number of columns
        mask9 = numpy.array([[1, 1, 0, 1, 0, 1],
                             [1, 1, 1, 0, 0, 0],
                             [0, 1, 1, 0, 1, 1],
                             [1, 0, 1, 1, 0, 1]])

        # Matrix with bad cells
        mask10 = [[1, 1, 0, 1, 0],
                  [1, 1, 1, 0, "snoopy"],
                  [1, 1, 0, 1, 1],
                  [1, 0, 1, 1, 0]]

        # Matrix with a bad row
        mask11 = [[1, 1, 0, 1, 0],
                  None,
                  [1, 1, 0, 1, 1],
                  [1, 0, 1, 1, 0]]

        # Array with non-numerical values
        mask12 = numpy.array([["a", "b", "c"],
                              ["e", "f", "g"]], "c")

        # Empty arrays
        mask13 = numpy.array([[]], "d")
        mask14 = []

        # Array of incorrect rank
        mask15 = numpy.array([[[1, 1], [0, 1], [1, 1]],
                              [[1, 1], [0, 1], [1, 1]],
                              [[1, 1], [1, 1], [1, 0]]])

        # References that cannot be converted to a matrix of int
        mask16 = "snoopy"
        mask17 = {"a": [[1, 0], [1, 1]]}
        mask18 = [None]

        self.assertRaises(ValueError, treecluster, data, mask7)
        self.assertRaises(ValueError, treecluster, data, mask8)
        self.assertRaises(ValueError, treecluster, data, mask9)
        self.assertRaises(ValueError, treecluster, data, mask10)
        self.assertRaises(ValueError, treecluster, data, mask11)
        self.assertRaises(ValueError, treecluster, data, mask12)
        self.assertRaises(ValueError, treecluster, data, mask13)
        self.assertRaises(ValueError, treecluster, data, mask14)
        self.assertRaises(ValueError, treecluster, data, mask15)
        self.assertRaises(ValueError, treecluster, data, mask16)
        self.assertRaises(TypeError, treecluster, data, mask17)
        self.assertRaises(TypeError, treecluster, data, mask18)

    def test_kcluster_arguments(self):
        # Test if incorrect arguments are caught by the C code
        if TestCluster.module == "Bio.Cluster":
            from Bio.Cluster._cluster import kcluster, clustercentroids
        elif TestCluster.module == "Pycluster":
            from Pycluster._cluster import kcluster, clustercentroids

        nclusters = 3
        weight = numpy.array([1.0, 1.0, 1.0, 1.0, 1.0])
        data = numpy.array([[1.1, 2.2, 3.3, 4.4, 5.5],
                            [3.1, 3.2, 1.3, 2.4, 1.5],
                            [4.1, 2.2, 0.3, 5.4, 0.5],
                            [9.9, 2.0, 0.0, 5.0, 0.0]])
        mask = numpy.array([[1, 1, 1, 1, 1],
                            [1, 1, 1, 1, 1],
                            [1, 1, 1, 1, 1],
                            [1, 1, 1, 1, 1]], numpy.int32)
        nrows, ncols = data.shape
        clusterid = numpy.zeros(nrows, numpy.int32)

        message = "^data matrix is empty$"
        with self.assertRaisesRegex(ValueError, message):
            kcluster(data[:0, :], nclusters=nclusters,
                     mask=mask, weight=weight,
                     transpose=False, npass=100, method="a", dist="e")
        message = "^mask has incorrect rank 1 \\(expected 2\\)$"
        with self.assertRaisesRegex(ValueError, message):
            kcluster(data, nclusters=nclusters,
                     mask=numpy.zeros(3), weight=weight,
                     transpose=False, npass=100, method="a", dist="e")
        message = "^mask has incorrect dimensions 4 x 3 \\(expected 4 x 5\\)$"
        with self.assertRaisesRegex(ValueError, message):
            kcluster(data, nclusters=nclusters,
                     mask=numpy.zeros((4, 3), numpy.int32), weight=weight,
                     transpose=False, npass=100, method="a", dist="e",
                     clusterid=clusterid)
        message = "^incorrect rank 2 \\(expected 1\\)$"
        with self.assertRaisesRegex(ValueError, message):
            kcluster(data, nclusters=nclusters,
                     mask=mask, weight=numpy.zeros((2, 2)),
                     transpose=False, npass=100, method="a", dist="e")
        message = "^weight has incorrect size 3 \\(expected 5\\)$"
        with self.assertRaisesRegex(ValueError, message):
            kcluster(data, nclusters=nclusters,
                     mask=mask, weight=numpy.zeros(3),
                     transpose=False, npass=100, method="a", dist="e",
                     clusterid=clusterid)
        message = "^nclusters should be positive$"
        with self.assertRaisesRegex(ValueError, message):
            kcluster(data, nclusters=-1, mask=mask, weight=weight,
                     transpose=False, npass=100, method="a", dist="e",
                     clusterid=clusterid)
        message = "^more clusters than items to be clustered$"
        with self.assertRaisesRegex(ValueError, message):
            kcluster(data, nclusters=1234, mask=mask, weight=weight,
                     transpose=False, npass=100, method="a", dist="e",
                     clusterid=clusterid)
        message = "^incorrect size \\(3, expected 4\\)$"
        with self.assertRaisesRegex(ValueError, message):
            kcluster(data, nclusters=nclusters, mask=mask, weight=weight,
                     transpose=False, npass=0, method="a", dist="e",
                     clusterid=clusterid[:3])
        message = "^more clusters requested than found in clusterid$"
        with self.assertRaisesRegex(ValueError, message):
            kcluster(data, nclusters=nclusters, mask=mask, weight=weight,
                     transpose=False, npass=0, method="a", dist="e",
                     clusterid=clusterid)
        clusterid = numpy.array([0, -1, 2, 3], numpy.int32)
        message = "^negative cluster number found$"
        with self.assertRaisesRegex(ValueError, message):
            kcluster(data, nclusters=nclusters, mask=mask, weight=weight,
                     transpose=False, npass=0, method="a", dist="e",
                     clusterid=clusterid)
        clusterid = numpy.array([0, 0, 2, 3], numpy.int32)
        message = "^cluster 1 is empty$"
        with self.assertRaisesRegex(ValueError, message):
            kcluster(data, nclusters=nclusters, mask=mask, weight=weight,
                     transpose=False, npass=0, method="a", dist="e",
                     clusterid=clusterid)

    def test_kcluster(self):
        if TestCluster.module == "Bio.Cluster":
            from Bio.Cluster import kcluster, clustercentroids
        elif TestCluster.module == "Pycluster":
            from Pycluster import kcluster, clustercentroids

        nclusters = 3

        # First data set
        weight = numpy.array([1, 1, 1, 1, 1])
        data = numpy.array([[1.1, 2.2, 3.3, 4.4, 5.5],
                            [3.1, 3.2, 1.3, 2.4, 1.5],
                            [4.1, 2.2, 0.3, 5.4, 0.5],
                            [9.9, 2.0, 0.0, 5.0, 0.0]])
        mask = numpy.array([[1, 1, 1, 1, 1],
                            [1, 1, 1, 1, 1],
                            [1, 1, 1, 1, 1],
                            [1, 1, 1, 1, 1]], int)
        nrows, ncols = data.shape

        clusterid, error, nfound = kcluster(data, nclusters=nclusters,
                                            mask=mask, weight=weight,
                                            transpose=False, npass=100,
                                            method="a", dist="e")
        self.assertEqual(len(clusterid), len(data))

        correct = [0, 1, 1, 2]
        mapping = [clusterid[correct.index(i)] for i in range(nclusters)]
        for i in range(len(clusterid)):
            self.assertEqual(clusterid[i], mapping[correct[i]])

        cdata, cmask = clustercentroids(data, mask=mask, clusterid=clusterid,
                                        method="a", transpose=False)

        self.assertEqual(cdata.shape, (nclusters, ncols))
        self.assertEqual(cmask.shape, (nclusters, ncols))
        for value in cmask.flat:
            self.assertEqual(value, 1)

        correct = numpy.array([[1.1, 2.2, 3.3, 4.4, 5.5],
                               [3.6, 2.7, 0.8, 3.9, 1.0],
                               [9.9, 2.0, 0.0, 5.0, 0.0]])
        for i in range(nclusters):
            for j in range(ncols):
                self.assertAlmostEqual(cdata[mapping[i], j], correct[i, j])

        # First data set, using transpose=True
        weight = numpy.array([1, 1, 1, 1])
        clusterid, error, nfound = kcluster(data, nclusters=nclusters,
                                            mask=mask, weight=weight,
                                            transpose=True, npass=100,
                                            method="a", dist="e")
        self.assertEqual(len(clusterid), ncols)

        correct = [0, 1, 1, 2, 1]
        mapping = [clusterid[correct.index(i)] for i in range(nclusters)]
        for i in range(len(clusterid)):
            self.assertEqual(clusterid[i], mapping[correct[i]])

        cdata, cmask = clustercentroids(data, mask=mask, clusterid=clusterid,
                                        method="a", transpose=True)

        self.assertEqual(cdata.shape, (nrows, nclusters))
        self.assertEqual(cmask.shape, (nrows, nclusters))
        for value in cmask.flat:
            self.assertEqual(value, 1)

        correct = numpy.array([[1.1, 3.6666666667, 4.4],
                               [3.1, 2.0000000000, 2.4],
                               [4.1, 1.0000000000, 5.4],
                               [9.9, 0.6666666667, 5.0]])
        for i in range(nrows):
            for j in range(nclusters):
                self.assertAlmostEqual(cdata[i, mapping[j]], correct[i, j])

        # Second data set
        weight = numpy.array([1, 1])
        data = numpy.array([[1.1, 1.2],
                            [1.4, 1.3],
                            [1.1, 1.5],
                            [2.0, 1.5],
                            [1.7, 1.9],
                            [1.7, 1.9],
                            [5.7, 5.9],
                            [5.7, 5.9],
                            [3.1, 3.3],
                            [5.4, 5.3],
                            [5.1, 5.5],
                            [5.0, 5.5],
                            [5.1, 5.2]])
        mask = numpy.array([[1, 1],
                            [1, 1],
                            [1, 1],
                            [1, 1],
                            [1, 1],
                            [1, 1],
                            [1, 1],
                            [1, 1],
                            [1, 1],
                            [1, 1],
                            [1, 1],
                            [1, 1],
                            [1, 1]], int)
        nrows, ncols = data.shape

        clusterid, error, nfound = kcluster(data, nclusters=3, mask=mask,
                                            weight=weight, transpose=False,
                                            npass=100, method="a", dist="e")
        self.assertEqual(len(clusterid), len(data))

        correct = [0, 0, 0, 0, 0, 0, 1, 1, 2, 1, 1, 1, 1]
        mapping = [clusterid[correct.index(i)] for i in range(nclusters)]
        for i in range(len(clusterid)):
            self.assertEqual(clusterid[i], mapping[correct[i]])

        cdata, cmask = clustercentroids(data, mask=mask, clusterid=clusterid,
                                        method="a", transpose=False)

        self.assertEqual(cdata.shape, (nclusters, ncols))
        self.assertEqual(cmask.shape, (nclusters, ncols))
        for value in cmask.flat:
            self.assertEqual(value, 1)

        correct = numpy.array([[1.5000000, 1.55],
                               [5.3333333, 5.55],
                               [3.1000000, 3.30]])
        for i in range(nclusters):
            for j in range(ncols):
                self.assertAlmostEqual(cdata[mapping[i], j], correct[i, j])

    def test_clusterdistance_arguments(self):
        # Test if incorrect arguments are caught by the C code
        if TestCluster.module == "Bio.Cluster":
            from Bio.Cluster._cluster import clusterdistance
        elif TestCluster.module == "Pycluster":
            from Pycluster._cluster import clusterdistance

        # First data set
        weight = numpy.array([1.0, 1.0, 1.0, 1.0, 1.0])
        data = numpy.array([[1.1, 2.2, 3.3, 4.4, 5.5],
                            [3.1, 3.2, 1.3, 2.4, 1.5],
                            [4.1, 2.2, 0.3, 5.4, 0.5],
                            [9.9, 2.0, 0.0, 5.0, 0.0]])
        mask = numpy.array([[1, 1, 1, 1, 1],
                            [1, 1, 1, 1, 1],
                            [1, 1, 1, 1, 1],
                            [1, 1, 1, 1, 1]], numpy.int32)

        # Cluster assignments
        c1 = numpy.array([0], numpy.int32)
        c2 = numpy.array([1, 2], numpy.int32)
        c3 = numpy.array([3], numpy.int32)

        message = "^data is None$"
        with self.assertRaisesRegex(RuntimeError, message):
            clusterdistance(data=None, mask=mask, weight=weight,
                            index1=c1, index2=c2, dist="e",
                            method="a", transpose=False)
        message = "^data matrix has unexpected format.$"
        with self.assertRaisesRegex(RuntimeError, message):
            clusterdistance(data=[None], mask=mask, weight=weight,
                            index1=c1, index2=c2, dist="e",
                            method="a", transpose=False)
        message = "^data matrix has incorrect rank 1 \\(expected 2\\)$"
        with self.assertRaisesRegex(RuntimeError, message):
            clusterdistance(data=numpy.zeros(3), mask=mask, weight=weight,
                            index1=c1, index2=c2, dist="e",
                            method="a", transpose=False)
        message = "^data matrix has incorrect data type$"
        with self.assertRaisesRegex(RuntimeError, message):
            clusterdistance(data=numpy.zeros((3, 3), dtype=numpy.int16),
                            mask=mask, weight=weight,
                            index1=c1, index2=c2, dist="e",
                            method="a", transpose=False)
        message = "^data matrix is empty$"
        with self.assertRaisesRegex(ValueError, message):
            clusterdistance(data=data[:0], mask=mask, weight=weight,
                            index1=c1, index2=c2, dist="e",
                            method="a", transpose=False)
        message = "^data is not contiguous$"
        with self.assertRaisesRegex(RuntimeError, message):
            clusterdistance(data=data[:, ::2], mask=mask, weight=weight,
                            index1=c1, index2=c2, dist="e",
                            method="a", transpose=False)
        message = "^mask has unexpected format.$"
        with self.assertRaisesRegex(RuntimeError, message):
            clusterdistance(data=data, mask=[None], weight=weight,
                            index1=c1, index2=c2, dist="e",
                            method="a", transpose=False)
        message = "^mask has incorrect rank 1 \\(expected 2\\)$"
        with self.assertRaisesRegex(ValueError, message):
            clusterdistance(data=data, mask=numpy.zeros(3), weight=weight,
                            index1=c1, index2=c2, dist="e",
                            method="a", transpose=False)
        message = "^mask has incorrect data type$"
        with self.assertRaisesRegex(RuntimeError, message):
            clusterdistance(data=data,
                            mask=numpy.ones((2, 2), dtype=numpy.int16),
                            weight=weight, index1=c1, index2=c2, dist="e",
                            method="a", transpose=False)
        message = "^mask is not contiguous$"
        with self.assertRaisesRegex(RuntimeError, message):
            clusterdistance(data=data, mask=mask[:, ::2], weight=weight,
                            index1=c1, index2=c2, dist="e",
                            method="a", transpose=False)
        message = "^unexpected format.$"
        with self.assertRaisesRegex(RuntimeError, message):
            clusterdistance(data=data, mask=mask, weight="nothing",
                            index1=c1, index2=c2, dist="e",
                            method="a", transpose=False)
        message = "^incorrect rank 2 \\(expected 1\\)$"
        with self.assertRaisesRegex(ValueError, message):
            clusterdistance(data=data, mask=mask, weight=numpy.zeros((2, 2)),
                            index1=c1, index2=c2, dist="e",
                            method="a", transpose=False)
        message = "^array has incorrect data type$"
        with self.assertRaisesRegex(RuntimeError, message):
            clusterdistance(data=data, mask=mask,
                            weight=numpy.ones(3, dtype=numpy.int16),
                            index1=c1, index2=c2, dist="e",
                            method="a", transpose=False)
        message = "^unexpected format.$"
        with self.assertRaisesRegex(RuntimeError, message):
            clusterdistance(data=data, mask=mask, weight=weight,
                            index1=None, index2=c2, dist="e",
                            method="a", transpose=False)
        message = "^incorrect rank 2 \\(expected 1\\)$"
        with self.assertRaisesRegex(ValueError, message):
            clusterdistance(data=data, mask=mask, weight=weight,
                            index1=numpy.zeros((2, 2)), index2=c2, dist="e",
                            method="a", transpose=False)
        message = "^argument has incorrect data type$"
        with self.assertRaisesRegex(RuntimeError, message):
            clusterdistance(data=data, mask=mask, weight=weight,
                            index1=numpy.zeros(2, numpy.int16), index2=c2,
                            dist="e", method="a", transpose=False)
        message = "^unexpected format.$"
        with self.assertRaisesRegex(RuntimeError, message):
            clusterdistance(data=data, mask=mask, weight=weight,
                            index1=c1, index2=None, dist="e",
                            method="a", transpose=False)
        message = "^incorrect rank 2 \\(expected 1\\)$"
        with self.assertRaisesRegex(ValueError, message):
            clusterdistance(data=data, mask=mask, weight=weight,
                            index1=c1, index2=numpy.zeros((2, 2)), dist="e",
                            method="a", transpose=False)
        message = "^argument has incorrect data type$"
        with self.assertRaisesRegex(RuntimeError, message):
            clusterdistance(data=data, mask=mask, weight=weight,
                            index1=c1, index2=numpy.zeros(2, numpy.int16),
                            dist="e", method="a", transpose=False)

    def test_clusterdistance(self):
        if TestCluster.module == "Bio.Cluster":
            from Bio.Cluster import clusterdistance
        elif TestCluster.module == "Pycluster":
            from Pycluster import clusterdistance

        # First data set
        weight = numpy.array([1, 1, 1, 1, 1])
        data = numpy.array([[1.1, 2.2, 3.3, 4.4, 5.5],
                            [3.1, 3.2, 1.3, 2.4, 1.5],
                            [4.1, 2.2, 0.3, 5.4, 0.5],
                            [9.9, 2.0, 0.0, 5.0, 0.0]])
        mask = numpy.array([[1, 1, 1, 1, 1],
                            [1, 1, 1, 1, 1],
                            [1, 1, 1, 1, 1],
                            [1, 1, 1, 1, 1]], int)

        # Cluster assignments
        c1 = [0]
        c2 = [1, 2]
        c3 = [3]

        distance = clusterdistance(data, mask=mask, weight=weight,
                                   index1=c1, index2=c2, dist="e",
                                   method="a", transpose=False)
        self.assertAlmostEqual(distance, 6.650, places=3)
        distance = clusterdistance(data, mask=mask, weight=weight,
                                   index1=c1, index2=c3, dist="e",
                                   method="a", transpose=False)
        self.assertAlmostEqual(distance, 23.796, places=3)
        distance = clusterdistance(data, mask=mask, weight=weight,
                                   index1=c2, index2=c3, dist="e",
                                   method="a", transpose=False)
        self.assertAlmostEqual(distance, 8.606, places=3)

        # First data set, using transpose=True
        weight = numpy.array([1, 1, 1, 1])

        # Cluster assignments
        c1 = [0, 2]
        c2 = [1, 4]
        c3 = [3]

        distance = clusterdistance(data, mask=mask, weight=weight,
                                   index1=c1, index2=c2, dist="e",
                                   method="a", transpose=True)
        self.assertAlmostEqual(distance, 4.7675, places=3)
        distance = clusterdistance(data, mask=mask, weight=weight,
                                   index1=c1, index2=c3, dist="e",
                                   method="a", transpose=True)
        self.assertAlmostEqual(distance, 3.780625, places=3)
        distance = clusterdistance(data, mask=mask, weight=weight,
                                   index1=c2, index2=c3, dist="e",
                                   method="a", transpose=True)
        self.assertAlmostEqual(distance, 8.176875, places=3)

        # Second data set
        weight = numpy.array([1, 1])
        data = numpy.array([[1.1, 1.2],
                            [1.4, 1.3],
                            [1.1, 1.5],
                            [2.0, 1.5],
                            [1.7, 1.9],
                            [1.7, 1.9],
                            [5.7, 5.9],
                            [5.7, 5.9],
                            [3.1, 3.3],
                            [5.4, 5.3],
                            [5.1, 5.5],
                            [5.0, 5.5],
                            [5.1, 5.2]])
        mask = numpy.array([[1, 1],
                            [1, 1],
                            [1, 1],
                            [1, 1],
                            [1, 1],
                            [1, 1],
                            [1, 1],
                            [1, 1],
                            [1, 1],
                            [1, 1],
                            [1, 1],
                            [1, 1],
                            [1, 1]], int)

        # Cluster assignments
        c1 = [0, 1, 2, 3]
        c2 = [4, 5, 6, 7]
        c3 = [8]

        distance = clusterdistance(data, mask=mask, weight=weight,
                                   index1=c1, index2=c2, dist="e",
                                   method="a", transpose=False)
        self.assertAlmostEqual(distance, 5.833, places=3)
        distance = clusterdistance(data, mask=mask, weight=weight,
                                   index1=c1, index2=c3, dist="e",
                                   method="a", transpose=False)
        self.assertAlmostEqual(distance, 3.298, places=3)
        distance = clusterdistance(data, mask=mask, weight=weight,
                                   index1=c2, index2=c3, dist="e",
                                   method="a", transpose=False)
        self.assertAlmostEqual(distance, 0.360, places=3)

    def test_treecluster_arguments(self):
        # Test if incorrect arguments are caught by the C code
        if TestCluster.module == "Bio.Cluster":
            from Bio.Cluster._cluster import treecluster, Tree
        elif TestCluster.module == "Pycluster":
            from Pycluster._cluster import treecluster, Tree
        weight = numpy.array([1.0, 1.0, 1.0, 1.0, 1.0])
        data = numpy.array([[1.1, 2.2, 3.3, 4.4, 5.5],
                            [3.1, 3.2, 1.3, 2.4, 1.5],
                            [4.1, 2.2, 0.3, 5.4, 0.5],
                            [9.7, 2.0, 0.0, 5.0, 0.0]])
        mask = numpy.array([[1, 1, 1, 1, 1],
                            [1, 1, 1, 1, 1],
                            [1, 1, 1, 1, 1],
                            [1, 1, 1, 1, 1]], numpy.int32)

        message = "^argument 1 must be _cluster.Tree, not None$"
        with self.assertRaisesRegex(TypeError, message):
            treecluster(None, data=data, mask=mask, weight=weight,
                        transpose=False, method="a", dist="e",
                        distancematrix=None)
        tree = Tree()
        message = "^neither data nor distancematrix was given$"
        with self.assertRaisesRegex(ValueError, message):
            treecluster(tree, data=None, mask=mask, weight=weight,
                        transpose=False, method="a", dist="e",
                        distancematrix=None)
        message = "^data matrix has unexpected format.$"
        with self.assertRaisesRegex(RuntimeError, message):
            treecluster(tree, data=[], mask=mask, weight=weight,
                        transpose=False, method="a", dist="e",
                        distancematrix=None)
        message = "^data matrix has incorrect data type$"
        with self.assertRaisesRegex(RuntimeError, message):
            treecluster(tree, data=numpy.zeros((3, 3), numpy.int32), mask=mask,
                        weight=weight, transpose=False, method="a", dist="e",
                        distancematrix=None)
        message = "^data matrix has incorrect rank 1 \\(expected 2\\)$"
        with self.assertRaisesRegex(RuntimeError, message):
            treecluster(tree, data=numpy.zeros(3), mask=mask, weight=weight,
                        transpose=False, method="a", dist="e",
                        distancematrix=None)

    def test_tree_arguments(self):
        # Test if incorrect arguments are caught by the C code
        if TestCluster.module == "Bio.Cluster":
            from Bio.Cluster._cluster import Node, Tree
        elif TestCluster.module == "Pycluster":
            from Pycluster._cluster import Node, Tree

        nodes = [Node(1, 2, 0.2), Node(0, -1, 0.5), Node(3, -2, 0.6)]
        indices = numpy.zeros(4, numpy.int32)
        tree = Tree(nodes)
        message = "^requested number of clusters should be positive$"
        with self.assertRaisesRegex(ValueError, message):
            tree.cut(indices, -5)
        message = "^more clusters requested than items available$"
        with self.assertRaisesRegex(ValueError, message):
            tree.cut(indices, +5)
        message = "^unexpected format.$"
        with self.assertRaisesRegex(RuntimeError, message):
            tree.sort(indices, "nothing")
        message = "^incorrect rank 2 \\(expected 1\\)$"
        with self.assertRaisesRegex(ValueError, message):
            tree.sort(indices, numpy.zeros((5, 5)))
        message = "^order array has incorrect size 2 \\(expected 4\\)$"
        with self.assertRaisesRegex(ValueError, message):
            tree.sort(indices, numpy.zeros(2))
        message = "^order array has incorrect size 6 \\(expected 4\\)$"
        with self.assertRaisesRegex(ValueError, message):
            tree.sort(indices, numpy.zeros(6))

    def test_tree(self):
        if TestCluster.module == "Bio.Cluster":
            from Bio.Cluster import Node, Tree
        elif TestCluster.module == "Pycluster":
            from Pycluster import Node, Tree

        node = Node(2, 3)
        self.assertEqual(node.left, 2)
        self.assertEqual(node.right, 3)
        self.assertAlmostEqual(node.distance, 0.0, places=3)
        node.left = 6
        node.right = 2
        node.distance = 0.73
        self.assertEqual(node.left, 6)
        self.assertEqual(node.right, 2)
        self.assertAlmostEqual(node.distance, 0.73, places=3)
        nodes = [Node(1, 2, 0.2), Node(0, 3, 0.5), Node(-2, 4, 0.6), Node(-1, -3, 0.9)]
        try:
            tree = Tree(nodes)
        except Exception:
            self.fail("failed to construct tree from nodes")
        nodes = [Node(1, 2, 0.2), Node(0, 2, 0.5)]
        self.assertRaises(ValueError, Tree, nodes)
        nodes = [Node(1, 2, 0.2), Node(0, -1, 0.5)]
        tree = Tree(nodes)
        self.assertEqual(tree[0].left, 1)
        self.assertEqual(tree[0].right, 2)
        self.assertAlmostEqual(tree[0].distance, 0.2)
        self.assertEqual(tree[1].left, 0)
        self.assertEqual(tree[1].right, -1)
        self.assertAlmostEqual(tree[1].distance, 0.5)
        tree = Tree([Node(1, 2, 0.1), Node(0, -1, 0.5), Node(-2, 3, 0.9)])
        nodes = tree[:]
        nodes[0] = Node(0, 1, 0.2)
        nodes[1].left = 2
        tree = Tree(nodes)
        self.assertEqual(tree[0].left, 0)
        self.assertEqual(tree[0].right, 1)
        self.assertAlmostEqual(tree[0].distance, 0.2)
        self.assertEqual(tree[1].left, 2)
        self.assertEqual(tree[1].right, -1)
        self.assertAlmostEqual(tree[1].distance, 0.5)
        self.assertEqual(tree[2].left, -2)
        self.assertEqual(tree[2].right, 3)
        self.assertAlmostEqual(tree[2].distance, 0.9)

    def test_treecluster(self):
        if TestCluster.module == "Bio.Cluster":
            from Bio.Cluster import treecluster
        elif TestCluster.module == "Pycluster":
            from Pycluster import treecluster

        # First data set
        weight1 = [1, 1, 1, 1, 1]
        data1 = numpy.array([[1.1, 2.2, 3.3, 4.4, 5.5],
                             [3.1, 3.2, 1.3, 2.4, 1.5],
                             [4.1, 2.2, 0.3, 5.4, 0.5],
                             [9.7, 2.0, 0.0, 5.0, 0.0]])
        mask1 = numpy.array([[1, 1, 1, 1, 1],
                             [1, 1, 1, 1, 1],
                             [1, 1, 1, 1, 1],
                             [1, 1, 1, 1, 1]], int)

        # Pairwise average-linkage clustering
        tree = treecluster(data=data1, mask=mask1, weight=weight1,
                           transpose=False, method="a", dist="e")
        self.assertEqual(len(tree), len(data1) - 1)
        self.assertEqual(tree[0].left, 2)
        self.assertEqual(tree[0].right, 1)
        self.assertAlmostEqual(tree[0].distance, 2.600, places=3)
        self.assertEqual(tree[1].left, -1)
        self.assertEqual(tree[1].right, 0)
        self.assertAlmostEqual(tree[1].distance, 7.300, places=3)
        self.assertEqual(tree[2].left, 3)
        self.assertEqual(tree[2].right, -2)
        self.assertAlmostEqual(tree[2].distance, 13.540, places=3)
        indices = tree.cut(nclusters=1)
        self.assertEqual(len(indices), len(data1))
        self.assertEqual(indices[0], 0)
        self.assertEqual(indices[1], 0)
        self.assertEqual(indices[2], 0)
        self.assertEqual(indices[3], 0)
        indices = tree.cut(nclusters=2)
        self.assertEqual(len(indices), len(data1))
        self.assertEqual(indices[0], 1)
        self.assertEqual(indices[1], 1)
        self.assertEqual(indices[2], 1)
        self.assertEqual(indices[3], 0)
        indices = tree.cut(nclusters=3)
        self.assertEqual(len(indices), len(data1))
        self.assertEqual(indices[0], 2)
        self.assertEqual(indices[1], 1)
        self.assertEqual(indices[2], 1)
        self.assertEqual(indices[3], 0)
        indices = tree.cut(nclusters=4)
        self.assertEqual(len(indices), len(data1))
        self.assertEqual(indices[0], 3)
        self.assertEqual(indices[1], 2)
        self.assertEqual(indices[2], 1)
        self.assertEqual(indices[3], 0)
        indices = tree.sort([0, 1, 2, 3])
        self.assertEqual(len(indices), len(data1))
        self.assertEqual(indices[0], 0)
        self.assertEqual(indices[1], 1)
        self.assertEqual(indices[2], 2)
        self.assertEqual(indices[3], 3)
        indices = tree.sort([0, 3, 2, 1])
        self.assertEqual(len(indices), len(data1))
        self.assertEqual(indices[0], 3)
        self.assertEqual(indices[1], 0)
        self.assertEqual(indices[2], 2)
        self.assertEqual(indices[3], 1)

        # Pairwise single-linkage clustering
        tree = treecluster(data=data1, mask=mask1, weight=weight1,
                           transpose=False, method="s", dist="e")
        self.assertEqual(len(tree), len(data1) - 1)
        self.assertEqual(tree[0].left, 1)
        self.assertEqual(tree[0].right, 2)
        self.assertAlmostEqual(tree[0].distance, 2.600, places=3)
        self.assertEqual(tree[1].left, 0)
        self.assertEqual(tree[1].right, -1)
        self.assertAlmostEqual(tree[1].distance, 5.800, places=3)
        self.assertEqual(tree[2].left, -2)
        self.assertEqual(tree[2].right, 3)
        self.assertAlmostEqual(tree[2].distance, 6.380, places=3)
        indices = tree.cut(nclusters=1)
        self.assertEqual(len(indices), len(data1))
        self.assertEqual(indices[0], 0)
        self.assertEqual(indices[1], 0)
        self.assertEqual(indices[2], 0)
        self.assertEqual(indices[3], 0)
        indices = tree.cut(nclusters=2)
        self.assertEqual(len(indices), len(data1))
        self.assertEqual(indices[0], 0)
        self.assertEqual(indices[1], 0)
        self.assertEqual(indices[2], 0)
        self.assertEqual(indices[3], 1)
        indices = tree.cut(nclusters=3)
        self.assertEqual(len(indices), len(data1))
        self.assertEqual(indices[0], 0)
        self.assertEqual(indices[1], 1)
        self.assertEqual(indices[2], 1)
        self.assertEqual(indices[3], 2)
        indices = tree.cut(nclusters=4)
        self.assertEqual(len(indices), len(data1))
        self.assertEqual(indices[0], 0)
        self.assertEqual(indices[1], 1)
        self.assertEqual(indices[2], 2)
        self.assertEqual(indices[3], 3)
        indices = tree.sort([0, 1, 2, 3])
        self.assertEqual(len(indices), len(data1))
        self.assertEqual(indices[0], 0)
        self.assertEqual(indices[1], 1)
        self.assertEqual(indices[2], 2)
        self.assertEqual(indices[3], 3)
        indices = tree.sort([0, 3, 2, 1])
        self.assertEqual(len(indices), len(data1))
        self.assertEqual(indices[0], 3)
        self.assertEqual(indices[1], 0)
        self.assertEqual(indices[2], 2)
        self.assertEqual(indices[3], 1)

        # Pairwise centroid-linkage clustering
        tree = treecluster(data=data1, mask=mask1, weight=weight1,
                           transpose=False, method="c", dist="e")
        self.assertEqual(len(tree), len(data1) - 1)
        self.assertEqual(tree[0].left, 1)
        self.assertEqual(tree[0].right, 2)
        self.assertAlmostEqual(tree[0].distance, 2.600, places=3)
        self.assertEqual(tree[1].left, 0)
        self.assertEqual(tree[1].right, -1)
        self.assertAlmostEqual(tree[1].distance, 6.650, places=3)
        self.assertEqual(tree[2].left, -2)
        self.assertEqual(tree[2].right, 3)
        self.assertAlmostEqual(tree[2].distance, 11.629, places=3)
        indices = tree.sort([0, 1, 2, 3])
        self.assertEqual(len(indices), len(data1))
        self.assertEqual(indices[0], 0)
        self.assertEqual(indices[1], 1)
        self.assertEqual(indices[2], 2)
        self.assertEqual(indices[3], 3)
        indices = tree.cut(nclusters=1)
        self.assertEqual(len(indices), len(data1))
        self.assertEqual(indices[0], 0)
        self.assertEqual(indices[1], 0)
        self.assertEqual(indices[2], 0)
        self.assertEqual(indices[3], 0)
        indices = tree.cut(nclusters=2)
        self.assertEqual(len(indices), len(data1))
        self.assertEqual(indices[0], 0)
        self.assertEqual(indices[1], 0)
        self.assertEqual(indices[2], 0)
        self.assertEqual(indices[3], 1)
        indices = tree.cut(nclusters=3)
        self.assertEqual(len(indices), len(data1))
        self.assertEqual(indices[0], 0)
        self.assertEqual(indices[1], 1)
        self.assertEqual(indices[2], 1)
        self.assertEqual(indices[3], 2)
        indices = tree.cut(nclusters=4)
        self.assertEqual(len(indices), len(data1))
        self.assertEqual(indices[0], 0)
        self.assertEqual(indices[1], 1)
        self.assertEqual(indices[2], 2)
        self.assertEqual(indices[3], 3)
        indices = tree.sort([0, 3, 2, 1])
        self.assertEqual(len(indices), len(data1))
        self.assertEqual(indices[0], 3)
        self.assertEqual(indices[1], 0)
        self.assertEqual(indices[2], 2)
        self.assertEqual(indices[3], 1)

        # Pairwise maximum-linkage clustering
        tree = treecluster(data=data1, mask=mask1, weight=weight1,
                           transpose=False, method="m", dist="e")
        self.assertEqual(len(tree), len(data1) - 1)
        self.assertEqual(tree[0].left, 2)
        self.assertEqual(tree[0].right, 1)
        self.assertAlmostEqual(tree[0].distance, 2.600, places=3)
        self.assertEqual(tree[1].left, -1)
        self.assertEqual(tree[1].right, 0)
        self.assertAlmostEqual(tree[1].distance, 8.800, places=3)
        self.assertEqual(tree[2].left, 3)
        self.assertEqual(tree[2].right, -2)
        self.assertAlmostEqual(tree[2].distance, 23.100, places=3)
        indices = tree.cut(nclusters=1)
        self.assertEqual(len(indices), len(data1))
        self.assertEqual(indices[0], 0)
        self.assertEqual(indices[1], 0)
        self.assertEqual(indices[2], 0)
        self.assertEqual(indices[3], 0)
        indices = tree.cut(nclusters=2)
        self.assertEqual(len(indices), len(data1))
        self.assertEqual(indices[0], 1)
        self.assertEqual(indices[1], 1)
        self.assertEqual(indices[2], 1)
        self.assertEqual(indices[3], 0)
        indices = tree.cut(nclusters=3)
        self.assertEqual(len(indices), len(data1))
        self.assertEqual(indices[0], 2)
        self.assertEqual(indices[1], 1)
        self.assertEqual(indices[2], 1)
        self.assertEqual(indices[3], 0)
        indices = tree.cut(nclusters=4)
        self.assertEqual(len(indices), len(data1))
        self.assertEqual(indices[0], 3)
        self.assertEqual(indices[1], 2)
        self.assertEqual(indices[2], 1)
        self.assertEqual(indices[3], 0)
        indices = tree.sort([0, 1, 2, 3])
        self.assertEqual(len(indices), len(data1))
        self.assertEqual(indices[0], 0)
        self.assertEqual(indices[1], 1)
        self.assertEqual(indices[2], 2)
        self.assertEqual(indices[3], 3)
        indices = tree.sort([0, 3, 2, 1])
        self.assertEqual(len(indices), len(data1))
        self.assertEqual(indices[0], 3)
        self.assertEqual(indices[1], 0)
        self.assertEqual(indices[2], 2)
        self.assertEqual(indices[3], 1)

        # First data set, using transpose=True
        weight1 = [1, 1, 1, 1]
        data1 = numpy.array([[1.1, 2.2, 3.3, 4.4, 5.5],
                             [3.1, 3.2, 1.3, 2.4, 1.5],
                             [4.1, 2.2, 0.3, 5.4, 0.5],
                             [9.7, 2.0, 0.0, 5.0, 0.0]])
        mask1 = numpy.array([[1, 1, 1, 1, 1],
                             [1, 1, 1, 1, 1],
                             [1, 1, 1, 1, 1],
                             [1, 1, 1, 1, 1]], int)
        nrows, ncols = data1.shape

        # Pairwise average-linkage clustering
        tree = treecluster(data=data1, mask=mask1, weight=weight1,
                           transpose=True, method="a", dist="e")
        self.assertEqual(len(tree), ncols - 1)
        self.assertEqual(tree[0].left, 4)
        self.assertEqual(tree[0].right, 2)
        self.assertAlmostEqual(tree[0].distance, 1.230, places=3)
        self.assertEqual(tree[1].left, -1)
        self.assertEqual(tree[1].right, 1)
        self.assertAlmostEqual(tree[1].distance, 4.1375, places=3)
        self.assertEqual(tree[2].left, 3)
        self.assertEqual(tree[2].right, 0)
        self.assertAlmostEqual(tree[2].distance, 8.790, places=3)
        self.assertEqual(tree[3].left, -2)
        self.assertEqual(tree[3].right, -3)
        self.assertAlmostEqual(tree[3].distance, 18.2867, places=3)
        indices = tree.cut(nclusters=1)
        self.assertEqual(len(indices), ncols)
        self.assertEqual(indices[0], 0)
        self.assertEqual(indices[1], 0)
        self.assertEqual(indices[2], 0)
        self.assertEqual(indices[3], 0)
        self.assertEqual(indices[4], 0)
        indices = tree.cut(nclusters=2)
        self.assertEqual(len(indices), ncols)
        self.assertEqual(indices[0], 1)
        self.assertEqual(indices[1], 0)
        self.assertEqual(indices[2], 0)
        self.assertEqual(indices[3], 1)
        self.assertEqual(indices[4], 0)
        indices = tree.cut(nclusters=3)
        self.assertEqual(len(indices), ncols)
        self.assertEqual(indices[0], 2)
        self.assertEqual(indices[1], 0)
        self.assertEqual(indices[2], 0)
        self.assertEqual(indices[3], 1)
        self.assertEqual(indices[4], 0)
        indices = tree.cut(nclusters=4)
        self.assertEqual(len(indices), ncols)
        self.assertEqual(indices[0], 3)
        self.assertEqual(indices[1], 1)
        self.assertEqual(indices[2], 0)
        self.assertEqual(indices[3], 2)
        self.assertEqual(indices[4], 0)
        indices = tree.sort([0, 1, 2, 3, 4])
        self.assertEqual(len(indices), ncols)
        self.assertEqual(indices[0], 0)
        self.assertEqual(indices[1], 3)
        self.assertEqual(indices[2], 1)
        self.assertEqual(indices[3], 2)
        self.assertEqual(indices[4], 4)
        indices = tree.sort([0, 4, 3, 2, 1])
        self.assertEqual(len(indices), ncols)
        self.assertEqual(indices[0], 0)
        self.assertEqual(indices[1], 3)
        self.assertEqual(indices[2], 4)
        self.assertEqual(indices[3], 2)
        self.assertEqual(indices[4], 1)

        # Pairwise single-linkage clustering
        tree = treecluster(data=data1, mask=mask1, weight=weight1,
                           transpose=True, method="s", dist="e")
        self.assertEqual(len(tree), ncols - 1)
        self.assertEqual(tree[0].left, 2)
        self.assertEqual(tree[0].right, 4)
        self.assertAlmostEqual(tree[0].distance, 1.230, places=3)
        self.assertEqual(tree[1].left, 1)
        self.assertEqual(tree[1].right, -1)
        self.assertAlmostEqual(tree[1].distance, 3.1075, places=3)
        self.assertEqual(tree[2].left, 3)
        self.assertEqual(tree[2].right, -2)
        self.assertAlmostEqual(tree[2].distance, 6.180, places=3)
        self.assertEqual(tree[3].left, 0)
        self.assertEqual(tree[3].right, -3)
        self.assertAlmostEqual(tree[3].distance, 8.790, places=3)
        indices = tree.cut(nclusters=1)
        self.assertEqual(len(indices), ncols)
        self.assertEqual(indices[0], 0)
        self.assertEqual(indices[1], 0)
        self.assertEqual(indices[2], 0)
        self.assertEqual(indices[3], 0)
        self.assertEqual(indices[4], 0)
        indices = tree.cut(nclusters=2)
        self.assertEqual(len(indices), ncols)
        self.assertEqual(indices[0], 0)
        self.assertEqual(indices[1], 1)
        self.assertEqual(indices[2], 1)
        self.assertEqual(indices[3], 1)
        self.assertEqual(indices[4], 1)
        indices = tree.cut(nclusters=3)
        self.assertEqual(len(indices), ncols)
        self.assertEqual(indices[0], 0)
        self.assertEqual(indices[1], 2)
        self.assertEqual(indices[2], 2)
        self.assertEqual(indices[3], 1)
        self.assertEqual(indices[4], 2)
        indices = tree.cut(nclusters=4)
        self.assertEqual(len(indices), ncols)
        self.assertEqual(indices[0], 0)
        self.assertEqual(indices[1], 2)
        self.assertEqual(indices[2], 3)
        self.assertEqual(indices[3], 1)
        self.assertEqual(indices[4], 3)
        indices = tree.sort([0, 1, 2, 3, 4])
        self.assertEqual(len(indices), ncols)
        self.assertEqual(indices[0], 0)
        self.assertEqual(indices[1], 1)
        self.assertEqual(indices[2], 2)
        self.assertEqual(indices[3], 4)
        self.assertEqual(indices[4], 3)
        indices = tree.sort([0, 4, 3, 2, 1])
        self.assertEqual(len(indices), ncols)
        self.assertEqual(indices[0], 0)
        self.assertEqual(indices[1], 3)
        self.assertEqual(indices[2], 4)
        self.assertEqual(indices[3], 2)
        self.assertEqual(indices[4], 1)

        # Pairwise centroid-linkage clustering
        tree = treecluster(data=data1, mask=mask1, weight=weight1,
                           transpose=True, method="c", dist="e")
        self.assertEqual(len(tree), ncols - 1)
        self.assertEqual(tree[0].left, 2)
        self.assertEqual(tree[0].right, 4)
        self.assertAlmostEqual(tree[0].distance, 1.23, places=3)
        self.assertEqual(tree[1].left, 1)
        self.assertEqual(tree[1].right, -1)
        self.assertAlmostEqual(tree[1].distance, 3.83, places=3)
        self.assertEqual(tree[2].left, 0)
        self.assertEqual(tree[2].right, 3)
        self.assertAlmostEqual(tree[2].distance, 8.79, places=3)
        self.assertEqual(tree[3].left, -3)
        self.assertEqual(tree[3].right, -2)
        self.assertAlmostEqual(tree[3].distance, 15.0331, places=3)
        indices = tree.sort([0, 1, 2, 3, 4])
        self.assertEqual(len(indices), ncols)
        self.assertEqual(indices[0], 0)
        self.assertEqual(indices[1], 3)
        self.assertEqual(indices[2], 1)
        self.assertEqual(indices[3], 2)
        self.assertEqual(indices[4], 4)
        indices = tree.cut(nclusters=1)
        self.assertEqual(len(indices), ncols)
        self.assertEqual(indices[0], 0)
        self.assertEqual(indices[1], 0)
        self.assertEqual(indices[2], 0)
        self.assertEqual(indices[3], 0)
        self.assertEqual(indices[4], 0)
        indices = tree.cut(nclusters=2)
        self.assertEqual(len(indices), ncols)
        self.assertEqual(indices[0], 0)
        self.assertEqual(indices[1], 1)
        self.assertEqual(indices[2], 1)
        self.assertEqual(indices[3], 0)
        self.assertEqual(indices[4], 1)
        indices = tree.cut(nclusters=3)
        self.assertEqual(len(indices), ncols)
        self.assertEqual(indices[0], 0)
        self.assertEqual(indices[1], 2)
        self.assertEqual(indices[2], 2)
        self.assertEqual(indices[3], 1)
        self.assertEqual(indices[4], 2)
        indices = tree.cut(nclusters=4)
        self.assertEqual(len(indices), ncols)
        self.assertEqual(indices[0], 0)
        self.assertEqual(indices[1], 2)
        self.assertEqual(indices[2], 3)
        self.assertEqual(indices[3], 1)
        self.assertEqual(indices[4], 3)
        indices = tree.sort([0, 4, 3, 2, 1])
        self.assertEqual(len(indices), ncols)
        self.assertEqual(indices[0], 0)
        self.assertEqual(indices[1], 3)
        self.assertEqual(indices[2], 4)
        self.assertEqual(indices[3], 2)
        self.assertEqual(indices[4], 1)

        # Pairwise maximum-linkage clustering
        tree = treecluster(data=data1, mask=mask1, weight=weight1,
                           transpose=True, method="m", dist="e")
        self.assertEqual(len(tree), ncols - 1)
        self.assertEqual(tree[0].left, 4)
        self.assertEqual(tree[0].right, 2)
        self.assertAlmostEqual(tree[0].distance, 1.230, places=3)
        self.assertEqual(tree[1].left, -1)
        self.assertEqual(tree[1].right, 1)
        self.assertAlmostEqual(tree[1].distance, 5.1675, places=3)
        self.assertEqual(tree[2].left, 3)
        self.assertEqual(tree[2].right, 0)
        self.assertAlmostEqual(tree[2].distance, 8.790, places=3)
        self.assertEqual(tree[3].left, -2)
        self.assertEqual(tree[3].right, -3)
        self.assertAlmostEqual(tree[3].distance, 32.2425, places=3)
        indices = tree.cut(nclusters=1)
        self.assertEqual(len(indices), ncols)
        self.assertEqual(indices[0], 0)
        self.assertEqual(indices[1], 0)
        self.assertEqual(indices[2], 0)
        self.assertEqual(indices[3], 0)
        self.assertEqual(indices[4], 0)
        indices = tree.cut(nclusters=2)
        self.assertEqual(len(indices), ncols)
        self.assertEqual(indices[0], 1)
        self.assertEqual(indices[1], 0)
        self.assertEqual(indices[2], 0)
        self.assertEqual(indices[3], 1)
        self.assertEqual(indices[4], 0)
        indices = tree.cut(nclusters=3)
        self.assertEqual(len(indices), ncols)
        self.assertEqual(indices[0], 2)
        self.assertEqual(indices[1], 0)
        self.assertEqual(indices[2], 0)
        self.assertEqual(indices[3], 1)
        self.assertEqual(indices[4], 0)
        indices = tree.cut(nclusters=4)
        self.assertEqual(len(indices), ncols)
        self.assertEqual(indices[0], 3)
        self.assertEqual(indices[1], 1)
        self.assertEqual(indices[2], 0)
        self.assertEqual(indices[3], 2)
        self.assertEqual(indices[4], 0)
        indices = tree.sort([0, 1, 2, 3, 4])
        self.assertEqual(len(indices), ncols)
        self.assertEqual(indices[0], 0)
        self.assertEqual(indices[1], 3)
        self.assertEqual(indices[2], 1)
        self.assertEqual(indices[3], 2)
        self.assertEqual(indices[4], 4)
        indices = tree.sort([0, 4, 3, 2, 1])
        self.assertEqual(len(indices), ncols)
        self.assertEqual(indices[0], 0)
        self.assertEqual(indices[1], 3)
        self.assertEqual(indices[2], 4)
        self.assertEqual(indices[3], 2)
        self.assertEqual(indices[4], 1)

        # Second data set
        weight2 = [1, 1]
        data2 = numpy.array([[0.8223, 0.9295],
                             [1.4365, 1.3223],
                             [1.1623, 1.5364],
                             [2.1826, 1.1934],
                             [1.7763, 1.9352],
                             [1.7215, 1.9912],
                             [2.1812, 5.9935],
                             [5.3290, 5.9452],
                             [3.1491, 3.3454],
                             [5.1923, 5.3156],
                             [4.7735, 5.4012],
                             [5.1297, 5.5645],
                             [5.3934, 5.1823]])
        mask2 = numpy.array([[1, 1],
                             [1, 1],
                             [1, 1],
                             [1, 1],
                             [1, 1],
                             [1, 1],
                             [1, 1],
                             [1, 1],
                             [1, 1],
                             [1, 1],
                             [1, 1],
                             [1, 1],
                             [1, 1]], int)

        # Test second data set
        # Pairwise average-linkage clustering
        tree = treecluster(data=data2, mask=mask2, weight=weight2,
                           transpose=False, method="a", dist="e")
        self.assertEqual(len(tree), len(data2) - 1)
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
        indices = tree.cut(nclusters=1)
        self.assertEqual(len(indices), len(data2))
        self.assertEqual(indices[0], 0)
        self.assertEqual(indices[1], 0)
        self.assertEqual(indices[2], 0)
        self.assertEqual(indices[3], 0)
        self.assertEqual(indices[4], 0)
        self.assertEqual(indices[5], 0)
        self.assertEqual(indices[6], 0)
        self.assertEqual(indices[7], 0)
        self.assertEqual(indices[8], 0)
        self.assertEqual(indices[9], 0)
        self.assertEqual(indices[10], 0)
        self.assertEqual(indices[11], 0)
        self.assertEqual(indices[12], 0)
        indices = tree.cut(nclusters=2)
        self.assertEqual(len(indices), len(data2))
        self.assertEqual(indices[0], 1)
        self.assertEqual(indices[1], 1)
        self.assertEqual(indices[2], 1)
        self.assertEqual(indices[3], 1)
        self.assertEqual(indices[4], 1)
        self.assertEqual(indices[5], 1)
        self.assertEqual(indices[6], 0)
        self.assertEqual(indices[7], 0)
        self.assertEqual(indices[8], 1)
        self.assertEqual(indices[9], 0)
        self.assertEqual(indices[10], 0)
        self.assertEqual(indices[11], 0)
        self.assertEqual(indices[12], 0)
        indices = tree.cut(nclusters=3)
        self.assertEqual(len(indices), len(data2))
        self.assertEqual(indices[0], 2)
        self.assertEqual(indices[1], 2)
        self.assertEqual(indices[2], 2)
        self.assertEqual(indices[3], 2)
        self.assertEqual(indices[4], 2)
        self.assertEqual(indices[5], 2)
        self.assertEqual(indices[6], 1)
        self.assertEqual(indices[7], 0)
        self.assertEqual(indices[8], 2)
        self.assertEqual(indices[9], 0)
        self.assertEqual(indices[10], 0)
        self.assertEqual(indices[11], 0)
        self.assertEqual(indices[12], 0)
        indices = tree.cut(nclusters=4)
        self.assertEqual(len(indices), len(data2))
        self.assertEqual(indices[0], 3)
        self.assertEqual(indices[1], 3)
        self.assertEqual(indices[2], 3)
        self.assertEqual(indices[3], 3)
        self.assertEqual(indices[4], 3)
        self.assertEqual(indices[5], 3)
        self.assertEqual(indices[6], 1)
        self.assertEqual(indices[7], 0)
        self.assertEqual(indices[8], 2)
        self.assertEqual(indices[9], 0)
        self.assertEqual(indices[10], 0)
        self.assertEqual(indices[11], 0)
        self.assertEqual(indices[12], 0)
        indices = tree.cut(nclusters=5)
        self.assertEqual(len(indices), len(data2))
        self.assertEqual(indices[0], 4)
        self.assertEqual(indices[1], 4)
        self.assertEqual(indices[2], 4)
        self.assertEqual(indices[3], 3)
        self.assertEqual(indices[4], 3)
        self.assertEqual(indices[5], 3)
        self.assertEqual(indices[6], 1)
        self.assertEqual(indices[7], 0)
        self.assertEqual(indices[8], 2)
        self.assertEqual(indices[9], 0)
        self.assertEqual(indices[10], 0)
        self.assertEqual(indices[11], 0)
        self.assertEqual(indices[12], 0)
        indices = tree.sort()
        self.assertEqual(len(indices), len(data2))
        self.assertEqual(indices[0], 7)
        self.assertEqual(indices[1], 11)
        self.assertEqual(indices[2], 9)
        self.assertEqual(indices[3], 12)
        self.assertEqual(indices[4], 10)
        self.assertEqual(indices[5], 6)
        self.assertEqual(indices[6], 8)
        self.assertEqual(indices[7], 5)
        self.assertEqual(indices[8], 4)
        self.assertEqual(indices[9], 3)
        self.assertEqual(indices[10], 2)
        self.assertEqual(indices[11], 1)
        self.assertEqual(indices[12], 0)

        # Pairwise single-linkage clustering
        tree = treecluster(data=data2, mask=mask2, weight=weight2,
                           transpose=False, method="s", dist="e")
        self.assertEqual(len(tree), len(data2) - 1)
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
        indices = tree.sort()
        self.assertEqual(len(indices), len(data2))
        self.assertEqual(indices[0], 6)
        self.assertEqual(indices[1], 3)
        self.assertEqual(indices[2], 0)
        self.assertEqual(indices[3], 1)
        self.assertEqual(indices[4], 2)
        self.assertEqual(indices[5], 4)
        self.assertEqual(indices[6], 5)
        self.assertEqual(indices[7], 8)
        self.assertEqual(indices[8], 7)
        self.assertEqual(indices[9], 10)
        self.assertEqual(indices[10], 11)
        self.assertEqual(indices[11], 9)
        self.assertEqual(indices[12], 12)

        # Pairwise centroid-linkage clustering
        tree = treecluster(data=data2, mask=mask2, weight=weight2,
                           transpose=False, method="c", dist="e")
        self.assertEqual(len(tree), len(data2) - 1)
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
        indices = tree.sort()
        self.assertEqual(len(indices), len(data2))
        self.assertEqual(indices[0], 0)
        self.assertEqual(indices[1], 1)
        self.assertEqual(indices[2], 2)
        self.assertEqual(indices[3], 3)
        self.assertEqual(indices[4], 4)
        self.assertEqual(indices[5], 5)
        self.assertEqual(indices[6], 8)
        self.assertEqual(indices[7], 6)
        self.assertEqual(indices[8], 10)
        self.assertEqual(indices[9], 12)
        self.assertEqual(indices[10], 9)
        self.assertEqual(indices[11], 11)
        self.assertEqual(indices[12], 7)

        # Pairwise maximum-linkage clustering
        tree = treecluster(data=data2, mask=mask2, weight=weight2,
                           transpose=False, method="m", dist="e")
        self.assertEqual(len(tree), len(data2) - 1)
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
        indices = tree.sort()
        self.assertEqual(len(indices), len(data2))
        self.assertEqual(indices[0], 8)
        self.assertEqual(indices[1], 6)
        self.assertEqual(indices[2], 9)
        self.assertEqual(indices[3], 12)
        self.assertEqual(indices[4], 11)
        self.assertEqual(indices[5], 10)
        self.assertEqual(indices[6], 7)
        self.assertEqual(indices[7], 5)
        self.assertEqual(indices[8], 4)
        self.assertEqual(indices[9], 3)
        self.assertEqual(indices[10], 2)
        self.assertEqual(indices[11], 1)
        self.assertEqual(indices[12], 0)

    def test_somcluster_arguments(self):
        # Test if incorrect arguments are caught by the C code
        if TestCluster.module == "Bio.Cluster":
            from Bio.Cluster._cluster import somcluster
        elif TestCluster.module == "Pycluster":
            from Pycluster._cluster import somcluster

        weight = numpy.array([1.0, 1.0, 1.0, 1.0, 1.0])
        data = numpy.array([[1.1, 2.2, 3.3, 4.4, 5.5],
                            [3.1, 3.2, 1.3, 2.4, 1.5],
                            [4.1, 2.2, 0.3, 5.4, 0.5],
                            [9.9, 2.0, 0.0, 5.0, 0.0]])
        mask = numpy.array([[1, 1, 1, 1, 1],
                            [1, 1, 1, 1, 1],
                            [1, 1, 1, 1, 1],
                            [1, 1, 1, 1, 1]], numpy.int32)
        nitems, ndata = data.shape
        nxgrid, nygrid = 10, 10
        clusterids = numpy.ones((nitems, 2), numpy.int32)
        celldata = numpy.zeros((nxgrid, nygrid, ndata), dtype="d")

        message = "^unexpected format.$"
        with self.assertRaisesRegex(RuntimeError, message):
            somcluster(clusterids=None, celldata=celldata,
                       data=data, mask=mask, weight=weight,
                       transpose=False, inittau=0.02, niter=100, dist="e")
        message = "^incorrect rank 1 \\(expected 2\\)$"
        with self.assertRaisesRegex(ValueError, message):
            somcluster(clusterids=numpy.ones(nitems, numpy.int32),
                       celldata=celldata, data=data, mask=mask, weight=weight,
                       transpose=False, inittau=0.02, niter=100, dist="e")
        message = "^argument has incorrect data type$"
        with self.assertRaisesRegex(RuntimeError, message):
            somcluster(clusterids=numpy.ones((nitems, 2), numpy.int16),
                       celldata=celldata, data=data, mask=mask, weight=weight,
                       transpose=False, inittau=0.02, niter=100, dist="e")
        message = "^array has 3 columns \\(expected 2\\)$"
        with self.assertRaisesRegex(ValueError, message):
            somcluster(clusterids=numpy.ones((nitems, 3), numpy.int32),
                       celldata=celldata, data=data, mask=mask, weight=weight,
                       transpose=False, inittau=0.02, niter=100, dist="e")
        message = "^celldata array has unexpected format.$"
        with self.assertRaisesRegex(RuntimeError, message):
            somcluster(clusterids=clusterids, celldata=None,
                       data=data, mask=mask, weight=weight,
                       transpose=False, inittau=0.02, niter=100, dist="e")
        message = "^celldata array has incorrect data type$"
        with self.assertRaisesRegex(RuntimeError, message):
            somcluster(clusterids=clusterids,
                       celldata=numpy.zeros((nxgrid, nygrid, ndata),
                                            dtype=numpy.int32),
                       data=data, mask=mask, weight=weight,
                       transpose=False, inittau=0.02, niter=100, dist="e")
        message = "^data is None$"
        with self.assertRaisesRegex(RuntimeError, message):
            somcluster(clusterids=clusterids, celldata=celldata,
                       data=None, mask=mask, weight=weight,
                       transpose=False, inittau=0.02, niter=100, dist="e")
        message = "^data matrix has unexpected format.$"
        with self.assertRaisesRegex(RuntimeError, message):
            somcluster(clusterids=clusterids, celldata=celldata,
                       data=[None], mask=mask, weight=weight,
                       transpose=False, inittau=0.02, niter=100, dist="e")
        message = "^data matrix has incorrect data type$"
        with self.assertRaisesRegex(RuntimeError, message):
            somcluster(clusterids=clusterids, celldata=celldata,
                       data=numpy.zeros((4, 5), dtype=numpy.int16),
                       mask=mask, weight=weight,
                       transpose=False, inittau=0.02, niter=100, dist="e")
        message = "^data matrix has incorrect rank 1 \\(expected 2\\)$"
        with self.assertRaisesRegex(RuntimeError, message):
            somcluster(clusterids=clusterids, celldata=celldata,
                       data=numpy.zeros(4), mask=mask, weight=weight,
                       transpose=False, inittau=0.02, niter=100, dist="e")
        message = "^data matrix is empty$"
        with self.assertRaisesRegex(ValueError, message):
            somcluster(clusterids=clusterids, celldata=celldata,
                       data=data[:0], mask=mask, weight=weight,
                       transpose=False, inittau=0.02, niter=100, dist="e")
        message = "^data is not contiguous$"
        with self.assertRaisesRegex(RuntimeError, message):
            somcluster(clusterids=clusterids, celldata=celldata,
                       data=data[:, ::2], mask=mask, weight=weight,
                       transpose=False, inittau=0.02, niter=100, dist="e")
        message = "^mask is None$"
        with self.assertRaisesRegex(RuntimeError, message):
            somcluster(clusterids=clusterids, celldata=celldata,
                       data=data, mask=None, weight=weight,
                       transpose=False, inittau=0.02, niter=100, dist="e")
        message = "^mask has unexpected format.$"
        with self.assertRaisesRegex(RuntimeError, message):
            somcluster(clusterids=clusterids, celldata=celldata,
                       data=data, mask=[None], weight=weight,
                       transpose=False, inittau=0.02, niter=100, dist="e")
        message = "^mask has incorrect rank 1 \\(expected 2\\)$"
        with self.assertRaisesRegex(ValueError, message):
            somcluster(clusterids=clusterids, celldata=celldata,
                       data=data, mask=numpy.array([1, 1, 1]), weight=weight,
                       transpose=False, inittau=0.02, niter=100, dist="e")
        message = "^mask has incorrect data type$"
        with self.assertRaisesRegex(RuntimeError, message):
            somcluster(clusterids=clusterids, celldata=celldata, data=data,
                       mask=numpy.array([[1, 1], [1, 1]], dtype=numpy.int16),
                       weight=weight, transpose=False, inittau=0.02, niter=100,
                       dist="e")
        message = "^mask is not contiguous$"
        with self.assertRaisesRegex(RuntimeError, message):
            somcluster(clusterids=clusterids, celldata=celldata,
                       data=data, mask=mask[:, ::2], weight=weight,
                       transpose=False, inittau=0.02, niter=100, dist="e")
        message = "^unexpected format.$"
        with self.assertRaisesRegex(RuntimeError, message):
            somcluster(clusterids=clusterids, celldata=celldata,
                       data=data, mask=mask, weight=None,
                       transpose=False, inittau=0.02, niter=100, dist="e")
        message = "^incorrect rank 3 \\(expected 1\\)$"
        with self.assertRaisesRegex(ValueError, message):
            somcluster(clusterids=clusterids, celldata=celldata,
                       data=data, mask=mask, weight=numpy.zeros((2, 2, 2)),
                       transpose=False, inittau=0.02, niter=100, dist="e")
        message = "^array has incorrect data type$"
        with self.assertRaisesRegex(RuntimeError, message):
            somcluster(clusterids=clusterids, celldata=celldata,
                       data=data, mask=mask,
                       weight=numpy.array([1, 1, 1], dtype=numpy.int16),
                       transpose=False, inittau=0.02, niter=100, dist="e")
        message = "^dist should be a string$"
        with self.assertRaisesRegex(ValueError, message):
            somcluster(clusterids=clusterids, celldata=celldata,
                       data=data, mask=mask, weight=weight,
                       transpose=False, inittau=0.02, niter=100, dist=5)
        message = "^dist should be a single character$"
        with self.assertRaisesRegex(ValueError, message):
            somcluster(clusterids=clusterids, celldata=celldata,
                       data=data, mask=mask, weight=weight,
                       transpose=False, inittau=0.02, niter=100,
                       dist="Pearson")
        message = "^unknown dist function specified \\(should be one of 'ebcauxsk'\\)$"
        with self.assertRaisesRegex(ValueError, message):
            somcluster(clusterids=clusterids, celldata=celldata,
                       data=data, mask=mask, weight=weight,
                       transpose=False, inittau=0.02, niter=100, dist="X")

    def test_somcluster(self):
        if TestCluster.module == "Bio.Cluster":
            from Bio.Cluster import somcluster
        elif TestCluster.module == "Pycluster":
            from Pycluster import somcluster

        # First data set
        weight = [1, 1, 1, 1, 1]
        data = numpy.array([[1.1, 2.2, 3.3, 4.4, 5.5],
                            [3.1, 3.2, 1.3, 2.4, 1.5],
                            [4.1, 2.2, 0.3, 5.4, 0.5],
                            [9.9, 2.0, 0.0, 5.0, 0.0]])
        mask = numpy.array([[1, 1, 1, 1, 1],
                            [1, 1, 1, 1, 1],
                            [1, 1, 1, 1, 1],
                            [1, 1, 1, 1, 1]], int)
        nrows, ncols = data.shape

        clusterid, celldata = somcluster(data=data, mask=mask, weight=weight,
                                         transpose=False, nxgrid=10, nygrid=10,
                                         inittau=0.02, niter=100, dist="e")
        self.assertEqual(len(clusterid), nrows)
        self.assertEqual(len(clusterid[0]), 2)

        # First data set, using transpose=True
        weight = [1, 1, 1, 1]
        clusterid, celldata = somcluster(data=data, mask=mask, weight=weight,
                                         transpose=True, nxgrid=10, nygrid=10,
                                         inittau=0.02, niter=100, dist="e")
        self.assertEqual(len(clusterid), ncols)
        self.assertEqual(len(clusterid[0]), 2)

        # Second data set
        weight = [1, 1]
        data = numpy.array([[1.1, 1.2],
                            [1.4, 1.3],
                            [1.1, 1.5],
                            [2.0, 1.5],
                            [1.7, 1.9],
                            [1.7, 1.9],
                            [5.7, 5.9],
                            [5.7, 5.9],
                            [3.1, 3.3],
                            [5.4, 5.3],
                            [5.1, 5.5],
                            [5.0, 5.5],
                            [5.1, 5.2]])
        mask = numpy.array([[1, 1],
                            [1, 1],
                            [1, 1],
                            [1, 1],
                            [1, 1],
                            [1, 1],
                            [1, 1],
                            [1, 1],
                            [1, 1],
                            [1, 1],
                            [1, 1],
                            [1, 1],
                            [1, 1]], int)

        clusterid, celldata = somcluster(data=data, mask=mask, weight=weight,
                                         transpose=False, nxgrid=10, nygrid=10,
                                         inittau=0.02, niter=100, dist="e")
        self.assertEqual(len(clusterid), len(data))
        self.assertEqual(len(clusterid[0]), 2)

    def test_distancematrix_arguments(self):
        # Test if incorrect arguments are caught by the C code
        if TestCluster.module == "Bio.Cluster":
            from Bio.Cluster._cluster import distancematrix
        elif TestCluster.module == "Pycluster":
            from Pycluster._cluster import distancematrix

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
                            [1, 1, 1]], numpy.int32)
        weight = numpy.array([2.0, 1.0, 0.5])
        message = "^data matrix is empty$"
        with self.assertRaisesRegex(ValueError, message):
            distancematrix(data[:0, :], mask=mask, weight=weight)
        message = "^mask has incorrect rank 1 \\(expected 2\\)$"
        with self.assertRaisesRegex(ValueError, message):
            distancematrix(data, mask=numpy.zeros(3), weight=weight)
        message = "^mask has incorrect dimensions \\(4 x 3, expected 9 x 3\\)$"
        with self.assertRaisesRegex(ValueError, message):
            distancematrix(data, mask=mask[:4, :], weight=weight,
                           transpose=False, dist="c",
                           distancematrix=[])
        message = "^incorrect rank 2 \\(expected 1\\)$"
        with self.assertRaisesRegex(ValueError, message):
            distancematrix(data, mask=mask, weight=numpy.zeros((2, 2)))
        message = "^weight has incorrect size 4 \\(expected 3\\)$"
        with self.assertRaisesRegex(ValueError, message):
            distancematrix(data, mask=mask, weight=numpy.zeros(4),
                           transpose=False, dist="c", distancematrix=[])

    def test_kmedoids_arguments(self):
        # Test if incorrect arguments are caught by the C code
        if TestCluster.module == "Bio.Cluster":
            from Bio.Cluster._cluster import distancematrix, kmedoids
        elif TestCluster.module == "Pycluster":
            from Pycluster._cluster import distancematrix, kmedoids

        clusterid = numpy.zeros(10, numpy.int32)
        message = "^failed to parse row 0.$"
        with self.assertRaisesRegex(RuntimeError, message):
            kmedoids([None])
        message = "^more clusters requested than items to be clustered$"
        with self.assertRaisesRegex(ValueError, message):
            kmedoids([], nclusters=2, npass=1000, clusterid=clusterid)
        message = "^distance matrix is not square.$"
        with self.assertRaisesRegex(ValueError, message):
            kmedoids(numpy.zeros((2, 3)), npass=1000)
        message = "^distance matrix has incorrect rank 3 \\(expected 1 or 2\\)$"
        with self.assertRaisesRegex(ValueError, message):
            kmedoids(numpy.zeros((2, 3, 4)), npass=1000)

    def test_distancematrix_kmedoids(self):
        if TestCluster.module == "Bio.Cluster":
            from Bio.Cluster import distancematrix, kmedoids
        elif TestCluster.module == "Pycluster":
            from Pycluster import distancematrix, kmedoids

        # transpose=False
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

        # check if default weights can be used
        matrix = distancematrix(data, mask=mask)
        self.assertEqual(len(matrix), 9)
        for i in range(3):
            self.assertEqual(len(matrix[i]), i)

        self.assertAlmostEqual(matrix[1][0], 1.687, places=3)

        self.assertAlmostEqual(matrix[2][0], 21.365, places=3)
        self.assertAlmostEqual(matrix[2][1], 38.560, places=3)

        self.assertAlmostEqual(matrix[3][0], 4.900, places=3)
        self.assertAlmostEqual(matrix[3][1], 7.793, places=3)
        self.assertAlmostEqual(matrix[3][2], 22.490, places=3)

        self.assertAlmostEqual(matrix[4][0], 3.687, places=3)
        self.assertAlmostEqual(matrix[4][1], 6.367, places=3)
        self.assertAlmostEqual(matrix[4][2], 22.025, places=3)
        self.assertAlmostEqual(matrix[4][3], 0.087, places=3)

        self.assertAlmostEqual(matrix[5][0], 0.040, places=3)
        self.assertAlmostEqual(matrix[5][1], 2.890, places=3)
        self.assertAlmostEqual(matrix[5][2], 34.810, places=3)
        self.assertAlmostEqual(matrix[5][3], 0.640, places=3)
        self.assertAlmostEqual(matrix[5][4], 0.490, places=3)

        self.assertAlmostEqual(matrix[6][0], 1.557, places=3)
        self.assertAlmostEqual(matrix[6][1], 0.990, places=3)
        self.assertAlmostEqual(matrix[6][2], 34.065, places=3)
        self.assertAlmostEqual(matrix[6][3], 3.937, places=3)
        self.assertAlmostEqual(matrix[6][4], 3.017, places=3)
        self.assertAlmostEqual(matrix[6][5], 3.610, places=3)

        self.assertAlmostEqual(matrix[7][0], 14.005, places=3)
        self.assertAlmostEqual(matrix[7][1], 9.050, places=3)
        self.assertAlmostEqual(matrix[7][2], 65.610, places=3)
        self.assertAlmostEqual(matrix[7][3], 30.465, places=3)
        self.assertAlmostEqual(matrix[7][4], 27.380, places=3)
        self.assertAlmostEqual(matrix[7][5], 0.000, places=3)
        self.assertAlmostEqual(matrix[7][6], 16.385, places=3)

        self.assertAlmostEqual(matrix[8][0], 14.167, places=3)
        self.assertAlmostEqual(matrix[8][1], 25.553, places=3)
        self.assertAlmostEqual(matrix[8][2], 0.010, places=3)
        self.assertAlmostEqual(matrix[8][3], 17.187, places=3)
        self.assertAlmostEqual(matrix[8][4], 16.380, places=3)
        self.assertAlmostEqual(matrix[8][5], 33.640, places=3)
        self.assertAlmostEqual(matrix[8][6], 22.497, places=3)
        self.assertAlmostEqual(matrix[8][7], 36.745, places=3)

        # transpose=True
        weight = numpy.array([2.0, 1.0, 0.5, 0.1, 0.9, 3.0, 2.0, 1.5, 0.2])
        matrix = distancematrix(data, mask=mask, weight=weight, transpose=True)
        self.assertEqual(len(matrix), 3)
        for i in range(3):
            self.assertEqual(len(matrix[i]), i)

        self.assertAlmostEqual(matrix[1][0], 3.080323, places=3)
        self.assertAlmostEqual(matrix[2][0], 9.324416, places=3)
        self.assertAlmostEqual(matrix[2][1], 11.569701, places=3)

        clusterid, error, nfound = kmedoids(matrix, npass=1000)
        self.assertEqual(clusterid[0], 0)
        self.assertEqual(clusterid[1], 0)
        self.assertEqual(clusterid[2], 2)
        self.assertAlmostEqual(error, 3.08032258, places=3)

        # check if default weights can be used
        matrix = distancematrix(data, mask=mask, transpose=True)
        self.assertEqual(len(matrix), 3)
        for i in range(3):
            self.assertEqual(len(matrix[i]), i)

        self.assertAlmostEqual(matrix[1][0], 10.47166667, places=3)
        self.assertAlmostEqual(matrix[2][0], 8.61571429, places=3)
        self.assertAlmostEqual(matrix[2][1], 21.24428571, places=3)

    def test_pca_arguments(self):
        if TestCluster.module == "Bio.Cluster":
            from Bio.Cluster._cluster import pca
        elif TestCluster.module == "Pycluster":
            from Pycluster._cluster import pca

        data = numpy.zeros((4, 2))
        columnmean = numpy.zeros(2)
        pc = numpy.zeros((2, 2), dtype="d")
        coordinates = numpy.zeros((4, 2), dtype="d")
        eigenvalues = numpy.zeros(2, dtype="d")

        message = "^data matrix has unexpected format.$"
        with self.assertRaisesRegex(RuntimeError, message):
            pca([None], columnmean, coordinates, pc, eigenvalues)
        message = "^data matrix has incorrect rank 1 \\(expected 2\\)$"
        with self.assertRaisesRegex(RuntimeError, message):
            pca(numpy.zeros(3), columnmean, coordinates, pc, eigenvalues)
        message = "^data matrix has incorrect data type$"
        with self.assertRaisesRegex(RuntimeError, message):
            pca(numpy.zeros((3, 3), dtype=numpy.int16),
                columnmean, coordinates, pc, eigenvalues)
        message = "^unexpected format.$"
        with self.assertRaisesRegex(RuntimeError, message):
            pca(data, "nothing", coordinates, pc, eigenvalues)
        message = "^incorrect rank 2 \\(expected 1\\)$"
        with self.assertRaisesRegex(ValueError, message):
            pca(data, numpy.zeros((2, 2)), coordinates, pc, eigenvalues)
        message = "^array has incorrect data type$"
        with self.assertRaisesRegex(RuntimeError, message):
            pca(data, numpy.ones(3, dtype=numpy.int16),
                coordinates, pc, eigenvalues)
        message = "^data matrix has unexpected format.$"
        with self.assertRaisesRegex(RuntimeError, message):
            pca(data, columnmean, [None], pc, eigenvalues)
        message = "^data matrix has incorrect rank 1 \\(expected 2\\)$"
        with self.assertRaisesRegex(RuntimeError, message):
            pca(data, columnmean, numpy.zeros(3), pc, eigenvalues)
        message = "^data matrix has incorrect data type$"
        with self.assertRaisesRegex(RuntimeError, message):
            pca(data, columnmean, numpy.zeros((3, 3), dtype=numpy.int16),
                pc, eigenvalues)
        message = "^data matrix has unexpected format.$"
        with self.assertRaisesRegex(RuntimeError, message):
            pca(data, columnmean, coordinates, [None], eigenvalues)
        message = "^data matrix has incorrect rank 1 \\(expected 2\\)$"
        with self.assertRaisesRegex(RuntimeError, message):
            pca(data, columnmean, coordinates, numpy.zeros(3), eigenvalues)
        message = "^data matrix has incorrect data type$"
        with self.assertRaisesRegex(RuntimeError, message):
            pca(data, columnmean, coordinates,
                numpy.zeros((3, 3), dtype=numpy.int16), eigenvalues)
        message = "^unexpected format.$"
        with self.assertRaisesRegex(RuntimeError, message):
            pca(data, columnmean, coordinates, pc, "nothing")
        message = "^incorrect rank 2 \\(expected 1\\)$"
        with self.assertRaisesRegex(ValueError, message):
            pca(data, columnmean, coordinates, pc, numpy.zeros((2, 2)))
        message = "^array has incorrect data type$"
        with self.assertRaisesRegex(RuntimeError, message):
            pca(data, columnmean, coordinates, pc,
                numpy.ones(3, dtype=numpy.int16))

    def test_pca(self):
        if TestCluster.module == "Bio.Cluster":
            from Bio.Cluster import pca
        elif TestCluster.module == "Pycluster":
            from Pycluster import pca

        data = numpy.array([[3.1, 1.2],
                            [1.4, 1.3],
                            [1.1, 1.5],
                            [2.0, 1.5],
                            [1.7, 1.9],
                            [1.7, 1.9],
                            [5.7, 5.9],
                            [5.7, 5.9],
                            [3.1, 3.3],
                            [5.4, 5.3],
                            [5.1, 5.5],
                            [5.0, 5.5],
                            [5.1, 5.2],
                            ])

        mean, coordinates, pc, eigenvalues = pca(data)
        self.assertAlmostEqual(mean[0], 3.5461538461538464)
        self.assertAlmostEqual(mean[1], 3.5307692307692311)
        self.assertAlmostEqual(coordinates[0, 0], 2.0323189722653883)
        self.assertAlmostEqual(coordinates[0, 1], 1.2252420399694917)
        self.assertAlmostEqual(coordinates[1, 0], 3.0936985166252251)
        self.assertAlmostEqual(coordinates[1, 1], -0.10647619705157851)
        self.assertAlmostEqual(coordinates[2, 0], 3.1453186907749426)
        self.assertAlmostEqual(coordinates[2, 1], -0.46331699855941139)
        self.assertAlmostEqual(coordinates[3, 0], 2.5440202962223761)
        self.assertAlmostEqual(coordinates[3, 1], 0.20633980959571077)
        self.assertAlmostEqual(coordinates[4, 0], 2.4468278463376221)
        self.assertAlmostEqual(coordinates[4, 1], -0.28412285736824866)
        self.assertAlmostEqual(coordinates[5, 0], 2.4468278463376221)
        self.assertAlmostEqual(coordinates[5, 1], -0.28412285736824866)
        self.assertAlmostEqual(coordinates[6, 0], -3.2018619434743254)
        self.assertAlmostEqual(coordinates[6, 1], 0.019692314198662915)
        self.assertAlmostEqual(coordinates[7, 0], -3.2018619434743254)
        self.assertAlmostEqual(coordinates[7, 1], 0.019692314198662915)
        self.assertAlmostEqual(coordinates[8, 0], 0.46978641990344067)
        self.assertAlmostEqual(coordinates[8, 1], -0.17778754731982949)
        self.assertAlmostEqual(coordinates[9, 0], -2.5549912731867215)
        self.assertAlmostEqual(coordinates[9, 1], 0.19733897451533403)
        self.assertAlmostEqual(coordinates[10, 0], -2.5033710990370044)
        self.assertAlmostEqual(coordinates[10, 1], -0.15950182699250004)
        self.assertAlmostEqual(coordinates[11, 0], -2.4365601663089413)
        self.assertAlmostEqual(coordinates[11, 1], -0.23390813900973562)
        self.assertAlmostEqual(coordinates[12, 0], -2.2801521629852974)
        self.assertAlmostEqual(coordinates[12, 1], 0.0409309711916888)
        self.assertAlmostEqual(pc[0, 0], -0.66810932728062988)
        self.assertAlmostEqual(pc[0, 1], -0.74406312017235743)
        self.assertAlmostEqual(pc[1, 0], 0.74406312017235743)
        self.assertAlmostEqual(pc[1, 1], -0.66810932728062988)
        self.assertAlmostEqual(eigenvalues[0], 9.3110471246032844)
        self.assertAlmostEqual(eigenvalues[1], 1.4437456297481428)

        data = numpy.array([[2.3, 4.5, 1.2, 6.7, 5.3, 7.1],
                            [1.3, 6.5, 2.2, 5.7, 6.2, 9.1],
                            [3.2, 7.2, 3.2, 7.4, 7.3, 8.9],
                            [4.2, 5.2, 9.2, 4.4, 6.3, 7.2]])
        mean, coordinates, pc, eigenvalues = pca(data)
        self.assertAlmostEqual(mean[0], 2.7500)
        self.assertAlmostEqual(mean[1], 5.8500)
        self.assertAlmostEqual(mean[2], 3.9500)
        self.assertAlmostEqual(mean[3], 6.0500)
        self.assertAlmostEqual(mean[4], 6.2750)
        self.assertAlmostEqual(mean[5], 8.0750)
        self.assertAlmostEqual(coordinates[0, 0], 2.6460846688406905)
        self.assertAlmostEqual(coordinates[0, 1], -2.1421701432732418)
        self.assertAlmostEqual(coordinates[0, 2], -0.56620932754145858)
        self.assertAlmostEqual(coordinates[0, 3], 0.0)
        self.assertAlmostEqual(coordinates[1, 0], 2.0644120899917544)
        self.assertAlmostEqual(coordinates[1, 1], 0.55542108669180323)
        self.assertAlmostEqual(coordinates[1, 2], 1.4818772348457117)
        self.assertAlmostEqual(coordinates[1, 3], 0.0)
        self.assertAlmostEqual(coordinates[2, 0], 1.0686641862092987)
        self.assertAlmostEqual(coordinates[2, 1], 1.9994412069101073)
        self.assertAlmostEqual(coordinates[2, 2], -1.000720598980291)
        self.assertAlmostEqual(coordinates[2, 3], 0.0)
        self.assertAlmostEqual(coordinates[3, 0], -5.77916094504174)
        self.assertAlmostEqual(coordinates[3, 1], -0.41269215032867046)
        self.assertAlmostEqual(coordinates[3, 2], 0.085052691676038017)
        self.assertAlmostEqual(coordinates[3, 3], 0.0)
        self.assertAlmostEqual(pc[0, 0], -0.26379660005997291)
        self.assertAlmostEqual(pc[0, 1], 0.064814972617134495)
        self.assertAlmostEqual(pc[0, 2], -0.91763310094893846)
        self.assertAlmostEqual(pc[0, 3], 0.26145408875373249)
        self.assertAlmostEqual(pc[1, 0], 0.05073770520434398)
        self.assertAlmostEqual(pc[1, 1], 0.68616983388698793)
        self.assertAlmostEqual(pc[1, 2], 0.13819106187213354)
        self.assertAlmostEqual(pc[1, 3], 0.19782544121828985)
        self.assertAlmostEqual(pc[2, 0], -0.63000893660095947)
        self.assertAlmostEqual(pc[2, 1], 0.091155993862151397)
        self.assertAlmostEqual(pc[2, 2], 0.045630391256086845)
        self.assertAlmostEqual(pc[2, 3], -0.67456694780914772)
        # As the last eigenvalue is zero, the corresponding eigenvector is
        # strongly affected by roundoff error, and is not being tested here.
        # For PCA, this doesn't matter since all data have a zero coefficient
        # along this eigenvector.
        self.assertAlmostEqual(eigenvalues[0], 6.7678878332578778)
        self.assertAlmostEqual(eigenvalues[1], 3.0108911400291856)
        self.assertAlmostEqual(eigenvalues[2], 1.8775592718563467)
        self.assertAlmostEqual(eigenvalues[3], 0.0)


if __name__ == "__main__":
    TestCluster.module = "Bio.Cluster"
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
