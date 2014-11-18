#!/usr/bin/env python
"""Tests for testing that clustering works as expected."""

# Copyright 2015 by DV Klopfenstein.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This test and test infrastructure was created to help verify this issue:
#
#   https://github.com/biopython/biopython/issues/424#issuecomment-62823368

import unittest

import sys
import copy
import collections as cx
from operator import itemgetter

try:
    from StringIO import StringIO # Python 2
    # Can't use StringIO, quoting the documentation,
    #   "Unlike the StringIO module, this module is not able to accept
    #    Unicode strings that cannot be encoded as plain ASCII strings."
    # Therefore can't use from Bio._py3k import StringIO
except ImportError:
    from io import StringIO # Python 3

try:
    import numpy as np
except ImportError:
    from Bio import MissingExternalDependencyError
    raise MissingExternalDependencyError(
        "Install NumPy if you want to use Bio.KDTree.")

try:
    from Bio import Cluster
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        "If you want to use Bio.Cluster, "
        "install NumPy first and then "
        "reinstall Biopython")

#############################################################################
# User-Created Tests: Edit this area to add more tests
#############################################################################
class TestClusterSmall(unittest.TestCase):
    """Unit tests; hierarchical clustering of small sets of data."""

    # Test input (Fails for Issue #440; method=('m', 's')
    i440_in = (" A B  C D     E     G H I K L     M N       O P Q    R")
    # Some expected clusters for input data, i440_in:
    i440_exp = [
        "A B",
        "C D",
        "A B  C D",
        # To (mock) "PASS" method='m': COMMENT OUT FOLLOWING LINE OUT
        "G H I K L",
        "M N",
        "O P Q"]

    def test_treecluster_i440_m(self):
        """Uncomment to see failure for issue biopython#424.

        ISSUE: Is it expected that cluster "G H I K L" is broken into
               "G H" and "I K L"?:

        ACTUAL CLUSTERS (method=m), 2nd FORMAT:
          Index 10s: 0         1         2         3         4         5
          Index  1s: 012345678901234567890123456789012345678901234567890123
          Array    :  A B  C D     E     G H I K L     M N       O P Q    R
            Level  1  <--------------------------------------------------->
            Level  2  <--------------------> <---------------------------->
            Level  3  <------>     <-------> <----------->       <-------->
            Level  4  <->  <->           <-> <--->     <->       <--->
            Level  5                         <->                   <->
         
        """
        #   'm' FAIL pairwise maximum- (or complete-) linkage clustering
        # Uncomment to see failure for issue biopython#424.
        #self.doit(TestClusterSmall.i440_in, TestClusterSmall.i440_exp, 'm')
        pass

    def test_treecluster_i440_s(self):
        """Uncomment to see failure for issue biopython#424.

        ISSUE: Is it expected that there are nested clusters? For example:

        1. Cluster "O    Q" is a child cluster of "O P Q".
           Would it bemore expected that "O P" or "P Q" would be a child
           cluster of "O P Q"?

        2. Is it expected that cluster "G H I K L" and cluster "A B C D E M N"
           are both on the same hierarchical level(3), when "G H I K L"
           is between "A B C D E" and "M N"?:

        HIERARCHY DIAGRAM (method=s), 1st FORMAT:
           i Hier   NodeID n.left n.right n.dist leaf-values
         --- ------ ------ ------ ------- ------ -----------
           0 -         -15 ( -14,   -11):   8.00 A B C D E G H I K L M N O P Q R
           1 --        -14 (  -9,   -13):   6.00 A B C D E G H I K L M N
           2 ---        -9 (   G,    -4):   2.00 G H I K L
           3 ----       -4 (   K,    -3):   2.00 H I K L
           4 -----      -3 (   I,    -2):   2.00 H I L
           5 ------     -2 (   H,     L):   2.00 H L
           6 ---       -13 (   E,   -12):   6.00 A B C D E M N
           7 ----      -12 ( -10,    -5):   6.00 A B C D M N
           8 -----     -10 (  -1,    -8):   3.00 A B C D
           9 ------     -1 (   A,     B):   2.00 A B
          10 ------     -8 (   C,     D):   2.00 C D
          11 -----      -5 (   M,     N):   2.00 M N
          12 --        -11 (  -7,     R):   5.00 O P Q R
          13 ---        -7 (   P,    -6):   2.00 O P Q
          14 ----       -6 (   O,     Q):   2.00 O Q
      
        HIERARCHY DIAGRAM (method=s), 2nd FORMAT:
          Index 10s: 0         1         2         3         4         5
          Index  1s: 012345678901234567890123456789012345678901234567890123
          Array    :  A B  C D     E     G H I K L     M N       O P Q    R
            Level  1  <--------------------------------------------------->
            Level  2  <---------------------------------->       <-------->
            Level  3  <------------------<------->------->       <--->
            Level  4  <--------------------<----->------->       <--->
            Level  5  <------>             <----->     <->
            Level  6  <->  <->             <----->

        """
        #   's' FAIL pairwise single-linkage clustering
        # Uncomment to see failure for issue biopython#424.
        #self.doit(TestClusterSmall.i440_in, TestClusterSmall.i440_exp, 's')
        pass

    def test_treecluster_i440_a(self):
        """This is a user-created test issue biopython#424."""
        #   method: 'a' pairwise average-linkage clustering
        self.doit(TestClusterSmall.i440_in, TestClusterSmall.i440_exp, 'a')
        pass


    def doit(self, usr_input, expected_clusters, method):
        test_obj = ClusterTestHelper(usr_input, method)
        test_obj.chk_cluster(expected_clusters)
        test_obj.tree_test_obj.prt_hier()
        test_obj.prt_ascii_art_levels()


#############################################################################
# TEST INFRASTRUCTURE: Do not edit when creating additional tests
#############################################################################
class ClustersStrInput:
    """Class for storing and printing clusters in a linear format."""

    def __init__(self, clusters_str):
        """Store a string representing clusters and its helper functions."""
        self.clu_str = clusters_str # e.g. "ABC   DE F  GH"
        self.names, self.locs = self.get_names_locs()

    def get_names_locs(self):
        """Given "A B   C  D", returns names and locations.

            names = ('A', 'B', 'C', 'D')
            locs  = ( 0 ,  1 ,  2 ,  3 )
        """
        names_locs = [(Name, Loc)
                      for Loc, Name in enumerate(self.clu_str)
                      if Name is not " "]
        return zip(*names_locs)

    def prt_arrays(self, prt=sys.stdout):
        """Prints the variables."""
        rng = range(len(self.names))
        index = ['  '.join('{0:>2}'.format(E) for E in rng)]
        names = [', '.join('{0:>2}'.format(E) for E in self.names)]
        locas = [', '.join('{0:>2}'.format(E) for E in self.locs)]
        prt.write('\n  VARIABLES:\n')
        prt.write('    index :   {0}\n'.format(index))
        prt.write('    names =   {0}\n'.format(names))
        prt.write('    locs  =   {0}\n'.format(locas))

    def prt_ascii_art(self, prt=sys.stdout, title="ASCII Art"):
        """Prints the User-Input in a human-readable way."""
        rng = range(len(self.clu_str))
        i10s = ''.join('{0}'.format(I/10) if I%10 == 0 else ' ' for I in rng)
        i01s = ''.join('{0}'.format(I%10) for I in rng)
        prt.write('\n  {0}:\n'.format(title))
        prt.write('    Index 10s: {0}\n'.format(i10s))
        prt.write('    Index  1s: {0}\n'.format(i01s))
        prt.write('    Array    : {0}\n'.format(self.clu_str))

    def __str__(self):
        sout = StringIO()
        self.prt_ascii_art(sout)
        self.prt_arrays(sout)
        ret_str = sout.getvalue()
        sout.close()
        return ret_str


class ClusterTree:
    """Used to create, store and print Biopython Cluster treecluster results."""

    def __init__(self, dist_matrix, method, names):
        #       meth  P/F
        # default 'm' FAIL pairwise maximum- (or complete-) linkage clustering
        #         's' NO!  pairwise single-linkage clustering
        #         'a' No   pairwise average-linkage clustering
        self.tree = Cluster.treecluster(distancematrix=dist_matrix, 
            dist='e', # Euclidean distance is the Bio.Cluster.treecluster default 
            method=method)
        # Reformat tree for easier printing in various formats
        # The following 3 values are set by set_hier_list():
        self.max_level = None
        self.hier_list = []
        self.actual_clusters = [] # Clusters found by treecluster
        self.set_hier_info(names)

    def set_hier_info(self, names):
        """Get the hierarchy."""
        # Assume top-level node is the last node in the treecluster list
        self._set_hier_list(names, -1*len(self.tree), level=1)
        self.max_level = max(elem[0] for elem in self.hier_list)
        # I        0   1      2    3  4       5  6        7     8
        # elem Level, idx, NodeID, L, L.left, R, R.right, D, num_leafs
        for elem in self.hier_list:
            clu_set = frozenset(self.get_all_locs_or_names(names, elem[4], elem[6]))
            self.actual_clusters.append(clu_set)
            elem.append(len(clu_set))

    def prt(self, prt=sys.stdout):
        """Print Tree without printing node names."""
        self.prt_w_names(None, prt)

    def prt_w_names(self, names, prt=sys.stdout):
        """Print Tree using node names."""
        prt.write('\n  TREE ARRAY FROM Bio.Cluster.treecluster:\n')
        prt.write('               NodeID (left, right): distance\n')
        for idx, elem in enumerate(self.tree):
            ClusterTree.prt_node(idx, elem, names, prt)

    @classmethod
    def prt_node(cls, idx, node, names, prt):
        """Print a node from Bio.Cluster.treecluster and additional info."""
        left = ClusterTree.get_node_str(node.left, names)
        right = ClusterTree.get_node_str(node.right, names)
        prt.write('    tree[{0:2}] = {1:6} ({2:>4}, {3:>5}): {4:>6}\n'.format(
            idx, -idx-1, left, right, node.distance))

    @classmethod
    def get_node_str(cls, lhs_rhs, names):
        """If available, return leaf node name."""
        return lhs_rhs if lhs_rhs < 0 or names is None else names[lhs_rhs]

    def prt_hier(self, prt=sys.stdout, title="HIERARCHY DIAGRAM"):
        """Print the hierarchy.

           How to read this hierarchy report:

             HIERARCHY DIAGRAM, 1st FORMAT:
                0 -        -15 ( -14,   -13):   51.0
                1 --       -14 ( -12,   -11):   28.0
                2 ---      -12 (  -8,    -6):   12.0
                3 ----      -8 (  -5,     L):    4.0
                4 -----     -5 (   K,     I):    2.0
                5 ----      -6 (   M,     N):    2.0
                6 ---      -11 (  -7,     R):    9.0
                7 ----      -7 (   O,    -4):    4.0
                8 -----     -4 (   P,     Q):    2.0
                9 --       -13 ( -10,    -9):   21.0
               10 ---      -10 (  -3,     E):    8.0
               11 ----      -3 (   H,     G):    2.0
               12 ---       -9 (  -2,    -1):    7.0
               13 ----      -2 (   D,     C):    2.0
               14 ----      -1 (   B,     A):    2.0

           The TOP-LEVEL NODE is id = -15:
                0 -        -15 ( -14,   -13):   51.0

           The TOP-LEVEL NODE CONTAINS TWO CLUSTER NODES:
                0 -        -15 ( -14,   -13):   51.0
                1 --       -14 ( -12,   -11):   28.0
                9 --       -13 ( -10,    -9):   21.0

           PATH FROM TOP TO LEAF NODES(A,B and C,D) on the LEFT-MOST CLUSTER:
                0 -        -15 ( -14,   -13):   51.0
                9 --       -13 ( -10,    -9):   21.0
               12 ---       -9 (  -2,    -1):    7.0
               13 ----      -2 (   D,     C):    2.0
               14 ----      -1 (   B,     A):    2.0

        """
        prt.write('\n  {0}:\n'.format(title))
        max_lev = self.max_level
        hier = '{0:{1}}'.format("Hier", max_lev)
        prt.write('     i {0} NodeID n.left n.right n.dist leaf-values\n'.format(hier))
        prt.write('   --- {0} ------ ------ ------- ------ -----------\n'.format('-'*max_lev))
        for idx, elem in enumerate(self.hier_list):
            hier_dashes = '-'*elem[0]
            hier = '{0:{1}}'.format(hier_dashes, max_lev)
            vals = ' '.join(sorted(self.actual_clusters[idx]))
            pat = '    {0:2} {1} {2:6} ({3:>4}, {4:>5}): {5:>6.2f} {6}\n'
            prt.write(pat.format(idx, hier, elem[2], elem[3], elem[5], elem[7], vals))

    def _set_hier_list(self, names, node_num, level):
        """Recursively called get_hierarchy."""
        # If this is not a leaf node
        if node_num < 0:
            idx = -1*node_num - 1
            node = self.tree[idx]
            # Print the intermediate node
            left = self.get_node_str(node.left, names)
            right = self.get_node_str(node.right, names)
            # Level, idx, NodeID, L, L.left, R, R.right, D
            elems = [level, idx, node_num, left,
                     node.left, right, node.right, node.distance]
            self.hier_list.append(elems)
            self._set_hier_list(names, node.left, level+1)
            self._set_hier_list(names, node.right, level+1)

    def __str__(self):
        sout = StringIO()
        self.prt(sout)
        self.prt_hier(sout)
        ret_str = sout.getvalue()
        sout.close()
        return ret_str

    def get_all_locs_or_names(self, all_locs, left, right):
        """Return either all the locations or all the names in a cluster."""
        cur_locs = []
        self._get_all_locs_or_names(cur_locs, all_locs, left)
        self._get_all_locs_or_names(cur_locs, all_locs, right)
        return cur_locs

    def _get_all_locs_or_names(self, cur_locs, all_locs, nodeid_locidx):
        """Internal method used by public method."""
        # If this is a leaf node, append cur_loc list
        if nodeid_locidx >= 0:
            cur_locs.append(all_locs[nodeid_locidx])
            return
        # If this is not intermediate node, keep drilling
        node = self.tree[-1*nodeid_locidx - 1]
        self._get_all_locs_or_names(cur_locs, all_locs, node.left)
        self._get_all_locs_or_names(cur_locs, all_locs, node.right)


class ClusterTestHelper:
    """Class used to take items in a string and generate clusters."""

    def __init__(self, elem_str, method='m'):
        self.method = method
        self.usr_input_obj = ClustersStrInput(elem_str)
        self.dist_matrix_orig = self._mk_dist_matrix()
        warn = "**WARNING: Bio.Cluster.treecluster " \
               "CHANGES THE DISTANCE MATRIX, SO MAKE A COPY\n"
        #sys.stdout.write(warn) # This is expected behavior
        self.dist_matrix = copy.deepcopy(self.dist_matrix_orig)
        #print self.dist_matrix_orig.dist_matrix
        #print self.dist_matrix
        self.tree_test_obj = ClusterTree(
            self.dist_matrix, method, self.usr_input_obj.names)
        #print self.dist_matrix
        #print self.dist_matrix_orig.dist_matrix
        self.level_pts = cx.defaultdict(list)
        self.set_cluster_pts()
        self.chk_cluster_overlap()

    def prt_dist_matrix(self, names, prt=sys.stdout):
        """Prints the distance matrix in a human-readable format."""
        cnt = 1
        cnt_stop = 1
        prt.write('\n  dist_matrix = [\n')
        prt.write('   # {0}\n  '.format(names[0]))
        for elem in self.dist_matrix_orig:
            prt.write('  {0:4.1f}, '.format(elem))
            if cnt == cnt_stop:
                prt.write(' # {0}\n  '.format(names[cnt_stop]))
                cnt = 0
                cnt_stop += 1
            cnt += 1
        prt.write('  ]\n')

    def prt_ascii_art_levels(self, prt=sys.stdout, title="HIERARCHY DIAGRAM"):
        """Prints ASCII Art representing a hierarchy diagram.

           GIVEN INPUT:  "A B  C D     E     G H I K L     M N      O P Q    R"

           GIVEN THIS TREE FROM Bio.Cluster.treecluster, THE TREE IS PRINTED:

             HIERARCHY DIAGRAM, 2nd FORMAT:
               Index 10s: 0         1         2         3         4         5
               Index  1s: 01234567890123456789012345678901234567890123456789012
               Array    :  A B  C D     E     G H I K L     M N      O P Q    R
                 Level  1  ####################################################
                 Level  2  ###################### #############################
                 Level  3  ########     ######### #############      ##########
                 Level  4  ###  ###           ### #####     ###      #####
                 Level  5                         ###                  ###

           THE ABOVE TREE FORMAT IS EASY TO CREATE usr_input_obj A SMALL AMOUNT
           OF CODE AND IS EQUIVALENT TO THIS MORE FAMILIAR FORMAT:


             #                                    |
             # Level 1               +------------a-------------+
             #                       |                          |
             # Level 2        +------b-----+            +-------c--------+
             #                |            |            |                |
             # Level 3     +--d-+      +---e--+     +---f---+         +--g---+
             #             |    |      |      |     |       |         |      |
             # Level 4    +h+  +i+     |     +k+  +-l+     +m+      +-n+     |
             #            | |  | |     |     | |  |  |     | |      |  |     |
             # Level 5    | |  | |     |     | | +o+ |     | |      | +p+    |
             #            | |  | |     |     | | | | |     | |      | | |    |
             do_cluster(" A B  C D     E     G H I K L     M N      O P Q    R")

        """
        # Sort hierarchy list by Hierarchy Level where 1 is the Top-Level
        self.usr_input_obj.prt_ascii_art(prt, title)
        num_locs = max(self.usr_input_obj.locs) + 1
        for level, clu_pts in sorted(self.level_pts.items()):
            ret_str = [' ']*num_locs
            for pts in clu_pts:
                clu_min = pts[0]
                clu_max = pts[1]
                ret_str[clu_min] = '<'
                ret_str[clu_max] = '>'
                for i in range(clu_min+1, clu_max):
                    if ret_str[i] == ' ':
                        ret_str[i] = '-'
            prt.write('      Level {0:2} {1}\n'.format(level, ''.join(ret_str)))

    def set_cluster_pts(self):
        """Creates lists of cluster endpoints for each level of hierarchy."""
        locs = self.usr_input_obj.locs
        for elem in self.tree_test_obj.hier_list:
            cur_level_locs = \
                self.tree_test_obj.get_all_locs_or_names(
                    locs, elem[4], elem[6])
            loc_min = min(cur_level_locs)
            loc_max = max(cur_level_locs)
            self.level_pts[elem[0]].append((loc_min, loc_max))

    def chk_cluster_overlap(self):
        """Checks that clusters on a single level do not overlap."""
        # Check for overlaps on each hierarchical level.
        overlap_levels = set()
        for level, ranges in sorted(self.level_pts.items()):
            ctr = cx.Counter()
            for start, stop in ranges:
                    for val in range(start, stop+1):
                        ctr[val] += 1
            if ctr.most_common()[0][1] > 1:
                overlap_levels.add(level)
        if overlap_levels:
            raise Exception('\n'.join([
                str(self),
                "CLUSTERS OVERLAP ON HIERARCHICAL LEVEL({0})".format(
                    ', '.join(map(str,overlap_levels)))]))

    def get_cluster_list(self):
        """Get a list of clusters to be used for self-checking tests."""
        # Sort hierarchy list by Hierarchy Level where 1 is the Top-Level
        cluster_list = []
        names = self.usr_input_obj.names
        for elem in self.tree_test_obj.hier_list:
            cluster_list.append(
                frozenset(self.tree_test_obj.get_all_locs_or_names(
                    names, elem[4], elem[6])))
        return cluster_list

    def __str__(self, title="HIERARCHY DIAGRAM"):
        sout = StringIO()
        sout.write('\nHUMAN INPUT:\n')
        sout.write(str(self.usr_input_obj))
        sout.write('\nINPUT TO Bio.Cluster.treecluster:\n')
        self.prt_dist_matrix(self.usr_input_obj.names, sout)
        sout.write('\nOUTPUT FROM Bio.Cluster.treecluster:\n')
        self.tree_test_obj.prt_w_names(self.usr_input_obj.names, sout)
        msg = "{0} (method={1})".format(title, self.method)
        self.tree_test_obj.prt_hier(sout, "{0}, 1st FORMAT".format(msg))
        self.prt_ascii_art_levels(sout, "{0}, 2nd FORMAT".format(msg))
        ret_str = sout.getvalue()
        sout.close()
        return ret_str

    def _mk_dist_matrix(self):
        """Given a set of locations, create a distance matrix."""
        locs = self.usr_input_obj.locs
        num_locs = len(locs)
        # Initialize distance matrix
        dist_matrix = np.zeros(((num_locs*(num_locs-1))/2,), dtype=np.float64)
        idx = 0
        for i in range(1, num_locs):
            for j in range(i):
                dist_matrix[idx] = abs(locs[j]-locs[i])
                idx += 1
        return dist_matrix

    def chk_cluster(self, exp_cluster_list):
        """Check that there exists a cluster with only the elements listed.

           For example, given the input:

             " A B  C D     E     G H I K L     M N      O P Q    R"

           It is expected that somewhere in the list of clusters generated by
           Bio.Cluster.treecluster, there will be a cluster that contains
           ONLY "G H I K L".

           If there is no cluster generated which contains "G H I K L"",
           then throw an Exception.
        """
        # Get the list of actual clusters generated by Bio.Cluster.treecluster
        actual_clusters = self.tree_test_obj.actual_clusters
        node_ids = [elem[2] for elem in self.tree_test_obj.hier_list]
        # Loop through the list of expected clusters
        for clu_str in exp_cluster_list:
            # Convert a string representing a cluster into a set
            expected = frozenset([Name for Name in clu_str if Name is not " "])
            test_pass = False
            for actual in actual_clusters:
                if actual == expected:
                    test_pass = True
                    break
            if test_pass is False:
                err_msg = ''.join([
                    self.__str__("ACTUAL CLUSTERS"),
                    "\n**FAILURE: ",
                    "EXPECTED CLUSTER({0}) NOT FOUND".format(clu_str),
                    " IN ANY ACTUAL CLUSTERS RETURNED BY ",
                    "Bio.Cluster.treecluster\n"])
                raise Exception(err_msg)
                self.tree_test_obj.prt_hier(sout, "ACTUAL CLUSTERS, 1st FORMAT")
                self.prt_ascii_art_levels(sout, "ACTUAL CLUSTERS, 2nd FORMAT")
            else:
                msg = '  Sub-test Passed.  Cluster found({0})\n'.format(clu_str)
                sys.stdout.write(msg)
        sys.stdout.write('  TEST PASSED.\n')

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))


