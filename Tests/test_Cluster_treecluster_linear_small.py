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
    """Unit tests for testing the clustering of small sets of data."""

    def test_treecluster_0(self):
        """This is a user-created test which fails for issue biopython#424."""
        usr_input = (" A B  C D     E     G H I K L     M N      O P Q    R")
        expected_clusters = [
            "A B",
            "C D",
            "A B  C D",
            # To (mock) "PASS": COMMENT OUT FOLLOWING LINE OUT
            "G H I K L",
            "M N",
            "O P Q"]

        # Run Hierarchical Tree Clustering
        test_obj = ClusterTestHelper(usr_input)
        # Check results
        test_obj.chk_cluster(expected_clusters)


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
        prt.write('    index :   {0}  \n'.format(index))
        prt.write('    names = [ {0} ]\n'.format(names))
        prt.write('    locs  = [ {0} ]\n'.format(locas))

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

    def __init__(self, dist_matrix, names):
        self.tree = Cluster.treecluster(distancematrix=dist_matrix)
        # Reformat tree for easier printing in various formats
        self.max_level = None # get_hier_list will set this value
        self.hier_list = self.get_hier_list(names)

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
        for idx, elem in enumerate(self.hier_list):
            hier_dashes = '-'*elem[0]
            hier = '{0:{1}}'.format(hier_dashes, max_lev)
            pat = '    {0:2} {1} {2:6} ({3:>4}, {4:>5}): {5:>6}\n'
            prt.write(pat.format(idx, hier, elem[2], elem[3], elem[5], elem[7]))

    def get_hier_list(self, names):
        """Get the hierarchy."""
        hier_list = []
        # Assume top-level node is the last node in the treecluster list
        self._get_hier_list(hier_list, names, -1*len(self.tree), level=1)
        self.max_level = max(elem[0] for elem in hier_list)
        return hier_list

    def _get_hier_list(self, hier_list, names, node_num, level):
        """Recursively called get_hierarchy."""
        # If this is not a leaf node
        if node_num < 0:
            idx = -1*node_num - 1
            node = self.tree[idx]
            # Print the intermediate node
            left = self.get_node_str(node.left, names)
            right = self.get_node_str(node.right, names)
            elems = [level, idx, node_num, left,
                     node.left, right, node.right, node.distance]
            hier_list.append(elems)
            self._get_hier_list(hier_list, names, node.left, level+1)
            self._get_hier_list(hier_list, names, node.right, level+1)

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

    def __init__(self, elem_str):
        self.usr_input_obj = ClustersStrInput(elem_str)
        self.dist_matrix_orig = self._mk_dist_matrix()
        warn = "**WARNING: Bio.Cluster.treecluster \
                CHANGES THE DISTANCE MATRIX, SO MAKE A COPY\n"
        sys.stdout.write(warn)
        self.dist_matrix = copy.deepcopy(self.dist_matrix_orig)
        #print self.dist_matrix_orig.dist_matrix
        #print self.dist_matrix
        self.tree_test_obj = ClusterTree(
            self.dist_matrix, self.usr_input_obj.names)
        #print self.dist_matrix
        #print self.dist_matrix_orig.dist_matrix

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
        level = cx.defaultdict(list)
        for elem in sorted(self.tree_test_obj.hier_list, key=itemgetter(0)):
            level[elem[0]].append((elem[1:]))
        num_locs = max(self.usr_input_obj.locs) + 1
        locs = self.usr_input_obj.locs
        self.usr_input_obj.prt_ascii_art(prt, title)
        for level, elem in sorted(level.items(), key=itemgetter(0)):
            ret_str = [' ']*num_locs
            for pts in elem:
                cur_level_locs = \
                    self.tree_test_obj.get_all_locs_or_names(
                        locs, pts[3], pts[5])
                loc_min = min(cur_level_locs)
                loc_max = max(cur_level_locs)
                for i in range(loc_min, loc_max+1):
                    ret_str[i] = '#'
            prt.write('      Level {0:2} {1}\n'.format(level, ''.join(ret_str)))


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

    def __str__(self):
        sout = StringIO()
        sout.write('\nHUMAN INPUT:\n')
        sout.write(str(self.usr_input_obj))
        sout.write('\nINPUT TO Bio.Cluster.treecluster:\n')
        self.prt_dist_matrix(self.usr_input_obj.names, sout)
        sout.write('\nOUTPUT FROM Bio.Cluster.treecluster:\n')
        self.tree_test_obj.prt_w_names(self.usr_input_obj.names, sout)
        self.tree_test_obj.prt_hier(sout, "HIERARCHY DIAGRAM, 1st FORMAT")
        self.prt_ascii_art_levels(sout, "HIERARCHY DIAGRAM, 2nd FORMAT")
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
        actual_clusters = self.get_cluster_list()
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
                pat = \
                    '  Bio.Cluster.treecluster NodeID({0:>3}) CLUSTER HAS: {1}'
                act_str = [
                    pat.format(num, ' '.join(clu))
                    for num, clu in zip(node_ids, actual_clusters)]

                err_msg = ''.join([
                    str(self),
                    '\nACTUAL CLUSTERS RETURNED BY ',
                    'Bio.Cluster.treecluster:\n',
                    '\n'.join(act_str),
                    "\n**FAILURE: ",
                    "EXPECTED CLUSTER({0}) NOT FOUND".format(clu_str),
                    " IN ANY ACTUAL CLUSTERS RETURNED BY ",
                    "Bio.Cluster.treecluster\n"])
                raise Exception(err_msg)
            else:
                msg = '  Sub-test Passed.  Cluster found({0})\n'.format(clu_str)
                sys.stdout.write(msg)
        sys.stdout.write('  TEST PASSED.\n')

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))


