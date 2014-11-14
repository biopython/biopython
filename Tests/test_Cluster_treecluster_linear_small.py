#!/usr/bin/env python

# Copyright 2015 by DV Klopfenstein.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This test and test infrastructure was created to help verify this issue:
# 
#   https://github.com/biopython/biopython/issues/424#issuecomment-62823368

import sys
import cStringIO
import copy
import numpy as np
import collections as cx
from operator import itemgetter

from Bio import Cluster


#############################################################################
# User-Created Tests: Edit this area to add more tests
#############################################################################
def run_all_tests():
  test_treecluster_0()



def test_treecluster_0():
  """This is a user-created test which fails for issue biopython#424."""
  INPUT = (" A B  C D     E     G H I K L     M N      O P Q    R")
  EXP_CLUSTERS = [
    "A B",
    "C D",
    "A B C D",
    "G H I K L", # To PASS: COMMENT THIS LINE OUT to make test PASS for issue biopython#424
    "M N",
    "O P Q"]

  A = ClusterTestHelper( INPUT ) # Run Hierarchical Tree Clustering
  A.chk_cluster( EXP_CLUSTERS )  # Check results




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
    nl = [ (Name, Loc) for Loc, Name in enumerate(self.clu_str) if Name is not " "]
    return zip(*nl)

  def prt_arrays(self, PRT=sys.stdout):
    L = len(self.names)
    PRT.write('\n  VARIABLES:\n')
    PRT.write('    index :   {}  \n'.format('  '.join('{:>2}'.format(E) for E in range(L))))
    PRT.write('    names = [ {} ]\n'.format(', '.join('{:>2}'.format(E) for E in self.names)))
    PRT.write('    locs  = [ {} ]\n'.format(', '.join('{:>2}'.format(E) for E in self.locs)))

  def prt_ascii_art(self, PRT=sys.stdout, title="ASCII Art"):
    L = len(self.clu_str)
    PRT.write('\n  {}:\n'.format(title))
    PRT.write('    Index 10s: {}\n'.format(''.join('{}'.format(I/10) if I%10==0 else ' ' for I in range(L))))
    PRT.write('    Index  1s: {}\n'.format(''.join('{}'.format(I%10) for I in range(L))))
    PRT.write('    Array    : {}\n'.format(self.clu_str))

  def __str__(self):
    SOUT = cStringIO.StringIO()
    self.prt_ascii_art(SOUT)
    self.prt_arrays(SOUT)
    Str = SOUT.getvalue()
    SOUT.close()
    return Str


class DistMatrixHelper:
  """Class used to create, print and store distance matrixes."""

  def __init__(self, locs):
    self.dist_matrix = self._mk_dist_matrix(locs)

  def _mk_dist_matrix(self, D):
    """Given a set of locations, create a distance matrix."""
    L = len(D)
    # Initialize distance matrix
    dist_matrix = np.zeros( ((L*(L-1))/2,), dtype=np.float64 )
    idx = 0
    for i in range(1,L):
      for j in range(i):
        dist_matrix[idx] = abs(D[j]-D[i])
        idx += 1
    return dist_matrix


class ClusterTree:
  """Used to create, store and print Biopython Cluster treecluster results."""

  def __init__(self, dist_matrix, names):
    self.tree = Cluster.treecluster(distancematrix = dist_matrix)
    # Reformat tree for easier printing in various formats
    self.max_level = None # get_hier_list will set this value
    self.hier_list = self.get_hier_list(names)

  def prt(self, PRT=sys.stdout):
    prt_w_names(None, PRT)

  def prt_w_names(self, names, PRT=sys.stdout):
    PRT.write('\n  TREE ARRAY FROM Bio.Cluster.treecluster:\n')
    PRT.write('               NodeID (left, right): distance\n')
    for I, E in enumerate(self.tree):
      self.prt_node(I, E, names, PRT)

  def str_node(self, I, node, names):
    SOUT = cStringIO.StringIO()
    self.prt_node(I, node, names, SOUT)
    Str = SOUT.getvalue()
    SOUT.close()
    return Str

  def prt_node(self, I, node, names, PRT):
    L = self.get_node_str(node.left,  names)
    R = self.get_node_str(node.right, names)
    PRT.write('    tree[{:2}] = {:6} ({:>4}, {:>5}): {:>6}\n'.format(I, -I-1, L, R, node.distance))

  def get_node_str(self, LR, names):
    return LR if LR  < 0 or names is None else names[LR]

  def prt_hier(self, names, PRT=sys.stdout, title="HIERARCHY DIAGRAM"):
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
    PRT.write('\n  {}:\n'.format(title))
    for I, E in enumerate(self.hier_list):
      PRT.write('    {:2} {:{}} {:6} ({:>4}, {:>5}): {:>6}\n'.format(
        I, '-'*E[0], self.max_level, E[2], E[3], E[5], E[7]))

  def get_hier_list(self, names):
    """Get the hierarchy."""
    hier_list = []
    # Assuming top-level node is the last node in the tree list from treecluster...
    self._get_hier_list(hier_list, names, -1*len(self.tree), level=1)
    self.max_level = max(E[0] for E in hier_list)
    return hier_list

  def _get_hier_list(self, hier_list, names, node_num, level):
    """Recursively called get_hierarchy."""
    # If this is not a leaf node
    if node_num < 0: 
      I = -1*node_num - 1
      node = self.tree[I]
      # Print the intermediate node
      L = self.get_node_str(node.left,  names)
      R = self.get_node_str(node.right, names)
      hier_list.append([level, I, node_num, L, node.left, R, node.right, node.distance])
      self._get_hier_list(hier_list, names, node.left,  level+1)
      self._get_hier_list(hier_list, names, node.right, level+1)

  def __str__(self):
    SOUT = cStringIO.StringIO()
    self.prt(SOUT)
    self.prt_hier(SOUT)
    Str = SOUT.getvalue()
    SOUT.close()
    return Str

  def get_all_locs_or_names(self, all_locs, L, R):
    cur_locs = []
    self._get_all_locs_or_names(cur_locs, all_locs, L)
    self._get_all_locs_or_names(cur_locs, all_locs, R)
    return cur_locs

  def _get_all_locs_or_names(self, cur_locs, all_locs, nodeid_locidx):
    # If this is a leaf node, append cur_loc list
    if nodeid_locidx >= 0: 
      cur_locs.append( all_locs[nodeid_locidx] )
      return
    # If this is not intermediate node, keep drilling
    node = self.tree[-1*nodeid_locidx - 1]
    self._get_all_locs_or_names(cur_locs, all_locs, node.left)
    self._get_all_locs_or_names(cur_locs, all_locs, node.right)


class ClusterTestHelper:
  """Class used to take items in a string and generate clusters."""

  def __init__(self, elem_str):
    self.IN = ClustersStrInput(elem_str)
    self.DM = DistMatrixHelper(self.IN.locs)
    sys.stdout.write("**WARNING: Bio.Cluster.treecluster CHANGES THE DISTANCE MATRIX, SO MAKE A COPY\n")
    self.dist_matrix = copy.deepcopy(self.DM.dist_matrix)
    #print self.DM.dist_matrix
    #print self.dist_matrix
    self.TR = ClusterTree(self.dist_matrix, self.IN.names)
    #print self.dist_matrix
    #print self.DM.dist_matrix

  def prt_dist_matrix(self, names, PRT=sys.stdout):
    cnt = 1
    cnt_stop  = 1
    PRT.write('\n  dist_matrix = [\n')
    PRT.write('   # {}\n  '.format(names[0]))
    for R in self.DM.dist_matrix:
      PRT.write('  {:4.1f}, '.format(R))
      if cnt == cnt_stop:
        PRT.write(' # {}\n  '.format(names[cnt_stop]))
        cnt = 0
        cnt_stop += 1
      cnt += 1
    PRT.write('  ]\n')

  def prt_ascii_art_levels(self, PRT=sys.stdout, title="HIERARCHY DIAGRAM"):
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

       THE ABOVE TREE FORMAT IS EASY TO CREATE IN A SMALL AMOUNT OF CODE
       AND IS EQUIVALENT TO THIS MORE FAMILIAR FORMAT:


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
    Lev = cx.defaultdict(list)
    for E in sorted(self.TR.hier_list, key=itemgetter(0)):
      Lev[E[0]].append((E[1:]))
    L = max(self.IN.locs) + 1
    locs = self.IN.locs
    self.IN.prt_ascii_art(PRT, title)
    for level, E in sorted(Lev.items(), key=itemgetter(0)):
      Str = [' ']*L
      for P in E:
        cur_level_locs = self.TR.get_all_locs_or_names(locs, P[3], P[5])
        loc_min = min(cur_level_locs)
        loc_max = max(cur_level_locs)
        for i in range(loc_min, loc_max+1):
          Str[i] = '#'
      PRT.write('      Level {:2} {}\n'.format(level, ''.join(Str)))
      

  def get_cluster_list(self):
    """Get a list of clusters to be used for self-checking tests."""
    # Sort hierarchy list by Hierarchy Level where 1 is the Top-Level
    cluster_list = []
    names = self.IN.names
    for E in self.TR.hier_list:
      cluster_list.append(frozenset(self.TR.get_all_locs_or_names(names, E[4], E[6])))
    return cluster_list
      
  def __str__(self):
    SOUT = cStringIO.StringIO()
    SOUT.write('\nHUMAN INPUT:\n')
    SOUT.write(str(self.IN))
    SOUT.write('\nINPUT TO Bio.Cluster.treecluster:\n')
    self.prt_dist_matrix(self.IN.names, SOUT)
    SOUT.write('\nOUTPUT FROM Bio.Cluster.treecluster:\n')
    self.TR.prt_w_names(self.IN.names, SOUT)
    self.TR.prt_hier(self.IN.names, SOUT, "HIERARCHY DIAGRAM, 1st FORMAT")
    self.prt_ascii_art_levels(SOUT,       "HIERARCHY DIAGRAM, 2nd FORMAT")
    Str = SOUT.getvalue()
    SOUT.close()
    return Str

  def chk_cluster(self, exp_cluster_list):
    """Check that there exists a cluster with only the specific elements listed.

       For example, given the input:

         " A B  C D     E     G H I K L     M N      O P Q    R"

       It is expected that somewhere in the list of clusters generated by Bio.Cluster.treecluster,
       there will be a cluster that contains ONLY "G H I K L".

       If there is no cluster at all generated which contains "G H I K L"", 
       then throw an Exception.
    """
    # Get the list of actual clusters generated by Bio.Cluster.treecluster
    act_clus = self.get_cluster_list()
    # Loop through the list of expected clusters
    for clu_str in exp_cluster_list:
      # Convert a string representing a cluster into a set 
      EXP = frozenset([ Name for Name in clu_str if Name is not " "])
      test_pass = False
      for ACT in act_clus:
        if ACT == EXP:
          test_pass = True
          break
      if test_pass is False:
        sys.stdout.write(str(self))
        sys.stdout.write('\n  LIST OF CLUSTERS RETURNED BY Bio.Cluster.treecluster:\n')
        for C in act_clus:
          sys.stdout.write('    ACTUAL CLUSTER: {}\n'.format(' '.join(C)))
        sys.stdout.write('\n')
        raise Exception("CLUSTER({}) NOT FOUND IN RETURNED VALUE FROM Bio.Cluster.treecluster".format(clu_str))
      else:
        sys.stdout.write('  Sub-test Passed.  Cluster found({})\n'.format(clu_str))
    sys.stdout.write('  TEST PASSED.\n'.format(clu_str))
    
if __name__ == '__main__':
  run_all_tests()


