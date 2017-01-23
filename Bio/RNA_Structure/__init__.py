# -*- coding: utf-8 -*-
# Copyright 2017 by Joanna Zbijewska, Agata Gruszczyńska, Michał Karlicki.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
# Note that:
# 1. RMSD scripts and license is available separately.
#    We added it in file: calculate_rmsd and README2.
#
# 2. RNAstructure license is available separately.
#    Please consult rna.urmc.rochester.edu .

import os
from subprocess import Popen, PIPE
import os
import sys


def rmsd_calculation(x, y):
    """ Enter two PDB format files with paths to generate rmsd value between them.
        That function uses script described below.
        Calculate Root-mean-square deviation (RMSD) between structure A and B, in XYZ
        or PDB format. The order of the atoms *must* be the same for both structures.

        citation:

         - Kabsch algorithm:
           Kabsch W., 1976, A solution for the best rotation to relate two sets of
           vectors, Acta Crystallographica, A32:922-923, doi:10.1107/S0567739476001873

         - Quaternion algorithm:
           Michael W. Walker and Lejun Shao and Richard A. Volz, 1991, Estimating 3-D
           location parameters using dual number quaternions, CVGIP: Image Understanding,
           54:358-367, doi: 10.1016/1049-9660(91)90036-o

         - Implementation:
           Calculate RMSD for two XYZ structures, GitHub,
           http://github.com/charnley/rmsd """
    current_path = os.path.dirname(os.path.realpath(__file__))
    cmd = "python "+current_path+"/calculate_rmsd.py "+x+" "+y
    p = Popen(cmd, shell=True, stdout=PIPE)
    out, err = p.communicate()
    splited_output = out.split("\n")
    return(splited_output)
