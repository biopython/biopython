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

"""ADPB"""
import os
from subprocess import Popen, PIPE


def rmsd_calculation(x, y):
    """ Enter two pdb files with paths to generate rmsd value between them """
    cmd = "python calculate_rmsd "+x+" "+y
    p = Popen(cmd, shell=True, stdout=PIPE)
    out, err = p.communicate()
    splited_output = out.split("\n")
    return(splited_output)
