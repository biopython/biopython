"""ADPB"""
print "module prettytable is needed"

#from RNA_structure import API_NDB
#from RNA_structure import API_RNA_STRAND
import os
#import subprocess
#import RNA_API
#import RBP_score
from subprocess import Popen, PIPE


def rmsd_calculation(x, y):
    """ Enter two pdb files with paths to generate rmsd value between them """
    cmd = "python calculate_rmsd "+x+" "+y
    p = Popen(cmd, shell=True, stdout=PIPE)
    out, err = p.communicate()
    return (out,err)
