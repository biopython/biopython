#!/usr/bin/python

# Copyright (C) 2002, Thomas Hamelryck (thamelry@vub.ac.be)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from optparse import OptionParser

from Bio.PDB import *
import sys

__doc__="""
This program calculates solvent exposure for all amino
acids in a PDB file using 5 different methods:

    -DSSP (DSSP needs to be installed)
    -Residue depth (MSMS needs to be installed)
    -Coordination number (ie. number of CA atoms within a sphere)
    -HSEalpha half sphere exposure
    -HSEbeta half sphere exposure

A PDB file can be written out with the exposure in the B factor field.
See --help for all options.
"""

if len(sys.argv)==1:
    print __doc__
    sys.exit()

# Get the user's options
parser=OptionParser(usage="usage: %prog [options] <PDB file>")

parser.add_option("-t", "--type", dest="exp", 
                  help="exposure type (HSEAU, HSEAD, HSEBU, HSEBD, CN, DSSPr, DSSPa, RD, RDa)",
                  default="HSEb")

parser.add_option("-o", "--out", dest="outfile", 
                  help="output to PDB file (B factor=exposure)")

parser.add_option("-r", "--radius", dest="radius", type="float", 
                  help="sphere radius (default 13.0 A)", 
                  default=13.0)

parser.add_option("-m", "--model", dest="model", type="int", 
                  help="model number (default 0)", 
                  default=0)

(options, args)=parser.parse_args()

pdbfile=args[0]

# Get the structure
p=PDBParser()
s=p.get_structure('X', pdbfile)

# First model by default
m=s[options.model]

RADIUS=options.radius

# d=dictionary of exposures
# k=position in ntuple containing the desired exposure

format="%4i"

options.exp=options.exp.upper()

if options.exp[0]=="H" and options.exp[3]=="A":
    hse=HSExposureCA(m, RADIUS)
    if options.exp[-1]=="D":
        k='EXP_HSE_A_D'
    else:
        k='EXP_HSE_A_U'
elif options.exp[0]=="H" and options.exp[3]=="B":
    hse=HSExposureCB(m, RADIUS)
    #hse.write_pymol_script()
    if options.exp[-1]=="D":
        k='EXP_HSE_B_U'
    else:
        k='EXP_HSE_B_D'
elif options.exp=="CN":
    hse=ExposureCN(m, RADIUS)
    k='EXP_CN'
elif options.exp=="ANGLE":
    hse=HSExposureCA(m, RADIUS)
    k='EXP_CB_PCB_ANGLE'
    format="%4.1f"
elif options.exp=="DSSPR":
    d=DSSP(m, pdbfile)
    k='EXP_DSSP_RASA'
    format="%.4f"
elif options.exp=="DSSPA":
    d=DSSP(m, pdbfile)
    k='EXP_DSSP_ASA'
elif options.exp=="RD":
    d=ResidueDepth(m, pdbfile)
    k='EXP_RD'
    format="%4.1f"
elif options.exp=="RDA":
    d=ResidueDepth(m, pdbfile)
    k='EXP_RD_CA'
    format="%4.1f"
else:
    print "ERROR: Unknown option."
    sys.exit()

residue_list=Selection.unfold_entities(m, 'R')

for r in residue_list:

    if r.xtra.has_key(k):

        exposure=r.xtra[k]

        if options.exp=="DSSPR":
            # to 0=exposed, 1=buried
            exposure=1-exposure

        # Print info
        hetflag, resseq, icode=r.get_id()

        if icode==' ':
            icode='_'

        resname=r.get_resname()

        print (("%s %4i %c\t"+format) % (resname, resseq, icode, exposure)) 
    else:
        exposure=0.0

    for atom in r.get_iterator():
        atom.set_bfactor(exposure)

if options.outfile:
    io=PDBIO()
    io.set_structure(s)
    io.save(options.outfile)

