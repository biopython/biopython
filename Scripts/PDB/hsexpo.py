#!/usr/bin/python

# Copyright (C) 2002, Thomas Hamelryck (thamelry@vub.ac.be)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Calculates solvent exposure for a PDB file using one of 5 different methods.

    -DSSP (DSSP needs to be installed)
    -Residue depth (MSMS needs to be installed)
    -Coordination number (ie. number of CA atoms within a sphere)
    -HSEalpha half sphere exposure
    -HSEbeta half sphere exposure

A PDB file can be written out with the exposure in the B factor field.
See --help for all args.
"""


import argparse
import sys

from Bio.PDB import (
    DSSP,
    ExposureCN,
    HSExposureCA,
    HSExposureCB,
    PDBParser,
    PDBIO,
    ResidueDepth,
    Selection,
)

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("pdbfile", help="Input structure in PDB format.")
ap.add_argument(
    "-t",
    "--type",
    dest="exp",
    choices=["HSEAU", "HSEAD", "HSEBU", "HSEBD", "CN", "DSSPr", "DSSPa", "RD", "RDa"],
    help="Exposure Type",
    default="HSEb",
)
ap.add_argument(
    "-o", "--out", dest="outfile", help="output to PDB file (B factor=exposure)"
)
ap.add_argument(
    "-r",
    "--radius",
    dest="radius",
    type=float,
    help="sphere radius (default 13.0 A)",
    default=13.0,
)
ap.add_argument(
    "-m",
    "--model",
    dest="model",
    type=int,
    help="PDB model number (default 0)",
    default=0,
)
ap.add_argument("--dssp", help="Path to the DSSP executable", default="dssp")
ap.add_argument("--msms", help="Path to the MSMS executable", default=None)
args = ap.parse_args()

# Get the structure
p = PDBParser()
s = p.get_structure("X", args.pdbfile)

# First model by default
m = s[args.model]

RADIUS = args.radius

# d=dictionary of exposures
# k=position in ntuple containing the desired exposure

format = "%4i"

args.exp = args.exp.upper()

if args.exp[0] == "H" and args.exp[3] == "A":
    hse = HSExposureCA(m, RADIUS)
    if args.exp[-1] == "D":
        k = "EXP_HSE_A_D"
    else:
        k = "EXP_HSE_A_U"
elif args.exp[0] == "H" and args.exp[3] == "B":
    hse = HSExposureCB(m, RADIUS)
    # hse.write_pymol_script()
    if args.exp[-1] == "D":
        k = "EXP_HSE_B_U"
    else:
        k = "EXP_HSE_B_D"
elif args.exp == "CN":
    hse = ExposureCN(m, RADIUS)
    k = "EXP_CN"
elif args.exp == "ANGLE":
    hse = HSExposureCA(m, RADIUS)
    k = "EXP_CB_PCB_ANGLE"
    format = "%4.1f"
elif args.exp == "DSSPR":
    d = DSSP(m, args.pdbfile, dssp=args.dssp)
    k = "EXP_DSSP_RASA"
    format = "%.4f"
elif args.exp == "DSSPA":
    d = DSSP(m, args.pdbfile, dssp=args.dssp)
    k = "EXP_DSSP_ASA"
elif args.exp == "RD":
    d = ResidueDepth(m, args.pdbfile, msms_exec=args.msms)
    k = "EXP_RD"
    format = "%4.1f"
elif args.exp == "RDA":
    d = ResidueDepth(m, args.pdbfile, msms_exec=args.msms)
    k = "EXP_RD_CA"
    format = "%4.1f"
else:
    print("ERROR: Unknown option.")
    sys.exit()

residue_list = Selection.unfold_entities(m, "R")

for r in residue_list:

    if k in r.xtra:

        exposure = r.xtra[k]

        if args.exp == "DSSPR":
            # to 0=exposed, 1=buried
            exposure = 1 - exposure

        # Print info
        hetflag, resseq, icode = r.get_id()

        if icode == " ":
            icode = "_"

        resname = r.get_resname()

        print(("%s %4i %c\t" + format) % (resname, resseq, icode, exposure))
    else:
        exposure = 0.0

    for atom in r.get_iterator():
        atom.set_bfactor(exposure)

if args.outfile:
    io = PDBIO()
    io.set_structure(s)
    io.save(args.outfile)
