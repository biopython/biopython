# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Classes that deal with macromolecular crystal structures.

Includes: PDB and mmCIF parsers, a Structure class, a module to keep a local
copy of the PDB up-to-date, selective IO of PDB files, etc.

Author: Thomas Hamelryck.  Additional code by Kristian Rother.
"""

# Get a Structure object from a PDB file
from .PDBParser import PDBParser

from .MMCIFParser import MMCIFParser
from .MMCIFParser import FastMMCIFParser

# Download from the PDB
from .PDBList import PDBList

# Parse PDB header directly
from .parse_pdb_header import parse_pdb_header

# Find connected polypeptides in a Structure
from .Polypeptide import PPBuilder, CaPPBuilder, is_aa, standard_aa_names

# This is also useful :-)
from Bio.Data.SCOPData import protein_letters_3to1

# IO of PDB files (including flexible selective output)
from .PDBIO import PDBIO, Select
from .mmcifio import MMCIFIO

# Some methods to eg. get a list of Residues
# from a list of Atoms.
from . import Selection

# Superimpose atom sets
from .Superimposer import Superimposer

# 3D vector class
from .vectors import Vector, calc_angle, calc_dihedral, refmat, rotmat, rotaxis
from .vectors import vector_to_axis, m2rotaxis, rotaxis2m

# Alignment module
from .StructureAlignment import StructureAlignment

# DSSP handle
# (secondary structure and solvent accessible area calculation)
from .DSSP import DSSP, make_dssp_dict

# Residue depth:
# distance of residue atoms from solvent accessible surface
from .ResidueDepth import ResidueDepth, get_surface

# Calculation of Half Sphere Solvent Exposure
from .HSExposure import HSExposureCA, HSExposureCB, ExposureCN

# Kolodny et al.'s backbone libraries
from .FragmentMapper import FragmentMapper

# Write out chain(start-end) to PDB file
from .Dice import extract

# Fast atom neighbor search
# Depends on kdtrees C module
try:
    from .NeighborSearch import NeighborSearch
except ImportError:
    pass
