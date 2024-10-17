# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Classes that deal with macromolecular crystal structures.

Includes: PDB and mmCIF parsers, a Structure class, a module to keep a local
copy of the PDB up-to-date, selective IO of PDB files, etc.

Original Author: Thomas Hamelryck.
Contributions by:
- Peter Cock
- Joe Greener
- Rob Miller
- Lenna X. Peterson
- Joao Rodrigues
- Kristian Rother
- Eric Talevich
- and many others.
"""

try:
    import numpy as np
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Please install NumPy if you want to use Bio.PDB. See http://www.numpy.org/"
    ) from None

# Get a Structure object from a PDB file
# Some methods to eg. get a list of Residues
# from a list of Atoms.
from . import Selection

# CEAlign structural alignment
from .cealign import CEAligner

# Write out chain(start-end) to PDB file
from .Dice import extract

# DSSP handle
# (secondary structure and solvent accessible area calculation)
from .DSSP import DSSP
from .DSSP import make_dssp_dict

# Kolodny et al.'s backbone libraries
from .FragmentMapper import FragmentMapper
from .HSExposure import ExposureCN

# Calculation of Half Sphere Solvent Exposure
from .HSExposure import HSExposureCA
from .HSExposure import HSExposureCB
from .mmcifio import MMCIFIO
from .MMCIFParser import FastMMCIFParser
from .MMCIFParser import MMCIFParser

# Parse PDB header directly
from .parse_pdb_header import parse_pdb_header

# IO of PDB files (including flexible selective output)
from .PDBIO import PDBIO
from .PDBIO import Select

# Download from the PDB
from .PDBList import PDBList
from .PDBMLParser import PDBMLParser
from .PDBParser import PDBParser
from .Polypeptide import CaPPBuilder
from .Polypeptide import is_aa
from .Polypeptide import is_nucleic

# Find connected polypeptides in a Structure
from .Polypeptide import PPBuilder
from .Polypeptide import standard_aa_names
from .ResidueDepth import get_surface

# Residue depth:
# distance of residue atoms from solvent accessible surface
from .ResidueDepth import ResidueDepth

# Alignment module
from .StructureAlignment import StructureAlignment

# Superimpose atom sets
from .Superimposer import Superimposer
from .vectors import calc_angle
from .vectors import calc_dihedral
from .vectors import m2rotaxis
from .vectors import refmat
from .vectors import rotaxis
from .vectors import rotaxis2m
from .vectors import rotmat

# 3D vector class
from .vectors import Vector
from .vectors import vector_to_axis

# Fast atom neighbor search
# Depends on kdtrees C module
try:
    from .NeighborSearch import NeighborSearch
except ImportError:
    pass

# Native Shrake-Rupley algorithm for SASA calculations.
# Depends on kdtrees C module
try:
    from .SASA import ShrakeRupley
except ImportError:
    pass
