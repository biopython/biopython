# Copyright (C) 2002,2020 Thomas Hamelryck (thamelry@binf.ku.dk)
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
- Simon Duerr
- Joe Greener
- Rob Miller
- Lenna X. Peterson
- Joao Rodrigues
- Kristian Rother
- Eric Talevich
- and many others.
"""

import gzip
import contextlib

from pathlib import Path

# Get a Structure object from a PDB file using available parsers
from .PDBParser import PDBParser
from .MMCIFParser import MMCIFParser
from .MMCIFParser import FastMMCIFParser

try:
    from .mmtf import MMTFParser
except ImportError:
    pass

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

# Native Shrake-Rupley algorithm for SASA calculations.
# Depends on kdtrees C module
try:
    from .SASA import ShrakeRupley
except ImportError:
    pass

_FormatToParser = {"pdb": PDBParser, "cif": MMCIFParser}


def _check_format(fmt, parseropt):
    """Return a parser that exists."""
    # TODO should a more helpful message be printed if parseropt
    # does not exist for parser and TypeError is raised?
    if fmt in _FormatToParser:
        return _FormatToParser[fmt](**parseropt)
    else:
        raise ValueError("Unknown format '%s'" % fmt)


def _as_handle_gzip(handle):
    """Open str, bytes or Path like object - gzipped or not."""
    try:
        with open(handle, "rb") as fp:
            magic = fp.read(2)
            fp.seek(0)
            if magic == b"\x1f\x8b":
                return gzip.open(handle, mode="rt")
            else:
                return open(handle)
    except TypeError:
        # should already be handle
        return handle


def read(handle, fmt, **kwargs):
    """Parse a file or handle using one of the available parsers to generate a structure.

    Arguments:
     - handle    - handle to the file, or the filename as a string
     - fmt       - string describing the file format. can be [pdb, cif, mmtf]

    If you have the file name in a string 'filename', use:

    >>> from Bio import PDB
    >>> filename = "../Tests/PDB/1A8O.pdb"
    >>> structure = PDB.read(filename, fmt="pdb")

    The function also supports file handles and gzip files

    >>> from Bio import PDB
    >>> structure = PDB.read("../Tests/PDB/1A8O.pdb.gz", fmt="pdb")

    The read function supports all arguments of the individual parsers.

    >>> from Bio import PDB
    >>> structure = PDB.read("../Tests/PDB/1A8O.pdb", fmt="pdb", QUIET=True, PERMISSIVE=True)
    """
    # set id to None internally for compatibility with parsers
    id = None
    # Try and give helpful error messages:
    if not isinstance(fmt, str):
        raise TypeError("Need a string for the file format")
    if not fmt:
        raise ValueError("Format required")

    # need to use different parser for mmtf versus pdb,cif
    if fmt.lower() == "mmtf":
        # MMTF needs a seperate parser as get_structure does only take 1 positional argument and does not support file handles
        if isinstance(handle, str):
            parser = MMTFParser()
            return parser.get_structure(handle)
        else:
            raise TypeError("'%s' parser only accepts str not handle" % fmt.lower())
    else:
        with _as_handle_gzip(handle) as fp:
            parser = _check_format(fmt.lower(), kwargs)
            try:
                structure = parser.get_structure(id, fp)
            except ValueError:
                raise ValueError(
                    "Could not parse structure using %s. Did you choose the correct parser (fmt) ?"
                    % fmt.lower()
                )
            # set id from header if given
            if structure.header["idcode"] != "":
                structure.id = structure.header["idcode"]
            return structure
