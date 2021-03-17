# Copyright (C) 2002,2020 Thomas Hamelryck (thamelry@binf.ku.dk), Simon Duerr
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
- Simon Duerr
- and many others.
"""

import gzip
from pathlib import Path
from mimetypes import guess_type
from functools import partial


# Get a Structure object from a PDB file using available parsers
from .PDBParser import PDBParser
from .MMCIFParser import MMCIFParser
from .MMCIFParser import FastMMCIFParser
from .mmtf import MMTFParser

from Bio.File import as_handle

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

_FormatToParser = {"pdb": PDBParser(), "cif": MMCIFParser()}


def read(handle, format, id=None):
    """Parse a file or handle using one of the available parsers to generate a structure.

    Arguments:
     - handle    - handle to the file, or the filename as a string
     - format    - string describing the file format. can be [pdb, cif, mmtf]
     - id        - optional id of the structure, if None header information or then base name is used

    If you have the file name in a string 'filename', use:

    >>> from Bio import PDB
    >>> filename = "1A8O.pdb"
    >>> PDB.read(filename, format="pdb")

    The function also supports file handles and gzip files

    >>> from Bio import PDB
    >>> PDB.read("1A8O.pdb.gz", format="pdb")
    """
    # TODO *kwargs for the parser to allow QUIET or permissive mode

    # Try and give helpful error messages:
    if not isinstance(format, str):
        raise TypeError("Need a string for the file format (lower case)")
    if not format:
        raise ValueError("Format required (lower case string)")
    if format != format.lower():
        raise ValueError("Format string '%s' should be lower case" % format)

    # need to use different parser for mmtf versus pdb,cif
    if format == "mmtf":
        # TODO: Raise warning if ID is set. MMTFParser does not support passing an id
        # MMTF needs a seperate parser as get_structure does only take 1 positional argument and does not support file handles
        if isinstance(handle, str):
            parser = MMTFParser()
            return parser.get_structure(handle)
        else:
            raise TypeError("'%s' parser only accepts str not handle" % format)

    elif format in _FormatToParser:  # Map the file format to a parser
        # if we are provided a path to a file, check if gzipped,
        # else open handle directly
        if isinstance(handle, str):
            encoding = guess_type(handle)[1]  # uses file extension
            _open = partial(gzip.open, mode="rt") if encoding == "gzip" else as_handle
        else:
            _open = as_handle

        with _open(handle) as fp:
            parser = _FormatToParser[format]
            structure = parser.get_structure(id, fp)
            # set id from header if given, else use basename
            if structure.header["idcode"] != "" and id is None:
                structure.id = structure.header["idcode"]
            elif isinstance(handle, str) and id is None:
                structure.id = Path(handle).stem
            # TODO set id if filehandle without header is passed
            return structure
    else:
        raise ValueError("Unknown format '%s'" % format)
