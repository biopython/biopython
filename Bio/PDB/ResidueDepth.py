# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# Copyright (C) 2017, Joao Rodrigues (j.p.g.l.m.rodrigues@gmail.com)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Calculation of residue depth using command line tool MSMS.

This module uses Michel Sanner's MSMS program for the surface calculation.
See: http://mgltools.scripps.edu/packages/MSMS

Residue depth is the average distance of the atoms of a residue from
the solvent accessible surface.

Residue Depth::

    from Bio.PDB.ResidueDepth import ResidueDepth
    from Bio.PDB.PDBParser import PDBParser
    parser = PDBParser()
    structure = parser.get_structure("1a8o", "Tests/PDB/1A8O.pdb")
    model = structure[0]
    rd = ResidueDepth(model)
    print(rd['A',(' ', 152, ' ')])

Direct MSMS interface, typical use::

    from Bio.PDB.ResidueDepth import get_surface
    surface = get_surface(model)

The surface is a Numeric array with all the surface vertices.

Distance to surface::

    from Bio.PDB.ResidueDepth import min_dist
    coord = (1.113, 35.393,  9.268)
    dist = min_dist(coord, surface)

where coord is the coord of an atom within the volume bound by
the surface (ie. atom depth).

To calculate the residue depth (average atom depth of the atoms
in a residue)::

    from Bio.PDB.ResidueDepth import residue_depth
    chain = model['A']
    res152 = chain[152]
    rd = residue_depth(res152, surface)

"""


import os
import tempfile
import warnings
import subprocess

import numpy

from Bio.PDB import PDBParser
from Bio.PDB import Selection
from Bio.PDB.AbstractPropertyMap import AbstractPropertyMap
from Bio.PDB.Polypeptide import is_aa

from Bio import BiopythonWarning

# Table 1: Atom Type to radius
_atomic_radii = {
    #   atom num dist  Rexplicit Runited-atom
    1: (0.57, 1.40, 1.40),
    2: (0.66, 1.40, 1.60),
    3: (0.57, 1.40, 1.40),
    4: (0.70, 1.54, 1.70),
    5: (0.70, 1.54, 1.80),
    6: (0.70, 1.54, 2.00),
    7: (0.77, 1.74, 2.00),
    8: (0.77, 1.74, 2.00),
    9: (0.77, 1.74, 2.00),
    10: (0.67, 1.74, 1.74),
    11: (0.70, 1.74, 1.86),
    12: (1.04, 1.80, 1.85),
    13: (1.04, 1.80, 1.80),  # P, S, and LonePairs
    14: (0.70, 1.54, 1.54),  # non-protonated nitrogens
    15: (0.37, 1.20, 1.20),  # H, D  hydrogen and deuterium
    16: (0.70, 0.00, 1.50),  # obsolete entry, purpose unknown
    17: (3.50, 5.00, 5.00),  # pseudoatom - big ball
    18: (1.74, 1.97, 1.97),  # Ca calcium
    19: (1.25, 1.40, 1.40),  # Zn zinc    (traditional radius)
    20: (1.17, 1.40, 1.40),  # Cu copper  (traditional radius)
    21: (1.45, 1.30, 1.30),  # Fe heme iron
    22: (1.41, 1.49, 1.49),  # Cd cadmium
    23: (0.01, 0.01, 0.01),  # pseudoatom - tiny dot
    24: (0.37, 1.20, 0.00),  # hydrogen vanishing if united-atoms
    25: (1.16, 1.24, 1.24),  # Fe not in heme
    26: (1.36, 1.60, 1.60),  # Mg magnesium
    27: (1.17, 1.24, 1.24),  # Mn manganese
    28: (1.16, 1.25, 1.25),  # Co cobalt
    29: (1.17, 2.15, 2.15),  # Se selenium
    30: (3.00, 3.00, 3.00),  # obsolete entry, original purpose unknown
    31: (1.15, 1.15, 1.15),  # Yb ytterbium +3 ion --- wild guess only
    38: (0.95, 1.80, 1.80),  # obsolete entry, original purpose unknown
}

# Table 2: Resname/Aname to Atom Type
# MSMS uses an awk/gawk pattern matching strategy that we cannot replicate
# We will take advantage of our parser to help us in the mapping.


def _get_atom_radius(atom, rtype="united"):
    """Translate an atom object to an atomic radius defined in MSMS (PRIVATE).

    Uses information from the parent residue and the atom object to define
    the atom type.

    Returns the radius (float) according to the selected type:
     - explicit (reads hydrogens)
     - united (default)

    """
    if rtype == "explicit":
        typekey = 1
    elif rtype == "united":
        typekey = 2
    else:
        raise ValueError(
            "Radius type (%r) not understood. Must be 'explicit' or 'united'" % rtype
        )

    resname = atom.parent.resname
    het_atm = atom.parent.id[0]

    at_name = atom.name
    at_elem = atom.element

    # Hydrogens
    if at_elem == "H" or at_elem == "D":
        return _atomic_radii[15][typekey]
    # HETATMs
    elif het_atm == "W" and at_elem == "O":
        return _atomic_radii[2][typekey]
    elif het_atm != " " and at_elem == "CA":
        return _atomic_radii[18][typekey]
    elif het_atm != " " and at_elem == "CD":
        return _atomic_radii[22][typekey]
    elif resname == "ACE" and at_name == "CA":
        return _atomic_radii[9][typekey]
    # Main chain atoms
    elif at_name == "N":
        return _atomic_radii[4][typekey]
    elif at_name == "CA":
        return _atomic_radii[7][typekey]
    elif at_name == "C":
        return _atomic_radii[10][typekey]
    elif at_name == "O":
        return _atomic_radii[1][typekey]
    elif at_name == "P":
        return _atomic_radii[13][typekey]
    # CB atoms
    elif at_name == "CB" and resname == "ALA":
        return _atomic_radii[9][typekey]
    elif at_name == "CB" and resname in {"ILE", "THR", "VAL"}:
        return _atomic_radii[7][typekey]
    elif at_name == "CB":
        return _atomic_radii[8][typekey]
    # CG atoms
    elif at_name == "CG" and resname in {
        "ASN",
        "ASP",
        "ASX",
        "HIS",
        "HIP",
        "HIE",
        "HID",
        "HISN",
        "HISL",
        "LEU",
        "PHE",
        "TRP",
        "TYR",
    }:
        return _atomic_radii[10][typekey]
    elif at_name == "CG" and resname == "LEU":
        return _atomic_radii[7][typekey]
    elif at_name == "CG":
        return _atomic_radii[8][typekey]
    # General amino acids in alphabetical order
    elif resname == "GLN" and at_elem == "O":
        return _atomic_radii[3][typekey]
    elif resname == "ACE" and at_name == "CH3":
        return _atomic_radii[9][typekey]
    elif resname == "ARG" and at_name == "CD":
        return _atomic_radii[8][typekey]
    elif resname == "ARG" and at_name in {"NE", "RE"}:
        return _atomic_radii[4][typekey]
    elif resname == "ARG" and at_name == "CZ":
        return _atomic_radii[10][typekey]
    elif resname == "ARG" and at_name.startswith(("NH", "RH")):
        return _atomic_radii[5][typekey]
    elif resname == "ASN" and at_name == "OD1":
        return _atomic_radii[1][typekey]
    elif resname == "ASN" and at_name == "ND2":
        return _atomic_radii[5][typekey]
    elif resname == "ASN" and at_name.startswith("AD"):
        return _atomic_radii[3][typekey]
    elif resname == "ASP" and at_name.startswith(("OD", "ED")):
        return _atomic_radii[3][typekey]
    elif resname == "ASX" and at_name.startswith("OD1"):
        return _atomic_radii[1][typekey]
    elif resname == "ASX" and at_name == "ND2":
        return _atomic_radii[3][typekey]
    elif resname == "ASX" and at_name.startswith(("OD", "AD")):
        return _atomic_radii[3][typekey]
    elif resname in {"CYS", "CYX", "CYM"} and at_name == "SG":
        return _atomic_radii[13][typekey]
    elif resname in {"CYS", "MET"} and at_name.startswith("LP"):
        return _atomic_radii[13][typekey]
    elif resname == "CUH" and at_name == "SG":
        return _atomic_radii[12][typekey]
    elif resname == "GLU" and at_name.startswith(("OE", "EE")):
        return _atomic_radii[3][typekey]
    elif resname in {"GLU", "GLN", "GLX"} and at_name == "CD":
        return _atomic_radii[10][typekey]
    elif resname == "GLN" and at_name == "OE1":
        return _atomic_radii[1][typekey]
    elif resname == "GLN" and at_name == "NE2":
        return _atomic_radii[5][typekey]
    elif resname in {"GLN", "GLX"} and at_name.startswith("AE"):
        return _atomic_radii[3][typekey]
    # Histdines and friends
    # There are 4 kinds of HIS rings: HIS (no protons), HID (proton on Delta),
    #   HIE (proton on epsilon), and HIP (protons on both)
    # Protonated nitrogens are numbered 4, else 14
    # HIS is treated here as the same as HIE
    #
    # HISL is a deprotonated HIS (the L means liganded)
    elif resname in {"HIS", "HID", "HIE", "HIP", "HISL"} and at_name in {"CE1", "CD2"}:
        return _atomic_radii[11][typekey]
    elif resname in {"HIS", "HID", "HIE", "HISL"} and at_name == "ND1":
        return _atomic_radii[14][typekey]
    elif resname in {"HID", "HIP"} and at_name in {"ND1", "RD1"}:
        return _atomic_radii[4][typekey]
    elif resname in {"HIS", "HIE", "HIP"} and at_name in {"NE2", "RE2"}:
        return _atomic_radii[4][typekey]
    elif resname in {"HID", "HISL"} and at_name in {"NE2", "RE2"}:
        return _atomic_radii[14][typekey]
    elif resname in {"HIS", "HID", "HIP", "HISL"} and at_name.startswith(("AD", "AE")):
        return _atomic_radii[4][typekey]
    # More amino acids
    elif resname == "ILE" and at_name == "CG1":
        return _atomic_radii[8][typekey]
    elif resname == "ILE" and at_name == "CG2":
        return _atomic_radii[9][typekey]
    elif resname == "ILE" and at_name in {"CD", "CD1"}:
        return _atomic_radii[9][typekey]
    elif resname == "LEU" and at_name.startswith("CD"):
        return _atomic_radii[9][typekey]
    elif resname == "LYS" and at_name in {"CG", "CD", "CE"}:
        return _atomic_radii[8][typekey]
    elif resname == "LYS" and at_name in {"NZ", "KZ"}:
        return _atomic_radii[6][typekey]
    elif resname == "MET" and at_name == "SD":
        return _atomic_radii[13][typekey]
    elif resname == "MET" and at_name == "CE":
        return _atomic_radii[9][typekey]
    elif resname == "PHE" and at_name.startswith(("CD", "CE", "CZ")):
        return _atomic_radii[11][typekey]
    elif resname == "PRO" and at_name in {"CG", "CD"}:
        return _atomic_radii[8][typekey]
    elif resname == "CSO" and at_name in {"SE", "SEG"}:
        return _atomic_radii[9][typekey]
    elif resname == "CSO" and at_name.startswith("OD"):
        return _atomic_radii[3][typekey]
    elif resname == "SER" and at_name == "OG":
        return _atomic_radii[2][typekey]
    elif resname == "THR" and at_name == "OG1":
        return _atomic_radii[2][typekey]
    elif resname == "THR" and at_name == "CG2":
        return _atomic_radii[9][typekey]
    elif resname == "TRP" and at_name == "CD1":
        return _atomic_radii[11][typekey]
    elif resname == "TRP" and at_name in {"CD2", "CE2"}:
        return _atomic_radii[10][typekey]
    elif resname == "TRP" and at_name == "NE1":
        return _atomic_radii[4][typekey]
    elif resname == "TRP" and at_name in {"CE3", "CZ2", "CZ3", "CH2"}:
        return _atomic_radii[11][typekey]
    elif resname == "TYR" and at_name in {"CD1", "CD2", "CE1", "CE2"}:
        return _atomic_radii[11][typekey]
    elif resname == "TYR" and at_name == "CZ":
        return _atomic_radii[10][typekey]
    elif resname == "TYR" and at_name == "OH":
        return _atomic_radii[2][typekey]
    elif resname == "VAL" and at_name in {"CG1", "CG2"}:
        return _atomic_radii[9][typekey]
    elif at_name in {"CD", "CD"}:
        return _atomic_radii[8][typekey]
    # Co-factors, and other weirdos
    elif (
        resname in {"FS3", "FS4"}
        and at_name.startswith("FE")
        and at_name.endswith(("1", "2", "3", "4", "5", "6", "7"))
    ):
        return _atomic_radii[21][typekey]
    elif (
        resname in {"FS3", "FS4"}
        and at_name.startswith("S")
        and at_name.endswith(("1", "2", "3", "4", "5", "6", "7"))
    ):
        return _atomic_radii[13][typekey]
    elif resname == "FS3" and at_name == "OXO":
        return _atomic_radii[1][typekey]
    elif resname == "FEO" and at_name in {"FE1", "FE2"}:
        return _atomic_radii[21][typekey]
    elif resname == "HEM" and at_name in {"O1", "O2"}:
        return _atomic_radii[1][typekey]
    elif resname == "HEM" and at_name == "FE":
        return _atomic_radii[21][typekey]
    elif resname == "HEM" and at_name in {
        "CHA",
        "CHB",
        "CHC",
        "CHD",
        "CAB",
        "CAC",
        "CBB",
        "CBC",
    }:
        return _atomic_radii[11][typekey]
    elif resname == "HEM" and at_name in {
        "NA",
        "NB",
        "NC",
        "ND",
        "N A",
        "N B",
        "N C",
        "N D",
    }:
        return _atomic_radii[14][typekey]
    elif resname == "HEM" and at_name in {
        "C1A",
        "C1B",
        "C1C",
        "C1D",
        "C2A",
        "C2B",
        "C2C",
        "C2D",
        "C3A",
        "C3B",
        "C3C",
        "C3D",
        "C4A",
        "C4B",
        "C4C",
        "C4D",
        "CGA",
        "CGD",
    }:
        return _atomic_radii[10][typekey]
    elif resname == "HEM" and at_name in {"CMA", "CMB", "CMC", "CMD"}:
        return _atomic_radii[9][typekey]
    elif resname == "HEM" and at_name == "OH2":
        return _atomic_radii[2][typekey]
    elif resname == "AZI" and at_name in {"N1", "N2", "N3"}:
        return _atomic_radii[14][typekey]
    elif resname == "MPD" and at_name in {"C1", "C5", "C6"}:
        return _atomic_radii[9][typekey]
    elif resname == "MPD" and at_name == "C2":
        return _atomic_radii[10][typekey]
    elif resname == "MPD" and at_name == "C3":
        return _atomic_radii[8][typekey]
    elif resname == "MPD" and at_name == "C4":
        return _atomic_radii[7][typekey]
    elif resname == "MPD" and at_name in {"O7", "O8"}:
        return _atomic_radii[2][typekey]
    elif resname in {"SO4", "SUL"} and at_name == "S":
        return _atomic_radii[13][typekey]
    elif resname in {"SO4", "SUL", "PO4", "PHO"} and at_name in {
        "O1",
        "O2",
        "O3",
        "O4",
    }:
        return _atomic_radii[3][typekey]
    elif resname == "PC " and at_name in {"O1", "O2", "O3", "O4"}:
        return _atomic_radii[3][typekey]
    elif resname == "PC " and at_name == "P1":
        return _atomic_radii[13][typekey]
    elif resname == "PC " and at_name in {"C1", "C2"}:
        return _atomic_radii[8][typekey]
    elif resname == "PC " and at_name in {"C3", "C4", "C5"}:
        return _atomic_radii[9][typekey]
    elif resname == "PC " and at_name == "N1":
        return _atomic_radii[14][typekey]
    elif resname == "BIG" and at_name == "BAL":
        return _atomic_radii[17][typekey]
    elif resname in {"POI", "DOT"} and at_name in {"POI", "DOT"}:
        return _atomic_radii[23][typekey]
    elif resname == "FMN" and at_name in {"N1", "N5", "N10"}:
        return _atomic_radii[4][typekey]
    elif resname == "FMN" and at_name in {
        "C2",
        "C4",
        "C7",
        "C8",
        "C10",
        "C4A",
        "C5A",
        "C9A",
    }:
        return _atomic_radii[10][typekey]
    elif resname == "FMN" and at_name in {"O2", "O4"}:
        return _atomic_radii[1][typekey]
    elif resname == "FMN" and at_name == "N3":
        return _atomic_radii[14][typekey]
    elif resname == "FMN" and at_name in {"C6", "C9"}:
        return _atomic_radii[11][typekey]
    elif resname == "FMN" and at_name in {"C7M", "C8M"}:
        return _atomic_radii[9][typekey]
    elif resname == "FMN" and at_name.startswith(("C1", "C2", "C3", "C4", "C5")):
        return _atomic_radii[8][typekey]
    elif resname == "FMN" and at_name.startswith(("O2", "O3", "O4")):
        return _atomic_radii[2][typekey]
    elif resname == "FMN" and at_name.startswith("O5"):
        return _atomic_radii[3][typekey]
    elif resname == "FMN" and at_name in {"OP1", "OP2", "OP3"}:
        return _atomic_radii[3][typekey]
    elif resname in {"ALK", "MYR"} and at_name == "OT1":
        return _atomic_radii[3][typekey]
    elif resname in {"ALK", "MYR"} and at_name == "C01":
        return _atomic_radii[10][typekey]
    elif resname == "ALK" and at_name == "C16":
        return _atomic_radii[9][typekey]
    elif resname == "MYR" and at_name == "C14":
        return _atomic_radii[9][typekey]
    elif resname in {"ALK", "MYR"} and at_name.startswith("C"):
        return _atomic_radii[8][typekey]
    # Metals
    elif at_elem == "CU":
        return _atomic_radii[20][typekey]
    elif at_elem == "ZN":
        return _atomic_radii[19][typekey]
    elif at_elem == "MN":
        return _atomic_radii[27][typekey]
    elif at_elem == "FE":
        return _atomic_radii[25][typekey]
    elif at_elem == "MG":
        return _atomic_radii[26][typekey]
    elif at_elem == "CO":
        return _atomic_radii[28][typekey]
    elif at_elem == "SE":
        return _atomic_radii[29][typekey]
    elif at_elem == "YB":
        return _atomic_radii[31][typekey]
    # Others
    elif at_name == "SEG":
        return _atomic_radii[9][typekey]
    elif at_name == "OXT":
        return _atomic_radii[3][typekey]
    # Catch-alls
    elif at_name.startswith(("OT", "E")):
        return _atomic_radii[3][typekey]
    elif at_name.startswith("S"):
        return _atomic_radii[13][typekey]
    elif at_name.startswith("C"):
        return _atomic_radii[7][typekey]
    elif at_name.startswith("A"):
        return _atomic_radii[11][typekey]
    elif at_name.startswith("O"):
        return _atomic_radii[1][typekey]
    elif at_name.startswith(("N", "R")):
        return _atomic_radii[4][typekey]
    elif at_name.startswith("K"):
        return _atomic_radii[6][typekey]
    elif at_name in {"PA", "PB", "PC", "PD"}:
        return _atomic_radii[13][typekey]
    elif at_name.startswith("P"):
        return _atomic_radii[13][typekey]
    elif resname in {"FAD", "NAD", "AMX", "APU"} and at_name.startswith("O"):
        return _atomic_radii[1][typekey]
    elif resname in {"FAD", "NAD", "AMX", "APU"} and at_name.startswith("N"):
        return _atomic_radii[4][typekey]
    elif resname in {"FAD", "NAD", "AMX", "APU"} and at_name.startswith("C"):
        return _atomic_radii[7][typekey]
    elif resname in {"FAD", "NAD", "AMX", "APU"} and at_name.startswith("P"):
        return _atomic_radii[13][typekey]
    elif resname in {"FAD", "NAD", "AMX", "APU"} and at_name.startswith("H"):
        return _atomic_radii[15][typekey]
    else:
        warnings.warn(f"{at_name}:{resname} not in radii library.", BiopythonWarning)
        return 0.01


def _read_vertex_array(filename):
    """Read the vertex list into a Numeric array (PRIVATE)."""
    with open(filename) as fp:
        vertex_list = []
        for l in fp:
            sl = l.split()
            if len(sl) != 9:
                # skip header
                continue
            vl = [float(x) for x in sl[0:3]]
            vertex_list.append(vl)
    return numpy.array(vertex_list)


def get_surface(model, MSMS="msms"):
    """Represent molecular surface as a vertex list array.

    Return a Numpy array that represents the vertex list of the
    molecular surface.

    Arguments:
     - MSMS - msms executable (used as argument to subprocess.call)

    """
    # Replace pdb_to_xyzr
    # Make x,y,z,radius file
    atom_list = Selection.unfold_entities(model, "A")

    xyz_tmp = tempfile.mktemp()
    with open(xyz_tmp, "w") as pdb_to_xyzr:
        for atom in atom_list:
            x, y, z = atom.coord
            radius = _get_atom_radius(atom, rtype="united")
            pdb_to_xyzr.write(f"{x:6.3f}\t{y:6.3f}\t{z:6.3f}\t{radius:1.2f}\n")

    # make surface
    surface_tmp = tempfile.mktemp()
    MSMS = MSMS + " -probe_radius 1.5 -if %s -of %s > " + tempfile.mktemp()
    make_surface = MSMS % (xyz_tmp, surface_tmp)
    subprocess.call(make_surface, shell=True)
    surface_file = surface_tmp + ".vert"
    if not os.path.isfile(surface_file):
        raise RuntimeError(
            "Failed to generate surface file using command:\n%s" % make_surface
        )

    # read surface vertices from vertex file
    surface = _read_vertex_array(surface_file)
    return surface


def min_dist(coord, surface):
    """Return minimum distance between coord and surface."""
    d = surface - coord
    d2 = numpy.sum(d * d, 1)
    return numpy.sqrt(min(d2))


def residue_depth(residue, surface):
    """Residue depth as average depth of all its atoms.

    Return average distance to surface for all atoms in a residue,
    ie. the residue depth.
    """
    atom_list = residue.get_unpacked_list()
    length = len(atom_list)
    d = 0
    for atom in atom_list:
        coord = atom.get_coord()
        d = d + min_dist(coord, surface)
    return d / length


def ca_depth(residue, surface):
    """Return CA depth."""
    if not residue.has_id("CA"):
        return None
    ca = residue["CA"]
    coord = ca.get_coord()
    return min_dist(coord, surface)


class ResidueDepth(AbstractPropertyMap):
    """Calculate residue and CA depth for all residues."""

    def __init__(self, model, msms_exec=None):
        """Initialize the class."""
        if msms_exec is None:
            msms_exec = "msms"

        depth_dict = {}
        depth_list = []
        depth_keys = []
        # get_residue
        residue_list = Selection.unfold_entities(model, "R")
        # make surface from PDB file using MSMS
        surface = get_surface(model, MSMS=msms_exec)
        # calculate rdepth for each residue
        for residue in residue_list:
            if not is_aa(residue):
                continue
            rd = residue_depth(residue, surface)
            ca_rd = ca_depth(residue, surface)
            # Get the key
            res_id = residue.get_id()
            chain_id = residue.get_parent().get_id()
            depth_dict[(chain_id, res_id)] = (rd, ca_rd)
            depth_list.append((residue, (rd, ca_rd)))
            depth_keys.append((chain_id, res_id))
            # Update xtra information
            residue.xtra["EXP_RD"] = rd
            residue.xtra["EXP_RD_CA"] = ca_rd
        AbstractPropertyMap.__init__(self, depth_dict, depth_keys, depth_list)
