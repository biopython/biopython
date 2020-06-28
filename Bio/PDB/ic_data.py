# Copyright 2019 by Robert T. Miller.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Per residue backbone and sidechain hedra and dihedra definitions.

Listed in order of output for internal coordinates (.pic) output file.
Require sufficient overlap to link all defined dihedra.  Entries in these
tables without corresponding atom coordinates are ignored.

<http://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/formuleAA/>
for naming of individual atoms
"""

# Backbone hedra and dihedra - within residue, no next or prev.
ic_data_backbone = (
    ("N", "CA", "C", "O"),  # locate backbone O
    ("O", "C", "CA", "CB"),  # locate CB
    ("CA", "C", "O"),
    ("CB", "CA", "C"),
    ("CA", "C", "OXT"),  # OXT if present
    ("N", "CA", "C", "OXT"),
    ("H", "N", "CA"),  # amide proton if present
    ("C", "CA", "N", "H"),
    ("HA", "CA", "C"),  # CA proton
    ("O", "C", "CA", "HA"),
    ("HA2", "CA", "C"),  # gly CA proton
    ("O", "C", "CA", "HA2"),
    ("HA3", "CA", "C"),  # gly CA proton
    ("O", "C", "CA", "HA3"),
    ("N", "CA", "CB"),
    ("N", "CA", "CB", "HB"),  # CB protons
    ("N", "CA", "CB", "HB1"),
    ("N", "CA", "CB", "HB2"),
    ("N", "CA", "CB", "HB3"),
    ("CA", "CB", "HB"),
    ("CA", "CB", "HB1"),
    ("CA", "CB", "HB2"),
    ("CA", "CB", "HB3"),
    ("H1", "N", "CA"),  # chain start protons
    ("H2", "N", "CA"),
    ("H3", "N", "CA"),
    ("C", "CA", "N", "H1"),
    ("C", "CA", "N", "H2"),
    ("C", "CA", "N", "H3"),
)

# Sidechain hedra and dihedra.
ic_data_sidechains = {
    "V": (
        ("CA", "CB", "CG1"),
        ("N", "CA", "CB", "CG1", "chi1"),  # chi1
        ("CA", "CB", "CG2"),
        ("N", "CA", "CB", "CG2"),
        ("CB", "CG1", "HG11"),
        ("CB", "CG1", "HG12"),
        ("CB", "CG1", "HG13"),
        ("CB", "CG2", "HG21"),
        ("CB", "CG2", "HG22"),
        ("CB", "CG2", "HG23"),
        ("CA", "CB", "CG1", "HG11"),
        ("CA", "CB", "CG1", "HG12"),
        ("CA", "CB", "CG1", "HG13"),
        ("CA", "CB", "CG2", "HG21"),
        ("CA", "CB", "CG2", "HG22"),
        ("CA", "CB", "CG2", "HG23"),
    ),
    "L": (
        ("CA", "CB", "CG"),
        ("N", "CA", "CB", "CG", "chi1"),  # chi1
        ("CB", "CG", "CD1"),
        ("CA", "CB", "CG", "CD1", "chi2"),  # chi2
        ("CB", "CG", "CD2"),
        ("CA", "CB", "CG", "CD2"),
        ("CB", "CG", "HG"),
        ("CA", "CB", "CG", "HG"),
        ("CG", "CD1", "HD11"),
        ("CG", "CD1", "HD12"),
        ("CG", "CD1", "HD13"),
        ("CG", "CD2", "HD21"),
        ("CG", "CD2", "HD22"),
        ("CG", "CD2", "HD23"),
        ("CB", "CG", "CD1", "HD11"),
        ("CB", "CG", "CD1", "HD12"),
        ("CB", "CG", "CD1", "HD13"),
        ("CB", "CG", "CD2", "HD21"),
        ("CB", "CG", "CD2", "HD22"),
        ("CB", "CG", "CD2", "HD23"),
    ),
    "I": (
        ("CA", "CB", "CG1"),
        ("N", "CA", "CB", "CG1", "chi1"),  # chi1
        ("CB", "CG1", "CD1"),
        ("CA", "CB", "CG1", "CD1", "chi2"),  # chi2
        ("CA", "CB", "CG2"),
        ("N", "CA", "CB", "CG2"),
        ("CB", "CG1", "HG12"),
        ("CB", "CG1", "HG13"),
        ("CB", "CG2", "HG21"),
        ("CB", "CG2", "HG22"),
        ("CB", "CG2", "HG23"),
        ("CA", "CB", "CG1", "HG12"),
        ("CA", "CB", "CG1", "HG13"),
        ("CA", "CB", "CG2", "HG21"),
        ("CA", "CB", "CG2", "HG22"),
        ("CA", "CB", "CG2", "HG23"),
        ("CG1", "CD1", "HD11"),
        ("CG1", "CD1", "HD12"),
        ("CG1", "CD1", "HD13"),
        ("CB", "CG1", "CD1", "HD11"),
        ("CB", "CG1", "CD1", "HD12"),
        ("CB", "CG1", "CD1", "HD13"),
    ),
    "M": (
        ("CA", "CB", "CG"),
        ("N", "CA", "CB", "CG", "chi1"),  # chi1
        ("CB", "CG", "SD"),
        ("CA", "CB", "CG", "SD", "chi2"),  # chi2
        ("CG", "SD", "CE"),
        ("CB", "CG", "SD", "CE", "chi3"),  # chi3
        ("CB", "CG", "HG2"),
        ("CB", "CG", "HG3"),
        ("CA", "CB", "CG", "HG2"),
        ("CA", "CB", "CG", "HG3"),
        ("SD", "CE", "HE1"),
        ("SD", "CE", "HE2"),
        ("SD", "CE", "HE3"),
        ("CG", "SD", "CE", "HE1"),
        ("CG", "SD", "CE", "HE2"),
        ("CG", "SD", "CE", "HE3"),
    ),
    "F": (
        ("CA", "CB", "CG"),
        ("N", "CA", "CB", "CG", "chi1"),  # chi1
        ("CB", "CG", "CD1"),
        ("CA", "CB", "CG", "CD1", "chi2"),  # chi2
        ("CG", "CD1", "CE1"),
        ("CB", "CG", "CD1", "CE1"),
        ("CD1", "CE1", "CZ"),
        ("CG", "CD1", "CE1", "CZ"),
        ("CB", "CG", "CD2"),
        ("CA", "CB", "CG", "CD2"),
        ("CG", "CD2", "CE2"),
        ("CB", "CG", "CD2", "CE2"),
        ("CG", "CD1", "HD1"),
        ("CB", "CG", "CD1", "HD1"),
        ("CG", "CD2", "HD2"),
        ("CB", "CG", "CD2", "HD2"),
        ("CD1", "CE1", "HE1"),
        ("CG", "CD1", "CE1", "HE1"),
        ("CD2", "CE2", "HE2"),
        ("CG", "CD2", "CE2", "HE2"),
        ("CE1", "CZ", "HZ"),
        ("CD1", "CE1", "CZ", "HZ"),
    ),
    "P": (
        ("CA", "CB", "CG"),
        ("N", "CA", "CB", "CG", "chi1"),  # chi1
        ("CB", "CG", "CD"),
        ("CA", "CB", "CG", "CD", "chi2"),  # chi2
        ("CB", "CG", "HG2"),
        ("CB", "CG", "HG3"),
        ("CA", "CB", "CG", "HG2"),
        ("CA", "CB", "CG", "HG3"),
        ("CG", "CD", "HD2"),
        ("CG", "CD", "HD3"),
        ("CB", "CG", "CD", "HD2"),
        ("CB", "CG", "CD", "HD3"),
    ),
    "S": (
        ("CA", "CB", "OG"),
        ("N", "CA", "CB", "OG", "chi1"),  # chi1
        ("CB", "OG", "HG"),
        ("CA", "CB", "OG", "HG"),
    ),
    "T": (
        ("CA", "CB", "OG1"),
        ("N", "CA", "CB", "OG1", "chi1"),  # chi1
        ("CA", "CB", "CG2"),
        ("N", "CA", "CB", "CG2"),
        ("CB", "OG1", "HG1"),
        ("CA", "CB", "OG1", "HG1"),
        ("CB", "CG2", "HG21"),
        ("CB", "CG2", "HG22"),
        ("CB", "CG2", "HG23"),
        ("CA", "CB", "CG2", "HG21"),
        ("CA", "CB", "CG2", "HG22"),
        ("CA", "CB", "CG2", "HG23"),
    ),
    "C": (
        ("CA", "CB", "SG"),
        ("N", "CA", "CB", "SG", "chi1"),  # chi1
        ("CB", "SG", "HG"),
        ("CA", "CB", "SG", "HG"),
    ),
    "N": (
        ("CA", "CB", "CG"),
        ("N", "CA", "CB", "CG", "chi1"),  # chi1
        ("CB", "CG", "OD1"),
        ("CA", "CB", "CG", "OD1", "chi2"),  # chi2
        ("CB", "CG", "ND2"),
        ("CA", "CB", "CG", "ND2"),
        ("CG", "ND2", "HD21"),
        ("CG", "ND2", "HD22"),
        ("CB", "CG", "ND2", "HD21"),
        ("CB", "CG", "ND2", "HD22"),
    ),
    "Q": (
        ("CA", "CB", "CG"),
        ("N", "CA", "CB", "CG", "chi1"),  # chi1
        ("CB", "CG", "CD"),
        ("CA", "CB", "CG", "CD", "chi2"),  # chi2
        ("CG", "CD", "OE1"),
        ("CB", "CG", "CD", "OE1", "chi3"),  # chi3
        ("CG", "CD", "NE2"),
        ("CB", "CG", "CD", "NE2"),
        ("CB", "CG", "HG2"),
        ("CB", "CG", "HG3"),
        ("CA", "CB", "CG", "HG2"),
        ("CA", "CB", "CG", "HG3"),
        ("CD", "NE2", "HE21"),
        ("CD", "NE2", "HE22"),
        ("CG", "CD", "NE2", "HE21"),
        ("CG", "CD", "NE2", "HE22"),
    ),
    "Y": (
        ("CA", "CB", "CG"),
        ("N", "CA", "CB", "CG", "chi1"),  # chi1
        ("CB", "CG", "CD1"),
        ("CA", "CB", "CG", "CD1", "chi2"),  # chi2
        ("CG", "CD1", "CE1"),
        ("CB", "CG", "CD1", "CE1"),
        ("CD1", "CE1", "CZ"),
        ("CG", "CD1", "CE1", "CZ"),
        ("CE1", "CZ", "OH"),
        ("CD1", "CE1", "CZ", "OH"),
        ("CB", "CG", "CD2"),
        ("CA", "CB", "CG", "CD2"),
        ("CG", "CD2", "CE2"),
        ("CB", "CG", "CD2", "CE2"),
        ("CG", "CD1", "HD1"),
        ("CB", "CG", "CD1", "HD1"),
        ("CG", "CD2", "HD2"),
        ("CB", "CG", "CD2", "HD2"),
        ("CD1", "CE1", "HE1"),
        ("CG", "CD1", "CE1", "HE1"),
        ("CD2", "CE2", "HE2"),
        ("CG", "CD2", "CE2", "HE2"),
        ("CZ", "OH", "HH"),
        ("CE1", "CZ", "OH", "HH"),
    ),
    "W": (
        ("CA", "CB", "CG"),
        ("N", "CA", "CB", "CG", "chi1"),  # chi1
        ("CB", "CG", "CD1"),
        ("CA", "CB", "CG", "CD1", "chi2"),  # chi2
        ("CG", "CD1", "NE1"),
        ("CB", "CG", "CD1", "NE1"),
        ("CB", "CG", "CD2"),
        ("CA", "CB", "CG", "CD2"),
        ("CG", "CD2", "CE2"),
        ("CB", "CG", "CD2", "CE2"),
        ("CD2", "CE2", "CZ2"),
        ("CG", "CD2", "CE2", "CZ2"),
        ("CE2", "CZ2", "CH2"),
        ("CD2", "CE2", "CZ2", "CH2"),
        ("CE2", "CZ2", "HZ2"),
        ("CD2", "CE2", "CZ2", "HZ2"),
        ("CG", "CD2", "CE3"),
        ("CB", "CG", "CD2", "CE3"),
        ("CZ2", "CH2", "CZ3"),
        ("CE2", "CZ2", "CH2", "CZ3"),
        ("CG", "CD1", "HD1"),
        ("CB", "CG", "CD1", "HD1"),
        ("CD1", "NE1", "HE1"),
        ("CG", "CD1", "NE1", "HE1"),
        ("CD2", "CE3", "HE3"),
        ("CG", "CD2", "CE3", "HE3"),
        ("CH2", "CZ3", "HZ3"),
        ("CZ2", "CH2", "CZ3", "HZ3"),
        ("CZ2", "CH2", "HH2"),
        ("CE2", "CZ2", "CH2", "HH2"),
    ),
    "D": (
        ("CA", "CB", "CG"),
        ("N", "CA", "CB", "CG", "chi1"),  # chi1
        ("CB", "CG", "OD1"),
        ("CA", "CB", "CG", "OD1", "chi2"),  # chi2
        ("CB", "CG", "OD2"),
        ("CA", "CB", "CG", "OD2"),
    ),
    "E": (
        ("CA", "CB", "CG"),
        ("N", "CA", "CB", "CG", "chi1"),  # chi1
        ("CB", "CG", "CD"),
        ("CA", "CB", "CG", "CD", "chi2"),  # chi2
        ("CG", "CD", "OE1"),
        ("CB", "CG", "CD", "OE1", "chi3"),  # chi3
        ("CG", "CD", "OE2"),
        ("CB", "CG", "CD", "OE2"),
        ("CB", "CG", "HG2"),
        ("CB", "CG", "HG3"),
        ("CA", "CB", "CG", "HG2"),
        ("CA", "CB", "CG", "HG3"),
    ),
    "H": (
        ("CA", "CB", "CG"),
        ("N", "CA", "CB", "CG", "chi1"),  # chi1
        ("CB", "CG", "ND1"),
        ("CA", "CB", "CG", "ND1", "chi2"),  # chi2
        ("CG", "ND1", "CE1"),
        ("CB", "CG", "ND1", "CE1"),
        ("CB", "CG", "CD2"),
        ("CA", "CB", "CG", "CD2"),
        ("CG", "CD2", "NE2"),
        ("CB", "CG", "CD2", "NE2"),
        ("CG", "ND1", "HD1"),
        ("CB", "CG", "ND1", "HD1"),
        ("CG", "CD2", "HD2"),
        ("CB", "CG", "CD2", "HD2"),
        ("ND1", "CE1", "HE1"),
        ("CG", "ND1", "CE1", "HE1"),
        ("CD2", "NE2", "HE2"),
        ("CG", "CD2", "NE2", "HE2"),
    ),
    "K": (
        ("CA", "CB", "CG"),
        ("N", "CA", "CB", "CG", "chi1"),  # chi1
        ("CB", "CG", "CD"),
        ("CA", "CB", "CG", "CD", "chi2"),  # chi2
        ("CG", "CD", "CE"),
        ("CB", "CG", "CD", "CE", "chi3"),  # chi3
        ("CD", "CE", "NZ"),
        ("CG", "CD", "CE", "NZ", "chi4"),  # chi4
        ("CB", "CG", "HG2"),
        ("CB", "CG", "HG3"),
        ("CA", "CB", "CG", "HG2"),
        ("CA", "CB", "CG", "HG3"),
        ("CG", "CD", "HD2"),
        ("CG", "CD", "HD3"),
        ("CB", "CG", "CD", "HD2"),
        ("CB", "CG", "CD", "HD3"),
        ("CD", "CE", "HE2"),
        ("CD", "CE", "HE3"),
        ("CG", "CD", "CE", "HE2"),
        ("CG", "CD", "CE", "HE3"),
        ("CE", "NZ", "HZ1"),
        ("CE", "NZ", "HZ2"),
        ("CE", "NZ", "HZ3"),
        ("CD", "CE", "NZ", "HZ1"),
        ("CD", "CE", "NZ", "HZ2"),
        ("CD", "CE", "NZ", "HZ3"),
    ),
    "R": (
        ("CA", "CB", "CG"),
        ("N", "CA", "CB", "CG", "chi1"),  # chi1
        ("CB", "CG", "CD"),
        ("CA", "CB", "CG", "CD", "chi2"),  # chi2
        ("CG", "CD", "NE"),
        ("CB", "CG", "CD", "NE", "chi3"),  # chi3
        ("CD", "NE", "CZ"),
        ("CG", "CD", "NE", "CZ", "chi4"),  # chi4
        ("NE", "CZ", "NH1"),
        ("CD", "NE", "CZ", "NH1", "chi5"),  # chi5
        ("NE", "CZ", "NH2"),
        ("CD", "NE", "CZ", "NH2"),
        ("CB", "CG", "HG2"),
        ("CB", "CG", "HG3"),
        ("CA", "CB", "CG", "HG2"),
        ("CA", "CB", "CG", "HG3"),
        ("CG", "CD", "HD2"),
        ("CG", "CD", "HD3"),
        ("CB", "CG", "CD", "HD2"),
        ("CB", "CG", "CD", "HD3"),
        ("CD", "NE", "HE"),
        ("CG", "CD", "NE", "HE"),
        ("CZ", "NH1", "HH11"),
        ("CZ", "NH1", "HH12"),
        ("NE", "CZ", "NH1", "HH11"),
        ("NE", "CZ", "NH1", "HH12"),
        ("CZ", "NH2", "HH21"),
        ("CZ", "NH2", "HH22"),
        ("NE", "CZ", "NH2", "HH21"),
        ("NE", "CZ", "NH2", "HH22"),
    ),
}

# Additional sidechain entries for explicit bonds.

# OpenSCAD output requires specification of bonds to be rendered as cylinders.
# These entries define hedra and dihedra to explicitly cover all bonds in rings,
# otherwise the entries above only capture atoms.

ic_data_sidechain_extras = {
    "F": (("CE1", "CZ", "CE2"), ("CD1", "CE1", "CZ", "CE2")),
    "P": (("CG", "CD", "N"), ("CB", "CG", "CD", "N")),
    "Y": (("CE1", "CZ", "CE2"), ("CD1", "CE1", "CZ", "CE2")),
    "W": (
        ("CD2", "CE3", "CZ3"),
        ("CG", "CD2", "CE3", "CZ3"),
        ("CD1", "NE1", "CE2"),
        ("CG", "CD1", "NE1", "CE2"),
    ),
    "H": (("ND1", "CE1", "NE2"), ("CG", "ND1", "CE1", "NE2")),
}

# Covalent radii for OpenSCAD output.

# Covalent radii from Heyrovska, Raji : 'Atomic Structures of all the Twenty
# Essential Amino Acids and a Tripeptide, with Bond Lengths as Sums of Atomic
# Covalent Radii' <https://arxiv.org/pdf/0804.2488.pdf>
# Adding Ores between Osb and Odb for Asp and Glu, Nres between Nsb and Ndb
# for Arg, as PDB does not specify

covalent_radii = {
    "Csb": 0.77,
    "Cres": 0.72,
    "Cdb": 0.67,
    "Osb": 0.67,
    "Ores": 0.635,
    "Odb": 0.60,
    "Nsb": 0.70,
    "Nres": 0.66,
    "Ndb": 0.62,
    "Hsb": 0.37,
    "Ssb": 1.04,
}

# Atom classes based on Heyrovska, Raji covalent radii paper.

residue_atom_bond_state = {
    "X": {
        "N": "Nsb",
        "CA": "Csb",
        "C": "Cdb",
        "O": "Odb",
        "OXT": "Osb",
        "CB": "Csb",
        "H": "Hsb",
    },
    "V": {"CG1": "Csb", "CG2": "Csb"},
    "L": {"CG": "Csb", "CD1": "Csb", "CD2": "Csb"},
    "I": {"CG1": "Csb", "CG2": "Csb", "CD1": "Csb"},
    "M": {"CG": "Csb", "SD": "Ssb", "CE": "Csb"},
    "F": {
        "CG": "Cdb",
        "CD1": "Cres",
        "CD2": "Cres",
        "CE1": "Cdb",
        "CE2": "Cdb",
        "CZ": "Cres",
    },
    "P": {"CG": "Csb", "CD": "Csb"},
    "S": {"OG": "Osb"},
    "T": {"OG1": "Osb", "CG2": "Csb"},
    "C": {"SG": "Ssb"},
    "N": {"CG": "Csb", "OD1": "Odb", "ND2": "Ndb"},
    "Q": {"CG": "Csb", "CD": "Csb", "OE1": "Odb", "NE2": "Ndb"},
    "Y": {
        "CG": "Cdb",
        "CD1": "Cres",
        "CD2": "Cres",
        "CE1": "Cdb",
        "CE2": "Cdb",
        "CZ": "Cres",
        "OH": "Osb",
    },
    "W": {
        "CG": "Cdb",
        "CD1": "Cdb",
        "CD2": "Cres",
        "NE1": "Nsb",
        "CE2": "Cdb",
        "CE3": "Cdb",
        "CZ2": "Cres",
        "CZ3": "Cres",
        "CH2": "Cdb",
    },
    "D": {"CG": "Csb", "OD1": "Ores", "OD2": "Ores"},
    "E": {"CG": "Csb", "CD": "Csb", "OE1": "Ores", "OE2": "Ores"},
    "H": {"CG": "Cdb", "CD2": "Cdb", "ND1": "Nsb", "CE1": "Cdb", "NE2": "Ndb"},
    "K": {"CG": "Csb", "CD": "Csb", "CE": "Csb", "NZ": "Nsb"},
    "R": {
        "CG": "Csb",
        "CD": "Csb",
        "NE": "Nsb",
        "CZ": "Cdb",
        "NH1": "Nres",
        "NH2": "Nres",
    },
}


# atomic weights of C,O,N,H,S
atomic_weight = {"C": 12.0107, "O": 15.9994, "N": 14.0067, "H": 1.0079, "S": 32.065}

# electronegativity values for C,O,N,H,S
electronegativity = {"C": 2.55, "O": 3.44, "N": 3.04, "H": 2.20, "S": 2.58}
