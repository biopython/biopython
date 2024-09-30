"""
A tool that introduces a point mutation to a given residue based on the Dunbrack rotamer library.
The tool uses sample residues stored in small PDB files, and uses a reduced rotamer library
to save space (only rotamers with higher than 5% probability are present).
The tool makes in place modification on a PDB structure object.
"""

from Bio.PDB.Atom import Atom
from Bio.PDB import Polypeptide
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Selection import unfold_entities
from Bio import SVDSuperimposer
from collections import defaultdict
import numpy as np
import os
import gzip


DATA_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "rotamers")
DUNBRACK_FILTERED = f"{DATA_DIR}/dunbrack_reduced.gz"
VW_RADII = {
    "ALA": {"N": 1.7, "CA": 2.0, "C": 1.7, "O": 1.4, "CB": 2.0},
    "CYS": {"N": 1.7, "CA": 2.0, "C": 1.7, "O": 1.4, "CB": 2.0, "SG": 1.8},
    "ASP": {
        "N": 1.7,
        "CA": 2.0,
        "C": 1.7,
        "O": 1.4,
        "CB": 2.0,
        "CG": 1.7,
        "OD1": 1.5,
        "OD2": 1.5,
    },
    "GLU": {
        "N": 1.7,
        "CA": 2.0,
        "C": 1.7,
        "O": 1.4,
        "CB": 2.0,
        "CG": 2.0,
        "CD": 1.7,
        "OE1": 1.5,
        "OE2": 1.5,
    },
    "PHE": {
        "N": 1.7,
        "CA": 2.0,
        "C": 1.7,
        "O": 1.4,
        "CB": 2.0,
        "CG": 1.7,
        "CD1": 1.9,
        "CD2": 1.9,
        "CE1": 1.9,
        "CE2": 1.9,
        "CZ": 1.9,
    },
    "GLY": {"N": 1.7, "CA": 2.0, "C": 1.7, "O": 1.4},
    "HIS": {
        "N": 1.7,
        "CA": 2.0,
        "C": 1.7,
        "O": 1.4,
        "CB": 2.0,
        "CG": 1.7,
        "ND1": 1.7,
        "CD2": 1.9,
        "CE1": 1.9,
        "NE2": 1.7,
    },
    "ILE": {
        "N": 1.7,
        "CA": 2.0,
        "C": 1.7,
        "O": 1.4,
        "CB": 2.0,
        "CG1": 2.0,
        "CG2": 2.0,
        "CD1": 2.0,
    },
    "LYS": {
        "N": 1.7,
        "CA": 2.0,
        "C": 1.7,
        "O": 1.4,
        "CB": 2.0,
        "CG": 2.0,
        "CD": 2.0,
        "CE": 2.0,
        "NZ": 2.0,
    },
    "LEU": {
        "N": 1.7,
        "CA": 2.0,
        "C": 1.7,
        "O": 1.4,
        "CB": 2.0,
        "CG": 2.0,
        "CD1": 2.0,
        "CD2": 2.0,
    },
    "MET": {
        "N": 1.7,
        "CA": 2.0,
        "C": 1.7,
        "O": 1.4,
        "CB": 2.0,
        "CG": 2.0,
        "SD": 1.8,
        "CE": 2.0,
    },
    "ASN": {
        "N": 1.7,
        "CA": 2.0,
        "C": 1.7,
        "O": 1.4,
        "CB": 2.0,
        "CG": 1.7,
        "OD1": 1.6,
        "ND2": 1.6,
    },
    "PRO": {"N": 1.7, "CA": 2.0, "C": 1.7, "O": 1.4, "CB": 2.0, "CG": 2.0, "CD": 2.0},
    "GLN": {
        "N": 1.7,
        "CA": 2.0,
        "C": 1.7,
        "O": 1.4,
        "CB": 2.0,
        "CG": 2.0,
        "CD": 1.7,
        "OE1": 1.6,
        "NE2": 1.6,
    },
    "ARG": {
        "N": 1.7,
        "CA": 2.0,
        "C": 1.7,
        "O": 1.4,
        "CB": 2.0,
        "CG": 2.0,
        "CD": 2.0,
        "NE": 1.7,
        "CZ": 2.0,
        "NH1": 2.0,
        "NH2": 2.0,
    },
    "SER": {"N": 1.7, "CA": 2.0, "C": 1.7, "O": 1.4, "CB": 2.0, "OG": 1.6},
    "THR": {"N": 1.7, "CA": 2.0, "C": 1.7, "O": 1.4, "CB": 2.0, "OG1": 1.6, "CG2": 2.0},
    "VAL": {"N": 1.7, "CA": 2.0, "C": 1.7, "O": 1.4, "CB": 2.0, "CG1": 2.0, "CG2": 2.0},
    "TRP": {
        "N": 1.7,
        "CA": 2.0,
        "C": 1.7,
        "O": 1.4,
        "CB": 2.0,
        "CG": 1.7,
        "CD1": 1.9,
        "CD2": 1.7,
        "NE1": 1.7,
        "CE2": 1.7,
        "CE3": 1.9,
        "CZ2": 1.9,
        "CZ3": 1.9,
        "CH2": 1.9,
    },
    "TYR": {
        "N": 1.7,
        "CA": 2.0,
        "C": 1.7,
        "O": 1.4,
        "CB": 2.0,
        "CG": 1.7,
        "CD1": 1.9,
        "CD2": 1.9,
        "CE1": 1.9,
        "CE2": 1.9,
        "CZ": 1.7,
        "OH": 1.6,
    },
}

CHI_ANGLES = {
    "CHI1": {
        "CYS": {"axis": ["CA", "CB"], "ref_plane": ["N", "CA", "CB", "SG"]},
        "ASP": {"axis": ["CA", "CB"], "ref_plane": ["N", "CA", "CB", "CG"]},
        "SER": {"axis": ["CA", "CB"], "ref_plane": ["N", "CA", "CB", "OG"]},
        "GLN": {"axis": ["CA", "CB"], "ref_plane": ["N", "CA", "CB", "CG"]},
        "LYS": {"axis": ["CA", "CB"], "ref_plane": ["N", "CA", "CB", "CG"]},
        "ILE": {"axis": ["CA", "CB"], "ref_plane": ["N", "CA", "CB", "CG1"]},
        "PRO": {"axis": ["CA", "CB"], "ref_plane": ["N", "CA", "CB", "CG"]},
        "THR": {"axis": ["CA", "CB"], "ref_plane": ["N", "CA", "CB", "OG1"]},
        "PHE": {"axis": ["CA", "CB"], "ref_plane": ["N", "CA", "CB", "CG"]},
        "ASN": {"axis": ["CA", "CB"], "ref_plane": ["N", "CA", "CB", "CG"]},
        "HIS": {"axis": ["CA", "CB"], "ref_plane": ["N", "CA", "CB", "CG"]},
        "LEU": {"axis": ["CA", "CB"], "ref_plane": ["N", "CA", "CB", "CG"]},
        "ARG": {"axis": ["CA", "CB"], "ref_plane": ["N", "CA", "CB", "CG"]},
        "TRP": {"axis": ["CA", "CB"], "ref_plane": ["N", "CA", "CB", "CG"]},
        "VAL": {"axis": ["CA", "CB"], "ref_plane": ["N", "CA", "CB", "CG1"]},
        "GLU": {"axis": ["CA", "CB"], "ref_plane": ["N", "CA", "CB", "CG"]},
        "TYR": {"axis": ["CA", "CB"], "ref_plane": ["N", "CA", "CB", "CG"]},
        "MET": {"axis": ["CA", "CB"], "ref_plane": ["N", "CA", "CB", "CG"]},
    },
    "CHI2": {
        "ASP": {"axis": ["CB", "CG"], "ref_plane": ["CA", "CB", "CG", "OD1"]},
        "GLN": {"axis": ["CB", "CG"], "ref_plane": ["CA", "CB", "CG", "CD"]},
        "LYS": {"axis": ["CB", "CG"], "ref_plane": ["CA", "CB", "CG", "CD"]},
        "ILE": {"axis": ["CB", "CG1"], "ref_plane": ["CA", "CB", "CG1", "CD1"]},
        "PRO": {"axis": ["CB", "CG"], "ref_plane": ["CA", "CB", "CG", "CD"]},
        "PHE": {"axis": ["CB", "CG"], "ref_plane": ["CA", "CB", "CG", "CD1"]},
        "ASN": {"axis": ["CB", "CG"], "ref_plane": ["CA", "CB", "CG", "OD1"]},
        "HIS": {"axis": ["CB", "CG"], "ref_plane": ["CA", "CB", "CG", "ND1"]},
        "LEU": {"axis": ["CB", "CG"], "ref_plane": ["CA", "CB", "CG", "CD1"]},
        "ARG": {"axis": ["CB", "CG"], "ref_plane": ["CA", "CB", "CG", "CD"]},
        "TRP": {"axis": ["CB", "CG"], "ref_plane": ["CA", "CB", "CG", "CD1"]},
        "GLU": {"axis": ["CB", "CG"], "ref_plane": ["CA", "CB", "CG", "CD"]},
        "TYR": {"axis": ["CB", "CG"], "ref_plane": ["CA", "CB", "CG", "CD1"]},
        "MET": {"axis": ["CB", "CG"], "ref_plane": ["CA", "CB", "CG", "SD"]},
    },
    "CHI3": {
        "ARG": {"axis": ["CG", "CD"], "ref_plane": ["CB", "CG", "CD", "NE"]},
        "GLN": {"axis": ["CG", "CD"], "ref_plane": ["CB", "CG", "CD", "OE1"]},
        "GLU": {"axis": ["CG", "CD"], "ref_plane": ["CB", "CG", "CD", "OE1"]},
        "LYS": {"axis": ["CG", "CD"], "ref_plane": ["CB", "CG", "CD", "CE"]},
        "MET": {"axis": ["CG", "SD"], "ref_plane": ["CB", "CG", "SD", "CE"]},
    },
    "CHI4": {
        "ARG": {"axis": ["CD", "NE"], "ref_plane": ["CG", "CD", "NE", "CZ"]},
        "LYS": {"axis": ["CG", "CE"], "ref_plane": ["CG", "CD", "CE", "NZ"]},
    },
}

RESIDUE_ORDER = {
    "CYS": ["N", "CA", "C", "O", "CB", "SG"],
    "ASP": ["N", "CA", "C", "O", "CB", "CG", "OD1", "OD2"],
    "SER": ["N", "CA", "C", "O", "CB", "OG"],
    "GLN": ["N", "CA", "C", "O", "CB", "CG", "CD", "NE2", "OE1"],
    "LYS": ["N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"],
    "ILE": ["N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"],
    "PRO": ["N", "CA", "C", "O", "CB", "CG", "CD"],
    "THR": ["N", "CA", "C", "O", "CB", "CG2", "OG1"],
    "PHE": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "ASN": ["N", "CA", "C", "O", "CB", "CG", "ND2", "OD1"],
    "GLY": ["N", "CA", "C", "O"],
    "HIS": ["N", "CA", "C", "O", "CB", "CG", "CD2", "ND1", "CE1", "NE2"],
    "LEU": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2"],
    "ARG": ["N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"],
    "TRP": [
        "N",
        "CA",
        "C",
        "O",
        "CB",
        "CG",
        "CD1",
        "CD2",
        "CE2",
        "CE3",
        "NE1",
        "CZ2",
        "CZ3",
        "CH2",
    ],
    "ALA": ["N", "CA", "C", "O", "CB"],
    "VAL": ["N", "CA", "C", "O", "CB", "CG1", "CG2"],
    "GLU": ["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2"],
    "TYR": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"],
    "MET": ["N", "CA", "C", "O", "CB", "CG", "SD", "CE"],
}


def load_rotamers(rotamer_loc=DUNBRACK_FILTERED):
    """
    Load the Dunbrack rotamer library
    """
    _dunbrack = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    with gzip.open(rotamer_loc, "rt") as fn:
        for line in fn:
            if line.startswith("#"):
                continue
            amino_acid, angle1, angle2, prob, chi1, chi2, chi3, chi4 = line.split()
            _dunbrack[amino_acid][int(angle1)][int(angle2)].append(
                {
                    "prob": float(prob),
                    "CHI1": float(chi1),
                    "CHI2": float(chi2),
                    "CHI3": float(chi3),
                    "CHI4": float(chi4),
                }
            )
    return _dunbrack


def rotation_matrix(axis, theta):
    """
    Calculates a rotational matrix around an axis with an angle theta
    """
    axis = np.asarray(axis)
    axis = axis / np.linalg.norm(axis)

    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)

    outer_product = np.outer(axis, axis)

    cross_product_matrix = np.array(
        [[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]]
    )

    rotation_mat = (
        cos_theta * np.eye(3)
        + sin_theta * cross_product_matrix
        + (1 - cos_theta) * outer_product
    )

    return rotation_mat


def dihedral_from_vectors(v1, v2, v3, v4):
    """Praxeolitic formula
    1 sqrt, 1 cross product"""
    b0 = -1.0 * (v2 - v1)
    b1 = v3 - v2
    b2 = v4 - v3

    b1 /= np.linalg.norm(b1)

    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1

    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.arctan2(y, x)


def distance(x, y):
    """
    Calculates the euclidian distance in 3D space
    """
    return np.sqrt((x[0] - y[0]) ** 2 + (x[1] - y[1]) ** 2 + (x[2] - y[2]) ** 2)


def read_sample_residue(residue_name):
    """
    Reads a sample PDB residue
    """
    sample_residue = {}
    with open(f"{DATA_DIR}/{residue_name.upper()}.pdb") as fn:
        for line in fn:
            sample_residue[line[12:16].strip()] = np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            )
    return sample_residue


def is_backbone(atom):
    """
    Checks whether a given residue is backbone or not
    """
    return atom.get_id() in ["C", "N", "CA", "O"]


def select_best_rotamer_based_on_clashes(
    pdb_object,
    chain,
    res_num,
    mutate_to,
    sample_residue,
    rotamers,
    skip_own_chain=False,
):
    """
    Calculates the best rotamer based on VdW energy
    """
    best_rotamer = None
    lowest_energy = float("inf")
    for rotamer in rotamers:
        vdw_energy = 0
        # Introduce the rotamer
        for angle in ["CHI1", "CHI2", "CHI3", "CHI4"]:
            if mutate_to not in CHI_ANGLES[angle]:
                continue
            dihedral_start = dihedral_from_vectors(
                *[sample_residue[x] for x in CHI_ANGLES[angle][mutate_to]["ref_plane"]]
            )
            rotation_angle = dihedral_start - np.deg2rad(rotamer[angle])
            axis = CHI_ANGLES[angle][mutate_to]["axis"]
            # print(angle)
            for atom in RESIDUE_ORDER[mutate_to][
                RESIDUE_ORDER[mutate_to].index(axis[1]) + 1 :
            ]:
                sample_residue[atom] = (
                    np.dot(
                        rotation_matrix(
                            sample_residue[axis[0]] - sample_residue[axis[1]],
                            rotation_angle,
                        ),
                        sample_residue[atom] - sample_residue[axis[1]],
                    )
                    + sample_residue[axis[1]]
                )

        for rotamer_atom, rotamer_vector in sample_residue.items():
            if vdw_energy > lowest_energy:  # Skip pointless rotamers
                break
            atoms = unfold_entities(pdb_object[0], "A")
            ns = NeighborSearch(atoms)
            close_atoms = ns.search(rotamer_vector, 5)  # 5 Angstrom radii
            for close_atom in close_atoms:
                close_residue = close_atom.get_parent()
                chain_if_close_atom = close_residue.get_parent()
                if (
                    close_atom.get_parent().get_id()[1] == res_num
                    and chain_if_close_atom.get_id() == chain
                ):
                    # print("skipping_own_atom")
                    continue

                if skip_own_chain and chain_if_close_atom.get_id() == chain:
                    # print("Skipping own chain")
                    continue
                if abs(close_residue.get_id()[1] - res_num) == 1 and is_backbone(
                    close_atom
                ):
                    continue
                dist = distance(close_atom.coord, rotamer_vector)
                if dist > 6:
                    continue
                try:
                    vdw_radi = (
                        VW_RADII[close_atom.get_parent().get_resname()][
                            close_atom.get_id()
                        ]
                        + VW_RADII[mutate_to][rotamer_atom]
                    )
                except KeyError:
                    continue
                attractive_force = (vdw_radi / dist) ** 6
                vdw_energy += attractive_force**2 - attractive_force
        if vdw_energy < lowest_energy:
            lowest_energy = vdw_energy
            best_rotamer = rotamer
    return best_rotamer


def mutate(
    pdb_obj,
    chain,
    res_num,
    mutate_to,
    rotamer_lib=None,
    mutation_type="best",
):
    """
    Mutates a given residue based on the best fitting rotamers
    """
    Polypeptide.Polypeptide(
        pdb_obj[0][chain]
    ).get_phi_psi_list()  # This generates the xtra attribute for PHI and PSI angles for residues
    try:
        _residue = pdb_obj[0][chain][res_num]
    except KeyError:
        raise KeyError(f"Residue {res_num} not found in chain {chain}!")
    for atom in list(_residue.get_atoms()):  # Create a copy to remove from
        if not is_backbone(atom):
            atom.parent.detach_child(atom.id)
    phi, psi = (
        round(np.rad2deg(y), -1) if y else 0
        for y in [_residue.xtra[x] for x in ["PHI", "PSI"]]
    )
    # GET_ATR IS WRONG
    sample_residue = read_sample_residue(mutate_to)
    starting_points = np.asmatrix(
        [sample_residue["N"], sample_residue["CA"], sample_residue["C"]]
    )
    end_points = np.asmatrix(
        [_residue["N"].coord, _residue["CA"].coord, _residue["C"].coord]
    )

    sup = SVDSuperimposer.SVDSuperimposer()
    sup.set(end_points, starting_points)
    sup.run()
    rot, tran = sup.get_rotran()

    for atom, coords in sample_residue.items():
        sample_residue[atom] = np.squeeze(np.asarray(np.dot(coords, rot) + tran))

    if mutate_to not in ["ALA", "GLY"]:
        if not rotamer_lib:
            rotamer_lib = load_rotamers()
        if mutation_type == "first":
            selected_rotamer = sorted(
                rotamer_lib[mutate_to][phi][psi], key=lambda x: x["prob"], reverse=True
            )[0]
        elif mutation_type == "random":
            p = np.array([x["prob"] for x in rotamer_lib[mutate_to][phi][psi]])
            p /= p.sum()
            selected_rotamer = np.random.choice(rotamer_lib[mutate_to][phi][psi], p=p)
        elif mutation_type == "best":
            selected_rotamer = select_best_rotamer_based_on_clashes(
                pdb_obj,
                chain,
                res_num,
                mutate_to,
                sample_residue,
                rotamer_lib[mutate_to][phi][psi],
            )
        elif mutation_type == "bestother":
            selected_rotamer = select_best_rotamer_based_on_clashes(
                pdb_obj,
                chain,
                res_num,
                mutate_to,
                sample_residue,
                rotamer_lib[mutate_to][phi][psi],
                skip_own_chain=True,
            )
        else:
            raise ValueError(
                f"Unknown mutation type {mutation_type}. Possible choices are 'first', 'random', 'best', 'bestother'"
            )

        # Introduce the rotamer
        for angle in ["CHI1", "CHI2", "CHI3", "CHI4"]:
            if mutate_to not in CHI_ANGLES[angle]:
                continue
            dihedral_start = dihedral_from_vectors(
                *[sample_residue[x] for x in CHI_ANGLES[angle][mutate_to]["ref_plane"]]
            )
            rotation_angle = dihedral_start - np.deg2rad(selected_rotamer[angle])
            axis = CHI_ANGLES[angle][mutate_to]["axis"]
            # print(angle)
            for atom in RESIDUE_ORDER[mutate_to][
                RESIDUE_ORDER[mutate_to].index(axis[1]) + 1 :
            ]:
                sample_residue[atom] = (
                    np.dot(
                        rotation_matrix(
                            sample_residue[axis[0]] - sample_residue[axis[1]],
                            rotation_angle,
                        ),
                        sample_residue[atom] - sample_residue[axis[1]],
                    )
                    + sample_residue[axis[1]]
                )
    for atom, coord in sample_residue.items():
        if atom not in ["C", "N", "CA", "O"]:
            new_atom = Atom(
                name=atom,
                element=atom[0],
                fullname="{}{}".format(
                    " " * (4 - len(atom)), atom
                ),  # for writing the structure, should be 4-char long
                coord=np.asarray(coord),
                bfactor=1.0,
                altloc=" ",
                occupancy=1.0,
                serial_number=9999,  # does not matter much, only for writing the struct.
            )
            _residue.add(new_atom)
    _residue.resname = mutate_to
    return
