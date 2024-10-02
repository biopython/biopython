"""
A tool that introduces a point mutation to a given residue based on the Dunbrack rotamer library.
The tool uses sample residues stored in small PDB files, and uses a reduced rotamer library
to save space (only rotamers with higher than 5% probability are present).
The tool makes in place modification on a PDB structure object.
"""

from Bio.Entrez.Parser import ValidationError
from Bio.PDB import Polypeptide
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Selection import unfold_entities
from Bio.PDB.vectors import calc_dihedral, Vector, rotaxis2m
from Bio import SVDSuperimposer
from collections import defaultdict
from Bio.PDB.ic_data import (
    CHI_ANGLES,
    RESIDUE_ORDER,
    VW_RADII,
    residue_atom_bond_state,
    ic_data_backbone,
)
from Bio.SeqUtils import seq3
import numpy as np
import os
import gzip

DATA_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "rotamers")
DUNBRACK_FILTERED = f"{DATA_DIR}/dunbrack_reduced.gz"


class Rotamer:
    """
    Object to store rotamers
    """

    def __init__(self, resname, psi, phi, chi1, chi2, chi3, chi4, probability):
        """
        Init the rotamer with given diherdral angles
        """
        self.psi = float(psi)
        self.phi = float(phi)
        self.chi1 = float(chi1)
        self.chi2 = float(chi2)
        self.chi3 = float(chi3)
        self.chi4 = float(chi4)
        self.resname = resname
        self.probability = float(probability)
        self.sample_residue = False

    @property
    def sample_residue(self):
        """Store a sample residue to use for the mutation tool in the Residue class"""
        if self._sample_residue:
            return self._sample_residue
        raise ValueError("Sample residue is not set, use RotamerSampling to set it!")

    @sample_residue.setter
    def sample_residue(self, sample_residue):
        """Property setter for the sample residue"""
        self._sample_residue = sample_residue


class RotamerSampling:
    """
    Object to sample rotamers based on certain criterion to use for mutation
    """

    def __init__(self, structure_object, rotamer_loc=DUNBRACK_FILTERED):
        """
        Load the Dunbrack rotamer library
        :param rotamer_loc: Location of the rotamer file to use
        :returns : Nested defaultdict of the torsion angles for each psi phi angle in the file
        """
        self.structure_object = structure_object
        self._rot_lib = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        with gzip.open(rotamer_loc, "rt") as fn:
            for line in fn:
                if line.startswith("#"):
                    continue
                amino_acid, angle1, angle2, prob, chi1, chi2, chi3, chi4 = line.split()
                rotamer = Rotamer(
                    amino_acid, angle1, angle2, chi1, chi2, chi3, chi4, prob
                )
                self._rot_lib[amino_acid][int(angle1)][int(angle2)].append(rotamer)

    def select_based_on_clashes(
        self,
        residue,
        mutate_to,
        sample_residue,
        rotamers,
        skip_own_chain=False,
    ):
        """
        Calculates the best rotamer based on VdW energy
        :param residue: Bio.PDB Residue object
        :param mutate_to: Three letter amino acid code in uppercase letters
        :param rotamers: Rotamers to use for selection. Obtained from rotamer_lib
        :param skip_own_chain: Boolean to skip the parent chain of the selected residue. Default: False
        :returns: Best rotamer from rotamers based on simplified Van der Waals energy.
        """
        best_rotamer = None
        lowest_energy = float("inf")
        for rotamer in rotamers:
            vdw_energy = 0
            # Introduce the rotamer
            for angle in ["CHI1", "CHI2", "CHI3", "CHI4"]:
                if mutate_to not in CHI_ANGLES[angle]:
                    continue
                dihedral_start = calc_dihedral(
                    *[
                        Vector(sample_residue[x])
                        for x in CHI_ANGLES[angle][mutate_to]["ref_plane"]
                    ]
                )
                rotation_angle = dihedral_start - np.deg2rad(
                    getattr(rotamer, angle.lower())
                )
                axis = CHI_ANGLES[angle][mutate_to]["axis"]
                for atom in RESIDUE_ORDER[mutate_to][
                    RESIDUE_ORDER[mutate_to].index(axis[1]) + 1 :
                ]:
                    sample_residue[atom] = (
                        np.dot(
                            rotaxis2m(
                                rotation_angle,
                                Vector(
                                    sample_residue[axis[0]] - sample_residue[axis[1]]
                                ),
                            ),
                            sample_residue[atom] - sample_residue[axis[1]],
                        )
                        + sample_residue[axis[1]]
                    )

            for rotamer_atom, rotamer_vector in sample_residue.items():
                if vdw_energy > lowest_energy:  # Skip pointless rotamers
                    break
                atoms = unfold_entities(self.structure_object, "A")
                ns = NeighborSearch(atoms)
                close_atoms = ns.search(rotamer_vector, 5)  # 5 Angstrom radii
                for close_atom in close_atoms:
                    close_residue = close_atom.get_parent()
                    chain_if_close_atom = close_residue.get_parent()
                    if (
                        close_atom.get_parent().get_id()[1] == residue.get_id()[1]
                        and chain_if_close_atom.get_id()
                        == residue.get_parent().get_id()
                    ):
                        continue

                    if (
                        skip_own_chain
                        and chain_if_close_atom.get_id()
                        == residue.get_parent().get_id()
                    ):
                        continue
                    if (
                        abs(close_residue.get_id()[1] - residue.get_id()[1]) == 1
                        and close_atom.is_backbone
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

    @staticmethod
    def read_sample_residue(residue_name):
        """
        Reads a sample PDB residue
        :param residue_name 3 letter name of a residue
        :returns: Formatted residue to use for mutation
        """
        sample_residue = {}
        with open(os.path.join(DATA_DIR, "sample_residues.pdb")) as fn:
            for line in fn:
                if line[17:20].strip() != residue_name.upper():
                    continue
                sample_residue[line[12:16].strip()] = np.array(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                )
        return sample_residue

    def sample(self, residue, mutate_to, method="best"):
        """
        Generate sample from a rotamer library
        :param residue: Bio.PBD.Residue object
        :param mutate_to: Three letter amino acid code in uppercase letters
        :param method: How to select the best rotamer in string format. Default: best
        best: Select the best rotamer based on simplified Van der Waals energy
        random: Select a random rotamer based on its probability distribution
        bestother: Same as best, but ignore clashes within the chain of the modified residue
        """
        sample_residue = self.read_sample_residue(mutate_to)
        starting_points = np.asmatrix(
            [sample_residue["N"], sample_residue["CA"], sample_residue["C"]]
        )
        end_points = np.asmatrix(
            [residue["N"].coord, residue["CA"].coord, residue["C"].coord]
        )

        sup = SVDSuperimposer.SVDSuperimposer()
        sup.set(end_points, starting_points)
        sup.run()
        rot, tran = sup.get_rotran()

        for atom, coords in sample_residue.items():
            sample_residue[atom] = np.squeeze(np.asarray(np.dot(coords, rot) + tran))
        Polypeptide.Polypeptide(
            residue.get_parent()
        ).get_phi_psi_list()  # This generates the xtra attribute for PHI and PSI angles for residues

        for atom in list(residue.get_atoms()):  # Create a copy to remove from
            if not atom.is_backbone:
                atom.parent.detach_child(atom.id)
        phi, psi = (
            round(np.rad2deg(y), -1) if y else 0
            for y in [residue.xtra[x] for x in ["PHI", "PSI"]]
        )
        if method == "first":
            selected_rotamer = sorted(
                self._rot_lib[mutate_to][phi][psi],
                key=lambda x: x["prob"],
                reverse=True,
            )[0]
        elif method == "random":
            p = np.array([x["prob"] for x in self._rot_lib[mutate_to][phi][psi]])
            p /= p.sum()
            selected_rotamer = np.random.choice(self._rot_lib[mutate_to][phi][psi], p=p)
        elif method == "best":
            selected_rotamer = self.select_based_on_clashes(
                residue,
                mutate_to,
                sample_residue,
                self._rot_lib[mutate_to][phi][psi],
            )
        elif method == "bestother":
            selected_rotamer = self.select_based_on_clashes(
                residue,
                mutate_to,
                sample_residue,
                self._rot_lib[mutate_to][phi][psi],
                skip_own_chain=True,
            )
        else:
            raise ValueError(
                f"Unknown mutation type {method}. Possible choices are 'first', 'random', 'best', 'bestother'"
            )
        selected_rotamer.sample_residue = sample_residue
        return selected_rotamer


def distance(x, y):
    """
    Calculates the euclidian distance in 3D space
    :param x, y: Points in 3D space in a form of vectors
    :returns: Distance of points
    """
    return np.sqrt((x[0] - y[0]) ** 2 + (x[1] - y[1]) ** 2 + (x[2] - y[2]) ** 2)
