# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""The structure class, representing a macromolecular structure."""

from typing import NamedTuple

from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Entity import Entity


class DisulfideBond(NamedTuple):
    """Stores the relevant information of a SSBOND line."""

    atom1: Atom
    atom2: Atom
    distance: float
    serial_number: str
    insertion_code1: str
    symmetry_operator1: int
    insertion_code2: str
    symmetry_operator2: int


class Link(NamedTuple):
    """Stores the relevant information of a LINK line."""

    atom1: Atom
    atom2: Atom
    distance: float
    alternate_location1: str
    insertion_code1: str
    symmetry_operator1: int
    alternate_location2: str
    insertion_code2: str
    symmetry_operator2: int


class Structure(Entity):
    """The Structure class contains a collection of Model instances."""

    def __init__(self, id):
        """Initialize the class."""
        self.level = "S"
        Entity.__init__(self, id)

    def __repr__(self):
        """Return the structure identifier."""
        return f"<Structure id={self.get_id()}>"

    def get_models(self):
        """Return models."""
        yield from self

    def get_chains(self):
        """Return chains from models."""
        for m in self.get_models():
            yield from m

    def get_residues(self):
        """Return residues from chains."""
        for c in self.get_chains():
            yield from c

    def get_atoms(self):
        """Return atoms from residue."""
        for r in self.get_residues():
            yield from r

    def atom_to_internal_coordinates(self, verbose: bool = False) -> None:
        """Create/update internal coordinates from Atom X,Y,Z coordinates.

        Internal coordinates are bond length, angle and dihedral angles.

        :param verbose bool: default False
            describe runtime problems

        """
        for chn in self.get_chains():
            chn.atom_to_internal_coordinates(verbose)

    def internal_to_atom_coordinates(self, verbose: bool = False) -> None:
        """Create/update atom coordinates from internal coordinates.

        :param verbose bool: default False
            describe runtime problems

        :raises Exception: if any chain does not have .internal_coord attribute
        """
        for chn in self.get_chains():
            chn.internal_to_atom_coordinates(verbose)

    def _find_residue_by_id(self, chain: Chain, id_: int):
        """
        Return the residue with `id_` from the `chain`.

        This is a workaround required due to Biopython's residue indexing e.g.
        ('H_GLC', 100, 'A') instead of 100. For more info please refer to the
        section "What is a residue id?" at
        https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ
        """
        for residue in chain.get_residues():
            if residue.id[1] == id_:
                return residue

    def get_links(self) -> list[Link]:
        """
        Return intramolecular links of the Structure.

        Each link consists of two Atom objects and their associated distance.
        """
        raw_links = self.header["links"]
        links = []
        for link in raw_links:
            for model in self.get_models():
                chain1 = model[link.res1.chain]
                residue1 = self._find_residue_by_id(chain1, link.res1.id)
                atom1 = residue1[link.res1.atom]

                chain2 = model[link.res2.chain]
                residue1 = self._find_residue_by_id(chain2, link.res2.id)
                atom2 = residue1[link.res2.atom]

                links.append(
                    Link(
                        atom1=atom1,
                        atom2=atom2,
                        alternate_location1=link.alternate_location1,
                        alternate_location2=link.alternate_location2,
                        symmetry_operator1=link.symmetry_operator1,
                        symmetry_operator2=link.symmetry_operator2,
                        insertion_code1=link.insertion_code1,
                        insertion_code2=link.insertion_code2,
                        distance=link.distance,
                    )
                )

        return links

    def get_disulfide_bonds(self) -> list[DisulfideBond]:
        """
        Return disulfide bonds as DisulfideBond objects.

        Each bond consists of two Atom objects and their associated distance.
        """
        raw_ss_bonds = self.header["ss_bonds"]
        disulfide_bonds = []
        for ss_bond in raw_ss_bonds:
            for model in self.get_models():

                residue1 = model[ss_bond.res1.chain][ss_bond.res1.id]
                residue2 = model[ss_bond.res2.chain][ss_bond.res2.id]
                atoms = [
                    atom
                    for i, residue in enumerate([residue1, residue2])
                    for atom in residue.get_atoms()
                    if "S" in atom.name  # sulfur atoms could be S, SD, SG, etc.
                ]

                disulfide_bonds.append(
                    DisulfideBond(
                        atom1=atoms[0],
                        atom2=atoms[1],
                        serial_number=ss_bond.serial_number,
                        symmetry_operator1=ss_bond.symmetry_operator1,
                        symmetry_operator2=ss_bond.symmetry_operator2,
                        insertion_code1=ss_bond.insertion_code1,
                        insertion_code2=ss_bond.insertion_code2,
                        distance=ss_bond.distance,
                    )
                )

        return disulfide_bonds
