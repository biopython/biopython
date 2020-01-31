# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""The structure class, representing a macromolecular structure."""

from Bio.PDB.Entity import Entity
from Bio.PDB.internal_coords import IC_Chain


class Structure(Entity):
    """The Structure class contains a collection of Model instances."""

    def __init__(self, id):
        """Initialize the class."""
        self.level = "S"
        Entity.__init__(self, id)

    def __repr__(self):
        """Return the structure identifier."""
        return "<Structure id=%s>" % self.get_id()

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

        :raises Exception: if any chain does not have .pic attribute
        """
        for chn in self.get_chains():
            chn.internal_to_atom_coordinates(verbose)
