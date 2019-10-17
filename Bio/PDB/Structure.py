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
        for m in self:
            yield m

    def get_chains(self):
        """Return chains from models."""
        for m in self.get_models():
            for c in m:
                yield c

    def get_residues(self):
        """Return residues from chains."""
        for c in self.get_chains():
            for r in c:
                yield r

    def get_atoms(self):
        """Return atoms from residue."""
        for r in self.get_residues():
            for a in r:
                yield a

    def atom_to_internal_coordinates(self, allBonds=False):
        """Create/update internal coordinates from Atom X,Y,Z coordinates.

        Internal coordinates are bond length, angle and dihedral angles.

        :param allBonds bool: default False
        include hedra and dihedra for bonds around sidechain rings.
        (not required to capture all atoms)
        """
        for chn in self.get_chains():
            chn.atom_to_internal_coordinates(allBonds)

    def internal_to_atom_coordinates(self):
        """Create/update atom coordinates from internal coordinates.

        :raises Exception: if any chain does not have .pic attribute
        """
        for chn in self.get_chains():
            chn.internal_to_atom_coordinates()
    
