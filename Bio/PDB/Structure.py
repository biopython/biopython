# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""The structure class, representing a macromolecular structure."""

from Bio.PDB.Entity import Entity


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
