# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.  

"""The structure class, representing a macromolecular structure."""

from Bio.PDB.Entity import Entity


class Structure(Entity):
    """
    The Structure class contains a collection of Model instances.
    """
    def __init__(self, id):
        self.level="S"
        Entity.__init__(self, id)

    # Special methods

    def __repr__(self):
        return "<Structure id=%s>" % self.get_id()

    # Private methods

    def _sort(self, m1, m2):
        """Sort models.

        This sorting function sorts the Model instances in the Structure instance.
        The sorting is done based on the model id, which is a simple int that 
        reflects the order of the models in the PDB file.

        Arguments:
        o m1, m2 - Model instances
        """
        return cmp(m1.get_id(), m2.get_id())

    # Public 

    def get_chains(self):
        for m in self:
            for c in m:
                yield c

    def get_residues(self):
        for c in self.get_chains():
            for r in c:
                yield r

    def get_atoms(self):
        for r in self.get_residues():
            for a in r:
                yield a
        

