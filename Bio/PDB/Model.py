# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Model class, used in Structure objects."""

from Bio.PDB.Entity import Entity


class Model(Entity):
    """The object representing a model in a structure.

    In a structure derived from an X-ray crystallography experiment,
    only a single model will be present (with some exceptions). NMR
    structures normally contain many different models.
    """

    def __init__(self, id, serial_num=None):
        """Initialize.

        Arguments:
         - id - int
         - serial_num - int

        """
        self.level = "M"
        if serial_num is None:
            self.serial_num = id
        else:
            self.serial_num = serial_num

        Entity.__init__(self, id)

    def __repr__(self):
        return "<Model id=%s>" % self.get_id()

    def get_chains(self):
        for c in self:
            yield c

    def get_residues(self):
        for c in self.get_chains():
            for r in c:
                yield r

    def get_atoms(self):
        for r in self.get_residues():
            for a in r:
                yield a
