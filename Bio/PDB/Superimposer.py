# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Superimpose two structures."""


import numpy as np

from Bio.SVDSuperimposer import SVDSuperimposer
from Bio.PDB.PDBExceptions import PDBException


class Superimposer:
    """Rotate/translate one set of atoms on top of another to minimize RMSD."""

    def __init__(self):
        """Initialize the class."""
        self.rotran = None
        self.rms = None

    def set_atoms(self, fixed, moving):
        """Prepare translation/rotation to minimize RMSD between atoms.

        Put (translate/rotate) the atoms in fixed on the atoms in
        moving, in such a way that the RMSD is minimized.

        :param fixed: list of (fixed) atoms
        :param moving: list of (moving) atoms
        :type fixed,moving: [L{Atom}, L{Atom},...]
        """
        if not len(fixed) == len(moving):
            raise PDBException("Fixed and moving atom lists differ in size")
        length = len(fixed)
        fixed_coord = np.zeros((length, 3))
        moving_coord = np.zeros((length, 3))
        for i in range(length):
            fixed_coord[i] = fixed[i].get_coord()
            moving_coord[i] = moving[i].get_coord()
        sup = SVDSuperimposer()
        sup.set(fixed_coord, moving_coord)
        sup.run()
        self.rms = sup.get_rms()
        self.rotran = sup.get_rotran()

    def apply(self, atom_list):
        """Rotate/translate a list of atoms."""
        if self.rotran is None:
            raise PDBException("No transformation has been calculated yet")
        rot, tran = self.rotran
        rot = rot.astype("f")
        tran = tran.astype("f")
        for atom in atom_list:
            atom.transform(rot, tran)
