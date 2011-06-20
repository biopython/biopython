# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Superimpose two structures."""

import numpy

from Bio.SVDSuperimposer import SVDSuperimposer
from Bio.PDB.PDBExceptions import PDBException


class Superimposer(object):
    """
    Rotate/translate one set of atoms on top of another,
    thereby minimizing the RMSD.
    """
    def __init__(self):
        self.rotran=None
        self.rms=None

    def set_atoms(self, fixed, moving):
        """
        Put (translate/rotate) the atoms in fixed on the atoms in 
        moving, in such a way that the RMSD is minimized.

        @param fixed: list of (fixed) atoms
        @param moving: list of (moving) atoms 
        @type fixed,moving: [L{Atom}, L{Atom},...]
        """
        if not (len(fixed)==len(moving)):
            raise PDBException("Fixed and moving atom lists differ in size")
        l=len(fixed)
        fixed_coord=numpy.zeros((l, 3))
        moving_coord=numpy.zeros((l, 3))
        for i in range(0, len(fixed)):
            fixed_coord[i]=fixed[i].get_coord()
            moving_coord[i]=moving[i].get_coord()
        sup=SVDSuperimposer()
        sup.set(fixed_coord, moving_coord)
        sup.run()
        self.rms=sup.get_rms()
        self.rotran=sup.get_rotran()

    def apply(self, atom_list):
        """
        Rotate/translate a list of atoms.
        """
        if self.rotran is None:
            raise PDBException("No transformation has been calculated yet")
        rot, tran=self.rotran
        rot=rot.astype('f')
        tran=tran.astype('f')
        for atom in atom_list:
            atom.transform(rot, tran)


if __name__=="__main__":
    import sys

    from Bio.PDB import PDBParser, Selection

    p=PDBParser()
    s1=p.get_structure("FIXED", sys.argv[1])
    fixed=Selection.unfold_entities(s1, "A")

    s2=p.get_structure("MOVING", sys.argv[1])
    moving=Selection.unfold_entities(s2, "A")

    rot=numpy.identity(3).astype('f')
    tran=numpy.array((1.0, 2.0, 3.0), 'f')

    for atom in moving:
        atom.transform(rot, tran)
    
    sup=Superimposer()

    sup.set_atoms(fixed, moving)

    print sup.rotran
    print sup.rms

    sup.apply(moving)
