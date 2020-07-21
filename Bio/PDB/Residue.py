# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Residue class, used by Structure objects."""

# My Stuff
import warnings
from Bio import BiopythonDeprecationWarning
from Bio.PDB.PDBExceptions import PDBConstructionException
from Bio.PDB.Entity import Entity, DisorderedEntityWrapper


_atom_name_dict = {}
_atom_name_dict["N"] = 1
_atom_name_dict["CA"] = 2
_atom_name_dict["C"] = 3
_atom_name_dict["O"] = 4


class Residue(Entity):
    """Represents a residue. A Residue object stores atoms."""

    def __init__(self, id, resname, segid):
        """Initialize the class."""
        self.level = "R"
        self.disordered = 0
        self.resname = resname
        self.segid = segid
        self.internal_coord = None
        Entity.__init__(self, id)

    def __repr__(self):
        """Return the residue full id."""
        resname = self.get_resname()
        hetflag, resseq, icode = self.get_id()
        full_id = (resname, hetflag, resseq, icode)
        return "<Residue %s het=%s resseq=%s icode=%s>" % full_id

    def add(self, atom):
        """Add an Atom object.

        Checks for adding duplicate atoms, and raises a
        PDBConstructionException if so.
        """
        atom_id = atom.get_id()
        if self.has_id(atom_id):
            raise PDBConstructionException(
                "Atom %s defined twice in residue %s" % (atom_id, self)
            )
        Entity.add(self, atom)

    def sort(self):
        """Sort child atoms.

        Atoms N, CA, C, O always come first, thereafter alphabetically
        by name, with any alternative location specifier for disordered
        atoms (altloc) as a tie-breaker.
        """
        warnings.warn(
            "The custom sort() method will be removed in the future in favour of rich "
            "comparison methods. Use the built-in sorted() function instead.",
            BiopythonDeprecationWarning,
        )
        self.child_list.sort()

    def flag_disordered(self):
        """Set the disordered flag."""
        self.disordered = 1

    def is_disordered(self):
        """Return 1 if the residue contains disordered atoms."""
        return self.disordered

    def get_resname(self):
        """Return the residue name."""
        return self.resname

    def get_unpacked_list(self):
        """Return the list of all atoms, unpack DisorderedAtoms."""
        atom_list = self.get_list()
        undisordered_atom_list = []
        for atom in atom_list:
            if atom.is_disordered():
                undisordered_atom_list += atom.disordered_get_list()
            else:
                undisordered_atom_list.append(atom)
        return undisordered_atom_list

    def get_segid(self):
        """Return the segment identifier."""
        return self.segid

    def get_atoms(self):
        """Return atoms."""
        yield from self

    def get_atom(self):
        """Return atom."""
        warnings.warn(
            "`get_atom` has been deprecated and we intend to remove it"
            " in a future release of Biopython. Please use `get_atoms` instead.",
            BiopythonDeprecationWarning,
        )
        yield from self


class DisorderedResidue(DisorderedEntityWrapper):
    """DisorderedResidue is a wrapper around two or more Residue objects.

    It is used to represent point mutations (e.g. there is a Ser 60 and a Cys 60
    residue, each with 50 % occupancy).
    """

    def __init__(self, id):
        """Initialize the class."""
        DisorderedEntityWrapper.__init__(self, id)

    def __repr__(self):
        """Return disordered residue full identifier."""
        resname = self.get_resname()
        hetflag, resseq, icode = self.get_id()
        full_id = (resname, hetflag, resseq, icode)
        return "<DisorderedResidue %s het=%s resseq=%i icode=%s>" % full_id

    def add(self, atom):
        """Add atom to residue."""
        residue = self.disordered_get()
        if not atom.is_disordered() == 2:
            # Atoms in disordered residues should have non-blank
            # altlocs, and are thus represented by DisorderedAtom objects.
            resname = residue.get_resname()
            het, resseq, icode = residue.get_id()
            # add atom anyway, if PDBParser ignores exception the atom will be part of the residue
            residue.add(atom)
            raise PDBConstructionException(
                "Blank altlocs in duplicate residue %s (%s, %i, %s)"
                % (resname, het, resseq, icode)
            )
        residue.add(atom)

    def sort(self):
        """Sort the atoms in the child Residue objects."""
        for residue in self.disordered_get_list():
            residue.sort()

    def disordered_add(self, residue):
        """Add a residue object and use its resname as key.

        Arguments:
         - residue - Residue object

        """
        resname = residue.get_resname()
        # add chain parent to residue
        chain = self.get_parent()
        residue.set_parent(chain)
        assert not self.disordered_has_id(resname)
        self[resname] = residue
        self.disordered_select(resname)
