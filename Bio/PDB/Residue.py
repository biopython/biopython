# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.  

# My Stuff
from Bio.PDB.PDBExceptions import PDBConstructionException
from Bio.PDB.Entity import Entity, DisorderedEntityWrapper

# GSOC 2011 - ExtendedResidue
from Bio.Data import IUPACData
from Bio.SCOP.Raf import to_one_letter_code

"""Residue class, used by Structure objects."""


_atom_name_dict={}
_atom_name_dict["N"]=1
_atom_name_dict["CA"]=2
_atom_name_dict["C"]=3
_atom_name_dict["O"]=4


class Residue(Entity):
    """
    Represents a residue. A Residue object stores atoms.
    """
    def __init__(self, id, resname, segid):
        self.level="R"
        self.disordered=0
        self.extended=0
        self.resname=resname
        self.segid=segid
        Entity.__init__(self, id)

    # Special methods

    def __repr__(self):
        resname=self.get_resname()
        hetflag, resseq, icode=self.get_id()
        full_id=(resname, hetflag, resseq, icode)
        return "<Residue %s het=%s resseq=%s icode=%s>" % full_id

    # Private methods

    def _sort(self, a1, a2):
        """Sort the Atom objects.

        Atoms are sorted alphabetically according to their name, 
        but N, CA, C, O always come first.

        Arguments:
        o a1, a2 - Atom objects
        """
        name1=a1.get_name()
        name2=a2.get_name()
        if name1==name2:
            return(cmp(a1.get_altloc(), a2.get_altloc()))
        if name1 in _atom_name_dict:
            index1=_atom_name_dict[name1]
        else:
            index1=None
        if name2 in _atom_name_dict:
            index2=_atom_name_dict[name2]
        else:
            index2=None
        if index1 and index2:
            return cmp(index1, index2)
        if index1:
            return -1
        if index2:
            return 1
        return cmp(name1, name2)

    # Public methods

    def add(self, atom):
        """Add an Atom object.

        Checks for adding duplicate atoms, and raises a
        PDBConstructionException if so.
        """
        atom_id=atom.get_id()
        if self.has_id(atom_id):
            raise PDBConstructionException( \
                "Atom %s defined twice in residue %s" % (atom_id, self))
        Entity.add(self, atom)

    def sort(self):
        self.child_list.sort(self._sort)

    def flag_disordered(self):
        "Set the disordered flag."
        self.disordered=1

    def is_disordered(self):
        "Return 1 if the residue contains disordered atoms."
        return self.disordered

    def get_resname(self):
        return self.resname

    def get_unpacked_list(self):
        """
        Returns the list of all atoms, unpack DisorderedAtoms."
        """
        atom_list=self.get_list()
        undisordered_atom_list=[]
        for atom in atom_list:
            if atom.is_disordered():
                undisordered_atom_list=(undisordered_atom_list+ atom.disordered_get_list())
            else:
                undisordered_atom_list.append(atom)
        return undisordered_atom_list       

    def get_segid(self):
        return self.segid

    def extend(self):
        self.__class__=ExtendedResidue
        self.__init__()
        return self

class DisorderedResidue(DisorderedEntityWrapper):
    """
    DisorderedResidue is a wrapper around two or more Residue objects. It is
    used to represent point mutations (e.g. there is a Ser 60 and a Cys 60 residue,
    each with 50 % occupancy).
    """
    def __init__(self, id):
        DisorderedEntityWrapper.__init__(self, id)

    def __repr__(self):
        resname=self.get_resname()
        hetflag, resseq, icode=self.get_id()
        full_id=(resname, hetflag, resseq, icode)
        return "<DisorderedResidue %s het=%s resseq=%i icode=%s>" % full_id

    def add(self, atom):
        residue=self.disordered_get()
        if not atom.is_disordered()==2:
            # Atoms in disordered residues should have non-blank
            # altlocs, and are thus represented by DisorderedAtom objects.
            resname=residue.get_resname()
            het, resseq, icode=residue.get_id() 
            # add atom anyway, if PDBParser ignores exception the atom will be part of the residue
            residue.add(atom)
            raise PDBConstructionException( \
                "Blank altlocs in duplicate residue %s (%s, %i, %s)" \
                % (resname, het, resseq, icode) )
        residue.add(atom)

    def sort(self):
        "Sort the atoms in the child Residue objects."
        for residue in self.disordered_get_list():
            residue.sort() 

    def disordered_add(self, residue):
        """Add a residue object and use its resname as key.

        Arguments:
        o residue - Residue object
        """
        resname=residue.get_resname()
        # add chain parent to residue
        chain=self.get_parent()
        residue.set_parent(chain)
        assert(not self.disordered_has_id(resname))
        self[resname]=residue
        self.disordered_select(resname)

class ExtendedResidue(Residue):
    """
    ExtendedResidue adds some physical and chemical information
    to a Residue object.

    - Hydrophobicity
    - etc etc etc

    """

    def __init__(self):
#        Residue.__init__(residue.id, residue.resname, residue.segid)
        self.extended=1
        self.set_hydrophobicity()
        self.set_charge()
        self.set_mass()
        
    def __repr__(self):
        resname=self.get_resname()
        hetflag, resseq, icode=self.get_id()
        full_id=(resname, hetflag, resseq, icode)
        return "<ExtendedResidue %s het=%s resseq=%i icode=%s>" % full_id

    def set_hydrophobicity(self, scale='consensus'):
        """
        Sets the hydrophobicity scale to use:
        - Kyte and Doolittle  (KD)
        - Sweet and Eisenberg (OHM)
        - Consensus of the previous two. (consensus) [default].
        
        Returns scale name.
        """

        set_scale=getattr(IUPACData, "protein_hydropathy_"+scale)
        
        self.hydrophobicity=set_scale[to_one_letter_code[self.resname]]
        self.hydrophobicity_scale = scale

        return scale

    def set_charge(self, pH=7.0):
        """
        Sets the charge depending on the pH (by default 7.0)
        and the pKa of the side chain        
        """

        try: 
            set_pka=getattr(IUPACData, "protein_pka_side_chain")
            pka=set_pka[to_one_letter_code[self.resname]]
        except IndexError:
            warnings.warn("pka not found for residue %s. See Data/IUPACData for a list of supported residues" %self.resname )
            self.charge = None
            return self.charge

        if pka and pka < pH:
            charge=1
        elif pka and pka > pH:
            charge=-1
        else:
            charge=0

        self.charge=charge
        
        return charge

    def set_mass(self):
        """
        Sets the mass of a residue from its atoms, if missing information, takes theorical weight
        """
        try:
            mass=sum([atom.mass for atom in self])
            self.mass=mass
        except:
            warnings.warn("An mass atom is missing. Using the theorical residue mass")
            theo_mass=getattr(IUPACData, "protein_weights")
            if (to_one_letter_code[self.resname] in theo_mass):
                self.mass=theo_mass[to_one_letter_code[self.resname]]
            else:
                warnings.warn("Residue %s has no theorical mass. Mass sets to None" %self.resname)
                self.mass=None
        return mass
