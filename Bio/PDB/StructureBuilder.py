# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Consumer class that builds a Structure object.

This is used by the PDBParser and MMCIFparser classes.
"""

import warnings

# SMCRA hierarchy
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue, DisorderedResidue
from Bio.PDB.Atom import Atom, DisorderedAtom

from Bio.PDB.PDBExceptions import PDBConstructionException
from Bio.PDB.PDBExceptions import PDBConstructionWarning


class StructureBuilder:
    """Deals with constructing the Structure object.

    The StructureBuilder class is used by the PDBParser classes to
    translate a file to a Structure object.
    """

    def __init__(self):
        """Initialize the class."""
        self.line_counter = 0
        self.header = {}

    def _is_completely_disordered(self, residue):
        """Return 1 if all atoms in the residue have a non blank altloc (PRIVATE)."""
        atom_list = residue.get_unpacked_list()
        for atom in atom_list:
            altloc = atom.get_altloc()
            if altloc == " ":
                return 0
        return 1

    # Public methods called by the Parser classes

    def set_header(self, header):
        """Set header."""
        self.header = header

    def set_line_counter(self, line_counter):
        """Tracks line in the PDB file that is being parsed.

        Arguments:
         - line_counter - int

        """
        self.line_counter = line_counter

    def init_structure(self, structure_id):
        """Initialize a new Structure object with given id.

        Arguments:
         - id - string

        """
        self.structure = Structure(structure_id)

    def init_model(self, model_id, serial_num=None):
        """Create a new Model object with given id.

        Arguments:
         - id - int
         - serial_num - int

        """
        self.model = Model(model_id, serial_num)
        self.structure.add(self.model)

    def init_chain(self, chain_id):
        """Create a new Chain object with given id.

        Arguments:
         - chain_id - string

        """
        if self.model.has_id(chain_id):
            self.chain = self.model[chain_id]
            warnings.warn(
                "WARNING: Chain %s is discontinuous at line %i."
                % (chain_id, self.line_counter),
                PDBConstructionWarning,
            )
        else:
            self.chain = Chain(chain_id)
            self.model.add(self.chain)

    def init_seg(self, segid):
        """Flag a change in segid.

        Arguments:
         - segid - string

        """
        self.segid = segid

    def init_residue(self, resname, field, resseq, icode):
        """Create a new Residue object.

        Arguments:
         - resname - string, e.g. "ASN"
         - field - hetero flag, "W" for waters, "H" for
           hetero residues, otherwise blank.
         - resseq - int, sequence identifier
         - icode - string, insertion code

        """
        if field != " ":
            if field == "H":
                # The hetero field consists of H_ + the residue name (e.g. H_FUC)
                field = "H_" + resname
        res_id = (field, resseq, icode)
        if field == " ":
            if self.chain.has_id(res_id):
                # There already is a residue with the id (field, resseq, icode).
                # This only makes sense in the case of a point mutation.
                warnings.warn(
                    "WARNING: Residue ('%s', %i, '%s') redefined at line %i."
                    % (field, resseq, icode, self.line_counter),
                    PDBConstructionWarning,
                )
                duplicate_residue = self.chain[res_id]
                if duplicate_residue.is_disordered() == 2:
                    # The residue in the chain is a DisorderedResidue object.
                    # So just add the last Residue object.
                    if duplicate_residue.disordered_has_id(resname):
                        # The residue was already made
                        self.residue = duplicate_residue
                        duplicate_residue.disordered_select(resname)
                    else:
                        # Make a new residue and add it to the already
                        # present DisorderedResidue
                        new_residue = Residue(res_id, resname, self.segid)
                        duplicate_residue.disordered_add(new_residue)
                        self.residue = duplicate_residue
                        return
                else:
                    if resname == duplicate_residue.resname:
                        warnings.warn(
                            "WARNING: Residue ('%s', %i, '%s','%s') already defined "
                            "with the same name at line  %i."
                            % (field, resseq, icode, resname, self.line_counter),
                            PDBConstructionWarning,
                        )
                        self.residue = duplicate_residue
                        return
                    # Make a new DisorderedResidue object and put all
                    # the Residue objects with the id (field, resseq, icode) in it.
                    # These residues each should have non-blank altlocs for all their atoms.
                    # If not, the PDB file probably contains an error.
                    if not self._is_completely_disordered(duplicate_residue):
                        # if this exception is ignored, a residue will be missing
                        self.residue = None
                        raise PDBConstructionException(
                            "Blank altlocs in duplicate residue %s ('%s', %i, '%s')"
                            % (resname, field, resseq, icode)
                        )
                    self.chain.detach_child(res_id)
                    new_residue = Residue(res_id, resname, self.segid)
                    disordered_residue = DisorderedResidue(res_id)
                    self.chain.add(disordered_residue)
                    disordered_residue.disordered_add(duplicate_residue)
                    disordered_residue.disordered_add(new_residue)
                    self.residue = disordered_residue
                    return
        self.residue = Residue(res_id, resname, self.segid)
        self.chain.add(self.residue)

    def init_atom(
        self,
        name,
        coord,
        b_factor,
        occupancy,
        altloc,
        fullname,
        serial_number=None,
        element=None,
        pqr_charge=None,
        radius=None,
        is_pqr=False,
    ):
        """Create a new Atom object.

        Arguments:
         - name - string, atom name, e.g. CA, spaces should be stripped
         - coord - Numeric array (Float0, size 3), atomic coordinates
         - b_factor - float, B factor
         - occupancy - float
         - altloc - string, alternative location specifier
         - fullname - string, atom name including spaces, e.g. " CA "
         - element - string, upper case, e.g. "HG" for mercury
         - pqr_charge - float, atom charge (PQR format)
         - radius - float, atom radius (PQR format)
         - is_pqr - boolean, flag to specify if a .pqr file is being parsed

        """
        residue = self.residue
        # if residue is None, an exception was generated during
        # the construction of the residue
        if residue is None:
            return
        # First check if this atom is already present in the residue.
        # If it is, it might be due to the fact that the two atoms have atom
        # names that differ only in spaces (e.g. "CA.." and ".CA.",
        # where the dots are spaces). If that is so, use all spaces
        # in the atom name of the current atom.
        if residue.has_id(name):
            duplicate_atom = residue[name]
            # atom name with spaces of duplicate atom
            duplicate_fullname = duplicate_atom.get_fullname()
            if duplicate_fullname != fullname:
                # name of current atom now includes spaces
                name = fullname
                warnings.warn(
                    "Atom names %r and %r differ only in spaces at line %i."
                    % (duplicate_fullname, fullname, self.line_counter),
                    PDBConstructionWarning,
                )
        if not is_pqr:
            self.atom = Atom(
                name,
                coord,
                b_factor,
                occupancy,
                altloc,
                fullname,
                serial_number,
                element,
            )
        elif is_pqr:
            self.atom = Atom(
                name,
                coord,
                None,
                None,
                altloc,
                fullname,
                serial_number,
                element,
                pqr_charge,
                radius,
            )
        if altloc != " ":
            # The atom is disordered
            if residue.has_id(name):
                # Residue already contains this atom
                duplicate_atom = residue[name]
                if duplicate_atom.is_disordered() == 2:
                    duplicate_atom.disordered_add(self.atom)
                else:
                    # This is an error in the PDB file:
                    # a disordered atom is found with a blank altloc
                    # Detach the duplicate atom, and put it in a
                    # DisorderedAtom object together with the current
                    # atom.
                    residue.detach_child(name)
                    disordered_atom = DisorderedAtom(name)
                    residue.add(disordered_atom)
                    disordered_atom.disordered_add(self.atom)
                    disordered_atom.disordered_add(duplicate_atom)
                    residue.flag_disordered()
                    warnings.warn(
                        "WARNING: disordered atom found with blank altloc before "
                        "line %i.\n" % self.line_counter,
                        PDBConstructionWarning,
                    )
            else:
                # The residue does not contain this disordered atom
                # so we create a new one.
                disordered_atom = DisorderedAtom(name)
                residue.add(disordered_atom)
                # Add the real atom to the disordered atom, and the
                # disordered atom to the residue
                disordered_atom.disordered_add(self.atom)
                residue.flag_disordered()
        else:
            # The atom is not disordered
            residue.add(self.atom)

    def set_anisou(self, anisou_array):
        """Set anisotropic B factor of current Atom."""
        self.atom.set_anisou(anisou_array)

    def set_siguij(self, siguij_array):
        """Set standard deviation of anisotropic B factor of current Atom."""
        self.atom.set_siguij(siguij_array)

    def set_sigatm(self, sigatm_array):
        """Set standard deviation of atom position of current Atom."""
        self.atom.set_sigatm(sigatm_array)

    def get_structure(self):
        """Return the structure."""
        # first sort everything
        # self.structure.sort()
        # Add the header dict
        self.structure.header = self.header
        return self.structure

    def set_symmetry(self, spacegroup, cell):
        """Set symmetry."""
        pass
