# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Output of PDB files."""


import warnings

# Exceptions and Warnings
from Bio import BiopythonWarning
from Bio.PDB.PDBExceptions import PDBIOException

# To allow saving of chains, residues, etc..
from Bio.PDB.StructureBuilder import StructureBuilder

# Allowed Elements
from Bio.Data.IUPACData import atom_weights


_ATOM_FORMAT_STRING = (
    "%s%5i %-4s%c%3s %c%4i%c   %8.3f%8.3f%8.3f%s%6.2f      %4s%2s%2s\n"
)
_PQR_ATOM_FORMAT_STRING = (
    "%s%5i %-4s%c%3s %c%4i%c   %8.3f%8.3f%8.3f %7s  %6s      %2s\n"
)

_TER_FORMAT_STRING = (
    "TER   %5i      %3s %c%4i%c                                                      \n"
)


class Select:
    """Select everything for PDB output (for use as a base class).

    Default selection (everything) during writing - can be used as base class
    to implement selective output. This selects which entities will be written out.
    """

    def __repr__(self):
        """Represent the output as a string for debugging."""
        return "<Select all>"

    def accept_model(self, model):
        """Overload this to reject models for output."""
        return 1

    def accept_chain(self, chain):
        """Overload this to reject chains for output."""
        return 1

    def accept_residue(self, residue):
        """Overload this to reject residues for output."""
        return 1

    def accept_atom(self, atom):
        """Overload this to reject atoms for output."""
        return 1


_select = Select()


class StructureIO:
    """Base class to derive structure file format writers from."""

    def __init__(self):
        """Initialise."""

    def set_structure(self, pdb_object):
        """Check what the user is providing and build a structure."""
        # The idea here is to build missing upstream components of
        # the SMCRA object representation. E.g., if the user provides
        # a Residue, build Structure/Model/Chain.

        if pdb_object.level == "S":
            structure = pdb_object
        else:  # Not a Structure
            sb = StructureBuilder()
            sb.init_structure("pdb")
            sb.init_seg(" ")

            if pdb_object.level == "M":
                sb.structure.add(pdb_object.copy())
                self.structure = sb.structure
            else:  # Not a Model
                sb.init_model(0)

                if pdb_object.level == "C":
                    sb.structure[0].add(pdb_object.copy())
                else:  # Not a Chain
                    chain_id = "A"  # default
                    sb.init_chain(chain_id)

                    if pdb_object.level == "R":  # Residue
                        # Residue extracted from a larger structure?
                        if pdb_object.parent is not None:
                            og_chain_id = pdb_object.parent.id
                            sb.structure[0][chain_id].id = og_chain_id
                            chain_id = og_chain_id

                        sb.structure[0][chain_id].add(pdb_object.copy())

                    else:  # Atom
                        sb.init_residue("DUM", " ", 1, " ")  # Dummy residue
                        sb.structure[0][chain_id].child_list[0].add(pdb_object.copy())

                        # Fix chain identifier if Atom has grandparents.
                        try:
                            og_chain_id = pdb_object.parent.parent.id
                        except AttributeError:  # pdb_object.parent == None
                            pass
                        else:
                            sb.structure[0][chain_id].id = og_chain_id

            # Return structure
            structure = sb.structure
        self.structure = structure


class PDBIO(StructureIO):
    """Write a Structure object (or a subset of a Structure object) as a PDB or PQR file.

    Examples
    --------
    >>> from Bio.PDB import PDBParser
    >>> from Bio.PDB.PDBIO import PDBIO
    >>> parser = PDBParser()
    >>> structure = parser.get_structure("1a8o", "PDB/1A8O.pdb")
    >>> io=PDBIO()
    >>> io.set_structure(structure)
    >>> io.save("bio-pdb-pdbio-out.pdb")
    >>> import os
    >>> os.remove("bio-pdb-pdbio-out.pdb")  # tidy up


    """

    def __init__(self, use_model_flag=0, is_pqr=False):
        """Create the PDBIO object.

        :param use_model_flag: if 1, force use of the MODEL record in output.
        :type use_model_flag: int
        :param is_pqr: if True, build PQR file. Otherwise build PDB file.
        :type is_pqr: Boolean
        """
        self.use_model_flag = use_model_flag
        self.is_pqr = is_pqr

    # private methods

    def _get_atom_line(
        self,
        atom,
        hetfield,
        segid,
        atom_number,
        resname,
        resseq,
        icode,
        chain_id,
        charge="  ",
    ):
        """Return an ATOM PDB string (PRIVATE)."""
        if hetfield != " ":
            record_type = "HETATM"
        else:
            record_type = "ATOM  "

        # Atom properties

        # Check if the atom serial number is an integer
        # Not always the case for structures built from
        # mmCIF files.
        try:
            atom_number = int(atom_number)
        except ValueError:
            raise ValueError(
                f"{atom_number!r} is not a number."
                "Atom serial numbers must be numerical"
                " If you are converting from an mmCIF"
                " structure, try using"
                " preserve_atom_numbering=False"
            )

        if atom_number > 99999:
            raise ValueError(
                f"Atom serial number ('{atom_number}') exceeds PDB format limit."
            )

        # Check if the element is valid, unknown (X), or blank
        if atom.element:
            element = atom.element.strip().upper()
            if element.capitalize() not in atom_weights and element != "X":
                raise ValueError(f"Unrecognised element {atom.element}")
            element = element.rjust(2)
        else:
            element = "  "

        # Format atom name
        # Pad if:
        #     - smaller than 4 characters
        # AND - is not C, N, O, S, H, F, P, ..., one letter elements
        # AND - first character is NOT numeric (funky hydrogen naming rules)
        name = atom.fullname.strip()
        if len(name) < 4 and name[:1].isalpha() and len(element.strip()) < 2:
            name = " " + name

        altloc = atom.altloc
        x, y, z = atom.coord

        # Write PDB format line
        if not self.is_pqr:
            bfactor = atom.bfactor
            try:
                occupancy = f"{atom.occupancy:6.2f}"
            except (TypeError, ValueError):
                if atom.occupancy is None:
                    occupancy = " " * 6
                    warnings.warn(
                        f"Missing occupancy in atom {atom.full_id!r} written as blank",
                        BiopythonWarning,
                    )
                else:
                    raise ValueError(
                        f"Invalid occupancy value: {atom.occupancy!r}"
                    ) from None

            args = (
                record_type,
                atom_number,
                name,
                altloc,
                resname,
                chain_id,
                resseq,
                icode,
                x,
                y,
                z,
                occupancy,
                bfactor,
                segid,
                element,
                charge,
            )
            return _ATOM_FORMAT_STRING % args

        # Write PQR format line
        else:
            try:
                pqr_charge = f"{atom.pqr_charge:7.4f}"
            except (TypeError, ValueError):
                if atom.pqr_charge is None:
                    pqr_charge = " " * 7
                    warnings.warn(
                        f"Missing PQR charge in atom {atom.full_id} written as blank",
                        BiopythonWarning,
                    )
                else:
                    raise ValueError(
                        f"Invalid PQR charge value: {atom.pqr_charge!r}"
                    ) from None

            try:
                radius = f"{atom.radius:6.4f}"
            except (TypeError, ValueError):
                if atom.radius is None:
                    radius = " " * 6
                    warnings.warn(
                        f"Missing radius in atom {atom.full_id} written as blank",
                        BiopythonWarning,
                    )
                else:
                    raise ValueError(f"Invalid radius value: {atom.radius}") from None

            args = (
                record_type,
                atom_number,
                name,
                altloc,
                resname,
                chain_id,
                resseq,
                icode,
                x,
                y,
                z,
                pqr_charge,
                radius,
                element,
            )

            return _PQR_ATOM_FORMAT_STRING % args

    # Public methods
    def save(self, file, select=_select, write_end=True, preserve_atom_numbering=False):
        """Save structure to a file.

        :param file: output file
        :type file: string or filehandle

        :param select: selects which entities will be written.
        :type select: object

        Typically select is a subclass of L{Select}, it should
        have the following methods:

         - accept_model(model)
         - accept_chain(chain)
         - accept_residue(residue)
         - accept_atom(atom)

        These methods should return 1 if the entity is to be
        written out, 0 otherwise.

        Typically select is a subclass of L{Select}.
        """
        if isinstance(file, str):
            fhandle = open(file, "w")
        else:
            # filehandle, I hope :-)
            fhandle = file

        get_atom_line = self._get_atom_line

        # multiple models?
        if len(self.structure) > 1 or self.use_model_flag:
            model_flag = 1
        else:
            model_flag = 0

        for model in self.structure.get_list():
            if not select.accept_model(model):
                continue
            # necessary for ENDMDL
            # do not write ENDMDL if no residues were written
            # for this model
            model_residues_written = 0
            if not preserve_atom_numbering:
                atom_number = 1
            if model_flag:
                fhandle.write(f"MODEL      {model.serial_num}\n")

            for chain in model.get_list():
                if not select.accept_chain(chain):
                    continue
                chain_id = chain.id
                if len(chain_id) > 1:
                    e = f"Chain id ('{chain_id}') exceeds PDB format limit."
                    raise PDBIOException(e)

                # necessary for TER
                # do not write TER if no residues were written
                # for this chain
                chain_residues_written = 0

                for residue in chain.get_unpacked_list():
                    if not select.accept_residue(residue):
                        continue
                    hetfield, resseq, icode = residue.id
                    resname = residue.resname
                    segid = residue.segid
                    resid = residue.id[1]
                    if resid > 9999:
                        e = f"Residue number ('{resid}') exceeds PDB format limit."
                        raise PDBIOException(e)

                    for atom in residue.get_unpacked_list():
                        if not select.accept_atom(atom):
                            continue
                        chain_residues_written = 1
                        model_residues_written = 1
                        if preserve_atom_numbering:
                            atom_number = atom.serial_number

                        try:
                            s = get_atom_line(
                                atom,
                                hetfield,
                                segid,
                                atom_number,
                                resname,
                                resseq,
                                icode,
                                chain_id,
                            )
                        except Exception as err:
                            # catch and re-raise with more information
                            raise PDBIOException(
                                f"Error when writing atom {atom.full_id}"
                            ) from err
                        else:
                            fhandle.write(s)
                            # inconsequential if preserve_atom_numbering is True
                            atom_number += 1

                if chain_residues_written:
                    fhandle.write(
                        _TER_FORMAT_STRING
                        % (atom_number, resname, chain_id, resseq, icode)
                    )

            if model_flag and model_residues_written:
                fhandle.write("ENDMDL\n")
        if write_end:
            fhandle.write("END   \n")

        if isinstance(file, str):
            fhandle.close()
