# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Output of PDB files."""


from itertools import groupby
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

_SSBOND_FORMAT_STRING = (
    "SSBOND %3i %3s %1s %4i%1s   %3s %1s %4i%1s                       %6i %6i %5.2f\n"
)

_LINK_FORMAT_STRING = "LINK        %-4s%1s%3s %1s%4i%1s               %-4s%1s%3s %1s%4i%1s  %6i %6i %5.2f\n"

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
        pass

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

        # Ensure chain id isn't longer than 1 character
        if len(chain_id) > 1:
            raise ValueError(f"Chain ID must be of length 1: {chain_id!r} > 1")

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

        element, name = self._extract_element_and_name(atom)

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

    @staticmethod
    def _extract_element_and_name(atom):
        """Extract element and atom name as strings."""
        # Check if the element is valid, unknown (X), or blank
        if atom.element:
            element = atom.element.strip().upper()
            if element.capitalize() not in atom_weights and element != "X":
                raise ValueError(f"Unrecognised element {atom.element}")
            element = element.rjust(2)
        else:
            element = "  "

        name = atom.fullname.strip()
        if len(name) < 4 and name[:1].isalpha() and len(element.strip()) < 2:
            name = " " + name

        return element, name

    def _get_link_line(self, link):
        """Return a LINK PDB string (PRIVATE)."""
        _, atom_name1 = self._extract_element_and_name(link.atom1)
        altloc1 = link.alternate_location1
        res_name1 = link.atom1.parent.resname
        chain_id1 = link.atom1.parent.parent.id
        res_id1 = link.atom1.parent.id[1]
        icode1 = link.insertion_code1
        _, atom_name2 = self._extract_element_and_name(link.atom2)
        altloc2 = link.alternate_location2
        res_name2 = link.atom2.parent.resname
        chain_id2 = link.atom2.parent.parent.id
        res_id2 = link.atom2.parent.id[1]
        icode2 = link.insertion_code2
        symmetry_operator1 = link.symmetry_operator1
        symmetry_operator2 = link.symmetry_operator2
        distance = link.distance

        args = (
            atom_name1,
            altloc1,
            res_name1,
            chain_id1,
            res_id1,
            icode1,
            atom_name2,
            altloc2,
            res_name2,
            chain_id2,
            res_id2,
            icode2,
            symmetry_operator1,
            symmetry_operator2,
            distance,
        )
        return _LINK_FORMAT_STRING % args

    def _get_ssbond_line(self, ssbond):
        """Return a SSBOND PDB string (PRIVATE)."""
        serial_number = ssbond.serial_number
        res_name1 = ssbond.atom1.parent.resname
        chain_id1 = ssbond.atom1.parent.parent.id
        res_id1 = ssbond.atom1.parent.id[1]
        icode1 = ssbond.insertion_code1
        res_name2 = ssbond.atom2.parent.resname
        chain_id2 = ssbond.atom2.parent.parent.id
        res_id2 = ssbond.atom2.parent.id[1]
        icode2 = ssbond.insertion_code2
        symmetry_operator1 = ssbond.symmetry_operator1
        symmetry_operator2 = ssbond.symmetry_operator2
        distance = ssbond.distance

        args = (
            serial_number,
            res_name1,
            chain_id1,
            res_id1,
            icode1,
            res_name2,
            chain_id2,
            res_id2,
            icode2,
            symmetry_operator1,
            symmetry_operator2,
            distance,
        )
        return _SSBOND_FORMAT_STRING % args

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

        with fhandle:
            get_atom_line = self._get_atom_line

            # multiple models?
            if len(self.structure) > 1 or self.use_model_flag:
                model_flag = 1
            else:
                model_flag = 0

            # only apply SSBONDs once independent of the number of models present
            ssbond_groups = groupby(
                self.structure.get_disulfide_bonds(),
                lambda x: (x.atom1.parent.id[1], x.atom2.parent.id[1]),
            )
            for _, ssbonds in ssbond_groups:
                ssbond = next(ssbonds)
                try:
                    line = self._get_ssbond_line(ssbond)
                except Exception as err:
                    # catch and re-raise with more information
                    raise PDBIOException(f"Error when writing SSBOND {ssbond}") from err
                else:
                    fhandle.write(line)

            # only apply LINKs once independent of the number of models present
            link_groups = groupby(
                self.structure.get_links(),
                lambda x: (
                    x.atom1.id,
                    x.atom1.parent.id[1],
                    x.atom2.id,
                    x.atom2.parent.id[1],
                ),
            )
            for _, links in link_groups:
                link = next(links)
                try:
                    line = self._get_link_line(link)
                except Exception as err:
                    # catch and re-raise with more information
                    raise PDBIOException(f"Error when writing LINK {link}") from err
                else:
                    fhandle.write(line)

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
                        for atom in residue.get_unpacked_list():
                            if select.accept_atom(atom):
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
