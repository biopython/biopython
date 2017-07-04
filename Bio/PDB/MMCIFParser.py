# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""mmCIF parsers."""

from __future__ import print_function

import numpy
import warnings

from Bio.File import as_handle
from Bio._py3k import range

from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.PDB.PDBExceptions import PDBConstructionException
from Bio.PDB.PDBExceptions import PDBConstructionWarning


class MMCIFParser(object):
    """Parse a mmCIF file and return a Structure object."""

    def __init__(self, structure_builder=None, QUIET=False):
        """Create a PDBParser object.

        The mmCIF parser calls a number of standard methods in an aggregated
        StructureBuilder object. Normally this object is instanciated by the
        MMCIParser object itself, but if the user provides his/her own
        StructureBuilder object, the latter is used instead.

        Arguments:
         - structure_builder - an optional user implemented StructureBuilder class.
         - QUIET - Evaluated as a Boolean. If true, warnings issued in constructing
           the SMCRA data will be suppressed. If false (DEFAULT), they will be shown.
           These warnings might be indicative of problems in the mmCIF file!

        """
        if structure_builder is not None:
            self._structure_builder = structure_builder
        else:
            self._structure_builder = StructureBuilder()
        # self.header = None
        # self.trailer = None
        self.line_counter = 0
        self.build_structure = None
        self.QUIET = bool(QUIET)

    # Public methods

    def get_structure(self, structure_id, filename):
        """Return the structure.

        Arguments:
         - structure_id - string, the id that will be used for the structure
         - filename - name of the mmCIF file OR an open filehandle

        """
        with warnings.catch_warnings():
            if self.QUIET:
                warnings.filterwarnings("ignore", category=PDBConstructionWarning)
            self._mmcif_dict = MMCIF2Dict(filename)
            self._build_structure(structure_id)

        return self._structure_builder.get_structure()

    # Private methods

    def _build_structure(self, structure_id):
        mmcif_dict = self._mmcif_dict
        atom_id_list = mmcif_dict["_atom_site.label_atom_id"]
        residue_id_list = mmcif_dict["_atom_site.label_comp_id"]
        try:
            element_list = mmcif_dict["_atom_site.type_symbol"]
        except KeyError:
            element_list = None
        seq_id_list = mmcif_dict["_atom_site.label_seq_id"]
        chain_id_list = mmcif_dict["_atom_site.auth_asym_id"]
        x_list = [float(x) for x in mmcif_dict["_atom_site.Cartn_x"]]
        y_list = [float(x) for x in mmcif_dict["_atom_site.Cartn_y"]]
        z_list = [float(x) for x in mmcif_dict["_atom_site.Cartn_z"]]
        alt_list = mmcif_dict["_atom_site.label_alt_id"]
        icode_list = mmcif_dict["_atom_site.pdbx_PDB_ins_code"]
        b_factor_list = mmcif_dict["_atom_site.B_iso_or_equiv"]
        occupancy_list = mmcif_dict["_atom_site.occupancy"]
        fieldname_list = mmcif_dict["_atom_site.group_PDB"]
        try:
            serial_list = [int(n) for n in mmcif_dict["_atom_site.pdbx_PDB_model_num"]]
        except KeyError:
            # No model number column
            serial_list = None
        except ValueError:
            # Invalid model number (malformed file)
            raise PDBConstructionException("Invalid model number")
        try:
            aniso_u11 = mmcif_dict["_atom_site.aniso_U[1][1]"]
            aniso_u12 = mmcif_dict["_atom_site.aniso_U[1][2]"]
            aniso_u13 = mmcif_dict["_atom_site.aniso_U[1][3]"]
            aniso_u22 = mmcif_dict["_atom_site.aniso_U[2][2]"]
            aniso_u23 = mmcif_dict["_atom_site.aniso_U[2][3]"]
            aniso_u33 = mmcif_dict["_atom_site.aniso_U[3][3]"]
            aniso_flag = 1
        except KeyError:
            # no anisotropic B factors
            aniso_flag = 0
        # if auth_seq_id is present, we use this.
        # Otherwise label_seq_id is used.
        if "_atom_site.auth_seq_id" in mmcif_dict:
            seq_id_list = mmcif_dict["_atom_site.auth_seq_id"]
        else:
            seq_id_list = mmcif_dict["_atom_site.label_seq_id"]
        # Now loop over atoms and build the structure
        current_chain_id = None
        current_residue_id = None
        current_resname = None
        structure_builder = self._structure_builder
        structure_builder.init_structure(structure_id)
        structure_builder.init_seg(" ")
        # Historically, Biopython PDB parser uses model_id to mean array index
        # so serial_id means the Model ID specified in the file
        current_model_id = -1
        current_serial_id = -1
        for i in range(0, len(atom_id_list)):

            # set the line_counter for 'ATOM' lines only and not
            # as a global line counter found in the PDBParser()
            # this number should match the '_atom_site.id' index in the MMCIF
            structure_builder.set_line_counter(i)

            x = x_list[i]
            y = y_list[i]
            z = z_list[i]
            resname = residue_id_list[i]
            chainid = chain_id_list[i]
            altloc = alt_list[i]
            if altloc == ".":
                altloc = " "
            int_resseq = int(seq_id_list[i])
            icode = icode_list[i]
            if icode == "?":
                icode = " "
            name = atom_id_list[i]
            # occupancy & B factor
            try:
                tempfactor = float(b_factor_list[i])
            except ValueError:
                raise PDBConstructionException("Invalid or missing B factor")
            try:
                occupancy = float(occupancy_list[i])
            except ValueError:
                raise PDBConstructionException("Invalid or missing occupancy")
            fieldname = fieldname_list[i]
            if fieldname == "HETATM":
                if resname == "HOH" or resname == "WAT":
                    hetatm_flag = "W"
                else:
                    hetatm_flag = "H"
            else:
                hetatm_flag = " "

            resseq = (hetatm_flag, int_resseq, icode)

            if serial_list is not None:
                # model column exists; use it
                serial_id = serial_list[i]
                if current_serial_id != serial_id:
                    # if serial changes, update it and start new model
                    current_serial_id = serial_id
                    current_model_id += 1
                    structure_builder.init_model(current_model_id, current_serial_id)
                    current_chain_id = None
                    current_residue_id = None
                    current_resname = None
            else:
                # no explicit model column; initialize single model
                structure_builder.init_model(current_model_id)

            if current_chain_id != chainid:
                current_chain_id = chainid
                structure_builder.init_chain(current_chain_id)
                current_residue_id = None
                current_resname = None

            if current_residue_id != resseq or current_resname != resname:
                current_residue_id = resseq
                current_resname = resname
                structure_builder.init_residue(resname, hetatm_flag, int_resseq, icode)

            coord = numpy.array((x, y, z), 'f')
            element = element_list[i] if element_list else None
            structure_builder.init_atom(name, coord, tempfactor, occupancy, altloc,
                name, element=element)
            if aniso_flag == 1:
                u = (aniso_u11[i], aniso_u12[i], aniso_u13[i],
                    aniso_u22[i], aniso_u23[i], aniso_u33[i])
                mapped_anisou = [float(x) for x in u]
                anisou_array = numpy.array(mapped_anisou, 'f')
                structure_builder.set_anisou(anisou_array)
        # Now try to set the cell
        try:
            a = float(mmcif_dict["_cell.length_a"])
            b = float(mmcif_dict["_cell.length_b"])
            c = float(mmcif_dict["_cell.length_c"])
            alpha = float(mmcif_dict["_cell.angle_alpha"])
            beta = float(mmcif_dict["_cell.angle_beta"])
            gamma = float(mmcif_dict["_cell.angle_gamma"])
            cell = numpy.array((a, b, c, alpha, beta, gamma), 'f')
            spacegroup = mmcif_dict["_symmetry.space_group_name_H-M"]
            spacegroup = spacegroup[1:-1]  # get rid of quotes!!
            if spacegroup is None:
                raise Exception
            structure_builder.set_symmetry(spacegroup, cell)
        except Exception:
            pass    # no cell found, so just ignore


class FastMMCIFParser(object):
    """Parse an MMCIF file and return a Structure object."""

    def __init__(self, structure_builder=None, QUIET=False):
        """Create a FastMMCIFParser object.

        The mmCIF parser calls a number of standard methods in an aggregated
        StructureBuilder object. Normally this object is instanciated by the
        parser object itself, but if the user provides his/her own
        StructureBuilder object, the latter is used instead.

        The main difference between this class and the regular MMCIFParser is
        that only 'ATOM' and 'HETATM' lines are parsed here. Use if you are
        interested only in coordinate information.

        Arguments:
         - structure_builder - an optional user implemented StructureBuilder class.
         - QUIET - Evaluated as a Boolean. If true, warnings issued in constructing
           the SMCRA data will be suppressed. If false (DEFAULT), they will be shown.
           These warnings might be indicative of problems in the mmCIF file!

        """
        if structure_builder is not None:
            self._structure_builder = structure_builder
        else:
            self._structure_builder = StructureBuilder()

        self.line_counter = 0
        self.build_structure = None
        self.QUIET = bool(QUIET)

    # Public methods

    def get_structure(self, structure_id, filename):
        """Return the structure.

        Arguments:
         - structure_id - string, the id that will be used for the structure
         - filename - name of the mmCIF file OR an open filehandle

        """
        with warnings.catch_warnings():
            if self.QUIET:
                warnings.filterwarnings("ignore", category=PDBConstructionWarning)
            with as_handle(filename) as handle:
                self._build_structure(structure_id, handle)

        return self._structure_builder.get_structure()

    # Private methods

    def _build_structure(self, structure_id, filehandle):

        # Read only _atom_site. and atom_site_anisotrop entries
        read_atom, read_aniso = False, False
        _fields, _records = [], []
        _anisof, _anisors = [], []
        for line in filehandle:
            if line.startswith('_atom_site.'):
                read_atom = True
                _fields.append(line.strip())
            elif line.startswith('_atom_site_anisotrop.'):
                read_aniso = True
                _anisof.append(line.strip())
            elif read_atom and line.startswith('#'):
                read_atom = False
            elif read_aniso and line.startswith('#'):
                read_aniso = False
            elif read_atom:
                _records.append(line.strip())
            elif read_aniso:
                _anisors.append(line.strip())

        # Dumping the shlex module here since this particular
        # category should be rather straightforward.
        # Quite a performance boost..
        _record_tbl = zip(*map(str.split, _records))
        _anisob_tbl = zip(*map(str.split, _anisors))

        mmcif_dict = dict(zip(_fields, _record_tbl))
        mmcif_dict.update(dict(zip(_anisof, _anisob_tbl)))

        # Build structure object
        atom_id_list = mmcif_dict["_atom_site.label_atom_id"]
        residue_id_list = mmcif_dict["_atom_site.label_comp_id"]

        try:
            element_list = mmcif_dict["_atom_site.type_symbol"]
        except KeyError:
            element_list = None

        seq_id_list = mmcif_dict["_atom_site.label_seq_id"]
        chain_id_list = mmcif_dict["_atom_site.auth_asym_id"]

        x_list = [float(x) for x in mmcif_dict["_atom_site.Cartn_x"]]
        y_list = [float(x) for x in mmcif_dict["_atom_site.Cartn_y"]]
        z_list = [float(x) for x in mmcif_dict["_atom_site.Cartn_z"]]
        alt_list = mmcif_dict["_atom_site.label_alt_id"]
        icode_list = mmcif_dict["_atom_site.pdbx_PDB_ins_code"]
        b_factor_list = mmcif_dict["_atom_site.B_iso_or_equiv"]
        occupancy_list = mmcif_dict["_atom_site.occupancy"]
        fieldname_list = mmcif_dict["_atom_site.group_PDB"]

        try:
            serial_list = [int(n) for n in mmcif_dict["_atom_site.pdbx_PDB_model_num"]]
        except KeyError:
            # No model number column
            serial_list = None
        except ValueError:
            # Invalid model number (malformed file)
            raise PDBConstructionException("Invalid model number")

        try:
            aniso_u11 = mmcif_dict["_atom_site.aniso_U[1][1]"]
            aniso_u12 = mmcif_dict["_atom_site.aniso_U[1][2]"]
            aniso_u13 = mmcif_dict["_atom_site.aniso_U[1][3]"]
            aniso_u22 = mmcif_dict["_atom_site.aniso_U[2][2]"]
            aniso_u23 = mmcif_dict["_atom_site.aniso_U[2][3]"]
            aniso_u33 = mmcif_dict["_atom_site.aniso_U[3][3]"]
            aniso_flag = 1
        except KeyError:
            # no anisotropic B factors
            aniso_flag = 0

        # if auth_seq_id is present, we use this.
        # Otherwise label_seq_id is used.
        if "_atom_site.auth_seq_id" in mmcif_dict:
            seq_id_list = mmcif_dict["_atom_site.auth_seq_id"]
        else:
            seq_id_list = mmcif_dict["_atom_site.label_seq_id"]

        # Now loop over atoms and build the structure
        current_chain_id = None
        current_residue_id = None
        current_resname = None
        structure_builder = self._structure_builder
        structure_builder.init_structure(structure_id)
        structure_builder.init_seg(" ")

        # Historically, Biopython PDB parser uses model_id to mean array index
        # so serial_id means the Model ID specified in the file
        current_model_id = -1
        current_serial_id = -1
        for i in range(0, len(atom_id_list)):

            # set the line_counter for 'ATOM' lines only and not
            # as a global line counter found in the PDBParser()
            # this number should match the '_atom_site.id' index in the MMCIF
            structure_builder.set_line_counter(i)

            x = x_list[i]
            y = y_list[i]
            z = z_list[i]
            resname = residue_id_list[i]
            chainid = chain_id_list[i]
            altloc = alt_list[i]
            if altloc == ".":
                altloc = " "
            int_resseq = int(seq_id_list[i])
            icode = icode_list[i]
            if icode == "?":
                icode = " "
            name = atom_id_list[i].strip('"')  # Remove occasional " from quoted atom names (e.g. xNA)

            # occupancy & B factor
            try:
                tempfactor = float(b_factor_list[i])
            except ValueError:
                raise PDBConstructionException("Invalid or missing B factor")

            try:
                occupancy = float(occupancy_list[i])
            except ValueError:
                raise PDBConstructionException("Invalid or missing occupancy")

            fieldname = fieldname_list[i]
            if fieldname == "HETATM":
                hetatm_flag = "H"
            else:
                hetatm_flag = " "

            resseq = (hetatm_flag, int_resseq, icode)

            if serial_list is not None:
                # model column exists; use it
                serial_id = serial_list[i]
                if current_serial_id != serial_id:
                    # if serial changes, update it and start new model
                    current_serial_id = serial_id
                    current_model_id += 1
                    structure_builder.init_model(current_model_id, current_serial_id)
                    current_chain_id = None
                    current_residue_id = None
                    current_resname = None
            else:
                # no explicit model column; initialize single model
                structure_builder.init_model(current_model_id)

            if current_chain_id != chainid:
                current_chain_id = chainid
                structure_builder.init_chain(current_chain_id)
                current_residue_id = None
                current_resname = None

            if current_residue_id != resseq or current_resname != resname:
                current_residue_id = resseq
                current_resname = resname
                structure_builder.init_residue(resname, hetatm_flag, int_resseq, icode)

            coord = numpy.array((x, y, z), 'f')
            element = element_list[i] if element_list else None
            structure_builder.init_atom(name, coord, tempfactor, occupancy, altloc,
                name, element=element)
            if aniso_flag == 1:
                u = (aniso_u11[i], aniso_u12[i], aniso_u13[i],
                    aniso_u22[i], aniso_u23[i], aniso_u33[i])
                mapped_anisou = [float(x) for x in u]
                anisou_array = numpy.array(mapped_anisou, 'f')
                structure_builder.set_anisou(anisou_array)
