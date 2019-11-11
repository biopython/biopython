# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Parser for PQR files."""

from __future__ import print_function

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        "Install NumPy if you want to use the PQR parser.")

from Bio.PDB.PDBExceptions import PDBConstructionException
from Bio.PDB.PDBParser import PDBParser


class PQRParser(PDBParser):
    """Create a PQRParser object."""

    def _parse_coordinates(self, coords_trailer):
        """Parse the atomic data in the PQR file (PRIVATE)."""
        local_line_counter = 0
        structure_builder = self.structure_builder
        current_model_id = 0
        # Flag we have an open model
        model_open = 0
        current_chain_id = None
        current_segid = None
        current_residue_id = None
        current_resname = None
        for i in range(0, len(coords_trailer)):
            line = coords_trailer[i].rstrip("\n")
            record_type = line[0:6]
            global_line_counter = self.line_counter + local_line_counter + 1
            structure_builder.set_line_counter(global_line_counter)
            if record_type == "ATOM  " or record_type == "HETATM":
                # Initialize the Model - there was no explicit MODEL record
                if not model_open:
                    structure_builder.init_model(current_model_id)
                    current_model_id += 1
                    model_open = 1
                fullname = line[12:16]
                # get rid of whitespace in atom names
                split_list = fullname.split()
                if len(split_list) != 1:
                    # atom name has internal spaces, e.g. " N B ", so
                    # we do not strip spaces
                    name = fullname
                else:
                    # atom name is like " CA ", so we can strip spaces
                    name = split_list[0]
                altloc = line[16]
                resname = line[17:20]
                chainid = line[21]
                try:
                    serial_number = int(line[6:11])
                except Exception:
                    serial_number = 0
                resseq = int(line[22:26].split()[0])  # sequence identifier
                icode = line[26]  # insertion code
                if record_type == "HETATM":  # hetero atom flag
                    if resname == "HOH" or resname == "WAT":
                        hetero_flag = "W"
                    else:
                        hetero_flag = "H"
                else:
                    hetero_flag = " "
                residue_id = (hetero_flag, resseq, icode)
                # atomic coordinates
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except Exception:
                    # Should we allow parsing to continue in permissive mode?
                    # If so, what coordinates should we default to?  Easier to abort!
                    raise PDBConstructionException("Invalid or missing coordinate(s) at line %i."
                                                   % global_line_counter)
                coord = numpy.array((x, y, z), "f")

                # charge & radius
                try:
                    charge = float(line[54:62])
                except Exception:
                    self._handle_PDB_exception("Invalid or missing charge",
                                               global_line_counter)
                    charge = None  # Rather than arbitrary zero or one
                try:
                    radius = float(line[62:70])
                except Exception:
                    self._handle_PDB_exception("Invalid or missing radius",
                                               global_line_counter)
                    radius = None
                if radius is not None and radius < 0:
                    # In permissive mode raise fatal exception.
                    message = "Negative atom radius"
                    self._handle_PDB_exception(message, global_line_counter)
                    radius = None
                segid = line[72:76]
                element = line[76:78].strip().upper()
                if current_segid != segid:
                    current_segid = segid
                    structure_builder.init_seg(current_segid)
                if current_chain_id != chainid:
                    current_chain_id = chainid
                    structure_builder.init_chain(current_chain_id)
                    current_residue_id = residue_id
                    current_resname = resname
                    try:
                        structure_builder.init_residue(resname, hetero_flag, resseq, icode)
                    except PDBConstructionException as message:
                        self._handle_PDB_exception(message, global_line_counter)
                elif current_residue_id != residue_id or current_resname != resname:
                    current_residue_id = residue_id
                    current_resname = resname
                    try:
                        structure_builder.init_residue(resname, hetero_flag, resseq, icode)
                    except PDBConstructionException as message:
                        self._handle_PDB_exception(message, global_line_counter)
                # init atom
                try:
                    structure_builder.init_atom_from_pqr(name, coord, charge, radius, altloc,
                                                         fullname, serial_number, element)
                except PDBConstructionException as message:
                    self._handle_PDB_exception(message, global_line_counter)
            elif record_type == "ANISOU":
                anisou = [float(x) for x in (line[28:35], line[35:42], line[43:49],
                                             line[49:56], line[56:63], line[63:70])]
                # U's are scaled by 10^4
                anisou_array = (numpy.array(anisou, "f") / 10000.0).astype("f")
                structure_builder.set_anisou(anisou_array)
            elif record_type == "MODEL ":
                try:
                    serial_num = int(line[10:14])
                except Exception:
                    self._handle_PDB_exception("Invalid or missing model serial number",
                                               global_line_counter)
                    serial_num = 0
                structure_builder.init_model(current_model_id, serial_num)
                current_model_id += 1
                model_open = 1
                current_chain_id = None
                current_residue_id = None
            elif record_type == "END   " or record_type == "CONECT":
                # End of atomic data, return the trailer
                self.line_counter += local_line_counter
                return coords_trailer[local_line_counter:]
            elif record_type == "ENDMDL":
                model_open = 0
                current_chain_id = None
                current_residue_id = None
            elif record_type == "SIGUIJ":
                # standard deviation of anisotropic B factor
                siguij = [float(x) for x in (line[28:35], line[35:42], line[42:49],
                                             line[49:56], line[56:63], line[63:70])]
                # U sigma's are scaled by 10^4
                siguij_array = (numpy.array(siguij, "f") / 10000.0).astype("f")
                structure_builder.set_siguij(siguij_array)
            elif record_type == "SIGATM":
                # standard deviation of atomic positions
                sigatm = [float(x) for x in (line[30:38], line[38:45], line[46:54],
                                             line[54:60], line[60:66])]
                sigatm_array = numpy.array(sigatm, "f")
                structure_builder.set_sigatm(sigatm_array)
            local_line_counter += 1
        # EOF (does not end in END or CONECT)
        self.line_counter = self.line_counter + local_line_counter
        return []
