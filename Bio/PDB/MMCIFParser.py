# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""mmCIF parser"""

from __future__ import print_function

from string import ascii_letters

import numpy

from Bio._py3k import range

from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.PDB.PDBExceptions import PDBConstructionException


class MMCIFParser(object):
    def get_structure(self, structure_id, filename):
        self._mmcif_dict=MMCIF2Dict(filename)
        self._structure_builder=StructureBuilder()
        self._build_structure(structure_id)
        return self._structure_builder.get_structure()

    def _build_structure(self, structure_id):
        mmcif_dict=self._mmcif_dict
        atom_id_list=mmcif_dict["_atom_site.label_atom_id"]
        residue_id_list=mmcif_dict["_atom_site.label_comp_id"]
        try:
            element_list = mmcif_dict["_atom_site.type_symbol"]
        except KeyError:
            element_list = None
        seq_id_list=mmcif_dict["_atom_site.label_seq_id"]
        chain_id_list=mmcif_dict["_atom_site.label_asym_id"]
        x_list = [float(x) for x in mmcif_dict["_atom_site.Cartn_x"]]
        y_list = [float(x) for x in mmcif_dict["_atom_site.Cartn_y"]]
        z_list = [float(x) for x in mmcif_dict["_atom_site.Cartn_z"]]
        alt_list=mmcif_dict["_atom_site.label_alt_id"]
        b_factor_list=mmcif_dict["_atom_site.B_iso_or_equiv"]
        occupancy_list=mmcif_dict["_atom_site.occupancy"]
        fieldname_list=mmcif_dict["_atom_site.group_PDB"]
        try:
            serial_list = [int(n) for n in mmcif_dict["_atom_site.pdbx_PDB_model_num"]]
        except KeyError:
            # No model number column
            serial_list = None
        except ValueError:
            # Invalid model number (malformed file)
            raise PDBConstructionException("Invalid model number")
        try:
            aniso_u11=mmcif_dict["_atom_site.aniso_U[1][1]"]
            aniso_u12=mmcif_dict["_atom_site.aniso_U[1][2]"]
            aniso_u13=mmcif_dict["_atom_site.aniso_U[1][3]"]
            aniso_u22=mmcif_dict["_atom_site.aniso_U[2][2]"]
            aniso_u23=mmcif_dict["_atom_site.aniso_U[2][3]"]
            aniso_u33=mmcif_dict["_atom_site.aniso_U[3][3]"]
            aniso_flag=1
        except KeyError:
            # no anisotropic B factors
            aniso_flag=0
        # if auth_seq_id is present, we use this.
        # Otherwise label_seq_id is used.
        if "_atom_site.auth_seq_id" in mmcif_dict:
            seq_id_list=mmcif_dict["_atom_site.auth_seq_id"]
        else:
            seq_id_list=mmcif_dict["_atom_site.label_seq_id"]
        # Now loop over atoms and build the structure
        current_chain_id=None
        current_residue_id=None
        structure_builder=self._structure_builder
        structure_builder.init_structure(structure_id)
        structure_builder.init_seg(" ")
        # Historically, Biopython PDB parser uses model_id to mean array index
        # so serial_id means the Model ID specified in the file
        current_model_id = 0
        current_serial_id = 0
        for i in range(0, len(atom_id_list)):
            x=x_list[i]
            y=y_list[i]
            z=z_list[i]
            resname=residue_id_list[i]
            chainid=chain_id_list[i]
            altloc=alt_list[i]
            if altloc==".":
                altloc=" "
            resseq=seq_id_list[i]
            name=atom_id_list[i]
            # occupancy & B factor
            try:
                tempfactor=float(b_factor_list[i])
            except ValueError:
                raise PDBConstructionException("Invalid or missing B factor")
            try:
                occupancy=float(occupancy_list[i])
            except ValueError:
                raise PDBConstructionException("Invalid or missing occupancy")
            fieldname=fieldname_list[i]
            if fieldname=="HETATM":
                hetatm_flag="H"
            else:
                hetatm_flag=" "
            if serial_list is not None:
                # model column exists; use it
                serial_id = serial_list[i]
                if current_serial_id != serial_id:
                    # if serial changes, update it and start new model
                    current_serial_id = serial_id
                    structure_builder.init_model(current_model_id, current_serial_id)
                    current_model_id += 1
            else:
                # no explicit model column; initialize single model
                structure_builder.init_model(current_model_id)
            if current_chain_id!=chainid:
                current_chain_id=chainid
                structure_builder.init_chain(current_chain_id)
                current_residue_id=resseq
                icode, int_resseq=self._get_icode(resseq)
                structure_builder.init_residue(resname, hetatm_flag, int_resseq,
                    icode)
            elif current_residue_id!=resseq:
                current_residue_id=resseq
                icode, int_resseq=self._get_icode(resseq)
                structure_builder.init_residue(resname, hetatm_flag, int_resseq,
                    icode)
            coord=numpy.array((x, y, z), 'f')
            element = element_list[i] if element_list else None
            structure_builder.init_atom(name, coord, tempfactor, occupancy, altloc,
                name, element=element)
            if aniso_flag==1:
                u=(aniso_u11[i], aniso_u12[i], aniso_u13[i],
                    aniso_u22[i], aniso_u23[i], aniso_u33[i])
                mapped_anisou = [float(x) for x in u]
                anisou_array=numpy.array(mapped_anisou, 'f')
                structure_builder.set_anisou(anisou_array)
        # Now try to set the cell
        try:
            a=float(mmcif_dict["_cell.length_a"])
            b=float(mmcif_dict["_cell.length_b"])
            c=float(mmcif_dict["_cell.length_c"])
            alpha=float(mmcif_dict["_cell.angle_alpha"])
            beta=float(mmcif_dict["_cell.angle_beta"])
            gamma=float(mmcif_dict["_cell.angle_gamma"])
            cell=numpy.array((a, b, c, alpha, beta, gamma), 'f')
            spacegroup=mmcif_dict["_symmetry.space_group_name_H-M"]
            spacegroup=spacegroup[1:-1] # get rid of quotes!!
            if spacegroup is None:
                raise Exception
            structure_builder.set_symmetry(spacegroup, cell)
        except:
            pass    # no cell found, so just ignore

    def _get_icode(self, resseq):
        """Tries to return the icode. In MMCIF files this is just part of
        resseq! In PDB files, it's a separate field."""
        last_resseq_char=resseq[-1]
        if last_resseq_char in ascii_letters:
            icode=last_resseq_char
            int_resseq=int(resseq[0:-1])
        else:
            icode=" "
            int_resseq=int(resseq)
        return icode, int_resseq


if __name__=="__main__":
    import sys

    if len(sys.argv) != 2:
        print("Usage: python MMCIFparser.py filename")
        raise SystemExit
    filename=sys.argv[1]

    p=MMCIFParser()

    structure=p.get_structure("test", filename)

    for model in structure.get_list():
        print(model)
        for chain in model.get_list():
            print(chain)
            print("Found %d residues." % len(chain.get_list()))
