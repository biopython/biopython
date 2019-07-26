# Copyright 2019 Joe Greener. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Write a MMTF file."""

from Bio._py3k import basestring
from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.PDB.PDBIO import Select
from mmtf.api.mmtf_writer import MMTFEncoder

_select = Select()

class MMTFIO(object):
    """Write a Structure object as a MMTF file.

    Examples
    --------
        >>> from Bio.PDB import MMCIFParser
        >>> from Bio.PDB.mmtf import MMTFIO
        >>> parser = MMCIFParser()
        >>> structure = parser.get_structure("1a8o", "PDB/1A8O.cif")
        >>> io=MMTFIO()
        >>> io.set_structure(structure)
        >>> io.save("bio-pdb-mmtf-out.mmtf")
        >>> import os
        >>> os.remove("bio-pdb-mmtf-out.mmtf")  # tidy up

    """

    def __init__(self):
        """Initialise."""
        pass

    def set_structure(self, pdb_object):
        """Check what object the user is providing and build a structure."""
        # This is duplicated from the PDBIO class
        if pdb_object.level == "S":
            structure = pdb_object
        else:
            sb = StructureBuilder()
            sb.init_structure('pdb')
            sb.init_seg(' ')
            # Build parts as necessary
            if pdb_object.level == "M":
                sb.structure.add(pdb_object)
                self.structure = sb.structure
            else:
                sb.init_model(0)
                if pdb_object.level == "C":
                    sb.structure[0].add(pdb_object)
                else:
                    sb.init_chain('A')
                    if pdb_object.level == "R":
                        try:
                            parent_id = pdb_object.parent.id
                            sb.structure[0]['A'].id = parent_id
                        except ValueError:
                            pass
                        sb.structure[0]['A'].add(pdb_object)
                    else:
                        # Atom
                        sb.init_residue('DUM', ' ', 1, ' ')
                        try:
                            parent_id = pdb_object.parent.parent.id
                            sb.structure[0]['A'].id = parent_id
                        except ValueError:
                            pass
                        sb.structure[0]['A'].child_list[0].add(pdb_object)

            # Return structure
            structure = sb.structure
        self.structure = structure

    def save(self, filepath, select=_select):
        """Save the structure to a file.

        :param filepath: output file
        :type filepath: string

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
        """
        # Similar to the PDBIO save method, we check if the filepath is a
        # string for a filepath or an open file handle
        if not isinstance(filepath, basestring):
            raise ValueError("Writing to a file handle is not supported for MMTF, filepath must be a string")
        if hasattr(self, "structure"):
            self._save_structure(filepath, select)
        else:
            raise ValueError("Use set_structure to set a structure to write out")

    def _save_structure(self, filepath, select):
        encoder = MMTFEncoder()
        encoder.init_structure(
            total_num_bonds=0,
            total_num_atoms=len(list(self.structure.get_atoms())),
            total_num_groups=len(list(self.structure.get_residues())),
            total_num_chains=len(list(self.structure.get_chains())),
            total_num_models=len(list(self.structure.get_models())),
            structure_id=self.structure.id
        )

        encoder.set_xtal_info(
            space_group="",
            unit_cell=[0, 0, 0, 0, 0, 0]
        )

        encoder.set_header_info(
            r_free=0.0,
            r_work=0.0,
            resolution=self.structure.header["resolution"],
            title=self.structure.header["name"],
            deposition_date=self.structure.header["deposition_date"],
            release_date=self.structure.header["release_date"],
            experimental_methods=self.structure.header["structure_method"]
        )

        for mi, model in enumerate(self.structure.get_models()):
            if not select.accept_model(model):
                continue

            encoder.set_model_info(
                model_id=mi,
                chain_count=len(list(model.get_chains()))
            )
            for ci, chain in enumerate(model.get_chains()):
                if not select.accept_chain(chain):
                    continue

                encoder.set_chain_info(
                    chain_id=chain.get_id(),
                    chain_name=chain.get_id(),
                    num_groups=len(list(chain.get_residues()))
                )

                prev_residue_type = ""
                prev_resname = ""

                for residue in chain.get_residues():
                    if not select.accept_residue(residue):
                        continue

                    hetfield, resseq, icode = residue.get_id()
                    if hetfield == " ":
                        residue_type = "ATOM"
                        entity_type = "polymer"
                    elif hetfield == "W":
                        residue_type = "HETATM"
                        entity_type = "water"
                    else:
                        residue_type = "HETATM"
                        entity_type = "non-polymer"
                    resname = residue.get_resname()
                    # Check if the molecule changes within the chain
                    # This will always increment for the first residue in a
                    # chain due to the starting values above
                    if residue_type != prev_residue_type or (residue_type == "HETATM" and resname != prev_resname):
                        encoder.set_entity_info(
                            chain_indices=[ci],
                            sequence="",
                            description="",
                            entity_type=entity_type
                        )
                    prev_residue_type = residue_type
                    prev_resname = resname

                    encoder.set_group_info(
                        group_name=resname,
                        group_number=residue.id[1],
                        insertion_code=residue.id[2],
                        group_type="",
                        atom_count=len(list(residue.get_atoms())),
                        bond_count=0,
                        single_letter_code="",
                        sequence_index=0,
                        secondary_structure_type=0
                    )
                    for atom in residue.get_atoms():
                        if select.accept_atom(atom):
                            encoder.set_atom_info(
                                atom_name=atom.name,
                                serial_number=atom.serial_number,
                                alternative_location_id=atom.altloc,
                                x=atom.coord[0],
                                y=atom.coord[1],
                                z=atom.coord[2],
                                occupancy=atom.occupancy,
                                temperature_factor=atom.bfactor,
                                element=atom.element,
                                charge=0
                            )

        encoder.finalize_structure()
        encoder.write_file(filepath)
