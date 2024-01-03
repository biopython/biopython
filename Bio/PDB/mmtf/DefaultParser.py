# Copyright 2016 Anthony Bradley.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Code handle loading mmtf-python into Biopython's structures."""

from Bio.PDB.StructureBuilder import StructureBuilder
import numpy as np


class StructureDecoder:
    """Class to pass the data from mmtf-python into a Biopython data structure."""

    def __init__(self):
        """Initialize the class."""
        self.this_type = ""

    def init_structure(
        self,
        total_num_bonds,
        total_num_atoms,
        total_num_groups,
        total_num_chains,
        total_num_models,
        structure_id,
    ):
        """Initialize the structure object.

        :param total_num_bonds: the number of bonds in the structure
        :param total_num_atoms: the number of atoms in the structure
        :param total_num_groups: the number of groups in the structure
        :param total_num_chains: the number of chains in the structure
        :param total_num_models: the number of models in the structure
        :param structure_id: the id of the structure (e.g. PDB id)

        """
        self.structure_builder = StructureBuilder()
        self.structure_builder.init_structure(structure_id=structure_id)
        self.chain_index_to_type_map = {}
        self.chain_index_to_seq_map = {}
        self.chain_index_to_description_map = {}
        self.chain_counter = 0

    def set_atom_info(
        self,
        atom_name,
        serial_number,
        alternative_location_id,
        x,
        y,
        z,
        occupancy,
        temperature_factor,
        element,
        charge,
    ):
        """Create an atom object an set the information.

        :param atom_name: the atom name, e.g. CA for this atom
        :param serial_number: the serial id of the atom (e.g. 1)
        :param alternative_location_id: the alternative location id for the atom, if present
        :param x: the x coordinate of the atom
        :param y: the y coordinate of the atom
        :param z: the z coordinate of the atom
        :param occupancy: the occupancy of the atom
        :param temperature_factor: the temperature factor of the atom
        :param element: the element of the atom, e.g. C for carbon. According to IUPAC. Calcium  is Ca
        :param charge: the formal atomic charge of the atom

        """
        # MMTF uses "\x00" (the NUL character) to indicate to altloc, so convert
        # that to the space required by StructureBuilder
        if alternative_location_id == "\x00":
            alternative_location_id = " "

        # Atom_name is in twice - the full_name is with spaces
        self.structure_builder.init_atom(
            str(atom_name),
            np.array((x, y, z), "f"),
            temperature_factor,
            occupancy,
            alternative_location_id,
            str(atom_name),
            serial_number=serial_number,
            element=str(element).upper(),
        )

    def set_chain_info(self, chain_id, chain_name, num_groups):
        """Set the chain information.

        :param chain_id: the asym chain id from mmCIF
        :param chain_name: the auth chain id from mmCIF
        :param num_groups: the number of groups this chain has

        """
        # A Bradley - chose to use chain_name (auth_id) as it complies
        # with current Biopython. Chain_id might be better.
        self.structure_builder.init_chain(chain_id=chain_name)
        if self.chain_index_to_type_map[self.chain_counter] == "polymer":
            self.this_type = " "
        elif self.chain_index_to_type_map[self.chain_counter] == "non-polymer":
            self.this_type = "H"
        elif self.chain_index_to_type_map[self.chain_counter] == "water":
            self.this_type = "W"
        self.chain_counter += 1

    def set_entity_info(self, chain_indices, sequence, description, entity_type):
        """Set the entity level information for the structure.

        :param chain_indices: the indices of the chains for this entity
        :param sequence: the one letter code sequence for this entity
        :param description: the description for this entity
        :param entity_type: the entity type (polymer,non-polymer,water)

        """
        for chain_ind in chain_indices:
            self.chain_index_to_type_map[chain_ind] = entity_type
            self.chain_index_to_seq_map[chain_ind] = sequence
            self.chain_index_to_description_map[chain_ind] = description

    def set_group_info(
        self,
        group_name,
        group_number,
        insertion_code,
        group_type,
        atom_count,
        bond_count,
        single_letter_code,
        sequence_index,
        secondary_structure_type,
    ):
        """Set the information for a group.

        :param group_name: the name of this group, e.g. LYS
        :param group_number: the residue number of this group
        :param insertion_code: the insertion code for this group
        :param group_type: a string indicating the type of group (as found in the chemcomp dictionary.
            Empty string if none available.
        :param atom_count: the number of atoms in the group
        :param bond_count: the number of unique bonds in the group
        :param single_letter_code: the single letter code of the group
        :param sequence_index: the index of this group in the sequence defined by the entity
        :param secondary_structure_type: the type of secondary structure used
            (types are according to DSSP and number to type mappings are defined in the specification)

        """
        # MMTF uses a NUL character to indicate a blank insertion code, but
        # StructureBuilder expects a space instead.
        if insertion_code == "\x00":
            insertion_code = " "

        self.structure_builder.init_seg(" ")
        self.structure_builder.init_residue(
            group_name, self.this_type, group_number, insertion_code
        )

    def set_model_info(self, model_id, chain_count):
        """Set the information for a model.

        :param model_id: the index for the model
        :param chain_count: the number of chains in the model

        """
        self.structure_builder.init_model(model_id)

    def set_xtal_info(self, space_group, unit_cell):
        """Set the crystallographic information for the structure.

        :param space_group: the space group name, e.g. "P 21 21 21"
        :param unit_cell: an array of length 6 with the unit cell parameters in order: a, b, c, alpha, beta, gamma

        """
        self.structure_builder.set_symmetry(space_group, unit_cell)

    def set_header_info(
        self,
        r_free,
        r_work,
        resolution,
        title,
        deposition_date,
        release_date,
        experimnetal_methods,
    ):
        """Set the header information.

        :param r_free: the measured R-Free for the structure
        :param r_work: the measure R-Work for the structure
        :param resolution: the resolution of the structure
        :param title: the title of the structure
        :param deposition_date: the deposition date of the structure
        :param release_date: the release date of the structure
        :param experimnetal_methods: the list of experimental methods in the structure

        """

    def set_bio_assembly_trans(
        self, bio_assembly_index, input_chain_indices, input_transform
    ):
        """Set the Bioassembly transformation information. A single bioassembly can have multiple transforms.

        :param bio_assembly_index: the integer index of the bioassembly
        :param input_chain_indices: the list of integer indices for the chains of this bioassembly
        :param input_transform: the list of doubles for  the transform of this bioassmbly transform.

        """

    def finalize_structure(self):
        """Any functions needed to cleanup the structure."""

    def set_group_bond(self, atom_index_one, atom_index_two, bond_order):
        """Add bonds within a group.

        :param atom_index_one: the integer atom index (in the group) of the first partner in the bond
        :param atom_index_two: the integer atom index (in the group) of the second partner in the bond
        :param bond_order: the integer bond order

        """

    def set_inter_group_bond(self, atom_index_one, atom_index_two, bond_order):
        """Add bonds between groups.

        :param atom_index_one: the integer atom index (in the structure) of the first partner in the bond
        :param atom_index_two: the integer atom index (in the structure) of the second partner in the bond
        :param bond_order: the bond order

        """
