from Bio.PDB.StructureBuilder import StructureBuilder

class StructureDecoder():

    def __init__(self):
        self.this_type = ""

    def init_structure(self, total_num_bonds, total_num_atoms, total_num_groups, total_num_chains, total_num_models, structure_id):
        self.structure_bulder = StructureBuilder()
        self.structure_bulder.init_structure(structure_id=structure_id)
        self.chain_index_to_type_map = {}
        self.chain_index_to_seq_map = {}
        self.chain_index_to_description_map = {}
        self.chain_counter = 0


    def set_atom_info(self, atom_name, serial_number, alternative_location_id, x, y, z, occupancy, temperature_factor, element, charge):
        # Atom_name is in twice - the full_name is with spaces
        self.structure_bulder.init_atom(str(atom_name), [x,y,z], temperature_factor, occupancy, alternative_location_id, str(atom_name),
                  serial_number=serial_number, element=str(element).upper())

    def set_chain_info(self, chain_id, chain_name, group_count):
        self.structure_bulder.init_chain(chain_id=chain_name)
        if self.chain_index_to_type_map[self.chain_counter] == "polymer":
            self.this_type = ""
        elif self.chain_index_to_type_map[self.chain_counter] == "non-polymer":
            self.this_type = "H"
        elif self.chain_index_to_type_map[self.chain_counter] == "water":
            self.this_type = "W"
        self.chain_counter+=1


    def set_entity_info(self, chain_indices, sequence, description, type_):
        for chain_ind in chain_indices:
            self.chain_index_to_type_map[chain_ind] = type_
            self.chain_index_to_seq_map[chain_ind] = sequence
            self.chain_index_to_description_map[chain_ind] = description

    def set_group_info(self, group_name, group_number, insertion_code, group_type, atom_count, bond_count, single_letter_code, sequence_index, secondary_structure_type):
        self.structure_bulder.init_seg(' ')
        self.structure_bulder.init_residue(group_name, self.this_type, group_number, insertion_code)

    def set_model_info(self, model_id, chain_count):
        self.structure_bulder.init_model(model_id)

    def set_xtal_info(self, space_group, unit_cell):
        self.structure_bulder.set_symmetry(space_group, unit_cell)

    def set_header_info(self, r_free, r_work, resolution, title, deposition_date, release_date, experimnetal_methods):
        pass


    def set_bio_assembly_trans(self, bio_assembly_index, input_chain_indices, input_transform):
        pass


    def finalize_structure(self):
        pass


    def set_group_bond(self, atom_index_one, atom_index_two, bond_order):
        pass
    def set_inter_group_bond(self, atom_index_one, atom_index_two, bond_order):
        pass


