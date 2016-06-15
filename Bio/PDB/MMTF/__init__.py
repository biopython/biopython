from mmtf import fetch,parse
from Bio.PDB.MMTF.DefaultParser import StructureDecoder


def get_from_decoded(decoder):
    structure_decoder = StructureDecoder()
    decoder.pass_data_on(structure_decoder)
    return structure_decoder.structure_bulder.get_structure()

class MMTFParser():

    def get_structure_from_url(self, pdb_id):
        """Get a structure from a URL - given a PDB id
        :param pdb_id the input PDB id
        :return the structure"""
        decoder = fetch(pdb_id)
        return get_from_decoded(decoder)

    def get_structure(self, file_path):
        """Get a structrue from a file - given a file path
        :param file_path the input file path
        :return the structure"""
        decoder = parse(file_path)
        return get_from_decoded(decoder)