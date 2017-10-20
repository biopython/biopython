"""Support for loading 3D structures stored in MMTF files."""
try:
    from mmtf import fetch, parse
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError("Install mmtf to use Bio.PDB.mmtf "
                                       "(e.g. pip install mmtf-python)")
from Bio.PDB.mmtf.DefaultParser import StructureDecoder


def get_from_decoded(decoder):
    structure_decoder = StructureDecoder()
    decoder.pass_data_on(structure_decoder)
    return structure_decoder.structure_bulder.get_structure()


class MMTFParser(object):
    """Class to get a Biopython structure from a URL or a filename."""

    @staticmethod
    def get_structure_from_url(pdb_id):
        """Get a structure from a URL - given a PDB id.

        :param pdb_id: the input PDB id
        :return: the structure

        """
        decoder = fetch(pdb_id)
        return get_from_decoded(decoder)

    @staticmethod
    def get_structure(file_path):
        """Get a structure from a file - given a file path.

        :param file_path: the input file path
        :return: the structure

        """
        decoder = parse(file_path)
        return get_from_decoded(decoder)
