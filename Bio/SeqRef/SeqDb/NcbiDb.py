from ._commons import _SeqDb


class NcbiNucleotideDb(_SeqDb):
    # https://www.ncbi.nlm.nih.gov/nuccore/CY009444.1/
    name = "NCBI Nucleotide"
    base_url = "https://www.ncbi.nlm.nih.gov/nuccore"

    @classmethod
    def make_identifier(cls, id_obj):
        identifier = id_obj.id
        if id_obj.version:
            identifier += "." + id_obj.version
        return identifier


class NcbiProteinDb(_SeqDb):
    # https://www.ncbi.nlm.nih.gov/protein/3LZG_L
    name = "NCBI Protein"
    base_url = "https://www.ncbi.nlm.nih.gov/protein"

    @classmethod
    def make_identifier(cls, id_obj):
        identifier = id_obj.id
        if id_obj.chain:
            identifier += "_" + id_obj.chain
        return identifier
