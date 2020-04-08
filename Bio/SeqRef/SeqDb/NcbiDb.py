from ._commons import _SeqDb


# https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=631490879&conwithfeat=on&withparts=on&hide-cdd=on
# https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=631490879&db=nuccore&report=genbank&conwithfeat=on&withparts=on&hide-cdd=on&retmode=html&withmarkup=on&tool=portal&log$=seqview&maxdownloadsize=1000000
# https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=gbwithparts&id=631490879&withparts=on
#
from copy import deepcopy


class NcbiNucleotideDb(_SeqDb):
    # https://www.ncbi.nlm.nih.gov/nuccore/CY009444.1/
    name = "NCBI Nucleotide"
    base_url = "https://www.ncbi.nlm.nih.gov/nuccore"
    fetch_url = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?save=file&log$=seqview&db=nuccore&report={format}&id={id}&withparts=on"
    fetch_file_format_map = deepcopy(_SeqDb.fetch_file_format_map)
    fetch_file_format_map.update({"gb": "gbwithparts"})

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
