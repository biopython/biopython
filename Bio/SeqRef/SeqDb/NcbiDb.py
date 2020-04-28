# Copyright 2020 by Tianyi Shi.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""This module contains public classes as APIs for several databases hosted
by NCBI. For more information see individual classes, as well as their parent,
_SeqDb.
"""

from ._SeqDb import _SeqDb
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
