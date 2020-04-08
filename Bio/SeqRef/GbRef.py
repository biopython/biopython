# Copyright 2020 by Tianyi Shi.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

from ._commons import SeqRef, SeqId
from .SeqDb import NcbiNucleotideDb, EbiDB


class GbId(SeqId):
    def __init__(self, id, version):
        if not version:
            split = id.split(".")
            if len(split) == 2:
                id, version = split
        self.id = id
        self.version = str(version)


class GbRef(SeqRef):
    databases = (NcbiNucleotideDb, EbiDB)

    def __init__(self, id, version=""):
        self.id = GbId(id, version)
        self.urls = self.get_urls()

    # https://www.ncbi.nlm.nih.gov/nuccore/CY073775.2/
