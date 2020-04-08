# Copyright 2020 by Tianyi Shi.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

from ._commons import SeqRef, SeqId

from .SeqDb import NcbiProteinDb, RcsbDb


class PdbId(SeqId):
    def __init__(self, id, chain):
        self.id = id
        self.chain = chain


class PdbRef(SeqRef):
    databases = (NcbiProteinDb, RcsbDb)
    # https://www.ncbi.nlm.nih.gov/protein/3LZG_L
    # http://www.rcsb.org/structure/3LZG

    def __init__(self, id, chain=""):
        self.id = PdbId(id, chain)
        self.urls = self.get_urls()

    def get_urls(self):
        return {db.name: db.make_entry_url(self.id) for db in self.databases}
