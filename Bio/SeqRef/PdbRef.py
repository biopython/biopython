# Copyright 2020 by Tianyi Shi.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

r"""Working with PDB references.

About: http://www.rcsb.org/pages/about-us/index

Contains one public class GbRef for representing PDB references.
"""

from ._SeqRef import SeqRef, _SeqId

from .SeqDb import NcbiProteinDb, RcsbDb


class _PdbId(_SeqId):
    def __init__(self, id, chain):
        self.id = id
        self.chain = chain

    def __str__(self):
        return self.id + "|" + self.chain if self.chain else self.id


class PdbRef(SeqRef):
    """NCBI GenBank reference."""

    name = "PDB"
    databases = (NcbiProteinDb, RcsbDb)
    # https://www.ncbi.nlm.nih.gov/protein/3LZG_L
    # http://www.rcsb.org/structure/3LZG

    def __init__(self, id, chain=""):
        """Initialize a PdbRef Object."""
        self.id = _PdbId(id, chain)
        self.urls = self.get_urls()
