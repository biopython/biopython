# Copyright 2020 by Tianyi Shi.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

r"""Working with NCBI's GenBank code and their equivalents, EMBL and DDBJ.

About: https://www.ncbi.nlm.nih.gov/genbank/
Examples: "CY073775.2", "CY073775"

Contains one public class GbRef for representing GenBank references.
"""

from ._SeqRef import SeqRef, _SeqId
from .SeqDb import NcbiNucleotideDb, EbiEnaDB, DdbjDb


class _GbId(_SeqId):
    def __init__(self, id, version):
        if not version:
            split = id.split(".")
            if len(split) == 2:
                id, version = split
        self.id = id
        self.version = str(version)

    def __str__(self):
        return self.id + "." + self.version if self.version else self.id


class GbRef(SeqRef):
    """NCBI's GenBank reference.
    """

    name = "GenBank"
    databases = (NcbiNucleotideDb, EbiEnaDB, DdbjDb)

    def __init__(self, id, version=""):
        """Initialize a GbRef object.

        Arguments:
            - id - accession code
            - version - (optional) version

        The version can either be specified explicitly, e.g. GbRef("CY073775", 2)
        or implicitly, e.g. GbRef("CY073775.2"), or it can be empty (then the
        latest version is used).
        """
        self.id = _GbId(id, version)
        self.urls = self.get_urls()


class EmblRef(GbRef):
    """EBI's EMBL reference (equivalent to GenBank).
    """

    pass


class DdbjRef(GbRef):
    """DDBJ's reference (equivalent to GenBank).
    """

    pass
