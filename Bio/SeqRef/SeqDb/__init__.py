# Copyright 2020 by Tianyi Shi.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

r"""Classes of sequence databases.

These classes know how to construct urls from given accession codes.
"""

from .NcbiDb import NcbiNucleotideDb, NcbiProteinDb
from .EbiDb import EbiEnaDB
from .RcsbDb import RcsbDb
from .DdbjDb import DdbjDb

reftype_to_default_db = {
    "gb": NcbiNucleotideDb,
    "ebi": EbiEnaDB,
    "rcsb": RcsbDb,
    "ddbj": DdbjDb,
    "dbj": DdbjDb,
}


def fetch_seq(accession_code, reftype, file_format=None, database=None):
    db = database if database else reftype_to_default_db[reftype]
    return db.fetch(accession_code, file_format)
