# Copyright 2020 by Tianyi Shi.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""This module contains one public class EbiEna as an API for RSCB PDB
(Protein Data Bank). For more information see _SeqDb.
"""

from ._SeqDb import _SeqDb


class RcsbDb(_SeqDb):
    """API for RCSB PDB (Protein Data Bank).

    RCSB stands for Research Collaboratory for Structural Bioinformatics.
    """

    name = "RCSB"
    base_url = "http://www.rcsb.org"
    entry_url = "http://www.rcsb.org/structure/{}"
