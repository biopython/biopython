# Copyright 2020 by Tianyi Shi.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""This module contains one public class EbiEna as an API for EBI-ENA
(European Nucleotide Archive). For more information see _SeqDb.
"""

from ._SeqDb import _SeqDb


class EbiEnaDB(_SeqDb):
    """API for DDBJ EBI-ENA (European Nucleotide Archive).
    """

    name = "EBI-ENA"
    base_url = "https://www.ebi.ac.uk/ena"
    entry_url = "old: https://www.ebi.ac.uk/ena/data/view/{0} new: https://www.ebi.ac.uk/ena/browser/view/{0}"
