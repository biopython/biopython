# Copyright 2020 by Tianyi Shi.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""This module contains one public class DdbjDb as an API for DDBJ
(DNA Data Bank of Japan). For more information see _SeqDb.
"""

from .NcbiDb import NcbiNucleotideDb


class DdbjDb(NcbiNucleotideDb):
    """API for DDBJ (DNA Data Bank of Japan).
    """

    name = "DDBJ"
    base_url = "http://getentry.ddbj.nig.ac.jp"
    entry_url = "http://getentry.ddbj.nig.ac.jp/getentry/na/{}"
    fetch_url = "http://getentry.ddbj.nig.ac.jp/getentry/na/{id}/?format={format}"
    fetch_file_format_map = {"fasta": "fasta", "ddbj": "flatfile"}
    fetch_file_format_default = "ddbj"
