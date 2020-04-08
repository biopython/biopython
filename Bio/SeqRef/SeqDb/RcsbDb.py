# Copyright 2020 by Tianyi Shi.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

from ._SeqDb import _SeqDb


class RcsbDb(_SeqDb):
    name = "RCSB"
    base_url = "http://www.rcsb.org"
    entry_url = "http://www.rcsb.org/structure/{}"
