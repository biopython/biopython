# Copyright 2020 by Tianyi Shi.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""This module contains public classes as APIs for several Uniprot databases.
For more information see individual classes, as well as their parent, _SeqDb.
"""

from ._SeqDb import _SeqDb


class UniprotDb(_SeqDb):
    pass
