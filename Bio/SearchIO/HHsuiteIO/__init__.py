# Copyright 2019 by Jens Thomas. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.SearchIO support for HHSUITE output formats.

This module adds support for parsing HHSUITE version 2 output.

More information about HHSUITE are available through these links:
- Github repository: https://github.com/soedinglab/hh-suite
- Wiki: https://github.com/soedinglab/hh-suite/wiki

"""

from .hhsuite2_text import Hhsuite2TextParser
