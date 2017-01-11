# Copyright 2007 by Tiago Antao.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""SimCoal2 execution module and support functions (DEPRECATED)."""

import os
import sys
import warnings

from Bio import BiopythonDeprecationWarning

warnings.warn("Bio.PopGen.SimCoal has been deprecated, and we intend to "
              " remove it in a future release of Biopython. If you would like"
              " to continue using it, please contact the Biopython developers"
              " via the mailing list.", BiopythonDeprecationWarning)


# This is a workaround to work with the test system
# In any case the problem is with the test system
for instance in sys.path:
    test_path = os.path.join(instance, 'Bio', 'PopGen', 'SimCoal', 'data')
    if os.access(test_path, os.F_OK):
        builtin_tpl_dir = test_path
        break
