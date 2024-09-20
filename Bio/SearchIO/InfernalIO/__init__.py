# Copyright 2024 by Samuel Prince. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.


from .infernal_tab import InfernalTabParser
from .infernal_tab import InfernalTabIndexer
from .infernal_text import InfernalTextParser
from .infernal_text import InfernalTextIndexer


# if not used as a module, run the doctest
if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
