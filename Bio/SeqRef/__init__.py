# Copyright 2020 by Tianyi Shi.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Object representation of sequence cross-references.

This subpackage has several classes in the form `xxxRef` (e.g. `PdbRef`,
`GbRef`) which represents different types of sequence references. They
are inherited from the _SeqDb class. These `xxxRef` objects are are meant
to be created by different parsers for files that may contain sequence
cross-references, and they are particularly useful in showing URLs to
the entries of a sequence in different databases.
"""

from .GbRef import GbRef, EmblRef, DdbjRef
from .PdbRef import PdbRef
