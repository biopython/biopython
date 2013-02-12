# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Some Bio.PDB-specific exceptions."""

from Bio import BiopythonWarning


# General error
class PDBException(Exception):
    pass


# The PDB file cannot be unambiguously represented in the SMCRA
# data structure
class PDBConstructionException(Exception):
    pass


class PDBConstructionWarning(BiopythonWarning):
    pass
