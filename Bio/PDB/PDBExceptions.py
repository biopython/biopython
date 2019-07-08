# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Some Bio.PDB-specific exceptions."""

from Bio import BiopythonWarning


# General error
class PDBException(Exception):
    """Define class PDBException."""

    pass


# The PDB file cannot be unambiguously represented in the SMCRA
# data structure
class PDBConstructionException(Exception):
    """Define class PDBConstructionException."""

    pass


class PDBConstructionWarning(BiopythonWarning):
    """Define class PDBConstructionWarning."""

    pass
