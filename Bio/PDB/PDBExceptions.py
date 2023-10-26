# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

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


# The SMCRA structure could not be written to file
class PDBIOException(Exception):
    """Define class PDBIOException."""

    pass
