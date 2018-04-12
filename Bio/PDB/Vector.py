# Copyright 2017 Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Deprecated in favor of Bio.PDB.vectors to solve name collision.

Historically the meaning of ``Bio.PDB.Vector`` has been ambiguous,
both a class ``Vector`` and the module which defines this class
and related functions (``Bio/PDB/Vector.py``).

The module has been renamed to ``Bio.PDB.vectors`` (lower case
in line with PEP8 guidelines).

Please use the following style imports in order to work on both
old and new versions of Biopython:

>>> from Bio.PDB import calc_angle
>>> from Bio.PDB import calc_dihedral
>>> from Bio.PDB import m2rotaxis
>>> from Bio.PDB import refmat
>>> from Bio.PDB import rotaxis2m
>>> from Bio.PDB import rotmat
>>> from Bio.PDB import vector_to_axis
>>> from Bio.PDB import Vector  # for the class

"""

import warnings

from .vectors import m2rotaxis, vector_to_axis, rotaxis2m
from .vectors import refmat, rotmat, calc_angle, calc_dihedral
from .vectors import Vector  # the class whose name clashed

from Bio import BiopythonDeprecationWarning

warnings.warn("The module Bio.PDB.Vector has been deprecated in "
              "favor of new module Bio.PDB.vectors to solve a "
              "name collision with the class Vector. For the "
              "class Vector, and vector functions like calc_angle, "
              "import from Bio.PDB instead.", BiopythonDeprecationWarning)
