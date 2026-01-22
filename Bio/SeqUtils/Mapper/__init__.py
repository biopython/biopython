# Copyright 2012 Lenna X. Peterson <arklenna@gmail.com>
# CoordinateMapper.py originally written by Reece Hart
# Older revisions may be found in this gist:
#     https://gist.github.com/3172753

# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Coordinate mapper for transformation of positions.

Transforms between genomic, CDS, and protein coordinates.
Includes methods for converting locations to and from HGVS conventions.
Genbank locations can be parsed with SeqIO.
"""

__all__ = [
    "MapPosition",
    "GenomePosition",
    "CDSPosition",
    "ProteinPosition",
    "CoordinateMapper",
]

import warnings

from Bio import BiopythonExperimentalWarning

from .MapPositions import MapPosition, GenomePosition, CDSPosition, ProteinPosition
from ._CoordinateMapper import CoordinateMapper

warnings.warn(
    "Bio.SeqUtils.Mapper is an experimental module which may undergo "
    "significant changes prior to its future official release.",
    BiopythonExperimentalWarning,
)
