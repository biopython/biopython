# Copyright 2003-2008 by Leighton Pritchard.  All rights reserved.
# Revisions copyright 2009 by Peter Cock.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
#
# Contact:       Leighton Pritchard, The James Hutton Institute,
#                Invergowrie, Dundee, Scotland, DD2 5DA, UK
#                Leighton.Pritchard@hutton.ac.uk
# #############################################################################

"""GenomeDiagram module integrated into Biopython."""

# Local imports, to make these classes available directly under the
# Bio.Graphics.GenomeDiagram namespace:

from ._Diagram import Diagram
from ._Track import Track
from ._FeatureSet import FeatureSet
from ._GraphSet import GraphSet
from ._CrossLink import CrossLink
from ._Colors import ColorTranslator
from ._Feature import Feature
from ._Graph import GraphData

__all__ = (
    "Diagram",
    "Track",
    "FeatureSet",
    "Feature",
    "GraphSet",
    "GraphData",
    "CrossLink",
    "ColorTranslator",
)
