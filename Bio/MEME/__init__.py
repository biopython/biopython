# Copyright 2005 by Jason A. Hackney.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
Parser for dealing with text output from the MEME motif search program (DEPRECATED).

Now superseded by Bio.Motif

Class mapping:
Bio.MEME.Motif.MEMEMotif -> Bio.Motif.Motif
Bio.MEME.parser.MEMEParser -> Bio.Motif.MEMEParser
Bio.MEME.parser.MASTParser -> Bio.Motif.MASTParser

Instead of using Bio.MEME.parser.MEMEParser class, you can use
Bio.Motif.parse(handle,format="MEME")

"""

import warnings
warnings.warn('Bio.MEME is deprecated. Please use Bio.Motif instead.',
              DeprecationWarning)
