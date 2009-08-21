# Copyright 2003 by Bartek Wilczynski.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
Parser and code for dealing with the standalone version of AlignAce, a motif search program (DEPRECATED).

It also includes a very simple interface to CompareACE tool.
http://atlas.med.harvard.edu/

Now superseded by Bio.Motif

Class mapping:
Bio.AlignAce.AlignAceStandalone -> Bio.Motif.Applications.AlignAceCommandline
Bio.AlignAce.CompareAceStandalone -> Bio.Motif.Applications.CompareAceCommandline
Bio.Alignace.Motif -> Bio.Motif.Motif
Bio.AlignAce.Parser.CompareAceParser -> Bio.Motif.CompareAceParser

Instead of using Bio.AlignAce.Parser.AlignAceParser class, you can use
Bio.Motif.parse(handle,format="AlignAce")

"""

import warnings
warnings.warn('Bio.AlignAce is deprecated. Please use Bio.Motif instead.',
              DeprecationWarning)
