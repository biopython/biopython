# Copyright 2000 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
This module is deprecated; its functions are now available from Bio.SCOP.

Provides code to access SCOP over the WWW.  The main SCOP web page
is available at:
http://scop.mrc-lmb.cam.ac.uk/scop/

Functions:
search       Access the main CGI script.
_open

"""

import warnings
warnings.warn("Bio.WWW.SCOP was deprecated. Its functionality is now available from Bio.SCOP.")


from Bio import SCOP
search = SCOP.search
_open = SCOP._open
