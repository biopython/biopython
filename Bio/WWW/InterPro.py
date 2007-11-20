# Copyright 2001 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


"""
This module is deprecated; its functions are now available from Bio.InterPro.

This module provides code to access resources at InterPro over the WWW.
http://www.ebi.ac.uk/interpro


Functions:
get_interpro_entry

"""


import warnings
warnings.warn("Bio.WWW.InterPro was deprecated. Its functionality is now available from Bio.InterPro.")


from Bio import InterPro
get_interpro_entry = InterPro.get_interpro_entry
