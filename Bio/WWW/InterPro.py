# Copyright 2001 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""InterPro.py

This module provides code to access resources at InterPro over the WWW.
http://www.ebi.ac.uk/interpro


Functions:
get_interpro_entry

"""
import time
import urllib

from Bio import File


def get_interpro_entry( id ):
    """get specified interpro entry"""
    handle = urllib.urlopen("http://www.ebi.ac.uk/interpro/IEntry?ac=" + id )

    # XXX need to check to see if the entry exists!
    return handle


