# Copyright 2013 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Dummy module to mimic Python 3 style urllib (PRIVATE).

We currently have lines like this under Python 2,

   from urllib import 
   handle = urlopen(request)

Under Python 3 after 2to3 fixing, this becomes:

   from urllib.request import urlopen
   handle = urlopen(request)

Instead, we can do this under either Python 2 or 3:

   from Bio._py3k.urllib.request import urlopen
   handle = urlopen(request)

Then, when we eventually drop Python 2, remove the Bio._py3k bit.
"""
