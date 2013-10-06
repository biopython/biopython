# Copyright 2010 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Dummy module to mimic Python 3 style builtins (PRIVATE).

We currently have lines like this under Python 2 in order
to use iterator based zip, map and filter:

    from future_builtins import zip

There is no similar option for range yet, other than:

    range = xrange

or:

    from __builtin__ import xrange as range

Under Python 3 this imports need to be removed. Also, deliberate
importing of built in functions like open changes from Python 2:

    from __builtin__ import open

to this under Python 3:

    from builtins import open

Instead, we can do this under either Python 2 or 3:

    from Bio._py3k.builtins import open
    from Bio._py3k.builtins import zip

Once we drop support for Python 2, the whole of Bio._py3k will
go away.
"""
import sys

if sys.version_info[0] >= 3:
    #Code for Python 3
    #Note if this is processed with 2to3 it will break!
    from builtins import open, zip, map, filter, range

    #Lots of our Python 2 code uses isinstance(x, basestring)
    #which after 2to3 becomes isinstance(x, str)
    basestring = str
else:
    #Code for Python 2
    from __builtin__ import open, basestring

    #Import Python3 like iterator functions:
    from future_builtins import zip, map, filter
    from __builtin__ import xrange as range
