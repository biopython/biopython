# Copyright 2000 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Collection of modules for dealing with biological data in Python.

The Biopython Project is an international association of developers 
of freely available Python tools for computational molecular biology.

http://biopython.org
"""

__version__ = "1.55+"

class MissingExternalDependencyError(Exception):
    pass

class BiopythonDeprecationWarning(UserWarning):
    """Code marked as deprecated is likely to be removed in a future version
of Biopython. To avoid removal of this code, please contact the Biopython
developers by sending an email to biopython-dev@biopython.org.
"""
    # Use this warning instead of Python's DeprecationWarning, since
    # those are ignored by default since Python 2.7.
