# Copyright 2000 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Collection of modules for dealing with biological data in Python.

The Biopython Project is an international association of developers 
of freely available Python tools for computational molecular biology.

http://biopython.org
"""

__docformat__ = "epytext en" #not just plaintext

__version__ = "1.59"

class MissingExternalDependencyError(Exception):
    """Missing an external dependency.

    Used for things like missing command line tools. Important for our unit
    tests to allow skipping tests with missing external dependencies.
    """
    pass

class MissingPythonDependencyError(MissingExternalDependencyError, ImportError):
    """Missing an external python dependency (subclass of ImportError).

    Used for missing Python modules (rather than just a typical ImportError).
    Important for our unit tests to allow skipping tests with missing external
    python dependencies, while also allowing the exception to be caught as an
    ImportError.
    """
    pass

class BiopythonWarning(Warning):
    """Biopython warning.
    
    Biopython should use this warning (or subclasses of it), making it easy to
    silence all our warning messages should you wish to:

        >>> import warnings
        >>> from Bio import BiopythonWarning
        >>> warnings.simplefilter('ignore', BiopythonWarning)

    Consult the warnings module documentation for more details.
    """
    pass

class BiopythonParserWarning(BiopythonWarning):
    """Biopython parser warning.
    
    Some in-valid data files cannot be parsed and will trigger an exception.
    Where a reasonable interpretation is possible, Biopython will issue this
    warning to indicate a potential problem. To silence these warnings, use:

        >>> import warnings
        >>> from Bio import BiopythonParserWarning
        >>> warnings.simplefilter('ignore', BiopythonParserWarning)

    Consult the warnings module documentation for more details.
    """
    pass

class BiopythonDeprecationWarning(BiopythonWarning):
    """Biopython deprecation warning.
    
    Biopython uses this warning instead of the built in DeprecationWarning
    since those are ignored by default since Python 2.7.

    To silence all our deprecation warning messages, use:

        >>> import warnings
        >>> from Bio import BiopythonDeprecationWarning
        >>> warnings.simplefilter('ignore', BiopythonDeprecationWarning)

    Code marked as deprecated is likely to be removed in a future version
    of Biopython. To avoid removal of this code, please contact the Biopython
    developers by sending an email to biopython-dev@biopython.org.
    """
    pass
