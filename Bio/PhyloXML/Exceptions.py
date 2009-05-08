# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Exception types specific to the PhyloXML module."""

class PhyloXMLError(Exception):
    """Exception raised when PhyloXML object construction cannot continue."""
    pass

class PhyloXMLWarning(Warning):
    """Warning for non-fatal parsing or construction issues."""
    pass
