# Copyright 2000 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Dom.py

This module provides code to work with the 'dom' files from SCOP.
http://scop.mrc-lmb.cam.ac.uk/scop/

Classes:
Domain         Holds information about one domain.
Iterator       Iterates over records in a 'dom' file.
DomainParser   Parse a domain record.

"""
import string
from types import *

import Location

class Domain:
    """Holds information for one SCOP domain.

    Members:
    sid          The SCOP ID of the entry, e.g. d1anu1
    pdbid        The PDB ID for the sequence, e.g. 1dan
    locations    A list of tuples (chain, start, end) describing the location.
    hierarchy    A string specifying where this domain is in the hierarchy.
    
    """
    def __init__(self):
        self.sid = ''
        self.pdbid = ''
        self.locations = []
        self.hierarchy = ''
        
    def __str__(self):
        s = []
        s.append(self.sid)
        s.append(self.pdbid)
        s.append(Location.str(self.locations))
        s.append(self.hierarchy)
        return string.join(s, "\t") + "\n"

class Iterator:
    """Iterates over a DOM file.

    Methods:
    next     Retrieve the next DOM record.
    
    """
    def __init__(self, handle, parser=None):
        """__init__(self, handle, parser=None)

        Create an object that iterates over a 'dom' file.  handle is a
        file-like object.  parser is an optional Parser object to change
        the results into another form.  If set to None, then the raw contents
        of the file will be returned.

        """
        if type(handle) is not FileType and type(handle) is not InstanceType:
            raise ValueError, "I expected a file handle or file-like object"
        self._handle = handle
        self._parser = parser

    def next(self):
        line = self._handle.readline()
        if not line:
            return None
        if self._parser is not None:
            return self._parser.parse(line)
        return line

class DomainParser:
    """Parses dom entries.

    Methods:
    parse        Parse an entry from a DOM file.
    
    """
    def parse(self, entry):
        """parse(self, entry) -> Domain"""
        entry = string.rstrip(entry)  # no trailing whitespace
        columns = string.split(entry, "\t")  # separate the tab-delineated cols
        if len(columns) != 4:
            raise error, "I don't understand the format of %s" % entry
        
        dom = Domain()
        dom.sid, dom.pdbid, locstr, dom.hierarchy = columns
        dom.locations = Location.parse(locstr)
        return dom
