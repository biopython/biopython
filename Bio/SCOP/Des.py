# Copyright 2001 by Gavin E. Crooks.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


""" Handle the SCOP DEScription file.

The file format is described in the scop
"release notes.":http://scop.berkeley.edu/release-notes-1.55.html 
The latest DES file can be found
"elsewhere at SCOP.":http://scop.mrc-lmb.cam.ac.uk/scop/parse/
  
"Release 1.55":http://scop.berkeley.edu/parse/des.cla.scop.txt_1.55 (July 2001)
"""

from types import *

class Record:
    """Holds information for one node in the SCOP hierarchy.

    sunid       -- SCOP unique identifiers

    nodetype    -- One of 'cl' (class), 'cf' (fold), 'sf' (superfamily),
                   'fa' (family), 'dm' (protein), 'sp' (species),
                   'px' (domain). Additional node types may be added.

    sccs        -- SCOP concise classification strings. e.g. b.1.2.1

    name        -- The SCOP ID (sid) for domains (e.g. d1anu1),
                   currently empty for other node types

    description --  e.g. "All beta proteins","Fibronectin type III", 
    
    """
    def __init__(self):
        self.sunid = ''
        self.nodetype = ''
        self.sccs = ''
        self.name = ''
        self.description =''
        
    def __str__(self):
        s = []
        s.append(self.sunid)
        s.append(self.nodetype)        
        s.append(self.sccs)        
        if self.name :
            s.append(self.name)
        else :
            s.append("-")
        s.append(self.description)        
        return "\t".join(map(str,s)) + "\n"

class Iterator:
    """Iterates over a DES file.
    """
    def __init__(self, handle, parser=None):
        """Create an object that iterates over a DES file.

        handle -- file-like object.

        parser -- an optional Parser object to chang the results into
                  another form.  If set to None, then the raw contents
                  of the file will be returned.
                  
        """
        if type(handle) is not FileType and type(handle) is not InstanceType:
            raise TypeError, "I expected a file handle or file-like object"
        self._handle = handle
        self._parser = parser

    def next(self):
        """Retrieve the next DES record."""
        while 1:
            line = self._handle.readline()
            if not line: return None
            if line[0] !='#':  break  # Not a comment line
        if self._parser is not None :    
            return self._parser.parse(line)
        return line
    
    def __iter__(self):
        return iter(self.next, None)


class Parser:
    """Parses DES records.
    
    Records consist of 5 tab deliminated fields,
    sunid, node type, sccs, node name, node description.
    """
    #For example ::
    #
    #21953   px      b.1.2.1 d1dan.1 1dan T:,U:91-106
    #48724   cl      b       -       All beta proteins
    #48725   cf      b.1     -       Immunoglobulin-like beta-sandwich
    #49265   sf      b.1.2   -       Fibronectin type III
    #49266   fa      b.1.2.1 -       Fibronectin type III

    def parse(self, entry):
        """Returns a Des Record """
        entry = entry.rstrip()  # no trailing whitespace
        columns = entry.split("\t")  # separate the tab-delineated cols
        if len(columns) != 5:
            raise ValueError, "I don't understand the format of %s" % entry
        
        rec = Record()
        rec.sunid, rec.nodetype, rec.sccs, rec.name, rec.description = columns
        if rec.name == '-' : rec.name =''
        rec.sunid = int(rec.sunid)
        return rec




