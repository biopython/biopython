# Copyright 2001 by Gavin E. Crooks.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

""" Handle the SCOP CLAssification file, which describes SCOP domains.

The file format is described in the scop
"release notes.":http://scop.berkeley.edu/release-notes-1.55.html 
The latest CLA file can be found
"elsewhere at SCOP.":http://scop.mrc-lmb.cam.ac.uk/scop/parse/
  
"Release 1.55":http://scop.berkeley.edu/parse/dir.cla.scop.txt_1.55 (July 2001)

"""


from types import *

from Residues import * 
from FileIndex import FileIndex


class Record:
    """Holds information for one SCOP domain

    sid         --  SCOP identifier. e.g. d1danl2

    residues    --  The domain definition as a Residues object

    sccs        --  SCOP concise classification strings.  e.g. b.1.2.1

    sunid       --  SCOP unique identifier for this domain

    hierarchy   --  A sequence of tuples (nodetype, sunid) describing the
                    location of this domain in the SCOP hierarchy.
                    See the Scop module for a description of nodetypes.
    """
    def __init__(self):
        self.sid = ''
        self.residues = None 
        self.sccs = ''
        self.sunid =''
        self.hierarchy = []
        
    def __str__(self):
        s = []
        s.append(self.sid)
        s += str(self.residues).split(" ")
        s.append(self.sccs)
        s.append(self.sunid)

        h=[]
        for ht in self.hierarchy:
             h.append("=".join(map(str,ht))) 
        s.append(",".join(h))
       
        return "\t".join(map(str,s)) + "\n"

class Iterator:
    """Iterates over a CLA file.
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
        """Retrieve the next CLA record."""    
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
    """Parses tab-deliminated CLA records.
    """
    def parse(self, entry):
        """Returns a Cla Record """        
        entry = entry.rstrip()        # no trailing whitespace
        columns = entry.split('\t')   # separate the tab-delineated cols
        if len(columns) != 6:
            raise SyntaxError, "I don't understand the format of %s" % entry
        
        rec = Record()
        rec.sid, pdbid, residues, rec.sccs, rec.sunid, hierarchy = columns
        rec.residues = Residues(residues)
        rec.residues.pdbid = pdbid
        rec.sunid = int(rec.sunid)
        
        h = []
        for ht in hierarchy.split(",") :
            h.append( ht.split('='))        
        for ht in h:
            ht[1] = int(ht[1])
        rec.hierarchy = h
        

        return rec


class Index(FileIndex):
    """A CLA file indexed by SCOP identifiers."""
    def __init__(self, filename, ) :
        iterator = lambda f : Iterator(f, Parser())
        key_gen  = lambda rec : rec.sid

        FileIndex.__init__(self, filename, iterator, key_gen)


       
    

    








