# Copyright 2001 by Gavin E. Crooks.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


""" Handle the SCOP HIErarchy files, which describe the SCOP hierarchy in
terms of SCOP unique identifiers (sunid).

The file format is described in the scop
"release notes.":http://scop.berkeley.edu/release-notes-1.55.html 
The latest HIE file can be found
"elsewhere at SCOP.":http://scop.mrc-lmb.cam.ac.uk/scop/parse/
  
"Release 1.55":http://scop.berkeley.edu/parse/dir.hie.scop.txt_1.55 (July 2001)
"""


from types import *

class Record:
    """Holds information for one node in the SCOP hierarchy.

    sunid      -- SCOP unique identifiers of this node

    parent     --  Parents sunid

    children   -- Sequence of childrens sunids
    """
    def __init__(self):
        self.sunid = ''
        self.parent = ''
        self.children = []
        
    def __str__(self):
        s = []
        s.append(self.sunid)

        if self.parent :
            s.append(self.parent)
        else:
            s.append('-')

        if self.children :
            s.append(",".join(self.children))
        else:
            s.append('-')

        return "\t".join(s) + "\n"


class Iterator:
    """Iterates over a HIE file.
    """
    def __init__(self, handle, parser=None):
        """Create an object that iterates over a HIE file.

        handle -- file-like object.

        parser -- an optional Parser object to change the results into
                  another form.  If set to None, then the raw contents
                  of the file will be returned.
                  
        """
        if type(handle) is not FileType and type(handle) is not InstanceType:
            raise TypeError, "I expected a file handle or file-like object"
        self._handle = handle
        self._parser = parser

    def next(self):
        """Retrieve the next HIE record."""
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
    """Parses HIE records.

    Records consist of 3 tab deliminated fields; node's sunid,
    parent's sunid, and a list of children's sunids.
    """
    #For example ::
    #
    #0       -       46456,48724,51349,53931,56572,56835,56992,57942
    #21953   49268   -
    #49267   49266   49268,49269
    def parse(self, entry):
        """Returns a Hie Record """
        entry = entry.rstrip()        # no trailing whitespace
        columns = entry.split('\t')   # separate the tab-delineated cols
        if len(columns) != 3:
            raise SyntaxError, "I don't understand the format of %s" % entry
        
        rec = Record()
        rec.sunid, rec.parent, children = columns

        if rec.sunid =='-' : rec.sunid =''
        if rec.parent =='-' : rec.parent =''

        if children =='-' :
            rec.children = ()
        else :
            rec.children = children.split(',')

        return rec


       
    

    








