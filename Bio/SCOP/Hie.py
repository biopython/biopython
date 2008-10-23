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


class Record:
    """Holds information for one node in the SCOP hierarchy.

    sunid      -- SCOP unique identifiers of this node

    parent     --  Parents sunid

    children   -- Sequence of childrens sunids
    """
    def __init__(self, line=None):
        self.sunid = ''
        self.parent = ''
        self.children = []
        if line:
            self._process(line)

    def _process(self, line):
        """Parses HIE records.

        Records consist of 3 tab deliminated fields; node's sunid,
        parent's sunid, and a list of children's sunids.
        """
        #For example ::
        #
        #0       -       46456,48724,51349,53931,56572,56835,56992,57942
        #21953   49268   -
        #49267   49266   49268,49269
        line = line.rstrip()        # no trailing whitespace
        columns = line.split('\t')   # separate the tab-delineated cols
        if len(columns) != 3:
            raise ValueError("I don't understand the format of %s" % line)
        
        sunid, parent, children = columns

        if sunid =='-':
            self.sunid = ''
        else:
            self.sunid = int(sunid)

        if parent=='-':
            self.parent = ''
        else:
            self.parent = int(parent)

        if children=='-':
            self.children = ()
        else:
            children = children.split(',')
            self.children = map(int, children)


    def __str__(self):
        s = []
        s.append(str(self.sunid))

        if self.parent:
            s.append(str(self.parent))
        else:
            if self.sunid != 0:
                s.append('0')
            else:
                s.append('-')
                

        if self.children :
            child_str = map(str, self.children)
            s.append(",".join(child_str))
        else:
            s.append('-')

        return "\t".join(s) + "\n"


def parse(handle):
    """Iterates over a HIE file, returning a Hie record for each line
    in the file.

    Arguments:

        handle -- file-like object.
    """
    for line in handle:
        yield Record(line)


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
        import warnings
        warnings.warn("Bio.SCOP.Hie.Iterator is deprecated. Please use Bio.SCOP.Hie.parse() instead.", DeprecationWarning)
        from types import FileType, InstanceType
        if type(handle) is not FileType and type(handle) is not InstanceType:
            raise TypeError("I expected a file handle or file-like object")
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
    def __init__(self):
        import warnings
        warnings.warn("""Bio.SCOP.Hie.Parser is deprecated.
        Instead of

        parser = Hie.Parser()
        record = parser.parse(entry)

        please use

        record = Hie.Record(entry)
        """, DeprecationWarning)

    def parse(self, entry):
        """Returns a Hie Record """
        return Record(entry)
