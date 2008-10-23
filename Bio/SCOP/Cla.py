# Copyright 2001 by Gavin E. Crooks.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

""" Handle the SCOP CLAssification file, which describes SCOP domains.

The file format is described in the scop
"release notes.":http://scop.mrc-lmb.cam.ac.uk/scop/release-notes.html
The latest CLA file can be found
"elsewhere at SCOP.":http://scop.mrc-lmb.cam.ac.uk/scop/parse/
  
"Release 1.73": http://scop.mrc-lmb.cam.ac.uk/scop/parse/dir.cla.scop.txt_1.73
(July 2008)

"""



from Residues import * 


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
    def __init__(self, line=None):
        self.sid = ''
        self.residues = None 
        self.sccs = ''
        self.sunid =''
        self.hierarchy = []
        if line:
            self._process(line)
        
    def _process(self, line):
        line = line.rstrip()         # no trailing whitespace
        columns = line.split('\t')   # separate the tab-delineated cols
        if len(columns) != 6:
            raise ValueError("I don't understand the format of %s" % line)
        
        self.sid, pdbid, residues, self.sccs, self.sunid, hierarchy = columns
        self.residues = Residues(residues)
        self.residues.pdbid = pdbid
        self.sunid = int(self.sunid)
        
        for ht in hierarchy.split(",") :
            key, value = ht.split('=')
            value = int(value)
            self.hierarchy.append([key, value])

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
        """Create an object that iterates over a CLA file.

        handle -- file-like object.

        parser -- an optional Parser object to chang the results into
                  another form.  If set to None, then the raw contents
                  of the file will be returned.

        """
        import warnings
        warnings.warn("Bio.SCOP.Cla.Iterator is deprecated. Please use Bio.SCOP.Cla.parse() instead.", DeprecationWarning)
        from types import FileType, InstanceType
        if type(handle) is not FileType and type(handle) is not InstanceType:
            raise TypeError("I expected a file handle or file-like object")
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
    def __init__(self):
        import warnings
        warnings.warn("""Bio.SCOP.Cla.Parser is deprecated.
        Instead of

        parser = Cla.Parser()
        record = parser.parse(entry)

        please use

        record = Cla.Record(entry)
        """, DeprecationWarning)

    def parse(self, entry):
        """Returns a Cla Record """        
        return Record(entry)


def parse(handle):
    """Iterates over a CLA file, returning a Cla record for each line
    in the file.

    Arguments:

        handle -- file-like object.
    """
    for line in handle:
        if line.startswith('#'):
            continue
        yield Record(line)


class Index(dict):
    """A CLA file indexed by SCOP identifiers, allowing rapid
       random access into a file."""
    def __init__(self, filename):
        """
        Arguments:
        
          filename  -- The file to index
        """
        dict.__init__(self)
        self.filename = filename
        f = open(self.filename)
        try:
            position = 0
            while True:
                line = f.readline()
                if not line: break
                record = Record(line)
                key = record.sid
                if key != None :
                    self[key] = position
                position = f.tell()
        finally:
            f.close()

    def __getitem__(self, key) :
        """ Return an item from the indexed file. """
        position = dict.__getitem__(self,key)

        f = open(self.filename)
        try:
            f.seek(position)
            line = f.readline()
            record = Record(line)
        finally:
            f.close()
        return record
