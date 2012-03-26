# Copyright 2001 by Gavin E. Crooks.  All rights reserved.
# Modifications Copyright 2010 Jeffrey Finkelstein. All rights reserved.
#
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


class Record(object):
    """Holds information for one SCOP domain.

    sid         --  SCOP identifier. e.g. d1danl2

    residues    --  The domain definition as a Residues object

    sccs        --  SCOP concise classification strings.  e.g. b.1.2.1

    sunid       --  SCOP unique identifier for this domain

    hierarchy   --  A dictionary, keys are nodetype, values are sunid,
                    describing the location of this domain in the SCOP
                    hierarchy. See the Scop module for a description of
                    nodetypes. This used to be a list of (key,value) tuples
                    in older versions of Biopython (see Bug 3109).
    """
    def __init__(self, line=None):
        self.sid = ''
        self.residues = None 
        self.sccs = ''
        self.sunid =''
        self.hierarchy = {}
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
        
        for ht in hierarchy.split(","):
            key, value = ht.split('=')
            self.hierarchy[key] = int(value)

    def __str__(self):
        s = []
        s.append(self.sid)
        s += str(self.residues).split(" ")
        s.append(self.sccs)
        s.append(self.sunid)

        s.append(','.join('='.join((key, str(value))) for key, value
                          in self.hierarchy.iteritems()))

        return "\t".join(map(str,s)) + "\n"


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
        f = open(self.filename, "rU")
        try:
            position = 0
            while True:
                line = f.readline()
                if not line: break
                if line.startswith('#'):
                    continue
                record = Record(line)
                key = record.sid
                if key != None:
                    self[key] = position
                position = f.tell()
        finally:
            f.close()

    def __getitem__(self, key):
        """ Return an item from the indexed file. """
        position = dict.__getitem__(self,key)

        f = open(self.filename, "rU")
        try:
            f.seek(position)
            line = f.readline()
            record = Record(line)
        finally:
            f.close()
        return record
