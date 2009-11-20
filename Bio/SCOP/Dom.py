# Copyright 2000 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# Gavin E. Crooks 2001-11-07:
#     Interface and comments modified to reflect changes to the SCOP
#     module, and to SCOP itself.

""" Handle the SCOP DOMain file.

The DOM file has been officially deprecated. For more information see
the SCOP"release notes.":http://scop.berkeley.edu/release-notes-1.55.html 
The DOM files for older releases can be found 
"elsewhere at SCOP.":http://scop.mrc-lmb.cam.ac.uk/scop/parse/
"""


from Residues import Residues

class Record:
    """Holds information for one SCOP domain.

    sid -- The SCOP ID of the entry, e.g. d1anu1

    residues -- The domain definition as a Residues object

    hierarchy -- A string specifying where this domain is in the hierarchy.
    """
    def __init__(self, line=None):
        self.sid = ''
        self.residues = []
        self.hierarchy = ''
        if line:
            self._process(line)

    def _process(self, line):
        """Parses DOM records.
    
        Records consist of 4 tab deliminated fields;
        sid, pdbid, residues, hierarchy
        """
        #For example ::
        #
        #d1sctg_ 1sct    g:      1.001.001.001.001.001
        #d1scth_ 1sct    h:      1.001.001.001.001.001
        #d1flp__ 1flp    -       1.001.001.001.001.002
        #d1moh__ 1moh    -       1.001.001.001.001.002

        line = line.rstrip()  # no trailing whitespace
        columns = line.split("\t")  # separate the tab-delineated cols
        if len(columns) != 4:
            raise ValueError("I don't understand the format of %s" % line)
        self.sid, pdbid, res, self.hierarchy = columns
        self.residues = Residues(res)
        self.residues.pdbid =pdbid


    def __str__(self):
        s = []
        s.append(self.sid)
        s.append(str(self.residues).replace(" ","\t") )
        s.append(self.hierarchy)
        return "\t".join(s) + "\n"


def parse(handle):
    """Iterates over a DOM file, returning a Dom record for each line
    in the file.

    Arguments:
        
        handle -- file-like object.
    """ 
    for line in handle:
        yield Record(line)

    
class Iterator:
    """Iterates over a DOM file.
    """
    def __init__(self, handle, parser=None):
        """Create an object that iterates over a DES file.

        handle -- file-like object.

        parser -- an optional Parser object to change the results into
                  another form.  If set to None, then the raw contents
                  of the file will be returned.
                  
        """
        import warnings
        warnings.warn("Bio.SCOP.Dom.Iterator is deprecated. Please use Bio.SCOP.Dom.parse() instead.", DeprecationWarning)
        from types import FileType, InstanceType
        if type(handle) is not FileType and type(handle) is not InstanceType:
            raise ValueError("I expected a file handle or file-like object")
        self._handle = handle
        self._parser = parser

    def next(self):
        line = self._handle.readline()
        if not line:
            return None
        if self._parser is not None:
            return self._parser.parse(line)
        return line

class Parser:
    def parse(self, entry):
        """Returns a Dom.Record """
        import warnings
        warnings.warn("""Bio.SCOP.Dom.Parser is deprecated.
        Instead of

        parser = Dom.Parser()
        record = parser.parse(entry)

        please use

        record = Dom.Record(entry)
        """, DeprecationWarning)
        return Record(entry)
