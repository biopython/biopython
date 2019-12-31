# Copyright 2000 by Jeffrey Chang.  All rights reserved.
# Revisions copyright 2001 by Gavin E. Crooks. All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Handle the SCOP DOMain file.

The DOM file has been officially deprecated. For more information see
the SCOP"release notes.":http://scop.berkeley.edu/release-notes-1.55.html
The DOM files for older releases can be found
"elsewhere at SCOP.":http://scop.mrc-lmb.cam.ac.uk/scop/parse/
"""

from .Residues import Residues


class Record:
    """Holds information for one SCOP domain.

    Attribues:
     - sid - The SCOP ID of the entry, e.g. d1anu1
     - residues - The domain definition as a Residues object
     - hierarchy - A string specifying where this domain is in the hierarchy.

    """

    def __init__(self, line=None):
        """Initialize the class."""
        self.sid = ""
        self.residues = []
        self.hierarchy = ""
        if line:
            self._process(line)

    def _process(self, line):
        """Parse DOM records (PRIVATE).

        Records consist of 4 tab deliminated fields;
        sid, pdbid, residues, hierarchy
        """
        # For example ::
        #
        # d1sctg_ 1sct    g:      1.001.001.001.001.001
        # d1scth_ 1sct    h:      1.001.001.001.001.001
        # d1flp__ 1flp    -       1.001.001.001.001.002
        # d1moh__ 1moh    -       1.001.001.001.001.002

        line = line.rstrip()  # no trailing whitespace
        columns = line.split("\t")  # separate the tab-delineated cols
        if len(columns) != 4:
            raise ValueError("I don't understand the format of %s" % line)
        self.sid, pdbid, res, self.hierarchy = columns
        self.residues = Residues(res)
        self.residues.pdbid = pdbid

    def __str__(self):
        """Represent the SCOP domain record as a tab-separated string."""
        s = []
        s.append(self.sid)
        s.append(str(self.residues).replace(" ", "\t"))
        s.append(self.hierarchy)
        return "\t".join(s) + "\n"


def parse(handle):
    """Iterate over a DOM file as a Dom record for each line.

    Arguments:
     - handle -- file-like object.

    """
    for line in handle:
        if line.startswith("#"):
            continue
        yield Record(line)
