# Copyright 2001 by Gavin E. Crooks.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Handle the SCOP DEScription file.

The file format is described in the scop
"release notes.":http://scop.berkeley.edu/release-notes-1.55.html
The latest DES file can be found
"elsewhere at SCOP.":http://scop.mrc-lmb.cam.ac.uk/scop/parse/

"Release 1.55":http://scop.berkeley.edu/parse/des.cla.scop.txt_1.55 (July 2001)
"""


class Record:
    """Holds information for one node in the SCOP hierarchy.

    Attributes:
     - sunid - SCOP unique identifiers
     - nodetype - One of 'cl' (class), 'cf' (fold), 'sf' (superfamily),
       'fa' (family), 'dm' (protein), 'sp' (species), 'px' (domain).
       Additional node types may be added.
     - sccs - SCOP concise classification strings. e.g. b.1.2.1
     - name - The SCOP ID (sid) for domains (e.g. d1anu1), currently empty for other node types
     - description - e.g. "All beta proteins","Fibronectin type III",

    """

    def __init__(self, line=None):
        """Initialize the class."""
        self.sunid = ""
        self.nodetype = ""
        self.sccs = ""
        self.name = ""
        self.description = ""
        if line:
            self._process(line)

    def _process(self, line):
        """Parse DES records (PRIVATE).

        Records consist of 5 tab deliminated fields,
        sunid, node type, sccs, node name, node description.
        """
        # For example ::
        #
        # 21953   px      b.1.2.1 d1dan.1 1dan T:,U:91-106
        # 48724   cl      b       -       All beta proteins
        # 48725   cf      b.1     -       Immunoglobulin-like beta-sandwich
        # 49265   sf      b.1.2   -       Fibronectin type III
        # 49266   fa      b.1.2.1 -       Fibronectin type III

        line = line.rstrip()  # no trailing whitespace
        columns = line.split("\t")  # separate the tab-delineated cols
        if len(columns) != 5:
            raise ValueError(f"I don't understand the format of {line}")

        sunid, self.nodetype, self.sccs, self.name, self.description = columns
        if self.name == "-":
            self.name = ""
        self.sunid = int(sunid)

    def __str__(self):
        """Represent the SCOP description record as a tab-separated string."""
        s = []
        s.append(self.sunid)
        s.append(self.nodetype)
        s.append(self.sccs)
        if self.name:
            s.append(self.name)
        else:
            s.append("-")
        s.append(self.description)
        return "\t".join(map(str, s)) + "\n"


def parse(handle):
    """Iterate over a DES file as a Des record for each line.

    Arguments:
     - handle - file-like object

    """
    for line in handle:
        if line.startswith("#"):
            continue
        yield Record(line)
