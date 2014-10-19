# Copyright 2001 by Gavin E. Crooks.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


"""Handle the SCOP HIErarchy files, which describe the SCOP hierarchy in
terms of SCOP unique identifiers (sunid).

The file format is described in the scop
"release notes.":http://scop.berkeley.edu/release-notes-1.55.html
The latest HIE file can be found
"elsewhere at SCOP.":http://scop.mrc-lmb.cam.ac.uk/scop/parse/

"Release 1.55":http://scop.berkeley.edu/parse/dir.hie.scop.txt_1.55 (July 2001)
"""


class Record(object):
    """Holds information for one node in the SCOP hierarchy.

    Attributes:

     - sunid - SCOP unique identifiers of this node
     - parent - Parents sunid
     - children - Sequence of childrens sunids
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
        # For example ::
        #
        # 0       -       46456,48724,51349,53931,56572,56835,56992,57942
        # 21953   49268   -
        # 49267   49266   49268,49269
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
            self.children = [int(x) for x in children]

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

        if self.children:
            s.append(",".join(str(x) for x in self.children))
        else:
            s.append('-')

        return "\t".join(s) + "\n"


def parse(handle):
    """Iterates over a HIE file as Hie records for each line.

    Arguments:

     - handle - file-like object.
    """
    for line in handle:
        if line.startswith('#'):
            continue
        yield Record(line)
