# Copyright 2000 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Lin.py

This module provides code to work with the 'lin' files from SCOP.
http://scop.mrc-lmb.cam.ac.uk/scop/

Classes:
Lin         Holds one record from a 'lin' file.
Iterator    Iterate over all records in a 'lin' file.
LinParser   Parse one record from a 'lin' file.

"""
import string
from types import *

from Bio import File
from Bio.ParserSupport import *
import Location

class Lin:
    """Holds information from a 'lin' file.

    Members:
    pdb          The PDB name for the sequence, e.g. 3sdh
    domain       The SCOP domain name for the sequence, e.g. d3sdha_
    hierarchy    The SCOP hierarchy, e.g. 1.001.001.001.001.001
    locations    A list of tuples (chain, start, end) describing the location.
    klass        SCOP class.
    fold         SCOP fold.
    superfamily  SCOP superfamily.
    family       SCOP family.
    protein      Description of the protein
    species      Species.
    
    """
    def __init__(self):
        self.pdb = ''
        self.domain = ''
        self.hierarchy = ''
        self.locations = []
        self.klass = ''
        self.fold = ''
        self.superfamily = ''
        self.family = ''
        self.protein = ''
        self.species = ''

    def __str__(self):
        lines = []
        lines.append("%s " % self.domain)
        lines.append("%s\t" % self.pdb)
        lines.append("%s\t" % Location.str(self.locations))
        lines.append("%s\t" % self.hierarchy)
        lines.append("Class: %s| " % self.klass)
        lines.append("Fold: %s| " % self.fold)
        lines.append("Superfamily: %s| " % self.superfamily)
        lines.append("Family: %s| " % self.family)
        lines.append("Protein: %s| " % self.protein)
        lines.append("Species: %s|" % self.species)
        return string.join(lines, "") + "\n"

class Iterator:
    """Iterates over a LIN file.

    Methods:
    next     Retrieve the next LIN record.
    
    """
    def __init__(self, handle, parser=None):
        """__init__(self, handle, parser=None)

        Create an object that iterates over a 'lin' file.  handle is a
        file-like object.  parser is an optional Parser object to change
        the results into another form.  If set to None, then the raw contents
        of the file will be returned.

        """
        if type(handle) is not FileType and type(handle) is not InstanceType:
            raise ValueError, "I expected a file handle or file-like object"
        self._handle = File.UndoHandle(handle)
        self._parser = parser

    def next(self):
        """S.next() -> obj or None"""
        # Return each record as a single line.
        lines = []
        while 1:
            line = self._handle.readline()
            if not line:
                break
            # If this is the first line I've read, make sure I'm at the
            # beginning of a record.
            if not lines:
                if string.find(line, 'Class:') < 0:
                    raise SyntaxError, "I don't see a record:\n%s" % line
            # If I'm at the beginning of the next record, then save the
            # line and stop reading.
            if lines and string.find(line, 'Class:') >= 0:
                self._handle.saveline(line)
                break
            # Strip off EOL characters.  XXX should optimize this.
            while line and line[-1] in ['\n', '\r']:
                line = line[:-1]
            lines.append(line)
        if not lines:
            return None
        line = string.join(lines, '')
        if self._parser is not None:
            return self._parser.parse(line)
        return line

class LinParser:
    """Parse information from a 'lin' file.

    Members:
    parse      Parse a record from a 'lin' file.
    
    """
    def parse(self, line):
        """S.parse(line) -> Lin

        data should be a string containing an entry from a 'lin' file.

        """
        # Sample record:
        # d1hbsg_ 1hbs    g:      1.001.001.001.001.016   Class: All alpha| \
        # Fold: Globin-like| Superfamily: Globin-like| Family: Globins| \
        # Protein: Hemoglobin, alpha-chain| Species: human (<i>Homo \
        # sapiens</i>)|
        rec = Lin()
        try:
            i = string.index(line, 'Class:')
        except ValueError:
            raise SyntaxError, "I couldn't find 'Class:' in\n%s" % line
        location, hierarchy = line[:i], line[i:]
        self._parse_info(location, rec)
        self._parse_hierarchy(hierarchy, rec)
        return rec

    def _parse_info(self, data, rec):
        cols = string.split(data)
        try:
            rec.domain, rec.pdb, locstr, rec.hierarchy = cols
        except ValueError:
            raise SyntaxError, "I expected 4 columns in %s" % data
        rec.locations = Location.parse(locstr)

    def _parse_hierarchy(self, data, rec):
        cols = string.split(data, '|')
        cols = map(string.strip, cols)
        if len(cols) != 7 or cols[-1] != '':
            raise SyntaxError, data
        labels = [('Class', 'klass'),
                  ('Fold', 'fold'),
                  ('Superfamily', 'superfamily'),
                  ('Family', 'family'),
                  ('Protein', 'protein'),
                  ('Species', 'species')
                  ]
        for i in range(len(labels)):
            label, attr = labels[i]
            value = self._clean(cols[i], label)
            setattr(rec, attr, value)

    def _clean(self, str, label):
        try:
            i = string.index(str, ':')
        except ValueError:
            raise SyntaxError, "I could not find the label in %s" % str
        l, data = str[:i], str[i+1:]
        if l != label:
            raise SyntaxError, "Expected %s but got %s" % (label, l)
        return string.lstrip(data)
