# Copyright 2001 by Gavin E. Crooks.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# Gavin E. Crooks 2001-10-10

"""ASTRAL RAF (Rapid Access Format) Sequence Maps.

The ASTRAL RAF Sequence Maps record the relationship between the PDB SEQRES
records (representing the sequence of the molecule used in an experiment) to
the ATOM records (representing the atoms experimentally observed).

This data is derived from the Protein Data Bank CIF files. Known errors in the
CIF files are corrected manually, with the original PDB file serving as the
final arbiter in case of discrepancies.

Residues are referenced by residue ID. This consists of a the PDB residue
sequence number (up to 4 digits) and an optional PDB insertion code (an
ascii alphabetic character, a-z, A-Z). e.g. "1", "10A", "1010b", "-1"

See "ASTRAL RAF Sequence Maps":http://astral.stanford.edu/raf.html

protein_letters_3to1 -- A mapping from the 3-letter amino acid codes found
                        in PDB files to 1-letter codes.  The 3-letter codes
                        include chemically modified residues.
"""

from __future__ import print_function
from Bio._py3k import basestring
from Bio._py3k import _universal_read_mode

from copy import copy

from Bio.Data.SCOPData import protein_letters_3to1

from Bio.SCOP.Residues import Residues

__docformat__ = "restructuredtext en"


def normalize_letters(one_letter_code):
    """Convert RAF one-letter amino acid codes into IUPAC standard codes.

    Letters are uppercased, and "." ("Unknown") is converted to "X".
    """
    if one_letter_code == '.':
        return 'X'
    else:
        return one_letter_code.upper()


class SeqMapIndex(dict):
    """An RAF file index.

    The RAF file itself is about 50 MB. This index provides rapid, random
    access of RAF records without having to load the entire file into memory.

    The index key is a concatenation of the  PDB ID and chain ID. e.g
    "2drcA", ``"155c_"``. RAF uses an underscore to indicate blank
    chain IDs.
    """

    def __init__(self, filename):
        """
        Arguments:

          filename  -- The file to index
        """
        dict.__init__(self)
        self.filename = filename

        with open(self.filename, _universal_read_mode) as f:
            position = 0
            while True:
                line = f.readline()
                if not line:
                    break
                key = line[0:5]
                if key is not None:
                    self[key] = position
                position = f.tell()

    def __getitem__(self, key):
        """ Return an item from the indexed file. """
        position = dict.__getitem__(self, key)

        with open(self.filename, _universal_read_mode) as f:
            f.seek(position)
            line = f.readline()
            record = SeqMap(line)
        return record

    def getSeqMap(self, residues):
        """Get the sequence map for a collection of residues.

        residues -- A Residues instance, or a string that can be converted into
                    a Residues instance.
        """
        if isinstance(residues, basestring):
            residues = Residues(residues)

        pdbid = residues.pdbid
        frags = residues.fragments
        if not frags:
            frags = (('_', '', ''),)  # All residues of unnamed chain

        seqMap = None
        for frag in frags:
            chainid = frag[0]
            if chainid in ['', '-', ' ', '_']:
                chainid = '_'
            id = pdbid + chainid

            sm = self[id]

            # Cut out fragment of interest
            start = 0
            end = len(sm.res)
            if frag[1]:
                start = int(sm.index(frag[1], chainid))
            if frag[2]:
                end = int(sm.index(frag[2], chainid)) + 1

            sm = sm[start:end]

            if seqMap is None:
                seqMap = sm
            else:
                seqMap += sm

        return seqMap


class SeqMap(object):
    """An ASTRAL RAF (Rapid Access Format) Sequence Map.

    This is a list like object; You can find the location of particular residues
    with index(), slice this SeqMap into fragments, and glue fragments back
    together with extend().

    pdbid -- The PDB 4 character ID

    pdb_datestamp -- From the PDB file

    version -- The RAF format version. e.g. 0.01

    flags -- RAF flags. (See release notes for more information.)

    res -- A list of Res objects, one for each residue in this sequence map
    """

    def __init__(self, line=None):
        self.pdbid = ''
        self.pdb_datestamp = ''
        self.version = ''
        self.flags = ''
        self.res = []
        if line:
            self._process(line)

    def _process(self, line):
        """Parses a RAF record into a SeqMap object.
        """
        header_len = 38

        line = line.rstrip()  # no trailing whitespace

        if len(line) < header_len:
            raise ValueError("Incomplete header: " + line)

        self.pdbid = line[0:4]
        chainid = line[4:5]

        self.version = line[6:10]

        # Raf format versions 0.01 and 0.02 are identical for practical purposes
        if(self.version != "0.01" and self.version != "0.02"):
            raise ValueError("Incompatible RAF version: " + self.version)

        self.pdb_datestamp = line[14:20]
        self.flags = line[21:27]

        for i in range(header_len, len(line), 7):
            f = line[i:i + 7]
            if len(f) != 7:
                raise ValueError("Corrupt Field: (" + f + ")")
            r = Res()
            r.chainid = chainid
            r.resid = f[0:5].strip()
            r.atom = normalize_letters(f[5:6])
            r.seqres = normalize_letters(f[6:7])

            self.res.append(r)

    def index(self, resid, chainid="_"):
        for i in range(0, len(self.res)):
            if self.res[i].resid == resid and self.res[i].chainid == chainid:
                return i
        raise KeyError("No such residue " + chainid + resid)

    def __getitem__(self, index):
        if not isinstance(index, slice):
            raise NotImplementedError
        s = copy(self)
        s.res = s.res[index]
        return s

    def append(self, res):
        """Append another Res object onto the list of residue mappings."""
        self.res.append(res)

    def extend(self, other):
        """Append another SeqMap onto the end of self.

        Both SeqMaps must have the same PDB ID, PDB datestamp and
        RAF version.  The RAF flags are erased if they are inconsistent. This
        may happen when fragments are taken from different chains.
        """
        if not isinstance(other, SeqMap):
            raise TypeError("Can only extend a SeqMap with a SeqMap.")
        if self.pdbid != other.pdbid:
            raise TypeError("Cannot add fragments from different proteins")
        if self.version != other.version:
            raise TypeError("Incompatible rafs")
        if self.pdb_datestamp != other.pdb_datestamp:
            raise TypeError("Different pdb dates!")
        if self.flags != other.flags:
            self.flags = ''
        self.res += other.res

    def __iadd__(self, other):
        self.extend(other)
        return self

    def __add__(self, other):
        s = copy(self)
        s.extend(other)
        return s

    def getAtoms(self, pdb_handle, out_handle):
        """Extract all relevant ATOM and HETATOM records from a PDB file.

        The PDB file is scanned for ATOM and HETATOM records. If the
        chain ID, residue ID (seqNum and iCode), and residue type match
        a residue in this sequence map, then the record is echoed to the
        output handle.

        This is typically used to find the coordinates of a domain, or other
        residue subset.

        pdb_handle -- A handle to the relevant PDB file.

        out_handle -- All output is written to this file like object.
        """
        # This code should be refactored when (if?) biopython gets a PDB parser

        # The set of residues that I have to find records for.
        resSet = {}
        for r in self.res:
            if r.atom == 'X':  # Unknown residue type
                continue
            chainid = r.chainid
            if chainid == '_':
                chainid = ' '
            resid = r.resid
            resSet[(chainid, resid)] = r

        resFound = {}
        for line in pdb_handle:
            if line.startswith("ATOM  ") or line.startswith("HETATM"):
                chainid = line[21:22]
                resid = line[22:27].strip()
                key = (chainid, resid)
                if key in resSet:
                    res = resSet[key]
                    atom_aa = res.atom
                    resName = line[17:20]
                    if resName in protein_letters_3to1:
                        if protein_letters_3to1[resName] == atom_aa:
                            out_handle.write(line)
                            resFound[key] = res

        if len(resSet) != len(resFound):
            # for k in resFound:
            #    del resSet[k]
            # print(resSet)

            raise RuntimeError('I could not find at least one ATOM or HETATM'
                   + ' record for each and every residue in this sequence map.')


class Res(object):
    """ A single residue mapping from a RAF record.

    chainid -- A single character chain ID.

    resid   -- The residue ID.

    atom    -- amino acid one-letter code from ATOM records.

    seqres  -- amino acid one-letter code from SEQRES records.
    """
    def __init__(self):
        self.chainid = ''
        self.resid = ''
        self.atom = ''
        self.seqres = ''


def parse(handle):
    """Iterates over a RAF file, returning a SeqMap object for each line
    in the file.

    Arguments:

        handle -- file-like object.
    """
    for line in handle:
        yield SeqMap(line)
