# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Classify protein backbone structure with Kolodny et al's fragment libraries.

It can be regarded as a form of objective secondary structure classification.
Only fragments of length 5 or 7 are supported (ie. there is a 'central'
residue).

Full reference:

Kolodny R, Koehl P, Guibas L, Levitt M.
Small libraries of protein fragments model native protein structures accurately.
J Mol Biol. 2002 323(2):297-307.

The definition files of the fragments can be obtained from:

http://github.com/csblab/fragments/

You need these files to use this module.

The following example uses the library with 10 fragments of length 5.
The library files can be found in directory 'fragment_data'.

    >>> from Bio.PDB.PDBParser import PDBParser
    >>> from Bio.PDB.FragmentMapper import FragmentMapper
    >>> parser = PDBParser()
    >>> structure = parser.get_structure("1a8o", "PDB/1A8O.pdb")
    >>> model = structure[0]
    >>> fm = FragmentMapper(model, lsize=10, flength=5, fdir="PDB")
    >>> chain = model['A']
    >>> res152 = chain[152]
    >>> res157 = chain[157]
    >>> res152 in fm # is res152 mapped? (fragment of a C-alpha polypeptide)
    False
    >>> res157 in fm # is res157 mapped? (fragment of a C-alpha polypeptide)
    True

"""


import numpy as np

from Bio.SVDSuperimposer import SVDSuperimposer

from Bio.PDB.PDBExceptions import PDBException
from Bio.PDB.Polypeptide import PPBuilder


# fragment file (lib_SIZE_z_LENGTH.txt)
# SIZE=number of fragments
# LENGTH=length of fragment (4,5,6,7)
_FRAGMENT_FILE = "lib_%s_z_%s.txt"


def _read_fragments(size, length, dir="."):
    """Read a fragment spec file (PRIVATE).

    Read a fragment spec file available from
    http://github.com/csblab/fragments/
    and return a list of Fragment objects.

    :param size: number of fragments in the library
    :type size: int

    :param length: length of the fragments
    :type length: int

    :param dir: directory where the fragment spec files can be found
    :type dir: string
    """
    filename = (dir + "/" + _FRAGMENT_FILE) % (size, length)
    with open(filename) as fp:
        flist = []
        # ID of fragment=rank in spec file
        fid = 0
        for line in fp:
            # skip comment and blank lines
            if line[0] == "*" or line[0] == "\n":
                continue
            sl = line.split()
            if sl[1] == "------":
                # Start of fragment definition
                f = Fragment(length, fid)
                flist.append(f)
                # increase fragment id (rank)
                fid += 1
                continue
            # Add CA coord to Fragment
            coord = np.array([float(x) for x in sl[0:3]])
            # XXX= dummy residue name
            f.add_residue("XXX", coord)
    return flist


class Fragment:
    """Represent a polypeptide C-alpha fragment."""

    def __init__(self, length, fid):
        """Initialize fragment object.

        :param length: length of the fragment
        :type length: int

        :param fid: id for the fragment
        :type fid: int
        """
        # nr of residues in fragment
        self.length = length
        # nr of residues added
        self.counter = 0
        self.resname_list = []
        # CA coordinate matrix
        self.coords_ca = np.zeros((length, 3), "d")
        self.fid = fid

    def get_resname_list(self):
        """Get residue list.

        :return: the residue names
        :rtype: [string, string,...]
        """
        return self.resname_list

    def get_id(self):
        """Get identifier for the fragment.

        :return: id for the fragment
        :rtype: int
        """
        return self.fid

    def get_coords(self):
        """Get the CA coordinates in the fragment.

        :return: the CA coords in the fragment
        :rtype: NumPy (Nx3) array
        """
        return self.coords_ca

    def add_residue(self, resname, ca_coord):
        """Add a residue.

        :param resname: residue name (eg. GLY).
        :type resname: string

        :param ca_coord: the c-alpha coordinates of the residues
        :type ca_coord: NumPy array with length 3
        """
        if self.counter >= self.length:
            raise PDBException("Fragment boundary exceeded.")
        self.resname_list.append(resname)
        self.coords_ca[self.counter] = ca_coord
        self.counter = self.counter + 1

    def __len__(self):
        """Return length of the fragment."""
        return self.length

    def __sub__(self, other):
        """Return rmsd between two fragments.

        :return: rmsd between fragments
        :rtype: float

        Examples
        --------
        This is an incomplete but illustrative example::

            rmsd = fragment1 - fragment2

        """
        sup = SVDSuperimposer()
        sup.set(self.coords_ca, other.coords_ca)
        sup.run()
        return sup.get_rms()

    def __repr__(self):
        """Represent the fragment object as a string.

        Returns <Fragment length=L id=ID> where L=length of fragment
        and ID the identifier (rank in the library).
        """
        return "<Fragment length=%i id=%i>" % (self.length, self.fid)


def _make_fragment_list(pp, length):
    """Dice up a peptide in fragments of length "length" (PRIVATE).

    :param pp: a list of residues (part of one peptide)
    :type pp: [L{Residue}, L{Residue}, ...]

    :param length: fragment length
    :type length: int
    """
    frag_list = []
    for i in range(len(pp) - length + 1):
        f = Fragment(length, -1)
        for j in range(length):
            residue = pp[i + j]
            resname = residue.get_resname()
            if residue.has_id("CA"):
                ca = residue["CA"]
            else:
                raise PDBException("CHAINBREAK")
            if ca.is_disordered():
                raise PDBException("CHAINBREAK")
            ca_coord = ca.get_coord()
            f.add_residue(resname, ca_coord)
        frag_list.append(f)
    return frag_list


def _map_fragment_list(flist, reflist):
    """Map flist fragments to closest entry in reflist (PRIVATE).

    Map all frgaments in flist to the closest (in RMSD) fragment in reflist.

    Returns a list of reflist indices.

    :param flist: list of protein fragments
    :type flist: [L{Fragment}, L{Fragment}, ...]

    :param reflist: list of reference (ie. library) fragments
    :type reflist: [L{Fragment}, L{Fragment}, ...]
    """
    mapped = []
    for f in flist:
        rank = []
        for i in range(len(reflist)):
            rf = reflist[i]
            rms = f - rf
            rank.append((rms, rf))
        rank.sort()
        fragment = rank[0][1]
        mapped.append(fragment)
    return mapped


class FragmentMapper:
    """Map polypeptides in a model to lists of representative fragments."""

    def __init__(self, model, lsize=20, flength=5, fdir="."):
        """Create instance of FragmentMapper.

        :param model: the model that will be mapped
        :type model: L{Model}

        :param lsize: number of fragments in the library
        :type lsize: int

        :param flength: length of fragments in the library
        :type flength: int

        :param fdir: directory where the definition files are
                     found (default=".")
        :type fdir: string
        """
        if flength == 5:
            self.edge = 2
        elif flength == 7:
            self.edge = 3
        else:
            raise PDBException("Fragment length should be 5 or 7.")
        self.flength = flength
        self.lsize = lsize
        self.reflist = _read_fragments(lsize, flength, fdir)
        self.model = model
        self.fd = self._map(self.model)

    def _map(self, model):
        """Map (PRIVATE).

        :param model: the model that will be mapped
        :type model: L{Model}
        """
        ppb = PPBuilder()
        ppl = ppb.build_peptides(model)
        fd = {}
        for pp in ppl:
            try:
                # make fragments
                flist = _make_fragment_list(pp, self.flength)
                # classify fragments
                mflist = _map_fragment_list(flist, self.reflist)
                for i in range(len(pp)):
                    res = pp[i]
                    if i < self.edge:
                        # start residues
                        continue
                    elif i >= (len(pp) - self.edge):
                        # end residues
                        continue
                    else:
                        # fragment
                        index = i - self.edge
                        assert index >= 0
                        fd[res] = mflist[index]
            except PDBException as why:
                if why == "CHAINBREAK":
                    # Funny polypeptide - skip
                    pass
                else:
                    raise PDBException(why) from None
        return fd

    def __contains__(self, res):
        """Check if the given residue is in any of the mapped fragments.

        :type res: L{Residue}
        """
        return res in self.fd

    def __getitem__(self, res):
        """Get an entry.

        :type res: L{Residue}

        :return: fragment classification
        :rtype: L{Fragment}
        """
        return self.fd[res]
