# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Polypeptide-related classes (construction and representation).

Simple example with multiple chains,

    >>> from Bio.PDB.PDBParser import PDBParser
    >>> from Bio.PDB.Polypeptide import PPBuilder
    >>> structure = PDBParser().get_structure('2BEG', 'PDB/2BEG.pdb')
    >>> ppb=PPBuilder()
    >>> for pp in ppb.build_peptides(structure):
    ...     print(pp.get_sequence())
    LVFFAEDVGSNKGAIIGLMVGGVVIA
    LVFFAEDVGSNKGAIIGLMVGGVVIA
    LVFFAEDVGSNKGAIIGLMVGGVVIA
    LVFFAEDVGSNKGAIIGLMVGGVVIA
    LVFFAEDVGSNKGAIIGLMVGGVVIA

Example with non-standard amino acids using HETATM lines in the PDB file,
in this case selenomethionine (MSE):

    >>> from Bio.PDB.PDBParser import PDBParser
    >>> from Bio.PDB.Polypeptide import PPBuilder
    >>> structure = PDBParser().get_structure('1A8O', 'PDB/1A8O.pdb')
    >>> ppb=PPBuilder()
    >>> for pp in ppb.build_peptides(structure):
    ...     print(pp.get_sequence())
    DIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNW
    TETLLVQNANPDCKTILKALGPGATLEE
    TACQG

If you want to, you can include non-standard amino acids in the peptides:

    >>> for pp in ppb.build_peptides(structure, aa_only=False):
    ...     print(pp.get_sequence())
    ...     print("%s %s" % (pp.get_sequence()[0], pp[0].get_resname()))
    ...     print("%s %s" % (pp.get_sequence()[-7], pp[-7].get_resname()))
    ...     print("%s %s" % (pp.get_sequence()[-6], pp[-6].get_resname()))
    MDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTILKALGPGATLEEMMTACQG
    M MSE
    M MSE
    M MSE

In this case the selenomethionines (the first and also seventh and sixth from
last residues) have been shown as M (methionine) by the get_sequence method.
"""

import warnings


from Bio.Data.PDBData import nucleic_letters_3to1
from Bio.Data.PDBData import nucleic_letters_3to1_extended
from Bio.Data.PDBData import protein_letters_3to1
from Bio.Data.PDBData import protein_letters_3to1_extended
from Bio.PDB.PDBExceptions import PDBException
from Bio.PDB.vectors import calc_dihedral, calc_angle
from Bio.Seq import Seq


# Sorted by 1-letter code
aa3, aa1 = zip(*sorted(protein_letters_3to1.items(), key=lambda x: x[1]))
standard_aa_names = aa3

d1_to_index = {}
dindex_to_1 = {}
d3_to_index = {}
dindex_to_3 = {}

# Create some lookup tables
for i in range(20):
    n1 = aa1[i]
    n3 = aa3[i]
    d1_to_index[n1] = i
    dindex_to_1[i] = n1
    d3_to_index[n3] = i
    dindex_to_3[i] = n3


def index_to_one(index):
    """Index to corresponding one letter amino acid name.

    >>> index_to_one(0)
    'A'
    >>> index_to_one(19)
    'Y'
    """
    return dindex_to_1[index]


def one_to_index(s):
    """One letter code to index.

    >>> one_to_index('A')
    0
    >>> one_to_index('Y')
    19
    """
    return d1_to_index[s]


def index_to_three(i):
    """Index to corresponding three letter amino acid name.

    >>> index_to_three(0)
    'ALA'
    >>> index_to_three(19)
    'TYR'
    """
    return dindex_to_3[i]


def three_to_index(s):
    """Three letter code to index.

    >>> three_to_index('ALA')
    0
    >>> three_to_index('TYR')
    19
    """
    return d3_to_index[s]


def is_aa(residue, standard=False):
    """Return True if residue object/string is an amino acid.

    :param residue: a L{Residue} object OR a three letter amino acid code
    :type residue: L{Residue} or string

    :param standard: flag to check for the 20 AA (default false)
    :type standard: boolean

    >>> is_aa('ALA')
    True

    Known three letter codes for modified amino acids are supported,

    >>> is_aa('FME')
    True
    >>> is_aa('FME', standard=True)
    False
    """
    if not isinstance(residue, str):
        residue = f"{residue.get_resname():<3s}"
    residue = residue.upper()
    if standard:
        return residue in protein_letters_3to1
    else:
        return residue in protein_letters_3to1_extended


def is_nucleic(residue, standard=False):
    """Return True if residue object/string is a nucleic acid.

    :param residue: a L{Residue} object OR a three letter code
    :type residue: L{Residue} or string

    :param standard: flag to check for the 8 (DNA + RNA) canonical bases.
        Default is False.
    :type standard: boolean

    >>> is_nucleic('DA ')
    True

    >>> is_nucleic('A  ')
    True

    Known three letter codes for modified nucleotides are supported,

    >>> is_nucleic('A2L')
    True
    >>> is_nucleic('A2L', standard=True)
    False
    """
    if not isinstance(residue, str):
        residue = f"{residue.get_resname():<3s}"
    residue = residue.upper()
    if standard:
        return residue in nucleic_letters_3to1
    else:
        return residue in nucleic_letters_3to1_extended


class Polypeptide(list):
    """A polypeptide is simply a list of L{Residue} objects."""

    def get_ca_list(self):
        """Get list of C-alpha atoms in the polypeptide.

        :return: the list of C-alpha atoms
        :rtype: [L{Atom}, L{Atom}, ...]
        """
        ca_list = []
        for res in self:
            ca = res["CA"]
            ca_list.append(ca)
        return ca_list

    def get_phi_psi_list(self):
        """Return the list of phi/psi dihedral angles."""
        ppl = []
        lng = len(self)
        for i in range(lng):
            res = self[i]
            try:
                n = res["N"].get_vector()
                ca = res["CA"].get_vector()
                c = res["C"].get_vector()
            except Exception:
                # Some atoms are missing
                # Phi/Psi cannot be calculated for this residue
                ppl.append((None, None))
                res.xtra["PHI"] = None
                res.xtra["PSI"] = None
                continue
            # Phi
            if i > 0:
                rp = self[i - 1]
                try:
                    cp = rp["C"].get_vector()
                    phi = calc_dihedral(cp, n, ca, c)
                except Exception:
                    phi = None
            else:
                # No phi for residue 0!
                phi = None
            # Psi
            if i < (lng - 1):
                rn = self[i + 1]
                try:
                    nn = rn["N"].get_vector()
                    psi = calc_dihedral(n, ca, c, nn)
                except Exception:
                    psi = None
            else:
                # No psi for last residue!
                psi = None
            ppl.append((phi, psi))
            # Add Phi/Psi to xtra dict of residue
            res.xtra["PHI"] = phi
            res.xtra["PSI"] = psi
        return ppl

    def get_tau_list(self):
        """List of tau torsions angles for all 4 consecutive Calpha atoms."""
        ca_list = self.get_ca_list()
        tau_list = []
        for i in range(len(ca_list) - 3):
            atom_list = (ca_list[i], ca_list[i + 1], ca_list[i + 2], ca_list[i + 3])
            v1, v2, v3, v4 = (a.get_vector() for a in atom_list)
            tau = calc_dihedral(v1, v2, v3, v4)
            tau_list.append(tau)
            # Put tau in xtra dict of residue
            res = ca_list[i + 2].get_parent()
            res.xtra["TAU"] = tau
        return tau_list

    def get_theta_list(self):
        """List of theta angles for all 3 consecutive Calpha atoms."""
        theta_list = []
        ca_list = self.get_ca_list()
        for i in range(len(ca_list) - 2):
            atom_list = (ca_list[i], ca_list[i + 1], ca_list[i + 2])
            v1, v2, v3 = (a.get_vector() for a in atom_list)
            theta = calc_angle(v1, v2, v3)
            theta_list.append(theta)
            # Put tau in xtra dict of residue
            res = ca_list[i + 1].get_parent()
            res.xtra["THETA"] = theta
        return theta_list

    def get_sequence(self):
        """Return the AA sequence as a Seq object.

        :return: polypeptide sequence
        :rtype: L{Seq}
        """
        s = "".join(
            protein_letters_3to1_extended.get(res.get_resname(), "X") for res in self
        )
        return Seq(s)

    def __repr__(self):
        """Return string representation of the polypeptide.

        Return <Polypeptide start=START end=END>, where START
        and END are sequence identifiers of the outer residues.
        """
        start = self[0].get_id()[1]
        end = self[-1].get_id()[1]
        return f"<Polypeptide start={start} end={end}>"


class _PPBuilder:
    """Base class to extract polypeptides.

    It checks if two consecutive residues in a chain are connected.
    The connectivity test is implemented by a subclass.

    This assumes you want both standard and non-standard amino acids.
    """

    def __init__(self, radius):
        """Initialize the base class.

        :param radius: distance
        :type radius: float
        """
        self.radius = radius

    def _accept(self, residue, standard_aa_only):
        """Check if the residue is an amino acid (PRIVATE)."""
        if is_aa(residue, standard=standard_aa_only):
            return True
        elif not standard_aa_only and "CA" in residue.child_dict:
            # It has an alpha carbon...
            # We probably need to update the hard coded list of
            # non-standard residues, see function is_aa for details.
            warnings.warn(
                "Assuming residue %s is an unknown modified amino acid"
                % residue.get_resname()
            )
            return True
        else:
            # not a standard AA so skip
            return False

    def build_peptides(self, entity, aa_only=1):
        """Build and return a list of Polypeptide objects.

        :param entity: polypeptides are searched for in this object
        :type entity: L{Structure}, L{Model} or L{Chain}

        :param aa_only: if 1, the residue needs to be a standard AA
        :type aa_only: int
        """
        is_connected = self._is_connected
        accept = self._accept
        level = entity.get_level()
        # Decide which entity we are dealing with
        if level == "S":
            model = entity[0]
            chain_list = model.get_list()
        elif level == "M":
            chain_list = entity.get_list()
        elif level == "C":
            chain_list = [entity]
        else:
            raise PDBException("Entity should be Structure, Model or Chain.")
        pp_list = []
        for chain in chain_list:
            chain_it = iter(chain)
            try:
                prev_res = next(chain_it)
                while not accept(prev_res, aa_only):
                    prev_res = next(chain_it)
            except StopIteration:
                # No interesting residues at all in this chain
                continue
            pp = None
            for next_res in chain_it:
                if (
                    accept(prev_res, aa_only)
                    and accept(next_res, aa_only)
                    and is_connected(prev_res, next_res)
                ):
                    if pp is None:
                        pp = Polypeptide()
                        pp.append(prev_res)
                        pp_list.append(pp)
                    pp.append(next_res)
                else:
                    # Either too far apart, or one of the residues is unwanted.
                    # End the current peptide
                    pp = None
                prev_res = next_res
        return pp_list


class CaPPBuilder(_PPBuilder):
    """Use CA--CA distance to find polypeptides."""

    def __init__(self, radius=4.3):
        """Initialize the class."""
        _PPBuilder.__init__(self, radius)

    def _is_connected(self, prev_res, next_res):
        for r in [prev_res, next_res]:
            if not r.has_id("CA"):
                return False
        n = next_res["CA"]
        p = prev_res["CA"]
        # Unpack disordered
        if n.is_disordered():
            nlist = n.disordered_get_list()
        else:
            nlist = [n]
        if p.is_disordered():
            plist = p.disordered_get_list()
        else:
            plist = [p]
        for nn in nlist:
            for pp in plist:
                if (nn - pp) < self.radius:
                    return True
        return False


class PPBuilder(_PPBuilder):
    """Use C--N distance to find polypeptides."""

    def __init__(self, radius=1.8):
        """Initialize the class."""
        _PPBuilder.__init__(self, radius)

    def _is_connected(self, prev_res, next_res):
        if not prev_res.has_id("C"):
            return False
        if not next_res.has_id("N"):
            return False
        test_dist = self._test_dist
        c = prev_res["C"]
        n = next_res["N"]
        # Test all disordered atom positions!
        if c.is_disordered():
            clist = c.disordered_get_list()
        else:
            clist = [c]
        if n.is_disordered():
            nlist = n.disordered_get_list()
        else:
            nlist = [n]
        for nn in nlist:
            for cc in clist:
                # To form a peptide bond, N and C must be
                # within radius and have the same altloc
                # identifier or one altloc blank
                n_altloc = nn.get_altloc()
                c_altloc = cc.get_altloc()
                if n_altloc == c_altloc or n_altloc == " " or c_altloc == " ":
                    if test_dist(nn, cc):
                        # Select the disordered atoms that
                        # are indeed bonded
                        if c.is_disordered():
                            c.disordered_select(c_altloc)
                        if n.is_disordered():
                            n.disordered_select(n_altloc)
                        return True
        return False

    def _test_dist(self, c, n):
        """Return 1 if distance between atoms<radius (PRIVATE)."""
        if (c - n) < self.radius:
            return 1
        else:
            return 0
