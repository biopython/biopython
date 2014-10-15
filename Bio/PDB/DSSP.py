# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Use the DSSP program to calculate secondary structure and accessibility.

You need to have a working version of DSSP (and a license, free for academic
use) in order to use this. For DSSP, see U{http://swift.cmbi.ru.nl/gv/dssp/}.

The DSSP codes for secondary structure used here are:

    - H        Alpha helix (4-12)
    - B        Isolated beta-bridge residue
    - E        Strand
    - G        3-10 helix
    - I        pi helix
    - T        Turn
    - S        Bend
    - -        None
"""

from __future__ import print_function

__docformat__ = "restructuredtext en"

import re
from Bio._py3k import StringIO
import subprocess

from Bio.Data import SCOPData

from Bio.PDB.AbstractPropertyMap import AbstractResiduePropertyMap
from Bio.PDB.PDBExceptions import PDBException
from Bio.PDB.PDBParser import PDBParser


# Match C in DSSP
_dssp_cys=re.compile('[a-z]')

# Maximal ASA of amino acids
# Values from Sander & Rost, (1994), Proteins, 20:216-226
# Used for relative accessibility
MAX_ACC={}
MAX_ACC["ALA"]=106.0
MAX_ACC["CYS"]=135.0
MAX_ACC["ASP"]=163.0
MAX_ACC["GLU"]=194.0
MAX_ACC["PHE"]=197.0
MAX_ACC["GLY"]=84.0
MAX_ACC["HIS"]=184.0
MAX_ACC["ILE"]=169.0
MAX_ACC["LYS"]=205.0
MAX_ACC["LEU"]=164.0
MAX_ACC["MET"]=188.0
MAX_ACC["ASN"]=157.0
MAX_ACC["PRO"]=136.0
MAX_ACC["GLN"]=198.0
MAX_ACC["ARG"]=248.0
MAX_ACC["SER"]=130.0
MAX_ACC["THR"]=142.0
MAX_ACC["VAL"]=142.0
MAX_ACC["TRP"]=227.0
MAX_ACC["TYR"]=222.0


def ss_to_index(ss):
    """
    Secondary structure symbol to index.
    H=0
    E=1
    C=2
    """
    if ss=='H':
        return 0
    if ss=='E':
        return 1
    if ss=='C':
        return 2
    assert 0


def dssp_dict_from_pdb_file(in_file, DSSP="dssp"):
    """
    Create a DSSP dictionary from a PDB file.

    Example:
    --------
    >>> dssp_dict=dssp_dict_from_pdb_file("1fat.pdb")
    >>> aa, ss, acc=dssp_dict[('A', 1)]

    ::

        @param in_file: pdb file
        @type in_file: string ::

        @param DSSP: DSSP executable (argument to os.system)
        @type DSSP: string ::

        @return: a dictionary that maps (chainid, resid) to
            amino acid type, secondary structure code and
            accessibility.
        @rtype: {}
    """
    # Using universal newlines is important on Python 3, this
    # gives unicode handles rather than bytes handles.
    p = subprocess.Popen([DSSP, in_file], universal_newlines=True,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    out_dict, keys = _make_dssp_dict(StringIO(out))
    return out_dict, keys


def make_dssp_dict(filename):
    """
    Return a DSSP dictionary that maps (chainid, resid) to
    aa, ss and accessibility, from a DSSP file. ::

    @param filename: the DSSP output file
    @type filename: string
    """
    with open(filename, "r") as handle:
        return _make_dssp_dict(handle)


def _make_dssp_dict(handle):
    """
    Return a DSSP dictionary that maps (chainid, resid) to
    aa, ss and accessibility, from an open DSSP file object. ::

    @param handle: the open DSSP output file handle
    @type handle: file
    """
    dssp = {}
    start = 0
    keys = []
    for l in handle.readlines():
        sl = l.split()
        if len(sl) < 2:
            continue
        if sl[1] == "RESIDUE":
            # Start parsing from here
            start = 1
            continue
        if not start:
            continue
        if l[9] == " ":
            # Skip -- missing residue
            continue
        resseq = int(l[5:10])
        icode = l[10]
        chainid = l[11]
        aa = l[13]
        ss = l[16]
        if ss == " ":
            ss = "-"
        try:
            acc = int(l[34:38])
            phi = float(l[103:109])
            psi = float(l[109:115])
        except ValueError as exc:
            # DSSP output breaks its own format when there are >9999
            # residues, since only 4 digits are allocated to the seq num
            # field.  See 3kic chain T res 321, 1vsy chain T res 6077.
            # Here, look for whitespace to figure out the number of extra
            # digits, and shift parsing the rest of the line by that amount.
            if l[34] != ' ':
                shift = l[34:].find(' ')
                acc = int((l[34+shift:38+shift]))
                phi = float(l[103+shift:109+shift])
                psi = float(l[109+shift:115+shift])
            else:
                raise ValueError(exc)
        res_id = (" ", resseq, icode)
        dssp[(chainid, res_id)] = (aa, ss, acc, phi, psi)
        keys.append((chainid, res_id))
    return dssp, keys


class DSSP(AbstractResiduePropertyMap):
    """
    Run DSSP on a pdb file, and provide a handle to the
    DSSP secondary structure and accessibility.

    **Note** that DSSP can only handle one model.

    Example:
    --------

    >>> p = PDBParser()
    >>> structure = p.get_structure("1MOT", "1MOT.pdb")
    >>> model = structure[0]
    >>> dssp = DSSP(model, "1MOT.pdb")
    >>> # DSSP data is accessed by a tuple (chain_id, res_id)
    >>> a_key = list(dssp)[2]
    >>> # residue object, secondary structure, solvent accessibility,
    >>> # relative accessiblity, phi, psi
    >>> dssp[a_key]
    (<Residue ALA het=  resseq=251 icode= >,
    'H',
    72,
    0.67924528301886788,
    -61.200000000000003,
    -42.399999999999999)
    """

    def __init__(self, model, pdb_file, dssp="dssp"):
        """
        ::

        @param model: the first model of the structure
        @type model: L{Model} ::

        @param pdb_file: a PDB file
        @type pdb_file: string ::

        @param dssp: the dssp executable (ie. the argument to os.system)
        @type dssp: string
        """
        # create DSSP dictionary
        dssp_dict, dssp_keys = dssp_dict_from_pdb_file(pdb_file, dssp)
        dssp_map = {}
        dssp_list = []

        def resid2code(res_id):
            """Serialize a residue's resseq and icode for easy comparison."""
            return '%s%s' % (res_id[1], res_id[2])

        # Now create a dictionary that maps Residue objects to
        # secondary structure and accessibility, and a list of
        # (residue, (secondary structure, accessibility)) tuples
        for key in dssp_keys:
            chain_id, res_id = key
            chain = model[chain_id]
            try:
                res = chain[res_id]
            except KeyError:
                # In DSSP, HET field is not considered in residue identifier.
                # Thus HETATM records may cause unnecessary exceptions.
                # (See 3jui chain A res 593.)
                # Try the lookup again with all HETATM other than water
                res_seq_icode = resid2code(res_id)
                for r in chain:
                    if r.id[0] not in (' ', 'W'):
                        # Compare resseq + icode
                        if resid2code(r.id) == res_seq_icode:
                            # Found a matching residue
                            res = r
                            break
                else:
                    raise KeyError(res_id)

            # For disordered residues of point mutations, BioPython uses the
            # last one as default, But DSSP takes the first one (alternative
            # location is blank, A or 1). See 1h9h chain E resi 22.
            # Here we select the res in which all atoms have altloc blank, A or
            # 1. If no such residues are found, simply use the first one appears
            # (as DSSP does).
            if res.is_disordered() == 2:
                for rk in res.disordered_get_id_list():
                    # All atoms in the disordered residue should have the same
                    # altloc, so it suffices to check the altloc of the first
                    # atom.
                    altloc = res.child_dict[rk].get_list()[0].get_altloc()
                    if altloc in tuple('A1 '):
                        res.disordered_select(rk)
                        break
                else:
                    # Simply select the first one
                    res.disordered_select(res.disordered_get_id_list()[0])

            # Sometimes point mutations are put into HETATM and ATOM with altloc
            # 'A' and 'B'.
            # See 3piu chain A residue 273:
            #   <Residue LLP het=H_LLP resseq=273 icode= >
            #   <Residue LYS het=  resseq=273 icode= >
            # DSSP uses the HETATM LLP as it has altloc 'A'
            # We check the altloc code here.
            elif res.is_disordered() == 1:
                # Check altloc of all atoms in the DisorderedResidue. If it
                # contains blank, A or 1, then use it.  Otherwise, look for HET
                # residues of the same seq+icode.  If not such HET residues are
                # found, just accept the current one.
                altlocs = set(a.get_altloc() for a in res.get_unpacked_list())
                if altlocs.isdisjoint('A1 '):
                    # Try again with all HETATM other than water
                    res_seq_icode = resid2code(res_id)
                    for r in chain:
                        if r.id[0] not in (' ', 'W'):
                            if resid2code(r.id) == res_seq_icode and \
                               r.get_list()[0].get_altloc() in tuple('A1 '):
                                res = r
                                break

            aa, ss, acc, phi, psi = dssp_dict[key]
            res.xtra["SS_DSSP"] = ss
            res.xtra["EXP_DSSP_ASA"] = acc
            res.xtra["PHI_DSSP"] = phi
            res.xtra["PSI_DSSP"] = psi
            # Relative accessibility
            resname = res.get_resname()
            try:
                rel_acc = acc/MAX_ACC[resname]
            except KeyError:
                # Invalid value for resname
                rel_acc = 'NA'
            else:
                if rel_acc > 1.0:
                    rel_acc = 1.0
            res.xtra["EXP_DSSP_RASA"] = rel_acc
            # Verify if AA in DSSP == AA in Structure
            # Something went wrong if this is not true!
            # NB: DSSP uses X often
            resname = SCOPData.protein_letters_3to1.get(resname, 'X')
            if resname == "C":
                # DSSP renames C in C-bridges to a,b,c,d,...
                # - we rename it back to 'C'
                if _dssp_cys.match(aa):
                    aa = 'C'
            # Take care of HETATM again
            if (resname != aa) and (res.id[0] == ' ' or aa != 'X'):
                raise PDBException("Structure/DSSP mismatch at %s" % res)
            dssp_map[key] = ((res, ss, acc, rel_acc, phi, psi))
            dssp_list.append((res, ss, acc, rel_acc, phi, psi))

        AbstractResiduePropertyMap.__init__(self, dssp_map, dssp_keys,
                dssp_list)


if __name__ == "__main__":
    import sys

    p = PDBParser()
    s = p.get_structure('X', sys.argv[1])
    model = s[0]
    d = DSSP(model, sys.argv[1])

    for r in d:
        print(r)
    print("Handled %i residues" % len(d))
    print(sorted(d))
    if ('A', 1) in d:
        print(d[('A', 1)])
        print(s[0]['A'][1].xtra)
    # Secondary structure
    print(''.join(item[1] for item in d))
