# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

r"""Use the DSSP program to calculate secondary structure and accessibility.

You need to have a working version of DSSP (and a license, free for academic
use) in order to use this. For DSSP, see http://swift.cmbi.ru.nl/gv/dssp/.

The following Accessible surface area (ASA) values can be used, defaulting
to the Sander and Rost values:

    Miller
        Miller et al. 1987 https://doi.org/10.1016/0022-2836(87)90038-6
    Sander
        Sander and Rost 1994 https://doi.org/10.1002/prot.340200303
    Wilke
        Tien et al. 2013 https://doi.org/10.1371/journal.pone.0080635

The DSSP codes for secondary structure used here are:

    =====     ====
    Code      Structure
    =====     ====
     H        Alpha helix (4-12)
     B        Isolated beta-bridge residue
     E        Strand
     G        3-10 helix
     I        Pi helix
     T        Turn
     S        Bend
     \-       None
    =====     ====

Usage
-----
The DSSP class can be used to run DSSP on a pdb file, and provides a
handle to the DSSP secondary structure and accessibility.

**Note** that DSSP can only handle one model, and will only run
calculations on the first model in the provided PDB file.

Examples
--------
>>> p = PDBParser()
>>> structure = p.get_structure("1MOT", "1mot.pdb")
>>> model = structure[0]
>>> dssp = DSSP(model, "1mot.pdb")

Note that the recent DSSP executable from the DSSP-2 package was
renamed from `dssp` to `mkdssp`. If using a recent DSSP release,
you may need to provide the name of your DSSP executable:

>>> dssp = DSSP(model, '1mot.pdb', dssp='mkdssp')

DSSP data is accessed by a tuple - (chain id, residue id):

>>> a_key = list(dssp.keys())[2]
>>> a_key
('A', (' ', 251, ' '))
>>> dssp[a_key]
(3, 'A', 'H', 0.7075471698113207, -61.2, -42.4,
 -2, -0.7, 4, -3.0, 1, -0.2, 5, -0.2)

The dssp data returned for a single residue is a tuple in the form:

    ============ ===
    Tuple Index  Value
    ============ ===
    0            DSSP index
    1            Amino acid
    2            Secondary structure
    3            Relative ASA
    4            Phi
    5            Psi
    6            NH-->O_1_relidx
    7            NH-->O_1_energy
    8            O-->NH_1_relidx
    9            O-->NH_1_energy
    10           NH-->O_2_relidx
    11           NH-->O_2_energy
    12           O-->NH_2_relidx
    13           O-->NH_2_energy
    ============ ===

"""

from __future__ import print_function

import re
from Bio._py3k import StringIO
import subprocess
import warnings

from Bio.PDB.AbstractPropertyMap import AbstractResiduePropertyMap
from Bio.PDB.PDBExceptions import PDBException
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import three_to_one

# Match C in DSSP
_dssp_cys = re.compile('[a-z]')

# Maximal ASA of amino acids
# Used for relative accessibility

residue_max_acc = {
    # Miller max acc: Miller et al. 1987 https://doi.org/10.1016/0022-2836(87)90038-6
    # Wilke: Tien et al. 2013 https://doi.org/10.1371/journal.pone.0080635
    # Sander: Sander & Rost 1994 https://doi.org/10.1002/prot.340200303
    'Miller': {
        'ALA': 113.0, 'ARG': 241.0, 'ASN': 158.0, 'ASP': 151.0,
        'CYS': 140.0, 'GLN': 189.0, 'GLU': 183.0, 'GLY': 85.0,
        'HIS': 194.0, 'ILE': 182.0, 'LEU': 180.0, 'LYS': 211.0,
        'MET': 204.0, 'PHE': 218.0, 'PRO': 143.0, 'SER': 122.0,
        'THR': 146.0, 'TRP': 259.0, 'TYR': 229.0, 'VAL': 160.0
    },
    'Wilke': {
        'ALA': 129.0, 'ARG': 274.0, 'ASN': 195.0, 'ASP': 193.0,
        'CYS': 167.0, 'GLN': 225.0, 'GLU': 223.0, 'GLY': 104.0,
        'HIS': 224.0, 'ILE': 197.0, 'LEU': 201.0, 'LYS': 236.0,
        'MET': 224.0, 'PHE': 240.0, 'PRO': 159.0, 'SER': 155.0,
        'THR': 172.0, 'TRP': 285.0, 'TYR': 263.0, 'VAL': 174.0
    },
    'Sander': {
        'ALA': 106.0, 'ARG': 248.0, 'ASN': 157.0, 'ASP': 163.0,
        'CYS': 135.0, 'GLN': 198.0, 'GLU': 194.0, 'GLY': 84.0,
        'HIS': 184.0, 'ILE': 169.0, 'LEU': 164.0, 'LYS': 205.0,
        'MET': 188.0, 'PHE': 197.0, 'PRO': 136.0, 'SER': 130.0,
        'THR': 142.0, 'TRP': 227.0, 'TYR': 222.0, 'VAL': 142.0
    }
}


def ss_to_index(ss):
    """Secondary structure symbol to index.

    H=0
    E=1
    C=2
    """
    if ss == 'H':
        return 0
    if ss == 'E':
        return 1
    if ss == 'C':
        return 2
    assert 0


def dssp_dict_from_pdb_file(in_file, DSSP="dssp"):
    """Create a DSSP dictionary from a PDB file.

    Parameters
    ----------
    in_file : string
        pdb file

    DSSP : string
        DSSP executable (argument to os.system)

    Returns
    -------
    (out_dict, keys) : tuple
        a dictionary that maps (chainid, resid) to
        amino acid type, secondary structure code and
        accessibility.

    Examples
    --------
    >>> dssp_dict=dssp_dict_from_pdb_file("1fat.pdb")
    >>> aa, ss, acc=dssp_dict[('A', 1)]

    """
    # Using universal newlines is important on Python 3, this
    # gives unicode handles rather than bytes handles.
    # Newer version of DSSP executable is named 'mkdssp',
    # and calling 'dssp' will hence not work in some operating systems
    # (Debian distribution of DSSP includes a symlink for 'dssp' argument)
    try:
        p = subprocess.Popen([DSSP, in_file], universal_newlines=True,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except OSError:  # TODO: Use FileNotFoundError once drop Python 2
        if DSSP == "mkdssp":
            raise
        p = subprocess.Popen(["mkdssp", in_file], universal_newlines=True,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()

    # Alert user for errors
    if err.strip():
        warnings.warn(err)
        if not out.strip():
            raise Exception('DSSP failed to produce an output')

    out_dict, keys = _make_dssp_dict(StringIO(out))
    return out_dict, keys


def make_dssp_dict(filename):
    """DSSP dictionary mapping identifiers to properties.

    Return a DSSP dictionary that maps (chainid, resid) to
    aa, ss and accessibility, from a DSSP file.

    Parameters
    ----------
    filename : string
        the DSSP output file

    """
    with open(filename, "r") as handle:
        return _make_dssp_dict(handle)


def _make_dssp_dict(handle):
    """Return a DSSP dictionary, used by mask_dssp_dict (PRIVATE).

    DSSP dictionary maps (chainid, resid) to an amino acid,
    secondary structure symbol, solvent accessibility value, and hydrogen bond
    information (relative dssp indices and hydrogen bond energies) from an open
    DSSP file object.

    Parameters
    ----------
    handle : file
        the open DSSP output file handle

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

        dssp_index = int(l[:5])
        resseq = int(l[5:10])
        icode = l[10]
        chainid = l[11]
        aa = l[13]
        ss = l[16]
        if ss == " ":
            ss = "-"
        try:
            NH_O_1_relidx = int(l[38:45])
            NH_O_1_energy = float(l[46:50])
            O_NH_1_relidx = int(l[50:56])
            O_NH_1_energy = float(l[57:61])
            NH_O_2_relidx = int(l[61:67])
            NH_O_2_energy = float(l[68:72])
            O_NH_2_relidx = int(l[72:78])
            O_NH_2_energy = float(l[79:83])

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

                NH_O_1_relidx = int(l[38 + shift:45 + shift])
                NH_O_1_energy = float(l[46 + shift:50 + shift])
                O_NH_1_relidx = int(l[50 + shift:56 + shift])
                O_NH_1_energy = float(l[57 + shift:61 + shift])
                NH_O_2_relidx = int(l[61 + shift:67 + shift])
                NH_O_2_energy = float(l[68 + shift:72 + shift])
                O_NH_2_relidx = int(l[72 + shift:78 + shift])
                O_NH_2_energy = float(l[79 + shift:83 + shift])

                acc = int((l[34 + shift:38 + shift]))
                phi = float(l[103 + shift:109 + shift])
                psi = float(l[109 + shift:115 + shift])
            else:
                raise ValueError(exc)
        res_id = (" ", resseq, icode)
        dssp[(chainid, res_id)] = (aa, ss, acc, phi, psi, dssp_index,
                NH_O_1_relidx, NH_O_1_energy, O_NH_1_relidx, O_NH_1_energy,
                NH_O_2_relidx, NH_O_2_energy, O_NH_2_relidx, O_NH_2_energy)
        keys.append((chainid, res_id))
    return dssp, keys


class DSSP(AbstractResiduePropertyMap):
    """Run DSSP and parse secondary structure and accessibility.

    Run DSSP on a pdb file, and provide a handle to the
    DSSP secondary structure and accessibility.

    **Note** that DSSP can only handle one model.

    Examples
    --------
    >>> p = PDBParser()
    >>> structure = p.get_structure("1MOT", "1mot.pdb")
    >>> model = structure[0]
    >>> dssp = DSSP(model, "1mot.pdb")
    >>> # DSSP data is accessed by a tuple (chain_id, res_id)
    >>> a_key = list(dssp.keys())[2]
    >>> # (dssp index, amino acid, secondary structure, relative ASA, phi, psi,
    >>> # NH_O_1_relidx, NH_O_1_energy, O_NH_1_relidx, O_NH_1_energy,
    >>> # NH_O_2_relidx, NH_O_2_energy, O_NH_2_relidx, O_NH_2_energy)
    >>> dssp[a_key]
    (3, 'A', 'H', 0.7075471698113207, -61.2, -42.4,
     -2, -0.7, 4, -3.0, 1, -0.2, 5, -0.2)

    """

    def __init__(self, model, in_file, dssp="dssp", acc_array="Sander", file_type='PDB'):
        """Create a DSSP object.

        Parameters
        ----------
        model : Model
            The first model of the structure
        in_file : string
            Either a PDB file or a DSSP file.
        dssp : string
            The dssp executable (ie. the argument to os.system)
        acc_array : string
            Accessible surface area (ASA) from either Miller et al. (1987),
            Sander & Rost (1994), or Wilke: Tien et al. 2013, as string
            Sander/Wilke/Miller. Defaults to Sander.
        file_type: string
            File type switch, either PDB or DSSP with PDB as default.

        """
        self.residue_max_acc = residue_max_acc[acc_array]

        # create DSSP dictionary
        file_type = file_type.upper()
        assert(file_type in ['PDB', 'DSSP'])
        # If the input file is a PDB file run DSSP and parse output:
        if file_type == 'PDB':
            # Newer versions of DSSP program call the binary 'mkdssp', so
            # calling 'dssp' will not work in some operating systems
            # (Debian distribution of DSSP includes a symlink for 'dssp' argument)
            try:
                dssp_dict, dssp_keys = dssp_dict_from_pdb_file(in_file, dssp)
            except OSError:  # TODO: Use FileNotFoundError once drop Python 2
                if dssp == 'dssp':
                    dssp = 'mkdssp'
                elif dssp == 'mkdssp':
                    dssp = 'dssp'
                else:
                    raise
            dssp_dict, dssp_keys = dssp_dict_from_pdb_file(in_file, dssp)
        # If the input file is a DSSP file just parse it directly:
        elif file_type == 'DSSP':
            dssp_dict, dssp_keys = make_dssp_dict(in_file)

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

            # For disordered residues of point mutations, Biopython uses the
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

            (aa, ss, acc, phi, psi, dssp_index,
                NH_O_1_relidx, NH_O_1_energy,
                O_NH_1_relidx, O_NH_1_energy,
                NH_O_2_relidx, NH_O_2_energy,
                O_NH_2_relidx, O_NH_2_energy) = dssp_dict[key]

            res.xtra["SS_DSSP"] = ss
            res.xtra["EXP_DSSP_ASA"] = acc
            res.xtra["PHI_DSSP"] = phi
            res.xtra["PSI_DSSP"] = psi
            res.xtra["DSSP_INDEX"] = dssp_index
            res.xtra["NH_O_1_RELIDX_DSSP"] = NH_O_1_relidx
            res.xtra["NH_O_1_ENERGY_DSSP"] = NH_O_1_energy
            res.xtra["O_NH_1_RELIDX_DSSP"] = O_NH_1_relidx
            res.xtra["O_NH_1_ENERGY_DSSP"] = O_NH_1_energy
            res.xtra["NH_O_2_RELIDX_DSSP"] = NH_O_2_relidx
            res.xtra["NH_O_2_ENERGY_DSSP"] = NH_O_2_energy
            res.xtra["O_NH_2_RELIDX_DSSP"] = O_NH_2_relidx
            res.xtra["O_NH_2_ENERGY_DSSP"] = O_NH_2_energy

            # Relative accessibility
            resname = res.get_resname()
            try:
                rel_acc = acc / self.residue_max_acc[resname]
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
            try:
                resname = three_to_one(resname)
            except KeyError:
                resname = 'X'
            if resname == "C":
                # DSSP renames C in C-bridges to a,b,c,d,...
                # - we rename it back to 'C'
                if _dssp_cys.match(aa):
                    aa = 'C'
            # Take care of HETATM again
            if (resname != aa) and (res.id[0] == ' ' or aa != 'X'):
                raise PDBException("Structure/DSSP mismatch at %s" % res)

            dssp_vals = (dssp_index, aa, ss, rel_acc, phi, psi,
                         NH_O_1_relidx, NH_O_1_energy,
                         O_NH_1_relidx, O_NH_1_energy,
                         NH_O_2_relidx, NH_O_2_energy,
                         O_NH_2_relidx, O_NH_2_energy)

            dssp_map[key] = dssp_vals
            dssp_list.append(dssp_vals)

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
