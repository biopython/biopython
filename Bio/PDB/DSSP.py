# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

r"""Use the DSSP program to calculate secondary structure and accessibility.

You need to have a working version of DSSP (and a license, free for academic
use) in order to use this. For DSSP, see https://swift.cmbi.umcn.nl/gv/dssp/.

The following Accessible surface area (ASA) values can be used, defaulting
to the Sander and Rost values:

    Ahmad
        Ahmad et al. 2003 https://doi.org/10.1002/prot.10328
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
The DSSP class can be used to run DSSP on a PDB or mmCIF file, and provides a
handle to the DSSP secondary structure and accessibility.

**Note** that DSSP can only handle one model, and will only run
calculations on the first model in the provided PDB file.

Examples
--------
Typical use::

    from Bio.PDB import PDBParser
    from Bio.PDB.DSSP import DSSP
    p = PDBParser()
    structure = p.get_structure("1MOT", "/local-pdb/1mot.pdb")
    model = structure[0]
    dssp = DSSP(model, "/local-pdb/1mot.pdb")

Note that the recent DSSP executable from the DSSP-2 package was
renamed from ``dssp`` to ``mkdssp``. If using a recent DSSP release,
you may need to provide the name of your DSSP executable::

    dssp = DSSP(model, '/local-pdb/1mot.pdb', dssp='mkdssp')

DSSP data is accessed by a tuple - (chain id, residue id)::

    a_key = list(dssp.keys())[2]
    dssp[a_key]

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


import re
import os
from io import StringIO
import subprocess
import warnings

from Bio.PDB.AbstractPropertyMap import AbstractResiduePropertyMap
from Bio.PDB.PDBExceptions import PDBException
from Bio.PDB.PDBParser import PDBParser
from Bio.Data.PDBData import protein_letters_3to1, residue_sasa_scales
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

# Match C in DSSP
_dssp_cys = re.compile("[a-z]")

# Maximal ASA of amino acids
# Used for relative accessibility
residue_max_acc = residue_sasa_scales


def version(version_string):
    """Parse semantic version scheme for easy comparison."""
    return tuple(map(int, (version_string.split("."))))


def ss_to_index(ss):
    """Secondary structure symbol to index.

    H=0
    E=1
    C=2
    """
    if ss == "H":
        return 0
    if ss == "E":
        return 1
    if ss == "C":
        return 2
    assert 0


def dssp_dict_from_pdb_file(in_file, DSSP="dssp", dssp_version="3.9.9"):
    """Create a DSSP dictionary from a PDB file.

    Parameters
    ----------
    in_file : string
        pdb file

    DSSP : string
        DSSP executable (argument to subprocess)

    dssp_version : string
        Version of DSSP excutable

    Returns
    -------
    (out_dict, keys) : tuple
        a dictionary that maps (chainid, resid) to
        amino acid type, secondary structure code and
        accessibility.

    Examples
    --------
    How dssp_dict_from_pdb_file could be used::

        from Bio.PDB.DSSP import dssp_dict_from_pdb_file
        dssp_tuple = dssp_dict_from_pdb_file("/local-pdb/1fat.pdb")
        dssp_dict = dssp_tuple[0]
        print(dssp_dict['A',(' ', 1, ' ')])

    """
    # Using universal newlines is important on Python 3, this
    # gives text handles rather than bytes handles.
    # Newer version of DSSP executable is named 'mkdssp',
    # and calling 'dssp' will hence not work in some operating systems
    # (Debian distribution of DSSP includes a symlink for 'dssp' argument)
    try:
        if version(dssp_version) < version("4.0.0"):
            DSSP_cmd = [DSSP, in_file]
        else:
            DSSP_cmd = [DSSP, "--output-format=dssp", in_file]
        p = subprocess.Popen(
            DSSP_cmd,
            universal_newlines=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except FileNotFoundError:
        if DSSP == "mkdssp":
            raise
        if version(dssp_version) < version("4.0.0"):
            DSSP_cmd = ["mkdssp", in_file]
        else:
            DSSP_cmd = ["mkdssp", "--output-format=dssp", in_file]
        p = subprocess.Popen(
            DSSP_cmd,
            universal_newlines=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    out, err = p.communicate()

    # Alert user for errors
    if err.strip():
        warnings.warn(err)
        if not out.strip():
            raise Exception("DSSP failed to produce an output")

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
    with open(filename) as handle:
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
    for line in handle:
        sl = line.split()
        if len(sl) < 2:
            continue
        if sl[1] == "RESIDUE":
            # Start parsing from here
            start = 1
            continue
        if not start:
            continue
        if line[9] == " ":
            # Skip -- missing residue
            continue

        dssp_index = int(line[:5])
        resseq = int(line[5:10])
        icode = line[10]
        chainid = line[11]
        aa = line[13]
        ss = line[16]
        if ss == " ":
            ss = "-"
        try:
            NH_O_1_relidx = int(line[38:45])
            NH_O_1_energy = float(line[46:50])
            O_NH_1_relidx = int(line[50:56])
            O_NH_1_energy = float(line[57:61])
            NH_O_2_relidx = int(line[61:67])
            NH_O_2_energy = float(line[68:72])
            O_NH_2_relidx = int(line[72:78])
            O_NH_2_energy = float(line[79:83])

            acc = int(line[34:38])
            phi = float(line[103:109])
            psi = float(line[109:115])
        except ValueError as exc:
            # DSSP output breaks its own format when there are >9999
            # residues, since only 4 digits are allocated to the seq num
            # field.  See 3kic chain T res 321, 1vsy chain T res 6077.
            # Here, look for whitespace to figure out the number of extra
            # digits, and shift parsing the rest of the line by that amount.
            if line[34] != " ":
                shift = line[34:].find(" ")

                NH_O_1_relidx = int(line[38 + shift : 45 + shift])
                NH_O_1_energy = float(line[46 + shift : 50 + shift])
                O_NH_1_relidx = int(line[50 + shift : 56 + shift])
                O_NH_1_energy = float(line[57 + shift : 61 + shift])
                NH_O_2_relidx = int(line[61 + shift : 67 + shift])
                NH_O_2_energy = float(line[68 + shift : 72 + shift])
                O_NH_2_relidx = int(line[72 + shift : 78 + shift])
                O_NH_2_energy = float(line[79 + shift : 83 + shift])

                acc = int(line[34 + shift : 38 + shift])
                phi = float(line[103 + shift : 109 + shift])
                psi = float(line[109 + shift : 115 + shift])
            else:
                raise ValueError(exc) from None
        res_id = (" ", resseq, icode)
        dssp[(chainid, res_id)] = (
            aa,
            ss,
            acc,
            phi,
            psi,
            dssp_index,
            NH_O_1_relidx,
            NH_O_1_energy,
            O_NH_1_relidx,
            O_NH_1_energy,
            NH_O_2_relidx,
            NH_O_2_energy,
            O_NH_2_relidx,
            O_NH_2_energy,
        )
        keys.append((chainid, res_id))
    return dssp, keys


class DSSP(AbstractResiduePropertyMap):
    """Run DSSP and parse secondary structure and accessibility.

    Run DSSP on a PDB/mmCIF file, and provide a handle to the
    DSSP secondary structure and accessibility.

    **Note** that DSSP can only handle one model.

    Examples
    --------
    How DSSP could be used::

        from Bio.PDB import PDBParser
        from Bio.PDB.DSSP import DSSP
        p = PDBParser()
        structure = p.get_structure("1MOT", "/local-pdb/1mot.pdb")
        model = structure[0]
        dssp = DSSP(model, "/local-pdb/1mot.pdb")
        # DSSP data is accessed by a tuple (chain_id, res_id)
        a_key = list(dssp.keys())[2]
        # (dssp index, amino acid, secondary structure, relative ASA, phi, psi,
        # NH_O_1_relidx, NH_O_1_energy, O_NH_1_relidx, O_NH_1_energy,
        # NH_O_2_relidx, NH_O_2_energy, O_NH_2_relidx, O_NH_2_energy)
        dssp[a_key]

    """

    def __init__(self, model, in_file, dssp="dssp", acc_array="Sander", file_type=""):
        """Create a DSSP object.

        Parameters
        ----------
        model : Model
            The first model of the structure
        in_file : string
            Either a PDB file or a DSSP file.
        dssp : string
            The dssp executable (ie. the argument to subprocess)
        acc_array : string
            Accessible surface area (ASA) from either Miller et al. (1987),
            Sander & Rost (1994), Wilke: Tien et al. 2013, or Ahmad et al.
            (2003) as string Sander/Wilke/Miller/Ahmad. Defaults to Sander.
        file_type: string
            File type switch: either PDB, MMCIF or DSSP. Inferred from the
            file extension by default.

        """
        self.residue_max_acc = residue_max_acc[acc_array]

        # create DSSP dictionary
        if file_type == "":
            file_type = os.path.splitext(in_file)[1][1:]
        file_type = file_type.upper()
        if file_type == "CIF":
            file_type = "MMCIF"
        assert file_type in [
            "PDB",
            "MMCIF",
            "DSSP",
        ], "File type must be PDB, mmCIF or DSSP"
        # If the input file is a PDB or mmCIF file run DSSP and parse output:
        if file_type == "PDB" or file_type == "MMCIF":
            # Newer versions of DSSP program call the binary 'mkdssp', so
            # calling 'dssp' will not work in some operating systems
            # (Debian distribution of DSSP includes a symlink for 'dssp' argument)
            try:
                version_string = subprocess.check_output(
                    [dssp, "--version"], universal_newlines=True
                )
                dssp_version = re.search(r"\s*([\d.]+)", version_string).group(1)
                dssp_dict, dssp_keys = dssp_dict_from_pdb_file(
                    in_file, dssp, dssp_version
                )
            except FileNotFoundError:
                if dssp == "dssp":
                    dssp = "mkdssp"
                elif dssp == "mkdssp":
                    dssp = "dssp"
                else:
                    raise
                version_string = subprocess.check_output(
                    [dssp, "--version"], universal_newlines=True
                )
                dssp_version = re.search(r"\s*([\d.]+)", version_string).group(1)
                dssp_dict, dssp_keys = dssp_dict_from_pdb_file(
                    in_file, dssp, dssp_version
                )
        # If the input file is a DSSP file just parse it directly:
        elif file_type == "DSSP":
            dssp_dict, dssp_keys = make_dssp_dict(in_file)

        dssp_map = {}
        dssp_list = []

        def resid2code(res_id):
            """Serialize a residue's resseq and icode for easy comparison."""
            return f"{res_id[1]}{res_id[2]}"

        # DSSP outputs label_asym_id from the mmCIF file as the chain ID
        # But MMCIFParser reads in the auth_asym_id
        # Here we create a dictionary to map label_asym_id to auth_asym_id
        # using the mmCIF file
        if file_type == "MMCIF" and version(dssp_version) < version("4.0.0"):
            mmcif_dict = MMCIF2Dict(in_file)
            mmcif_chain_dict = {}
            for i, c in enumerate(mmcif_dict["_atom_site.label_asym_id"]):
                if c not in mmcif_chain_dict:
                    mmcif_chain_dict[c] = mmcif_dict["_atom_site.auth_asym_id"][i]
            dssp_mapped_keys = []

        # Now create a dictionary that maps Residue objects to
        # secondary structure and accessibility, and a list of
        # (residue, (secondary structure, accessibility)) tuples
        for key in dssp_keys:
            chain_id, res_id = key
            if file_type == "MMCIF" and version(dssp_version) < version("4.0.0"):
                chain_id = mmcif_chain_dict[chain_id]
                dssp_mapped_keys.append((chain_id, res_id))
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
                    if r.id[0] not in (" ", "W"):
                        # Compare resseq + icode
                        if resid2code(r.id) == res_seq_icode:
                            # Found a matching residue
                            res = r
                            break
                else:
                    raise KeyError(res_id) from None

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
                    if altloc in tuple("A1 "):
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
                altlocs = {a.get_altloc() for a in res.get_unpacked_list()}
                if altlocs.isdisjoint("A1 "):
                    # Try again with all HETATM other than water
                    res_seq_icode = resid2code(res_id)
                    for r in chain:
                        if r.id[0] not in (" ", "W"):
                            if resid2code(r.id) == res_seq_icode and r.get_list()[
                                0
                            ].get_altloc() in tuple("A1 "):
                                res = r
                                break

            (
                aa,
                ss,
                acc,
                phi,
                psi,
                dssp_index,
                NH_O_1_relidx,
                NH_O_1_energy,
                O_NH_1_relidx,
                O_NH_1_energy,
                NH_O_2_relidx,
                NH_O_2_energy,
                O_NH_2_relidx,
                O_NH_2_energy,
            ) = dssp_dict[key]

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
                rel_acc = "NA"
            else:
                if rel_acc > 1.0:
                    rel_acc = 1.0
            res.xtra["EXP_DSSP_RASA"] = rel_acc
            # Verify if AA in DSSP == AA in Structure
            # Something went wrong if this is not true!
            # NB: DSSP uses X often
            resname = protein_letters_3to1.get(resname, "X")
            if resname == "C":
                # DSSP renames C in C-bridges to a,b,c,d,...
                # - we rename it back to 'C'
                if _dssp_cys.match(aa):
                    aa = "C"
            # Take care of HETATM again
            if (resname != aa) and (res.id[0] == " " or aa != "X"):
                raise PDBException(f"Structure/DSSP mismatch at {res}")

            dssp_vals = (
                dssp_index,
                aa,
                ss,
                rel_acc,
                phi,
                psi,
                NH_O_1_relidx,
                NH_O_1_energy,
                O_NH_1_relidx,
                O_NH_1_energy,
                NH_O_2_relidx,
                NH_O_2_energy,
                O_NH_2_relidx,
                O_NH_2_energy,
            )

            dssp_map[(chain_id, res_id)] = dssp_vals
            dssp_list.append(dssp_vals)

        if file_type == "MMCIF" and version(dssp_version) < version("4.0.0"):
            dssp_keys = dssp_mapped_keys
        AbstractResiduePropertyMap.__init__(self, dssp_map, dssp_keys, dssp_list)
