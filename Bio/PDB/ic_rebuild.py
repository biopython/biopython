# Copyright 2019 by Robert T. Miller.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Convert XYZ Structure to internal coordinates and back, test result."""

import re

from itertools import zip_longest

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install NumPy to build proteins from internal coordinates."
    )

from Bio.PDB.PDBExceptions import PDBException
from io import StringIO
from Bio.File import as_handle
from Bio.PDB.PDBIO import PDBIO

from Bio.PDB.Structure import Structure
from Bio.PDB.internal_coords import IC_Residue, IC_Chain
from Bio.PDB.PICIO import write_PIC, read_PIC, enumerate_atoms, pdb_date

# for typing
from typing import Dict, Union, Any
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue, DisorderedResidue
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.PDB.Entity import Entity


def structure_rebuild_test(entity, verbose: bool = False) -> Dict:
    """Test rebuild PDB structure from internal coordinates.

    :param entity: Biopython Structure, Model or Chain
        Structure to test
    :param verbose: bool
        print extra messages
    :returns: dict
        comparison dict from compare_residues()
    """
    sp = StringIO()
    entity.atom_to_internal_coordinates(verbose)
    write_PIC(entity, sp)
    sp.seek(0)
    pdb2 = read_PIC(sp)
    if verbose:
        report_IC(pdb2, verbose=True)
    pdb2.internal_to_atom_coordinates(verbose)
    r = compare_residues(entity, pdb2, verbose=verbose)
    return r


def report_IC(
    entity: Union[Structure, Model, Chain, Residue],
    reportDict: Dict[str, Any] = None,
    verbose: bool = False,
) -> Dict[str, Any]:
    """Generate dict with counts of ic data elements for each entity level.

    reportDict entries are:
        - idcode : PDB ID
        - hdr : PDB header lines
        - mdl : models
        - chn : chains
        - res : residue objects
        - res_e : residues with dihedra and/or hedra
        - dih : dihedra
        - hed : hedra

    :param Entity entity: Biopython PDB Entity object: S, M, C or R
    :raises PDBException: if entity level not S, M, C, or R
    :raises Exception: if entity does not have .level attribute
    :returns: dict with counts of IC data elements
    """
    if reportDict is None:
        reportDict = {
            "idcode": None,
            "hdr": 0,
            "mdl": 0,
            "chn": 0,
            "chn_ids": [],
            "res": 0,
            "res_e": 0,
            "dih": 0,
            "hed": 0,
        }
    try:
        if "A" == entity.level:
            raise PDBException("No IC output at Atom level")
        elif isinstance(entity, Residue) or isinstance(
            entity, DisorderedResidue
        ):  # "R" == entity.level:
            if entity.internal_coord:
                reportDict["res"] += 1
                dlen = len(entity.internal_coord.dihedra)
                hlen = len(entity.internal_coord.hedra)
                if 0 < dlen or 0 < hlen:
                    reportDict["res_e"] += 1
                    reportDict["dih"] += dlen
                    reportDict["hed"] += hlen

        elif isinstance(entity, Chain):  # "C" == entity.level:
            reportDict["chn"] += 1
            reportDict["chn_ids"].append(entity.id)
            for res in entity:
                reportDict = report_IC(res, reportDict)

        elif isinstance(entity, Model):  # "M" == entity.level:
            reportDict["mdl"] += 1
            for chn in entity:
                reportDict = report_IC(chn, reportDict)

        elif isinstance(entity, Structure):  # "S" == entity.level:
            if hasattr(entity, "header"):
                if reportDict["idcode"] is None:
                    reportDict["idcode"] = entity.header.get("idcode", None)

                hdr = entity.header.get("head", None)
                if hdr:
                    reportDict["hdr"] += 1
                nam = entity.header.get("name", None)
                if nam:
                    reportDict["hdr"] += 1
            for mdl in entity:
                reportDict = report_IC(mdl, reportDict)
        else:
            raise PDBException("Cannot identify level: " + str(entity.level))
    except KeyError:
        raise Exception(
            "write_PIC: argument is not a Biopython PDB Entity " + str(entity)
        )

    if verbose:
        print(
            "{} : {} models {} chains {} {} residue objects "
            "{} residues with {} dihedra {} hedra".format(
                reportDict["idcode"],
                reportDict["mdl"],
                reportDict["chn"],
                reportDict["chn_ids"],
                reportDict["res"],
                reportDict["res_e"],
                reportDict["dih"],
                reportDict["hed"],
            )
        )

    return reportDict


def IC_duplicate(entity) -> Structure:
    """Duplicate structure entity with IC data, no atom coordinates.

    Employs write_PIC(), read_PIC() with StringIO buffer.
    Calls atom_to_internal_coordinates() if needed.

    :param entity: Biopython PDB Entity (will fail for Atom)
    :returns: Biopython PDBStructure, no Atom objects
    """
    sp = StringIO()
    hasInternalCoords = False
    for res in entity.get_residues():
        if res.internal_coord:
            if len(res.internal_coord.hedra) > 0:
                hasInternalCoords = True
                break
    if not hasInternalCoords:
        if isinstance(entity, Residue):  # "R" == entity.level:
            # works better at chain level but leave option here
            if not res.internal_coord:
                res.internal_coord = IC_Residue(entity)
            res.internal_coord.atom_to_internal_coordinates()
        else:
            entity.atom_to_internal_coordinates()

    write_PIC(entity, sp)
    sp.seek(0)
    return read_PIC(sp)


def _cmp_atm(
    r0: Residue, r1: Residue, a0: Atom, a1: Atom, verbose: bool, cmpdict: Dict
) -> None:
    cmpdict["aCount"] += 1
    if a0 is None:
        if verbose:
            print(r1.get_full_id(), "None !=", a1.get_full_id(), a1.parent.resname)
    elif a1 is None:
        if verbose:
            print(r0.get_full_id(), a0.get_full_id(), a0.parent.resname, "!= None")
    else:
        if a0.get_full_id() == a1.get_full_id():
            cmpdict["aFullIdMatchCount"] += 1
        elif verbose:
            print(
                r0.get_full_id(),
                a0.get_full_id(),
                a0.parent.resname,
                "!=",
                a1.get_full_id(),
            )
        a0c = a0.get_coord()
        a1c = a1.get_coord()
        if numpy.allclose(a0c, a1c, rtol=1e-05, atol=1e-08):
            cmpdict["aCoordMatchCount"] += 1
        elif verbose:
            print(
                "atom coords disagree:",
                r0.get_full_id(),
                a0.get_full_id(),
                a1.get_full_id(),
                a0c,
                "!=",
                a1c,
            )


def _cmp_res(r0: Residue, r1: Residue, verbose: bool, cmpdict: Dict) -> None:
    r0id, r0fid, r1fid = r0.id, r0.full_id, r1.full_id
    chn = r0.parent.id
    if chn not in cmpdict["chains"]:
        cmpdict["chains"].append(chn)
    cmpdict["rCount"] += 1
    if r0fid == r1fid:
        cmpdict["rMatchCount"] += 1
    elif verbose:
        print(r0fid, "!=", r1fid)
    if " " == r0id[0] and not (" " == r0.resname[0] or 2 == len(r0.resname)):
        # skip water, DNA (' ' == [0] for pdb, 2 == len() for mmcif)
        cmpdict["residues"] += 1
        longer = r0 if len(r0.child_dict) >= len(r1.child_dict) else r1
        for ak in longer.child_dict:
            a0 = r0.child_dict.get(ak, None)
            if a0 is None:
                aknd = re.sub("D", "H", ak, count=1)
                a0 = r0.child_dict.get(aknd, None)
            a1 = r1.child_dict.get(ak, None)
            if a1 is None:
                aknd = re.sub("D", "H", ak, count=1)
                a1 = r1.child_dict.get(aknd, None)
            if (
                a0 is None
                or a1 is None
                or 0 == a0.is_disordered() == a1.is_disordered()
            ):
                _cmp_atm(r0, r1, a0, a1, verbose, cmpdict)
            elif 2 == a0.is_disordered() == a1.is_disordered():
                cmpdict["disAtmCount"] += 1
                for da0k in a0.child_dict:
                    _cmp_atm(
                        r0,
                        r1,
                        a0.child_dict.get(da0k, None),
                        a1.child_dict.get(da0k, None),
                        verbose,
                        cmpdict,
                    )
            else:
                if verbose:
                    print("disorder disagreement:", r0.get_full_id(), ak)
                cmpdict["aCount"] += 1


def compare_residues(
    e0: Union[Structure, Model, Chain],
    e1: Union[Structure, Model, Chain],
    verbose: bool = False,
) -> Dict[str, Any]:
    """Compare full IDs and atom coordinates for 2 Biopython PDB entities.

    Skip DNA and HETATMs.

    :param e0, e1: Biopython PDB Entity objects (S, M or C)
        Structures, Models or Chains to be compared
    :param verbose: Bool
        whether to print mismatch info, default False
    :returns: Dictionary
        Result counts for Residues, Full ID match Residues, Atoms,
        Full ID match atoms, and Coordinate match atoms; report string;
        error status (bool)
    """
    cmpdict: Dict[str, Any] = {}
    cmpdict["chains"] = []
    cmpdict["residues"] = 0
    cmpdict["rCount"] = 0
    cmpdict["rMatchCount"] = 0
    cmpdict["aCount"] = 0
    cmpdict["disAtmCount"] = 0
    cmpdict["aCoordMatchCount"] = 0
    cmpdict["aFullIdMatchCount"] = 0
    cmpdict["id0"] = e0.get_full_id()
    cmpdict["id1"] = e1.get_full_id()
    cmpdict["pass"] = None
    cmpdict["report"] = None

    for r0, r1 in zip_longest(e0.get_residues(), e1.get_residues()):
        if 2 == r0.is_disordered() == r1.is_disordered():
            for dr0, dr1 in zip_longest(r0.child_dict.values(), r1.child_dict.values()):
                _cmp_res(dr0, dr1, verbose, cmpdict)
        else:
            _cmp_res(r0, r1, verbose, cmpdict)

    if (
        cmpdict["rMatchCount"] == cmpdict["rCount"]
        and cmpdict["aCoordMatchCount"] == cmpdict["aCount"]
        and cmpdict["aFullIdMatchCount"] == cmpdict["aCount"]
    ):
        cmpdict["pass"] = True
    else:
        cmpdict["pass"] = False

    rstr = (
        "{}:{} {} -- {} of {} residue IDs match; {} residues {} atom coords, "
        "{} full IDs of {} atoms ({} disordered) match : {}".format(
            cmpdict["id0"],
            cmpdict["id1"],
            cmpdict["chains"],
            cmpdict["rMatchCount"],
            cmpdict["rCount"],
            cmpdict["residues"],
            cmpdict["aCoordMatchCount"],
            cmpdict["aFullIdMatchCount"],
            cmpdict["aCount"],
            cmpdict["disAtmCount"],
            "ERROR" if not cmpdict["pass"] else "ALL OK",
        )
    )
    if not cmpdict["pass"]:
        if cmpdict["rMatchCount"] != cmpdict["rCount"]:
            rstr += " -RESIDUE IDS-"
        if cmpdict["aCoordMatchCount"] != cmpdict["aFullIdMatchCount"]:
            rstr += " -COORDINATES-"
        if cmpdict["aFullIdMatchCount"] != cmpdict["aCount"]:
            rstr += " -ATOM IDS-"
    cmpdict["report"] = rstr

    return cmpdict


def write_PDB(
    entity: Structure, file: str, pdbid: str = None, chainid: str = None
) -> None:
    """Write PDB file with HEADER and TITLE."""
    enumerate_atoms(entity)
    with as_handle(file, "w") as fp:
        try:
            if "S" == entity.level:
                if hasattr(entity, "header"):
                    if not pdbid:
                        pdbid = entity.header.get("idcode", None)
                    hdr = entity.header.get("head", None)
                    dd = pdb_date(entity.header.get("deposition_date", None))

                    if hdr:
                        fp.write(
                            ("HEADER    {:40}{:8}   {:4}\n").format(
                                hdr.upper(), (dd or ""), (pdbid or "")
                            )
                        )
                    nam = entity.header.get("name", None)
                    if nam:
                        fp.write("TITLE     " + nam.upper() + "\n")
                io = PDBIO()
                io.set_structure(entity)
                io.save(fp, preserve_atom_numbering=True)

            else:
                raise PDBException("level not 'S': " + str(entity.level))
        except KeyError:
            raise Exception(
                "write_PIC: argument is not a Biopython PDB Entity " + str(entity)
            )
