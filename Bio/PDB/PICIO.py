# Copyright 2019 by Robert T. Miller.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""PICIO: read and write Protein Internal Coordinate (.pic) data files."""

import re

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install NumPy to build proteins from internal coordinates."
    )

from Bio.File import as_handle
from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.PDB.parse_pdb_header import _parse_pdb_header_list
from Bio.PDB.PDBExceptions import PDBException

from Bio.PDB.internal_coords import IC_Residue, IC_Chain, Edron, AtomKey

from typing import Dict, TextIO
from Bio.PDB.Structure import Structure


def read_PIC(file: TextIO, verbose: bool = False) -> Structure:
    """Load Protein Internal Coordinate (.pic) data from file.

    PIC file format:
        - comment lines start with #
        - (optional) PDB HEADER record
           - idcode and deposition date recommended but optional
           - deposition date in PDB format or as changed by Biopython
        - (optional) PDB TITLE record
        - repeat:
           - Biopython Residue Full ID - sets residue IDs of returned structure
           - (optional) PDB N, CA, C ATOM records for chain start
           - (optional) PIC Hedra records for residue
           - (optional) PIC Dihedra records for residue
           - (optional) BFAC records listing AtomKeys and b-factors

    An improvement would define relative positions for HOH (water) entries.

    N.B. dihedron (i-1)C-N-CA-CB is ignored in assembly if O exists.

    C-beta is by default placed using O-C-CA-CB, but O is missing
    in some PDB file residues, which means the sidechain cannot be
    placed.  The alternate CB path (i-1)C-N-CA-CB is provided to
    circumvent this, but if this is needed then it must be adjusted in
    conjunction with PHI ((i-1)C-N-CA-C) as they overlap.  (i-1)C-N-CA-CB is
    included by default in .pic files for consistency and informational
    (e.g. statistics gathering) purposes, as otherwise the dihedron would only
    appear in the few cases it is needed for.

    :param Bio.File file: file name or handle
    :param bool verbose: complain when lines not as expected
    :returns: Biopython Structure object, Residues with .internal_coord attributes
        but no coordinates except for chain start N, CA, C atoms if supplied,
        **OR** None on parse fail (silent unless verbose=True)

    """
    pdb_hdr_re = re.compile(
        r"^HEADER\s{4}(?P<cf>.{1,40})"
        r"(?:\s+(?P<dd>\d\d\d\d-\d\d-\d\d|\d\d-\w\w\w-\d\d))?"
        r"(?:\s+(?P<id>[0-9A-Z]{4}))?\s*$"
    )
    # ^\('(?P<pid>\w*)',\s(?P<mdl>\d+),\s'(?P<chn>\w)',\s\('(?P<het>\s|[\w-]+)',\s(?P<pos>\d+),\s'(?P<icode>\s|\w)'\)\)\s(?P<res>[A-Z]{3})\s(\[(?P<segid>[a-zA-z\s]{4})\])?\s*$
    pdb_ttl_re = re.compile(r"^TITLE\s{5}(?P<ttl>.+)\s*$")
    biop_id_re = re.compile(
        r"^\('(?P<pid>[^\s]*)',\s(?P<mdl>\d+),\s"
        r"'(?P<chn>\s|\w)',\s\('(?P<het>\s|[\w\s-]+)"
        r"',\s(?P<pos>-?\d+),\s'(?P<icode>\s|\w)'\)\)"
        r"\s+(?P<res>[\w]{1,3})"
        r"(\s\[(?P<segid>[a-zA-z\s]+)\])?"
        r"\s*$"
    )
    pdb_atm_re = re.compile(
        r"^ATOM\s\s(?:\s*(?P<ser>\d+))\s(?P<atm>[\w\s]{4})"
        r"(?P<alc>\w|\s)(?P<res>[\w]{3})\s(?P<chn>.)"
        r"(?P<pos>[\s\-\d]{4})(?P<icode>[A-Za-z\s])\s\s\s"
        r"(?P<x>[\s\-\d\.]{8})(?P<y>[\s\-\d\.]{8})"
        r"(?P<z>[\s\-\d\.]{8})(?P<occ>[\s\d\.]{6})"
        r"(?P<tfac>[\s\d\.]{6})\s{6}"
        r"(?P<segid>[a-zA-z\s]{4})(?P<elm>.{2})"
        r"(?P<chg>.{2})?\s*$"
    )
    bfac_re = re.compile(
        r"^BFAC:\s([^\s]+\s+[\-\d\.]+)"
        r"\s*([^\s]+\s+[\-\d\.]+)?"
        r"\s*([^\s]+\s+[\-\d\.]+)?"
        r"\s*([^\s]+\s+[\-\d\.]+)?"
        r"\s*([^\s]+\s+[\-\d\.]+)?"
    )
    bfac2_re = re.compile(r"([^\s]+)\s+([\-\d\.]+)")
    struct_builder = StructureBuilder()

    # init empty header dict
    # - could use to parse HEADER and TITLE lines except
    #   deposition_date format changed from original PDB header
    header_dict = _parse_pdb_header_list([])

    curr_SMCS = [None, None, None, None]  # struct model chain seg
    SMCS_init = [
        struct_builder.init_structure,
        struct_builder.init_model,
        struct_builder.init_chain,
        struct_builder.init_seg,
    ]

    sb_res = None

    with as_handle(file, mode="r") as handle:
        for aline in handle.readlines():
            if aline.startswith("#"):
                pass  # skip comment lines
            elif aline.startswith("HEADER "):
                m = pdb_hdr_re.match(aline)
                if m:
                    header_dict["head"] = m.group("cf")  # classification
                    header_dict["idcode"] = m.group("id")
                    header_dict["deposition_date"] = m.group("dd")
                elif verbose:
                    print("Reading pic file", file, "HEADER parse fail: ", aline)
            elif aline.startswith("TITLE "):
                m = pdb_ttl_re.match(aline)
                if m:
                    header_dict["name"] = m.group("ttl").strip()
                    # print('TTL: ', m.group('ttl').strip())
                elif verbose:
                    print("Reading pic file", file, "TITLE parse fail:, ", aline)
            elif aline.startswith("("):  # Biopython ID line for Residue
                m = biop_id_re.match(aline)
                if m:
                    # check SMCS = Structure, Model, Chain, SegID
                    segid = m.group(9)
                    if segid is None:
                        segid = "    "
                    this_SMCS = [m.group(1), int(m.group(2)), m.group(3), segid]
                    if curr_SMCS != this_SMCS:
                        # init new SMCS level as needed
                        for i in range(4):
                            if curr_SMCS[i] != this_SMCS[i]:
                                SMCS_init[i](this_SMCS[i])
                                curr_SMCS[i] = this_SMCS[i]
                                if 0 == i:
                                    # 0 = init structure so add header
                                    struct_builder.set_header(header_dict)
                                elif 1 == i:
                                    # new model means new chain and new segid
                                    curr_SMCS[2] = curr_SMCS[3] = None

                    struct_builder.init_residue(
                        m.group("res"),
                        m.group("het"),
                        int(m.group("pos")),
                        m.group("icode"),
                    )

                    sb_res = struct_builder.residue
                    if 2 == sb_res.is_disordered():
                        for r in sb_res.child_dict.values():
                            if not r.internal_coord:
                                sb_res = r
                                break
                    sb_res.internal_coord = IC_Residue(sb_res)
                    # print('res id:', m.groupdict())
                    # print(report_IC(struct_builder.get_structure()))
                else:
                    if verbose:
                        print(
                            "Reading pic file", file, "residue ID parse fail: ", aline
                        )
                    return None
            elif aline.startswith("ATOM "):
                m = pdb_atm_re.match(aline)
                if m:
                    if sb_res is None:
                        # ATOM without res spec already loaded, not a pic file
                        if verbose:
                            print(
                                "Reading pic file",
                                file,
                                "ATOM without residue configured:, ",
                                aline,
                            )
                        return None
                    if sb_res.resname != m.group("res") or sb_res.id[1] != int(
                        m.group("pos")
                    ):
                        if verbose:
                            print(
                                "Reading pic file",
                                file,
                                "ATOM not in configured residue (",
                                sb_res.resname,
                                str(sb_res.id),
                                "):",
                                aline,
                            )
                        return None
                    coord = numpy.array(
                        (float(m.group("x")), float(m.group("y")), float(m.group("z"))),
                        "f",
                    )
                    struct_builder.init_atom(
                        m.group("atm").strip(),
                        coord,
                        float(m.group("tfac")),
                        float(m.group("occ")),
                        m.group("alc"),
                        m.group("atm"),
                        int(m.group("ser")),
                        m.group("elm").strip(),
                    )

                    # print('atom: ', m.groupdict())
                # elif verbose:
                #     print("Reading pic file", file, "ATOM parse fail:", aline)
            elif aline.startswith("BFAC: "):
                m = bfac_re.match(aline)
                if m:
                    for bfac_pair in m.groups():
                        if bfac_pair is not None:
                            m2 = bfac2_re.match(bfac_pair)
                            if m2 and sb_res is not None and sb_res.internal_coord:
                                rp = sb_res.internal_coord
                                rp.bfactors[m2.group(1)] = float(m2.group(2))
                # else:
                #    print('Reading pic file', file, 'B-factor line fail: ', aline)
            else:
                m = Edron.edron_re.match(aline)
                if m and sb_res is not None:
                    sb_res.internal_coord.load_PIC(m.groupdict())
                elif m:
                    print(
                        "PIC file: ",
                        file,
                        " error: no residue info before reading (di/h)edron data: ",
                        aline,
                    )
                    return None
                elif aline.strip():
                    if verbose:
                        print("Reading PIC file", file, "parse fail on: .", aline, ".")
                    return None

    struct = struct_builder.get_structure()
    for chn in struct.get_chains():
        chnp = chn.internal_coord = IC_Chain(chn)
        # done in IC_Chain init : chnp.set_residues()
        chnp.link_residues()
        chnp.init_edra()

    # print(report_PIC(struct_builder.get_structure()))
    return struct


def _wpr(entity, fp, pdbid, chainid):
    if entity.internal_coord:
        if not chainid or not pdbid:
            chain = entity.parent
            if not chainid:
                chainid = chain.id
            if not pdbid:
                struct = chain.parent.parent
                pdbid = struct.header.get("idcode")

        fp.write(entity.internal_coord.write_PIC(pdbid, chainid))
    else:
        fp.write(IC_Residue._residue_string(entity))


def _enumerate_entity_atoms(entity):
    need = False
    for atm in entity.get_atoms():
        need = not atm.get_serial_number()
        break
    if need:
        anum = 1
        for res in entity.get_residues():
            if 2 == res.is_disordered():
                for r in res.child_dict.values():
                    for atm in r.get_unpacked_list():
                        atm.set_serial_number(anum)
                        anum += 1
            else:
                for atm in res.get_unpacked_list():
                    atm.set_serial_number(anum)
                    anum += 1


def enumerate_atoms(entity):
    """Ensure all atoms in entity have serial_number set."""
    while entity.get_parent():
        entity = entity.get_parent()  # get to top level

    if "S" == entity.level:
        for mdl in entity:  # each model starts with 1
            _enumerate_entity_atoms(mdl)
    else:  # only Chain or Residue, start with 1
        _enumerate_entity_atoms(entity)


def pdb_date(datestr: str) -> str:
    """Convert yyyy-mm-dd date to dd-month-yy."""
    if datestr:
        m = re.match(r"(\d{4})-(\d{2})-(\d{2})", datestr)
        if m:
            mo = [
                "XXX",
                "JAN",
                "FEB",
                "MAR",
                "APR",
                "MAY",
                "JUN",
                "JUL",
                "AUG",
                "SEP",
                "OCT",
                "NOV",
                "DEC",
            ][int(m.group(2))]
            datestr = m.group(3) + "-" + mo + "-" + m.group(1)[-2:]
    return datestr


def write_PIC(entity, file, pdbid=None, chainid=None):
    """Write Protein Internal Coordinates (PIC) to file.

    See read_PIC() for file format.  Recurses to lower entity levels (M, C, R).

    :param Entity entity: Biopython PDB Entity object: S, M, C or R
    :param Bio.File file: file name or handle
    :param str pdbid: PDB idcode, read from entity if not supplied
    :param char chainid: PDB Chain ID, set from C level entity.id if needed
    :raises PDBException: if entity level not S, M, C, or R
    :raises Exception: if entity does not have .level attribute
    """
    enumerate_atoms(entity)

    with as_handle(file, "w") as fp:
        try:
            if "A" == entity.level:
                raise PDBException("No PIC output at Atom level")
            elif "R" == entity.level:
                if 2 == entity.is_disordered():
                    for r in entity.child_dict.values():
                        _wpr(r, fp, pdbid, chainid)
                else:
                    _wpr(entity, fp, pdbid, chainid)
            elif "C" == entity.level:
                if not chainid:
                    chainid = entity.id
                for res in entity:
                    write_PIC(res, fp, pdbid, chainid)
                    pass
            elif "M" == entity.level:
                for chn in entity:
                    write_PIC(chn, fp, pdbid, chainid)
            elif "S" == entity.level:
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
                for mdl in entity:
                    write_PIC(mdl, fp, pdbid, chainid)
            else:
                raise PDBException("Cannot identify level: " + str(entity.level))
        except KeyError:
            raise Exception(
                "write_PIC: argument is not a Biopython PDB Entity " + str(entity)
            )
