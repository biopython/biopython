# Copyright 2012 by Eric Talevich.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Bio.SeqIO support for accessing sequences in PDB and mmCIF files."""

import collections
import shutil
import warnings

from Bio import BiopythonParserWarning
from Bio._py3k import StringIO
from Bio.Alphabet import generic_protein
from Bio.Data.SCOPData import protein_letters_3to1
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def AtomIterator(pdb_id, struct):
    """Return SeqRecords from Structure objects.

    Base function for sequence parsers that read structures Bio.PDB parsers.

    Once a parser from Bio.PDB has been used to load a structure into a
    Bio.PDB.Structure.Structure object, there is no difference in how the
    sequence parser interprets the residue sequence. The functions in this
    module may be used by SeqIO modules wishing to parse sequences from lists
    of residues.

    Calling funtions must pass a Bio.PDB.Structure.Structure object.


    See Bio.SeqIO.PdbIO.PdbAtomIterator and Bio.SeqIO.PdbIO.CifAtomIterator for
    details.
    """
    from Bio.SeqUtils import seq1

    def restype(residue):
        """Return a residue's type as a one-letter code.

        Non-standard residues (e.g. CSD, ANP) are returned as 'X'.
        """
        return seq1(residue.resname, custom_map=protein_letters_3to1)

    model = struct[0]
    for chn_id, chain in sorted(model.child_dict.items()):
        # HETATM mod. res. policy: remove mod if in sequence, else discard
        residues = [res for res in chain.get_unpacked_list()
                    if seq1(res.get_resname().upper(),
                            custom_map=protein_letters_3to1) != "X"]
        if not residues:
            continue
        # Identify missing residues in the structure
        # (fill the sequence with 'X' residues in these regions)
        gaps = []
        rnumbers = [r.id[1] for r in residues]
        for i, rnum in enumerate(rnumbers[:-1]):
            if rnumbers[i + 1] != rnum + 1:
                # It's a gap!
                gaps.append((i + 1, rnum, rnumbers[i + 1]))
        if gaps:
            res_out = []
            prev_idx = 0
            for i, pregap, postgap in gaps:
                if postgap > pregap:
                    gapsize = postgap - pregap - 1
                    res_out.extend(restype(x) for x in residues[prev_idx:i])
                    prev_idx = i
                    res_out.append('X' * gapsize)
                else:
                    warnings.warn("Ignoring out-of-order residues after a gap",
                                  BiopythonParserWarning)
                    # Keep the normal part, drop the out-of-order segment
                    # (presumably modified or hetatm residues, e.g. 3BEG)
                    res_out.extend(restype(x) for x in residues[prev_idx:i])
                    break
            else:
                # Last segment
                res_out.extend(restype(x) for x in residues[prev_idx:])
        else:
            # No gaps
            res_out = [restype(x) for x in residues]
        record_id = "%s:%s" % (pdb_id, chn_id)
        # ENH - model number in SeqRecord id if multiple models?
        # id = "Chain%s" % str(chain.id)
        # if len(structure) > 1 :
        #     id = ("Model%s|" % str(model.id)) + id

        record = SeqRecord(Seq(''.join(res_out), generic_protein),
                           id=record_id, description=record_id)

        record.annotations["model"] = model.id
        record.annotations["chain"] = chain.id

        record.annotations["start"] = int(rnumbers[0])
        record.annotations["end"] = int(rnumbers[-1])
        yield record


def PdbSeqresIterator(handle):
    """Return SeqRecord objects for each chain in a PDB file.

    The sequences are derived from the SEQRES lines in the
    PDB file header, not the atoms of the 3D structure.

    Specifically, these PDB records are handled: DBREF, SEQADV, SEQRES, MODRES

    See: http://www.wwpdb.org/documentation/format23/sect3.html

    This gets called internally via Bio.SeqIO for the SEQRES based interpretation
    of the PDB file format:

    >>> from Bio import SeqIO
    >>> for record in SeqIO.parse("PDB/1A8O.pdb", "pdb-seqres"):
    ...     print("Record id %s, chain %s" % (record.id, record.annotations["chain"]))
    ...     print(record.dbxrefs)
    ...
    Record id 1A8O:A, chain A
    ['UNP:P12497', 'UNP:POL_HV1N5']

    Equivalently,

    >>> with open("PDB/1A8O.pdb") as handle:
    ...     for record in PdbSeqresIterator(handle):
    ...         print("Record id %s, chain %s" % (record.id, record.annotations["chain"]))
    ...         print(record.dbxrefs)
    ...
    Record id 1A8O:A, chain A
    ['UNP:P12497', 'UNP:POL_HV1N5']

    Note the chain is recorded in the annotations dictionary, and any PDB DBREF
    lines are recorded in the database cross-references list.
    """
    # Late-binding import to avoid circular dependency on SeqIO in Bio.SeqUtils
    from Bio.SeqUtils import seq1

    chains = collections.defaultdict(list)
    metadata = collections.defaultdict(list)
    for line in handle:
        rec_name = line[0:6].strip()
        if rec_name == 'SEQRES':
            # NB: We only actually need chain ID and the residues here;
            # commented bits are placeholders from the wwPDB spec.
            # Serial number of the SEQRES record for the current chain.
            # Starts at 1 and increments by one each line.
            # Reset to 1 for each chain.
            # ser_num = int(line[8:10])
            # Chain identifier. This may be any single legal character,
            # including a blank which is used if there is only one chain.
            chn_id = line[11]
            # Number of residues in the chain (repeated on every record)
            # num_res = int(line[13:17])
            residues = [seq1(res, custom_map=protein_letters_3to1) for res in line[19:].split()]
            chains[chn_id].extend(residues)
        elif rec_name == 'DBREF':
            #  ID code of this entry (PDB ID)
            pdb_id = line[7:11]
            # Chain identifier.
            chn_id = line[12]
            # Initial sequence number of the PDB sequence segment.
            # seq_begin = int(line[14:18])
            # Initial insertion code of the PDB sequence segment.
            # icode_begin = line[18]
            # Ending sequence number of the PDB sequence segment.
            # seq_end = int(line[20:24])
            # Ending insertion code of the PDB sequence segment.
            # icode_end = line[24]
            # Sequence database name.
            database = line[26:32].strip()
            # Sequence database accession code.
            db_acc = line[33:41].strip()
            # Sequence database identification code.
            db_id_code = line[42:54].strip()
            # Initial sequence number of the database seqment.
            # db_seq_begin = int(line[55:60])
            # Insertion code of initial residue of the segment, if PDB is the
            # reference.
            # db_icode_begin = line[60]
            # Ending sequence number of the database segment.
            # db_seq_end = int(line[62:67])
            # Insertion code of the ending residue of the segment, if PDB is the
            # reference.
            # db_icode_end = line[67]
            metadata[chn_id].append({'pdb_id': pdb_id, 'database': database,
                                    'db_acc': db_acc, 'db_id_code': db_id_code})
        # ENH: 'SEQADV' 'MODRES'

    for chn_id, residues in sorted(chains.items()):
        record = SeqRecord(Seq(''.join(residues), generic_protein))
        record.annotations = {"chain": chn_id}
        if chn_id in metadata:
            m = metadata[chn_id][0]
            record.id = record.name = "%s:%s" % (m['pdb_id'], chn_id)
            record.description = ("%s:%s %s" % (m['database'],
                                                m['db_acc'],
                                                m['db_id_code']))
            for melem in metadata[chn_id]:
                record.dbxrefs.extend([
                    "%s:%s" % (melem['database'], melem['db_acc']),
                    "%s:%s" % (melem['database'], melem['db_id_code'])])
        else:
            record.id = chn_id
        yield record


def PdbAtomIterator(handle):
    """Return SeqRecord objects for each chain in a PDB file.

    The sequences are derived from the 3D structure (ATOM records), not the
    SEQRES lines in the PDB file header.

    Unrecognised three letter amino acid codes (e.g. "CSD") from HETATM entries
    are converted to "X" in the sequence.

    In addition to information from the PDB header (which is the same for all
    records), the following chain specific information is placed in the
    annotation:

    record.annotations["residues"] = List of residue ID strings
    record.annotations["chain"] = Chain ID (typically A, B ,...)
    record.annotations["model"] = Model ID (typically zero)

    Where amino acids are missing from the structure, as indicated by residue
    numbering, the sequence is filled in with 'X' characters to match the size
    of the missing region, and  None is included as the corresponding entry in
    the list record.annotations["residues"].

    This function uses the Bio.PDB module to do most of the hard work. The
    annotation information could be improved but this extra parsing should be
    done in parse_pdb_header, not this module.

    This gets called internally via Bio.SeqIO for the atom based interpretation
    of the PDB file format:

    >>> from Bio import SeqIO
    >>> for record in SeqIO.parse("PDB/1A8O.pdb", "pdb-atom"):
    ...     print("Record id %s, chain %s" % (record.id, record.annotations["chain"]))
    ...
    Record id 1A8O:A, chain A

    Equivalently,

    >>> with open("PDB/1A8O.pdb") as handle:
    ...     for record in PdbAtomIterator(handle):
    ...         print("Record id %s, chain %s" % (record.id, record.annotations["chain"]))
    ...
    Record id 1A8O:A, chain A

    """
    # TODO - Add record.annotations to the doctest, esp the residues (not working?)

    # Only import PDB when needed, to avoid/delay NumPy dependency in SeqIO
    from Bio.PDB import PDBParser

    # Deduce the PDB ID from the PDB header
    # ENH: or filename?
    from Bio.File import UndoHandle
    undo_handle = UndoHandle(handle)
    firstline = undo_handle.peekline()

    # check if file is empty
    if firstline == '':
        raise ValueError("Empty file.")

    if firstline.startswith("HEADER"):
        pdb_id = firstline[62:66]
    else:
        warnings.warn("First line is not a 'HEADER'; can't determine PDB ID. "
                      "Line: %r" % firstline, BiopythonParserWarning)
        pdb_id = '????'

    struct = PDBParser().get_structure(pdb_id, undo_handle)
    for record in AtomIterator(pdb_id, struct):
        # The PDB header was loaded as a dictionary, so let's reuse it all
        record.annotations.update(struct.header)

        # ENH - add letter annotations -- per-residue info, e.g. numbers

        yield record


# We can't be sure that we have the enum module in python 2.7, so
# we will refer to the field names directly in the code.
PDBX_POLY_SEQ_SCHEME_FIELDS = (
    "_pdbx_poly_seq_scheme.asym_id",     # Chain ID
    "_pdbx_poly_seq_scheme.mon_id")      # Residue type

STRUCT_REF_FIELDS = (
    "_struct_ref.id",                    # ID of this reference
    "_struct_ref.db_name",               # Name of the database
    "_struct_ref.db_code",               # Code for this entity
    "_struct_ref.pdbx_db_accession")     # DB accession ID of ref

STRUCT_REF_SEQ_FIELDS = (
    "_struct_ref_seq.ref_id",            # Pointer to _struct_ref
    "_struct_ref_seq.pdbx_PDB_id_code",  # PDB ID of this structure
    "_struct_ref_seq.pdbx_strand_id")    # Chain ID of the reference


def CifSeqresIterator(handle):
    """Return SeqRecord objects for each chain in an mmCIF file.

    The sequences are derived from the _entity_poly_seq entries in the mmCIF
    file, not the atoms of the 3D structure.

    Specifically, these mmCIF records are handled: _pdbx_poly_seq_scheme and
    _struct_ref_seq. The _pdbx_poly_seq records contain sequence information,
    and the _struct_ref_seq records contain database cross-references.

    See:
    http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v40.dic/Categories/pdbx_poly_seq_scheme.html
    and
    http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Categories/struct_ref_seq.html

    This gets called internally via Bio.SeqIO for the sequence-based
    interpretation of the mmCIF file format:

    >>> from Bio import SeqIO
    >>> for record in SeqIO.parse("PDB/1A8O.cif", "cif-seqres"):
    ...     print("Record id %s, chain %s" % (record.id, record.annotations["chain"]))
    ...     print(record.dbxrefs)
    ...
    Record id 1A8O:A, chain A
    ['UNP:P12497', 'UNP:POL_HV1N5']

    Equivalently,

    >>> with open("PDB/1A8O.cif") as handle:
    ...     for record in CifSeqresIterator(handle):
    ...         print("Record id %s, chain %s" % (record.id, record.annotations["chain"]))
    ...         print(record.dbxrefs)
    ...
    Record id 1A8O:A, chain A
    ['UNP:P12497', 'UNP:POL_HV1N5']

    Note the chain is recorded in the annotations dictionary, and any mmCIF
    _struct_ref_seq entries are recorded in the database cross-references list.
    """
    # Late-binding import to avoid circular dependency on SeqIO in Bio.SeqUtils
    from Bio.SeqUtils import seq1

    # Only import PDB when needed, to avoid/delay NumPy dependency in SeqIO
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict

    chains = collections.defaultdict(list)
    metadata = collections.defaultdict(list)
    records = MMCIF2Dict(handle)

    # Explicitly convert records to list (See #1533).
    # If an item is not present, use an empty list
    for field in (
            PDBX_POLY_SEQ_SCHEME_FIELDS
            + STRUCT_REF_SEQ_FIELDS
            + STRUCT_REF_FIELDS):
        if field not in records:
            records[field] = []
        elif not isinstance(records[field], list):
            records[field] = [records[field]]

    for asym_id, mon_id in zip(records["_pdbx_poly_seq_scheme.asym_id"],
                               records["_pdbx_poly_seq_scheme.mon_id"]):
        mon_id_1l = seq1(mon_id, custom_map=protein_letters_3to1)
        chains[asym_id].append(mon_id_1l)

    # Build a dict of _struct_ref records, indexed by the id field:
    struct_refs = {}
    for fields in zip(records["_struct_ref.id"],
                      records["_struct_ref.db_name"],
                      records["_struct_ref.db_code"],
                      records["_struct_ref.pdbx_db_accession"]):
        ref_id, db_name, db_code, db_acc = fields
        struct_refs[ref_id] = {
            "database": db_name,
            "db_id_code": db_code,
            "db_acc": db_acc}

    # Look through _struct_ref_seq records, look up the corresponding
    # _struct_ref and add an entry to the metadata list for this chain.
    for fields in zip(records["_struct_ref_seq.ref_id"],
                      records["_struct_ref_seq.pdbx_PDB_id_code"],
                      records["_struct_ref_seq.pdbx_strand_id"]):
        ref_id, pdb_id, chain_id = fields
        struct_ref = struct_refs[ref_id]

        # The names here mirror those in PdbIO
        metadata[chain_id].append({'pdb_id': pdb_id})
        metadata[chain_id][-1].update(struct_ref)

    for chn_id, residues in sorted(chains.items()):
        record = SeqRecord(Seq(''.join(residues), generic_protein))
        record.annotations = {"chain": chn_id}
        if chn_id in metadata:
            m = metadata[chn_id][0]
            record.id = record.name = "%s:%s" % (m['pdb_id'], chn_id)
            record.description = ("%s:%s %s" % (m['database'],
                                                m['db_acc'],
                                                m['db_id_code']))
            for melem in metadata[chn_id]:
                record.dbxrefs.extend([
                    "%s:%s" % (melem['database'], melem['db_acc']),
                    "%s:%s" % (melem['database'], melem['db_id_code'])])
        else:
            record.id = chn_id
        yield record


def CifAtomIterator(handle):
    """Return SeqRecord objects for each chain in a PDB file.

    The sequences are derived from the 3D structure (ATOM records), not the
    SEQRES lines in the PDB file header.

    Unrecognised three letter amino acid codes (e.g. "CSD") from HETATM entries
    are converted to "X" in the sequence.

    In addition to information from the PDB header (which is the same for all
    records), the following chain specific information is placed in the
    annotation:

    record.annotations["residues"] = List of residue ID strings
    record.annotations["chain"] = Chain ID (typically A, B ,...)
    record.annotations["model"] = Model ID (typically zero)

    Where amino acids are missing from the structure, as indicated by residue
    numbering, the sequence is filled in with 'X' characters to match the size
    of the missing region, and  None is included as the corresponding entry in
    the list record.annotations["residues"].

    This function uses the Bio.PDB module to do most of the hard work. The
    annotation information could be improved but this extra parsing should be
    done in parse_pdb_header, not this module.

    This gets called internally via Bio.SeqIO for the atom based interpretation
    of the PDB file format:

    >>> from Bio import SeqIO
    >>> for record in SeqIO.parse("PDB/1A8O.cif", "cif-atom"):
    ...     print("Record id %s, chain %s" % (record.id, record.annotations["chain"]))
    ...
    Record id 1A8O:A, chain A

    Equivalently,

    >>> with open("PDB/1A8O.cif") as handle:
    ...     for record in CifAtomIterator(handle):
    ...         print("Record id %s, chain %s" % (record.id, record.annotations["chain"]))
    ...
    Record id 1A8O:A, chain A

    """
    # TODO - Add record.annotations to the doctest, esp the residues (not working?)

    # Only import parser when needed, to avoid/delay NumPy dependency in SeqIO
    from Bio.PDB.MMCIFParser import MMCIFParser
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict

    # The PdbAtomIterator uses UndoHandle to peek at the first line and get the
    # PDB ID. The equivalent for mmCIF is the _entry.id field. AFAIK, the mmCIF
    # format does not constrain the order of fields, so we need to parse the
    # entire file using MMCIF2Dict. We copy the contents of the handle into a
    # StringIO buffer first, so that both MMCIF2Dict and MMCIFParser can
    # consume the handle.
    buffer = StringIO()
    shutil.copyfileobj(handle, buffer)

    # check if file is empty
    if len(buffer.getvalue()) == 0:
        raise ValueError("Empty file.")

    buffer.seek(0)
    mmcif_dict = MMCIF2Dict(buffer)
    if "_entry.id" in mmcif_dict:
        pdb_id = mmcif_dict["_entry.id"]
        if isinstance(pdb_id, list):
            pdb_id = pdb_id[0]
    else:
        warnings.warn("Could not find the '_entry.id' field; can't determine "
                      "PDB ID.", BiopythonParserWarning)
        pdb_id = '????'

    buffer.seek(0)
    struct = MMCIFParser().get_structure(pdb_id, buffer)
    for record in AtomIterator(pdb_id, struct):
        yield record


if __name__ == '__main__':
    from Bio._utils import run_doctest
    run_doctest(verbose=0)
