"""Base function for sequence parsers that read structures Bio.PDB parsers.

Once a parser from Bio.PDB has been used to load a structure into a
Bio.PDB.Structure.Structure object, there is no difference in how the sequence
parser interprets the residue sequence. The functions in this module may be
used by SeqIO modules wishing to parse sequences from lists of residues.

Calling funtions must pass a Bio.PDB.Structure.Structure object.

Note: This module lives in the Bio.SeqIO package, but is not registered as a
format in the Bio.SeqIO._FormatToIterator dictionary.
"""
import warnings

from Bio import BiopythonWarning
from Bio.Alphabet import generic_protein
from Bio.Data.SCOPData import protein_letters_3to1
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def AtomIterator(pdb_id, struct):
    """Return SeqRecords from Structure objects.

    See Bio.SeqIO.PdbIO.PdbAtomIterator and Bio.SeqIO.CifIO.CifAtomIterator for
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
                                  BiopythonWarning)
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
