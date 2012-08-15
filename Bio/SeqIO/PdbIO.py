# Copyright 2012 by Eric Talevich.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import collections

from Bio.Alphabet import generic_protein
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SCOP.three_to_one_dict import to_one_letter_code


def PdbSeqresIterator(handle):
    """Returns SeqRecord objects for each chain in a PDB file.

    The sequences are derived from the SEQRES lines in the
    PDB file header, not the atoms of the 3D structure.

    Specifically, these PDB records are handled: DBREF, SEQADV, SEQRES, MODRES

    See: http://www.wwpdb.org/documentation/format23/sect3.html
    """
    chains = collections.defaultdict(list)
    metadata = collections.defaultdict(list)
    for line in handle:
        rec_name = line[0:6].strip()
        if rec_name == 'SEQRES':
            # NB: We only actually need chain ID and the residues here
            # Serial number of the SEQRES record for the current chain.
            # Starts at 1 and increments by one each line.
            # Reset to 1 for each chain.
            # ser_num = int(line[8:10])
            # Chain identifier. This may be any single legal character,
            # including a blank which is used if there is only one chain.
            chn_id = line[11]
            # Number of residues in the chain (repeated on every record)
            # num_res = int(line[13:17])
            residues = [to_one_letter_code.get(res, 'X')
                        for res in line[19:].split()]
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

    for chn_id, residues in sorted(chains.iteritems()):
        record = SeqRecord(Seq(''.join(residues), generic_protein))
        record.annotations = {"chain": chn_id}
        if chn_id in metadata:
            m = metadata[chn_id][0]
            record.id = record.name = "%s:%s" % (m['pdb_id'], chn_id)
            record.description = ("%s:%s %s" % (m['database'],
                                                m['db_acc'],
                                                m['db_id_code']))
            # XXX Is this the best way to do this?
            for melem in metadata[chn_id]:
                record.dbxrefs.extend([
                    "%s:%s" % (melem['database'], melem['db_acc']),
                    "%s:%s" % (melem['database'], melem['db_id_code'])])
        else:
            record.id = chn_id
        yield record


if __name__ == '__main__':
    # Test
    import sys
    from Bio import SeqIO
    with open(sys.argv[1]) as handle:
        result = PdbSeqresIterator(handle)
        results = list(result)
        for r in results:
            print r
        SeqIO.write(results, sys.stdout, 'fasta')

