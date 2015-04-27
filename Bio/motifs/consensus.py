# Copyright 2015 by Gert Hulselmans.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Converting list of degenerate consensus sequences to list of motifs."""

from Bio import motifs
from Bio.Alphabet import IUPAC
from Bio.Data import IUPACData
from Bio.Seq import Seq


__docformat__ = "restructuredtext en"


def __get_unambigous_nulceotides_from_IUPAC_nucleotide_codes():
    """Create a dictionary with all ambigous DNA nucleotides as keys and with the corresponding list of unambigous DNA
    nucleotides as values. The list of unambigous nucleotides is repeated till the total number of unambigous DNA
    nucleotides is equal to 12.
    """

    unambigous_nucleotide_sequences_for_IUPAC_nucleotide_codes = dict()

    for ambigous_DNA_nucleotide, unambigous_DNA_nucleotides in IUPACData.ambiguous_dna_values.iteritems():
        # Get the unambigous DNA nucleotides for a ambigous DNA nucleotide and repeat the sequence till the total
        # number of unambigous DNA nucleotides is equal to 12.
        unambigous_nucleotide_sequences_for_IUPAC_nucleotide_codes[ambigous_DNA_nucleotide] \
            = sorted(unambigous_DNA_nucleotides) * (12 / len(unambigous_DNA_nucleotides))

    return unambigous_nucleotide_sequences_for_IUPAC_nucleotide_codes


unambigous_nucleotide_sequences_for_IUPAC_nucleotide_codes = __get_unambigous_nulceotides_from_IUPAC_nucleotide_codes()


class Record(list):
    """A Bio.motifs.consensus.Record represents a list of motifs made from a list of degenerate consensus sequences."""

    def __str__(self):
        return "\n".join(str(motif) for motif in self)


def read(handle):
    """Read degenerate consensus sequence(s) from a file and convert them to motifs.

    The consensus sequence file should contain two columns (separated by TABs):
      - column 1: motif name
      - column 2: degenerate consensus sequence
    """

    record = Record()
    for line in handle:
        line = line.strip()
        if line:
            motif_name, degenerate_consensus_sequence = line.split("\t")
            motif = create_motif_from_degenerate_consensus_sequence(motif_name, degenerate_consensus_sequence)
            record.append(motif)
    return record


def create_motif_from_degenerate_consensus_sequence(motif_name, degenerate_consensus_sequence):
    """Create motif from a degenerate consensus sequence.

    For example:

    >>> import Bio.motifs.consensus
    >>> m = Bio.motifs.consensus.create_motif_from_degenerate_consensus_sequence("consensus-motifID1", "ACGTRYSWKMBDHVN")
    >>> print(m.degenerate_consensus)
    ACGTRYSWKMBDHVN
    >>> print(m.counts)
            0      1      2      3      4      5      6      7      8      9     10     11     12     13     14
    A:  12.00   0.00   0.00   0.00   6.00   0.00   0.00   6.00   0.00   6.00   0.00   4.00   4.00   4.00   3.00
    C:   0.00  12.00   0.00   0.00   0.00   6.00   6.00   0.00   0.00   6.00   4.00   0.00   4.00   4.00   3.00
    G:   0.00   0.00  12.00   0.00   6.00   0.00   6.00   0.00   6.00   0.00   4.00   4.00   0.00   4.00   3.00
    T:   0.00   0.00   0.00  12.00   0.00   6.00   0.00   6.00   6.00   0.00   4.00   4.00   4.00   0.00   3.00
    <BLANKLINE>
    """

    empty_sequence = [" "] * len(degenerate_consensus_sequence)
    unambigous_DNA_list = [list(empty_sequence) for i in range(0, 12)]

    for nucleotide_position, ambigous_DNA_nucleotide in enumerate(degenerate_consensus_sequence):
        sequence_number = 0
        for unambigous_DNA_nucleotide in unambigous_nucleotide_sequences_for_IUPAC_nucleotide_codes[ambigous_DNA_nucleotide]:
            unambigous_DNA_list[sequence_number][nucleotide_position] = unambigous_DNA_nucleotide
            sequence_number += 1

    unambigous_DNA_sequences = [Seq("".join(unambigous_DNA), IUPAC.unambiguous_dna)
                                for unambigous_DNA in unambigous_DNA_list]
    instances = motifs.Instances(unambigous_DNA_sequences)

    motif = motifs.Motif(alphabet=IUPAC.unambiguous_dna, instances=instances)
    motif.name = motif_name

    return motif
