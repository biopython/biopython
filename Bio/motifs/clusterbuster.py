# Copyright 2015 by Gert Hulselmans.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Parse Cluster Buster position frequency matrix files."""

from Bio import motifs
from Bio.Alphabet import IUPAC


class Record(list):
    """Class to store the information in a Cluster Buster matrix table.

    The record inherits from a list containing the individual motifs.
    """

    def __str__(self):
        return "\n".join(str(motif) for motif in self)


def read(handle):
    """Read motifs in Cluster Buster position frequency matrix format from a file handle.

    Cluster Buster motif format: http://zlab.bu.edu/cluster-buster/help/cis-format.html
    """
    motif_nbr = 0
    record = Record()
    nucleotide_counts = {'A': [], 'C': [], 'G': [], 'T': []}
    motif_name = ""

    for line in handle:
        line = line.strip()
        if line:
            if line.startswith('>'):

                if motif_nbr != 0:
                    motif = motifs.Motif(alphabet=IUPAC.unambiguous_dna, counts=nucleotide_counts)
                    motif.name = motif_name
                    record.append(motif)

                motif_name = line[1:].strip()
                nucleotide_counts = {'A': [], 'C': [], 'G': [], 'T': []}
                motif_nbr += 1
            else:
                if line.startswith('#'):
                    continue

                matrix_columns = line.split()

                if len(matrix_columns) == 4:
                    [nucleotide_counts[nucleotide].append(float(nucleotide_count))
                     for nucleotide, nucleotide_count in zip(['A', 'C', 'G', 'T'], matrix_columns)]

    motif = motifs.Motif(alphabet=IUPAC.unambiguous_dna, counts=nucleotide_counts)
    motif.name = motif_name
    record.append(motif)

    return record


def write(motifs):
    """Return the representation of motifs in Cluster Buster position frequency matrix format."""
    lines = []
    for m in motifs:
        line = '>{0}\n'.format(m.name)
        lines.append(line)
        for ACGT_counts in zip(m.counts['A'], m.counts['C'], m.counts['G'], m.counts['T']):
            lines.append('{0:0.0f}\t{1:0.0f}\t{2:0.0f}\t{3:0.0f}\n'.format(*ACGT_counts))

    # Finished; glue the lines together.
    text = ''.join(lines)

    return text
