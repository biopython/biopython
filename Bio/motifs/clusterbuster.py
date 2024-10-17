# Copyright 2015 by Gert Hulselmans.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Parse Cluster Buster position frequency matrix files."""

from Bio import motifs


class Record(list):
    """Class to store the information in a Cluster Buster matrix table.

    The record inherits from a list containing the individual motifs.
    """

    def __str__(self):
        """Return a string representation of the motifs in the Record object."""
        return "\n".join(str(motif) for motif in self)


def read(handle):
    """Read motifs in Cluster Buster position frequency matrix format from a file handle.

    Cluster Buster motif format: https://bu.wenglab.org/cluster-buster/help/cis-format.html
    """
    motif_nbr = 0
    record = Record()
    nucleotide_counts = {"A": [], "C": [], "G": [], "T": []}
    motif_name = ""
    motif_gap = None
    motif_weight = None

    for line in handle:
        line = line.strip()
        if line:
            if line.startswith(">"):
                if motif_nbr != 0:
                    motif = motifs.Motif(alphabet="GATC", counts=nucleotide_counts)
                    motif.name = motif_name
                    motif.gap = motif_gap
                    motif.weight = motif_weight
                    record.append(motif)

                motif_name = line[1:].strip()
                nucleotide_counts = {"A": [], "C": [], "G": [], "T": []}
                motif_gap = None
                motif_weight = None
                motif_nbr += 1
            else:
                if line.startswith("# GAP"):
                    motif_gap = float(line.split()[2])
                    continue
                elif line.startswith("# WEIGHT"):
                    motif_weight = float(line.split()[2])
                    continue
                elif line.startswith("#"):
                    continue

                matrix_columns = line.split()

                if len(matrix_columns) == 4:
                    [
                        nucleotide_counts[nucleotide].append(float(nucleotide_count))
                        for nucleotide, nucleotide_count in zip(
                            ["A", "C", "G", "T"], matrix_columns
                        )
                    ]

    motif = motifs.Motif(alphabet="GATC", counts=nucleotide_counts)
    motif.name = motif_name
    motif.gap = motif_gap
    motif.weight = motif_weight
    record.append(motif)

    return record


def write(motifs, precision=0):
    """Return the representation of motifs in Cluster Buster position frequency matrix format.

    By default (`precision=0`) Cluster Buster position frequency matrices will be written
    with integer values.
    If a higher precision value is set, Cluster Buster position frequency matrices will be
    written as floats with `x` decimal places.
    """
    lines = []
    for m in motifs:
        lines.append(f">{m.name}\n")
        if m.weight:
            lines.append(f"# WEIGHT: {m.weight}\n")
        if m.gap:
            lines.append(f"# GAP: {m.gap}\n")
        for ACGT_counts in zip(
            m.counts["A"], m.counts["C"], m.counts["G"], m.counts["T"]
        ):
            lines.append(
                "{:0.0f}\t{:0.0f}\t{:0.0f}\t{:0.0f}\n".replace(
                    "0f", str(int(precision)) + "f"
                ).format(*ACGT_counts)
            )

    # Finished; glue the lines together.
    text = "".join(lines)

    return text
