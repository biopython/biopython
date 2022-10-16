# Copyright 2003 Yair Benita.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Methods for codon usage calculations."""

import math
import warnings
from .CodonUsageIndices import SharpEcoliIndex
from Bio import SeqIO  # To parse a FASTA file
from Bio import BiopythonDeprecationWarning


warnings.warn(
    "This module has been DEPRECATED. Please use the CodonAdaptationIndex "
    "class in Bio.SeqUtils instead. Note that this class has been updated to "
    "use modern Python, and may give slightly different results from the "
    "CodonAdaptationIndex class in Bio.SeqUtils.CodonUsage, as the code was "
    "changed to be consistent with the published paper by Sharp and Li. The "
    "code in the old CodonAdaptationIndex class in Bio.SeqUtils.CodonUsage was "
    "not changed.",
    BiopythonDeprecationWarning,
)


# Turn black code style off
# fmt: off

CodonsDict = {
    "TTT": 0, "TTC": 0, "TTA": 0, "TTG": 0,
    "CTT": 0, "CTC": 0, "CTA": 0, "CTG": 0,
    "ATT": 0, "ATC": 0, "ATA": 0, "ATG": 0,
    "GTT": 0, "GTC": 0, "GTA": 0, "GTG": 0,
    "TAT": 0, "TAC": 0, "TAA": 0, "TAG": 0,
    "CAT": 0, "CAC": 0, "CAA": 0, "CAG": 0,
    "AAT": 0, "AAC": 0, "AAA": 0, "AAG": 0,
    "GAT": 0, "GAC": 0, "GAA": 0, "GAG": 0,
    "TCT": 0, "TCC": 0, "TCA": 0, "TCG": 0,
    "CCT": 0, "CCC": 0, "CCA": 0, "CCG": 0,
    "ACT": 0, "ACC": 0, "ACA": 0, "ACG": 0,
    "GCT": 0, "GCC": 0, "GCA": 0, "GCG": 0,
    "TGT": 0, "TGC": 0, "TGA": 0, "TGG": 0,
    "CGT": 0, "CGC": 0, "CGA": 0, "CGG": 0,
    "AGT": 0, "AGC": 0, "AGA": 0, "AGG": 0,
    "GGT": 0, "GGC": 0, "GGA": 0, "GGG": 0}

# Turn black code style on
# fmt: on


# this dictionary shows which codons encode the same AA
SynonymousCodons = {
    "CYS": ["TGT", "TGC"],
    "ASP": ["GAT", "GAC"],
    "SER": ["TCT", "TCG", "TCA", "TCC", "AGC", "AGT"],
    "GLN": ["CAA", "CAG"],
    "MET": ["ATG"],
    "ASN": ["AAC", "AAT"],
    "PRO": ["CCT", "CCG", "CCA", "CCC"],
    "LYS": ["AAG", "AAA"],
    "STOP": ["TAG", "TGA", "TAA"],
    "THR": ["ACC", "ACA", "ACG", "ACT"],
    "PHE": ["TTT", "TTC"],
    "ALA": ["GCA", "GCC", "GCG", "GCT"],
    "GLY": ["GGT", "GGG", "GGA", "GGC"],
    "ILE": ["ATC", "ATA", "ATT"],
    "LEU": ["TTA", "TTG", "CTC", "CTT", "CTG", "CTA"],
    "HIS": ["CAT", "CAC"],
    "ARG": ["CGA", "CGC", "CGG", "CGT", "AGG", "AGA"],
    "TRP": ["TGG"],
    "VAL": ["GTA", "GTC", "GTG", "GTT"],
    "GLU": ["GAG", "GAA"],
    "TYR": ["TAT", "TAC"],
}


class CodonAdaptationIndex:
    """A codon adaptation index (CAI) implementation.

    Implements the codon adaptation index (CAI) described by Sharp and
    Li (Nucleic Acids Res. 1987 Feb 11;15(3):1281-95).

    NOTE - This implementation does not currently cope with alternative genetic
    codes: only the synonymous codons in the standard table are considered.
    """

    def __init__(self):
        """Initialize the class."""
        self.index = {}
        self.codon_count = {}

    # use this method with predefined CAI index
    def set_cai_index(self, index):
        """Set up an index to be used when calculating CAI for a gene.

        Just pass a dictionary similar to the SharpEcoliIndex in the
        CodonUsageIndices module.
        """
        self.index = index

    def generate_index(self, fasta_file):
        """Generate a codon usage index from a FASTA file of CDS sequences.

        Takes a location of a Fasta file containing CDS sequences
        (which must all have a whole number of codons) and generates a codon
        usage index.

        """
        # first make sure we're not overwriting an existing index:
        if self.index != {} or self.codon_count != {}:
            raise ValueError(
                "an index has already been set or a codon count "
                "has been done. Cannot overwrite either."
            )

        # count codon occurrences in the file.
        self._count_codons(fasta_file)

        # now to calculate the index we first need to sum the number of times
        # synonymous codons were used all together.
        for aa in SynonymousCodons:
            codons = SynonymousCodons[aa]

            count_max = max(self.codon_count[codon] for codon in codons)
            if count_max == 0:  # the residue does not occur at all
                for codon in codons:
                    self.index[codon] = None
            else:
                # now generate the index W=RCSUi/RCSUmax = COUNTi/COUNTmax:
                # see equation 2 in Sharp & Li 1987 NAR
                for codon in codons:
                    self.index[codon] = self.codon_count[codon] / count_max

    def cai_for_gene(self, dna_sequence):
        """Calculate the CAI (float) for the provided DNA sequence (string).

        This method uses the Index (either the one you set or the one you
        generated) and returns the CAI for the DNA sequence.
        """
        cai_value, cai_length = 0, 0

        # if no index is set or generated, the default SharpEcoliIndex will
        # be used.
        if self.index == {}:
            self.set_cai_index(SharpEcoliIndex)

        dna_sequence = dna_sequence.upper()

        for i in range(0, len(dna_sequence), 3):
            codon = dna_sequence[i : i + 3]
            if codon in self.index:
                # these two codons are always one, exclude them:
                if codon not in ["ATG", "TGG"]:
                    cai_value += math.log(self.index[codon])
                    cai_length += 1
            # some indices may not include stop codons:
            elif codon not in ["TGA", "TAA", "TAG"]:
                raise TypeError(f"illegal codon in sequence: {codon}.\n{self.index}")

        return math.exp(cai_value / (cai_length - 1.0))

    def _count_codons(self, fasta_file):
        with open(fasta_file) as handle:

            # make the codon dictionary local
            self.codon_count = CodonsDict.copy()

            # iterate over sequence and count all the codons in the FastaFile.
            for record in SeqIO.parse(handle, "fasta"):
                sequence = record.seq.upper()
                for i in range(0, len(sequence), 3):
                    codon = sequence[i : i + 3]
                    try:
                        self.codon_count[codon] += 1
                    except KeyError:
                        raise ValueError(
                            f"illegal codon '{codon}' in gene: {record.id}"
                        ) from None

    def __str__(self):
        lines = []
        for i in sorted(self.index):
            line = f"{i}\t{self.index[i]:.3f}"
            lines.append(line)
        return "\n".join(lines) + "\n"

    def print_index(self):
        """Print out the index used.

        This just gives the index when the objects is printed.

        >>> from Bio.SeqUtils.CodonUsage import CodonAdaptationIndex
        >>> from Bio.SeqUtils.CodonUsageIndices import SharpEcoliIndex
        >>> cai = CodonAdaptationIndex()
        >>> cai.print_index()
        >>> cai.set_cai_index({'GCA':None, 'GCC':1.0})
        >>> cai.print_index()  #doctest: +NORMALIZE_WHITESPACE
        GCA	None
        GCC	1.000
        >>> cai.set_cai_index(SharpEcoliIndex)
        >>> cai.print_index()  #doctest: +NORMALIZE_WHITESPACE
        AAA	1.000
        AAC	1.000
        AAG	0.253
        AAT	0.051
        ACA	0.076
        ACC	1.000
        ACG	0.099
        ACT	0.965
        AGA	0.004
        AGC	0.410
        AGG	0.002
        AGT	0.085
        ATA	0.003
        ATC	1.000
        ATG	1.000
        ATT	0.185
        CAA	0.124
        CAC	1.000
        CAG	1.000
        CAT	0.291
        CCA	0.135
        CCC	0.012
        CCG	1.000
        CCT	0.070
        CGA	0.004
        CGC	0.356
        CGG	0.004
        CGT	1.000
        CTA	0.007
        CTC	0.037
        CTG	1.000
        CTT	0.042
        GAA	1.000
        GAC	1.000
        GAG	0.259
        GAT	0.434
        GCA	0.586
        GCC	0.122
        GCG	0.424
        GCT	1.000
        GGA	0.010
        GGC	0.724
        GGG	0.019
        GGT	1.000
        GTA	0.495
        GTC	0.066
        GTG	0.221
        GTT	1.000
        TAC	1.000
        TAT	0.239
        TCA	0.077
        TCC	0.744
        TCG	0.017
        TCT	1.000
        TGC	1.000
        TGG	1.000
        TGT	0.500
        TTA	0.020
        TTC	1.000
        TTG	0.020
        TTT	0.296
        """

        warnings.warn(
            "The print_index method is deprecated; instead of "
            "self.print_index(), please use print(self).",
            BiopythonDeprecationWarning,
        )
        print(self)