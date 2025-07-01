# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Map residues of two structures to each other based on a FASTA alignment."""

from typing import Optional
import warnings

from Bio.Align import Alignment, MultipleSeqAlignment, PairwiseAligner
from Bio.Data import PDBData
from Bio.PDB import Selection
from Bio.PDB.Model import Model
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.Residue import Residue
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class StructureAlignment:
    """Class to align two structures based on an alignment of their sequences."""

    def __init__(
        self,
        fasta_align: Optional[MultipleSeqAlignment | Alignment] = None,
        m1: Model | None = None,
        m2: Model | None = None,
        si: int = 0,
        sj: int = 1,
        aligner: Optional[PairwiseAligner] = None,
    ) -> None:
        """Initialize.

        Attributes:
         - fasta_align - Alignment object / MSA object or None (if None, one will be generated automatically)
         - m1, m2 - two models (Bio.PDB.Model.Model objects). Their default values are set to None to maintain legacy code, but they CANNOT be None
         - si, sj - the sequences in the Alignment object that correspond to the structures
         - aligner - Optional aligner, mutually exclusive with fasta_align, allows for customization of automatic
           alignment, otherwise blastp defaults will be used

        """
        # Check that fasta_align and aligner are not both provided
        if fasta_align is not None and aligner is not None:
            raise ValueError(
                "fasta_align and aligner cannot be both provided, fasta_align uses precomputed alignments",
                "while the aligner is to allow a custom alignment tool to be used for automatically aligning",
                "input sequences",
            )

        self.aligner = aligner

        # Validate that models are provided
        if m1 is None or m2 is None:
            raise ValueError("Both m1 and m2 models must be provided")

        if fasta_align is None:
            fasta_align = self._generate_alignment_from_models(
                m1, m2
            )  # MultipleSeqAlignment
        else:
            # if fasta_align is explicitly provided, raise DeprecationWarning letting user
            # know that fasta_align is no longer required
            warnings.warn(
                "fasta_align is no longer a required argument (will be automatically computed "
                "if not provided), and the function signature will change in a future release",
                DeprecationWarning,
                stacklevel=2,
            )

        if isinstance(fasta_align, MultipleSeqAlignment):
            ncolumns = fasta_align.get_alignment_length()
        elif isinstance(fasta_align, Alignment):
            _, ncolumns = fasta_align.shape
        else:
            raise ValueError(
                f"Could not infer length from alignment object: {fasta_align}. "
                "Alignment must be of type MultipleSeqAlignment or Alignment."
            )

        # Get the residues in the models
        rl1 = Selection.unfold_entities(m1, "R")
        rl2 = Selection.unfold_entities(m2, "R")
        # Residue positions
        p1 = 0
        p2 = 0
        # Map equivalent residues to each other
        map12 = {}
        map21 = {}
        # List of residue pairs (None if -)
        duos = []
        for i in range(ncolumns):
            aa1 = fasta_align[si, i]
            aa2 = fasta_align[sj, i]

            if aa1 != "-":
                # Position in seq1 is not -
                while True:
                    # Loop until an aa is found
                    r1 = rl1[p1]
                    p1 = p1 + 1
                    if is_aa(r1):
                        break
                self._test_equivalence(r1, aa1)
            else:
                r1 = None
            if aa2 != "-":
                # Position in seq2 is not -
                while True:
                    # Loop until an aa is found
                    r2 = rl2[p2]
                    p2 = p2 + 1
                    if is_aa(r2):
                        break
                self._test_equivalence(r2, aa2)
            else:
                r2 = None
            if r1:
                # Map residue in seq1 to its equivalent in seq2
                map12[r1] = r2
            if r2:
                # Map residue in seq2 to its equivalent in seq1
                map21[r2] = r1
            # Append aligned pair (r is None if gap)
            duos.append((r1, r2))
        self.map12 = map12
        self.map21 = map21
        self.duos = duos

    def _generate_alignment_from_models(self, m1, m2) -> Alignment:
        """Generate a MultipleSeqAlignment from two protein models .

        Uses Bio.Align.PairwiseAligner to create a proper sequence alignment
        with gaps, rather than simply placing sequences side-by-side.
        This ensures that homologous residues are properly aligned.
        """
        seq1_record = self._extract_sequence_from_model(m1, "structure1")
        seq2_record = self._extract_sequence_from_model(m2, "structure2")

        if self.aligner is not None:
            assert isinstance(
                self.aligner, PairwiseAligner
            ), f"custom aligner must be a PairwiseAligner object, not {type(self.aligner)}"
        else:
            aligner = PairwiseAligner("blastp")

        alignments = aligner.align(seq1_record.seq, seq2_record.seq)
        best_alignment = alignments[0]

        return best_alignment

    def _extract_sequence_from_model(self, model, seq_id):
        """Extract amino acid sequence from a protein model."""
        residues = Selection.unfold_entities(model, "R")
        sequence = ""

        for residue in residues:
            if isinstance(residue, Residue) and is_aa(residue):
                resname = residue.get_resname()
                aa_code = PDBData.protein_letters_3to1_extended.get(resname, "X")
                sequence += aa_code

        return SeqRecord(Seq(sequence), id=seq_id)

    def _test_equivalence(self, r1, aa1):
        """Test if aa in sequence fits aa in structure (PRIVATE)."""
        resname = r1.get_resname()
        resname = PDBData.protein_letters_3to1_extended[resname]
        assert aa1 == resname

    def get_maps(self):
        """Map residues between the structures.

        Return two dictionaries that map a residue in one structure to
        the equivealent residue in the other structure.
        """
        return self.map12, self.map21

    def get_iterator(self):
        """Create an iterator over all residue pairs."""
        for i in range(len(self.duos)):
            yield self.duos[i]
