# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Map residues of two structures to each other based on a FASTA alignment."""

from typing import Optional

from Bio.Align import Alignment, MultipleSeqAlignment
from Bio.Data import PDBData
from Bio.PDB import Selection
from Bio.PDB.Model import Model
from Bio.PDB.Polypeptide import is_aa


class StructureAlignment:
    """Class to align two structures based on an alignment of their sequences."""

    def __init__(
        self,
        fasta_align: Optional[MultipleSeqAlignment | Alignment] = None,
        m1: Model | None = None,
        m2: Model | None = None,
        si: int = 0,
        sj: int = 1,
    ) -> None:
        """Initialize.

        Attributes:
         - fasta_align - Alignment object / MSA object or None (if None, one will be generated automatically)
         - m1, m2 - two models (Bio.PDB.Model.Model objects). Their default values are set to None to maintain legacy code, but they CANNOT be None
         - si, sj - the sequences in the Alignment object that correspond to the structures

        """
        if fasta_align is None:
            raise NotImplementedError(
                "fasta_align cannot be None for now to pass precommit"
            )
            # add code for auto generation
        elif isinstance(fasta_align, MultipleSeqAlignment):
            ncolumns = fasta_align.get_alignment_length()
        elif isinstance(fasta_align, Alignment):
            nrows, ncolumns = fasta_align.shape
        else:
            raise ValueError(
                "Invalid alignment object, must either be an Alignment object, "
                "MultipleSeqAlignment object, or None"
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
            column = fasta_align[:, i]
            aa1 = column[si]
            aa2 = column[sj]
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
