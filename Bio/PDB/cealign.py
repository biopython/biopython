# Copyright (C) 2022, Joao Rodrigues (j.p.g.l.m.rodrigues@gmail.com)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Protein Structural Alignment using Combinatorial Extension.

Python code written by Joao Rodrigues. C++ code and Python/C++ interface
adapted from open-source Pymol and originally written by Jason Vertrees.
The original license and notices are available in `cealign` folder.

Reference
---------
Shindyalov, I.N., Bourne P.E. (1998).
"Protein structure alignment by incremental combinatorial extension (CE)
of the optimal path". Protein Engineering. 11 (9): 739–747. PMID 9796821.
"""

import numpy as np

from Bio.PDB.ccealign import run_cealign
from Bio.PDB.PDBExceptions import PDBException
from Bio.PDB.qcprot import QCPSuperimposer

_RESID_SORTER = lambda r: r.id[1]  # noqa: E731


class CEAligner:
    """Protein Structure Alignment by Combinatorial Extension."""

    def __init__(self, window_size=8, max_gap=30):
        """Superimpose one set of atoms onto another using structural data.

        Structures are superimposed using guide atoms, CA and C4', for protein
        and nucleic acid molecules respectively.

        Parameters
        ----------
        window_size : float, optional
            CE algorithm parameter. Used to define paths when building the
            CE similarity matrix. Default is 8.
        max_gap : float, optional
            CE algorithm parameter. Maximum gap size. Default is 30.
        """
        assert window_size > 0, "window_size must be greater than 0"
        assert max_gap >= 0, "max_gap must be positive (or zero)"

        self.window_size = window_size
        self.max_gap = max_gap

        self.rms = None
        self._rigid_motion = None
        self.refcoord = None
        self._coord = None
        self._superimposer = QCPSuperimposer()

    def get_guide_coord_from_structure(self, structure):
        """Return the coordinates of guide atoms in the structure.

        We use guide atoms (C-alpha and C4' atoms) since it is much faster than
        using all atoms in the calculation without a significant loss in
        accuracy.
        """
        coords = []
        # CE algorithm is sensitive to atom ordering. To reproduce Pymol
        # results, sort atoms by chain and then residue number.
        for chain in sorted(structure.get_chains()):
            for resid in sorted(chain, key=_RESID_SORTER):
                if "CA" in resid:
                    coords.append(resid["CA"].coord.tolist())
                elif "C4'" in resid:
                    coords.append(resid["C4'"].coord.tolist())
        if not coords:
            msg = f"Structure {structure.id} does not have any guide atoms."
            raise PDBException(msg)
        return coords

    def set_reference(self, structure):
        """Define a reference structure onto which all others will be aligned."""
        self.refcoord = self.get_guide_coord_from_structure(structure)

        if len(self.refcoord) < self.window_size * 2:
            n_atoms = len(self.refcoord)
            msg = (
                f"Too few atoms in the reference structure ({n_atoms}). "
                "Try reducing the window_size parameter."
            )
            raise PDBException(msg)

    def align(self, structure, transform=True, *, final_optimization=True):
        """Align the input structure onto the reference structure.

        Parameters
        ----------
        transform: bool, optional
            If True (default), apply the rotation/translation that minimizes
            the RMSD between the two structures to the input structure. If
            False, the structure is not modified but the optimal RMSD will
            still be calculated.
        final_optimization: bool, optional
            If True (default), apply additional optimization to statistically
            significant alignments.
        """
        self.rms = None  # clear before aligning
        self._rigid_motion = None

        coord = self.get_guide_coord_from_structure(structure)
        self._coord = coord

        if len(coord) < self.window_size * 2:
            n_atoms = len(coord)
            msg = (
                f"Too few atoms in the mobile structure ({n_atoms}). "
                "Try reducing the window_size parameter."
            )
            raise PDBException(msg)

        # Run CEAlign
        # CEAlign returns the best N paths, sorted descending by length,
        # where each path is a pair of lists with aligned atom indices.
        alignments = run_cealign(self.refcoord, coord, self.window_size, self.max_gap)
        longest_alignments = [
            alignment
            for alignment in alignments
            if alignment.length == alignments[0].length
        ]

        # Iterate over paths and find the one that gives the lowest
        # corresponding RMSD. Use QCP to align the molecules.
        self.rms = float("inf")
        superimposer = self._superimposer
        best_alignment = None
        for alignment in longest_alignments:
            idxA, idxB = alignment.path

            coordsA = np.array([self.refcoord[i] for i in idxA])
            coordsB = np.array([coord[i] for i in idxB])

            superimposer.set(coordsA, coordsB)
            superimposer.run()
            if superimposer.rms < self.rms:
                best_alignment = alignment
                self.rms = superimposer.rms
                self._rigid_motion = (superimposer.rot, superimposer.tran)

        if best_alignment is None:
            raise RuntimeError("Failed to find a suitable alignment.")

        # Gap optimization
        if final_optimization and best_alignment.z_score >= 3.5:
            self._optimize(best_alignment)

        if transform:
            # Transform all atoms
            rotmtx, trvec = self._rigid_motion
            for chain in structure.get_chains():
                for resid in chain.get_unpacked_list():
                    for atom in resid.get_unpacked_list():
                        atom.transform(rotmtx, trvec)

    def _optimize(self, alignment):
        best_path = alignment.path
        superimposer = self._superimposer

        for ab_index in [0, 1]:
            for index in range(1, len(best_path[ab_index]) - 1):
                best_shift = 0
                left = best_path[ab_index][index - 1]
                center = best_path[ab_index][index]
                right = best_path[ab_index][index + 1]

                for shift in range(
                    max(-self.window_size // 2, left - center + 1),
                    min(self.window_size // 2 + 1, right - center),
                ):
                    best_path[ab_index][index] += shift
                    idxA, idxB = best_path

                    coordsA = np.array([self.refcoord[i] for i in idxA])
                    coordsB = np.array([self._coord[i] for i in idxB])

                    superimposer.set(coordsA, coordsB)
                    superimposer.run()
                    best_path[ab_index][index] -= shift

                    if superimposer.rms < self.rms:
                        best_shift = shift
                        self.rms = superimposer.rms
                        self._rigid_motion = (superimposer.rot, superimposer.tran)

                best_path[ab_index][index] += best_shift
