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

    def align(self, structure, transform=True):
        """Align the input structure onto the reference structure.

        Parameters
        ----------
        transform: bool, optional
            If True (default), apply the rotation/translation that minimizes
            the RMSD between the two structures to the input structure. If
            False, the structure is not modified but the optimal RMSD will
            still be calculated.
        """
        self.rms = None  # clear before aligning

        coord = self.get_guide_coord_from_structure(structure)

        if len(coord) < self.window_size * 2:
            n_atoms = len(coord)
            msg = (
                f"Too few atoms in the mobile structure ({n_atoms}). "
                "Try reducing the window_size parameter."
            )
            raise PDBException(msg)

        # Run CEAlign
        # CEAlign returns the best N paths, where each path is a pair of lists
        # with aligned atom indices. Paths are not guaranteed to be unique.
        paths = run_cealign(self.refcoord, coord, self.window_size, self.max_gap)
        unique_paths = {(tuple(pA), tuple(pB)) for pA, pB in paths}

        # Iterate over unique paths and find the one that gives the lowest
        # corresponding RMSD. Use QCP to align the molecules.
        best_rmsd, best_u = 1e6, None
        for u_path in unique_paths:
            idxA, idxB = u_path

            coordsA = np.array([self.refcoord[i] for i in idxA])
            coordsB = np.array([coord[i] for i in idxB])

            aln = QCPSuperimposer()
            aln.set(coordsA, coordsB)
            aln.run()
            if aln.rms < best_rmsd:
                best_rmsd = aln.rms
                best_u = (aln.rot, aln.tran)

        if best_u is None:
            raise RuntimeError("Failed to find a suitable alignment.")

        if transform:
            # Transform all atoms
            rotmtx, trvec = best_u
            for chain in structure.get_chains():
                for resid in chain.get_unpacked_list():
                    for atom in resid.get_unpacked_list():
                        atom.transform(rotmtx, trvec)

        self.rms = best_rmsd
