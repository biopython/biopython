# Copyright (C) 2020, Joao Rodrigues (j.p.g.l.m.rodrigues@gmail.com)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Calculation of solvent accessible surface areas for Bio.PDB entities.

Uses the "rolling ball" algorithm developed by Shrake & Rupley algorithm,
which uses a sphere (of equal radius to a solvent molecule) to probe the
surface of the molecule.

Reference:
Shrake, A; Rupley, JA. (1973).
"Environment and exposure to solvent of protein atoms. Lysozyme and insulin".
J Mol Biol
"""

import collections
import math

import numpy as np

from Bio.PDB.kdtrees import KDTree

# vdW radii taken from:
# https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
#
# References:
# A. Bondi (1964). "van der Waals Volumes and Radii".
ATOMIC_RADII = collections.defaultdict(lambda x: 2.0)
ATOMIC_RADII.update(
    {
        "H": 1.200,
        "C": 1.700,
        "N": 1.550,
        "O": 1.520,
        "S": 1.800,
        "P": 1.800,
    }
)


class ShrakeRupley:
    """Calculates SASAs using the Shrake-Rupley algorithm."""

    def __init__(self, probe_radius=1.40, n_points=100):
        """Initialize the class.

        :param probe_radius: radius of the probe in A.
        :type probe_radius: float

        :param n_points: resolution of the surface of each atom.
        :type n_points: int
        """
        if probe_radius <= 0.0:
            raise ValueError(
                f"Probe radius must be a positive number: {probe_radius} <= 0"
            )

        self.probe_radius = float(probe_radius)

        if n_points < 1:
            raise ValueError(
                f"Number of sphere points must be larger than 1: {n_points}"
            )
        self.n_points = n_points

        # Pre-compute reference sphere
        self._sphere = self._compute_sphere()

    def _compute_sphere(self):
        """Return the 3D coordinates of n points on a sphere.

        Uses the golden spiral algorithm to place points 'evenly'
        on the sphere surface. We compute this once and then move
        the sphere to the centroid of each atom as we compute the
        ASAs.
        """
        n = self.n_points

        dl = np.pi * (3 - 5 ** 0.5)
        dz = 2.0 / n

        longitude = 0
        z = 1 - dz / 2

        coords = np.zeros((n, 3), dtype=np.float32)
        for k in range(n):
            r = (1 - z * z) ** 0.5
            coords[k, 0] = math.cos(longitude) * r
            coords[k, 1] = math.sin(longitude) * r
            coords[k, 2] = z
            z -= dz
            longitude += dl

        return coords

    def compute(self, entity):
        """Calculate atomic surface accessibility values for a given entity.

        :param entity: input entity.
        :type entity: Bio.PDB.Entity.
        """
        # Get atoms onto list for lookup
        atoms = list(entity.get_atoms())
        n_atoms = len(atoms)

        # Get coordinates as a numpy array
        # We trust DisorderedAtom and friends to pick representatives.
        coords = np.array([a.coord for a in atoms], dtype=np.float64)

        # Pre-compute atom neighbors using KDTree
        kdt = KDTree(coords, 10)

        # Pre-compute radius * probe table
        radii = np.array([ATOMIC_RADII[a.element] for a in atoms], dtype=np.float64)
        radii += self.probe_radius
        twice_maxradii = max(radii) * 2

        # Calculate ASAs
        asa_array = np.zeros((n_atoms, 1), dtype=np.int)
        for i in range(n_atoms):

            r_i = radii[i]

            # Move sphere to atom
            s_on_i = (np.array(self._sphere, copy=True) * r_i) + coords[i]
            available_set = set(range(self.n_points))

            # KDtree for sphere points
            kdt_sphere = KDTree(s_on_i, 10)

            # Iterate over neighbors of atom i
            for jj in kdt.search(coords[i], twice_maxradii):
                j = jj.index
                if i == j:
                    continue

                if jj.radius < (r_i + radii[j]):
                    # Remove overlapping points on sphere from available set
                    available_set -= {
                        pt.index for pt in kdt_sphere.search(coords[j], radii[j])
                    }

            asa_array[i] = len(available_set)  # update counts

        # Convert accessible point count to surface area in A**2
        f = radii * radii * (4 * np.pi / self.n_points)
        asa_array = asa_array * f[:, np.newaxis]

        # Set atom attribute
        for i, atom in enumerate(atoms):
            atom.sasa = asa_array[i, 0]
            # atom.xtra["SASA"] = asa_array[i, 0]

        return asa_array
