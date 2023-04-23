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
    Shrake, A; Rupley, JA. (1973). J Mol Biol
    "Environment and exposure to solvent of protein atoms. Lysozyme and insulin".
"""

import collections
import math

import numpy as np

from Bio.PDB.kdtrees import KDTree

__all__ = ["ShrakeRupley"]

_ENTITY_HIERARCHY = {
    "A": 0,
    "R": 1,
    "C": 2,
    "M": 3,
    "S": 4,
}

# vdW radii taken from:
# https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
#
# Radii for CL, K, NA, etc are _not_ ionic radii.
#
# References:
# A. Bondi (1964). "van der Waals Volumes and Radii".
# M. Mantina, A.C. et al., J. Phys. Chem. 2009, 113, 5806.
ATOMIC_RADII = collections.defaultdict(lambda: 2.0)
ATOMIC_RADII.update(
    {
        "H": 1.200,
        "HE": 1.400,
        "C": 1.700,
        "N": 1.550,
        "O": 1.520,
        "F": 1.470,
        "NA": 2.270,
        "MG": 1.730,
        "P": 1.800,
        "S": 1.800,
        "CL": 1.750,
        "K": 2.750,
        "CA": 2.310,
        "NI": 1.630,
        "CU": 1.400,
        "ZN": 1.390,
        "SE": 1.900,
        "BR": 1.850,
        "CD": 1.580,
        "I": 1.980,
        "HG": 1.550,
    }
)


class ShrakeRupley:
    """Calculates SASAs using the Shrake-Rupley algorithm."""

    def __init__(self, probe_radius=1.40, n_points=100, radii_dict=None):
        """Initialize the class.

        :param probe_radius: radius of the probe in A. Default is 1.40, roughly
            the radius of a water molecule.
        :type probe_radius: float

        :param n_points: resolution of the surface of each atom. Default is 100.
            A higher number of points results in more precise measurements, but
            slows down the calculation.
        :type n_points: int

        :param radii_dict: user-provided dictionary of atomic radii to use in
            the calculation. Values will replace/complement those in the
            default ATOMIC_RADII dictionary.
        :type radii_dict: dict

        >>> sr = ShrakeRupley()
        >>> sr = ShrakeRupley(n_points=960)
        >>> sr = ShrakeRupley(radii_dict={"O": 3.1415})
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

        # Update radii list with user provided lists.
        self.radii_dict = ATOMIC_RADII.copy()
        if radii_dict is not None:
            self.radii_dict.update(radii_dict)

        # Pre-compute reference sphere
        self._sphere = self._compute_sphere()

    def _compute_sphere(self):
        """Return the 3D coordinates of n points on a sphere.

        Uses the golden spiral algorithm to place points 'evenly' on the sphere
        surface. We compute this once and then move the sphere to the centroid
        of each atom as we compute the ASAs.
        """
        n = self.n_points

        dl = np.pi * (3 - 5**0.5)
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

    def compute(self, entity, level="A"):
        """Calculate surface accessibility surface area for an entity.

        The resulting atomic surface accessibility values are attached to the
        .sasa attribute of each entity (or atom), depending on the level. For
        example, if level="R", all residues will have a .sasa attribute. Atoms
        will always be assigned a .sasa attribute with their individual values.

        :param entity: input entity.
        :type entity: Bio.PDB.Entity, e.g. Residue, Chain, ...

        :param level: the level at which ASA values are assigned, which can be
            one of "A" (Atom), "R" (Residue), "C" (Chain), "M" (Model), or
            "S" (Structure). The ASA value of an entity is the sum of all ASA
            values of its children. Defaults to "A".
        :type entity: Bio.PDB.Entity

        >>> from Bio.PDB import PDBParser
        >>> from Bio.PDB.SASA import ShrakeRupley
        >>> p = PDBParser(QUIET=1)
        >>> # This assumes you have a local copy of 1LCD.pdb in a directory called "PDB"
        >>> struct = p.get_structure("1LCD", "PDB/1LCD.pdb")
        >>> sr = ShrakeRupley()
        >>> sr.compute(struct, level="S")
        >>> print(round(struct.sasa, 2))
        7053.43
        >>> print(round(struct[0]["A"][11]["OE1"].sasa, 2))
        9.64
        """
        is_valid = hasattr(entity, "level") and entity.level in {"R", "C", "M", "S"}
        if not is_valid:
            raise ValueError(
                f"Invalid entity type '{type(entity)}'. "
                "Must be Residue, Chain, Model, or Structure"
            )

        if level not in _ENTITY_HIERARCHY:
            raise ValueError(f"Invalid level '{level}'. Must be A, R, C, M, or S.")
        elif _ENTITY_HIERARCHY[level] > _ENTITY_HIERARCHY[entity.level]:
            raise ValueError(
                f"Level '{level}' must be equal or smaller than input entity: {entity.level}"
            )

        # Get atoms onto list for lookup
        atoms = list(entity.get_atoms())
        n_atoms = len(atoms)
        if not n_atoms:
            raise ValueError("Entity has no child atoms.")

        # Get coordinates as a numpy array
        # We trust DisorderedAtom and friends to pick representatives.
        coords = np.array([a.coord for a in atoms], dtype=np.float64)

        # Pre-compute atom neighbors using KDTree
        kdt = KDTree(coords, 10)

        # Pre-compute radius * probe table
        radii_dict = self.radii_dict
        radii = np.array([radii_dict[a.element] for a in atoms], dtype=np.float64)
        radii += self.probe_radius
        twice_maxradii = np.max(radii) * 2

        # Calculate ASAs
        asa_array = np.zeros((n_atoms, 1), dtype=np.int64)
        ptset = set(range(self.n_points))
        for i in range(n_atoms):

            r_i = radii[i]

            # Move sphere to atom
            s_on_i = (np.array(self._sphere, copy=True) * r_i) + coords[i]
            available_set = ptset.copy()

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

        # Set atom .sasa
        for i, atom in enumerate(atoms):
            atom.sasa = asa_array[i, 0]

        # Aggregate values per entity level if necessary
        if level != "A":
            entities = set(atoms)
            target = _ENTITY_HIERARCHY[level]
            for _ in range(target):
                entities = {e.parent for e in entities}

            atomdict = {a.full_id: idx for idx, a in enumerate(atoms)}
            for e in entities:
                e_atoms = [atomdict[a.full_id] for a in e.get_atoms()]
                e.sasa = asa_array[e_atoms].sum()
