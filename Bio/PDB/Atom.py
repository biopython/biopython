# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Atom class, used in Structure objects."""

import copy
import sys
import warnings

import numpy as np

from Bio.PDB.Entity import DisorderedEntityWrapper
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.PDB.vectors import Vector
from Bio.Data import IUPACData


class Atom:
    """Define Atom class.

    The Atom object stores atom name (both with and without spaces),
    coordinates, B factor, occupancy, alternative location specifier
    and (optionally) anisotropic B factor and standard deviations of
    B factor and positions.

    In the case of PQR files, B factor and occupancy are replaced by
    atomic charge and radius.
    """

    def __init__(
        self,
        name,
        coord,
        bfactor,
        occupancy,
        altloc,
        fullname,
        serial_number,
        element=None,
        pqr_charge=None,
        radius=None,
    ):
        """Initialize Atom object.

        :param name: atom name (eg. "CA"). Note that spaces are normally stripped.
        :type name: string

        :param coord: atomic coordinates (x,y,z)
        :type coord: NumPy array (Float0, length 3)

        :param bfactor: isotropic B factor
        :type bfactor: number

        :param occupancy: occupancy (0.0-1.0)
        :type occupancy: number

        :param altloc: alternative location specifier for disordered atoms
        :type altloc: string

        :param fullname: full atom name, including spaces, e.g. " CA ". Normally
                         these spaces are stripped from the atom name.
        :type fullname: string

        :param element: atom element, e.g. "C" for Carbon, "HG" for mercury,
        :type element: uppercase string (or None if unknown)

        :param pqr_charge: atom charge
        :type pqr_charge: number

        :param radius: atom radius
        :type radius: number
        """
        self.level = "A"
        # Reference to the residue
        self.parent = None
        # the atomic data
        self.name = name  # eg. CA, spaces are removed from atom name
        self.fullname = fullname  # e.g. " CA ", spaces included
        self.coord = coord
        self.bfactor = bfactor
        self.occupancy = occupancy
        self.altloc = altloc
        self.full_id = None  # (structure id, model id, chain id, residue id, atom id)
        self.id = name  # id of atom is the atom name (e.g. "CA")
        self.disordered_flag = 0
        self.anisou_array = None
        self.siguij_array = None
        self.sigatm_array = None
        self.serial_number = serial_number
        # Dictionary that keeps additional properties
        self.xtra = {}
        assert not element or element == element.upper(), element
        self.element = self._assign_element(element)
        self.mass = self._assign_atom_mass()
        self.pqr_charge = pqr_charge
        self.radius = radius

        # For atom sorting (protein backbone atoms first)
        self._sorting_keys = {"N": 0, "CA": 1, "C": 2, "O": 3}

    # Sorting Methods
    # standard across different objects and allows direct comparison
    def __eq__(self, other):
        """Test equality."""
        if isinstance(other, Atom):
            return self.full_id[1:] == other.full_id[1:]
        else:
            return NotImplemented

    def __ne__(self, other):
        """Test inequality."""
        if isinstance(other, Atom):
            return self.full_id[1:] != other.full_id[1:]
        else:
            return NotImplemented

    def __gt__(self, other):
        """Test greater than."""
        if isinstance(other, Atom):
            if self.parent != other.parent:
                return self.parent > other.parent
            order_s = self._sorting_keys.get(self.name, 4)
            order_o = self._sorting_keys.get(other.name, 4)
            if order_s != order_o:
                return order_s > order_o
            elif self.name != other.name:
                return self.name > other.name
            else:
                return self.altloc > other.altloc
        else:
            return NotImplemented

    def __ge__(self, other):
        """Test greater or equal."""
        if isinstance(other, Atom):
            if self.parent != other.parent:
                return self.parent >= other.parent
            order_s = self._sorting_keys.get(self.name, 4)
            order_o = self._sorting_keys.get(other.name, 4)
            if order_s != order_o:
                return order_s >= order_o
            elif self.name != other.name:
                return self.name >= other.name
            else:
                return self.altloc >= other.altloc
        else:
            return NotImplemented

    def __lt__(self, other):
        """Test less than."""
        if isinstance(other, Atom):
            if self.parent != other.parent:
                return self.parent < other.parent
            order_s = self._sorting_keys.get(self.name, 4)
            order_o = self._sorting_keys.get(other.name, 4)
            if order_s != order_o:
                return order_s < order_o
            elif self.name != other.name:
                return self.name < other.name
            else:
                return self.altloc < other.altloc
        else:
            return NotImplemented

    def __le__(self, other):
        """Test less or equal."""
        if isinstance(other, Atom):
            if self.parent != other.parent:
                return self.parent <= other.parent
            order_s = self._sorting_keys.get(self.name, 4)
            order_o = self._sorting_keys.get(other.name, 4)
            if order_s != order_o:
                return order_s <= order_o
            elif self.name != other.name:
                return self.name <= other.name
            else:
                return self.altloc <= other.altloc
        else:
            return NotImplemented

    # Hash method to allow uniqueness (set)
    def __hash__(self):
        """Return atom full identifier."""
        return hash(self.get_full_id())

    def _assign_element(self, element):
        """Guess element from atom name if not recognised (PRIVATE).

        There is little documentation about extracting/encoding element
        information in atom names, but some conventions seem to prevail:

            - C, N, O, S, H, P, F atom names start with a blank space (e.g. " CA ")
              unless the name is 4 characters long (e.g. HE21 in glutamine). In both
              these cases, the element is the first character.

            - Inorganic elements do not have a blank space (e.g. "CA  " for calcium)
              but one must check the full name to differentiate between e.g. helium
              ("HE  ") and long-name hydrogens (e.g. "HE21").

            - Atoms with unknown or ambiguous elements are marked with 'X', e.g.
              PDB 4cpa. If we fail to identify an element, we should mark it as
              such.

        """
        if not element or element.capitalize() not in IUPACData.atom_weights:
            if self.fullname[0].isalpha() and not self.fullname[2:].isdigit():
                putative_element = self.name.strip()
            else:
                # Hs may have digit in [0]
                if self.name[0].isdigit():
                    putative_element = self.name[1]
                else:
                    putative_element = self.name[0]

            if putative_element.capitalize() in IUPACData.atom_weights:
                msg = "Used element %r for Atom (name=%s) with given element %r" % (
                    putative_element,
                    self.name,
                    element,
                )
                element = putative_element
            else:
                msg = (
                    "Could not assign element %r for Atom (name=%s) with given element %r"
                    % (putative_element, self.name, element)
                )
                element = "X"  # mark as unknown/ambiguous
            warnings.warn(msg, PDBConstructionWarning)

        return element

    def _assign_atom_mass(self):
        """Return atom weight (PRIVATE)."""
        try:
            return IUPACData.atom_weights[self.element.capitalize()]
        except (AttributeError, KeyError):
            return float("NaN")

    # Special methods

    def __repr__(self):
        """Print Atom object as <Atom atom_name>."""
        return f"<Atom {self.get_id()}>"

    def __sub__(self, other):
        """Calculate distance between two atoms.

        :param other: the other atom
        :type other: L{Atom}

        Examples
        --------
        This is an incomplete but illustrative example::

            distance = atom1 - atom2

        """
        diff = self.coord - other.coord
        return np.sqrt(np.dot(diff, diff))

    # set methods

    def set_serial_number(self, n):
        """Set serial number."""
        self.serial_number = n

    def set_bfactor(self, bfactor):
        """Set isotroptic B factor."""
        self.bfactor = bfactor

    def set_coord(self, coord):
        """Set coordinates."""
        self.coord = coord

    def set_altloc(self, altloc):
        """Set alternative location specifier."""
        self.altloc = altloc

    def set_occupancy(self, occupancy):
        """Set occupancy."""
        self.occupancy = occupancy

    def set_sigatm(self, sigatm_array):
        """Set standard deviation of atomic parameters.

        The standard deviation of atomic parameters consists
        of 3 positional, 1 B factor and 1 occupancy standard
        deviation.

        :param sigatm_array: standard deviations of atomic parameters.
        :type sigatm_array: NumPy array (length 5)
        """
        self.sigatm_array = sigatm_array

    def set_siguij(self, siguij_array):
        """Set standard deviations of anisotropic temperature factors.

        :param siguij_array: standard deviations of anisotropic temperature factors.
        :type siguij_array: NumPy array (length 6)
        """
        self.siguij_array = siguij_array

    def set_anisou(self, anisou_array):
        """Set anisotropic B factor.

        :param anisou_array: anisotropic B factor.
        :type anisou_array: NumPy array (length 6)
        """
        self.anisou_array = anisou_array

    def set_charge(self, pqr_charge):
        """Set charge."""
        self.pqr_charge = pqr_charge

    def set_radius(self, radius):
        """Set radius."""
        self.radius = radius

    # Public methods

    def flag_disorder(self):
        """Set the disordered flag to 1.

        The disordered flag indicates whether the atom is disordered or not.
        """
        self.disordered_flag = 1

    def is_disordered(self):
        """Return the disordered flag (1 if disordered, 0 otherwise)."""
        return self.disordered_flag

    def set_parent(self, parent):
        """Set the parent residue.

        Arguments:
         - parent - Residue object

        """
        self.parent = parent
        self.full_id = self.get_full_id()

    def detach_parent(self):
        """Remove reference to parent."""
        self.parent = None

    def get_sigatm(self):
        """Return standard deviation of atomic parameters."""
        return self.sigatm_array

    def get_siguij(self):
        """Return standard deviations of anisotropic temperature factors."""
        return self.siguij_array

    def get_anisou(self):
        """Return anisotropic B factor."""
        return self.anisou_array

    def get_parent(self):
        """Return parent residue."""
        return self.parent

    def get_serial_number(self):
        """Return the serial number."""
        return self.serial_number

    def get_name(self):
        """Return atom name."""
        return self.name

    def get_id(self):
        """Return the id of the atom (which is its atom name)."""
        return self.id

    def get_full_id(self):
        """Return the full id of the atom.

        The full id of an atom is a tuple used to uniquely identify
        the atom and consists of the following elements:
        (structure id, model id, chain id, residue id, atom name, altloc)
        """
        try:
            return self.parent.get_full_id() + ((self.name, self.altloc),)
        except AttributeError:
            return (None, None, None, None, self.name, self.altloc)

    def get_coord(self):
        """Return atomic coordinates."""
        return self.coord

    def get_bfactor(self):
        """Return B factor."""
        return self.bfactor

    def get_occupancy(self):
        """Return occupancy."""
        return self.occupancy

    def get_fullname(self):
        """Return the atom name, including leading and trailing spaces."""
        return self.fullname

    def get_altloc(self):
        """Return alternative location specifier."""
        return self.altloc

    def get_level(self):
        """Return level."""
        return self.level

    def get_charge(self):
        """Return charge."""
        return self.pqr_charge

    def get_radius(self):
        """Return radius."""
        return self.radius

    def transform(self, rot, tran):
        """Apply rotation and translation to the atomic coordinates.

        :param rot: A right multiplying rotation matrix
        :type rot: 3x3 NumPy array

        :param tran: the translation vector
        :type tran: size 3 NumPy array

        Examples
        --------
        This is an incomplete but illustrative example::

            from numpy import pi, array
            from Bio.PDB.vectors import Vector, rotmat
            rotation = rotmat(pi, Vector(1, 0, 0))
            translation = array((0, 0, 1), 'f')
            atom.transform(rotation, translation)

        """
        self.coord = np.dot(self.coord, rot) + tran

    def get_vector(self):
        """Return coordinates as Vector.

        :return: coordinates as 3D vector
        :rtype: Bio.PDB.Vector class
        """
        x, y, z = self.coord
        return Vector(x, y, z)

    def copy(self):
        """Create a copy of the Atom.

        Parent information is lost.
        """
        # Do a shallow copy then explicitly copy what needs to be deeper.
        shallow = copy.copy(self)
        shallow.detach_parent()
        shallow.set_coord(copy.copy(self.get_coord()))
        shallow.xtra = self.xtra.copy()
        return shallow


class DisorderedAtom(DisorderedEntityWrapper):
    """Contains all Atom objects that represent the same disordered atom.

    One of these atoms is "selected" and all method calls not caught
    by DisorderedAtom are forwarded to the selected Atom object. In that way, a
    DisorderedAtom behaves exactly like a normal Atom. By default, the selected
    Atom object represents the Atom object with the highest occupancy, but a
    different Atom object can be selected by using the disordered_select(altloc)
    method.
    """

    def __init__(self, id):
        """Create DisorderedAtom.

        Arguments:
         - id - string, atom name

        """
        # TODO - make this a private attribute?
        self.last_occupancy = -sys.maxsize
        DisorderedEntityWrapper.__init__(self, id)

    # Special methods
    # Override parent class __iter__ method
    def __iter__(self):
        """Iterate through disordered atoms."""
        yield from self.disordered_get_list()

    def __repr__(self):
        """Return disordered atom identifier."""
        if self.child_dict:
            return f"<DisorderedAtom {self.get_id()}>"
        else:
            return f"<Empty DisorderedAtom {self.get_id()}>"

    # This is a separate method from Entity.center_of_mass since DisorderedAtoms
    # will be unpacked by Residue.get_unpacked_list(). Here we allow for a very
    # specific use case that is much simpler than the general implementation.
    def center_of_mass(self):
        """Return the center of mass of the DisorderedAtom as a numpy array.

        Assumes all child atoms have the same mass (same element).
        """
        children = self.disordered_get_list()

        if not children:
            raise ValueError(f"{self} does not have children")

        coords = np.asarray([a.coord for a in children], dtype=np.float32)
        return np.average(coords, axis=0, weights=None)

    def disordered_get_list(self):
        """Return list of atom instances.

        Sorts children by altloc (empty, then alphabetical).
        """
        return sorted(self.child_dict.values(), key=lambda a: ord(a.altloc))

    def disordered_add(self, atom):
        """Add a disordered atom."""
        # Add atom to dict, use altloc as key
        atom.flag_disorder()
        # set the residue parent of the added atom
        residue = self.get_parent()
        atom.set_parent(residue)
        altloc = atom.get_altloc()
        occupancy = atom.get_occupancy()
        self[altloc] = atom
        if occupancy > self.last_occupancy:
            self.last_occupancy = occupancy
            self.disordered_select(altloc)

    def disordered_remove(self, altloc):
        """Remove a child atom altloc from the DisorderedAtom.

        Arguments:
         - altloc - name of the altloc to remove, as a string.

        """
        # Get child altloc
        atom = self.child_dict[altloc]
        is_selected = self.selected_child is atom

        # Detach
        del self.child_dict[altloc]
        atom.detach_parent()

        if is_selected and self.child_dict:  # pick next highest occupancy
            child = sorted(self.child_dict.values(), key=lambda a: a.occupancy)[-1]
            self.disordered_select(child.altloc)
        elif not self.child_dict:
            self.selected_child = None
            self.last_occupancy = -sys.maxsize

    def transform(self, rot, tran):
        """Apply rotation and translation to all children.

        See the documentation of Atom.transform for details.
        """
        for child in self:
            child.coord = np.dot(child.coord, rot) + tran
